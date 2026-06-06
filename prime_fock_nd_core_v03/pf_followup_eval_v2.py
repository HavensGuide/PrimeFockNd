"""
pf_followup_eval_v2.py
=======================

Patch / extension for pf_followup_eval.py.

Adds the missing checklist items:
  - alternative nonlinear eval for PCA/random projection baselines
  - future volatility / abs-change / squared-change / dispersion targets
  - future correlation-change targets
  - regime ranking metrics: AUC, precision@top-decile, top-vs-bottom spread
  - kernel / representation alignment metrics: CKA, Gram error, pairwise-distance rank corr
  - wave conservation diagnostics: diffusion mass/positivity, unitary norm, damped norm decay
  - graph-change forecasting targets from cached graph bundles

This script reuses pf_bench_cache and the graph/descriptor builders from
pf_followup_eval.py and pf_better_benchmark.py. It does not require sklearn.

Run examples:
  python pf_followup_eval_v2.py --preset lite --datasets plotly_stock_returns --only-cached-graphs --max-pf-configs 2
  python pf_followup_eval_v2.py --preset lite --only-cached-graphs
  python pf_followup_eval_v2.py --preset markets --only-cached-graphs
"""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np
import pandas as pd

import pf_better_benchmark as base
import pf_followup_eval as fup


# -----------------------------------------------------------------------------
# Metrics and target helpers
# -----------------------------------------------------------------------------

def _as_2d(y: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    return np.nan_to_num(y, nan=0.0, posinf=0.0, neginf=0.0)


def rankdata_average(x: np.ndarray) -> np.ndarray:
    """Average ranks, pure NumPy. Good enough for diagnostics."""
    x = np.asarray(x, dtype=float)
    order = np.argsort(x)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(len(x), dtype=float)
    # Tie handling: average equal-value ranks.
    xs = x[order]
    i = 0
    while i < len(x):
        j = i + 1
        while j < len(x) and xs[j] == xs[i]:
            j += 1
        if j - i > 1:
            avg = 0.5 * (i + j - 1)
            ranks[order[i:j]] = avg
        i = j
    return ranks


def spearman_corr(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size < 3 or np.std(x) <= 1e-12 or np.std(y) <= 1e-12:
        return np.nan
    rx = rankdata_average(x)
    ry = rankdata_average(y)
    return float(np.corrcoef(rx, ry)[0, 1])


def binary_auc(scores: np.ndarray, labels: np.ndarray) -> float:
    scores = np.asarray(scores, dtype=float).ravel()
    labels = np.asarray(labels, dtype=float).ravel() > 0.5
    mask = np.isfinite(scores) & np.isfinite(labels)
    scores = scores[mask]
    labels = labels[mask]
    n_pos = int(labels.sum())
    n_neg = int((~labels).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = rankdata_average(scores)
    sum_pos = ranks[labels].sum()
    auc = (sum_pos - n_pos * (n_pos - 1) / 2.0) / (n_pos * n_neg)
    return float(auc)


def ranking_metrics(y_true: np.ndarray, score: np.ndarray, top_frac: float = 0.10) -> dict[str, float]:
    y = np.asarray(y_true, dtype=float).ravel()
    s = np.asarray(score, dtype=float).ravel()
    mask = np.isfinite(y) & np.isfinite(s)
    y = y[mask]
    s = s[mask]
    if y.size < 10:
        return {"rank_spearman": np.nan, "precision_top_decile": np.nan, "top_bottom_target_spread": np.nan, "auc_top_quartile": np.nan}
    n_top = max(1, int(np.ceil(top_frac * len(y))))
    order = np.argsort(s)
    top = order[-n_top:]
    bottom = order[:n_top]
    # Regime label: top quartile of realized target.
    q75 = np.quantile(y, 0.75)
    labels = (y >= q75).astype(float)
    return {
        "rank_spearman": spearman_corr(y, s),
        "precision_top_decile": float(labels[top].mean()),
        "top_bottom_target_spread": float(y[top].mean() - y[bottom].mean()),
        "auc_top_quartile": binary_auc(s, labels),
    }


def centered_kernel_alignment(X: np.ndarray, Z: np.ndarray) -> float:
    X = np.nan_to_num(np.asarray(X, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    Z = np.nan_to_num(np.asarray(Z, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    if X.shape[0] < 3 or Z.shape[0] != X.shape[0]:
        return np.nan
    Kx = X @ X.T
    Kz = Z @ Z.T
    H = np.eye(Kx.shape[0]) - np.ones_like(Kx) / Kx.shape[0]
    Kxc = H @ Kx @ H
    Kzc = H @ Kz @ H
    denom = np.linalg.norm(Kxc, "fro") * np.linalg.norm(Kzc, "fro")
    if denom <= 1e-12:
        return np.nan
    return float(np.sum(Kxc * Kzc) / denom)


def gram_fro_error(X: np.ndarray, Z: np.ndarray) -> float:
    X = np.nan_to_num(np.asarray(X, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    Z = np.nan_to_num(np.asarray(Z, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    if X.shape[0] < 3 or Z.shape[0] != X.shape[0]:
        return np.nan
    Kx = X @ X.T
    Kz = Z @ Z.T
    # Normalize scale before comparing.
    nx = np.linalg.norm(Kx, "fro") + 1e-12
    nz = np.linalg.norm(Kz, "fro") + 1e-12
    return float(np.linalg.norm(Kx / nx - Kz / nz, "fro"))


def pairwise_distance_rank_corr(X: np.ndarray, Z: np.ndarray, max_pairs: int = 20000, seed: int = 0) -> float:
    X = np.nan_to_num(np.asarray(X, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    Z = np.nan_to_num(np.asarray(Z, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    n = X.shape[0]
    if n < 4 or Z.shape[0] != n:
        return np.nan
    rng = np.random.default_rng(seed)
    total = n * (n - 1) // 2
    m = min(max_pairs, total)
    if total <= max_pairs:
        ii, jj = np.triu_indices(n, k=1)
    else:
        ii = rng.integers(0, n, size=m)
        jj = rng.integers(0, n, size=m)
        mask = ii != jj
        ii = ii[mask]
        jj = jj[mask]
    dx = np.linalg.norm(X[ii] - X[jj], axis=1)
    dz = np.linalg.norm(Z[ii] - Z[jj], axis=1)
    return spearman_corr(dx, dz)


def representation_alignment_row(dataset: str, window_id: int, horizon: int, descriptor: str, eval_method: str, k: int,
                                 X_test: np.ndarray, Z_test: np.ndarray, extra: dict[str, Any], seed: int) -> dict[str, Any]:
    return {
        "dataset": dataset,
        "window": int(window_id),
        "horizon": int(horizon),
        "descriptor": descriptor,
        "eval_method": eval_method,
        "k": int(k),
        "feature_dim_in": int(Z_test.shape[1]) if Z_test.ndim == 2 else 1,
        "alignment_cka": centered_kernel_alignment(X_test, Z_test),
        "alignment_gram_fro_error": gram_fro_error(X_test, Z_test),
        "alignment_pairwise_spearman": pairwise_distance_rank_corr(X_test, Z_test, seed=seed),
        **extra,
    }


def corr_mat(M: np.ndarray) -> np.ndarray:
    if M.shape[0] < 3:
        return np.zeros((M.shape[1], M.shape[1]))
    C = np.corrcoef(M, rowvar=False)
    return np.nan_to_num(C, nan=0.0, posinf=0.0, neginf=0.0)


def make_scalar_targets_for_indices(Xstd_all: np.ndarray, pair_indices: np.ndarray, horizon: int, vol_window: int = 5, corr_window: int = 20) -> dict[str, np.ndarray]:
    """Targets computed from a standardized full sequence. Targets are aligned to X_t rows."""
    X = np.asarray(Xstd_all, dtype=float)
    out = {
        "future_realized_volatility": [],
        "future_cross_sectional_dispersion": [],
        "future_corr_change": [],
    }
    n = X.shape[0]
    for t in pair_indices:
        t = int(t)
        future_start = min(n, t + 1)
        future_end = min(n, t + 1 + vol_window)
        future = X[future_start:future_end]
        if future.shape[0] < 2:
            future = X[max(0, t):min(n, t + 2)]
        out["future_realized_volatility"].append(float(np.sqrt(np.mean(future ** 2))) if future.size else 0.0)
        y_idx = min(n - 1, t + horizon)
        out["future_cross_sectional_dispersion"].append(float(np.std(X[y_idx])))
        past = X[max(0, t - corr_window + 1):t + 1]
        fut = X[t + 1:min(n, t + 1 + corr_window)]
        if past.shape[0] >= 3 and fut.shape[0] >= 3 and X.shape[1] >= 2:
            C0 = corr_mat(past)
            C1 = corr_mat(fut)
            out["future_corr_change"].append(float(np.linalg.norm(C1 - C0, ord="fro") / X.shape[1]))
        else:
            out["future_corr_change"].append(0.0)
    return {k: np.asarray(v, dtype=float)[:, None] for k, v in out.items()}


def graph_change_target(graph: Any, horizon: int) -> np.ndarray:
    A = np.asarray(graph.adjacency, dtype=float)
    n = A.shape[0]
    y = np.zeros(n, dtype=float)
    for s in range(n):
        j = min(n - 1, s + horizon)
        denom = A.shape[1] + 1e-12
        y[s] = float(np.linalg.norm(A[j] - A[s], ord="fro") / denom)
    return y[:, None]


# -----------------------------------------------------------------------------
# Evaluation stack for arbitrary targets
# -----------------------------------------------------------------------------

def eval_prediction_row(dataset: str, window_id: int, horizon: int, descriptor: str, eval_method: str, k: int,
                        Z_train: np.ndarray, Z_test: np.ndarray, X_train: np.ndarray, X_test: np.ndarray,
                        Y_train: np.ndarray, Y_test: np.ndarray, cfg: base.BenchConfig,
                        target_name: str, extra: dict[str, Any]) -> dict[str, Any]:
    Y_train = _as_2d(Y_train)
    Y_test = _as_2d(Y_test)
    Y_pred = base.ridge_fit_predict(Z_train, Y_train, Z_test, alpha=cfg.ridge_alpha)
    X_rec = base.ridge_fit_predict(Z_train, X_train, Z_test, alpha=cfg.ridge_alpha)
    mf = base.regression_metrics_matrix(Y_test, Y_pred)
    mr = base.regression_metrics_matrix(X_test, X_rec)
    row = {
        "dataset": dataset,
        "window": int(window_id),
        "horizon": int(horizon),
        "descriptor": descriptor,
        "eval_method": eval_method,
        "target_name": target_name,
        "k": int(k),
        "feature_dim_in": int(Z_train.shape[1]),
        "forecast_r2": mf["r2"],
        "forecast_rmse": mf["rmse"],
        "forecast_mae": mf["mae"],
        "forecast_corr": mf["mean_channel_corr"],
        "reconstruction_r2": mr["r2"],
        "reconstruction_rmse": mr["rmse"],
        "corr_matrix_error": base.corr_matrix_error(X_test, X_rec),
        "compression_ratio": float(k / X_train.shape[1]) if X_train.shape[1] else np.nan,
        **extra,
    }
    # Ranking/regime metrics for scalar targets.
    if Y_test.shape[1] == 1:
        row.update(ranking_metrics(Y_test[:, 0], Y_pred[:, 0]))
    return row


def evaluate_block_v2(dataset: str, window_id: int, horizon: int, descriptor_name: str,
                      F_train: np.ndarray, F_test: np.ndarray,
                      X_train: np.ndarray, X_test: np.ndarray,
                      target_map_train: dict[str, np.ndarray], target_map_test: dict[str, np.ndarray],
                      cfg: base.BenchConfig, k_values: Sequence[int], seed: int,
                      extra: dict[str, Any], include_alignment: bool = True) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    align_rows: list[dict[str, Any]] = []
    F_train = np.nan_to_num(np.asarray(F_train, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    F_test = np.nan_to_num(np.asarray(F_test, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)

    def add_evals(method: str, k: int, Ztr: np.ndarray, Zte: np.ndarray) -> None:
        for target_name, Ytr in target_map_train.items():
            Yte = target_map_test[target_name]
            rows.append(eval_prediction_row(dataset, window_id, horizon, descriptor_name, method, k,
                                            Ztr, Zte, X_train, X_test, Ytr, Yte, cfg, target_name, extra))
        if include_alignment:
            align_rows.append(representation_alignment_row(dataset, window_id, horizon, descriptor_name, method, k,
                                                           X_test, Zte, extra, seed=seed + 37 * int(k)))

    # direct ridge, no PCA.
    Ztr, Zte = fup.standardize_features(F_train, F_test)
    add_evals("direct_ridge", Ztr.shape[1], Ztr, Zte)

    for k in k_values:
        k_eff = min(max(1, int(k)), F_train.shape[1], max(1, F_train.shape[0] - 1))
        Ztr, Zte = base.svd_project_train_test(F_train, F_test, k_eff)
        add_evals("svd_ridge", k_eff, Ztr, Zte)

        # top-k is target-specific because selection uses Y_train.
        for target_name, Ytr_target in target_map_train.items():
            Ztr_top, Zte_top, idx = fup.select_topcorr(F_train, F_test, _as_2d(Ytr_target), k_eff)
            rows.append(eval_prediction_row(dataset, window_id, horizon, descriptor_name, "topcorr_ridge", k_eff,
                                            Ztr_top, Zte_top, X_train, X_test, Ytr_target, target_map_test[target_name],
                                            cfg, target_name, {**extra, "selected_features": int(len(idx))}))
            if include_alignment:
                align_rows.append(representation_alignment_row(dataset, window_id, horizon, descriptor_name,
                                                               f"topcorr_ridge__{target_name}", k_eff,
                                                               X_test, Zte_top, {**extra, "selected_features": int(len(idx))},
                                                               seed=seed + 101 * k_eff))

        Ztr, Zte = fup.random_tanh_features(F_train, F_test, width=k_eff, seed=seed + 1009 * k_eff)
        add_evals("elm_tanh_ridge", k_eff, Ztr, Zte)

        Ztr, Zte = fup.rff_rbf_features(F_train, F_test, width=k_eff, seed=seed + 2003 * k_eff)
        add_evals("rff_rbf_ridge", k_eff, Ztr, Zte)

    return rows, align_rows


def pca_random_baseline_blocks(X_train: np.ndarray, X_test: np.ndarray, cfg: base.BenchConfig, seed: int) -> list[tuple[str, np.ndarray, np.ndarray, dict[str, Any]]]:
    blocks = [("raw_lag", X_train, X_test, {"descriptor_family": "raw_lag"})]
    for k in cfg.k_values:
        k_eff = min(max(1, int(k)), X_train.shape[1], max(1, X_train.shape[0] - 1))
        Ptr, Pte = base.svd_project_train_test(X_train, X_test, k_eff)
        blocks.append((f"baseline_pca_k{k_eff}", Ptr, Pte, {"descriptor_family": "baseline_pca", "baseline_k": k_eff}))
        Rtr, Rte = base.random_project_train_test(X_train, X_test, k_eff, seed + 31 * k_eff)
        blocks.append((f"baseline_random_projection_k{k_eff}", Rtr, Rte, {"descriptor_family": "baseline_random_projection", "baseline_k": k_eff}))
    return blocks


# -----------------------------------------------------------------------------
# Conservation diagnostics
# -----------------------------------------------------------------------------

def wave_conservation_diagnostics(bundle: dict[str, Any], X: np.ndarray, cfg: base.BenchConfig, split: str, gamma_values: Sequence[float]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    graph = bundle["prop_graph_train"] if split == "train" else bundle["prop_graph_test"]
    diff = base.GraphDiffusionDynamics(generator=cfg.graph_generator, normalize_output=True)
    X = np.asarray(X, dtype=float)
    for tau in cfg.tau_values:
        diff_mass_err = []
        diff_min = []
        for s in range(X.shape[0]):
            rho0, _ = base.density_and_phase_from_state(X, s, cfg.use_phase_from_delta)
            op = diff.prepare(graph, s)
            rho = diff.evolve(op, rho0, float(tau))
            diff_mass_err.append(abs(float(np.sum(rho)) - 1.0))
            diff_min.append(float(np.min(rho)))
        rows.append({
            "split": split,
            "diagnostic": "diffusion_mass_positivity",
            "tau": float(tau),
            "gamma": 0.0,
            "mean_mass_error": float(np.mean(diff_mass_err)),
            "max_mass_error": float(np.max(diff_mass_err)),
            "min_density": float(np.min(diff_min)),
            "graph_key": bundle.get("key"),
            **bundle.get("key_obj", {}),
        })

    for gamma in gamma_values:
        wave = base.GraphWaveDynamics(generator=cfg.graph_generator, gamma=float(gamma))
        for tau in cfg.tau_values:
            norm_err = []
            norms = []
            residual_energy = []
            for s in range(X.shape[0]):
                rho0, psi0 = base.density_and_phase_from_state(X, s, cfg.use_phase_from_delta)
                op = wave.prepare(graph, s)
                psi = wave.evolve(op, psi0, float(tau))
                rho_w = np.abs(psi) ** 2
                nrm = float(np.sum(rho_w))
                norms.append(nrm)
                if gamma == 0.0:
                    norm_err.append(abs(nrm - 1.0))
                else:
                    # Damped wave should generally not amplify norm above initial norm.
                    norm_err.append(max(0.0, nrm - 1.0))
                residual_energy.append(float(np.sum(np.abs(rho_w - rho0))))
            rows.append({
                "split": split,
                "diagnostic": "unitary_norm" if gamma == 0.0 else "damped_norm_decay",
                "tau": float(tau),
                "gamma": float(gamma),
                "mean_norm": float(np.mean(norms)),
                "std_norm": float(np.std(norms)),
                "mean_norm_error_or_amplification": float(np.mean(norm_err)),
                "max_norm_error_or_amplification": float(np.max(norm_err)),
                "mean_wave_state_change_energy": float(np.mean(residual_energy)),
                "graph_key": bundle.get("key"),
                **bundle.get("key_obj", {}),
            })
    return rows


# -----------------------------------------------------------------------------
# Dataset runner
# -----------------------------------------------------------------------------

def run_dataset_v2(dataset_name: str, X_df: pd.DataFrame, cfg: base.BenchConfig, dirs: dict[str, Path],
                   only_cached_graphs: bool, max_pf_configs: Optional[int], target_suite: str,
                   diagnostics: bool) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    X_raw = X_df.to_numpy(dtype=float)
    item_names = list(map(str, X_df.columns))
    pred_rows: list[dict[str, Any]] = []
    align_rows: list[dict[str, Any]] = []
    diag_rows: list[dict[str, Any]] = []

    for horizon in cfg.horizons:
        X0_raw, Y_raw = base.supervised_pairs(X_raw, horizon)
        windows = base.walk_forward_windows(len(X0_raw), cfg.train_frac, cfg.test_frac, cfg.n_windows)
        print(f"[WINDOWS] {dataset_name} horizon={horizon}: {len(windows)} windows")

        for w_id, (train_idx, test_idx) in enumerate(windows):
            X_train_raw = X0_raw[train_idx]
            X_test_raw = X0_raw[test_idx]
            Y_train_raw = Y_raw[train_idx]
            Y_test_raw = Y_raw[test_idx]

            X_train, X_test, mu, sd = base.standardize_train_test(X_train_raw, X_test_raw)
            Y_train = np.nan_to_num((Y_train_raw - mu) / sd, nan=0.0, posinf=0.0, neginf=0.0)
            Y_test = np.nan_to_num((Y_test_raw - mu) / sd, nan=0.0, posinf=0.0, neginf=0.0)
            Xstd_all = np.nan_to_num((X_raw - mu) / sd, nan=0.0, posinf=0.0, neginf=0.0)

            target_train = {
                "next_state": Y_train,
                "residual_vs_persistence": Y_train - X_train,
            }
            target_test = {
                "next_state": Y_test,
                "residual_vs_persistence": Y_test - X_test,
            }
            if target_suite in {"expanded", "full"}:
                target_train.update({
                    "future_abs_change": np.abs(Y_train - X_train),
                    "future_squared_change": (Y_train - X_train) ** 2,
                })
                target_test.update({
                    "future_abs_change": np.abs(Y_test - X_test),
                    "future_squared_change": (Y_test - X_test) ** 2,
                })
                aux_train = make_scalar_targets_for_indices(Xstd_all, train_idx, horizon)
                aux_test = make_scalar_targets_for_indices(Xstd_all, test_idx, horizon)
                target_train.update(aux_train)
                target_test.update(aux_test)

            # Persistence baseline row for next-state only.
            m = base.regression_metrics_matrix(Y_test, X_test)
            pred_rows.append({
                "dataset": dataset_name, "window": w_id, "horizon": horizon,
                "descriptor": "baseline_persistence", "eval_method": "persistence", "target_name": "next_state",
                "k": X_train.shape[1], "feature_dim_in": X_train.shape[1],
                "forecast_r2": m["r2"], "forecast_rmse": m["rmse"], "forecast_mae": m["mae"],
                "forecast_corr": m["mean_channel_corr"], "reconstruction_r2": 1.0,
                "reconstruction_rmse": 0.0, "corr_matrix_error": 0.0,
                "compression_ratio": 1.0, "descriptor_family": "baseline_persistence",
            })

            # Baselines: raw, PCA, random projection, all with alt/nonlinear evaluators and all targets.
            for desc_name, Ftr, Fte, extra in pca_random_baseline_blocks(X_train, X_test, cfg, cfg.seed + w_id):
                rows, aligns = evaluate_block_v2(dataset_name, w_id, horizon, desc_name, Ftr, Fte,
                                                 X_train, X_test, target_train, target_test,
                                                 cfg, cfg.k_values, cfg.seed + 123 * w_id,
                                                 extra=extra, include_alignment=True)
                pred_rows.extend(rows)
                align_rows.extend(aligns)

            pf_count = 0
            for kernel in cfg.kernels:
                for weight_scheme in cfg.weight_schemes:
                    for include_item in cfg.item_identity_options:
                        for use_stability in cfg.stability_options:
                            if max_pf_configs is not None and pf_count >= max_pf_configs:
                                break
                            path = fup.graph_bundle_path_for(X_train, X_test, item_names, dataset_name, w_id, horizon,
                                                             kernel, weight_scheme, include_item, use_stability, dirs)
                            if only_cached_graphs and not path.exists():
                                print(f"[SKIP] missing cached graph: {path.name}")
                                continue
                            try:
                                bundle = base.build_or_load_pf_bundle(
                                    X_train=X_train, X_test=X_test, item_names=item_names,
                                    dataset=dataset_name, window_id=w_id, horizon=horizon,
                                    kernel=kernel, weight_scheme=weight_scheme,
                                    include_item_identity=include_item, use_stability=use_stability,
                                    cfg=cfg, dirs=dirs,
                                )
                                pf_count += 1
                                combo_extra = {
                                    "kernel": kernel,
                                    "weight_scheme": weight_scheme,
                                    "item_identity": bool(include_item),
                                    "stability": bool(use_stability),
                                }
                                # Graph-change target is bundle-specific, so add it here.
                                if target_suite == "full":
                                    target_train_pf = dict(target_train)
                                    target_test_pf = dict(target_test)
                                    target_train_pf["future_graph_change"] = graph_change_target(bundle["prop_graph_train"], horizon)
                                    target_test_pf["future_graph_change"] = graph_change_target(bundle["prop_graph_test"], horizon)
                                else:
                                    target_train_pf = target_train
                                    target_test_pf = target_test

                                blocks = fup.descriptor_blocks_for_bundle(bundle, X_train, X_test, item_names, cfg, dirs)
                                for desc_name, Ftr, Fte, desc_extra in blocks:
                                    rows, aligns = evaluate_block_v2(dataset_name, w_id, horizon, desc_name, Ftr, Fte,
                                                                     X_train, X_test, target_train_pf, target_test_pf,
                                                                     cfg, cfg.k_values,
                                                                     seed=cfg.seed + 10000 * w_id + 97 * pf_count,
                                                                     extra={**combo_extra, **desc_extra},
                                                                     include_alignment=True)
                                    pred_rows.extend(rows)
                                    align_rows.extend(aligns)
                                if diagnostics:
                                    gammas = [0.0] if not cfg.damped_gammas else [0.0, *list(cfg.damped_gammas)]
                                    diag_rows.extend(wave_conservation_diagnostics(bundle, X_train, cfg, "train", gammas))
                                    diag_rows.extend(wave_conservation_diagnostics(bundle, X_test, cfg, "test", gammas))
                            except Exception as e:
                                print("[WARN] v2 config failed:", {
                                    "dataset": dataset_name,
                                    "window": w_id,
                                    "kernel": kernel,
                                    "weight_scheme": weight_scheme,
                                    "item_identity": include_item,
                                    "stability": use_stability,
                                    "error": repr(e),
                                })
                        if max_pf_configs is not None and pf_count >= max_pf_configs:
                            break
                    if max_pf_configs is not None and pf_count >= max_pf_configs:
                        break
                if max_pf_configs is not None and pf_count >= max_pf_configs:
                    break

    return pd.DataFrame(pred_rows), pd.DataFrame(align_rows), pd.DataFrame(diag_rows)


def aggregate_prediction(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    group_cols = ["dataset", "descriptor", "eval_method", "target_name", "k"]
    optional = ["kernel", "weight_scheme", "item_identity", "stability", "descriptor_family", "gamma", "baseline_k"]
    group_cols += [c for c in optional if c in df.columns]
    metrics = [
        "forecast_r2", "forecast_rmse", "forecast_mae", "forecast_corr",
        "reconstruction_r2", "reconstruction_rmse", "corr_matrix_error",
        "rank_spearman", "precision_top_decile", "top_bottom_target_spread", "auc_top_quartile",
    ]
    metrics = [c for c in metrics if c in df.columns]
    out = df.groupby(group_cols, dropna=False)[metrics].agg(["mean", "std", "median", "count"])
    out.columns = ["_".join(c).strip("_") for c in out.columns.to_flat_index()]
    out = out.reset_index()
    if "forecast_r2_mean" in out.columns:
        out = out.sort_values(["dataset", "target_name", "forecast_r2_mean"], ascending=[True, True, False])
    return out


def aggregate_alignment(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    group_cols = ["dataset", "descriptor", "eval_method", "k"]
    optional = ["kernel", "weight_scheme", "item_identity", "stability", "descriptor_family", "gamma", "baseline_k"]
    group_cols += [c for c in optional if c in df.columns]
    metrics = ["alignment_cka", "alignment_gram_fro_error", "alignment_pairwise_spearman"]
    out = df.groupby(group_cols, dropna=False)[metrics].agg(["mean", "std", "median", "count"])
    out.columns = ["_".join(c).strip("_") for c in out.columns.to_flat_index()]
    out = out.reset_index()
    if "alignment_cka_mean" in out.columns:
        out = out.sort_values(["dataset", "alignment_cka_mean"], ascending=[True, False])
    return out


def run_v2(cfg: base.BenchConfig, only_cached_graphs: bool, max_pf_configs: Optional[int], target_suite: str,
           diagnostics: bool) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    dirs = base.ensure_dirs(cfg)
    all_pred = []
    all_align = []
    all_diag = []
    for dataset_name in cfg.datasets:
        print("\n" + "=" * 100)
        print(f"[FOLLOWUP V2 DATASET] {dataset_name}")
        X_df = base.cached_dataset(dataset_name, cfg, dirs)
        print(f"[SHAPE] {dataset_name}: {X_df.shape}")
        pred, align, diag = run_dataset_v2(dataset_name, X_df, cfg, dirs, only_cached_graphs, max_pf_configs, target_suite, diagnostics)
        pred_path = dirs["results"] / f"followup_v2_prediction_{dataset_name}_{cfg.preset}_{target_suite}.csv"
        align_path = dirs["results"] / f"followup_v2_alignment_{dataset_name}_{cfg.preset}_{target_suite}.csv"
        diag_path = dirs["results"] / f"followup_v2_diagnostics_{dataset_name}_{cfg.preset}_{target_suite}.csv"
        pred.to_csv(pred_path, index=False)
        align.to_csv(align_path, index=False)
        diag.to_csv(diag_path, index=False)
        print(f"[SAVE] {pred_path}")
        print(f"[SAVE] {align_path}")
        print(f"[SAVE] {diag_path}")
        all_pred.append(pred)
        all_align.append(align)
        all_diag.append(diag)

    pred_all = pd.concat(all_pred, ignore_index=True) if all_pred else pd.DataFrame()
    align_all = pd.concat(all_align, ignore_index=True) if all_align else pd.DataFrame()
    diag_all = pd.concat(all_diag, ignore_index=True) if all_diag else pd.DataFrame()

    pred_all_path = dirs["results"] / f"followup_v2_prediction_combined_{cfg.preset}_{target_suite}.csv"
    align_all_path = dirs["results"] / f"followup_v2_alignment_combined_{cfg.preset}_{target_suite}.csv"
    diag_all_path = dirs["results"] / f"followup_v2_diagnostics_combined_{cfg.preset}_{target_suite}.csv"
    pred_all.to_csv(pred_all_path, index=False)
    align_all.to_csv(align_all_path, index=False)
    diag_all.to_csv(diag_all_path, index=False)

    pred_agg = aggregate_prediction(pred_all)
    align_agg = aggregate_alignment(align_all)
    pred_agg_path = dirs["results"] / f"followup_v2_prediction_aggregate_{cfg.preset}_{target_suite}.csv"
    align_agg_path = dirs["results"] / f"followup_v2_alignment_aggregate_{cfg.preset}_{target_suite}.csv"
    pred_agg.to_csv(pred_agg_path, index=False)
    align_agg.to_csv(align_agg_path, index=False)

    print("\n[SAVE]", pred_all_path)
    print("[SAVE]", align_all_path)
    print("[SAVE]", diag_all_path)
    print("[SAVE]", pred_agg_path)
    print("[SAVE]", align_agg_path)

    if not pred_agg.empty:
        print("\nTop prediction rows by dataset/target:")
        print(pred_agg.groupby(["dataset", "target_name"], group_keys=False).head(12).to_string(index=False))
    if not align_agg.empty:
        print("\nTop alignment rows by dataset:")
        print(align_agg.groupby("dataset", group_keys=False).head(12).to_string(index=False))
    if not diag_all.empty:
        print("\nDiagnostics preview:")
        print(diag_all.head(20).to_string(index=False))

    return pred_all, align_all, diag_all


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--preset", default="lite", choices=["lite", "markets", "full"])
    p.add_argument("--cache-dir", default=None)
    p.add_argument("--datasets", nargs="*", default=None,
                   choices=["plotly_stock_returns", "uci_household_power", "uci_appliances"])
    p.add_argument("--only-cached-graphs", action="store_true")
    p.add_argument("--max-pf-configs", type=int, default=None)
    p.add_argument("--target-suite", default="expanded", choices=["basic", "expanded", "full"],
                   help="basic: next/residual. expanded: volatility/corr/regime targets. full: expanded + graph-change target.")
    p.add_argument("--no-diagnostics", action="store_true")
    p.add_argument("--no-unitary", action="store_true")
    p.add_argument("--no-diffusion", action="store_true")
    p.add_argument("--no-phase", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    cfg = base.config_for_preset(args.preset)
    if args.cache_dir is not None:
        cfg = base.BenchConfig(**{**asdict(cfg), "cache_dir": args.cache_dir})
    if args.datasets is not None:
        cfg = base.BenchConfig(**{**asdict(cfg), "datasets": tuple(args.datasets)})
    if args.no_unitary:
        cfg = base.BenchConfig(**{**asdict(cfg), "include_unitary": False})
    if args.no_diffusion:
        cfg = base.BenchConfig(**{**asdict(cfg), "include_diffusion": False})
    if args.no_phase:
        cfg = base.BenchConfig(**{**asdict(cfg), "use_phase_from_delta": False})

    print("[FOLLOWUP V2 CONFIG]")
    print(json.dumps(asdict(cfg), indent=2, default=str))
    print("[ONLY CACHED GRAPHS]", args.only_cached_graphs)
    print("[MAX PF CONFIGS]", args.max_pf_configs)
    print("[TARGET SUITE]", args.target_suite)
    print("[DIAGNOSTICS]", not args.no_diagnostics)
    run_v2(cfg, args.only_cached_graphs, args.max_pf_configs, args.target_suite, diagnostics=not args.no_diagnostics)


if __name__ == "__main__":
    main()
