"""
pf_flux_continuum_eval.py
==========================

Cached-graph follow-up focused on flux direction / continuum-mechanics style tests.

It reuses graph bundles already built by pf_better_benchmark.py and descriptor helpers
from pf_followup_eval.py / pf_followup_eval_v2.py. It does NOT rebuild graphs when
--only-cached-graphs is passed.

What it tests
-------------
1. Wave graph-current / fringe flux direction features
   J_ij = A_ij |psi_i| |psi_j| sin(theta_j - theta_i)

2. Diffusion flux / pressure features
   F_ij = A_ij (rho_j - rho_i)

3. Raw continuum / graph-gradient features from the current standardized state
   G_ij = A_ij (x_j - x_i)

4. Conditional flux features:
   momentum alignment, acceleration alignment, pressure/concentration, strain/vorticity proxies.

5. Targets better suited to markets than raw next returns:
   - future_residual_return
   - future_rank_return
   - future_acceleration
   - future_abs_acceleration
   - future_squared_acceleration
   - future_acceleration_norm
   - future_acceleration_dispersion
   - future_dispersion_change
   - future_abs_return_mean
   - future_squared_return_mean
   - next_state / residual_vs_persistence for comparison

Run examples
------------
python pf_flux_continuum_eval.py --preset lite --datasets plotly_stock_returns --only-cached-graphs --max-pf-configs 2
python pf_flux_continuum_eval.py --preset lite --datasets plotly_stock_returns --only-cached-graphs --target-suite flux
python pf_flux_continuum_eval.py --preset lite --only-cached-graphs --target-suite flux

Notes
-----
This script is intentionally conservative: all supervised feature selection uses train-only data;
no future values are used as input features. Alignment-with-future is NOT used as a feature.
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
import pf_followup_eval_v2 as v2


EPS = 1e-12


# -----------------------------------------------------------------------------
# Graph helpers
# -----------------------------------------------------------------------------

def dense_adjacency_at(graph: Any, s: int) -> np.ndarray:
    if hasattr(graph, "adjacency"):
        A = np.asarray(graph.adjacency[s], dtype=float)
        A = 0.5 * (A + A.T)
        np.fill_diagonal(A, 0.0)
        return np.nan_to_num(A, nan=0.0, posinf=0.0, neginf=0.0)
    if hasattr(graph, "edges"):
        n = int(graph.n_items)
        A = np.zeros((n, n), dtype=float)
        part = graph.edges[graph.edges["slice_idx"] == s]
        for row in part.itertuples(index=False):
            i = int(row.source_idx); j = int(row.target_idx); w = float(row.weight)
            A[i, j] += w; A[j, i] += w
        return A
    raise TypeError("graph bundle graph must expose adjacency or edges")


def wrap_phase_diff(dtheta: np.ndarray) -> np.ndarray:
    return (dtheta + np.pi) % (2.0 * np.pi) - np.pi


def safe_skew(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float).ravel()
    if x.size == 0:
        return 0.0
    m = float(np.mean(x)); s = float(np.std(x))
    if s <= EPS:
        return 0.0
    return float(np.mean(((x - m) / s) ** 3))


def safe_kurtosis(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float).ravel()
    if x.size == 0:
        return 0.0
    m = float(np.mean(x)); s = float(np.std(x))
    if s <= EPS:
        return 0.0
    return float(np.mean(((x - m) / s) ** 4) - 3.0)


def summarize_node_field(x: np.ndarray, prefix: str) -> list[float]:
    x = np.nan_to_num(np.asarray(x, dtype=float).ravel(), nan=0.0, posinf=0.0, neginf=0.0)
    if x.size == 0:
        return [0.0] * 10
    ax = np.abs(x)
    return [
        float(np.mean(x)),
        float(np.std(x)),
        float(np.min(x)),
        float(np.max(x)),
        float(np.mean(ax)),
        float(np.max(ax)),
        float(np.linalg.norm(x) / np.sqrt(max(1, x.size))),
        float(np.quantile(x, 0.10)),
        float(np.quantile(x, 0.90)),
        safe_skew(x),
    ]


def node_summary_names(prefix: str) -> list[str]:
    return [
        f"{prefix}_mean", f"{prefix}_std", f"{prefix}_min", f"{prefix}_max",
        f"{prefix}_abs_mean", f"{prefix}_abs_max", f"{prefix}_rms",
        f"{prefix}_q10", f"{prefix}_q90", f"{prefix}_skew",
    ]


def graph_current_summaries(A: np.ndarray, J: np.ndarray, phase: Optional[np.ndarray] = None, amp: Optional[np.ndarray] = None) -> tuple[list[float], list[str]]:
    deg = A.sum(axis=1)
    absJ = np.abs(J)
    pos = np.maximum(J, 0.0)
    neg = np.maximum(-J, 0.0)
    net = J.sum(axis=1)
    incoming = pos.sum(axis=0) + neg.sum(axis=1)  # positive flow into node under orientation convention
    outgoing = pos.sum(axis=1) + neg.sum(axis=0)
    abs_flow = absJ.sum(axis=1)
    pressure = abs_flow / (deg + EPS)
    imbalance = incoming - outgoing
    concentration = float((np.max(abs_flow) + EPS) / (np.sum(abs_flow) + EPS))
    gini_like = float(np.mean(np.abs(abs_flow[:, None] - abs_flow[None, :])) / (2.0 * (np.mean(abs_flow) + EPS))) if abs_flow.size else 0.0
    names: list[str] = []
    vals: list[float] = []
    for name, field in [
        ("net_flux", net),
        ("abs_flux", abs_flow),
        ("incoming_flux", incoming),
        ("outgoing_flux", outgoing),
        ("flux_imbalance", imbalance),
        ("pressure", pressure),
    ]:
        vals.extend(summarize_node_field(field, name)); names.extend(node_summary_names(name))
    vals.extend([
        float(np.sum(pos)),
        float(np.sum(neg)),
        float(np.sum(absJ)),
        float(np.linalg.norm(J, ord="fro")),
        float(concentration),
        float(gini_like),
        float(np.trace(J @ J.T) / (np.sum(A > 0) + EPS)),
    ])
    names.extend([
        "positive_current_mass", "negative_current_mass", "absolute_current_mass",
        "current_fro_energy", "current_concentration", "current_gini_like",
        "edge_current_energy_per_edge",
    ])
    if phase is not None:
        dtheta = wrap_phase_diff(phase[None, :] - phase[:, None])
        grad = (A * np.abs(dtheta)).sum(axis=1) / (deg + EPS)
        vals.extend(summarize_node_field(grad, "phase_gradient")); names.extend(node_summary_names("phase_gradient"))
        vals.append(float(np.abs(np.mean(np.exp(1j * phase)))))
        names.append("phase_coherence")
        if amp is not None:
            w = np.asarray(amp, dtype=float)
            vals.append(float(np.abs(np.sum(w * np.exp(1j * phase))) / (np.sum(np.abs(w)) + EPS)))
            names.append("amp_weighted_phase_coherence")
    return vals, names


# -----------------------------------------------------------------------------
# Flux / continuum descriptors
# -----------------------------------------------------------------------------

def wave_flux_descriptor_for_graph(graph: Any, X: np.ndarray, cfg: base.BenchConfig, gamma: float) -> tuple[np.ndarray, list[str]]:
    wave = base.GraphWaveDynamics(generator=cfg.graph_generator, gamma=float(gamma))
    rows: list[list[float]] = []
    names: Optional[list[str]] = None
    for s in range(X.shape[0]):
        A = dense_adjacency_at(graph, s)
        rho0, psi0 = base.density_and_phase_from_state(X, s, cfg.use_phase_from_delta)
        op = wave.prepare(graph, s)
        # Use the first tau value for flux descriptors to keep feature budget stable.
        psi = wave.evolve(op, psi0, float(cfg.tau_values[0]))
        amp = np.abs(psi)
        phase = np.angle(psi)
        dtheta = wrap_phase_diff(phase[None, :] - phase[:, None])
        J = A * (amp[:, None] * amp[None, :]) * np.sin(dtheta)
        vals, nm = graph_current_summaries(A, J, phase=phase, amp=amp)
        # Momentum/flow conditional features use only current/past state.
        x = np.asarray(X[s], dtype=float)
        if s > 0:
            mom = np.asarray(X[s] - X[s - 1], dtype=float)
        else:
            mom = np.zeros_like(x)
        net = J.sum(axis=1)
        abs_flow = np.abs(J).sum(axis=1)
        vals.extend([
            float(np.dot(net, x) / ((np.linalg.norm(net) * np.linalg.norm(x)) + EPS)),
            float(np.dot(net, mom) / ((np.linalg.norm(net) * np.linalg.norm(mom)) + EPS)),
            float(np.dot(abs_flow, np.abs(mom)) / ((np.linalg.norm(abs_flow) * np.linalg.norm(np.abs(mom))) + EPS)),
            float(np.std(net)),
            float(np.std(abs_flow)),
        ])
        nm.extend(["flux_state_cosine", "flux_momentum_cosine", "pressure_abs_momentum_cosine", "net_flux_dispersion", "pressure_dispersion"])
        rows.append(vals)
        names = nm
    return np.asarray(rows, dtype=float), list(names or [])


def diffusion_flux_descriptor_for_graph(graph: Any, X: np.ndarray, cfg: base.BenchConfig) -> tuple[np.ndarray, list[str]]:
    diff = base.GraphDiffusionDynamics(generator=cfg.graph_generator, normalize_output=True)
    rows: list[list[float]] = []
    names: Optional[list[str]] = None
    for s in range(X.shape[0]):
        A = dense_adjacency_at(graph, s)
        rho0, _ = base.density_and_phase_from_state(X, s, cfg.use_phase_from_delta)
        op = diff.prepare(graph, s)
        rho = diff.evolve(op, rho0, float(cfg.tau_values[0]))
        # Flow from lower density to higher density under simple graph-gradient convention.
        F = A * (rho[None, :] - rho[:, None])
        vals, nm = graph_current_summaries(A, F)
        x = np.asarray(X[s], dtype=float)
        net = F.sum(axis=1)
        vals.append(float(np.dot(net, x) / ((np.linalg.norm(net) * np.linalg.norm(x)) + EPS)))
        nm.append("diffusion_flux_state_cosine")
        rows.append(vals)
        names = nm
    return np.asarray(rows, dtype=float), list(names or [])


def raw_continuum_descriptor_for_graph(graph: Any, X: np.ndarray) -> tuple[np.ndarray, list[str]]:
    rows: list[list[float]] = []
    names: Optional[list[str]] = None
    for s in range(X.shape[0]):
        A = dense_adjacency_at(graph, s)
        x = np.asarray(X[s], dtype=float)
        if s > 0:
            mom = np.asarray(X[s] - X[s - 1], dtype=float)
        else:
            mom = np.zeros_like(x)
        # Graph-gradient current on observed field and on momentum field.
        Gx = A * (x[None, :] - x[:, None])
        Gm = A * (mom[None, :] - mom[:, None])
        vals_x, nm_x = graph_current_summaries(A, Gx)
        vals_m, nm_m = graph_current_summaries(A, Gm)
        vals = vals_x + vals_m
        nm = [f"stategrad_{n}" for n in nm_x] + [f"momgrad_{n}" for n in nm_m]
        # Continuum-ish scalar summaries.
        deg = A.sum(axis=1)
        lap_x = (A * (x[:, None] - x[None, :])).sum(axis=1)
        lap_m = (A * (mom[:, None] - mom[None, :])).sum(axis=1)
        local_strain = np.abs(Gx).sum(axis=1) / (deg + EPS)
        local_mom_strain = np.abs(Gm).sum(axis=1) / (deg + EPS)
        for name, field in [("lap_state", lap_x), ("lap_momentum", lap_m), ("local_strain", local_strain), ("local_mom_strain", local_mom_strain)]:
            vals.extend(summarize_node_field(field, name)); nm.extend(node_summary_names(name))
        vals.extend([
            float(np.dot(lap_x, mom) / ((np.linalg.norm(lap_x) * np.linalg.norm(mom)) + EPS)),
            float(np.dot(local_strain, np.abs(mom)) / ((np.linalg.norm(local_strain) * np.linalg.norm(np.abs(mom))) + EPS)),
        ])
        nm.extend(["lap_state_momentum_cosine", "strain_abs_momentum_cosine"])
        rows.append(vals)
        names = nm
    return np.asarray(rows, dtype=float), list(names or [])


def build_flux_blocks(bundle: dict[str, Any], X_train: np.ndarray, X_test: np.ndarray, cfg: base.BenchConfig, dirs: dict[str, Path]) -> list[tuple[str, np.ndarray, np.ndarray, dict[str, Any]]]:
    blocks: list[tuple[str, np.ndarray, np.ndarray, dict[str, Any]]] = []
    key = str(bundle.get("key", "nokey"))
    desc_dir = dirs["descriptors"] / "flux_continuum"
    desc_dir.mkdir(parents=True, exist_ok=True)

    def cached(name: str, builder):
        path = desc_dir / f"{name}_{key}.npz"
        if path.exists():
            data = np.load(path, allow_pickle=True)
            return data["train"], data["test"]
        Ftr, names_tr = builder("train")
        Fte, names_te = builder("test")
        np.savez_compressed(path, train=Ftr, test=Fte, names=np.asarray(names_tr, dtype=object))
        return Ftr, Fte

    Ftr, Fte = cached("raw_continuum_v1", lambda split: raw_continuum_descriptor_for_graph(bundle["prop_graph_train"] if split == "train" else bundle["prop_graph_test"], X_train if split == "train" else X_test))
    blocks.append(("pf_raw_continuum_flux", Ftr, Fte, {"descriptor_family": "raw_continuum_flux"}))

    Ftr, Fte = cached("diffusion_flux_v1", lambda split: diffusion_flux_descriptor_for_graph(bundle["prop_graph_train"] if split == "train" else bundle["prop_graph_test"], X_train if split == "train" else X_test, cfg))
    blocks.append(("pf_diffusion_flux", Ftr, Fte, {"descriptor_family": "diffusion_flux", "gamma": 0.0}))

    # Keep prior arithmetic descriptor available to test arithmetic+flux mixtures.
    Atr, Ate = fup.build_arithmetic_direct_descriptor(bundle)
    blocks.append(("pf_arithmetic_direct", Atr, Ate, {"descriptor_family": "arithmetic_direct"}))
    blocks.append(("pf_arithmetic_plus_raw_continuum_flux", np.hstack([Atr, Ftr]), np.hstack([Ate, Fte]), {"descriptor_family": "arithmetic+raw_continuum_flux"}))

    gammas = [0.0, *list(cfg.damped_gammas)]
    for gamma in gammas:
        label = "unitary" if float(gamma) == 0.0 else f"damped_g{gamma}"
        Wtr, Wte = cached(f"wave_flux_{label}_v1", lambda split, gamma=gamma: wave_flux_descriptor_for_graph(bundle["prop_graph_train"] if split == "train" else bundle["prop_graph_test"], X_train if split == "train" else X_test, cfg, float(gamma)))
        dname = "pf_unitary_wave_flux" if float(gamma) == 0.0 else f"pf_damped_wave_flux_g{gamma}"
        fam = "unitary_wave_flux" if float(gamma) == 0.0 else "damped_wave_flux"
        blocks.append((dname, Wtr, Wte, {"descriptor_family": fam, "gamma": float(gamma)}))
        blocks.append((f"pf_arithmetic_plus_{dname[3:]}", np.hstack([Atr, Wtr]), np.hstack([Ate, Wte]), {"descriptor_family": "arithmetic+" + fam, "gamma": float(gamma)}))
        blocks.append((f"pf_continuum_plus_{dname[3:]}", np.hstack([Ftr, Wtr]), np.hstack([Fte, Wte]), {"descriptor_family": "continuum+" + fam, "gamma": float(gamma)}))
    return blocks


# -----------------------------------------------------------------------------
# Target helpers
# -----------------------------------------------------------------------------

def rank_rows(M: np.ndarray) -> np.ndarray:
    M = np.asarray(M, dtype=float)
    out = np.zeros_like(M, dtype=float)
    for i in range(M.shape[0]):
        r = v2.rankdata_average(M[i])
        if M.shape[1] > 1:
            r = (r - np.mean(r)) / (np.std(r) + EPS)
        out[i] = r
    return out


def scalar_col(x: np.ndarray) -> np.ndarray:
    return np.asarray(x, dtype=float).reshape(-1, 1)


def make_flux_target_maps(X_train: np.ndarray, X_test: np.ndarray, Y_train: np.ndarray, Y_test: np.ndarray) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    dtr = Y_train - X_train
    dte = Y_test - X_test
    # Cross-sectional residual return removes contemporaneous market mode from the future vector.
    res_tr = Y_train - Y_train.mean(axis=1, keepdims=True)
    res_te = Y_test - Y_test.mean(axis=1, keepdims=True)
    target_train: dict[str, np.ndarray] = {
        "next_state": Y_train,
        "residual_vs_persistence": dtr,
        "future_residual_return": res_tr,
        "future_rank_return": rank_rows(Y_train),
        "future_acceleration": dtr,
        "future_abs_acceleration": np.abs(dtr),
        "future_squared_acceleration": dtr ** 2,
        "future_acceleration_norm": scalar_col(np.linalg.norm(dtr, axis=1) / np.sqrt(max(1, dtr.shape[1]))),
        "future_acceleration_dispersion": scalar_col(np.std(dtr, axis=1)),
        "future_dispersion_change": scalar_col(np.std(Y_train, axis=1) - np.std(X_train, axis=1)),
        "future_abs_return_mean": scalar_col(np.mean(np.abs(Y_train), axis=1)),
        "future_squared_return_mean": scalar_col(np.mean(Y_train ** 2, axis=1)),
    }
    target_test: dict[str, np.ndarray] = {
        "next_state": Y_test,
        "residual_vs_persistence": dte,
        "future_residual_return": res_te,
        "future_rank_return": rank_rows(Y_test),
        "future_acceleration": dte,
        "future_abs_acceleration": np.abs(dte),
        "future_squared_acceleration": dte ** 2,
        "future_acceleration_norm": scalar_col(np.linalg.norm(dte, axis=1) / np.sqrt(max(1, dte.shape[1]))),
        "future_acceleration_dispersion": scalar_col(np.std(dte, axis=1)),
        "future_dispersion_change": scalar_col(np.std(Y_test, axis=1) - np.std(X_test, axis=1)),
        "future_abs_return_mean": scalar_col(np.mean(np.abs(Y_test), axis=1)),
        "future_squared_return_mean": scalar_col(np.mean(Y_test ** 2, axis=1)),
    }
    return target_train, target_test


# -----------------------------------------------------------------------------
# Runner
# -----------------------------------------------------------------------------

def aggregate_prediction(df: pd.DataFrame) -> pd.DataFrame:
    return v2.aggregate_prediction(df)


def run_flux_dataset(dataset_name: str, X_df: pd.DataFrame, cfg: base.BenchConfig, dirs: dict[str, Path], only_cached_graphs: bool, max_pf_configs: Optional[int]) -> pd.DataFrame:
    X_raw = X_df.to_numpy(dtype=float)
    item_names = list(map(str, X_df.columns))
    rows: list[dict[str, Any]] = []
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
            target_train, target_test = make_flux_target_maps(X_train, X_test, Y_train, Y_test)

            # Baselines with the same evaluator stack.
            for desc_name, Ftr, Fte, extra in v2.pca_random_baseline_blocks(X_train, X_test, cfg, cfg.seed + w_id):
                r, _ = v2.evaluate_block_v2(dataset_name, w_id, horizon, desc_name, Ftr, Fte, X_train, X_test, target_train, target_test, cfg, cfg.k_values, cfg.seed + 117 * w_id, extra=extra, include_alignment=False)
                rows.extend(r)

            pf_count = 0
            for kernel in cfg.kernels:
                for weight_scheme in cfg.weight_schemes:
                    for include_item in cfg.item_identity_options:
                        for use_stability in cfg.stability_options:
                            if max_pf_configs is not None and pf_count >= max_pf_configs:
                                break
                            path = fup.graph_bundle_path_for(X_train, X_test, item_names, dataset_name, w_id, horizon, kernel, weight_scheme, include_item, use_stability, dirs)
                            if only_cached_graphs and not path.exists():
                                print(f"[SKIP] missing cached graph: {path.name}")
                                continue
                            try:
                                bundle = base.build_or_load_pf_bundle(X_train=X_train, X_test=X_test, item_names=item_names, dataset=dataset_name, window_id=w_id, horizon=horizon, kernel=kernel, weight_scheme=weight_scheme, include_item_identity=include_item, use_stability=use_stability, cfg=cfg, dirs=dirs)
                                pf_count += 1
                                combo_extra = {"kernel": kernel, "weight_scheme": weight_scheme, "item_identity": bool(include_item), "stability": bool(use_stability)}
                                blocks = build_flux_blocks(bundle, X_train, X_test, cfg, dirs)
                                for desc_name, Ftr, Fte, extra in blocks:
                                    r, _ = v2.evaluate_block_v2(dataset_name, w_id, horizon, desc_name, Ftr, Fte, X_train, X_test, target_train, target_test, cfg, cfg.k_values, cfg.seed + 10000 * w_id + 97 * pf_count, extra={**combo_extra, **extra}, include_alignment=False)
                                    rows.extend(r)
                            except Exception as e:
                                print("[WARN] flux config failed:", {"dataset": dataset_name, "window": w_id, "kernel": kernel, "weight_scheme": weight_scheme, "item_identity": include_item, "stability": use_stability, "error": repr(e)})
                        if max_pf_configs is not None and pf_count >= max_pf_configs:
                            break
                    if max_pf_configs is not None and pf_count >= max_pf_configs:
                        break
                if max_pf_configs is not None and pf_count >= max_pf_configs:
                    break
    return pd.DataFrame(rows)


def run_flux(cfg: base.BenchConfig, only_cached_graphs: bool, max_pf_configs: Optional[int]) -> pd.DataFrame:
    dirs = base.ensure_dirs(cfg)
    all_rows: list[pd.DataFrame] = []
    for dataset_name in cfg.datasets:
        print("\n" + "=" * 100)
        print(f"[FLUX/CONTINUUM DATASET] {dataset_name}")
        X_df = base.cached_dataset(dataset_name, cfg, dirs)
        print(f"[SHAPE] {dataset_name}: {X_df.shape}")
        df = run_flux_dataset(dataset_name, X_df, cfg, dirs, only_cached_graphs, max_pf_configs)
        out_path = dirs["results"] / f"flux_continuum_prediction_{dataset_name}_{cfg.preset}.csv"
        df.to_csv(out_path, index=False)
        print(f"[SAVE] {out_path}")
        all_rows.append(df)
    combined = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    combined_path = dirs["results"] / f"flux_continuum_prediction_combined_{cfg.preset}.csv"
    combined.to_csv(combined_path, index=False)
    agg = aggregate_prediction(combined)
    agg_path = dirs["results"] / f"flux_continuum_prediction_aggregate_{cfg.preset}.csv"
    agg.to_csv(agg_path, index=False)
    print(f"[SAVE] {combined_path}")
    print(f"[SAVE] {agg_path}")
    if not agg.empty:
        focus = [
            "future_residual_return", "future_rank_return", "future_acceleration",
            "future_acceleration_norm", "future_acceleration_dispersion",
            "future_dispersion_change", "future_abs_return_mean", "future_squared_return_mean",
        ]
        print("\nTop flux/continuum rows by dataset/target:")
        show = agg[agg["target_name"].isin(focus)] if "target_name" in agg.columns else agg
        print(show.groupby(["dataset", "target_name"], group_keys=False).head(10).to_string(index=False))
    return combined


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--preset", default="lite", choices=["lite", "markets", "full"])
    p.add_argument("--cache-dir", default=None)
    p.add_argument("--datasets", nargs="*", default=None, choices=["plotly_stock_returns", "uci_household_power", "uci_appliances"])
    p.add_argument("--only-cached-graphs", action="store_true")
    p.add_argument("--max-pf-configs", type=int, default=None)
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
    print("[FLUX/CONTINUUM CONFIG]")
    print(json.dumps(asdict(cfg), indent=2, default=str))
    print("[ONLY CACHED GRAPHS]", args.only_cached_graphs)
    print("[MAX PF CONFIGS]", args.max_pf_configs)
    run_flux(cfg, args.only_cached_graphs, args.max_pf_configs)


if __name__ == "__main__":
    main()
