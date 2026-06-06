"""
pf_followup_eval.py
====================

Follow-up evaluator for pf_better_benchmark.py.

Purpose
-------
Reuse cached datasets / graph bundles from pf_better_benchmark and run deeper
post-hoc tests without refetching data or reconstructing graphs:

  1. Direct arithmetic feature descriptors from cached valuation tensors.
  2. Fringe / phase-current / orientation descriptors for unitary and damped waves.
  3. Alternative evaluation methods for both old and new descriptors:
       - direct ridge, no PCA
       - SVD/PCA + ridge, original-style redo
       - supervised top-k correlation feature selection + ridge
       - random tanh/ELM nonlinear ridge
       - random Fourier/RBF nonlinear ridge
       - residual prediction over persistence baseline

This script is deliberately pure NumPy/Pandas: no sklearn dependency.

Expected layout
---------------
    project/
        prime_fock/
            ...
        pf_better_benchmark.py
        pf_followup_eval.py
        pf_bench_cache/
            datasets/
            graphs/
            descriptors/
            results/

Typical commands
----------------
    python pf_followup_eval.py --preset lite --datasets plotly_stock_returns
    python pf_followup_eval.py --preset markets
    python pf_followup_eval.py --preset lite --only-cached-graphs

Notes
-----
If graph bundles do not exist, this script can build them by calling the same
builder in pf_better_benchmark. Use --only-cached-graphs to avoid that and only
evaluate already-trained/cached graph bundles.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any, Literal, Optional, Sequence

import numpy as np
import pandas as pd

# Import the original benchmark as a module. Keep this file next to it.
import pf_better_benchmark as base


DescriptorName = Literal[
    "static",
    "arithmetic_direct",
    "diffusion",
    "unitary_wave",
    "damped_wave",
    "unitary_orientation",
    "damped_orientation",
    "static_plus_arithmetic",
    "arithmetic_plus_diffusion",
    "arithmetic_plus_unitary_orientation",
    "arithmetic_plus_damped_orientation",
    "static_plus_arithmetic_plus_damped_orientation",
]


# -----------------------------------------------------------------------------
# Small numerical helpers
# -----------------------------------------------------------------------------

def standardize_features(F_train: np.ndarray, F_test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    F_train = np.asarray(F_train, dtype=float)
    F_test = np.asarray(F_test, dtype=float)
    mu = np.nanmean(F_train, axis=0)
    sd = np.nanstd(F_train, axis=0)
    sd = np.where(sd <= 1e-12, 1.0, sd)
    Ztr = np.nan_to_num((F_train - mu) / sd, nan=0.0, posinf=0.0, neginf=0.0)
    Zte = np.nan_to_num((F_test - mu) / sd, nan=0.0, posinf=0.0, neginf=0.0)
    return Ztr, Zte


def safe_corr_scores(F_train: np.ndarray, Y_train: np.ndarray) -> np.ndarray:
    """
    Multi-output supervised feature score using train-only absolute correlation.
    score_j = mean_c |corr(F_j, Y_c)|
    """
    F, _ = standardize_features(F_train, F_train)
    Y, _ = standardize_features(Y_train, Y_train)
    n = max(1, F.shape[0] - 1)
    C = (F.T @ Y) / float(n)
    C = np.nan_to_num(C, nan=0.0, posinf=0.0, neginf=0.0)
    return np.mean(np.abs(C), axis=1)


def select_topcorr(F_train: np.ndarray, F_test: np.ndarray, Y_train: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    scores = safe_corr_scores(F_train, Y_train)
    k_eff = min(max(1, int(k)), F_train.shape[1])
    idx = np.argsort(scores)[-k_eff:]
    idx = idx[np.argsort(scores[idx])[::-1]]
    return F_train[:, idx], F_test[:, idx], idx


def random_tanh_features(F_train: np.ndarray, F_test: np.ndarray, width: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    Xtr, Xte = standardize_features(F_train, F_test)
    rng = np.random.default_rng(seed)
    d = max(1, Xtr.shape[1])
    W = rng.normal(0.0, 1.0 / np.sqrt(d), size=(d, width))
    b = rng.uniform(-np.pi, np.pi, size=width)
    return np.tanh(Xtr @ W + b), np.tanh(Xte @ W + b)


def rff_rbf_features(F_train: np.ndarray, F_test: np.ndarray, width: int, seed: int, gamma: Optional[float] = None) -> tuple[np.ndarray, np.ndarray]:
    Xtr, Xte = standardize_features(F_train, F_test)
    rng = np.random.default_rng(seed)
    d = max(1, Xtr.shape[1])
    if gamma is None:
        # For standardized inputs, this is a conservative default.
        gamma = 1.0 / float(d)
    W = rng.normal(0.0, np.sqrt(2.0 * gamma), size=(d, width))
    b = rng.uniform(0.0, 2.0 * np.pi, size=width)
    scale = np.sqrt(2.0 / float(width))
    return scale * np.cos(Xtr @ W + b), scale * np.cos(Xte @ W + b)


def eval_Z(
    dataset: str,
    window_id: int,
    horizon: int,
    descriptor: str,
    eval_method: str,
    k: int,
    Z_train: np.ndarray,
    Z_test: np.ndarray,
    X_train: np.ndarray,
    X_test: np.ndarray,
    Y_train: np.ndarray,
    Y_test: np.ndarray,
    cfg: base.BenchConfig,
    target_mode: Literal["next_state", "residual_vs_persistence"] = "next_state",
    extra: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    """Evaluate a representation with either normal next-state target or residual-over-persistence target."""
    if target_mode == "next_state":
        Yfit = Y_train
        Y_pred = base.ridge_fit_predict(Z_train, Yfit, Z_test, alpha=cfg.ridge_alpha)
    elif target_mode == "residual_vs_persistence":
        # Train target is the part persistence does not explain.
        Yfit = Y_train - X_train
        residual_pred = base.ridge_fit_predict(Z_train, Yfit, Z_test, alpha=cfg.ridge_alpha)
        Y_pred = X_test + residual_pred
    else:
        raise ValueError(target_mode)

    # Reconstruction is always same-time X from descriptor, useful for compression diagnostics.
    X_rec = base.ridge_fit_predict(Z_train, X_train, Z_test, alpha=cfg.ridge_alpha)
    mf = base.regression_metrics_matrix(Y_test, Y_pred)
    mr = base.regression_metrics_matrix(X_test, X_rec)

    row = {
        "dataset": dataset,
        "window": int(window_id),
        "horizon": int(horizon),
        "descriptor": descriptor,
        "eval_method": eval_method,
        "target_mode": target_mode,
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
    }
    if extra:
        row.update(extra)
    return row


# -----------------------------------------------------------------------------
# Direct arithmetic descriptors
# -----------------------------------------------------------------------------

def _resolve_weights(bundle: dict[str, Any], encoded_key: str) -> np.ndarray:
    encoded = bundle[encoded_key]
    if bundle.get("mode_weights") is None:
        return np.log(np.asarray(encoded.primes, dtype=float))
    return np.asarray(bundle["mode_weights"], dtype=float)


def _weighted_mass(E: np.ndarray, w: np.ndarray) -> np.ndarray:
    return E @ w


def _pairwise_arithmetic_mats(E: np.ndarray, w: np.ndarray, eps: float = 1e-12) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return gcd, lcm, distance, jaccard, coprime matrices for one slice E[I,P]."""
    # I is small in our benchmarks, so dense pairwise is fine here.
    mn = np.minimum(E[:, None, :], E[None, :, :])
    mx = np.maximum(E[:, None, :], E[None, :, :])
    gcd = np.tensordot(mn, w, axes=([2], [0]))
    lcm = np.tensordot(mx, w, axes=([2], [0]))
    dist = np.maximum(lcm - gcd, 0.0)
    jac = gcd / (lcm + eps)
    cop = (gcd <= eps).astype(float)
    for M in (gcd, lcm, dist, jac, cop):
        np.fill_diagonal(M, 0.0)
    return gcd, lcm, dist, jac, cop


def _row_summaries(M: np.ndarray, prefix: str) -> np.ndarray:
    I = M.shape[0]
    if I <= 1:
        vals = np.zeros((I, 7), dtype=float)
        return vals
    denom = max(1, I - 1)
    row_sum = M.sum(axis=1)
    row_mean = row_sum / denom
    row_max = M.max(axis=1)
    row_std = M.std(axis=1)
    # top-k pooled edge mass per node.
    k_top = min(3, max(1, I - 1))
    top = np.sort(M, axis=1)[:, -k_top:]
    top_mean = top.mean(axis=1)
    top_sum = top.sum(axis=1)
    concentration = row_max / (row_sum + 1e-12)
    return np.column_stack([row_sum, row_mean, row_max, row_std, top_mean, top_sum, concentration])


def _global_summaries(*arrays: np.ndarray) -> list[float]:
    out: list[float] = []
    for a in arrays:
        x = np.asarray(a, dtype=float).ravel()
        if x.size == 0:
            out.extend([0.0, 0.0, 0.0, 0.0])
        else:
            out.extend([float(np.mean(x)), float(np.std(x)), float(np.max(x)), float(np.sum(x))])
    return out


def arithmetic_direct_matrix(bundle: dict[str, Any], encoded_key: str, graph_key: str, stability_key: str) -> np.ndarray:
    """
    State-level direct arithmetic descriptor from cached valuation tensor.

    Uses explicit gcd/lcm/distance/jaccard/coprime summaries plus mass and
    stable-vs-raw graph degree summaries. This tests arithmetic as features,
    not merely as a graph kernel.
    """
    encoded = bundle[encoded_key]
    Eall = np.asarray(encoded.exponents, dtype=float)  # S,I,P
    w = _resolve_weights(bundle, encoded_key)
    graph = bundle[graph_key]
    stability = bundle[stability_key]
    A_raw = np.asarray(graph.adjacency, dtype=float)
    A_stable = np.asarray(stability.stable_adjacency, dtype=float)

    rows = []
    for s in range(Eall.shape[0]):
        E = Eall[s]
        mass = _weighted_mass(E, w)
        support = (E > 0).sum(axis=1).astype(float)
        exponent_sum = E.sum(axis=1)

        gcd, lcm, dist, jac, cop = _pairwise_arithmetic_mats(E, w)
        node_parts = [
            mass[:, None],
            support[:, None],
            exponent_sum[:, None],
            _row_summaries(gcd, "gcd"),
            _row_summaries(lcm, "lcm"),
            _row_summaries(dist, "dist"),
            _row_summaries(jac, "jac"),
            _row_summaries(cop, "cop"),
            _row_summaries(A_raw[s], "rawA"),
            _row_summaries(A_stable[s], "stableA"),
        ]
        node_feat = np.hstack(node_parts).ravel()

        stable_degree = A_stable[s].sum(axis=1)
        raw_degree = A_raw[s].sum(axis=1)
        stable_fraction = stable_degree / (raw_degree + 1e-12)

        global_feat = _global_summaries(
            mass,
            support,
            gcd,
            lcm,
            dist,
            jac,
            cop,
            raw_degree,
            stable_degree,
            stable_fraction,
        )
        rows.append(np.concatenate([node_feat, np.asarray(global_feat, dtype=float)]))

    return np.nan_to_num(np.asarray(rows, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)


def build_arithmetic_direct_descriptor(bundle: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    return (
        arithmetic_direct_matrix(bundle, "encoded_train", "graph_train", "stability_train"),
        arithmetic_direct_matrix(bundle, "encoded_test", "graph_test", "stability_test"),
    )


# -----------------------------------------------------------------------------
# Wave fringe / orientation descriptors
# -----------------------------------------------------------------------------

def _adjacency_for_slice(graph: Any, s: int) -> np.ndarray:
    if hasattr(graph, "adjacency"):
        A = np.asarray(graph.adjacency, dtype=float)[s]
        return 0.5 * (A + A.T)
    raise TypeError("This follow-up evaluator expects dense SharedSupportGraph cache objects.")


def _orientation_features(A: np.ndarray, psi: np.ndarray) -> np.ndarray:
    """
    Fringe orientation / graph current features from evolved complex wave.

    current_ij = W_ij * amp_i * amp_j * sin(theta_j - theta_i)

    This captures signed phase-flow direction over the graph. It is distinct
    from interference magnitude: it asks where phase/fringe structure is trying
    to push mass, not just how strong cancellation/reinforcement is.
    """
    A = np.asarray(A, dtype=float)
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    amp = np.abs(psi)
    phase = np.angle(psi)
    dtheta = np.angle(np.exp(1j * (phase[None, :] - phase[:, None])))

    amp_pair = amp[:, None] * amp[None, :]
    current = A * amp_pair * np.sin(dtheta)
    abs_current = np.abs(current)

    # Weighted wrapped phase gradient and neighbor phase alignment.
    deg = A.sum(axis=1) + 1e-12
    phase_grad_abs = (A * np.abs(dtheta)).sum(axis=1) / deg
    phase_alignment = (A * np.cos(dtheta)).sum(axis=1) / deg

    node_signed_flow = current.sum(axis=1)
    node_abs_flow = abs_current.sum(axis=1)
    node_out_pos = np.maximum(current, 0.0).sum(axis=1)
    node_in_neg = np.maximum(-current, 0.0).sum(axis=1)
    node_anisotropy = np.abs(node_signed_flow) / (node_abs_flow + 1e-12)

    # Circular/Kuramoto coherence. Weighted by amplitude as well as unweighted.
    coherence = np.abs(np.mean(np.exp(1j * phase)))
    amp_coherence = np.abs(np.sum(amp * np.exp(1j * phase)) / (amp.sum() + 1e-12))

    node_feat = np.column_stack([
        amp,
        phase,
        node_signed_flow,
        node_abs_flow,
        node_out_pos,
        node_in_neg,
        node_anisotropy,
        phase_grad_abs,
        phase_alignment,
    ]).ravel()

    global_feat = [
        float(np.linalg.norm(node_signed_flow)),
        float(np.sum(node_abs_flow)),
        float(np.mean(node_abs_flow)),
        float(np.std(node_abs_flow)),
        float(np.max(node_abs_flow)) if node_abs_flow.size else 0.0,
        float(np.mean(phase_grad_abs)),
        float(np.std(phase_grad_abs)),
        float(np.mean(phase_alignment)),
        float(np.std(phase_alignment)),
        float(coherence),
        float(amp_coherence),
        float(np.sum(np.maximum(current, 0.0))),
        float(np.sum(np.maximum(-current, 0.0))),
        float(np.linalg.norm(current, ord="fro")),
        float(np.linalg.norm(abs_current, ord="fro")),
    ]
    return np.nan_to_num(np.concatenate([node_feat, np.asarray(global_feat, dtype=float)]), nan=0.0, posinf=0.0, neginf=0.0)


def orientation_state_matrix(
    graph: Any,
    X: np.ndarray,
    tau_values: Sequence[float],
    gamma: float,
    cfg: base.BenchConfig,
) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    wave = base.GraphWaveDynamics(generator=cfg.graph_generator, gamma=gamma)
    rows = []
    for s in range(X.shape[0]):
        _, psi0 = base.density_and_phase_from_state(X, s, cfg.use_phase_from_delta)
        op = wave.prepare(graph, s)
        A = _adjacency_for_slice(graph, s)
        feat = []
        for tau in tau_values:
            psi = wave.evolve(op, psi0, float(tau))
            feat.extend(_orientation_features(A, psi).tolist())
        rows.append(feat)
    return np.nan_to_num(np.asarray(rows, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)


def build_orientation_descriptor(bundle: dict[str, Any], X_train: np.ndarray, X_test: np.ndarray, gamma: float, cfg: base.BenchConfig) -> tuple[np.ndarray, np.ndarray]:
    return (
        orientation_state_matrix(bundle["prop_graph_train"], X_train, cfg.tau_values, gamma, cfg),
        orientation_state_matrix(bundle["prop_graph_test"], X_test, cfg.tau_values, gamma, cfg),
    )


# -----------------------------------------------------------------------------
# Descriptor cache wrappers
# -----------------------------------------------------------------------------

def cached_followup_descriptor(
    name: str,
    bundle: dict[str, Any],
    X_train: np.ndarray,
    X_test: np.ndarray,
    cfg: base.BenchConfig,
    dirs: dict[str, Path],
    builder_fn,
) -> tuple[np.ndarray, np.ndarray]:
    key = base.stable_json_hash({
        "followup_version": 3,
        "name": name,
        "graph_key": bundle["key"],
        "train_hash": base.matrix_hash(X_train),
        "test_hash": base.matrix_hash(X_test),
        "tau": cfg.tau_values,
        "phase": cfg.use_phase_from_delta,
        "generator": cfg.graph_generator,
    })
    path = dirs["descriptors"] / f"followup_descriptor_{key}.npz"
    if path.exists():
        z = np.load(path)
        return z["train"], z["test"]
    F_train, F_test = builder_fn()
    np.savez_compressed(path, train=F_train, test=F_test)
    return F_train, F_test


def graph_bundle_path_for(
    X_train: np.ndarray,
    X_test: np.ndarray,
    item_names: Sequence[str],
    dataset: str,
    window_id: int,
    horizon: int,
    kernel: base.KernelName,
    weight_scheme: base.WeightScheme,
    include_item_identity: bool,
    use_stability: bool,
    dirs: dict[str, Path],
) -> Path:
    key_obj = {
        "dataset": dataset,
        "window": window_id,
        "horizon": horizon,
        "kernel": kernel,
        "weight_scheme": weight_scheme,
        "include_item_identity": include_item_identity,
        "use_stability": use_stability,
        "train_hash": base.matrix_hash(X_train),
        "test_hash": base.matrix_hash(X_test),
        "items": list(map(str, item_names)),
    }
    key = base.stable_json_hash(key_obj)
    return dirs["graphs"] / f"graph_bundle_{key}.pkl"


# -----------------------------------------------------------------------------
# Alternative evaluator stack
# -----------------------------------------------------------------------------

def evaluate_feature_block(
    dataset: str,
    window_id: int,
    horizon: int,
    descriptor_name: str,
    F_train: np.ndarray,
    F_test: np.ndarray,
    X_train: np.ndarray,
    X_test: np.ndarray,
    Y_train: np.ndarray,
    Y_test: np.ndarray,
    cfg: base.BenchConfig,
    k_values: Sequence[int],
    seed: int,
    extra: dict[str, Any],
    include_residual: bool = True,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    F_train = np.nan_to_num(np.asarray(F_train, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    F_test = np.nan_to_num(np.asarray(F_test, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)

    # 1. Direct ridge on full feature block. This avoids PCA entirely.
    Ztr, Zte = standardize_features(F_train, F_test)
    rows.append(eval_Z(dataset, window_id, horizon, descriptor_name, "direct_ridge", Ztr.shape[1],
                       Ztr, Zte, X_train, X_test, Y_train, Y_test, cfg, "next_state", extra))
    if include_residual:
        rows.append(eval_Z(dataset, window_id, horizon, descriptor_name, "direct_ridge", Ztr.shape[1],
                           Ztr, Zte, X_train, X_test, Y_train, Y_test, cfg, "residual_vs_persistence", extra))

    for k in k_values:
        k_eff = min(max(1, int(k)), F_train.shape[1], max(1, F_train.shape[0] - 1))

        # 2. Original-style SVD/PCA compression + ridge.
        Ztr, Zte = base.svd_project_train_test(F_train, F_test, k_eff)
        rows.append(eval_Z(dataset, window_id, horizon, descriptor_name, "svd_ridge", k_eff,
                           Ztr, Zte, X_train, X_test, Y_train, Y_test, cfg, "next_state", extra))
        if include_residual:
            rows.append(eval_Z(dataset, window_id, horizon, descriptor_name, "svd_ridge", k_eff,
                               Ztr, Zte, X_train, X_test, Y_train, Y_test, cfg, "residual_vs_persistence", extra))

        # 3. Supervised feature budget. Train-only top-k correlation.
        Ztr, Zte, idx = select_topcorr(F_train, F_test, Y_train, k_eff)
        rows.append(eval_Z(dataset, window_id, horizon, descriptor_name, "topcorr_ridge", k_eff,
                           Ztr, Zte, X_train, X_test, Y_train, Y_test, cfg, "next_state", {**extra, "selected_features": int(len(idx))}))
        if include_residual:
            rows.append(eval_Z(dataset, window_id, horizon, descriptor_name, "topcorr_ridge", k_eff,
                               Ztr, Zte, X_train, X_test, Y_train, Y_test, cfg, "residual_vs_persistence", {**extra, "selected_features": int(len(idx))}))

        # 4. Nonlinear ELM/tanh random features + ridge.
        Ztr, Zte = random_tanh_features(F_train, F_test, width=k_eff, seed=seed + 1009 * k_eff)
        rows.append(eval_Z(dataset, window_id, horizon, descriptor_name, "elm_tanh_ridge", k_eff,
                           Ztr, Zte, X_train, X_test, Y_train, Y_test, cfg, "next_state", extra))
        if include_residual:
            rows.append(eval_Z(dataset, window_id, horizon, descriptor_name, "elm_tanh_ridge", k_eff,
                               Ztr, Zte, X_train, X_test, Y_train, Y_test, cfg, "residual_vs_persistence", extra))

        # 5. Nonlinear random Fourier/RBF random features + ridge.
        Ztr, Zte = rff_rbf_features(F_train, F_test, width=k_eff, seed=seed + 2003 * k_eff)
        rows.append(eval_Z(dataset, window_id, horizon, descriptor_name, "rff_rbf_ridge", k_eff,
                           Ztr, Zte, X_train, X_test, Y_train, Y_test, cfg, "next_state", extra))
        if include_residual:
            rows.append(eval_Z(dataset, window_id, horizon, descriptor_name, "rff_rbf_ridge", k_eff,
                               Ztr, Zte, X_train, X_test, Y_train, Y_test, cfg, "residual_vs_persistence", extra))

    return rows


# -----------------------------------------------------------------------------
# Main follow-up run
# -----------------------------------------------------------------------------

def descriptor_blocks_for_bundle(
    bundle: dict[str, Any],
    X_train: np.ndarray,
    X_test: np.ndarray,
    item_names: Sequence[str],
    cfg: base.BenchConfig,
    dirs: dict[str, Path],
) -> list[tuple[str, np.ndarray, np.ndarray, dict[str, Any]]]:
    """Build/load old and new descriptor families for one cached graph bundle."""
    blocks: list[tuple[str, np.ndarray, np.ndarray, dict[str, Any]]] = []

    F_static_tr, F_static_te = cached_followup_descriptor(
        "static_via_original_emitter",
        bundle, X_train, X_test, cfg, dirs,
        lambda: base.build_static_descriptor(bundle, item_names),
    )
    blocks.append(("pf_static", F_static_tr, F_static_te, {"descriptor_family": "static"}))

    F_arith_tr, F_arith_te = cached_followup_descriptor(
        "arithmetic_direct_v3",
        bundle, X_train, X_test, cfg, dirs,
        lambda: build_arithmetic_direct_descriptor(bundle),
    )
    blocks.append(("pf_arithmetic_direct", F_arith_tr, F_arith_te, {"descriptor_family": "arithmetic_direct"}))
    blocks.append(("pf_static_plus_arithmetic", np.hstack([F_static_tr, F_arith_tr]), np.hstack([F_static_te, F_arith_te]), {"descriptor_family": "static+arithmetic_direct"}))

    if cfg.include_diffusion:
        Dtr, Dte = cached_followup_descriptor(
            "diffusion_original",
            bundle, X_train, X_test, cfg, dirs,
            lambda: (
                base.dynamics_state_matrix(bundle["prop_graph_train"], X_train, cfg.tau_values, "diffusion", 0.0, cfg),
                base.dynamics_state_matrix(bundle["prop_graph_test"], X_test, cfg.tau_values, "diffusion", 0.0, cfg),
            ),
        )
        blocks.append(("pf_diffusion", Dtr, Dte, {"descriptor_family": "diffusion", "gamma": 0.0}))
        blocks.append(("pf_arithmetic_plus_diffusion", np.hstack([F_arith_tr, Dtr]), np.hstack([F_arith_te, Dte]), {"descriptor_family": "arithmetic+diffusion", "gamma": 0.0}))

    if cfg.include_unitary:
        Utr, Ute = cached_followup_descriptor(
            "unitary_wave_original",
            bundle, X_train, X_test, cfg, dirs,
            lambda: (
                base.dynamics_state_matrix(bundle["prop_graph_train"], X_train, cfg.tau_values, "unitary", 0.0, cfg),
                base.dynamics_state_matrix(bundle["prop_graph_test"], X_test, cfg.tau_values, "unitary", 0.0, cfg),
            ),
        )
        Otr, Ote = cached_followup_descriptor(
            "unitary_orientation_v3",
            bundle, X_train, X_test, cfg, dirs,
            lambda: build_orientation_descriptor(bundle, X_train, X_test, 0.0, cfg),
        )
        blocks.append(("pf_unitary_wave", Utr, Ute, {"descriptor_family": "unitary_wave", "gamma": 0.0}))
        blocks.append(("pf_unitary_orientation", Otr, Ote, {"descriptor_family": "unitary_orientation", "gamma": 0.0}))
        blocks.append(("pf_unitary_wave_plus_orientation", np.hstack([Utr, Otr]), np.hstack([Ute, Ote]), {"descriptor_family": "unitary_wave+orientation", "gamma": 0.0}))
        blocks.append(("pf_arithmetic_plus_unitary_orientation", np.hstack([F_arith_tr, Otr]), np.hstack([F_arith_te, Ote]), {"descriptor_family": "arithmetic+unitary_orientation", "gamma": 0.0}))

    for gamma in cfg.damped_gammas:
        Gtr, Gte = cached_followup_descriptor(
            f"damped_wave_original_g{gamma}",
            bundle, X_train, X_test, cfg, dirs,
            lambda gamma=gamma: (
                base.dynamics_state_matrix(bundle["prop_graph_train"], X_train, cfg.tau_values, "damped", float(gamma), cfg),
                base.dynamics_state_matrix(bundle["prop_graph_test"], X_test, cfg.tau_values, "damped", float(gamma), cfg),
            ),
        )
        Otr, Ote = cached_followup_descriptor(
            f"damped_orientation_v3_g{gamma}",
            bundle, X_train, X_test, cfg, dirs,
            lambda gamma=gamma: build_orientation_descriptor(bundle, X_train, X_test, float(gamma), cfg),
        )
        blocks.append((f"pf_damped_wave_g{gamma}", Gtr, Gte, {"descriptor_family": "damped_wave", "gamma": float(gamma)}))
        blocks.append((f"pf_damped_orientation_g{gamma}", Otr, Ote, {"descriptor_family": "damped_orientation", "gamma": float(gamma)}))
        blocks.append((f"pf_damped_wave_plus_orientation_g{gamma}", np.hstack([Gtr, Otr]), np.hstack([Gte, Ote]), {"descriptor_family": "damped_wave+orientation", "gamma": float(gamma)}))
        blocks.append((f"pf_arithmetic_plus_damped_orientation_g{gamma}", np.hstack([F_arith_tr, Otr]), np.hstack([F_arith_te, Ote]), {"descriptor_family": "arithmetic+damped_orientation", "gamma": float(gamma)}))
        blocks.append((f"pf_static_plus_arithmetic_plus_damped_orientation_g{gamma}", np.hstack([F_static_tr, F_arith_tr, Otr]), np.hstack([F_static_te, F_arith_te, Ote]), {"descriptor_family": "static+arithmetic+damped_orientation", "gamma": float(gamma)}))

    return blocks


def run_followup_dataset(
    dataset_name: str,
    X_df: pd.DataFrame,
    cfg: base.BenchConfig,
    dirs: dict[str, Path],
    only_cached_graphs: bool,
    max_pf_configs: Optional[int] = None,
) -> pd.DataFrame:
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

            # Redo original non-PF baselines with alternative eval framing.
            # Baselines are included so result CSV can rank everything together.
            m = base.regression_metrics_matrix(Y_test, X_test)
            rows.append({
                "dataset": dataset_name, "window": w_id, "horizon": horizon,
                "descriptor": "baseline_persistence", "eval_method": "persistence",
                "target_mode": "next_state", "k": X_train.shape[1], "feature_dim_in": X_train.shape[1],
                "forecast_r2": m["r2"], "forecast_rmse": m["rmse"], "forecast_mae": m["mae"],
                "forecast_corr": m["mean_channel_corr"], "reconstruction_r2": 1.0,
                "reconstruction_rmse": 0.0, "corr_matrix_error": 0.0,
                "compression_ratio": 1.0,
            })

            # Raw full lag with nonlinear eval methods too.
            rows.extend(evaluate_feature_block(
                dataset_name, w_id, horizon, "raw_lag", X_train, X_test,
                X_train, X_test, Y_train, Y_test, cfg, cfg.k_values,
                seed=cfg.seed + 17 * w_id, extra={"descriptor_family": "raw_lag"}, include_residual=True,
            ))

            pf_count = 0
            for kernel in cfg.kernels:
                for weight_scheme in cfg.weight_schemes:
                    for include_item in cfg.item_identity_options:
                        for use_stability in cfg.stability_options:
                            if max_pf_configs is not None and pf_count >= max_pf_configs:
                                break
                            path = graph_bundle_path_for(
                                X_train, X_test, item_names, dataset_name, w_id, horizon,
                                kernel, weight_scheme, include_item, use_stability, dirs,
                            )
                            if only_cached_graphs and not path.exists():
                                print(f"[SKIP] missing cached graph: {path.name}")
                                continue
                            try:
                                bundle = base.build_or_load_pf_bundle(
                                    X_train=X_train,
                                    X_test=X_test,
                                    item_names=item_names,
                                    dataset=dataset_name,
                                    window_id=w_id,
                                    horizon=horizon,
                                    kernel=kernel,
                                    weight_scheme=weight_scheme,
                                    include_item_identity=include_item,
                                    use_stability=use_stability,
                                    cfg=cfg,
                                    dirs=dirs,
                                )
                                pf_count += 1
                                combo_extra = {
                                    "kernel": kernel,
                                    "weight_scheme": weight_scheme,
                                    "item_identity": bool(include_item),
                                    "stability": bool(use_stability),
                                }
                                blocks = descriptor_blocks_for_bundle(bundle, X_train, X_test, item_names, cfg, dirs)
                                for desc_name, Ftr, Fte, desc_extra in blocks:
                                    rows.extend(evaluate_feature_block(
                                        dataset_name, w_id, horizon, desc_name, Ftr, Fte,
                                        X_train, X_test, Y_train, Y_test, cfg, cfg.k_values,
                                        seed=cfg.seed + 10000 * w_id + 97 * pf_count,
                                        extra={**combo_extra, **desc_extra},
                                        include_residual=True,
                                    ))
                            except Exception as e:
                                print("[WARN] follow-up config failed:", {
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

    return pd.DataFrame(rows)


def aggregate_followup(df: pd.DataFrame) -> pd.DataFrame:
    group_cols = [
        "dataset", "descriptor", "eval_method", "target_mode", "k",
    ]
    optional = ["kernel", "weight_scheme", "item_identity", "stability", "descriptor_family", "gamma"]
    group_cols += [c for c in optional if c in df.columns]
    metric_cols = [
        "forecast_r2", "forecast_rmse", "forecast_mae", "forecast_corr",
        "reconstruction_r2", "reconstruction_rmse", "corr_matrix_error",
    ]
    metric_cols = [c for c in metric_cols if c in df.columns]
    out = df.groupby(group_cols, dropna=False)[metric_cols].agg(["mean", "std", "median", "count"])
    out.columns = ["_".join(c).strip("_") for c in out.columns.to_flat_index()]
    out = out.reset_index()
    if "forecast_r2_mean" in out.columns:
        out = out.sort_values(["dataset", "forecast_r2_mean"], ascending=[True, False])
    return out


def run_followup(cfg: base.BenchConfig, only_cached_graphs: bool, max_pf_configs: Optional[int]) -> pd.DataFrame:
    dirs = base.ensure_dirs(cfg)
    all_results = []
    for dataset_name in cfg.datasets:
        print("\n" + "=" * 100)
        print(f"[FOLLOWUP DATASET] {dataset_name}")
        X_df = base.cached_dataset(dataset_name, cfg, dirs)
        print(f"[SHAPE] {dataset_name}: {X_df.shape}")
        res = run_followup_dataset(dataset_name, X_df, cfg, dirs, only_cached_graphs, max_pf_configs=max_pf_configs)
        out_path = dirs["results"] / f"followup_results_{dataset_name}_{cfg.preset}.csv"
        res.to_csv(out_path, index=False)
        print(f"[SAVE] {out_path}")
        all_results.append(res)

    combined = pd.concat(all_results, ignore_index=True)
    combined_path = dirs["results"] / f"followup_results_combined_{cfg.preset}.csv"
    combined.to_csv(combined_path, index=False)
    agg = aggregate_followup(combined)
    agg_path = dirs["results"] / f"followup_results_aggregate_{cfg.preset}.csv"
    agg.to_csv(agg_path, index=False)

    print("\n[SAVE]", combined_path)
    print("[SAVE]", agg_path)
    print("\nTop rows by dataset:")
    print(agg.groupby("dataset", group_keys=False).head(30).to_string(index=False))
    return combined


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--preset", default="lite", choices=["lite", "markets", "full"])
    p.add_argument("--cache-dir", default=None)
    p.add_argument("--datasets", nargs="*", default=None,
                   choices=["plotly_stock_returns", "uci_household_power", "uci_appliances"])
    p.add_argument("--only-cached-graphs", action="store_true",
                   help="Skip PF configs whose graph bundles have not already been cached by pf_better_benchmark.py.")
    p.add_argument("--max-pf-configs", type=int, default=None,
                   help="Debug/smoke-test limit per window. Example: --max-pf-configs 2")
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

    print("[FOLLOWUP CONFIG]")
    print(json.dumps(asdict(cfg), indent=2, default=str))
    print(f"[ONLY CACHED GRAPHS] {args.only_cached_graphs}")
    print(f"[MAX PF CONFIGS] {args.max_pf_configs}")
    run_followup(cfg, only_cached_graphs=args.only_cached_graphs, max_pf_configs=args.max_pf_configs)


if __name__ == "__main__":
    main()
