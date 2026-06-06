"""
pf_targeted_flux_test.py
=========================

Targeted expansion around the promising flux/continuum results.

Purpose
-------
Run a narrower, more decisive test of the empirical lead:

    graph/continuum flux features -> future acceleration dispersion,
    acceleration norm, dispersion change, and second-order volatility proxies.

It reuses the existing pf_bench_cache layout and can either:
  - use only cached graph bundles, or
  - build missing targeted bundles and cache them.

It also writes a simple time-leakage audit showing the train/test index ranges and
confirming that every test window starts strictly after its train window.

Recommended first runs
----------------------
# Audit cache and run only existing graphs
python pf_targeted_flux_test.py --only-cached-graphs --max-pf-configs 2

# Build/test the targeted grid
python pf_targeted_flux_test.py

# Faster targeted grid, no stability ablation
python pf_targeted_flux_test.py --no-stability

# More windows after the first pass looks sane
python pf_targeted_flux_test.py --n-windows 5
"""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

import pf_better_benchmark as base
import pf_flux_continuum_eval as flux


FOCUS_TARGETS = [
    "future_acceleration_dispersion",
    "future_acceleration_norm",
    "future_dispersion_change",
    "future_abs_return_mean",
    "future_squared_return_mean",
    "future_residual_return",
    "future_rank_return",
]

FOCUS_DESCRIPTORS = [
    "pf_raw_continuum_flux",
    "pf_arithmetic_direct",
    "pf_arithmetic_plus_raw_continuum_flux",
    "pf_diffusion_flux",
    "pf_damped_wave_flux_g0.001",
    "pf_damped_wave_flux_g0.005",
    "pf_damped_wave_flux_g0.01",
    "pf_arithmetic_plus_damped_wave_flux_g0.001",
    "pf_arithmetic_plus_damped_wave_flux_g0.005",
    "pf_arithmetic_plus_damped_wave_flux_g0.01",
]


def targeted_config(args: argparse.Namespace) -> base.BenchConfig:
    # Keep preset='lite' so this reuses pf_bench_cache/datasets/*_lite.pkl and writes to familiar result paths.
    cfg = base.config_for_preset("lite")
    kernels = tuple(args.kernels) if args.kernels else (
        "log_gcd",
        "arithmetic_distance",
        "jaccard_gcd_lcm",
    )
    weights = tuple(args.weight_schemes) if args.weight_schemes else (
        "log_prime_weights",
        "shuffled_mode_weights",
        "uniform_weights",
    )
    stability = (False,) if args.no_stability else (False, True)
    cfg = base.BenchConfig(**{
        **asdict(cfg),
        "preset": "lite",
        "datasets": tuple(args.datasets),
        "k_values": tuple(args.k_values),
        "horizons": tuple(args.horizons),
        "n_windows": int(args.n_windows),
        "kernels": kernels,
        "weight_schemes": weights,
        "item_identity_options": (False,),
        "stability_options": stability,
        "tau_values": tuple(args.tau_values),
        "damped_gammas": tuple(args.damped_gammas),
        "include_unitary": bool(args.include_unitary),
        "include_diffusion": True,
        "use_phase_from_delta": not bool(args.no_phase),
        "max_rows_stocks": int(args.max_rows_stocks),
        "cache_dir": args.cache_dir,
        "ridge_alpha": float(args.ridge_alpha),
    })
    return cfg


def write_time_audit(cfg: base.BenchConfig, dirs: dict[str, Path]) -> pd.DataFrame:
    rows = []
    for dataset_name in cfg.datasets:
        X_df = base.cached_dataset(dataset_name, cfg, dirs)
        X_raw = X_df.to_numpy(dtype=float)
        for horizon in cfg.horizons:
            X0, Y = base.supervised_pairs(X_raw, horizon)
            windows = base.walk_forward_windows(len(X0), cfg.train_frac, cfg.test_frac, cfg.n_windows)
            for w_id, (train_idx, test_idx) in enumerate(windows):
                train_start = int(train_idx[0]) if len(train_idx) else -1
                train_end = int(train_idx[-1]) if len(train_idx) else -1
                test_start = int(test_idx[0]) if len(test_idx) else -1
                test_end = int(test_idx[-1]) if len(test_idx) else -1
                rows.append({
                    "dataset": dataset_name,
                    "horizon": int(horizon),
                    "window": int(w_id),
                    "n_rows_X0": int(len(X0)),
                    "n_train": int(len(train_idx)),
                    "n_test": int(len(test_idx)),
                    "train_feature_start_t": train_start,
                    "train_feature_end_t": train_end,
                    "train_target_start_t": train_start + int(horizon),
                    "train_target_end_t": train_end + int(horizon),
                    "test_feature_start_t": test_start,
                    "test_feature_end_t": test_end,
                    "test_target_start_t": test_start + int(horizon),
                    "test_target_end_t": test_end + int(horizon),
                    "train_before_test_features": bool(train_end < test_start),
                    "train_targets_end_before_test_features_start": bool((train_end + int(horizon)) < test_start),
                    "feature_target_horizon": int(horizon),
                    "standardization_fit_on": "train feature rows only",
                    "encoder_fit_on": "train feature rows only",
                    "test_graph_uses": "test feature slice X_t only; target Y_t+h is label only",
                })
    audit = pd.DataFrame(rows)
    out = dirs["results"] / f"targeted_flux_time_audit_{cfg.preset}.csv"
    audit.to_csv(out, index=False)
    print(f"[SAVE] {out}")
    if not audit.empty:
        ok = bool(audit["train_before_test_features"].all())
        target_ok = bool(audit["train_targets_end_before_test_features_start"].all())
        print("[TIME AUDIT] train feature windows before test feature windows:", ok)
        print("[TIME AUDIT] train target end before test feature start:", target_ok)
        print(audit.to_string(index=False))
    return audit


def summarize_focus(cfg: base.BenchConfig, combined: pd.DataFrame, dirs: dict[str, Path]) -> pd.DataFrame:
    if combined.empty:
        return combined
    agg = flux.aggregate_prediction(combined)
    focus = agg[agg["target_name"].isin(FOCUS_TARGETS)].copy()
    # Keep all descriptors, but surface focus descriptors first by sorting.
    focus["is_focus_descriptor"] = focus["descriptor"].isin(FOCUS_DESCRIPTORS)
    sort_cols = ["dataset", "target_name", "forecast_r2_mean"]
    focus = focus.sort_values(sort_cols, ascending=[True, True, False])
    out = dirs["results"] / f"targeted_flux_focus_summary_{cfg.preset}.csv"
    focus.to_csv(out, index=False)
    print(f"[SAVE] {out}")
    if not focus.empty:
        print("\n[FOCUS SUMMARY: top 12 per target]")
        show_cols = [
            "dataset", "target_name", "descriptor", "eval_method", "k", "kernel", "weight_scheme",
            "stability", "gamma", "forecast_r2_mean", "forecast_corr_mean", "rank_spearman_mean",
            "precision_top_decile_mean", "auc_top_quartile_mean",
        ]
        show_cols = [c for c in show_cols if c in focus.columns]
        print(focus.groupby(["dataset", "target_name"], group_keys=False).head(12)[show_cols].to_string(index=False))
    return focus


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default="pf_bench_cache")
    p.add_argument("--datasets", nargs="*", default=["plotly_stock_returns"], choices=["plotly_stock_returns", "uci_household_power", "uci_appliances"])
    p.add_argument("--only-cached-graphs", action="store_true", help="Skip missing graph bundles instead of building them.")
    p.add_argument("--max-pf-configs", type=int, default=None, help="Optional cap per window for smoke tests.")
    p.add_argument("--n-windows", type=int, default=3)
    p.add_argument("--max-rows-stocks", type=int, default=900)
    p.add_argument("--horizons", nargs="*", type=int, default=[1])
    p.add_argument("--k-values", nargs="*", type=int, default=[2, 4, 8])
    p.add_argument("--kernels", nargs="*", default=None, choices=["log_gcd", "log_lcm", "arithmetic_distance", "jaccard_gcd_lcm", "coprime_indicator"])
    p.add_argument("--weight-schemes", nargs="*", default=None, choices=["log_prime_weights", "uniform_weights", "random_monotone_weights", "shuffled_mode_weights"])
    p.add_argument("--no-stability", action="store_true")
    p.add_argument("--include-unitary", action="store_true", help="Include unitary wave flux blocks; off by default for speed.")
    p.add_argument("--no-phase", action="store_true")
    p.add_argument("--tau-values", nargs="*", type=float, default=[0.25])
    p.add_argument("--damped-gammas", nargs="*", type=float, default=[0.001, 0.005, 0.01])
    p.add_argument("--ridge-alpha", type=float, default=10.0)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    cfg = targeted_config(args)
    dirs = base.ensure_dirs(cfg)
    print("[TARGETED FLUX CONFIG]")
    print(json.dumps(asdict(cfg), indent=2, default=str))
    print("[ONLY CACHED GRAPHS]", bool(args.only_cached_graphs))
    print("[MAX PF CONFIGS]", args.max_pf_configs)
    write_time_audit(cfg, dirs)
    combined = flux.run_flux(cfg, only_cached_graphs=bool(args.only_cached_graphs), max_pf_configs=args.max_pf_configs)
    summarize_focus(cfg, combined, dirs)


if __name__ == "__main__":
    main()
