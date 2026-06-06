"""
pf_targeted_flux_test_purged.py
================================

Leak-safe / purged version of pf_targeted_flux_test.py.

This script runs the same targeted flux/continuum tests, but replaces the default
walk-forward splitter with a purged splitter:

    train feature: 0 ... end-1
    train target:  horizon ... end-1+horizon
    test feature:  end+embargo ... end+embargo+test_len-1

With the default embargo="auto", embargo=horizon, so the first test feature row
is strictly after the last train target row. This prevents the boundary overlap:

    last train target X[end] == first test feature X[end]

that appeared in the earlier audit for horizon=1.

Recommended runs
----------------
# First just audit the purged windows. This does not build graphs.
python pf_targeted_flux_test_purged.py --audit-only

# Build/test the purged targeted grid. Existing unpurged graph bundles will not be reused
# because the train/test hashes change.
python pf_targeted_flux_test_purged.py --no-stability

# Smoke build for one or two graph configs.
python pf_targeted_flux_test_purged.py --no-stability --max-pf-configs 2

# Broader validation after the smoke pass.
python pf_targeted_flux_test_purged.py --no-stability --n-windows 5
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


def walk_forward_windows_purged(
    n: int,
    train_frac: float,
    test_frac: float,
    n_windows: int,
    embargo: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Expanding walk-forward windows with an embargo gap after train features.

    The default base splitter uses test_start=end. That is safe for feature/feature
    separation but not for horizon labels, because train target at end-1+horizon can
    equal the first test feature at end. Here test_start=end+embargo.
    """
    embargo = max(0, int(embargo))
    test_len = max(10, int(round(n * test_frac)))
    train_len0 = max(20, int(round(n * train_frac)))
    windows: list[tuple[np.ndarray, np.ndarray]] = []

    # The train end is exclusive. Need end + embargo + test_len <= n.
    max_train_end = n - embargo - test_len
    if max_train_end <= 20:
        return []
    train_len0 = min(train_len0, max_train_end)

    if n_windows <= 1:
        train_ends = [train_len0]
    else:
        train_ends = np.linspace(train_len0, max_train_end, n_windows).astype(int).tolist()

    # Remove duplicates caused by integer linspace on short arrays.
    seen = set()
    for end in train_ends:
        end = int(end)
        if end in seen:
            continue
        seen.add(end)
        test_start = end + embargo
        test_end = min(n, test_start + test_len)
        train_idx = np.arange(0, end)
        test_idx = np.arange(test_start, test_end)
        if len(train_idx) >= 20 and len(test_idx) >= 5:
            windows.append((train_idx, test_idx))
    return windows


def install_purged_splitter(embargo: int) -> None:
    """Monkey-patch the shared base splitter used by pf_flux_continuum_eval."""
    def _splitter(n: int, train_frac: float, test_frac: float, n_windows: int):
        return walk_forward_windows_purged(n, train_frac, test_frac, n_windows, embargo=embargo)

    base.walk_forward_windows = _splitter  # type: ignore[assignment]
    flux.base.walk_forward_windows = _splitter  # same module object, explicit for clarity


def targeted_config(args: argparse.Namespace) -> base.BenchConfig:
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
    return base.BenchConfig(**{
        **asdict(cfg),
        "preset": "lite_purged" if args.separate_preset_cache else "lite",
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


def embargo_for_horizon(args: argparse.Namespace, horizon: int) -> int:
    if args.embargo == "auto":
        return int(horizon)
    return int(args.embargo)


def write_time_audit(cfg: base.BenchConfig, dirs: dict[str, Path], args: argparse.Namespace) -> pd.DataFrame:
    rows = []
    for dataset_name in cfg.datasets:
        X_df = base.cached_dataset(dataset_name, cfg, dirs)
        X_raw = X_df.to_numpy(dtype=float)
        for horizon in cfg.horizons:
            emb = embargo_for_horizon(args, int(horizon))
            X0, _ = base.supervised_pairs(X_raw, int(horizon))
            windows = walk_forward_windows_purged(len(X0), cfg.train_frac, cfg.test_frac, cfg.n_windows, embargo=emb)
            for w_id, (train_idx, test_idx) in enumerate(windows):
                train_start = int(train_idx[0]) if len(train_idx) else -1
                train_end = int(train_idx[-1]) if len(train_idx) else -1
                test_start = int(test_idx[0]) if len(test_idx) else -1
                test_end = int(test_idx[-1]) if len(test_idx) else -1
                train_target_end = train_end + int(horizon)
                rows.append({
                    "dataset": dataset_name,
                    "horizon": int(horizon),
                    "embargo": int(emb),
                    "window": int(w_id),
                    "n_rows_X0": int(len(X0)),
                    "n_train": int(len(train_idx)),
                    "n_test": int(len(test_idx)),
                    "train_feature_start_t": train_start,
                    "train_feature_end_t": train_end,
                    "train_target_start_t": train_start + int(horizon),
                    "train_target_end_t": train_target_end,
                    "test_feature_start_t": test_start,
                    "test_feature_end_t": test_end,
                    "test_target_start_t": test_start + int(horizon),
                    "test_target_end_t": test_end + int(horizon),
                    "train_before_test_features": bool(train_end < test_start),
                    "train_targets_end_before_test_features_start": bool(train_target_end < test_start),
                    "purge_gap_after_train_feature_end": int(test_start - train_end - 1),
                    "purge_gap_after_train_target_end": int(test_start - train_target_end - 1),
                    "feature_target_horizon": int(horizon),
                    "standardization_fit_on": "train feature rows only",
                    "encoder_fit_on": "train feature rows only",
                    "test_graph_uses": "test feature slice X_t only; target Y_t+h is label only",
                })
    audit = pd.DataFrame(rows)
    out = dirs["results"] / f"targeted_flux_time_audit_{cfg.preset}_embargo{args.embargo}.csv"
    audit.to_csv(out, index=False)
    print(f"[SAVE] {out}")
    if not audit.empty:
        ok = bool(audit["train_before_test_features"].all())
        target_ok = bool(audit["train_targets_end_before_test_features_start"].all())
        print("[TIME AUDIT] train feature windows before test feature windows:", ok)
        print("[TIME AUDIT] train target end before test feature start:", target_ok)
        print(audit.to_string(index=False))
    else:
        print("[TIME AUDIT] no valid windows produced")
    return audit


def summarize_focus(cfg: base.BenchConfig, combined: pd.DataFrame, dirs: dict[str, Path], args: argparse.Namespace) -> pd.DataFrame:
    if combined.empty:
        print("[FOCUS SUMMARY] no prediction rows")
        return combined
    agg = flux.aggregate_prediction(combined)
    focus = agg[agg["target_name"].isin(FOCUS_TARGETS)].copy()
    focus["is_focus_descriptor"] = focus["descriptor"].isin(FOCUS_DESCRIPTORS)
    focus = focus.sort_values(["dataset", "target_name", "forecast_r2_mean"], ascending=[True, True, False])
    out = dirs["results"] / f"targeted_flux_focus_summary_{cfg.preset}_embargo{args.embargo}.csv"
    focus.to_csv(out, index=False)
    print(f"[SAVE] {out}")
    if not focus.empty:
        print("\n[FOCUS SUMMARY: top 12 per target]")
        show_cols = [
            "dataset", "target_name", "descriptor", "eval_method", "k", "kernel", "weight_scheme",
            "stability", "gamma", "forecast_r2_mean", "forecast_r2_std", "forecast_r2_count",
            "forecast_corr_mean", "rank_spearman_mean", "precision_top_decile_mean", "auc_top_quartile_mean",
        ]
        show_cols = [c for c in show_cols if c in focus.columns]
        print(focus.groupby(["dataset", "target_name"], group_keys=False).head(12)[show_cols].to_string(index=False))
    return focus


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default="pf_bench_cache")
    p.add_argument("--datasets", nargs="*", default=["plotly_stock_returns"], choices=["plotly_stock_returns", "uci_household_power", "uci_appliances"])
    p.add_argument("--only-cached-graphs", action="store_true", help="Skip missing purged graph bundles instead of building them.")
    p.add_argument("--max-pf-configs", type=int, default=None, help="Optional cap per window for smoke tests.")
    p.add_argument("--audit-only", action="store_true", help="Only write the purged time audit; do not run flux models.")
    p.add_argument("--n-windows", type=int, default=3)
    p.add_argument("--max-rows-stocks", type=int, default=900)
    p.add_argument("--horizons", nargs="*", type=int, default=[1])
    p.add_argument("--embargo", default="auto", help="Rows between train feature end and test feature start. Use 'auto' for embargo=horizon.")
    p.add_argument("--k-values", nargs="*", type=int, default=[2, 4, 8])
    p.add_argument("--kernels", nargs="*", default=None, choices=["log_gcd", "log_lcm", "arithmetic_distance", "jaccard_gcd_lcm", "coprime_indicator"])
    p.add_argument("--weight-schemes", nargs="*", default=None, choices=["log_prime_weights", "uniform_weights", "random_monotone_weights", "shuffled_mode_weights"])
    p.add_argument("--no-stability", action="store_true")
    p.add_argument("--include-unitary", action="store_true", help="Include unitary wave flux blocks; off by default for speed.")
    p.add_argument("--no-phase", action="store_true")
    p.add_argument("--tau-values", nargs="*", type=float, default=[0.25])
    p.add_argument("--damped-gammas", nargs="*", type=float, default=[0.001, 0.005, 0.01])
    p.add_argument("--ridge-alpha", type=float, default=10.0)
    p.add_argument("--separate-preset-cache", action="store_true", help="Use dataset/result suffix lite_purged instead of lite. Graph hashes still isolate purged splits either way.")
    args = p.parse_args()
    if args.embargo != "auto":
        try:
            int(args.embargo)
        except ValueError as e:
            raise SystemExit("--embargo must be 'auto' or an integer") from e
    return args


def main() -> None:
    args = parse_args()
    # Use the smallest horizon for the monkey patch only when multiple horizons are present.
    # The audit itself uses per-horizon embargo; the runner calls one splitter signature without
    # horizon, so for multi-horizon runs use max requested embargo to remain conservative.
    if args.embargo == "auto":
        runner_embargo = max(int(h) for h in args.horizons)
    else:
        runner_embargo = int(args.embargo)
    install_purged_splitter(runner_embargo)

    cfg = targeted_config(args)
    dirs = base.ensure_dirs(cfg)
    print("[TARGETED FLUX PURGED CONFIG]")
    print(json.dumps({**asdict(cfg), "embargo": args.embargo, "runner_embargo": runner_embargo}, indent=2, default=str))
    print("[ONLY CACHED GRAPHS]", bool(args.only_cached_graphs))
    print("[MAX PF CONFIGS]", args.max_pf_configs)
    audit = write_time_audit(cfg, dirs, args)
    if args.audit_only:
        if not audit.empty and not bool(audit["train_targets_end_before_test_features_start"].all()):
            raise SystemExit("Purged audit failed: train target overlaps test feature start.")
        return
    if not audit.empty and not bool(audit["train_targets_end_before_test_features_start"].all()):
        raise SystemExit("Purged audit failed: refusing to run models.")

    combined = flux.run_flux(cfg, only_cached_graphs=bool(args.only_cached_graphs), max_pf_configs=args.max_pf_configs)
    summarize_focus(cfg, combined, dirs, args)


if __name__ == "__main__":
    main()
