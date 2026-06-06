"""
macro_graph_entry_exit_nn.py
================================

No-lookahead Prime-Fock + macro/credit graph + eigen/flux + multi-head NN
benchmark for entry/exit prediction.

Pipeline implemented:
  1. Load stock daily data
  2. Build normalized stock features
  3. Add macro/credit proxy nodes
  4. Build stock-stock and stock-macro Prime-Fock graph
  5. Extract arithmetic, field/flux, dynamics, and eigen features
  6. Create triple-barrier entry/exit labels
  7. Train MLP multi-head NN
  8. Run purged walk-forward validation
  9. Report trading metrics and ablations

Design principles:
  - No lookahead: rolling features use data available at date t; Prime-Fock binning
    is fitted inside each walk-forward fold on train rows only; scalers are fitted
    on train only; train/test split uses an embargo at least as large as the
    label horizon.
  - Macro/credit nodes are graph nodes. They are not prediction targets.
  - The NN predicts trade states, not raw forward returns.

Example smoke:
  python macro_graph_entry_exit_nn.py --smoke --epochs 3 --n-windows 1

Example real data, wide CSV with Date plus ticker columns:
  python macro_graph_entry_exit_nn.py ^
    --csv all_stocks_5yr.csv ^
    --date-col Date ^
    --max-symbols 80 ^
    --epochs 20 ^
    --out metrics.csv
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd


def _ensure_prime_fock_importable() -> None:
    """Allow running this script beside the Prime-Fock module files."""
    try:
        import prime_fock  # noqa: F401
        return
    except Exception:
        pass

    here = Path(__file__).resolve().parent
    init_py = here / "__init__.py"
    if init_py.exists():
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "prime_fock", str(init_py), submodule_search_locations=[str(here)]
        )
        if spec is None or spec.loader is None:
            raise ImportError("Could not construct import spec for local prime_fock package")
        module = importlib.util.module_from_spec(spec)
        sys.modules["prime_fock"] = module
        spec.loader.exec_module(module)
        return
    raise ImportError(
        "Could not import prime_fock. Run from project root or place this script beside the package files."
    )


_ensure_prime_fock_importable()

from prime_fock import (  # noqa: E402
    ArithmeticGraphBuilder,
    CategoricalSpec,
    DynamicsFeatureEmitter,
    GenericFeatureEmitter,
    GraphDiffusionDynamics,
    GraphWaveDynamics,
    NumericBinningSpec,
    PrimeFockNDConfig,
    PrimeFockNDEncoder,
    SharedSupportGraph,
    StabilityFilter,
    StabilityResult,
)

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
except Exception as exc:  # pragma: no cover
    torch = None
    nn = None
    F = None
    _TORCH_IMPORT_ERROR = exc
else:
    _TORCH_IMPORT_ERROR = None
    try:
        torch.set_num_threads(max(1, min(4, os.cpu_count() or 1)))
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

DEFAULT_MACRO_SYMBOLS = (
    "SPY", "QQQ", "IWM", "DIA",
    "TLT", "IEF", "SHY",
    "HYG", "LQD", "JNK",
    "GLD", "SLV",
    "USO", "UUP",
    "VIX", "^VIX", "VIXY",
)


@dataclass(frozen=True)
class RunConfig:
    csv: Optional[str] = None
    macro_csv: Optional[str] = None
    date_col: str = "Date"
    ticker_col: Optional[str] = None
    close_col: str = "Close"
    max_symbols: int = 80
    min_history: int = 260
    macro_symbols: tuple[str, ...] = DEFAULT_MACRO_SYMBOLS

    lookback_z: int = 252
    fast_lookback: int = 21
    corr_lookback: int = 63
    beta_lookback: int = 63
    winsor_z: float = 5.0

    label_horizon: int = 20
    exit_horizon: int = 10
    vol_lookback: int = 20
    take_profit_vol: float = 1.5
    stop_loss_vol: float = 1.0
    relative_to_market: bool = True
    market_symbol: str = "SPY"

    train_frac: float = 0.60
    test_frac: float = 0.20
    n_windows: int = 3
    embargo: Optional[int] = None
    min_train_dates: int = 252

    graph_kernel: str = "log_gcd"
    include_item_identity: bool = False
    use_stability: bool = True
    dynamics_taus: tuple[float, ...] = ()
    wave_gamma: float = 0.05
    top_eigs: int = 6

    hidden_dim: int = 256
    depth: int = 3
    dropout: float = 0.15
    lr: float = 1e-3
    weight_decay: float = 1e-4
    batch_size: int = 4096
    epochs: int = 20
    patience: int = 5
    seed: int = 0
    device: str = "auto"

    top_frac: float = 0.10
    cost_bps: float = 5.0
    out: str = "macro_graph_entry_exit_nn_metrics.csv"
    preds_out: str = "macro_graph_entry_exit_nn_predictions.csv"
    smoke: bool = False


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def seed_everything(seed: int) -> None:
    np.random.seed(seed)
    if torch is not None:
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)


def safe_float_array(x) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    arr[~np.isfinite(arr)] = np.nan
    return arr


def winsorize_df(df: pd.DataFrame, cols: Sequence[str], limit: float) -> pd.DataFrame:
    out = df.copy()
    for c in cols:
        out[c] = out[c].clip(-limit, limit)
    return out


def rolling_z_by_col(mat: pd.DataFrame, lookback: int, min_periods: Optional[int] = None) -> pd.DataFrame:
    minp = min_periods or max(20, lookback // 5)
    mean = mat.rolling(lookback, min_periods=minp).mean()
    std = mat.rolling(lookback, min_periods=minp).std(ddof=0)
    return (mat - mean) / std.replace(0, np.nan)


def rolling_vol(ret: pd.DataFrame, lookback: int) -> pd.DataFrame:
    return ret.rolling(lookback, min_periods=max(10, lookback // 2)).std(ddof=0)


def quantile_bins(values: np.ndarray, q: Sequence[float] = (0.1, 0.25, 0.5, 0.75, 0.9)) -> list[float]:
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return [-1.0, 0.0, 1.0]
    qs = np.quantile(v, q)
    qs = np.unique(qs)
    if qs.size == 0:
        return [-1.0, 0.0, 1.0]
    # NumericBinningSpec requires strictly increasing bins.
    out = []
    last = -np.inf
    for x in qs:
        x = float(x)
        if x > last + 1e-12:
            out.append(x)
            last = x
    if len(out) < 2:
        center = out[0] if out else 0.0
        out = [center - 1e-6, center + 1e-6]
    return out


def normalized_adjacency(A: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    A = np.asarray(A, dtype=float)
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    deg = A.sum(axis=1)
    inv = np.where(deg > eps, 1.0 / np.sqrt(deg + eps), 0.0)
    return (inv[:, None] * A) * inv[None, :]


def dense_raw_stability_result(graph: SharedSupportGraph) -> StabilityResult:
    A = np.asarray(graph.adjacency, dtype=float)
    spectral = np.zeros(A.shape[0], dtype=float)
    for s in range(A.shape[0]):
        if A.shape[1] > 0:
            spectral[s] = float(np.max(np.abs(np.linalg.eigvalsh((A[s] + A[s].T) / 2.0))))
    return StabilityResult(
        stable_adjacency=A.copy(),
        spectral_radius=spectral,
        local_stability=(A > 0).astype(float),
        slices=graph.slices,
        items=graph.items,
    )


def graph_with_adjacency(base: SharedSupportGraph, adjacency: np.ndarray, suffix: str) -> SharedSupportGraph:
    return SharedSupportGraph(
        adjacency=np.asarray(adjacency, dtype=float),
        active_modes=base.active_modes,
        slices=base.slices,
        items=base.items,
        kernel=f"{base.kernel}_{suffix}",
    )


# ---------------------------------------------------------------------------
# Data loading and daily feature panel
# ---------------------------------------------------------------------------

def generate_synthetic_prices(n_dates: int = 780, n_stocks: int = 48, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    dates = pd.bdate_range("2020-01-01", periods=n_dates)
    macro_names = ["SPY", "QQQ", "IWM", "TLT", "HYG", "LQD", "GLD", "USO", "UUP", "VIXY"]
    stock_names = [f"S{i:03d}" for i in range(n_stocks)]
    all_names = macro_names + stock_names

    # Latent macro factors with regime shifts.
    f_market = rng.normal(0, 0.008, n_dates)
    f_rates = rng.normal(0, 0.004, n_dates)
    f_credit = rng.normal(0, 0.005, n_dates)
    f_commod = rng.normal(0, 0.006, n_dates)
    for t in range(100, n_dates, 170):
        f_market[t:t+20] += rng.normal(-0.006, 0.012, min(20, n_dates-t))
        f_credit[t:t+20] += rng.normal(-0.006, 0.010, min(20, n_dates-t))

    rets = {}
    rets["SPY"] = f_market + rng.normal(0, 0.004, n_dates)
    rets["QQQ"] = 1.25 * f_market - 0.4 * f_rates + rng.normal(0, 0.006, n_dates)
    rets["IWM"] = 1.15 * f_market + 0.7 * f_credit + rng.normal(0, 0.007, n_dates)
    rets["TLT"] = -1.5 * f_rates - 0.15 * f_market + rng.normal(0, 0.006, n_dates)
    rets["HYG"] = 0.45 * f_market + 1.6 * f_credit + rng.normal(0, 0.004, n_dates)
    rets["LQD"] = 0.15 * f_market + 0.6 * f_credit - 0.4 * f_rates + rng.normal(0, 0.003, n_dates)
    rets["GLD"] = -0.2 * f_rates + 0.3 * f_commod + rng.normal(0, 0.006, n_dates)
    rets["USO"] = 1.5 * f_commod + rng.normal(0, 0.010, n_dates)
    rets["UUP"] = 0.6 * f_rates - 0.2 * f_commod + rng.normal(0, 0.003, n_dates)
    rets["VIXY"] = -1.8 * f_market - 0.5 * f_credit + rng.normal(0, 0.018, n_dates)

    sector_beta = rng.normal(0, 0.004, (6, n_dates))
    for k, name in enumerate(stock_names):
        sector = k % 6
        b_m = rng.uniform(0.5, 1.5)
        b_r = rng.normal(0, 0.5)
        b_c = rng.normal(0.2, 0.7)
        b_o = rng.normal(0, 0.4)
        alpha_pulse = np.zeros(n_dates)
        # synthetic entry/exit signal depends on macro stress and sector waves
        stress = -f_market - f_credit
        alpha_pulse[1:] = 0.002 * np.tanh(-40 * stress[:-1]) + 0.001 * np.sign(sector_beta[sector, :-1])
        rets[name] = (
            b_m * f_market + b_r * f_rates + b_c * f_credit + b_o * f_commod
            + sector_beta[sector] + alpha_pulse + rng.normal(0, 0.012, n_dates)
        )

    prices = pd.DataFrame(index=dates)
    for name in all_names:
        prices[name] = 100.0 * np.exp(np.cumsum(rets[name]))
    prices.index.name = "Date"
    return prices.reset_index()


def load_price_matrix(cfg: RunConfig) -> pd.DataFrame:
    if cfg.smoke or cfg.csv is None:
        n_stocks = max(12, min(48, cfg.max_symbols - 10 if cfg.max_symbols else 48))
        return generate_synthetic_prices(n_stocks=n_stocks, seed=cfg.seed)

    df = pd.read_csv(cfg.csv)
    if cfg.date_col not in df.columns:
        raise ValueError(f"date column {cfg.date_col!r} not found")
    df[cfg.date_col] = pd.to_datetime(df[cfg.date_col])

    # Long format: Date, ticker, Close.
    if cfg.ticker_col and cfg.ticker_col in df.columns and cfg.close_col in df.columns:
        wide = df.pivot_table(index=cfg.date_col, columns=cfg.ticker_col, values=cfg.close_col, aggfunc="last")
    else:
        numeric_cols = [c for c in df.columns if c != cfg.date_col and pd.api.types.is_numeric_dtype(df[c])]
        if not numeric_cols:
            raise ValueError("Could not infer wide price columns. Provide --ticker-col and --close-col for long data.")
        wide = df[[cfg.date_col] + numeric_cols].set_index(cfg.date_col).sort_index()

    wide = wide.sort_index()
    wide = wide.replace([np.inf, -np.inf], np.nan).dropna(axis=1, thresh=cfg.min_history)
    wide = wide.ffill().dropna(axis=1)
    if cfg.max_symbols and wide.shape[1] > cfg.max_symbols:
        # Prefer keeping macro symbols if available, then most complete/variable columns.
        macro = [c for c in wide.columns if str(c).upper() in {m.upper() for m in cfg.macro_symbols}]
        rest = [c for c in wide.columns if c not in macro]
        vol_rank = wide[rest].pct_change().std().sort_values(ascending=False).index.tolist()
        keep = (macro + vol_rank)[:cfg.max_symbols]
        wide = wide[keep]
    return wide.reset_index().rename(columns={wide.index.name or cfg.date_col: "Date"})


def add_ratio_nodes(price_wide: pd.DataFrame) -> pd.DataFrame:
    out = price_wide.copy()
    cols_upper = {str(c).upper(): c for c in out.columns}
    pairs = [("HYG", "LQD", "HYG_LQD"), ("QQQ", "SPY", "QQQ_SPY"), ("IWM", "SPY", "IWM_SPY"), ("TLT", "SPY", "TLT_SPY")]
    for a, b, name in pairs:
        if a in cols_upper and b in cols_upper and name not in out.columns:
            out[name] = out[cols_upper[a]] / out[cols_upper[b]].replace(0, np.nan)
    return out


def identify_nodes(price_wide: pd.DataFrame, cfg: RunConfig) -> tuple[list[str], list[str]]:
    symbols = [c for c in price_wide.columns if c != "Date"]
    macro_set = {m.upper() for m in cfg.macro_symbols}
    macro_nodes = []
    for s in symbols:
        su = str(s).upper()
        if su in macro_set or su in {"HYG_LQD", "QQQ_SPY", "IWM_SPY", "TLT_SPY"}:
            macro_nodes.append(s)
    stock_nodes = [s for s in symbols if s not in macro_nodes]
    if not stock_nodes:
        # In degenerate data, trade all non-ratio nodes except obvious VIX-style items.
        stock_nodes = [s for s in symbols if "VIX" not in str(s).upper() and "_" not in str(s)]
        macro_nodes = [s for s in symbols if s not in stock_nodes]
    return stock_nodes, macro_nodes


def build_node_feature_panel(price_df: pd.DataFrame, cfg: RunConfig) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    wide = price_df.copy()
    if "Date" not in wide.columns:
        wide = wide.rename(columns={wide.columns[0]: "Date"})
    wide["Date"] = pd.to_datetime(wide["Date"])
    wide = wide.sort_values("Date").drop_duplicates("Date").set_index("Date")
    wide = wide.apply(pd.to_numeric, errors="coerce").ffill()
    wide = add_ratio_nodes(wide)
    stock_nodes, macro_nodes = identify_nodes(wide.reset_index(), cfg)
    all_nodes = stock_nodes + [m for m in macro_nodes if m not in stock_nodes]
    wide = wide[all_nodes].dropna(axis=0, how="all").ffill().dropna(axis=1)
    stock_nodes = [s for s in stock_nodes if s in wide.columns]
    macro_nodes = [s for s in macro_nodes if s in wide.columns]
    all_nodes = stock_nodes + [m for m in macro_nodes if m not in stock_nodes]
    wide = wide[all_nodes]

    logp = np.log(wide.replace(0, np.nan))
    ret1 = logp.diff()
    ret5 = logp.diff(5)
    ret21 = logp.diff(21)
    vol21 = rolling_vol(ret1, cfg.fast_lookback)
    vol63 = rolling_vol(ret1, cfg.corr_lookback)

    feat_mats = {
        "raw_ret_1d_z": rolling_z_by_col(ret1, cfg.lookback_z),
        "raw_ret_5d_z": rolling_z_by_col(ret5, cfg.lookback_z),
        "raw_mom_21d_z": rolling_z_by_col(ret21, cfg.lookback_z),
        "raw_vol_21d_z": rolling_z_by_col(vol21, cfg.lookback_z),
        "raw_vol_63d_z": rolling_z_by_col(vol63, cfg.lookback_z),
        "raw_price_z": rolling_z_by_col(logp, cfg.lookback_z),
    }

    rows = []
    dates = list(wide.index)
    for t_idx, dt in enumerate(dates):
        for item in all_nodes:
            row = {
                "date": dt,
                "t": int(t_idx),
                "item": str(item),
                "node_type": "macro" if item in macro_nodes else "stock",
                "is_macro": float(item in macro_nodes),
            }
            for name, mat in feat_mats.items():
                row[name] = float(mat.loc[dt, item]) if dt in mat.index and item in mat.columns else np.nan
            # Category signs for Prime-Fock.
            r = row.get("raw_ret_1d_z", np.nan)
            m = row.get("raw_mom_21d_z", np.nan)
            row["ret_sign"] = "pos" if pd.notna(r) and r >= 0 else "neg"
            row["mom_sign"] = "pos" if pd.notna(m) and m >= 0 else "neg"
            rows.append(row)
    panel = pd.DataFrame(rows)
    num_cols = [c for c in panel.columns if c.startswith("raw_")]
    panel[num_cols] = panel[num_cols].replace([np.inf, -np.inf], np.nan).fillna(0.0)
    panel = winsorize_df(panel, num_cols, cfg.winsor_z)
    return panel, wide.reset_index(), stock_nodes, macro_nodes


# ---------------------------------------------------------------------------
# Labels
# ---------------------------------------------------------------------------

def create_triple_barrier_labels(price_wide: pd.DataFrame, stock_nodes: Sequence[str], cfg: RunConfig) -> pd.DataFrame:
    df = price_wide.copy()
    df["Date"] = pd.to_datetime(df["Date"])
    df = df.sort_values("Date").set_index("Date")
    prices = df[list(stock_nodes)].astype(float)
    logp = np.log(prices.replace(0, np.nan))
    ret1 = logp.diff()
    vol = rolling_vol(ret1, cfg.vol_lookback).shift(0)

    market = None
    if cfg.relative_to_market:
        if cfg.market_symbol in df.columns:
            market = np.log(df[cfg.market_symbol].replace(0, np.nan)).diff()
        elif "SPY" in df.columns:
            market = np.log(df["SPY"].replace(0, np.nan)).diff()

    rows = []
    dates = list(prices.index)
    n = len(dates)
    for t in range(n):
        dt = dates[t]
        for item in stock_nodes:
            if not np.isfinite(prices.iloc[t][item]) or t + 1 >= n:
                continue
            v_daily = vol.iloc[t][item]
            if not np.isfinite(v_daily) or v_daily <= 1e-8:
                v_daily = np.nanmedian(vol[item].iloc[max(0, t-252):t+1].to_numpy(dtype=float))
            if not np.isfinite(v_daily) or v_daily <= 1e-8:
                v_daily = 0.02
            # Barriers are cumulative-log-return barriers scaled by horizon volatility.
            upper = cfg.take_profit_vol * v_daily * math.sqrt(cfg.label_horizon)
            lower = -cfg.stop_loss_vol * v_daily * math.sqrt(cfg.label_horizon)
            end = min(n - 1, t + cfg.label_horizon)
            exit_end = min(n - 1, t + cfg.exit_horizon)
            path = logp[item].iloc[t+1:end+1].to_numpy(dtype=float) - float(logp[item].iloc[t])
            if path.size == 0:
                continue
            hit_up = np.flatnonzero(path >= upper)
            hit_dn = np.flatnonzero(path <= lower)
            first_up = int(hit_up[0]) if hit_up.size else None
            first_dn = int(hit_dn[0]) if hit_dn.size else None
            if first_up is not None and (first_dn is None or first_up < first_dn):
                tb_label = 1
            elif first_dn is not None and (first_up is None or first_dn < first_up):
                tb_label = -1
            else:
                tb_label = 0
            fwd_ret = float(path[-1])
            rel_ret = fwd_ret
            if market is not None:
                mpath = market.iloc[t+1:end+1].sum()
                if np.isfinite(mpath):
                    rel_ret = fwd_ret - float(mpath)
            mfe = float(np.nanmax(path))
            mae = float(np.nanmin(path))
            hold_quality = mfe / (abs(mae) + 1e-6)
            # Good entry must hit upper first and not badly underperform market.
            entry = 1.0 if (tb_label == 1 and (not cfg.relative_to_market or rel_ret > 0)) else 0.0
            exit_path = logp[item].iloc[t+1:exit_end+1].to_numpy(dtype=float) - float(logp[item].iloc[t])
            exit_risk = 1.0 if exit_path.size and np.nanmin(exit_path) <= lower * math.sqrt(cfg.exit_horizon / cfg.label_horizon) else 0.0
            rows.append({
                "date": dt,
                "t": int(t),
                "item": str(item),
                "entry_label": entry,
                "exit_label": exit_risk,
                "triple_barrier_label": int(tb_label),
                "fwd_return": fwd_ret,
                "rel_fwd_return": rel_ret,
                "mfe": mfe,
                "mae": mae,
                "hold_quality": float(np.clip(hold_quality, 0, 25)),
                "upper_barrier": upper,
                "lower_barrier": lower,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Prime-Fock graph and feature extraction per fold
# ---------------------------------------------------------------------------

def make_pf_config(train_panel: pd.DataFrame, include_item_identity: bool) -> PrimeFockNDConfig:
    raw_cols = [c for c in train_panel.columns if c.startswith("raw_")]
    specs = [NumericBinningSpec(name=c, bins=quantile_bins(train_panel[c].to_numpy())) for c in raw_cols]
    specs += [
        CategoricalSpec(name="node_type", categories=["stock", "macro"], unknown_policy="zero"),
        CategoricalSpec(name="ret_sign", categories=["neg", "pos"], unknown_policy="zero"),
        CategoricalSpec(name="mom_sign", categories=["neg", "pos"], unknown_policy="zero"),
    ]
    if include_item_identity:
        specs.append(CategoricalSpec(name="item", categories=None, unknown_policy="zero"))
    return PrimeFockNDConfig(index_cols=["t", "item"], item_col="item", slice_cols=["t"], specs=specs)


def build_pf_graph_features(
    train_panel: pd.DataFrame,
    apply_panel: pd.DataFrame,
    cfg: RunConfig,
) -> tuple[pd.DataFrame, SharedSupportGraph, object, object]:
    pf_config = make_pf_config(train_panel, cfg.include_item_identity)
    encoder = PrimeFockNDEncoder(pf_config).fit(train_panel)
    encoded = encoder.transform(apply_panel)
    graph = ArithmeticGraphBuilder(kernel=cfg.graph_kernel).build(encoded)
    if cfg.use_stability:
        stability = StabilityFilter().fit_transform(graph)
        prop_graph = graph_with_adjacency(graph, stability.stable_adjacency, "stable")
    else:
        stability = dense_raw_stability_result(graph)
        prop_graph = graph
    pf = GenericFeatureEmitter().transform(encoded, graph, stability)
    pf = _normalize_slice_col(pf)
    pf = pf.rename(columns={c: f"pf_{c}" for c in pf.columns if c not in {"t", "item", "slice"}})
    # Return the raw graph for eigen/flux dynamics; stability is still emitted as PF features.
    # Stability filtering can legitimately prune broad market graphs to near-zero,
    # which is useful as a feature but bad as the only geometry for spectral/flux readout.
    return pf, graph, encoded, stability


def _normalize_slice_col(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "slice" in out.columns:
        out["t"] = out["slice"].apply(lambda x: int(x[0]) if isinstance(x, tuple) else int(x))
    return out


def raw_feature_rows(panel: pd.DataFrame) -> pd.DataFrame:
    keep = ["t", "date", "item", "node_type", "is_macro"] + [c for c in panel.columns if c.startswith("raw_")]
    out = panel[keep].copy()
    return out.rename(columns={c: c if c.startswith("raw_") else c for c in out.columns})


def emit_eigen_flux_features(
    graph: SharedSupportGraph,
    panel: pd.DataFrame,
    macro_nodes: Sequence[str],
    top_eigs: int,
) -> pd.DataFrame:
    # Map t to node signals.
    signal_col = "raw_ret_1d_z"
    mom_col = "raw_mom_21d_z"
    panel_key = panel[["t", "item", signal_col, mom_col, "is_macro"]].copy()
    item_to_macro = set(map(str, macro_nodes))
    rows = []
    A_all = np.asarray(graph.adjacency, dtype=float)
    items = list(map(str, graph.items))
    item_to_idx = {it: i for i, it in enumerate(items)}
    macro_idx = [item_to_idx[m] for m in item_to_macro if m in item_to_idx]

    prev_vecs = None
    for s, slice_key in enumerate(graph.slices):
        t = int(slice_key[0]) if isinstance(slice_key, tuple) else int(slice_key)
        A = 0.5 * (A_all[s] + A_all[s].T)
        np.fill_diagonal(A, 0.0)
        deg = A.sum(axis=1)
        An = normalized_adjacency(A)
        try:
            vals, vecs = np.linalg.eigh(An)
            order = np.argsort(vals)[::-1]
            vals = vals[order]
            vecs = vecs[:, order]
        except np.linalg.LinAlgError:
            vals = np.zeros(A.shape[0])
            vecs = np.eye(A.shape[0])
        k = min(top_eigs, vecs.shape[1])
        if prev_vecs is not None:
            # align signs to previous eigenvectors for smoother deltas
            for j in range(min(k, prev_vecs.shape[1])):
                if float(np.dot(vecs[:, j], prev_vecs[:, j])) < 0:
                    vecs[:, j] *= -1
        gap = float(vals[0] - vals[1]) if len(vals) > 1 else float(vals[0]) if len(vals) else 0.0
        eff_rank = float((np.sum(np.abs(vals)) ** 2) / (np.sum(vals ** 2) + 1e-12)) if len(vals) else 0.0
        concentration = float(np.sum(vecs[:, 0] ** 4)) if k > 0 else 0.0
        sub = panel_key[panel_key["t"] == t].set_index("item")
        signal = np.zeros(len(items))
        mom = np.zeros(len(items))
        is_macro = np.zeros(len(items))
        for i, it in enumerate(items):
            if it in sub.index:
                signal[i] = float(sub.loc[it, signal_col])
                mom[i] = float(sub.loc[it, mom_col])
                is_macro[i] = float(sub.loc[it, "is_macro"])
        diff = signal[None, :] - signal[:, None]
        flux_mat = A * diff
        net_flux = flux_mat.sum(axis=1)
        abs_flux = np.abs(flux_mat).sum(axis=1)
        pos_flux = np.maximum(flux_mat, 0).sum(axis=1)
        neg_flux = np.minimum(flux_mat, 0).sum(axis=1)
        graph_grad_mom = (A * (mom[None, :] - mom[:, None])).sum(axis=1)
        macro_pressure = np.zeros(len(items))
        macro_abs_pressure = np.zeros(len(items))
        if macro_idx:
            macro_pressure = (A[:, macro_idx] * signal[macro_idx][None, :]).sum(axis=1)
            macro_abs_pressure = np.abs(A[:, macro_idx] * signal[macro_idx][None, :]).sum(axis=1)
        for i, it in enumerate(items):
            row = {
                "t": t,
                "item": it,
                "eig_degree": float(deg[i]),
                "eig_lambda1": float(vals[0]) if len(vals) else 0.0,
                "eig_lambda2": float(vals[1]) if len(vals) > 1 else 0.0,
                "eig_spectral_gap": gap,
                "eig_effective_rank": eff_rank,
                "eig_top_concentration": concentration,
                "flux_net": float(net_flux[i]),
                "flux_abs": float(abs_flux[i]),
                "flux_pos": float(pos_flux[i]),
                "flux_neg": float(neg_flux[i]),
                "flux_graph_grad_mom": float(graph_grad_mom[i]),
                "macro_pressure": float(macro_pressure[i]),
                "macro_abs_pressure": float(macro_abs_pressure[i]),
            }
            for j in range(k):
                row[f"eig_load_{j+1}"] = float(vecs[i, j])
                if prev_vecs is not None and j < prev_vecs.shape[1]:
                    row[f"eig_delta_load_{j+1}"] = float(vecs[i, j] - prev_vecs[i, j])
                else:
                    row[f"eig_delta_load_{j+1}"] = 0.0
            # Distances to macro nodes in eigenspace.
            if macro_idx and k > 0:
                coords = vecs[:, :k]
                for macro_name in list(item_to_macro)[:8]:
                    if macro_name in item_to_idx:
                        mi = item_to_idx[macro_name]
                        row[f"eig_dist_to_{macro_name}"] = float(np.linalg.norm(coords[i] - coords[mi]))
            rows.append(row)
        prev_vecs = vecs[:, :k].copy()
    return pd.DataFrame(rows)


def emit_dynamics_features(graph: SharedSupportGraph, panel: pd.DataFrame, cfg: RunConfig) -> pd.DataFrame:
    if not cfg.dynamics_taus:
        return pd.DataFrame(columns=["t", "item"])
    rows = []
    item_list = list(map(str, graph.items))
    panel_key = panel[["t", "item", "raw_ret_1d_z"]].copy()
    try:
        emitter = DynamicsFeatureEmitter(
            diffusion=GraphDiffusionDynamics(generator="laplacian", normalize_output=True),
            wave=GraphWaveDynamics(generator="laplacian", gamma=cfg.wave_gamma),
        )
        for s, slice_key in enumerate(graph.slices):
            t = int(slice_key[0]) if isinstance(slice_key, tuple) else int(slice_key)
            sub = panel_key[panel_key["t"] == t].set_index("item")
            rho = np.ones(len(item_list), dtype=float)
            for i, it in enumerate(item_list):
                if it in sub.index:
                    rho[i] = abs(float(sub.loc[it, "raw_ret_1d_z"])) + 1e-3
            rho = rho / max(rho.sum(), 1e-12)
            dyn = emitter.transform(graph, s, rho, times=np.asarray(cfg.dynamics_taus, dtype=float))
            dyn = _normalize_slice_col(dyn)
            rows.append(dyn)
    except Exception as exc:
        # Dynamics is optional; return empty rather than killing the whole benchmark.
        print(f"[WARN] dynamics feature extraction failed: {exc}")
        return pd.DataFrame(columns=["t", "item"])
    if not rows:
        return pd.DataFrame(columns=["t", "item"])
    df = pd.concat(rows, ignore_index=True)
    # Pivot tau dimension into columns.
    val_cols = [c for c in df.columns if c not in {"slice", "t", "item", "generator", "wave_gamma", "t"} and pd.api.types.is_numeric_dtype(df[c])]
    # Careful: dynamics uses column named 't' for tau originally then overwritten by slice t. Preserve tau as dyn_tau if present.
    # The emitter's tau column is overwritten by _normalize_slice_col. Recreate by grouping order not possible here, so skip tau pivot
    # by aggregating mean/max over requested taus.
    agg = df.groupby(["t", "item"])[val_cols].agg(["mean", "max", "min"]) if val_cols else pd.DataFrame()
    if agg.empty:
        return pd.DataFrame(columns=["t", "item"])
    agg.columns = [f"dyn_{a}_{b}" for a, b in agg.columns]
    return agg.reset_index()


def build_fold_features(
    train_panel: pd.DataFrame,
    apply_panel: pd.DataFrame,
    macro_nodes: Sequence[str],
    cfg: RunConfig,
) -> pd.DataFrame:
    base = raw_feature_rows(apply_panel)
    pf, graph, encoded, stability = build_pf_graph_features(train_panel, apply_panel, cfg)
    eigflux = emit_eigen_flux_features(graph, apply_panel, macro_nodes, cfg.top_eigs)
    dyn = emit_dynamics_features(graph, apply_panel, cfg)

    out = base.merge(pf.drop(columns=["slice"], errors="ignore"), on=["t", "item"], how="left")
    out = out.merge(eigflux, on=["t", "item"], how="left")
    if not dyn.empty:
        out = out.merge(dyn, on=["t", "item"], how="left")
    num_cols = [c for c in out.columns if c not in {"date", "item", "node_type"} and pd.api.types.is_numeric_dtype(out[c])]
    out[num_cols] = out[num_cols].replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return out


# ---------------------------------------------------------------------------
# Walk-forward and NN
# ---------------------------------------------------------------------------

def purged_windows(date_values: Sequence[int], cfg: RunConfig) -> list[dict[str, int]]:
    dates = np.array(sorted(set(map(int, date_values))))
    n = len(dates)
    train_len = max(cfg.min_train_dates, int(n * cfg.train_frac))
    test_len = max(20, int(n * cfg.test_frac))
    embargo = cfg.embargo if cfg.embargo is not None else cfg.label_horizon
    windows = []
    if train_len + embargo + test_len >= n:
        train_len = max(cfg.min_train_dates, int(n * 0.55))
        test_len = max(20, int(n * 0.15))
    max_start = n - train_len - embargo - test_len
    if max_start < 0:
        raise ValueError(f"Not enough dates ({n}) for train/test/embargo. Reduce min_train_dates/n_windows.")
    starts = np.linspace(0, max_start, cfg.n_windows).round().astype(int) if cfg.n_windows > 1 else np.array([0])
    for wid, start in enumerate(starts):
        tr0 = int(dates[start])
        tr1 = int(dates[start + train_len - 1])
        te0 = int(dates[start + train_len + embargo])
        te1 = int(dates[start + train_len + embargo + test_len - 1])
        windows.append({"window": wid, "train_start_t": tr0, "train_end_t": tr1, "test_start_t": te0, "test_end_t": te1, "embargo": embargo})
    return windows


class MultiHeadMLP(nn.Module):
    def __init__(self, n_in: int, hidden: int = 256, depth: int = 3, dropout: float = 0.15):
        super().__init__()
        layers = [nn.LayerNorm(n_in), nn.Linear(n_in, hidden), nn.GELU(), nn.Dropout(dropout)]
        for _ in range(max(0, depth - 1)):
            layers += [nn.Linear(hidden, hidden), nn.GELU(), nn.Dropout(dropout)]
        self.body = nn.Sequential(*layers)
        self.entry = nn.Linear(hidden, 1)
        self.exit = nn.Linear(hidden, 1)
        self.hold = nn.Linear(hidden, 1)

    def forward(self, x):
        h = self.body(x)
        return {
            "entry_logit": self.entry(h).squeeze(-1),
            "exit_logit": self.exit(h).squeeze(-1),
            "hold": self.hold(h).squeeze(-1),
        }


def feature_blocks(columns: Sequence[str]) -> dict[str, list[str]]:
    cols = list(columns)
    raw = [c for c in cols if c.startswith("raw_") or c in {"is_macro"}]
    pf = [c for c in cols if c.startswith("pf_")]
    flux = [c for c in cols if c.startswith("flux_") or c.startswith("macro_")]
    eig = [c for c in cols if c.startswith("eig_")]
    dyn = [c for c in cols if c.startswith("dyn_")]
    return {
        "raw_only": raw,
        "raw_pf": raw + pf,
        "raw_flux_macro": raw + flux,
        "raw_eig": raw + eig,
        "raw_pf_flux_eig": raw + pf + flux + eig,
        "all": raw + pf + flux + eig + dyn,
        "no_pf": raw + flux + eig + dyn,
    }


def standardize_train_test(train: pd.DataFrame, test: pd.DataFrame, cols: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    Xtr = train[list(cols)].to_numpy(dtype=float)
    Xte = test[list(cols)].to_numpy(dtype=float)
    mu = np.nanmean(Xtr, axis=0)
    sig = np.nanstd(Xtr, axis=0)
    sig = np.where(sig > 1e-8, sig, 1.0)
    Xtr = np.nan_to_num((Xtr - mu) / sig, nan=0.0, posinf=0.0, neginf=0.0)
    Xte = np.nan_to_num((Xte - mu) / sig, nan=0.0, posinf=0.0, neginf=0.0)
    return Xtr.astype(np.float32), Xte.astype(np.float32)


def train_model(train_df: pd.DataFrame, test_df: pd.DataFrame, feature_cols: Sequence[str], cfg: RunConfig):
    if torch is None:
        raise RuntimeError(f"PyTorch is required for NN training. Import error: {_TORCH_IMPORT_ERROR}")
    device = cfg.device
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    Xtr, Xte = standardize_train_test(train_df, test_df, feature_cols)
    y_entry = train_df["entry_label"].to_numpy(dtype=np.float32)
    y_exit = train_df["exit_label"].to_numpy(dtype=np.float32)
    y_hold = train_df["hold_quality"].to_numpy(dtype=np.float32)

    # Chronological validation split within train.
    ts = train_df["t"].to_numpy()
    cut_t = np.quantile(ts, 0.85)
    tr_idx = np.flatnonzero(ts <= cut_t)
    va_idx = np.flatnonzero(ts > cut_t)
    if va_idx.size < 100:
        rng = np.random.default_rng(cfg.seed)
        idx = np.arange(len(train_df))
        rng.shuffle(idx)
        split = int(0.85 * len(idx))
        tr_idx, va_idx = idx[:split], idx[split:]

    model = MultiHeadMLP(Xtr.shape[1], cfg.hidden_dim, cfg.depth, cfg.dropout).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=cfg.lr, weight_decay=cfg.weight_decay)

    Xt = torch.tensor(Xtr, device=device)
    ye = torch.tensor(y_entry, device=device)
    yx = torch.tensor(y_exit, device=device)
    yh = torch.tensor(y_hold, device=device)

    # class weights for imbalance
    pos_e = max(float(y_entry.mean()), 1e-4)
    pos_x = max(float(y_exit.mean()), 1e-4)
    w_entry = torch.tensor((1 - pos_e) / pos_e, device=device, dtype=torch.float32)
    w_exit = torch.tensor((1 - pos_x) / pos_x, device=device, dtype=torch.float32)

    best_state = None
    best_val = float("inf")
    stale = 0
    rng = np.random.default_rng(cfg.seed)
    for epoch in range(cfg.epochs):
        model.train()
        order = tr_idx.copy()
        rng.shuffle(order)
        losses = []
        for start in range(0, len(order), cfg.batch_size):
            batch = order[start:start+cfg.batch_size]
            out = model(Xt[batch])
            loss_entry = F.binary_cross_entropy_with_logits(out["entry_logit"], ye[batch], pos_weight=w_entry)
            loss_exit = F.binary_cross_entropy_with_logits(out["exit_logit"], yx[batch], pos_weight=w_exit)
            loss_hold = F.smooth_l1_loss(out["hold"], yh[batch])
            loss = loss_entry + 0.7 * loss_exit + 0.1 * loss_hold
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
            losses.append(float(loss.detach().cpu()))
        model.eval()
        with torch.no_grad():
            outv = model(Xt[va_idx])
            val = (
                F.binary_cross_entropy_with_logits(outv["entry_logit"], ye[va_idx], pos_weight=w_entry)
                + 0.7 * F.binary_cross_entropy_with_logits(outv["exit_logit"], yx[va_idx], pos_weight=w_exit)
                + 0.1 * F.smooth_l1_loss(outv["hold"], yh[va_idx])
            )
            val_f = float(val.detach().cpu())
        if val_f < best_val - 1e-5:
            best_val = val_f
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
        if stale >= cfg.patience:
            break
    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        pred = model(torch.tensor(Xte, device=device))
        entry_prob = torch.sigmoid(pred["entry_logit"]).detach().cpu().numpy()
        exit_prob = torch.sigmoid(pred["exit_logit"]).detach().cpu().numpy()
        hold_pred = pred["hold"].detach().cpu().numpy()
    return entry_prob, exit_prob, hold_pred, {"best_val_loss": best_val, "epochs_ran": epoch + 1, "device": device}


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def evaluate_predictions(df: pd.DataFrame, cfg: RunConfig) -> dict[str, float]:
    out = {}
    n = len(df)
    if n == 0:
        return {"n": 0}
    out["n"] = float(n)
    out["entry_base_rate"] = float(df["entry_label"].mean())
    out["exit_base_rate"] = float(df["exit_label"].mean())
    out["entry_auc_proxy"] = rank_auc_proxy(df["entry_label"].to_numpy(), df["entry_prob"].to_numpy())
    k = max(1, int(math.ceil(cfg.top_frac * n)))
    top = df.sort_values("entry_prob", ascending=False).head(k)
    out["topk_n"] = float(k)
    out["topk_entry_precision"] = float(top["entry_label"].mean())
    out["topk_avg_fwd_return_bps"] = float(top["fwd_return"].mean() * 10000.0 - cfg.cost_bps)
    out["topk_avg_rel_return_bps"] = float(top["rel_fwd_return"].mean() * 10000.0 - cfg.cost_bps)
    out["topk_avg_mfe_bps"] = float(top["mfe"].mean() * 10000.0)
    out["topk_avg_mae_bps"] = float(top["mae"].mean() * 10000.0)
    out["topk_hold_quality"] = float(top["hold_quality"].mean())
    # Exit signal quality among highest exit-prob rows.
    ex = df.sort_values("exit_prob", ascending=False).head(k)
    out["topk_exit_precision"] = float(ex["exit_label"].mean())
    out["topk_exit_next_return_bps"] = float(ex["fwd_return"].mean() * 10000.0)
    return out


def rank_auc_proxy(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=float)
    score = np.asarray(score, dtype=float)
    pos = score[y > 0.5]
    neg = score[y <= 0.5]
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    # Mann-Whitney AUC via ranks.
    ranks = pd.Series(np.concatenate([pos, neg])).rank(method="average").to_numpy()
    rpos = ranks[:pos.size].sum()
    auc = (rpos - pos.size * (pos.size + 1) / 2) / (pos.size * neg.size)
    return float(auc)


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run(cfg: RunConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    seed_everything(cfg.seed)
    t0 = time.time()
    price_df = load_price_matrix(cfg)
    panel, price_wide, stock_nodes, macro_nodes = build_node_feature_panel(price_df, cfg)
    labels = create_triple_barrier_labels(price_wide, stock_nodes, cfg)
    valid_t = sorted(set(labels["t"].astype(int)))
    windows = purged_windows(valid_t, cfg)
    print(f"[DATA] dates={price_wide.shape[0]} stock_nodes={len(stock_nodes)} macro_nodes={len(macro_nodes)} labels={len(labels)}")
    print(f"[WINDOWS] {windows}")

    all_metrics = []
    all_preds = []
    for w in windows:
        wid = w["window"]
        train_mask = (labels["t"] >= w["train_start_t"]) & (labels["t"] <= w["train_end_t"])
        test_mask = (labels["t"] >= w["test_start_t"]) & (labels["t"] <= w["test_end_t"])
        train_labels = labels.loc[train_mask].copy()
        test_labels = labels.loc[test_mask].copy()
        if train_labels.empty or test_labels.empty:
            print(f"[SKIP] window={wid} empty labels")
            continue
        train_panel = panel[(panel["t"] >= w["train_start_t"]) & (panel["t"] <= w["train_end_t"])].copy()
        # Need apply_panel for train and test separately, fitted on train_panel.
        apply_train_panel = train_panel
        apply_test_panel = panel[(panel["t"] >= w["test_start_t"]) & (panel["t"] <= w["test_end_t"])].copy()
        print(f"[BUILD] window={wid} train_dates={train_panel['t'].nunique()} test_dates={apply_test_panel['t'].nunique()}")
        train_feat = build_fold_features(train_panel, apply_train_panel, macro_nodes, cfg)
        test_feat = build_fold_features(train_panel, apply_test_panel, macro_nodes, cfg)
        train_df = train_labels.merge(train_feat, on=["t", "item"], how="inner", suffixes=("", "_feat"))
        test_df = test_labels.merge(test_feat, on=["t", "item"], how="inner", suffixes=("", "_feat"))
        # Train only on stock rows; labels already stock-only but guard.
        train_df = train_df[train_df["node_type"].eq("stock")].copy()
        test_df = test_df[test_df["node_type"].eq("stock")].copy()
        if train_df.empty or test_df.empty:
            print(f"[SKIP] window={wid} empty merged features")
            continue
        numeric_cols = [
            c for c in train_df.columns
            if pd.api.types.is_numeric_dtype(train_df[c])
            and c not in {
                "entry_label", "exit_label", "triple_barrier_label", "fwd_return", "rel_fwd_return",
                "mfe", "mae", "hold_quality", "upper_barrier", "lower_barrier", "t"
            }
        ]
        blocks = feature_blocks(numeric_cols)
        for ablation, cols in blocks.items():
            cols = [c for c in cols if c in train_df.columns and c in test_df.columns]
            # Remove constants based on train.
            cols = [c for c in cols if float(np.nanstd(train_df[c].to_numpy(dtype=float))) > 1e-12]
            if len(cols) == 0:
                continue
            print(f"[TRAIN] window={wid} ablation={ablation} features={len(cols)} train={len(train_df)} test={len(test_df)}")
            entry_p, exit_p, hold_p, train_info = train_model(train_df, test_df, cols, cfg)
            date_col = "date" if "date" in test_df.columns else ("date_x" if "date_x" in test_df.columns else None)
            pred_cols = ([date_col] if date_col else []) + ["t", "item", "entry_label", "exit_label", "fwd_return", "rel_fwd_return", "mfe", "mae", "hold_quality"]
            pred_df = test_df[pred_cols].copy()
            if date_col and date_col != "date":
                pred_df = pred_df.rename(columns={date_col: "date"})
            pred_df["entry_prob"] = entry_p
            pred_df["exit_prob"] = exit_p
            pred_df["hold_pred"] = hold_p
            pred_df["window"] = wid
            pred_df["ablation"] = ablation
            all_preds.append(pred_df)
            metrics = evaluate_predictions(pred_df, cfg)
            metrics.update({
                "window": wid,
                "ablation": ablation,
                "n_features": len(cols),
                "train_n": len(train_df),
                **train_info,
            })
            all_metrics.append(metrics)
            print(
                f"[METRIC] w={wid} {ablation} precision@top={metrics.get('topk_entry_precision', np.nan):.3f} "
                f"base={metrics.get('entry_base_rate', np.nan):.3f} rel_bps={metrics.get('topk_avg_rel_return_bps', np.nan):.1f}"
            )
    metrics_df = pd.DataFrame(all_metrics)
    preds_df = pd.concat(all_preds, ignore_index=True) if all_preds else pd.DataFrame()
    if not metrics_df.empty:
        # Add aggregate rows.
        agg_cols = [c for c in metrics_df.columns if pd.api.types.is_numeric_dtype(metrics_df[c]) and c != "window"]
        agg = metrics_df.groupby("ablation")[agg_cols].mean(numeric_only=True).reset_index()
        agg["window"] = -1
        agg["row_type"] = "mean_by_ablation"
        metrics_df["row_type"] = "window"
        for c in metrics_df.columns:
            if c not in agg.columns:
                agg[c] = np.nan
        metrics_df = pd.concat([metrics_df, agg[metrics_df.columns]], ignore_index=True, sort=False)
    Path(cfg.out).parent.mkdir(parents=True, exist_ok=True) if Path(cfg.out).parent != Path(".") else None
    metrics_df.to_csv(cfg.out, index=False)
    preds_df.to_csv(cfg.preds_out, index=False)
    print(f"[SAVE] metrics={cfg.out} predictions={cfg.preds_out} seconds={time.time()-t0:.2f}")
    return metrics_df, preds_df


def parse_args() -> RunConfig:
    p = argparse.ArgumentParser(description="Prime-Fock macro/credit graph NN entry/exit benchmark")
    p.add_argument("--csv", type=str, default=None)
    p.add_argument("--macro-csv", type=str, default=None)
    p.add_argument("--date-col", type=str, default="Date")
    p.add_argument("--ticker-col", type=str, default=None)
    p.add_argument("--close-col", type=str, default="Close")
    p.add_argument("--max-symbols", type=int, default=80)
    p.add_argument("--min-history", type=int, default=260)
    p.add_argument("--macro-symbols", type=str, default=",".join(DEFAULT_MACRO_SYMBOLS))
    p.add_argument("--lookback-z", type=int, default=252)
    p.add_argument("--fast-lookback", type=int, default=21)
    p.add_argument("--corr-lookback", type=int, default=63)
    p.add_argument("--winsor-z", type=float, default=5.0)
    p.add_argument("--label-horizon", type=int, default=20)
    p.add_argument("--exit-horizon", type=int, default=10)
    p.add_argument("--take-profit-vol", type=float, default=1.5)
    p.add_argument("--stop-loss-vol", type=float, default=1.0)
    p.add_argument("--train-frac", type=float, default=0.60)
    p.add_argument("--test-frac", type=float, default=0.20)
    p.add_argument("--n-windows", type=int, default=3)
    p.add_argument("--embargo", type=int, default=None)
    p.add_argument("--min-train-dates", type=int, default=252)
    p.add_argument("--graph-kernel", type=str, default="log_gcd")
    p.add_argument("--include-item-identity", action="store_true")
    p.add_argument("--no-stability", action="store_true")
    p.add_argument("--dynamics-taus", type=str, default="")
    p.add_argument("--wave-gamma", type=float, default=0.05)
    p.add_argument("--top-eigs", type=int, default=6)
    p.add_argument("--hidden-dim", type=int, default=256)
    p.add_argument("--depth", type=int, default=3)
    p.add_argument("--dropout", type=float, default=0.15)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--weight-decay", type=float, default=1e-4)
    p.add_argument("--batch-size", type=int, default=4096)
    p.add_argument("--epochs", type=int, default=20)
    p.add_argument("--patience", type=int, default=5)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--device", type=str, default="auto")
    p.add_argument("--top-frac", type=float, default=0.10)
    p.add_argument("--cost-bps", type=float, default=5.0)
    p.add_argument("--out", type=str, default="macro_graph_entry_exit_nn_metrics.csv")
    p.add_argument("--preds-out", type=str, default="macro_graph_entry_exit_nn_predictions.csv")
    p.add_argument("--smoke", action="store_true")
    a = p.parse_args()
    return RunConfig(
        csv=a.csv,
        macro_csv=a.macro_csv,
        date_col=a.date_col,
        ticker_col=a.ticker_col,
        close_col=a.close_col,
        max_symbols=a.max_symbols,
        min_history=a.min_history,
        macro_symbols=tuple([x.strip() for x in a.macro_symbols.split(",") if x.strip()]),
        lookback_z=a.lookback_z,
        fast_lookback=a.fast_lookback,
        corr_lookback=a.corr_lookback,
        winsor_z=a.winsor_z,
        label_horizon=a.label_horizon,
        exit_horizon=a.exit_horizon,
        take_profit_vol=a.take_profit_vol,
        stop_loss_vol=a.stop_loss_vol,
        train_frac=a.train_frac,
        test_frac=a.test_frac,
        n_windows=a.n_windows,
        embargo=a.embargo,
        min_train_dates=a.min_train_dates,
        graph_kernel=a.graph_kernel,
        include_item_identity=a.include_item_identity,
        use_stability=not a.no_stability,
        dynamics_taus=tuple(float(x) for x in a.dynamics_taus.split(",") if x.strip()),
        wave_gamma=a.wave_gamma,
        top_eigs=a.top_eigs,
        hidden_dim=a.hidden_dim,
        depth=a.depth,
        dropout=a.dropout,
        lr=a.lr,
        weight_decay=a.weight_decay,
        batch_size=a.batch_size,
        epochs=a.epochs,
        patience=a.patience,
        seed=a.seed,
        device=a.device,
        top_frac=a.top_frac,
        cost_bps=a.cost_bps,
        out=a.out,
        preds_out=a.preds_out,
        smoke=a.smoke,
    )


if __name__ == "__main__":
    cfg = parse_args()
    print("[CONFIG]", json.dumps(cfg.__dict__, default=str))
    run(cfg)
