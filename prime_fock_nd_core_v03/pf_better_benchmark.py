"""
pf_better_benchmark.py
======================

Second-generation Prime-Fock benchmark for compression/dynamics, not generic
static tabular prediction.

What it adds over the first bench:
  - dataset caching: no repeated fetches after first run
  - graph caching: encoded tensors, dense arithmetic graphs, stability results
    are saved and reloaded by deterministic cache keys
  - persistence and mean baselines
  - walk-forward chronological validation
  - arithmetic-kernel ablations: log_gcd, log_lcm, arithmetic_distance,
    jaccard_gcd_lcm, coprime_indicator
  - item identity on/off ablation
  - stability on/off ablation
  - diffusion, unitary wave, and damped wave feature descriptors
  - unitary gamma=0 and damped gamma>0 tested with same graph operator

Expected project layout:

    project/
        prime_fock/
            __init__.py
            ... package files ...
        pf_better_benchmark.py

Run:

    python pf_better_benchmark.py --preset lite
    python pf_better_benchmark.py --preset markets
    python pf_better_benchmark.py --preset full

Outputs go under:

    pf_bench_cache/
        datasets/
        graphs/
        descriptors/
        results/

The graph cache is deliberately pickle-based because the library graph/stability
objects are dataclasses. This is for local experimental iteration, not a stable
cross-version artifact format.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import pickle
import zipfile
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Literal, Optional, Sequence

import numpy as np
import pandas as pd

from prime_fock_nd import (
    ArithmeticGraphBuilder,
    CategoricalSpec,
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


WeightScheme = Literal[
    "log_prime_weights",
    "uniform_weights",
    "random_monotone_weights",
    "shuffled_mode_weights",
]

KernelName = Literal[
    "log_gcd",
    "log_lcm",
    "arithmetic_distance",
    "jaccard_gcd_lcm",
    "coprime_indicator",
]


# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class BenchConfig:
    cache_dir: str = "pf_bench_cache"
    preset: str = "lite"
    datasets: tuple[str, ...] = ("plotly_stock_returns", "uci_household_power", "uci_appliances")
    k_values: tuple[int, ...] = (2, 4, 8)
    horizons: tuple[int, ...] = (1,)
    train_frac: float = 0.60
    test_frac: float = 0.20
    n_windows: int = 2
    ridge_alpha: float = 10.0
    seed: int = 0
    max_rows_stocks: int = 1000
    max_rows_power: int = 2000
    max_rows_appliances: int = 2500
    max_symbols_stocks: int = 40
    kernels: tuple[KernelName, ...] = ("log_gcd", "arithmetic_distance", "jaccard_gcd_lcm")
    weight_schemes: tuple[WeightScheme, ...] = (
        "log_prime_weights",
        "uniform_weights",
        "shuffled_mode_weights",
    )
    item_identity_options: tuple[bool, ...] = (False, True)
    stability_options: tuple[bool, ...] = (False, True)
    tau_values: tuple[float, ...] = (0.25, 1.0)
    damped_gammas: tuple[float, ...] = (0.05, 0.20)
    include_unitary: bool = True
    include_diffusion: bool = True
    use_phase_from_delta: bool = True
    graph_generator: str = "laplacian"


def config_for_preset(preset: str) -> BenchConfig:
    """
    if preset == "lite":
        return BenchConfig(
            preset="lite",
            datasets=("plotly_stock_returns", "uci_household_power", "uci_appliances"),
            k_values=(2, 4, 8),
            n_windows=2,
            max_rows_stocks=900,
            max_rows_power=1200,
            max_rows_appliances=1500,
            kernels=("log_gcd", "arithmetic_distance"),
            weight_schemes=("log_prime_weights", "uniform_weights", "shuffled_mode_weights"),
            item_identity_options=(False, True),
            stability_options=(False, True),
            tau_values=(0.25, 1.0),
            damped_gammas=(0.10,),
        )
    """
    if preset == "lite":
        return BenchConfig(
            preset="lite",
            datasets=("plotly_stock_returns",),
            k_values=(2, 4),
            n_windows=1,
            max_rows_stocks=500,
            kernels=("log_gcd",),
            weight_schemes=("log_prime_weights", "shuffled_mode_weights"),
            item_identity_options=(False,),
            stability_options=(False,),
            tau_values=(0.25,),
            damped_gammas=(0.10,),
            include_unitary=False,
        )
    if preset == "markets":
        return BenchConfig(
            preset="markets",
            datasets=("plotly_stock_returns",),
            k_values=(2, 4, 8, 12, 16),
            n_windows=3,
            max_rows_stocks=1200,
            kernels=("log_gcd", "log_lcm", "arithmetic_distance", "jaccard_gcd_lcm"),
            weight_schemes=("log_prime_weights", "uniform_weights", "random_monotone_weights", "shuffled_mode_weights"),
            item_identity_options=(False, True),
            stability_options=(False, True),
            tau_values=(0.10, 0.25, 1.0),
            damped_gammas=(0.05, 0.20),
        )
    if preset == "full":
        return BenchConfig(
            preset="full",
            datasets=("plotly_stock_returns", "uci_household_power", "uci_appliances"),
            k_values=(2, 4, 8, 12, 16),
            n_windows=3,
            max_rows_stocks=1200,
            max_rows_power=3000,
            max_rows_appliances=4000,
            kernels=("log_gcd", "log_lcm", "arithmetic_distance", "jaccard_gcd_lcm", "coprime_indicator"),
            weight_schemes=("log_prime_weights", "uniform_weights", "random_monotone_weights", "shuffled_mode_weights"),
            item_identity_options=(False, True),
            stability_options=(False, True),
            tau_values=(0.10, 0.25, 0.50, 1.0),
            damped_gammas=(0.05, 0.20, 0.50),
        )
    raise ValueError("preset must be one of: lite, markets, full")


# -----------------------------------------------------------------------------
# Cache helpers
# -----------------------------------------------------------------------------

def ensure_dirs(cfg: BenchConfig) -> dict[str, Path]:
    root = Path(cfg.cache_dir)
    dirs = {
        "root": root,
        "datasets": root / "datasets",
        "graphs": root / "graphs",
        "descriptors": root / "descriptors",
        "results": root / "results",
    }
    for p in dirs.values():
        p.mkdir(parents=True, exist_ok=True)
    return dirs


def sha1_bytes(x: bytes) -> str:
    return hashlib.sha1(x).hexdigest()[:16]


def stable_json_hash(obj: Any) -> str:
    data = json.dumps(obj, sort_keys=True, default=str).encode("utf-8")
    return sha1_bytes(data)


def matrix_hash(X: np.ndarray) -> str:
    X = np.asarray(X, dtype=np.float64, order="C")
    h = hashlib.sha1()
    h.update(str(X.shape).encode())
    h.update(X.tobytes())
    return h.hexdigest()[:16]


def load_pickle(path: Path):
    with path.open("rb") as f:
        return pickle.load(f)


def save_pickle(path: Path, obj: Any) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("wb") as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    tmp.replace(path)


# -----------------------------------------------------------------------------
# Dataset loaders with cache
# -----------------------------------------------------------------------------

def cached_dataset(name: str, cfg: BenchConfig, dirs: dict[str, Path]) -> pd.DataFrame:
    path = dirs["datasets"] / f"{name}_{cfg.preset}.pkl"
    if path.exists():
        print(f"[CACHE] dataset {name}: {path}")
        return load_pickle(path)

    if name == "plotly_stock_returns":
        df = load_plotly_stock_panel(max_symbols=cfg.max_symbols_stocks).tail(cfg.max_rows_stocks)
    elif name == "uci_household_power":
        df = load_uci_household_power(resample_rule="1h", max_rows=cfg.max_rows_power)
    elif name == "uci_appliances":
        df = load_uci_appliances_energy(max_rows=cfg.max_rows_appliances)
    else:
        raise ValueError(f"unknown dataset: {name}")

    save_pickle(path, df)
    return df


def load_plotly_stock_panel(max_symbols: int = 40) -> pd.DataFrame:
    url = "https://raw.githubusercontent.com/plotly/datasets/master/all_stocks_5yr.csv"
    print("[LOAD] Plotly all_stocks_5yr.csv")
    df = pd.read_csv(url)
    required = {"date", "close", "Name"}
    if not required.issubset(df.columns):
        raise RuntimeError(f"Unexpected Plotly stock columns: {df.columns.tolist()}")
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df["close"] = pd.to_numeric(df["close"], errors="coerce")
    df = df.dropna(subset=["date", "close", "Name"])
    counts = df.groupby("Name")["date"].nunique().sort_values(ascending=False)
    symbols = list(counts.head(max_symbols).index)
    panel = (
        df[df["Name"].isin(symbols)]
        .pivot_table(index="date", columns="Name", values="close", aggfunc="last")
        .sort_index()
        .dropna(axis=0, how="any")
    )
    panel = np.log(panel).diff().dropna()
    panel = panel.replace([np.inf, -np.inf], np.nan).dropna(axis=1, how="any")
    print("[DATA] Plotly stock returns:", panel.shape)
    return panel


def load_uci_household_power(resample_rule: str = "1h", max_rows: Optional[int] = None) -> pd.DataFrame:
    url = "https://archive.ics.uci.edu/static/public/235/individual+household+electric+power+consumption.zip"
    print("[LOAD] UCI household power")
    raw = urllib.request.urlopen(url).read()
    zf = zipfile.ZipFile(io.BytesIO(raw))
    name = [n for n in zf.namelist() if n.endswith(".txt")][0]
    df = pd.read_csv(zf.open(name), sep=";", na_values=["?"], low_memory=False)
    df["datetime"] = pd.to_datetime(df["Date"] + " " + df["Time"], dayfirst=True, errors="coerce")
    df = df.drop(columns=["Date", "Time"]).set_index("datetime")
    df = df.apply(pd.to_numeric, errors="coerce").dropna()
    df = df.resample(resample_rule).mean().dropna()
    cols = [
        "Global_active_power",
        "Global_reactive_power",
        "Voltage",
        "Global_intensity",
        "Sub_metering_1",
        "Sub_metering_2",
        "Sub_metering_3",
    ]
    df = df[cols].dropna()
    if max_rows is not None:
        df = df.iloc[:max_rows]
    print("[DATA] Household power:", df.shape)
    return df


def load_uci_appliances_energy(max_rows: Optional[int] = None) -> pd.DataFrame:
    url = "https://archive.ics.uci.edu/static/public/374/appliances+energy+prediction.zip"
    print("[LOAD] UCI appliances energy")
    raw = urllib.request.urlopen(url).read()
    zf = zipfile.ZipFile(io.BytesIO(raw))
    csv_name = [n for n in zf.namelist() if n.endswith(".csv")][0]
    df = pd.read_csv(zf.open(csv_name))
    df["date"] = pd.to_datetime(df["date"])
    df = df.set_index("date")
    df = df.drop(columns=[c for c in ["rv1", "rv2"] if c in df.columns])
    df = df.apply(pd.to_numeric, errors="coerce").dropna()
    if max_rows is not None:
        df = df.iloc[:max_rows]
    print("[DATA] Appliances:", df.shape)
    return df


# -----------------------------------------------------------------------------
# Numerical baselines
# -----------------------------------------------------------------------------

def standardize_train_test(X_train: np.ndarray, X_test: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mu = np.nanmean(X_train, axis=0)
    sd = np.nanstd(X_train, axis=0)
    sd = np.where(sd <= 1e-12, 1.0, sd)
    Xtr = np.nan_to_num((X_train - mu) / sd, nan=0.0, posinf=0.0, neginf=0.0)
    Xte = np.nan_to_num((X_test - mu) / sd, nan=0.0, posinf=0.0, neginf=0.0)
    return Xtr, Xte, mu, sd


def ridge_fit_predict(Z_train: np.ndarray, Y_train: np.ndarray, Z_test: np.ndarray, alpha: float) -> np.ndarray:
    Z_train = np.asarray(Z_train, dtype=float)
    Z_test = np.asarray(Z_test, dtype=float)
    Y_train = np.asarray(Y_train, dtype=float)
    z_mu = Z_train.mean(axis=0)
    z_sd = Z_train.std(axis=0)
    z_sd = np.where(z_sd <= 1e-12, 1.0, z_sd)
    Ztr = (Z_train - z_mu) / z_sd
    Zte = (Z_test - z_mu) / z_sd
    Ztr = np.hstack([np.ones((Ztr.shape[0], 1)), Ztr])
    Zte = np.hstack([np.ones((Zte.shape[0], 1)), Zte])
    penalty = np.eye(Ztr.shape[1]) * alpha
    penalty[0, 0] = 0.0
    beta = np.linalg.solve(Ztr.T @ Ztr + penalty, Ztr.T @ Y_train)
    return Zte @ beta


def svd_project_train_test(X_train: np.ndarray, X_test: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    mu = X_train.mean(axis=0)
    Xc_train = X_train - mu
    Xc_test = X_test - mu
    _, _, vt = np.linalg.svd(Xc_train, full_matrices=False)
    k_eff = min(k, vt.shape[0])
    B = vt[:k_eff].T
    return Xc_train @ B, Xc_test @ B


def random_project_train_test(X_train: np.ndarray, X_test: np.ndarray, k: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    R = rng.normal(0.0, 1.0 / np.sqrt(k), size=(X_train.shape[1], k))
    return X_train @ R, X_test @ R


def regression_metrics_matrix(Y_true: np.ndarray, Y_pred: np.ndarray) -> dict[str, float]:
    Y_true = np.asarray(Y_true, dtype=float)
    Y_pred = np.asarray(Y_pred, dtype=float)
    err = Y_true - Y_pred
    denom = float(np.sum((Y_true - Y_true.mean(axis=0, keepdims=True)) ** 2))
    r2 = np.nan if denom <= 1e-12 else 1.0 - float(np.sum(err ** 2)) / denom
    corrs = []
    for j in range(Y_true.shape[1]):
        a = Y_true[:, j]
        b = Y_pred[:, j]
        if a.std() > 1e-12 and b.std() > 1e-12:
            corrs.append(float(np.corrcoef(a, b)[0, 1]))
    return {
        "r2": float(r2),
        "rmse": float(np.sqrt(np.mean(err ** 2))),
        "mae": float(np.mean(np.abs(err))),
        "mean_channel_corr": float(np.mean(corrs)) if corrs else np.nan,
        "median_channel_corr": float(np.median(corrs)) if corrs else np.nan,
        "pred_std": float(np.std(Y_pred)),
        "target_std": float(np.std(Y_true)),
    }


def corr_matrix_error(X_true: np.ndarray, X_pred: np.ndarray) -> float:
    if X_true.shape[1] < 2:
        return np.nan
    C0 = np.corrcoef(X_true, rowvar=False)
    C1 = np.corrcoef(X_pred, rowvar=False)
    C0 = np.nan_to_num(C0, nan=0.0, posinf=0.0, neginf=0.0)
    C1 = np.nan_to_num(C1, nan=0.0, posinf=0.0, neginf=0.0)
    return float(np.linalg.norm(C0 - C1, ord="fro") / C0.shape[0])


def evaluate_descriptor(dataset: str, window_id: int, horizon: int, method: str, k: int,
                        Z_train: np.ndarray, Z_test: np.ndarray,
                        X_train: np.ndarray, X_test: np.ndarray,
                        Y_train: np.ndarray, Y_test: np.ndarray,
                        cfg: BenchConfig, extra: Optional[dict[str, Any]] = None) -> dict[str, Any]:
    Y_pred = ridge_fit_predict(Z_train, Y_train, Z_test, alpha=cfg.ridge_alpha)
    X_rec = ridge_fit_predict(Z_train, X_train, Z_test, alpha=cfg.ridge_alpha)
    mf = regression_metrics_matrix(Y_test, Y_pred)
    mr = regression_metrics_matrix(X_test, X_rec)
    row = {
        "dataset": dataset,
        "window": window_id,
        "horizon": horizon,
        "method": method,
        "k": int(k),
        "forecast_r2": mf["r2"],
        "forecast_rmse": mf["rmse"],
        "forecast_mae": mf["mae"],
        "forecast_corr": mf["mean_channel_corr"],
        "reconstruction_r2": mr["r2"],
        "reconstruction_rmse": mr["rmse"],
        "corr_matrix_error": corr_matrix_error(X_test, X_rec),
        "compression_ratio": float(k / X_train.shape[1]) if X_train.shape[1] else np.nan,
    }
    if extra:
        row.update(extra)
    return row


# -----------------------------------------------------------------------------
# Windows and supervised pairs
# -----------------------------------------------------------------------------

def supervised_pairs(X: np.ndarray, horizon: int) -> tuple[np.ndarray, np.ndarray]:
    return X[:-horizon], X[horizon:]


def walk_forward_windows(n: int, train_frac: float, test_frac: float, n_windows: int) -> list[tuple[np.ndarray, np.ndarray]]:
    test_len = max(10, int(round(n * test_frac)))
    train_len0 = max(20, int(round(n * train_frac)))
    windows = []
    max_start = n - test_len
    if max_start <= train_len0:
        return [(np.arange(0, train_len0), np.arange(train_len0, min(n, train_len0 + test_len)))]
    if n_windows == 1:
        train_ends = [train_len0]
    else:
        train_ends = np.linspace(train_len0, max_start, n_windows).astype(int).tolist()
    for end in train_ends:
        train_idx = np.arange(0, end)
        test_idx = np.arange(end, min(n, end + test_len))
        if len(test_idx) >= 5:
            windows.append((train_idx, test_idx))
    return windows


# -----------------------------------------------------------------------------
# Prime-Fock encoding/graph caching
# -----------------------------------------------------------------------------

def quantile_bins(values: np.ndarray, qs=(0.05, 0.15, 0.30, 0.50, 0.70, 0.85, 0.95)) -> list[float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return [-1.0, 0.0, 1.0]
    out = sorted(set(float(x) for x in np.quantile(vals, qs)))
    if len(out) < 2:
        m = float(np.mean(vals))
        s = float(np.std(vals)) or 1.0
        out = [m - s, m, m + s]
    return out


def matrix_to_long_panel(X: np.ndarray, item_names: Sequence[str], start_t: int = 0) -> pd.DataFrame:
    X = np.asarray(X, dtype=float)
    rows = []
    prev = np.zeros(X.shape[1], dtype=float)
    for local_t in range(X.shape[0]):
        t = start_t + local_t
        for j, item in enumerate(item_names):
            value = float(X[local_t, j])
            delta = float(value - prev[j]) if local_t > 0 else 0.0
            rows.append({
                "t": int(t),
                "item": str(item),
                "value": value,
                "delta": delta,
                "abs_delta": abs(delta),
                "sign": "pos" if value >= 0 else "neg",
                "delta_sign": "pos" if delta >= 0 else "neg",
            })
        prev = X[local_t].copy()
    return pd.DataFrame(rows)


def make_pf_config(train_long: pd.DataFrame, include_item_identity: bool) -> PrimeFockNDConfig:
    specs = [
        NumericBinningSpec(name="value", bins=quantile_bins(train_long["value"].to_numpy())),
        NumericBinningSpec(name="delta", bins=quantile_bins(train_long["delta"].to_numpy())),
        NumericBinningSpec(name="abs_delta", bins=quantile_bins(train_long["abs_delta"].to_numpy())),
        CategoricalSpec(name="sign", categories=None, unknown_policy="zero"),
        CategoricalSpec(name="delta_sign", categories=None, unknown_policy="zero"),
    ]
    if include_item_identity:
        specs.append(CategoricalSpec(name="item", categories=None, unknown_policy="zero"))
    return PrimeFockNDConfig(index_cols=["t", "item"], item_col="item", slice_cols=["t"], specs=specs)


def make_mode_weights(scheme: WeightScheme, encoded_train, seed: int) -> Optional[np.ndarray]:
    rng = np.random.default_rng(seed)
    P = encoded_train.exponents.shape[2]
    logw = np.log(encoded_train.primes.astype(float))
    if scheme == "log_prime_weights":
        return None
    if scheme == "uniform_weights":
        return np.ones(P) * float(np.mean(logw))
    if scheme == "random_monotone_weights":
        raw = np.cumsum(rng.gamma(shape=2.0, scale=1.0, size=P))
        return raw * (float(np.mean(logw)) / float(np.mean(raw)))
    if scheme == "shuffled_mode_weights":
        w = logw.copy()
        rng.shuffle(w)
        return w
    raise ValueError(scheme)


def raw_stability_result(graph: SharedSupportGraph) -> StabilityResult:
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


def graph_with_adjacency(base: SharedSupportGraph, adjacency: np.ndarray, kernel_suffix: str) -> SharedSupportGraph:
    return SharedSupportGraph(
        adjacency=np.asarray(adjacency, dtype=float),
        active_modes=base.active_modes,
        slices=base.slices,
        items=base.items,
        kernel=f"{base.kernel}_{kernel_suffix}",
    )


def build_or_load_pf_bundle(
    X_train: np.ndarray,
    X_test: np.ndarray,
    item_names: Sequence[str],
    dataset: str,
    window_id: int,
    horizon: int,
    kernel: KernelName,
    weight_scheme: WeightScheme,
    include_item_identity: bool,
    use_stability: bool,
    cfg: BenchConfig,
    dirs: dict[str, Path],
) -> dict[str, Any]:
    key_obj = {
        "dataset": dataset,
        "window": window_id,
        "horizon": horizon,
        "kernel": kernel,
        "weight_scheme": weight_scheme,
        "include_item_identity": include_item_identity,
        "use_stability": use_stability,
        "train_hash": matrix_hash(X_train),
        "test_hash": matrix_hash(X_test),
        "items": list(map(str, item_names)),
    }
    key = stable_json_hash(key_obj)
    path = dirs["graphs"] / f"graph_bundle_{key}.pkl"
    if path.exists():
        print(f"[CACHE] graph bundle {key}")
        return load_pickle(path)

    print(f"[BUILD] graph bundle {key_obj}")
    train_long = matrix_to_long_panel(X_train, item_names=item_names, start_t=0)
    test_long = matrix_to_long_panel(X_test, item_names=item_names, start_t=len(X_train))
    pf_config = make_pf_config(train_long, include_item_identity=include_item_identity)
    encoder = PrimeFockNDEncoder(pf_config).fit(train_long)
    encoded_train = encoder.transform(train_long)
    encoded_test = encoder.transform(test_long)
    mode_weights = make_mode_weights(weight_scheme, encoded_train, cfg.seed + window_id)

    builder = ArithmeticGraphBuilder(mode_weights=mode_weights, kernel=kernel)
    graph_train = builder.build(encoded_train)
    graph_test = builder.build(encoded_test)

    if use_stability:
        stab_train = StabilityFilter().fit_transform(graph_train)
        stab_test = StabilityFilter().fit_transform(graph_test)
        prop_graph_train = graph_with_adjacency(graph_train, stab_train.stable_adjacency, "stable")
        prop_graph_test = graph_with_adjacency(graph_test, stab_test.stable_adjacency, "stable")
    else:
        stab_train = raw_stability_result(graph_train)
        stab_test = raw_stability_result(graph_test)
        prop_graph_train = graph_train
        prop_graph_test = graph_test

    bundle = {
        "key": key,
        "key_obj": key_obj,
        "encoder": encoder,
        "encoded_train": encoded_train,
        "encoded_test": encoded_test,
        "graph_train": graph_train,
        "graph_test": graph_test,
        "prop_graph_train": prop_graph_train,
        "prop_graph_test": prop_graph_test,
        "stability_train": stab_train,
        "stability_test": stab_test,
        "mode_weights": mode_weights,
    }
    save_pickle(path, bundle)
    return bundle


# -----------------------------------------------------------------------------
# Prime-Fock descriptors
# -----------------------------------------------------------------------------

def pf_features_to_state_matrix(features: pd.DataFrame, item_names: Sequence[str]) -> np.ndarray:
    numeric_cols = [
        c for c in features.columns
        if c not in {"slice", "item"} and pd.api.types.is_numeric_dtype(features[c])
    ]
    pieces = []
    for col in numeric_cols:
        piv = features.pivot_table(index="slice", columns="item", values=col, aggfunc="mean")
        for item in item_names:
            if item not in piv.columns:
                piv[item] = 0.0
        piv = piv[list(item_names)]
        piv.columns = [f"{col}__{item}" for item in item_names]
        pieces.append(piv)
    if not pieces:
        raise RuntimeError("GenericFeatureEmitter produced no numeric features")
    out = pd.concat(pieces, axis=1).sort_index()
    return out.replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy(dtype=float)


def build_static_descriptor(bundle: dict[str, Any], item_names: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    emitter = GenericFeatureEmitter()
    f_train = emitter.transform(bundle["encoded_train"], bundle["graph_train"], bundle["stability_train"])
    f_test = emitter.transform(bundle["encoded_test"], bundle["graph_test"], bundle["stability_test"])
    return pf_features_to_state_matrix(f_train, item_names), pf_features_to_state_matrix(f_test, item_names)


def density_and_phase_from_state(X: np.ndarray, row_idx: int, use_phase: bool) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(X[row_idx], dtype=float)
    rho = np.abs(x) + 1e-12
    rho = rho / rho.sum()
    if use_phase and row_idx > 0:
        dx = x - np.asarray(X[row_idx - 1], dtype=float)
        theta = np.where(dx >= 0, 0.0, np.pi)
    else:
        theta = np.zeros_like(x)
    psi = np.sqrt(rho).astype(complex) * np.exp(1j * theta)
    return rho, psi


def dynamics_state_matrix(
    graph: SharedSupportGraph,
    X: np.ndarray,
    tau_values: Sequence[float],
    mode: Literal["diffusion", "unitary", "damped", "all"],
    gamma: float,
    cfg: BenchConfig,
) -> np.ndarray:
    """
    State-level dynamics descriptor. For each slice/time row, concatenate
    per-item propagated densities/residuals plus global interference summaries.
    """
    X = np.asarray(X, dtype=float)
    diff = GraphDiffusionDynamics(generator=cfg.graph_generator, normalize_output=True)
    wave = GraphWaveDynamics(generator=cfg.graph_generator, gamma=gamma)
    rows = []
    for s in range(X.shape[0]):
        rho0, psi0 = density_and_phase_from_state(X, s, cfg.use_phase_from_delta)
        op_d = diff.prepare(graph, s)
        op_w = wave.prepare(graph, s)
        feat = []
        for tau in tau_values:
            rho_d = diff.evolve(op_d, rho0, float(tau))
            if mode in {"diffusion", "all"}:
                feat.extend(rho_d.tolist())
                p = np.clip(rho_d, 1e-12, None)
                feat.append(float(-(p * np.log(p)).sum()))
                feat.append(float((rho_d ** 2).sum()))
            if mode in {"unitary", "damped", "all"}:
                rho_w = wave.probability(op_w, psi0, float(tau))
                residual = rho_w - rho_d
                feat.extend(rho_w.tolist())
                feat.extend(residual.tolist())
                feat.append(float(np.sum(rho_w)))              # wave norm
                feat.append(float(np.sum(np.abs(residual))))   # total residual/interference
                feat.append(float(np.max(np.abs(residual))))   # peak residual
                try:
                    M = wave.interference_matrix(op_w, float(tau))
                    feat.append(float(np.linalg.norm(M, ord="fro")))
                    feat.append(float(np.mean(M)))
                    feat.append(float(np.max(M)))
                except Exception:
                    feat.extend([0.0, 0.0, 0.0])
        rows.append(feat)
    return np.nan_to_num(np.asarray(rows, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)


def cached_descriptor(
    name: str,
    bundle: dict[str, Any],
    X_train: np.ndarray,
    X_test: np.ndarray,
    item_names: Sequence[str],
    cfg: BenchConfig,
    dirs: dict[str, Path],
    builder_fn,
) -> tuple[np.ndarray, np.ndarray]:
    key = stable_json_hash({
        "name": name,
        "graph_key": bundle["key"],
        "train_hash": matrix_hash(X_train),
        "test_hash": matrix_hash(X_test),
        "tau": cfg.tau_values,
        "phase": cfg.use_phase_from_delta,
        "generator": cfg.graph_generator,
    })
    path = dirs["descriptors"] / f"descriptor_{key}.npz"
    if path.exists():
        z = np.load(path)
        return z["train"], z["test"]
    F_train, F_test = builder_fn()
    np.savez_compressed(path, train=F_train, test=F_test)
    return F_train, F_test


# -----------------------------------------------------------------------------
# Main benchmark
# -----------------------------------------------------------------------------

def run_dataset(dataset_name: str, X_df: pd.DataFrame, cfg: BenchConfig, dirs: dict[str, Path]) -> pd.DataFrame:
    X_raw = X_df.to_numpy(dtype=float)
    item_names = list(map(str, X_df.columns))
    rows: list[dict[str, Any]] = []

    for horizon in cfg.horizons:
        X0_raw, Y_raw = supervised_pairs(X_raw, horizon)
        windows = walk_forward_windows(len(X0_raw), cfg.train_frac, cfg.test_frac, cfg.n_windows)
        print(f"[WINDOWS] {dataset_name} horizon={horizon}: {len(windows)} windows")

        for w_id, (train_idx, test_idx) in enumerate(windows):
            X_train_raw = X0_raw[train_idx]
            X_test_raw = X0_raw[test_idx]
            Y_train_raw = Y_raw[train_idx]
            Y_test_raw = Y_raw[test_idx]

            X_train, X_test, mu, sd = standardize_train_test(X_train_raw, X_test_raw)
            Y_train = np.nan_to_num((Y_train_raw - mu) / sd, nan=0.0, posinf=0.0, neginf=0.0)
            Y_test = np.nan_to_num((Y_test_raw - mu) / sd, nan=0.0, posinf=0.0, neginf=0.0)

            # Mean baseline.
            Y_mean = np.broadcast_to(Y_train.mean(axis=0, keepdims=True), Y_test.shape)
            m = regression_metrics_matrix(Y_test, Y_mean)
            rows.append({
                "dataset": dataset_name, "window": w_id, "horizon": horizon,
                "method": "mean", "k": 0,
                "forecast_r2": m["r2"], "forecast_rmse": m["rmse"], "forecast_mae": m["mae"],
                "forecast_corr": m["mean_channel_corr"],
                "reconstruction_r2": np.nan, "reconstruction_rmse": np.nan,
                "corr_matrix_error": np.nan, "compression_ratio": 0.0,
            })

            # Persistence baseline: X[t+h] = X[t].
            m = regression_metrics_matrix(Y_test, X_test)
            rows.append({
                "dataset": dataset_name, "window": w_id, "horizon": horizon,
                "method": "persistence", "k": X_train.shape[1],
                "forecast_r2": m["r2"], "forecast_rmse": m["rmse"], "forecast_mae": m["mae"],
                "forecast_corr": m["mean_channel_corr"],
                "reconstruction_r2": 1.0, "reconstruction_rmse": 0.0,
                "corr_matrix_error": 0.0, "compression_ratio": 1.0,
            })

            # Raw lag ridge baseline.
            rows.append(evaluate_descriptor(dataset_name, w_id, horizon, "raw_lag_ridge", X_train.shape[1],
                                            X_train, X_test, X_train, X_test, Y_train, Y_test, cfg))

            for k in cfg.k_values:
                if k <= 0:
                    continue
                # PCA / random projection.
                Ztr, Zte = svd_project_train_test(X_train, X_test, k=k)
                rows.append(evaluate_descriptor(dataset_name, w_id, horizon, "pca", k, Ztr, Zte,
                                                X_train, X_test, Y_train, Y_test, cfg))
                Ztr, Zte = random_project_train_test(X_train, X_test, k=k, seed=cfg.seed + k + w_id)
                rows.append(evaluate_descriptor(dataset_name, w_id, horizon, "random_projection", k, Ztr, Zte,
                                                X_train, X_test, Y_train, Y_test, cfg))

            # Prime-Fock ablations.
            for kernel in cfg.kernels:
                for weight_scheme in cfg.weight_schemes:
                    for include_item in cfg.item_identity_options:
                        for use_stability in cfg.stability_options:
                            try:
                                bundle = build_or_load_pf_bundle(
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
                                extra = {
                                    "kernel": kernel,
                                    "weight_scheme": weight_scheme,
                                    "item_identity": include_item,
                                    "stability": use_stability,
                                }

                                F_static_tr, F_static_te = cached_descriptor(
                                    "static",
                                    bundle,
                                    X_train,
                                    X_test,
                                    item_names,
                                    cfg,
                                    dirs,
                                    lambda: build_static_descriptor(bundle, item_names),
                                )

                                descriptor_specs: list[tuple[str, np.ndarray, np.ndarray, dict[str, Any]]] = []
                                descriptor_specs.append(("pf_static", F_static_tr, F_static_te, {"dynamics": "static"}))

                                if cfg.include_diffusion:
                                    Dtr, Dte = cached_descriptor(
                                        "diffusion",
                                        bundle,
                                        X_train,
                                        X_test,
                                        item_names,
                                        cfg,
                                        dirs,
                                        lambda: (
                                            dynamics_state_matrix(bundle["prop_graph_train"], X_train, cfg.tau_values, "diffusion", 0.0, cfg),
                                            dynamics_state_matrix(bundle["prop_graph_test"], X_test, cfg.tau_values, "diffusion", 0.0, cfg),
                                        ),
                                    )
                                    descriptor_specs.append(("pf_diffusion", Dtr, Dte, {"dynamics": "diffusion", "gamma": 0.0}))
                                    descriptor_specs.append(("pf_static_plus_diffusion", np.hstack([F_static_tr, Dtr]), np.hstack([F_static_te, Dte]), {"dynamics": "static+diffusion", "gamma": 0.0}))

                                if cfg.include_unitary:
                                    Utr, Ute = cached_descriptor(
                                        "unitary_wave",
                                        bundle,
                                        X_train,
                                        X_test,
                                        item_names,
                                        cfg,
                                        dirs,
                                        lambda: (
                                            dynamics_state_matrix(bundle["prop_graph_train"], X_train, cfg.tau_values, "unitary", 0.0, cfg),
                                            dynamics_state_matrix(bundle["prop_graph_test"], X_test, cfg.tau_values, "unitary", 0.0, cfg),
                                        ),
                                    )
                                    descriptor_specs.append(("pf_unitary_wave", Utr, Ute, {"dynamics": "unitary_wave", "gamma": 0.0}))
                                    descriptor_specs.append(("pf_static_plus_unitary", np.hstack([F_static_tr, Utr]), np.hstack([F_static_te, Ute]), {"dynamics": "static+unitary_wave", "gamma": 0.0}))

                                for gamma in cfg.damped_gammas:
                                    Gtr, Gte = cached_descriptor(
                                        f"damped_wave_gamma_{gamma}",
                                        bundle,
                                        X_train,
                                        X_test,
                                        item_names,
                                        cfg,
                                        dirs,
                                        lambda gamma=gamma: (
                                            dynamics_state_matrix(bundle["prop_graph_train"], X_train, cfg.tau_values, "damped", float(gamma), cfg),
                                            dynamics_state_matrix(bundle["prop_graph_test"], X_test, cfg.tau_values, "damped", float(gamma), cfg),
                                        ),
                                    )
                                    descriptor_specs.append((f"pf_damped_wave_g{gamma}", Gtr, Gte, {"dynamics": "damped_wave", "gamma": float(gamma)}))
                                    descriptor_specs.append((f"pf_static_plus_damped_g{gamma}", np.hstack([F_static_tr, Gtr]), np.hstack([F_static_te, Gte]), {"dynamics": "static+damped_wave", "gamma": float(gamma)}))

                                for base_method, Ftr, Fte, dyn_extra in descriptor_specs:
                                    for k in cfg.k_values:
                                        if k <= 0:
                                            continue
                                        k_eff = min(k, Ftr.shape[1], max(1, Ftr.shape[0] - 1))
                                        Ztr, Zte = svd_project_train_test(Ftr, Fte, k=k_eff)
                                        method = f"{base_method}_svd"
                                        rows.append(evaluate_descriptor(
                                            dataset_name, w_id, horizon, method, k_eff,
                                            Ztr, Zte, X_train, X_test, Y_train, Y_test, cfg,
                                            extra={**extra, **dyn_extra, "raw_feature_dim": int(Ftr.shape[1])},
                                        ))

                            except Exception as e:
                                print("[WARN] PF config failed:", {
                                    "dataset": dataset_name,
                                    "window": w_id,
                                    "kernel": kernel,
                                    "weight_scheme": weight_scheme,
                                    "item_identity": include_item,
                                    "stability": use_stability,
                                    "error": repr(e),
                                })

    return pd.DataFrame(rows)


def aggregate_results(df: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["dataset", "method", "k"]
    optional = ["kernel", "weight_scheme", "item_identity", "stability", "dynamics", "gamma"]
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


def run_benchmark(cfg: BenchConfig) -> pd.DataFrame:
    dirs = ensure_dirs(cfg)
    all_results = []
    for dataset_name in cfg.datasets:
        print("\n" + "=" * 90)
        print(f"[DATASET] {dataset_name}")
        X_df = cached_dataset(dataset_name, cfg, dirs)
        print(f"[SHAPE] {dataset_name}: {X_df.shape}")
        res = run_dataset(dataset_name, X_df, cfg, dirs)
        out_path = dirs["results"] / f"results_{dataset_name}_{cfg.preset}.csv"
        res.to_csv(out_path, index=False)
        print(f"[SAVE] {out_path}")
        all_results.append(res)

    combined = pd.concat(all_results, ignore_index=True)
    combined_path = dirs["results"] / f"results_combined_{cfg.preset}.csv"
    combined.to_csv(combined_path, index=False)
    agg = aggregate_results(combined)
    agg_path = dirs["results"] / f"results_aggregate_{cfg.preset}.csv"
    agg.to_csv(agg_path, index=False)

    print("\n[SAVE]", combined_path)
    print("[SAVE]", agg_path)
    print("\nTop rows by dataset:")
    print(agg.groupby("dataset", group_keys=False).head(20).to_string(index=False))
    return combined


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--preset", default="lite", choices=["lite", "markets", "full"])
    p.add_argument("--cache-dir", default=None)
    p.add_argument("--datasets", nargs="*", default=None,
                   choices=["plotly_stock_returns", "uci_household_power", "uci_appliances"])
    p.add_argument("--no-unitary", action="store_true")
    p.add_argument("--no-diffusion", action="store_true")
    p.add_argument("--no-phase", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    cfg = config_for_preset(args.preset)
    if args.cache_dir is not None:
        cfg = BenchConfig(**{**asdict(cfg), "cache_dir": args.cache_dir})
    if args.datasets is not None:
        cfg = BenchConfig(**{**asdict(cfg), "datasets": tuple(args.datasets)})
    if args.no_unitary:
        cfg = BenchConfig(**{**asdict(cfg), "include_unitary": False})
    if args.no_diffusion:
        cfg = BenchConfig(**{**asdict(cfg), "include_diffusion": False})
    if args.no_phase:
        cfg = BenchConfig(**{**asdict(cfg), "use_phase_from_delta": False})

    print("[CONFIG]")
    print(json.dumps(asdict(cfg), indent=2, default=str))
    run_benchmark(cfg)


if __name__ == "__main__":
    main()
