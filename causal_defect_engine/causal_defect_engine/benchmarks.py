"""
benchmarks.py
=============

Small proof-of-work benchmarks for the Causal Defect Engine.

These are intentionally lightweight.  They are meant to verify that each module
runs and to provide rough timing/compression diagnostics, not to serve as a
formal performance suite.
"""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter
from typing import Any, Dict, Optional

import numpy as np

from .core import Defect, DefectState
from .memory import PleromaMemoryEngine, MemoryConfig
from .coarse_grain import coarse_grain_state, CoarseGrainConfig
from .optim import neutral_cycle_search, optimize_mps_ordering
from .state_geometry import CurvatureSparseStateEngine, StateGeometryConfig


@dataclass
class BenchmarkResult:
    name: str
    seconds: float
    metrics: Dict[str, Any]


def benchmark_memory(n: int = 256, d: int = 16, *, seed: int = 0) -> BenchmarkResult:
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, d))
    eng = PleromaMemoryEngine(MemoryConfig(dim=d, max_items=n, compression_rank=min(16, d)))
    t0 = perf_counter()
    for i, x in enumerate(X):
        eng.write(x, key=f"x{i}")
    recall = eng.recall(X[0], k=5)
    dt = perf_counter() - t0
    return BenchmarkResult("memory", dt, {"n": n, "d": d, "top_key": recall.keys[0] if recall.keys else None, "top_score": float(recall.scores[0]) if len(recall.scores) else None})


def benchmark_coarse_grain(n_defects: int = 128, d: int = 2, *, seed: int = 1) -> BenchmarkResult:
    rng = np.random.default_rng(seed)
    defects = [Defect("toy", charge=float(rng.choice([-1, 1])), location=rng.random(d), strength=1.0) for _ in range(n_defects)]
    omega = rng.standard_normal((64, 64))
    omega = (omega + omega.T) / 2
    state = DefectState(events=list(range(64)), omega=omega, defects=defects, valid=True)
    t0 = perf_counter()
    out, report = coarse_grain_state(state, CoarseGrainConfig(merge_radius=0.08, omega_rank=8))
    dt = perf_counter() - t0
    return BenchmarkResult("coarse_grain", dt, {"defects_in": report.n_defects_in, "defects_out": report.n_defects_out, "omega_rank": report.omega_rank, "energy": report.retained_energy})


def benchmark_gf2(rows: int = 96, cols: int = 64, *, seed: int = 2) -> BenchmarkResult:
    rng = np.random.default_rng(seed)
    M = (rng.random((rows, cols)) < 0.08).astype(np.uint8)
    t0 = perf_counter()
    ns = neutral_cycle_search(M, left=True)
    dt = perf_counter() - t0
    return BenchmarkResult("gf2_neutral_cycles", dt, {"rows": rows, "cols": cols, "basis_dim": int(len(ns))})


def benchmark_ordering(n: int = 48, *, seed: int = 3) -> BenchmarkResult:
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, 4))
    omega = np.abs(X @ X.T)
    t0 = perf_counter()
    res = optimize_mps_ordering(omega, random_state=seed, max_iter=200)
    dt = perf_counter() - t0
    return BenchmarkResult("mps_ordering", dt, {"n": n, "score": float(res.score), "success": res.success})


def benchmark_state_geometry(n: int = 32, *, seed: int = 4) -> BenchmarkResult:
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, 4))
    omega = X @ X.T
    omega = omega / (np.max(np.abs(omega)) + 1e-12)
    omega = 0.5 * (omega + omega.T)
    t0 = perf_counter()
    eng = CurvatureSparseStateEngine(StateGeometryConfig(connection_threshold=0.25, curvature_tol=0.05, max_curvature_defects=64))
    eng.fit(omega).estimate_generators()
    eng.build_connection()
    eng.detect_defects()
    rank = eng.estimate_defect_rank()
    dt = perf_counter() - t0
    return BenchmarkResult("state_geometry", dt, {"n": n, "edges": len(eng.geometry.edges), "curvature_defects": len(eng.geometry.curvature_defects), "r_def_store": rank["r_def_store"]})


def run_all() -> Dict[str, BenchmarkResult]:
    return {
        "memory": benchmark_memory(),
        "coarse_grain": benchmark_coarse_grain(),
        "gf2": benchmark_gf2(),
        "ordering": benchmark_ordering(),
        "state_geometry": benchmark_state_geometry(),
    }


__all__ = [
    "BenchmarkResult",
    "benchmark_memory",
    "benchmark_coarse_grain",
    "benchmark_gf2",
    "benchmark_ordering",
    "benchmark_state_geometry",
    "run_all",
]
