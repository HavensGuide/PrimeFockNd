"""
core.py
=======

Base abstractions and shared accelerated kernels for the Causal Defect Algebra
Engine. Backends expose selectable execution modes so users can compare legacy
all-pairs behavior against faster neighbor-local/lazy modes.

Common mode names
-----------------
legacy / exact
    Preserve the original all-pairs semantics. Useful for regression tests and
    small problems.
local / efficient
    Prefer k-nearest-neighbor/local calculations and lazy representations.
auto
    Pick a safe efficient mode when optional dependencies are present, otherwise
    fall back to exact/legacy behavior.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np

try:  # optional acceleration
    from numba import njit, prange  # type: ignore
    NUMBA_AVAILABLE = True
except Exception:  # pragma: no cover
    NUMBA_AVAILABLE = False

    def njit(*args, **kwargs):  # type: ignore
        def deco(fn):
            return fn
        if args and callable(args[0]):
            return args[0]
        return deco

    prange = range  # type: ignore

try:  # optional nearest-neighbor acceleration
    from scipy.spatial import cKDTree  # type: ignore
    SCIPY_KDTREE_AVAILABLE = True
except Exception:  # pragma: no cover
    cKDTree = None
    SCIPY_KDTREE_AVAILABLE = False

EPS = 1e-12
ArrayLike = Union[np.ndarray, Sequence[float], Sequence[Sequence[float]]]
Mode = str


def normalize_mode(mode: Optional[str], *, default: str = "auto") -> str:
    """Canonicalize user-facing execution modes."""
    if mode is None:
        mode = default
    m = str(mode).lower().strip()
    aliases = {
        "old": "legacy",
        "original": "legacy",
        "brute": "legacy",
        "allpairs": "legacy",
        "all_pairs": "legacy",
        "full": "exact",
        "fast": "efficient",
        "knn": "local",
        "neighbor": "local",
        "neighbors": "local",
    }
    return aliases.get(m, m)


@dataclass
class RuntimeConfig:
    """Shared runtime knobs used by all backends."""

    mode: str = "auto"
    use_numba: bool = True
    neighbor_backend: str = "auto"  # auto | scipy | brute
    chunk_size: int = 2048
    lazy: bool = True

    def resolved_mode(self) -> str:
        return normalize_mode(self.mode)


@dataclass
class Defect:
    """A lightweight defect record shared by all backends."""

    kind: str
    charge: float = 0.0
    location: Optional[np.ndarray] = None
    support: Optional[Tuple[int, ...]] = None
    strength: float = 1.0
    metadata: Dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, Any]:
        out = {
            "kind": self.kind,
            "charge": self.charge,
            "strength": self.strength,
            "metadata": dict(self.metadata),
        }
        if self.location is not None:
            out["location"] = np.asarray(self.location).tolist()
        if self.support is not None:
            out["support"] = tuple(self.support)
        return out


@dataclass
class DefectState:
    """Domain-independent container for a fitted defect-algebra state."""

    events: Any = None
    generators: Any = None
    connection: Any = None
    fusion_rules: Any = None
    projectors: Any = None
    omega: Any = None
    defects: List[Defect] = field(default_factory=list)
    observables: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)
    valid: bool = True

    def summary(self) -> Dict[str, Any]:
        return {
            "valid": self.valid,
            "n_defects": len(self.defects),
            "observables": list(self.observables.keys()),
            "metadata": dict(self.metadata),
        }


class CompositionRule(ABC):
    """Abstract composition law for defect generators."""

    @abstractmethod
    def identity(self) -> Any:
        """Return the vacuum/identity element."""

    @abstractmethod
    def compose(self, a: Any, b: Any) -> Any:
        """Compose two elements."""

    def invert(self, a: Any) -> Any:
        raise NotImplementedError(f"{self.__class__.__name__} has no inverse")

    def power(self, a: Any, k: int) -> Any:
        if k < 0:
            return self.power(self.invert(a), -k)
        out = self.identity()
        base = a
        while k:
            if k & 1:
                out = self.compose(out, base)
            base = self.compose(base, base)
            k >>= 1
        return out


class CausalDefectEngine(ABC):
    """Base class for all defect-algebra engines."""

    def __init__(
        self,
        *,
        mode: str = "auto",
        use_numba: bool = True,
        neighbor_backend: str = "auto",
        chunk_size: int = 2048,
        lazy: bool = True,
        random_state: Optional[int] = None,
    ):
        self.runtime = RuntimeConfig(
            mode=normalize_mode(mode),
            use_numba=bool(use_numba),
            neighbor_backend=neighbor_backend,
            chunk_size=int(chunk_size),
            lazy=bool(lazy),
        )
        self.use_numba = bool(use_numba and NUMBA_AVAILABLE)
        self.rng = np.random.default_rng(random_state)
        self.state = DefectState(valid=False, metadata={"numba": self.use_numba, "mode": self.runtime.mode})

    @abstractmethod
    def fit(self, events: Any, **kwargs) -> "CausalDefectEngine":
        """Fit/prepare backend event space."""

    @abstractmethod
    def estimate_generators(self, **kwargs) -> Any:
        """Estimate local primitive generators."""

    @abstractmethod
    def build_connection(self, **kwargs) -> Any:
        """Build transport/connection/composition data."""

    @abstractmethod
    def detect_defects(self, **kwargs) -> List[Defect]:
        """Detect defects, neutral cycles, or holonomy obstructions."""

    def set_mode(self, mode: str) -> "CausalDefectEngine":
        """Switch execution mode before re-running fit/run."""
        self.runtime.mode = normalize_mode(mode)
        self.state.metadata["mode"] = self.runtime.mode
        return self

    def compose(self, a: Any, b: Any) -> Any:
        if self.state.fusion_rules is None:
            raise RuntimeError("No fusion/composition rule has been installed")
        return self.state.fusion_rules.compose(a, b)

    def holonomy(self, path: Sequence[int]) -> Any:
        raise NotImplementedError(f"{self.__class__.__name__} does not implement holonomy")

    def project(self, sector: Any, x: Any) -> Any:
        if self.state.projectors is None:
            raise RuntimeError("No projectors installed")
        return self.state.projectors[sector](x)

    def observables(self) -> Dict[str, Any]:
        return dict(self.state.observables)

    def run(self, events: Any, **kwargs) -> DefectState:
        """Convenience full pipeline."""
        self.fit(events, **kwargs)
        self.estimate_generators(**kwargs)
        self.build_connection(**kwargs)
        self.detect_defects(**kwargs)
        return self.state


# ---------------------------------------------------------------------------
# Common Numba kernels and neighborhood helpers
# ---------------------------------------------------------------------------


@njit(cache=True)
def wrap_angle(x: float) -> float:
    """Wrap an angle to [-pi, pi]."""
    return (x + np.pi) % (2.0 * np.pi) - np.pi


@njit(cache=True, parallel=True)
def pairwise_sqdist_numba(X: np.ndarray) -> np.ndarray:
    """Pairwise squared Euclidean distances."""
    n, d = X.shape
    out = np.empty((n, n), dtype=np.float64)
    for i in prange(n):
        out[i, i] = 0.0
        for j in range(i + 1, n):
            s = 0.0
            for k in range(d):
                diff = X[i, k] - X[j, k]
                s += diff * diff
            out[i, j] = s
            out[j, i] = s
    return out


@njit(cache=True, parallel=True)
def knn_bruteforce_numba(X: np.ndarray, k: int) -> Tuple[np.ndarray, np.ndarray]:
    """Brute-force kNN with O(n^2 d) time but O(nk) output memory."""
    n, d = X.shape
    kk = min(k, n - 1)
    idx = np.empty((n, kk), dtype=np.int64)
    dst = np.empty((n, kk), dtype=np.float64)
    for i in prange(n):
        best_d = np.full(kk, np.inf, dtype=np.float64)
        best_i = np.full(kk, -1, dtype=np.int64)
        for j in range(n):
            if i == j:
                continue
            s = 0.0
            for c in range(d):
                diff = X[i, c] - X[j, c]
                s += diff * diff
            # insert into descending worst slot
            worst = 0
            for q in range(1, kk):
                if best_d[q] > best_d[worst]:
                    worst = q
            if s < best_d[worst]:
                best_d[worst] = s
                best_i[worst] = j
        # sort local neighbors
        order = np.argsort(best_d)
        for q in range(kk):
            idx[i, q] = best_i[order[q]]
            dst[i, q] = best_d[order[q]]
    return idx, dst


def pairwise_sqdist(X: np.ndarray, use_numba: bool = True) -> np.ndarray:
    X = np.asarray(X, dtype=np.float64)
    if use_numba and NUMBA_AVAILABLE:
        return pairwise_sqdist_numba(X)
    diff = X[:, None, :] - X[None, :, :]
    return np.einsum("ijk,ijk->ij", diff, diff)


def knn_indices(
    X: np.ndarray,
    k: int,
    *,
    backend: str = "auto",
    use_numba: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return k nearest-neighbor indices and squared distances.

    In efficient/local modes this avoids materializing the full n x n matrix if
    SciPy's cKDTree is available. Otherwise it falls back to brute force with
    O(nk) output memory.
    """
    X = np.asarray(X, dtype=np.float64)
    n = len(X)
    if n <= 1:
        return np.zeros((n, 0), dtype=np.int64), np.zeros((n, 0), dtype=np.float64)
    kk = min(int(k), n - 1)
    b = str(backend).lower()
    if b in {"auto", "scipy", "kdtree"} and SCIPY_KDTREE_AVAILABLE:
        tree = cKDTree(X)
        d, idx = tree.query(X, k=kk + 1)
        return idx[:, 1:].astype(np.int64), (d[:, 1:] ** 2).astype(np.float64)
    if use_numba and NUMBA_AVAILABLE:
        return knn_bruteforce_numba(X, kk)
    d2 = pairwise_sqdist(X, use_numba=False)
    np.fill_diagonal(d2, np.inf)
    idx = np.argpartition(d2, kk, axis=1)[:, :kk]
    dst = np.take_along_axis(d2, idx, axis=1)
    order = np.argsort(dst, axis=1)
    return np.take_along_axis(idx, order, axis=1).astype(np.int64), np.take_along_axis(dst, order, axis=1)


@njit(cache=True)
def winding_loop_numba(theta: np.ndarray, tri: np.ndarray) -> float:
    i, j, k = int(tri[0]), int(tri[1]), int(tri[2])
    return (
        wrap_angle(theta[j] - theta[i])
        + wrap_angle(theta[k] - theta[j])
        + wrap_angle(theta[i] - theta[k])
    ) / (2.0 * np.pi)


@njit(cache=True)
def bit_count_uint64(x: np.uint64) -> int:
    c = 0
    y = x
    while y:
        y &= y - np.uint64(1)
        c += 1
    return c
