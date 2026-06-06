"""
tensor_defect.py
================

Defect-guided tensorization utilities.

This module does not try to simulate arbitrary quantum states.  It implements a
practical bridge from a Causal Defect State to MPS-like compression:

    DefectState -> connected correlation graph/Omega -> learned ordering
    -> bond-dimension profile -> optional TT-SVD/MPS tensors.

The core idea is that MPS bond indices should represent active defect channels
across a cut, rather than arbitrary raw coordinates.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np

from .core import CausalDefectEngine, DefectState, Defect, normalize_mode, knn_indices, NUMBA_AVAILABLE, njit, prange, EPS
from .numba_kernels import corr_from_features_numba


@njit(cache=True, parallel=True)
def _corr_from_features_numba(X: np.ndarray) -> np.ndarray:
    n, d = X.shape
    means = np.zeros(d)
    for j in range(d):
        s = 0.0
        for i in range(n):
            s += X[i, j]
        means[j] = s / max(n, 1)
    stds = np.zeros(d)
    for j in range(d):
        s = 0.0
        for i in range(n):
            z = X[i, j] - means[j]
            s += z * z
        stds[j] = np.sqrt(s / max(n - 1, 1)) + 1e-12
    C = np.empty((d, d), dtype=np.float64)
    for a in prange(d):
        for b in range(d):
            s = 0.0
            for i in range(n):
                s += (X[i, a] - means[a]) * (X[i, b] - means[b])
            C[a, b] = s / max(n - 1, 1) / (stds[a] * stds[b])
    return C


@dataclass
class TensorDefectConfig:
    """Configuration for defect-guided tensorization."""

    mode: str = "auto"  # auto | exact | efficient | local | legacy
    max_bond: Optional[int] = None
    svd_tol: float = 1e-8
    corr_method: str = "auto"  # auto | omega | features | generators
    ordering_method: str = "greedy"  # greedy | spectral | natural | random
    use_numba: bool = True
    random_state: Optional[int] = None


@dataclass
class SimpleMPS:
    """Minimal MPS/TT container.

    tensors[k] has shape (left_bond, physical_dim_k, right_bond).  This object
    is deliberately small and dependency-free; it is compatible with standard
    TT-SVD output and can reconstruct modest tensors for testing.
    """

    tensors: List[np.ndarray]
    ordering: Optional[List[int]] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_sites(self) -> int:
        return len(self.tensors)

    @property
    def physical_dims(self) -> List[int]:
        return [int(t.shape[1]) for t in self.tensors]

    @property
    def bond_dims(self) -> List[int]:
        if not self.tensors:
            return []
        bonds = [int(self.tensors[0].shape[0])]
        bonds.extend(int(t.shape[2]) for t in self.tensors)
        return bonds

    def reconstruct(self) -> np.ndarray:
        """Reconstruct the dense tensor. Suitable for tests/small tensors."""
        if not self.tensors:
            return np.array([])
        acc = self.tensors[0]
        # acc shape: (1, d0, r0)
        for t in self.tensors[1:]:
            acc = np.tensordot(acc, t, axes=([-1], [0]))
            # shape accumulates as (1, d0, ..., dk, rk)
        arr = np.squeeze(acc, axis=(0, -1))
        return arr

    def entanglement_profile(self) -> List[float]:
        """Schmidt entropy proxy from singular values stored during TT-SVD."""
        svs = self.metadata.get("singular_values", [])
        out: List[float] = []
        for s in svs:
            ss = np.asarray(s, dtype=np.float64)
            if ss.size == 0:
                out.append(0.0)
                continue
            p = ss * ss
            p = p / (p.sum() + EPS)
            out.append(float(-(p * np.log(p + EPS)).sum()))
        return out


@dataclass
class TensorDefectState:
    sectors: List[Any]
    omega: np.ndarray
    ordering: List[int]
    bond_dims: List[int]
    correlation_graph: np.ndarray
    mps: Optional[SimpleMPS] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def summary(self) -> Dict[str, Any]:
        return {
            "n_sectors": len(self.sectors),
            "ordering": list(self.ordering),
            "bond_dims": list(self.bond_dims),
            "max_bond": int(max(self.bond_dims)) if self.bond_dims else 0,
            "has_mps": self.mps is not None,
            "metadata": dict(self.metadata),
        }


class TensorDefectEngine(CausalDefectEngine):
    """Backend that builds MPS-like compression data from defect correlations."""

    def __init__(self, config: Optional[TensorDefectConfig] = None, **kwargs: Any):
        self.config = config or TensorDefectConfig(**{k: v for k, v in kwargs.items() if k in TensorDefectConfig.__annotations__})
        super().__init__(
            mode=self.config.mode,
            use_numba=self.config.use_numba,
            random_state=self.config.random_state,
        )
        self.tensor_state: Optional[TensorDefectState] = None
        self._source: Any = None

    def fit(self, events: Any, **kwargs: Any) -> "TensorDefectEngine":
        """Fit from DefectState, dense feature matrix, or dense tensor.

        If ``events`` is a DefectState with ``omega`` present, that Omega is used
        as the defect-correlation field.  Otherwise a 2D array is interpreted as
        samples x features and converted to a connected correlation matrix.
        Arrays with ndim >= 3 are treated as dense tensors for TT-SVD, and a
        site correlation proxy is built by flattening each mode.
        """
        self._source = events
        omega = extract_omega(events, method=self.config.corr_method, use_numba=self.use_numba)
        sectors = infer_sectors(events, omega.shape[0])
        corr_graph = np.abs(omega).astype(np.float64)
        np.fill_diagonal(corr_graph, 0.0)
        ordering = learn_mps_ordering(corr_graph, method=self.config.ordering_method, rng=self.rng)
        bond_dims = estimate_bond_dimensions(omega, ordering, tol=self.config.svd_tol, max_bond=self.config.max_bond)
        self.tensor_state = TensorDefectState(
            sectors=sectors,
            omega=omega,
            ordering=ordering,
            bond_dims=bond_dims,
            correlation_graph=corr_graph,
            metadata={"mode": self.runtime.mode, "corr_method": self.config.corr_method, "ordering_method": self.config.ordering_method},
        )
        self.state = DefectState(
            events=events,
            generators=None,
            omega=omega,
            defects=[],
            observables={
                "ordering": ordering,
                "bond_dims": bond_dims,
                "max_bond": int(max(bond_dims)) if bond_dims else 0,
                "correlation_graph": corr_graph,
            },
            metadata={"backend": "tensor_defect", **self.tensor_state.metadata},
            valid=True,
        )
        return self

    def estimate_generators(self, **kwargs: Any) -> Any:
        if self.tensor_state is None:
            raise RuntimeError("call fit() first")
        # Sectors plus rows of Omega act as primitive correlation-channel generators.
        gens = self.tensor_state.omega[self.tensor_state.ordering]
        self.state.generators = gens
        return gens

    def build_connection(self, **kwargs: Any) -> Any:
        if self.tensor_state is None:
            raise RuntimeError("call fit() first")
        order = self.tensor_state.ordering
        G = self.tensor_state.correlation_graph
        # Chain connection uses adjacent correlations in learned order.
        conn: Dict[Tuple[int, int], float] = {}
        for a, b in zip(order[:-1], order[1:]):
            conn[(int(a), int(b))] = float(G[int(a), int(b)])
            conn[(int(b), int(a))] = float(G[int(a), int(b)])
        self.state.connection = conn
        return conn

    def detect_defects(self, **kwargs: Any) -> List[Defect]:
        """Detect high-rank cuts as tensor defects/bottlenecks."""
        if self.tensor_state is None:
            raise RuntimeError("call fit() first")
        b = self.tensor_state.bond_dims
        if not b:
            defects: List[Defect] = []
        else:
            threshold = kwargs.get("threshold", None)
            if threshold is None:
                threshold = max(2.0, float(np.median(b) + np.std(b)))
            defects = []
            for cut, chi in enumerate(b):
                if float(chi) >= float(threshold):
                    defects.append(
                        Defect(
                            kind="tensor_bottleneck",
                            charge=float(chi),
                            support=(cut,),
                            strength=float(chi) / (max(b) + EPS),
                            metadata={"cut": cut, "bond_dim": int(chi)},
                        )
                    )
        self.state.defects = defects
        return defects

    def to_mps(self, dense_tensor: Optional[np.ndarray] = None, **kwargs: Any) -> SimpleMPS:
        if dense_tensor is None:
            if isinstance(self._source, np.ndarray) and self._source.ndim >= 3:
                dense_tensor = np.asarray(self._source, dtype=np.float64)
            else:
                raise ValueError("dense_tensor is required unless fit() source was an ndim>=3 array")
        max_bond = kwargs.get("max_bond", self.config.max_bond)
        tol = kwargs.get("tol", self.config.svd_tol)
        order = kwargs.get("ordering", self.tensor_state.ordering if self.tensor_state else None)
        mps = dense_to_mps(np.asarray(dense_tensor, dtype=np.float64), max_bond=max_bond, tol=tol, ordering=order)
        if self.tensor_state is not None:
            self.tensor_state.mps = mps
            self.state.observables["mps_bond_dims"] = mps.bond_dims
            self.state.observables["mps_entropy"] = mps.entanglement_profile()
        return mps

    def holonomy(self, path: Sequence[int]) -> float:
        if self.state.connection is None:
            self.build_connection()
        total = 0.0
        for a, b in zip(path[:-1], path[1:]):
            total += float(self.state.connection.get((int(a), int(b)), 0.0))
        return total


def infer_sectors(source: Any, n: int) -> List[Any]:
    if isinstance(source, DefectState):
        if source.projectors is not None:
            try:
                keys = list(source.projectors.keys())
                if len(keys) == n:
                    return keys
            except Exception:
                pass
        if isinstance(source.observables, dict) and "sectors" in source.observables:
            sec = list(source.observables["sectors"])
            if len(sec) == n:
                return sec
    return list(range(n))


def extract_omega(source: Any, *, method: str = "auto", use_numba: bool = True) -> np.ndarray:
    """Extract or estimate a connected correlation/Omega matrix."""
    m = str(method).lower()
    if isinstance(source, DefectState):
        if source.omega is not None and m in {"auto", "omega"}:
            omega = np.asarray(source.omega, dtype=np.float64)
            return _symmetrize_square(omega)
        # Fall through to generators/events if requested/available.
        if source.generators is not None and m in {"auto", "features", "generators"}:
            return correlation_from_features(np.asarray(source.generators, dtype=np.float64), use_numba=use_numba)
        if source.events is not None:
            return extract_omega(source.events, method="auto", use_numba=use_numba)
        raise ValueError("DefectState does not contain omega, generators, or events usable for Omega extraction")

    arr = np.asarray(source, dtype=np.float64)
    if arr.ndim == 2:
        return correlation_from_features(arr, use_numba=use_numba)
    if arr.ndim >= 3:
        return mode_correlation_from_tensor(arr)
    if arr.ndim == 1:
        return correlation_from_features(arr[:, None], use_numba=use_numba)
    raise ValueError("cannot extract omega from source")


def _symmetrize_square(A: np.ndarray) -> np.ndarray:
    A = np.asarray(A, dtype=np.float64)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("omega must be a square 2D array")
    out = 0.5 * (A + A.T)
    np.fill_diagonal(out, np.diag(out))
    return out


def correlation_from_features(X: np.ndarray, *, use_numba: bool = True) -> np.ndarray:
    X = np.asarray(X, dtype=np.float64)
    if X.ndim != 2:
        raise ValueError("feature matrix must be 2D: samples x features")
    if X.shape[1] == 1:
        return np.ones((1, 1), dtype=np.float64)
    if use_numba and NUMBA_AVAILABLE:
        C = corr_from_features_numba(X)
    else:
        Xc = X - X.mean(axis=0, keepdims=True)
        std = Xc.std(axis=0, ddof=1) + EPS
        Xn = Xc / std[None, :]
        C = (Xn.T @ Xn) / max(len(X) - 1, 1)
    C = np.nan_to_num(C, nan=0.0, posinf=0.0, neginf=0.0)
    C = 0.5 * (C + C.T)
    np.fill_diagonal(C, 1.0)
    return C


def mode_correlation_from_tensor(T: np.ndarray) -> np.ndarray:
    """Correlation proxy between tensor modes/sites."""
    dims = T.shape
    n = len(dims)
    feats = []
    for axis in range(n):
        X = np.moveaxis(T, axis, 0).reshape(dims[axis], -1)
        # Use leading normalized singular profile as a mode descriptor.
        u, s, _ = np.linalg.svd(X, full_matrices=False)
        k = min(8, len(s))
        desc = np.zeros(8)
        if k:
            ss = s[:k] / (np.linalg.norm(s[:k]) + EPS)
            desc[:k] = ss
        feats.append(desc)
    return correlation_from_features(np.asarray(feats), use_numba=False)


def learn_mps_ordering(correlation_graph: np.ndarray, *, method: str = "greedy", rng: Optional[np.random.Generator] = None) -> List[int]:
    G = np.asarray(correlation_graph, dtype=np.float64)
    n = G.shape[0]
    if n == 0:
        return []
    method = str(method).lower()
    if method in {"natural", "identity", "legacy"}:
        return list(range(n))
    if method == "random":
        rr = rng or np.random.default_rng()
        return list(map(int, rr.permutation(n)))
    if method == "spectral":
        return spectral_ordering(G)
    if method in {"greedy", "auto", "efficient"}:
        return greedy_ordering(G)
    raise ValueError(f"unknown ordering_method: {method}")


def greedy_ordering(G: np.ndarray) -> List[int]:
    """Nearest-neighbor style ordering on correlation weights."""
    n = G.shape[0]
    if n <= 1:
        return list(range(n))
    # Start from strongest total-correlation sector.
    remaining = set(range(n))
    cur = int(np.argmax(G.sum(axis=1)))
    order = [cur]
    remaining.remove(cur)
    while remaining:
        nxt = max(remaining, key=lambda j: G[cur, j])
        order.append(int(nxt))
        remaining.remove(nxt)
        cur = int(nxt)
    return order


def spectral_ordering(G: np.ndarray) -> List[int]:
    """Fiedler-vector ordering of a weighted graph."""
    n = G.shape[0]
    if n <= 2:
        return list(range(n))
    W = np.maximum(G, 0.0)
    np.fill_diagonal(W, 0.0)
    D = np.diag(W.sum(axis=1))
    L = D - W
    try:
        vals, vecs = np.linalg.eigh(L)
        fiedler = vecs[:, 1]
        return list(map(int, np.argsort(fiedler)))
    except np.linalg.LinAlgError:
        return greedy_ordering(G)


def estimate_bond_dimensions(omega: np.ndarray, ordering: Sequence[int], *, tol: float = 1e-8, max_bond: Optional[int] = None) -> List[int]:
    """Estimate active defect channels across each MPS cut from Omega rank."""
    O = np.asarray(omega, dtype=np.float64)
    order = list(map(int, ordering))
    n = len(order)
    if n <= 1:
        return [1] if n == 1 else []
    dims: List[int] = []
    cap = np.inf if max_bond is None else int(max_bond)
    for cut in range(1, n):
        left = order[:cut]
        right = order[cut:]
        block = O[np.ix_(left, right)]
        if block.size == 0:
            r = 1
        else:
            s = np.linalg.svd(block, compute_uv=False)
            if s.size == 0:
                r = 1
            else:
                thresh = max(float(tol), float(tol) * float(s[0]))
                r = int(np.sum(s > thresh))
                r = max(1, r)
        dims.append(int(min(r, cap)))
    return dims


def dense_to_mps(tensor: np.ndarray, *, max_bond: Optional[int] = None, tol: float = 1e-10, ordering: Optional[Sequence[int]] = None) -> SimpleMPS:
    """TT-SVD conversion of a dense tensor to a SimpleMPS.

    Parameters
    ----------
    tensor:
        Dense tensor with one physical axis per site.
    max_bond:
        Optional maximum TT rank.
    tol:
        Relative singular-value threshold.  Singular values below
        ``tol * s[0]`` at each cut are discarded.
    ordering:
        Optional axis permutation before TT-SVD.  The returned MPS stores this
        ordering in metadata; reconstruction returns the permuted tensor.
    """
    X = np.asarray(tensor, dtype=np.float64)
    if X.ndim < 1:
        raise ValueError("tensor must have at least one axis")
    if ordering is not None:
        ordering = list(map(int, ordering))
        if len(ordering) == X.ndim:
            X = np.transpose(X, axes=ordering)
        else:
            # Defect ordering may refer to sectors rather than tensor axes; ignore safely.
            ordering = None
    dims = list(X.shape)
    n = len(dims)
    tensors: List[np.ndarray] = []
    singular_values: List[np.ndarray] = []
    rank_left = 1
    unfolding = X.reshape(1 * dims[0], -1)
    for site in range(n - 1):
        d = dims[site]
        unfolding = unfolding.reshape(rank_left * d, -1)
        U, S, Vt = np.linalg.svd(unfolding, full_matrices=False)
        if S.size == 0:
            r = 1
        else:
            thresh = float(tol) * float(S[0])
            r = int(np.sum(S > thresh))
            r = max(1, r)
        if max_bond is not None:
            r = min(r, int(max_bond))
        U = U[:, :r]
        S_keep = S[:r]
        Vt = Vt[:r, :]
        tensors.append(U.reshape(rank_left, d, r))
        singular_values.append(S_keep.copy())
        unfolding = (S_keep[:, None] * Vt)
        rank_left = r
        if site + 1 < n:
            next_dim = dims[site + 1]
            unfolding = unfolding.reshape(rank_left * next_dim, -1)
    tensors.append(unfolding.reshape(rank_left, dims[-1], 1))
    return SimpleMPS(tensors=tensors, ordering=None if ordering is None else list(ordering), metadata={"singular_values": singular_values, "tol": tol, "max_bond": max_bond})


def compress_mps(mps: SimpleMPS, *, max_bond: Optional[int] = None, tol: float = 1e-10) -> SimpleMPS:
    """Simple recompression by reconstructing and TT-SVDing. For small tests."""
    return dense_to_mps(mps.reconstruct(), max_bond=max_bond, tol=tol, ordering=None)


def connected_correlator(X: np.ndarray) -> np.ndarray:
    """Return connected covariance/correlation field from samples x channels."""
    return correlation_from_features(X, use_numba=False)


def defect_guided_mps_from_state(state: DefectState, *, max_bond: Optional[int] = None, tol: float = 1e-8, ordering_method: str = "greedy") -> TensorDefectState:
    """Convenience functional API."""
    omega = extract_omega(state)
    G = np.abs(omega)
    np.fill_diagonal(G, 0.0)
    order = learn_mps_ordering(G, method=ordering_method)
    bonds = estimate_bond_dimensions(omega, order, tol=tol, max_bond=max_bond)
    return TensorDefectState(sectors=infer_sectors(state, omega.shape[0]), omega=omega, ordering=order, bond_dims=bonds, correlation_graph=G, metadata={"source": "DefectState"})
