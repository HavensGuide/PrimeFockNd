"""
state_geometry.py
=================

Curvature-sparse defect geometry for compressed representations of
random-looking but generatively structured exponential states.

This module sits above ``tensor_defect.py`` and ``coarse_grain.py``.  It does
not claim that arbitrary Haar-random states are compressible.  Instead it gives
an explicit software object for the useful class:

    state = flat/default completion rule + sparse gluing/curvature defects

The representation tracks the distinction between:

    C_store  -- how many parameters/defects are stored,
    C_eval   -- rough cost to evaluate an amplitude/boundary configuration,
    C_learn  -- placeholder metadata for recognition/learning cost.

Supported construction paths in this v0.6 prototype:

* from a connected-correlation/Omega matrix,
* from a SimpleMPS or dense tensor via TT-SVD,
* from a simple Ising/Boltzmann-style generator,
* direct manual construction.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np

from .core import CausalDefectEngine, DefectState, Defect, EPS, NUMBA_AVAILABLE, njit
from .tensor_defect import (
    SimpleMPS,
    dense_to_mps,
    extract_omega,
    learn_mps_ordering,
    estimate_bond_dimensions,
)


@njit(cache=True)
def _ising_energy_numba(x: np.ndarray, h: np.ndarray, edges_i: np.ndarray, edges_j: np.ndarray, J: np.ndarray) -> float:
    # x is assumed encoded as 0/1; spins are -1/+1.
    e = 0.0
    for a in range(x.shape[0]):
        s = 1.0 if x[a] > 0 else -1.0
        e -= h[a] * s
    for k in range(J.shape[0]):
        si = 1.0 if x[edges_i[k]] > 0 else -1.0
        sj = 1.0 if x[edges_j[k]] > 0 else -1.0
        e -= J[k] * si * sj
    return e


@dataclass
class StateGeometryConfig:
    """Configuration for curvature-sparse state geometry."""

    mode: str = "auto"  # auto | exact | efficient | legacy
    local_dim: int = 2
    connection_threshold: float = 0.15
    curvature_tol: float = 0.05
    max_curvature_defects: Optional[int] = 512
    max_bond: Optional[int] = None
    svd_tol: float = 1e-8
    eval_rule: str = "auto"  # auto | mps | ising | holonomy_product | default
    max_dense_states: int = 1 << 20
    use_numba: bool = True
    random_state: Optional[int] = None


@dataclass
class DefectComplexity:
    """Rough complexity summary for a defect-geometric representation."""

    n_sites: int
    n_edges: int
    n_curvature_defects: int
    storage_parameters: int
    eval_cost_estimate: int
    learn_cost_estimate: Optional[int] = None
    omega_rank: Optional[int] = None
    max_bond: Optional[int] = None
    compression_ratio_vs_dense: Optional[float] = None
    notes: Dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "n_sites": self.n_sites,
            "n_edges": self.n_edges,
            "n_curvature_defects": self.n_curvature_defects,
            "storage_parameters": self.storage_parameters,
            "eval_cost_estimate": self.eval_cost_estimate,
            "learn_cost_estimate": self.learn_cost_estimate,
            "omega_rank": self.omega_rank,
            "max_bond": self.max_bond,
            "compression_ratio_vs_dense": self.compression_ratio_vs_dense,
            "notes": dict(self.notes),
        }


@dataclass
class CurvatureSparseStateGeometry:
    """Defect-geometric compressed state representation.

    Parameters
    ----------
    graph_nodes:
        Site/sector labels.
    edges:
        Adaptive gluing graph.  Edges are pairs of node indices.
    local_dims:
        Local physical/fiber dimension at each node.
    connections:
        Edge gluing maps.  For binary states, a map may be a ``2 x 2`` matrix
        with entry ``A[x_i, x_j]``.
    curvature_defects:
        Sparse non-flat loop/gluing obstructions.
    omega:
        Optional connected-correlation/defect information matrix.
    eval_rule:
        Semantics used by ``eval_amplitude``.
    mps:
        Optional MPS/TT backend.
    boundary_data:
        Generator-specific data, e.g. Ising fields/couplings or local factors.
    """

    graph_nodes: List[Any]
    edges: List[Tuple[int, int]] = field(default_factory=list)
    local_dims: List[int] = field(default_factory=list)
    connections: Dict[Tuple[int, int], np.ndarray] = field(default_factory=dict)
    curvature_defects: List[Defect] = field(default_factory=list)
    projectors: Optional[Any] = None
    boundary_data: Dict[str, Any] = field(default_factory=dict)
    renormalization: Optional[Any] = None
    eval_rule: str = "default"
    omega: Optional[np.ndarray] = None
    mps: Optional[SimpleMPS] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_sites(self) -> int:
        return len(self.graph_nodes)

    def summary(self) -> Dict[str, Any]:
        c = self.estimate_complexity()
        return {
            "n_sites": self.n_sites,
            "n_edges": len(self.edges),
            "n_curvature_defects": len(self.curvature_defects),
            "eval_rule": self.eval_rule,
            "has_mps": self.mps is not None,
            "complexity": c.as_dict(),
            "metadata": dict(self.metadata),
        }

    def storage_size(self) -> int:
        """Approximate number of stored scalar parameters."""
        size = 0
        size += len(self.graph_nodes)
        size += 2 * len(self.edges)
        for A in self.connections.values():
            size += int(np.asarray(A).size)
        size += 4 * len(self.curvature_defects)
        if self.omega is not None:
            # Count effective stored Omega size according to metadata when present.
            rank = self.metadata.get("omega_rank", None)
            n = int(np.asarray(self.omega).shape[0])
            size += int(n * rank) if rank is not None else int(np.asarray(self.omega).size)
        if self.mps is not None:
            size += int(sum(np.asarray(t).size for t in self.mps.tensors))
        for v in self.boundary_data.values():
            if isinstance(v, np.ndarray):
                size += int(v.size)
            elif isinstance(v, (list, tuple)):
                size += len(v)
            elif isinstance(v, dict):
                size += len(v)
            else:
                size += 1
        return int(size)

    def estimate_complexity(self) -> DefectComplexity:
        """Estimate storage/evaluation complexity for this representation."""
        dense = 1
        for d in (self.local_dims or [2] * self.n_sites):
            dense *= max(int(d), 1)
            if dense > 10**18:
                break
        storage = self.storage_size()
        if self.eval_rule == "mps" and self.mps is not None:
            eval_cost = int(sum(t.shape[0] * t.shape[1] * t.shape[2] for t in self.mps.tensors))
            max_bond = int(max(self.mps.bond_dims)) if self.mps.bond_dims else 1
        elif self.eval_rule == "ising":
            eval_cost = int(self.n_sites + len(self.edges))
            max_bond = self.metadata.get("max_bond")
        else:
            eval_cost = int(self.n_sites + len(self.edges) + len(self.curvature_defects))
            max_bond = self.metadata.get("max_bond")
        omega_rank = self.metadata.get("omega_rank")
        ratio = None
        if dense > 0 and storage > 0 and dense < 10**18:
            ratio = float(storage) / float(dense)
        return DefectComplexity(
            n_sites=self.n_sites,
            n_edges=len(self.edges),
            n_curvature_defects=len(self.curvature_defects),
            storage_parameters=storage,
            eval_cost_estimate=eval_cost,
            learn_cost_estimate=self.metadata.get("learn_cost_estimate"),
            omega_rank=None if omega_rank is None else int(omega_rank),
            max_bond=None if max_bond is None else int(max_bond),
            compression_ratio_vs_dense=ratio,
            notes={"dense_state_count": dense, "eval_rule": self.eval_rule},
        )

    def estimate_defect_rank(self, epsilon: float = 1e-8) -> Dict[str, Any]:
        """Estimate epsilon-defect rank via Omega rank, curvature count, and MPS bonds."""
        out: Dict[str, Any] = self.estimate_complexity().as_dict()
        if self.omega is not None:
            s = np.linalg.svd(np.asarray(self.omega, dtype=float), compute_uv=False)
            if s.size:
                thresh = max(float(epsilon), float(epsilon) * float(s[0]))
                out["r_def_omega"] = int(np.sum(s > thresh))
                out["omega_singular_values"] = s[: min(16, len(s))].tolist()
            else:
                out["r_def_omega"] = 0
        if self.mps is not None:
            out["r_def_mps_max_bond"] = int(max(self.mps.bond_dims)) if self.mps.bond_dims else 1
            out["mps_bond_dims"] = list(map(int, self.mps.bond_dims))
        out["r_def_curvature"] = int(len(self.curvature_defects))
        out["r_def_store"] = int(self.storage_size())
        return out

    def eval_amplitude(self, x: Sequence[int]) -> complex:
        """Evaluate amplitude/weight for boundary condition ``x``.

        This is intentionally semantic, not universal.  For MPS it contracts the
        selected tensor slices.  For Ising it returns ``exp(-beta E(x))`` or a
        phase if ``phase=True`` is stored in boundary data.  For the default
        holonomy product it multiplies edge gluing weights and sparse curvature
        insertions.
        """
        xx = np.asarray(x, dtype=np.int64).ravel()
        if len(xx) != self.n_sites:
            raise ValueError(f"expected boundary of length {self.n_sites}, got {len(xx)}")
        rule = self.eval_rule
        if rule == "mps":
            if self.mps is None:
                raise ValueError("eval_rule='mps' requires an mps")
            return _eval_mps_config(self.mps, xx)
        if rule == "ising":
            return self._eval_ising(xx)
        # default / holonomy_product: local factors times edge connections times defect insertions.
        amp: complex = complex(1.0)
        local = self.boundary_data.get("local_factors")
        if local is not None:
            for i, xi in enumerate(xx):
                fac = np.asarray(local[i]) if isinstance(local, (list, tuple)) else np.asarray(local)
                if fac.ndim == 2:
                    amp *= complex(fac[i, int(xi)])
                else:
                    amp *= complex(fac[int(xi)])
        for i, j in self.edges:
            A = self.connections.get((int(i), int(j)))
            if A is not None:
                AA = np.asarray(A)
                amp *= complex(AA[int(xx[i]) % AA.shape[0], int(xx[j]) % AA.shape[1]])
        for d in self.curvature_defects:
            # Sparse curvature insertion.  If parity metadata is provided, only
            # active on matching boundary parity; otherwise always active.
            parity = d.metadata.get("parity") if isinstance(d.metadata, dict) else None
            if parity is not None and d.support is not None:
                val = int(np.sum(xx[list(d.support)]) % 2)
                if val != int(parity):
                    continue
            amp *= complex(1.0 + d.charge * d.strength)
        return amp

    def _eval_ising(self, xx: np.ndarray) -> complex:
        h = np.asarray(self.boundary_data.get("h", np.zeros(self.n_sites)), dtype=float)
        beta = float(self.boundary_data.get("beta", 1.0))
        phase = bool(self.boundary_data.get("phase", False))
        J_items = self.boundary_data.get("J_items")
        if J_items is None:
            Jdict = self.boundary_data.get("J", {})
            J_items = [(int(i), int(j), float(w)) for (i, j), w in Jdict.items()]
        if len(J_items):
            ei = np.asarray([a for a, _, _ in J_items], dtype=np.int64)
            ej = np.asarray([b for _, b, _ in J_items], dtype=np.int64)
            jw = np.asarray([w for _, _, w in J_items], dtype=float)
            if NUMBA_AVAILABLE:
                e = _ising_energy_numba(xx, h, ei, ej, jw)
            else:
                spins = 2 * xx.astype(float) - 1.0
                e = -float(np.dot(h, spins))
                for i, j, w in J_items:
                    e -= float(w) * spins[i] * spins[j]
        else:
            spins = 2 * xx.astype(float) - 1.0
            e = -float(np.dot(h, spins))
        return complex(np.exp(1j * beta * e) if phase else np.exp(-beta * e))

    def sample_boundary(self, n_samples: int = 1, *, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Naive sampler for small geometries; exact enumeration when feasible.

        For large systems this falls back to independent uniform proposals.  It
        is meant for smoke tests and small toy models, not production sampling.
        """
        rr = rng or np.random.default_rng()
        dims = self.local_dims or [2] * self.n_sites
        total = int(np.prod(dims)) if self.n_sites <= 30 else 10**18
        if total <= 65536:
            configs = _enumerate_configs(dims)
            weights = np.asarray([abs(self.eval_amplitude(c)) ** 2 for c in configs], dtype=float)
            weights = weights / (weights.sum() + EPS)
            idx = rr.choice(len(configs), size=int(n_samples), replace=True, p=weights)
            return configs[idx]
        out = np.zeros((int(n_samples), self.n_sites), dtype=np.int64)
        for i, d in enumerate(dims):
            out[:, i] = rr.integers(0, int(d), size=int(n_samples))
        return out

    def to_mps(self, *, max_bond: Optional[int] = None, tol: float = 1e-10, max_states: Optional[int] = None) -> SimpleMPS:
        """Compile to a SimpleMPS when possible.

        If an MPS already exists, it is returned.  Otherwise this enumerates the
        dense tensor only for small local state spaces and TT-SVDs it.
        """
        if self.mps is not None:
            return self.mps
        dims = self.local_dims or [2] * self.n_sites
        total = int(np.prod(dims)) if self.n_sites <= 30 else 10**18
        limit = int(max_states if max_states is not None else self.metadata.get("max_dense_states", 1 << 20))
        if total > limit:
            raise ValueError(f"state space {total} exceeds max_states={limit}; cannot compile by enumeration")
        dense = np.zeros(tuple(dims), dtype=np.complex128)
        for cfg in _enumerate_configs(dims):
            dense[tuple(cfg)] = self.eval_amplitude(cfg)
        # dense_to_mps currently casts to float; keep real if possible, otherwise use magnitudes as a safe fallback.
        if np.max(np.abs(dense.imag)) < 1e-12:
            tensor = dense.real
        else:
            tensor = np.abs(dense)
            self.metadata["to_mps_complex_policy"] = "absolute_value"
        self.mps = dense_to_mps(tensor, max_bond=max_bond, tol=tol)
        self.eval_rule = "mps"
        return self.mps

    def compress_curvature(self, tol: float = 1e-8, max_defects: Optional[int] = None) -> "CurvatureSparseStateGeometry":
        """Drop weak curvature defects in-place and return self."""
        keep = [d for d in self.curvature_defects if abs(float(d.charge) * float(d.strength)) > float(tol)]
        keep.sort(key=lambda d: abs(float(d.charge) * float(d.strength)), reverse=True)
        if max_defects is not None:
            keep = keep[: int(max_defects)]
        self.curvature_defects = keep
        self.metadata["curvature_compressed_tol"] = float(tol)
        self.metadata["curvature_max_defects"] = max_defects
        return self

    def to_defect_state(self) -> DefectState:
        return DefectState(
            events=list(self.graph_nodes),
            connection=dict(self.connections),
            omega=self.omega,
            defects=list(self.curvature_defects),
            observables={
                "state_geometry_summary": self.summary(),
                "defect_rank": self.estimate_defect_rank(),
            },
            metadata={"backend": "state_geometry", **dict(self.metadata)},
            valid=True,
        )


class CurvatureSparseStateEngine(CausalDefectEngine):
    """CausalDefectEngine wrapper around CurvatureSparseStateGeometry."""

    def __init__(self, config: Optional[StateGeometryConfig] = None, **kwargs: Any):
        self.config = config or StateGeometryConfig(**{k: v for k, v in kwargs.items() if k in StateGeometryConfig.__annotations__})
        super().__init__(mode=self.config.mode, use_numba=self.config.use_numba, random_state=self.config.random_state)
        self.geometry: Optional[CurvatureSparseStateGeometry] = None
        self._source: Any = None

    def fit(self, events: Any, **kwargs: Any) -> "CurvatureSparseStateEngine":
        self._source = events
        cfg = self.config
        if isinstance(events, CurvatureSparseStateGeometry):
            geom = events
        elif isinstance(events, SimpleMPS):
            geom = geometry_from_mps(events, config=cfg)
        elif isinstance(events, Mapping) and str(events.get("type", "")).lower() in {"ising", "boltzmann", "phase_ising"}:
            geom = geometry_from_ising(events, config=cfg)
        elif isinstance(events, np.ndarray) and events.ndim >= 3:
            mps = dense_to_mps(np.asarray(events), max_bond=cfg.max_bond, tol=cfg.svd_tol)
            geom = geometry_from_mps(mps, config=cfg)
            geom.metadata["source"] = "dense_tensor"
        else:
            # A square 2D array that is already symmetric is treated as an
            # Omega/correlation matrix.  Other 2D arrays are samples x features.
            if isinstance(events, np.ndarray) and events.ndim == 2 and events.shape[0] == events.shape[1] and np.allclose(events, events.T, atol=1e-10):
                omega = np.asarray(events, dtype=float)
            else:
                omega = extract_omega(events, use_numba=self.use_numba)
            geom = geometry_from_omega(omega, config=cfg)
        self.geometry = geom
        self.state = geom.to_defect_state()
        return self

    def estimate_generators(self, **kwargs: Any) -> Any:
        if self.geometry is None:
            raise RuntimeError("call fit() first")
        gens = self.geometry.omega if self.geometry.omega is not None else np.eye(self.geometry.n_sites)
        self.state.generators = gens
        return gens

    def build_connection(self, **kwargs: Any) -> Any:
        if self.geometry is None:
            raise RuntimeError("call fit() first")
        self.state.connection = self.geometry.connections
        return self.state.connection

    def detect_defects(self, **kwargs: Any) -> List[Defect]:
        if self.geometry is None:
            raise RuntimeError("call fit() first")
        tol = kwargs.get("tol", self.config.curvature_tol)
        max_defects = kwargs.get("max_defects", self.config.max_curvature_defects)
        self.geometry.compress_curvature(tol=tol, max_defects=max_defects)
        self.state.defects = list(self.geometry.curvature_defects)
        return self.state.defects

    def eval_amplitude(self, x: Sequence[int]) -> complex:
        if self.geometry is None:
            raise RuntimeError("call fit() first")
        return self.geometry.eval_amplitude(x)

    def estimate_defect_rank(self, epsilon: float = 1e-8) -> Dict[str, Any]:
        if self.geometry is None:
            raise RuntimeError("call fit() first")
        return self.geometry.estimate_defect_rank(epsilon)

    def to_mps(self, **kwargs: Any) -> SimpleMPS:
        if self.geometry is None:
            raise RuntimeError("call fit() first")
        return self.geometry.to_mps(**kwargs)


# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------


def geometry_from_omega(omega: np.ndarray, *, config: Optional[StateGeometryConfig] = None) -> CurvatureSparseStateGeometry:
    """Build a curvature-sparse geometry from an Omega/correlation matrix."""
    cfg = config or StateGeometryConfig()
    O = np.asarray(omega, dtype=float)
    if O.ndim != 2 or O.shape[0] != O.shape[1]:
        raise ValueError("omega must be square")
    n = O.shape[0]
    W = np.abs(O.copy())
    np.fill_diagonal(W, 0.0)
    edges: List[Tuple[int, int]] = []
    connections: Dict[Tuple[int, int], np.ndarray] = {}
    for i in range(n):
        for j in range(i + 1, n):
            w = float(W[i, j])
            if w >= cfg.connection_threshold:
                edges.append((i, j))
                # Binary gluing map: flat identity-like background plus correlation coupling.
                A = np.array([[1.0 + w, 1.0 - w], [1.0 - w, 1.0 + w]], dtype=float)
                connections[(i, j)] = A
                connections[(j, i)] = A.T
    defects = _triangle_curvature_defects(W, edges, tol=cfg.curvature_tol, max_defects=cfg.max_curvature_defects)
    s = np.linalg.svd(O, compute_uv=False)
    omega_rank = int(np.sum(s > max(cfg.svd_tol, cfg.svd_tol * float(s[0] if len(s) else 1.0)))) if len(s) else 0
    ordering = learn_mps_ordering(W, method="greedy")
    bonds = estimate_bond_dimensions(O, ordering, tol=cfg.svd_tol, max_bond=cfg.max_bond)
    return CurvatureSparseStateGeometry(
        graph_nodes=list(range(n)),
        edges=edges,
        local_dims=[int(cfg.local_dim)] * n,
        connections=connections,
        curvature_defects=defects,
        omega=O,
        eval_rule="holonomy_product" if cfg.eval_rule == "auto" else cfg.eval_rule,
        metadata={
            "source": "omega",
            "omega_rank": omega_rank,
            "ordering": ordering,
            "bond_dims_est": bonds,
            "max_bond": int(max(bonds)) if bonds else 1,
            "connection_threshold": cfg.connection_threshold,
            "curvature_tol": cfg.curvature_tol,
            "max_dense_states": cfg.max_dense_states,
        },
    )


def geometry_from_mps(mps: SimpleMPS, *, config: Optional[StateGeometryConfig] = None) -> CurvatureSparseStateGeometry:
    """Wrap a SimpleMPS as a defect geometry.

    Bonds become gluing channels; large bond dimensions become bottleneck
    defects.  Evaluation is exact MPS slice contraction.
    """
    cfg = config or StateGeometryConfig()
    n = mps.n_sites
    dims = mps.physical_dims
    edges = [(i, i + 1) for i in range(max(0, n - 1))]
    conns: Dict[Tuple[int, int], np.ndarray] = {}
    # Store a lightweight scalar connection proportional to adjacent bond size.
    bonds = mps.bond_dims
    for i, j in edges:
        chi = float(bonds[i + 1] if i + 1 < len(bonds) else 1.0)
        A = np.ones((dims[i], dims[j]), dtype=float) * chi
        conns[(i, j)] = A
        conns[(j, i)] = A.T
    defects: List[Defect] = []
    if bonds:
        med = float(np.median(bonds))
        for cut, chi in enumerate(bonds):
            if chi > max(1.0, med):
                defects.append(Defect("bond_gluing_defect", charge=float(chi), support=(cut,), strength=float(chi) / (max(bonds) + EPS), metadata={"cut": cut, "bond_dim": int(chi)}))
    return CurvatureSparseStateGeometry(
        graph_nodes=list(range(n)),
        edges=edges,
        local_dims=dims,
        connections=conns,
        curvature_defects=defects,
        eval_rule="mps",
        mps=mps,
        metadata={
            "source": "mps",
            "max_bond": int(max(bonds)) if bonds else 1,
            "bond_dims": list(map(int, bonds)),
            "max_dense_states": cfg.max_dense_states,
        },
    )


def geometry_from_ising(spec: Mapping[str, Any], *, config: Optional[StateGeometryConfig] = None) -> CurvatureSparseStateGeometry:
    """Build from compact Ising/Boltzmann/phase generator.

    Spec keys:
        type: 'ising' | 'boltzmann' | 'phase_ising'
        n: number of spins
        h: optional fields length n
        J: dict {(i,j): weight} or list [(i,j,weight), ...]
        beta: inverse temperature / phase scale
        phase: bool; if True returns exp(i beta E), else exp(-beta E)
    """
    cfg = config or StateGeometryConfig()
    n = int(spec.get("n", len(spec.get("h", []))))
    if n <= 0:
        raise ValueError("Ising geometry requires n or h")
    h = np.asarray(spec.get("h", np.zeros(n)), dtype=float)
    if h.size != n:
        raise ValueError("h must have length n")
    rawJ = spec.get("J", {})
    Jdict: Dict[Tuple[int, int], float] = {}
    if isinstance(rawJ, Mapping):
        for k, v in rawJ.items():
            i, j = k
            Jdict[(int(i), int(j))] = float(v)
    else:
        for i, j, w in rawJ:
            Jdict[(int(i), int(j))] = float(w)
    edges = [(i, j) for i, j in Jdict.keys()]
    conns: Dict[Tuple[int, int], np.ndarray] = {}
    beta = float(spec.get("beta", 1.0))
    for (i, j), w in Jdict.items():
        # Local pair Boltzmann factor for binary states.
        A = np.zeros((2, 2), dtype=float)
        for xi in (0, 1):
            for xj in (0, 1):
                si = 1.0 if xi else -1.0
                sj = 1.0 if xj else -1.0
                A[xi, xj] = np.exp(beta * float(w) * si * sj)
        conns[(i, j)] = A
        conns[(j, i)] = A.T
    defects = [Defect("ising_interaction", charge=float(w), support=(i, j), strength=abs(float(w)), metadata={"edge": (i, j)}) for (i, j), w in Jdict.items() if abs(float(w)) > cfg.curvature_tol]
    return CurvatureSparseStateGeometry(
        graph_nodes=list(range(n)),
        edges=edges,
        local_dims=[2] * n,
        connections=conns,
        curvature_defects=defects[: cfg.max_curvature_defects],
        boundary_data={
            "h": h,
            "J": Jdict,
            "J_items": [(i, j, w) for (i, j), w in Jdict.items()],
            "beta": beta,
            "phase": bool(spec.get("phase", str(spec.get("type", "")).lower() == "phase_ising")),
        },
        eval_rule="ising",
        metadata={"source": "ising_generator", "max_dense_states": cfg.max_dense_states},
    )


def fit_curvature_sparse_geometry(source: Any, *, config: Optional[StateGeometryConfig] = None, **kwargs: Any) -> CurvatureSparseStateGeometry:
    """Convenience constructor mirroring CurvatureSparseStateEngine.fit."""
    cfg = config or StateGeometryConfig(**{k: v for k, v in kwargs.items() if k in StateGeometryConfig.__annotations__})
    eng = CurvatureSparseStateEngine(cfg)
    eng.fit(source)
    assert eng.geometry is not None
    return eng.geometry


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _eval_mps_config(mps: SimpleMPS, x: np.ndarray) -> complex:
    if len(x) != mps.n_sites:
        raise ValueError(f"expected {mps.n_sites} sites, got {len(x)}")
    vec = np.array([1.0], dtype=np.result_type(*[t.dtype for t in mps.tensors]))
    for site, xi in enumerate(x):
        T = np.asarray(mps.tensors[site])
        if xi < 0 or xi >= T.shape[1]:
            raise ValueError(f"site {site} index {xi} outside physical dim {T.shape[1]}")
        mat = T[:, int(xi), :]
        vec = vec @ mat
    return complex(vec.reshape(-1)[0])


def _enumerate_configs(dims: Sequence[int]) -> np.ndarray:
    grids = np.meshgrid(*[np.arange(int(d), dtype=np.int64) for d in dims], indexing="ij")
    return np.stack([g.ravel() for g in grids], axis=1)


def _triangle_curvature_defects(W: np.ndarray, edges: Sequence[Tuple[int, int]], *, tol: float, max_defects: Optional[int]) -> List[Defect]:
    n = W.shape[0]
    adj = np.zeros((n, n), dtype=bool)
    for i, j in edges:
        adj[int(i), int(j)] = True
        adj[int(j), int(i)] = True
    defects: List[Defect] = []
    # Triangle loop curvature proxy: product of signed/absolute couplings around closed loop.
    for i in range(n):
        nbrs = np.flatnonzero(adj[i])
        for a_idx in range(len(nbrs)):
            j = int(nbrs[a_idx])
            if j <= i:
                continue
            for k in nbrs[a_idx + 1 :]:
                k = int(k)
                if k <= j or not adj[j, k]:
                    continue
                curv = float(W[i, j] * W[j, k] * W[k, i])
                if abs(curv) > tol:
                    defects.append(Defect("curvature_loop", charge=curv, support=(i, j, k), strength=abs(curv), metadata={"loop": (i, j, k)}))
    defects.sort(key=lambda d: abs(d.charge) * d.strength, reverse=True)
    if max_defects is not None:
        defects = defects[: int(max_defects)]
    return defects


__all__ = [
    "StateGeometryConfig",
    "DefectComplexity",
    "CurvatureSparseStateGeometry",
    "CurvatureSparseStateEngine",
    "geometry_from_omega",
    "geometry_from_mps",
    "geometry_from_ising",
    "fit_curvature_sparse_geometry",
]
