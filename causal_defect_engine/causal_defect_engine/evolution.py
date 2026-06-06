r"""
evolution.py
============

Evolution mechanisms for Causal Defect Engine objects.

The package already had representation, memory, optimization, coarse graining,
and geometry compilation.  This module adds the missing dynamics layer:

    D_{t+dt} = R_epsilon \circ G \circ Phi_dt(D_t)

where Phi_dt is a raw evolution rule, G is optional gauge fixing, and R_epsilon
is optional recompression/coarse-graining.

The implementation is deliberately conservative:

* it operates on the common DefectState container and on
  CurvatureSparseStateGeometry when available,
* it keeps legacy/exact behavior possible through explicit configuration,
* it uses Numba kernels for edge diffusion, edge transport, curvature energy,
  and dense Euler steps when Numba is installed,
* it falls back to NumPy/Python on machines without Numba.

This is not meant to be a full PDE or quantum simulator.  It is the first
unified interface for evolving compressed defect geometries.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union
import copy

import numpy as np

from .core import Defect, DefectState, EPS, NUMBA_AVAILABLE, njit, normalize_mode

try:  # optional imports; keep this module importable if higher layers change
    from .coarse_grain import CoarseGrainConfig, coarse_grain_state
except Exception:  # pragma: no cover
    CoarseGrainConfig = None  # type: ignore
    coarse_grain_state = None  # type: ignore

try:
    from .geometry_compiler import CompilerConfig, gauge_fix_geometry, complexity_report
    from .state_geometry import CurvatureSparseStateGeometry
except Exception:  # pragma: no cover
    CompilerConfig = None  # type: ignore
    gauge_fix_geometry = None  # type: ignore
    complexity_report = None  # type: ignore
    CurvatureSparseStateGeometry = None  # type: ignore


# ---------------------------------------------------------------------------
# Numba kernels
# ---------------------------------------------------------------------------


@njit(cache=True)
def _edge_diffusion_step_numba(q: np.ndarray, ei: np.ndarray, ej: np.ndarray, w: np.ndarray, diffusion: float, dt: float) -> np.ndarray:
    out = q.copy()
    for k in range(ei.shape[0]):
        i = ei[k]
        j = ej[k]
        wij = w[k]
        dq = diffusion * dt * wij * (q[j] - q[i])
        out[i] += dq
        out[j] -= dq
    return out


@njit(cache=True)
def _edge_transport_features_numba(X: np.ndarray, ei: np.ndarray, ej: np.ndarray, w: np.ndarray, dt: float, damping: float) -> np.ndarray:
    out = (1.0 - damping * dt) * X.copy()
    n_edges = ei.shape[0]
    n_feat = X.shape[1]
    for k in range(n_edges):
        i = ei[k]
        j = ej[k]
        wij = w[k] * dt
        for a in range(n_feat):
            msg_ij = wij * X[i, a]
            msg_ji = wij * X[j, a]
            out[j, a] += msg_ij
            out[i, a] += msg_ji
    return out


@njit(cache=True)
def _dense_euler_numba(x: np.ndarray, A: np.ndarray, dt: float) -> np.ndarray:
    n = A.shape[0]
    out = x.copy()
    if x.ndim == 1:
        for i in range(n):
            acc = 0.0
            for j in range(n):
                acc += A[i, j] * x[j]
            out[i] += dt * acc
    else:
        m = x.shape[1]
        for i in range(n):
            for a in range(m):
                acc = 0.0
                for j in range(n):
                    acc += A[i, j] * x[j, a]
                out[i, a] += dt * acc
    return out


@njit(cache=True)
def _curvature_energy_numba(omega: np.ndarray) -> float:
    e = 0.0
    for i in range(omega.shape[0]):
        for j in range(omega.shape[1]):
            e += omega[i, j] * omega[i, j]
    return e


@njit(cache=True)
def _reaction_clip_numba(q: np.ndarray, growth: float, saturation: float, dt: float) -> np.ndarray:
    out = q.copy()
    for i in range(q.shape[0]):
        # simple stable cubic reaction: q' = q + dt*(growth*q - saturation*q^3)
        out[i] = q[i] + dt * (growth * q[i] - saturation * q[i] * q[i] * q[i])
    return out


# ---------------------------------------------------------------------------
# Config and reports
# ---------------------------------------------------------------------------


@dataclass
class EvolutionConfig:
    """Shared configuration for defect evolution rules."""

    mode: str = "auto"  # auto | efficient | exact | legacy
    dt: float = 1.0
    steps: int = 1
    use_numba: bool = True
    copy: bool = True
    renormalize: bool = False
    gauge_fix: bool = False
    coarse_grain: bool = False
    coarse_omega_rank: Optional[int] = None
    merge_radius: float = 0.1
    track_history: bool = True
    random_state: Optional[int] = None

    def resolved_mode(self) -> str:
        return normalize_mode(self.mode)


@dataclass
class EvolutionReport:
    """Diagnostics emitted by an evolution run."""

    rule: str
    steps: int
    dt: float
    energies: List[float] = field(default_factory=list)
    n_defects: List[int] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "rule": self.rule,
            "steps": self.steps,
            "dt": self.dt,
            "energies": list(self.energies),
            "n_defects": list(self.n_defects),
            "metadata": dict(self.metadata),
        }


StateLike = Union[DefectState, Any]


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def _clone_state_like(obj: StateLike) -> StateLike:
    if isinstance(obj, DefectState):
        return DefectState(
            events=obj.events,
            generators=None if obj.generators is None else np.array(obj.generators, copy=True) if isinstance(obj.generators, np.ndarray) else copy.deepcopy(obj.generators),
            connection=copy.deepcopy(obj.connection),
            fusion_rules=obj.fusion_rules,
            projectors=copy.deepcopy(obj.projectors),
            omega=None if obj.omega is None else np.array(obj.omega, copy=True) if isinstance(obj.omega, np.ndarray) else copy.deepcopy(obj.omega),
            defects=copy.deepcopy(obj.defects),
            observables=copy.deepcopy(obj.observables),
            metadata=copy.deepcopy(obj.metadata),
            valid=obj.valid,
        )
    return copy.deepcopy(obj)


def _is_geometry(obj: Any) -> bool:
    return CurvatureSparseStateGeometry is not None and isinstance(obj, CurvatureSparseStateGeometry)


def _state_energy(obj: StateLike) -> float:
    omega = getattr(obj, "omega", None)
    if omega is not None:
        O = np.asarray(omega, dtype=float)
        if O.ndim == 2:
            if NUMBA_AVAILABLE:
                return float(_curvature_energy_numba(O))
            return float(np.sum(O * O))
    defects = getattr(obj, "defects", None)
    if defects is None and _is_geometry(obj):
        defects = getattr(obj, "curvature_defects", [])
    if defects is not None:
        return float(sum((float(d.charge) ** 2 + float(d.strength) ** 2) for d in defects))
    return 0.0


def _n_defects(obj: StateLike) -> int:
    if isinstance(obj, DefectState):
        return len(obj.defects)
    if _is_geometry(obj):
        return len(getattr(obj, "curvature_defects", []))
    return 0


def _extract_edges(connection: Any = None, *, edges: Optional[Sequence[Tuple[int, int]]] = None, n: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Extract undirected edge arrays and weights from flexible connection data."""
    edge_list: List[Tuple[int, int]] = []
    weights: List[float] = []

    if edges is not None:
        for e in edges:
            if len(e) >= 2:
                edge_list.append((int(e[0]), int(e[1])))
                weights.append(1.0)

    if connection is not None:
        if isinstance(connection, np.ndarray):
            A = np.asarray(connection)
            if A.ndim == 2 and A.shape[0] == A.shape[1]:
                if n is None:
                    n = A.shape[0]
                ii, jj = np.nonzero(np.triu(np.abs(A) > EPS, k=1))
                for i, j in zip(ii, jj):
                    edge_list.append((int(i), int(j)))
                    weights.append(float(abs(A[i, j])))
        elif isinstance(connection, Mapping):
            if "edges" in connection:
                conn_edges = connection.get("edges")
                conn_w = connection.get("weights", None)
                for idx, e in enumerate(conn_edges):
                    if len(e) >= 2:
                        edge_list.append((int(e[0]), int(e[1])))
                        if conn_w is None:
                            weights.append(1.0)
                        else:
                            weights.append(float(np.asarray(conn_w).ravel()[idx]))
            elif "adjacency" in connection:
                A = np.asarray(connection["adjacency"])
                if n is None and A.ndim == 2:
                    n = A.shape[0]
                if A.ndim == 2:
                    ii, jj = np.nonzero(np.triu(np.abs(A) > EPS, k=1))
                    for i, j in zip(ii, jj):
                        edge_list.append((int(i), int(j)))
                        weights.append(float(abs(A[i, j])))
            else:
                # Mapping from edge tuple -> matrix/phase/weight.
                for key, val in connection.items():
                    if isinstance(key, tuple) and len(key) >= 2:
                        i, j = int(key[0]), int(key[1])
                        edge_list.append((i, j))
                        arr = np.asarray(val)
                        if arr.ndim == 0:
                            weights.append(float(abs(arr)))
                        else:
                            weights.append(float(np.linalg.norm(arr)))

    if n is None:
        max_idx = -1
        for i, j in edge_list:
            max_idx = max(max_idx, i, j)
        n = max_idx + 1 if max_idx >= 0 else 0

    if not edge_list:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64), np.empty(0, dtype=float), int(n or 0)

    # De-duplicate roughly, accumulating weights.
    acc: Dict[Tuple[int, int], float] = {}
    for (i, j), w in zip(edge_list, weights):
        if i == j:
            continue
        a, b = (i, j) if i < j else (j, i)
        acc[(a, b)] = acc.get((a, b), 0.0) + float(w)
    keys = sorted(acc.keys())
    ei = np.array([k[0] for k in keys], dtype=np.int64)
    ej = np.array([k[1] for k in keys], dtype=np.int64)
    ww = np.array([acc[k] for k in keys], dtype=float)
    if n is None or n == 0:
        n = int(max(ei.max(initial=-1), ej.max(initial=-1)) + 1)
    return ei, ej, ww, int(n)


def graph_laplacian_from_edges(n: int, edges: Sequence[Tuple[int, int]], weights: Optional[Sequence[float]] = None) -> np.ndarray:
    """Build a dense graph Laplacian for small/medium graphs."""
    L = np.zeros((int(n), int(n)), dtype=float)
    if weights is None:
        weights = [1.0] * len(edges)
    for (i, j), w in zip(edges, weights):
        i, j = int(i), int(j)
        wij = float(w)
        if i == j:
            continue
        L[i, i] += wij
        L[j, j] += wij
        L[i, j] -= wij
        L[j, i] -= wij
    return L


def defect_charges_to_vector(defects: Sequence[Defect], n: Optional[int] = None) -> np.ndarray:
    """Project defect charges onto node support vector."""
    if n is None:
        n = 0
        for d in defects:
            if d.support:
                n = max(n, max(d.support) + 1)
    q = np.zeros(int(n or 0), dtype=float)
    for d in defects:
        if d.support:
            supp = [int(s) for s in d.support if int(s) < len(q)]
            if supp:
                q[supp] += float(d.charge) / max(len(supp), 1)
    return q


def vector_to_defects(q: np.ndarray, *, kind: str = "diffused_defect", threshold: float = 1e-9) -> List[Defect]:
    """Convert a charge vector back into point-supported defects."""
    out: List[Defect] = []
    for i, val in enumerate(np.asarray(q, dtype=float).ravel()):
        if abs(val) > threshold:
            out.append(Defect(kind=kind, charge=float(val), support=(int(i),), strength=float(abs(val)), metadata={"source": "evolution"}))
    return out


# ---------------------------------------------------------------------------
# Base rule and concrete rules
# ---------------------------------------------------------------------------


class EvolutionRule:
    """Base class for defect evolution rules."""

    name = "evolution"

    def __init__(self, config: Optional[EvolutionConfig] = None, **kwargs: Any):
        cfg = config or EvolutionConfig()
        for k, v in kwargs.items():
            if hasattr(cfg, k):
                setattr(cfg, k, v)
        self.config = cfg
        self.use_numba = bool(cfg.use_numba and NUMBA_AVAILABLE)

    def step(self, state: StateLike, dt: Optional[float] = None) -> StateLike:
        raise NotImplementedError

    def energy(self, state: StateLike) -> float:
        return _state_energy(state)

    def evolve(self, state: StateLike, *, steps: Optional[int] = None, dt: Optional[float] = None) -> Tuple[StateLike, EvolutionReport]:
        cfg = self.config
        steps_i = int(cfg.steps if steps is None else steps)
        dt_f = float(cfg.dt if dt is None else dt)
        cur = _clone_state_like(state) if cfg.copy else state
        report = EvolutionReport(rule=self.name, steps=steps_i, dt=dt_f, metadata={"mode": cfg.resolved_mode(), "numba": self.use_numba})
        if cfg.track_history:
            report.energies.append(self.energy(cur))
            report.n_defects.append(_n_defects(cur))
        for _ in range(steps_i):
            cur = self.step(cur, dt=dt_f)
            cur = _postprocess_evolved(cur, cfg)
            if cfg.track_history:
                report.energies.append(self.energy(cur))
                report.n_defects.append(_n_defects(cur))
        return cur, report


class TransportEvolution(EvolutionRule):
    """Transport generators/features over a graph or connection.

    For ``DefectState`` objects with an array ``generators`` this performs
    message passing along ``state.connection``.  If ``connection`` is a dense
    square matrix, an Euler step ``X <- X + dt A X`` is used.  If connection is
    an edge list/dict, neighbor transport is used.

    For ``CurvatureSparseStateGeometry`` objects, the same operation transports
    a ``boundary_data['features']`` array when present; otherwise it records the
    operation in metadata without inventing features.
    """

    name = "transport"

    def __init__(self, config: Optional[EvolutionConfig] = None, *, damping: float = 0.0, connection: Any = None, **kwargs: Any):
        super().__init__(config, **kwargs)
        self.damping = float(damping)
        self.override_connection = connection

    def step(self, state: StateLike, dt: Optional[float] = None) -> StateLike:
        dt_f = float(self.config.dt if dt is None else dt)
        if isinstance(state, DefectState):
            X = state.generators
            if X is None:
                state.metadata["last_evolution"] = "transport:no_generators"
                return state
            Xarr = np.asarray(X, dtype=float)
            was_1d = Xarr.ndim == 1
            X2 = Xarr.reshape((-1, 1)) if was_1d else Xarr
            conn = self.override_connection if self.override_connection is not None else state.connection
            if isinstance(conn, np.ndarray) and conn.ndim == 2 and conn.shape[0] == X2.shape[0]:
                if self.use_numba:
                    out = _dense_euler_numba(X2, np.asarray(conn, dtype=float), dt_f)
                else:
                    out = X2 + dt_f * (np.asarray(conn, dtype=float) @ X2)
                if self.damping:
                    out *= max(0.0, 1.0 - self.damping * dt_f)
            else:
                ei, ej, w, _ = _extract_edges(conn, n=X2.shape[0])
                if len(ei) == 0:
                    out = X2.copy()
                elif self.use_numba:
                    out = _edge_transport_features_numba(X2, ei, ej, w, dt_f, self.damping)
                else:
                    out = (1.0 - self.damping * dt_f) * X2.copy()
                    for i, j, wij in zip(ei, ej, w):
                        out[j] += dt_f * wij * X2[i]
                        out[i] += dt_f * wij * X2[j]
            state.generators = out.ravel() if was_1d else out
            state.metadata["last_evolution"] = "transport"
            return state

        if _is_geometry(state):
            feats = state.boundary_data.get("features") if hasattr(state, "boundary_data") else None
            if feats is not None:
                Xarr = np.asarray(feats, dtype=float)
                was_1d = Xarr.ndim == 1
                X2 = Xarr.reshape((-1, 1)) if was_1d else Xarr
                ei, ej, w, _ = _extract_edges(state.connections, edges=state.edges, n=len(state.graph_nodes))
                if len(ei):
                    if self.use_numba:
                        out = _edge_transport_features_numba(X2, ei, ej, w, dt_f, self.damping)
                    else:
                        out = (1.0 - self.damping * dt_f) * X2.copy()
                        for i, j, wij in zip(ei, ej, w):
                            out[j] += dt_f * wij * X2[i]
                            out[i] += dt_f * wij * X2[j]
                    state.boundary_data["features"] = out.ravel() if was_1d else out
            state.metadata["last_evolution"] = "transport"
            return state
        return state


class CurvatureRelaxation(EvolutionRule):
    """Relax curvature/connection inconsistency by shrinking Omega and weak defects.

    This is a simple, stable proxy for gradient descent on

        E_curv = ||Omega||_F^2 + lambda * sum_a strength_a^2.

    The rule is intentionally local and cheap.  Stable/topological defects can
    be preserved by setting ``preserve_topological=True`` or marking defects
    with ``metadata['stable']=True``.
    """

    name = "curvature_relaxation"

    def __init__(
        self,
        config: Optional[EvolutionConfig] = None,
        *,
        eta: float = 0.1,
        defect_decay: float = 0.05,
        preserve_topological: bool = True,
        min_strength: float = 1e-10,
        **kwargs: Any,
    ):
        super().__init__(config, **kwargs)
        self.eta = float(eta)
        self.defect_decay = float(defect_decay)
        self.preserve_topological = bool(preserve_topological)
        self.min_strength = float(min_strength)

    def _relax_defects(self, defects: Sequence[Defect], dt_f: float) -> List[Defect]:
        out: List[Defect] = []
        factor = max(0.0, 1.0 - self.defect_decay * dt_f)
        for d in defects:
            stable = bool(d.metadata.get("stable", False)) or (self.preserve_topological and "topological" in d.kind)
            nd = copy.deepcopy(d)
            if not stable:
                nd.strength = float(nd.strength) * factor
                nd.charge = float(nd.charge) * factor
            if abs(nd.strength) >= self.min_strength or stable:
                out.append(nd)
        return out

    def step(self, state: StateLike, dt: Optional[float] = None) -> StateLike:
        dt_f = float(self.config.dt if dt is None else dt)
        factor = max(0.0, 1.0 - self.eta * dt_f)
        if isinstance(state, DefectState):
            if state.omega is not None:
                state.omega = np.asarray(state.omega, dtype=float) * factor
            state.defects = self._relax_defects(state.defects, dt_f)
            state.observables["curvature_energy"] = self.energy(state)
            state.metadata["last_evolution"] = "curvature_relaxation"
            return state
        if _is_geometry(state):
            if state.omega is not None:
                state.omega = np.asarray(state.omega, dtype=float) * factor
            state.curvature_defects = self._relax_defects(state.curvature_defects, dt_f)
            # shrink scalar/matrix connection deviations toward identity/flatness
            new_conn: Dict[Tuple[int, int], np.ndarray] = {}
            for e, A in state.connections.items():
                arr = np.asarray(A, dtype=float)
                if arr.ndim == 2 and arr.shape[0] == arr.shape[1]:
                    I = np.eye(arr.shape[0])
                    new_conn[e] = I + factor * (arr - I)
                else:
                    new_conn[e] = factor * arr
            state.connections = new_conn
            state.metadata["last_evolution"] = "curvature_relaxation"
            state.metadata["curvature_energy"] = self.energy(state)
            return state
        return state


class DefectDiffusion(EvolutionRule):
    """Diffuse defect charges over the support graph.

    Defect charges are projected to node supports, evolved by an edge diffusion
    step, then converted back to point-supported defects.  If a charge vector is
    already present under ``state.observables['defect_charges']`` it is used
    directly.
    """

    name = "defect_diffusion"

    def __init__(self, config: Optional[EvolutionConfig] = None, *, diffusion: float = 0.1, threshold: float = 1e-9, **kwargs: Any):
        super().__init__(config, **kwargs)
        self.diffusion = float(diffusion)
        self.threshold = float(threshold)

    def step(self, state: StateLike, dt: Optional[float] = None) -> StateLike:
        dt_f = float(self.config.dt if dt is None else dt)
        if isinstance(state, DefectState):
            conn = state.connection
            n_hint = None
            if state.generators is not None:
                try:
                    n_hint = len(state.generators)
                except Exception:
                    n_hint = None
            ei, ej, w, n = _extract_edges(conn, n=n_hint)
            if n == 0:
                state.metadata["last_evolution"] = "defect_diffusion:no_graph"
                return state
            if "defect_charges" in state.observables:
                q = np.asarray(state.observables["defect_charges"], dtype=float).ravel()
                if len(q) != n:
                    q = np.resize(q, n).astype(float)
            else:
                q = defect_charges_to_vector(state.defects, n=n)
            if len(ei):
                if self.use_numba:
                    q2 = _edge_diffusion_step_numba(q, ei, ej, w, self.diffusion, dt_f)
                else:
                    q2 = q.copy()
                    for i, j, wij in zip(ei, ej, w):
                        dq = self.diffusion * dt_f * wij * (q[j] - q[i])
                        q2[i] += dq
                        q2[j] -= dq
                q = q2
            state.observables["defect_charges"] = q
            state.defects = vector_to_defects(q, threshold=self.threshold)
            state.metadata["last_evolution"] = "defect_diffusion"
            return state

        if _is_geometry(state):
            ei, ej, w, n = _extract_edges(state.connections, edges=state.edges, n=len(state.graph_nodes))
            q = defect_charges_to_vector(state.curvature_defects, n=n)
            if len(ei):
                if self.use_numba:
                    q = _edge_diffusion_step_numba(q, ei, ej, w, self.diffusion, dt_f)
                else:
                    q2 = q.copy()
                    for i, j, wij in zip(ei, ej, w):
                        dq = self.diffusion * dt_f * wij * (q[j] - q[i])
                        q2[i] += dq
                        q2[j] -= dq
                    q = q2
            state.boundary_data["defect_charges"] = q
            state.curvature_defects = vector_to_defects(q, kind="diffused_curvature", threshold=self.threshold)
            state.metadata["last_evolution"] = "defect_diffusion"
            return state
        return state


class ReactionDiffusionEvolution(DefectDiffusion):
    """Defect diffusion plus a stable local reaction term.

    The default reaction is cubic saturation:

        dq/dt = growth*q - saturation*q^3.

    A custom reaction callable can be supplied.  It should accept ``q`` and
    return ``dq/dt`` with the same shape.
    """

    name = "reaction_diffusion"

    def __init__(
        self,
        config: Optional[EvolutionConfig] = None,
        *,
        diffusion: float = 0.1,
        growth: float = 0.0,
        saturation: float = 1.0,
        reaction: Optional[Callable[[np.ndarray], np.ndarray]] = None,
        threshold: float = 1e-9,
        **kwargs: Any,
    ):
        super().__init__(config, diffusion=diffusion, threshold=threshold, **kwargs)
        self.growth = float(growth)
        self.saturation = float(saturation)
        self.reaction = reaction

    def step(self, state: StateLike, dt: Optional[float] = None) -> StateLike:
        dt_f = float(self.config.dt if dt is None else dt)
        state = super().step(state, dt=dt_f)
        if isinstance(state, DefectState):
            q = np.asarray(state.observables.get("defect_charges", []), dtype=float).ravel()
            if q.size:
                if self.reaction is not None:
                    q = q + dt_f * np.asarray(self.reaction(q), dtype=float)
                elif self.use_numba:
                    q = _reaction_clip_numba(q, self.growth, self.saturation, dt_f)
                else:
                    q = q + dt_f * (self.growth * q - self.saturation * q ** 3)
                state.observables["defect_charges"] = q
                state.defects = vector_to_defects(q, threshold=self.threshold)
                state.metadata["last_evolution"] = "reaction_diffusion"
        elif _is_geometry(state):
            q = np.asarray(state.boundary_data.get("defect_charges", []), dtype=float).ravel()
            if q.size:
                if self.reaction is not None:
                    q = q + dt_f * np.asarray(self.reaction(q), dtype=float)
                elif self.use_numba:
                    q = _reaction_clip_numba(q, self.growth, self.saturation, dt_f)
                else:
                    q = q + dt_f * (self.growth * q - self.saturation * q ** 3)
                state.boundary_data["defect_charges"] = q
                state.curvature_defects = vector_to_defects(q, kind="reaction_curvature", threshold=self.threshold)
                state.metadata["last_evolution"] = "reaction_diffusion"
        return state


class HamiltonianEvolution(EvolutionRule):
    """Small/medium dense Hamiltonian evolution for state-like vectors.

    This is intentionally simple and meant as an integration point for MFVM-like
    backends.  It supports explicit Euler, exact eigendecomposition for small
    dense Hermitian matrices, and imaginary-time Euler with normalization.
    """

    name = "hamiltonian"

    def __init__(
        self,
        config: Optional[EvolutionConfig] = None,
        *,
        H: Optional[np.ndarray] = None,
        imaginary_time: bool = False,
        method: str = "auto",  # auto | euler | exact
        normalize: bool = True,
        max_exact_dim: int = 512,
        **kwargs: Any,
    ):
        super().__init__(config, **kwargs)
        self.H = None if H is None else np.asarray(H)
        self.imaginary_time = bool(imaginary_time)
        self.method = method
        self.normalize = bool(normalize)
        self.max_exact_dim = int(max_exact_dim)

    def _evolve_vector(self, psi: np.ndarray, dt_f: float) -> np.ndarray:
        if self.H is None:
            return psi
        H = np.asarray(self.H)
        method = self.method
        if method == "auto":
            method = "exact" if H.ndim == 2 and H.shape[0] <= self.max_exact_dim else "euler"
        if method == "exact" and H.ndim == 2 and H.shape[0] == H.shape[1]:
            vals, vecs = np.linalg.eigh(H.astype(complex))
            coeff = vecs.conj().T @ psi
            if self.imaginary_time:
                phase = np.exp(-dt_f * vals)
            else:
                phase = np.exp(-1j * dt_f * vals)
            out = vecs @ (phase * coeff)
        else:
            if self.imaginary_time:
                out = psi - dt_f * (H @ psi)
            else:
                out = psi - 1j * dt_f * (H @ psi)
        if self.normalize:
            norm = np.linalg.norm(out)
            if norm > EPS:
                out = out / norm
        return out

    def step(self, state: StateLike, dt: Optional[float] = None) -> StateLike:
        dt_f = float(self.config.dt if dt is None else dt)
        if isinstance(state, DefectState):
            psi = state.observables.get("state_vector", None)
            if psi is None and state.generators is not None:
                arr = np.asarray(state.generators)
                if arr.ndim == 1 and np.iscomplexobj(arr):
                    psi = arr
            if psi is None:
                state.metadata["last_evolution"] = "hamiltonian:no_state_vector"
                return state
            out = self._evolve_vector(np.asarray(psi, dtype=complex).ravel(), dt_f)
            state.observables["state_vector"] = out
            state.generators = out
            state.metadata["last_evolution"] = "hamiltonian"
            state.metadata["imaginary_time"] = self.imaginary_time
            return state
        if _is_geometry(state):
            psi = state.boundary_data.get("state_vector") if hasattr(state, "boundary_data") else None
            if psi is not None:
                state.boundary_data["state_vector"] = self._evolve_vector(np.asarray(psi, dtype=complex).ravel(), dt_f)
                state.metadata["last_evolution"] = "hamiltonian"
            return state
        return state


# ---------------------------------------------------------------------------
# Pipeline and postprocessing
# ---------------------------------------------------------------------------


def _postprocess_evolved(state: StateLike, cfg: EvolutionConfig) -> StateLike:
    """Apply optional gauge fixing and coarse graining after one step."""
    out = state
    if cfg.gauge_fix and _is_geometry(out) and gauge_fix_geometry is not None:
        try:
            out = gauge_fix_geometry(out)
        except Exception:
            out.metadata["gauge_fix_failed"] = True
    if cfg.coarse_grain and isinstance(out, DefectState) and coarse_grain_state is not None and CoarseGrainConfig is not None:
        try:
            cg_cfg = CoarseGrainConfig(merge_radius=cfg.merge_radius, omega_rank=cfg.coarse_omega_rank)
            out, report = coarse_grain_state(out, cg_cfg)
            out.metadata["last_coarse_grain_report"] = report.__dict__
        except Exception:
            out.metadata["coarse_grain_failed"] = True
    if cfg.coarse_grain and _is_geometry(out):
        # Geometry-level lightweight coarse step: drop tiny curvature defects.
        try:
            if hasattr(out, "compress_curvature"):
                out.compress_curvature()
        except Exception:
            pass
    return out


class EvolutionPipeline:
    """Composable evolution pipeline.

    A pipeline can mix transport, relaxation, diffusion, Hamiltonian steps, and
    user-supplied rules.  Each global step applies all rules in order, then the
    shared postprocessor.
    """

    def __init__(self, rules: Sequence[Union[EvolutionRule, Callable[[StateLike, float], StateLike]]], config: Optional[EvolutionConfig] = None):
        self.rules = list(rules)
        self.config = config or EvolutionConfig()

    def step(self, state: StateLike, dt: Optional[float] = None) -> StateLike:
        dt_f = float(self.config.dt if dt is None else dt)
        cur = state
        for rule in self.rules:
            if isinstance(rule, EvolutionRule):
                cur = rule.step(cur, dt=dt_f)
            else:
                cur = rule(cur, dt_f)
        return _postprocess_evolved(cur, self.config)

    def evolve(self, state: StateLike, *, steps: Optional[int] = None, dt: Optional[float] = None) -> Tuple[StateLike, EvolutionReport]:
        cfg = self.config
        steps_i = int(cfg.steps if steps is None else steps)
        dt_f = float(cfg.dt if dt is None else dt)
        cur = _clone_state_like(state) if cfg.copy else state
        report = EvolutionReport(rule="pipeline", steps=steps_i, dt=dt_f, metadata={"n_rules": len(self.rules), "mode": cfg.resolved_mode()})
        if cfg.track_history:
            report.energies.append(_state_energy(cur))
            report.n_defects.append(_n_defects(cur))
        for _ in range(steps_i):
            cur = self.step(cur, dt=dt_f)
            if cfg.track_history:
                report.energies.append(_state_energy(cur))
                report.n_defects.append(_n_defects(cur))
        return cur, report


# ---------------------------------------------------------------------------
# Convenience functions and VM integration
# ---------------------------------------------------------------------------


def evolve_state(state: StateLike, rule: EvolutionRule, *, steps: int = 1, dt: float = 1.0) -> Tuple[StateLike, EvolutionReport]:
    """Evolve ``state`` under a single rule."""
    return rule.evolve(state, steps=steps, dt=dt)


def transport_along_connection(state: StateLike, *, steps: int = 1, dt: float = 1.0, damping: float = 0.0, connection: Any = None) -> Tuple[StateLike, EvolutionReport]:
    rule = TransportEvolution(EvolutionConfig(steps=steps, dt=dt), damping=damping, connection=connection)
    return rule.evolve(state, steps=steps, dt=dt)


def relax_curvature(state: StateLike, *, steps: int = 1, dt: float = 1.0, eta: float = 0.1, defect_decay: float = 0.05) -> Tuple[StateLike, EvolutionReport]:
    rule = CurvatureRelaxation(EvolutionConfig(steps=steps, dt=dt), eta=eta, defect_decay=defect_decay)
    return rule.evolve(state, steps=steps, dt=dt)


def diffuse_defects(state: StateLike, *, steps: int = 1, dt: float = 1.0, diffusion: float = 0.1) -> Tuple[StateLike, EvolutionReport]:
    rule = DefectDiffusion(EvolutionConfig(steps=steps, dt=dt), diffusion=diffusion)
    return rule.evolve(state, steps=steps, dt=dt)


def evolve_and_recompress(
    state: StateLike,
    rule: Union[EvolutionRule, EvolutionPipeline],
    *,
    steps: int = 1,
    dt: float = 1.0,
    gauge_fix: bool = True,
    coarse_grain: bool = True,
    coarse_omega_rank: Optional[int] = None,
) -> Tuple[StateLike, EvolutionReport]:
    """Evolve with post-step gauge fixing and recompression enabled."""
    cfg = EvolutionConfig(steps=steps, dt=dt, gauge_fix=gauge_fix, coarse_grain=coarse_grain, coarse_omega_rank=coarse_omega_rank)
    if isinstance(rule, EvolutionPipeline):
        rule.config = cfg
        return rule.evolve(state, steps=steps, dt=dt)
    rule.config.gauge_fix = gauge_fix
    rule.config.coarse_grain = coarse_grain
    rule.config.coarse_omega_rank = coarse_omega_rank
    return rule.evolve(state, steps=steps, dt=dt)


def register_vm_ops(vm: Any) -> Any:
    """Register evolution operations on a DefectVM-like object."""

    def _update_vm_state(result: Tuple[StateLike, EvolutionReport]) -> Dict[str, Any]:
        new_state, report = result
        # DefectVM stores state by attribute.  Replace it when possible.
        try:
            vm.state = new_state
            vm.context["state"] = new_state
        except Exception:
            pass
        return {"state": new_state, "report": report.as_dict()}

    vm.register("evolve_transport", lambda **kwargs: _update_vm_state(transport_along_connection(vm.state, **kwargs)))
    vm.register("evolve_relax_curvature", lambda **kwargs: _update_vm_state(relax_curvature(vm.state, **kwargs)))
    vm.register("evolve_diffuse_defects", lambda **kwargs: _update_vm_state(diffuse_defects(vm.state, **kwargs)))
    vm.register("evolve", lambda rule, **kwargs: _update_vm_state(evolve_state(vm.state, rule, **kwargs)))
    return vm


__all__ = [
    "EvolutionConfig",
    "EvolutionReport",
    "EvolutionRule",
    "TransportEvolution",
    "CurvatureRelaxation",
    "DefectDiffusion",
    "ReactionDiffusionEvolution",
    "HamiltonianEvolution",
    "EvolutionPipeline",
    "graph_laplacian_from_edges",
    "defect_charges_to_vector",
    "vector_to_defects",
    "evolve_state",
    "transport_along_connection",
    "relax_curvature",
    "diffuse_defects",
    "evolve_and_recompress",
    "register_vm_ops",
]
