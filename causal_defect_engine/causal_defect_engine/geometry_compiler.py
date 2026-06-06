"""
geometry_compiler.py
====================

Advanced curvature-sparse defect-geometry compiler.

This module turns the v0.6 representation into a more explicit compiler/runtime
layer.  It can learn an adaptive topology from Omega/features, gauge-fix the
learned connection so a flat spanning-tree background is removed, store only
residual curvature/gluing defects, and compile to MPS/TT-style targets when a
small dense enumeration is acceptable.

The implementation is intentionally dependency-light and conservative:
- exact/legacy modes keep simple threshold graph construction,
- efficient/auto modes use sparse topology construction and optional Numba
  kernels for dense scans, connection construction, gauge residuals, and Ising
  batch helpers,
- all paths keep NumPy/Python fallbacks.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union, Callable

import numpy as np

from .core import Defect, DefectState, normalize_mode, NUMBA_AVAILABLE, EPS
from .tensor_defect import SimpleMPS, dense_to_mps, extract_omega, learn_mps_ordering, estimate_bond_dimensions
from .state_geometry import (
    StateGeometryConfig,
    CurvatureSparseStateGeometry,
    geometry_from_omega,
    geometry_from_mps,
    geometry_from_ising,
)
from .numba_kernels import (
    threshold_edges_numba,
    edge_connection_binary_numba,
    triangle_curvature_candidates_numba,
    tree_gauge_residual_numba,
)


@dataclass
class CompilerConfig(StateGeometryConfig):
    """Configuration for advanced defect-geometry compilation."""

    topology_method: str = "auto"  # auto | threshold | mst | mst_chords | knn | legacy
    max_edges: Optional[int] = None
    k_neighbors: int = 4
    chord_budget: Optional[int] = None
    gauge_group: str = "u1"  # u1 in this prototype
    gauge_fix: bool = True
    store_flat_tree: bool = False
    phase_from: str = "omega"  # omega | sign | weight
    compile_target: str = "auto"  # auto | geometry | mps
    max_treewidth_proxy: Optional[int] = None
    error_metric: str = "storage"  # storage | amplitude | observable placeholder
    confidence_samples: int = 0


@dataclass
class GeometryCompileReport:
    """Compiler diagnostics separating store/eval/sample/learn costs."""

    storage_complexity: int
    eval_complexity: int
    sample_complexity: Optional[int]
    learn_complexity_estimate: Optional[int]
    n_nodes: int
    n_edges: int
    n_tree_edges: int
    n_residual_edges: int
    n_curvature_defects: int
    curvature_sparsity: float
    omega_rank: Optional[int]
    effective_bond_rank: Optional[int]
    compression_ratio_vs_dense: Optional[float]
    mode: str
    backend: str
    notes: Dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "storage_complexity": self.storage_complexity,
            "eval_complexity": self.eval_complexity,
            "sample_complexity": self.sample_complexity,
            "learn_complexity_estimate": self.learn_complexity_estimate,
            "n_nodes": self.n_nodes,
            "n_edges": self.n_edges,
            "n_tree_edges": self.n_tree_edges,
            "n_residual_edges": self.n_residual_edges,
            "n_curvature_defects": self.n_curvature_defects,
            "curvature_sparsity": self.curvature_sparsity,
            "omega_rank": self.omega_rank,
            "effective_bond_rank": self.effective_bond_rank,
            "compression_ratio_vs_dense": self.compression_ratio_vs_dense,
            "mode": self.mode,
            "backend": self.backend,
            "notes": dict(self.notes),
        }


# ---------------------------------------------------------------------------
# Topology learning
# ---------------------------------------------------------------------------


def learn_adaptive_topology(
    omega: np.ndarray,
    *,
    config: Optional[CompilerConfig] = None,
    method: Optional[str] = None,
) -> Tuple[List[Tuple[int, int]], Dict[str, Any]]:
    """Learn a sparse gluing graph from an Omega/correlation matrix.

    Efficient default: maximum spanning forest backbone plus a bounded set of
    strongest residual chords.  This captures the minimum connected skeleton and
    then stores non-tree residual couplings as possible curvature sources.
    """
    cfg = config or CompilerConfig()
    mode = normalize_mode(cfg.mode)
    O = np.asarray(omega, dtype=float)
    if O.ndim != 2 or O.shape[0] != O.shape[1]:
        raise ValueError("omega must be square")
    n = O.shape[0]
    W = np.abs(O.copy())
    np.fill_diagonal(W, 0.0)
    meth = normalize_mode(method or cfg.topology_method)
    if meth == "auto":
        meth = "mst_chords" if mode in {"auto", "efficient", "local"} else "threshold"
    if meth == "legacy":
        meth = "threshold"

    max_edges = cfg.max_edges
    if max_edges is None:
        # Sparse default: connected tree + O(kn) chords.
        max_edges = max(n - 1, (n * max(int(cfg.k_neighbors), 1)) // 2)

    if meth == "threshold":
        if cfg.use_numba and NUMBA_AVAILABLE and n * n > 256:
            raw = threshold_edges_numba(W, float(cfg.connection_threshold), int(max_edges))
            edges = [(int(i), int(j)) for i, j in raw if int(i) >= 0]
        else:
            edges = [(i, j) for i in range(n) for j in range(i + 1, n) if W[i, j] >= cfg.connection_threshold]
            edges = sorted(edges, key=lambda e: W[e[0], e[1]], reverse=True)[: int(max_edges)]
        return edges, {"topology_method": meth, "max_edges": max_edges}

    if meth in {"mst", "mst_chords"}:
        tree = _maximum_spanning_tree(W)
        used = set(tuple(sorted(e)) for e in tree)
        edges = list(tree)
        if meth == "mst_chords":
            budget = cfg.chord_budget
            if budget is None:
                budget = max(0, int(max_edges) - len(edges))
            candidates = []
            for i in range(n):
                for j in range(i + 1, n):
                    if (i, j) not in used and W[i, j] >= cfg.connection_threshold:
                        candidates.append((float(W[i, j]), i, j))
            candidates.sort(reverse=True)
            for _, i, j in candidates[: int(budget)]:
                edges.append((i, j))
        return edges, {"topology_method": meth, "max_edges": max_edges, "tree_edges": len(tree)}

    if meth == "knn":
        k = max(1, int(cfg.k_neighbors))
        edge_set = set()
        for i in range(n):
            idx = np.argsort(-W[i])[:k]
            for j in idx:
                if i != int(j) and W[i, int(j)] >= cfg.connection_threshold:
                    edge_set.add(tuple(sorted((i, int(j)))))
        edges = sorted(edge_set, key=lambda e: W[e[0], e[1]], reverse=True)[: int(max_edges)]
        return edges, {"topology_method": meth, "k_neighbors": k, "max_edges": max_edges}

    raise ValueError(f"unknown topology method {method!r}")


def _maximum_spanning_tree(W: np.ndarray) -> List[Tuple[int, int]]:
    n = W.shape[0]
    if n <= 1:
        return []
    parent = list(range(n))
    rank = [0] * n

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a: int, b: int) -> bool:
        ra, rb = find(a), find(b)
        if ra == rb:
            return False
        if rank[ra] < rank[rb]:
            parent[ra] = rb
        elif rank[ra] > rank[rb]:
            parent[rb] = ra
        else:
            parent[rb] = ra
            rank[ra] += 1
        return True

    candidates = [(float(W[i, j]), i, j) for i in range(n) for j in range(i + 1, n) if W[i, j] > 0]
    candidates.sort(reverse=True)
    edges: List[Tuple[int, int]] = []
    for _, i, j in candidates:
        if union(i, j):
            edges.append((i, j))
            if len(edges) == n - 1:
                break
    return edges


def build_binary_connections(omega: np.ndarray, edges: Sequence[Tuple[int, int]], *, use_numba: bool = True) -> Dict[Tuple[int, int], np.ndarray]:
    O = np.asarray(omega, dtype=float)
    edge_arr = np.asarray(edges, dtype=np.int64)
    conns: Dict[Tuple[int, int], np.ndarray] = {}
    if len(edges) and use_numba and NUMBA_AVAILABLE:
        mats = edge_connection_binary_numba(np.abs(O), edge_arr)
        for idx, (i, j) in enumerate(edges):
            A = np.asarray(mats[idx], dtype=float)
            conns[(int(i), int(j))] = A
            conns[(int(j), int(i))] = A.T
    else:
        for i, j in edges:
            w = abs(float(O[i, j]))
            A = np.array([[1.0 + w, 1.0 - w], [1.0 - w, 1.0 + w]], dtype=float)
            conns[(int(i), int(j))] = A
            conns[(int(j), int(i))] = A.T
    return conns


def detect_curvature_loops(
    omega: np.ndarray,
    edges: Sequence[Tuple[int, int]],
    *,
    tol: float = 0.05,
    max_defects: Optional[int] = 512,
    use_numba: bool = True,
) -> List[Defect]:
    """Detect sparse triangle curvature defects in the learned topology."""
    W = np.abs(np.asarray(omega, dtype=float))
    edge_arr = np.asarray(edges, dtype=np.int64)
    cap = int(max_defects if max_defects is not None else max(1024, len(edges) * 4 + 1))
    defects: List[Defect] = []
    if len(edges) and use_numba and NUMBA_AVAILABLE and W.shape[0] > 8:
        loops, charges = triangle_curvature_candidates_numba(W, edge_arr, float(tol), cap)
        for loop, charge in zip(loops, charges):
            if int(loop[0]) < 0:
                continue
            c = float(charge)
            defects.append(Defect("curvature_loop", charge=c, support=tuple(map(int, loop)), strength=abs(c), metadata={"loop": tuple(map(int, loop)), "detector": "numba_triangle"}))
    else:
        adj = np.zeros_like(W, dtype=bool)
        for i, j in edges:
            adj[int(i), int(j)] = True; adj[int(j), int(i)] = True
        n = W.shape[0]
        for i in range(n):
            nbrs = np.flatnonzero(adj[i])
            for a, j in enumerate(nbrs):
                j = int(j)
                if j <= i:
                    continue
                for k in nbrs[a + 1:]:
                    k = int(k)
                    if k <= j or not adj[j, k]:
                        continue
                    c = float(W[i, j] * W[j, k] * W[k, i])
                    if abs(c) > tol:
                        defects.append(Defect("curvature_loop", charge=c, support=(i, j, k), strength=abs(c), metadata={"loop": (i, j, k), "detector": "python_triangle"}))
    defects.sort(key=lambda d: abs(float(d.charge)) * float(d.strength), reverse=True)
    if max_defects is not None:
        defects = defects[: int(max_defects)]
    return defects


# ---------------------------------------------------------------------------
# Gauge fixing
# ---------------------------------------------------------------------------


def edge_phases_from_omega(omega: np.ndarray, edges: Sequence[Tuple[int, int]], *, phase_from: str = "omega") -> np.ndarray:
    O = np.asarray(omega)
    phases = np.zeros(len(edges), dtype=float)
    for k, (i, j) in enumerate(edges):
        val = O[int(i), int(j)]
        if np.iscomplexobj(O):
            phases[k] = float(np.angle(val))
        elif phase_from == "sign":
            phases[k] = 0.0 if float(val) >= 0 else np.pi
        elif phase_from == "weight":
            phases[k] = float(val)
        else:
            # real-valued Omega has no phase, but sign carries a useful Z2/U1 proxy.
            phases[k] = 0.0 if float(val) >= 0 else np.pi
    return (phases + np.pi) % (2.0 * np.pi) - np.pi


def gauge_fix_geometry(geometry: CurvatureSparseStateGeometry, *, config: Optional[CompilerConfig] = None) -> CurvatureSparseStateGeometry:
    """Gauge-fix a U(1)-proxy connection along a spanning forest.

    Tree residual phases are flattened to zero; non-tree residual phases are
    stored as ``residual_holonomy`` defects when nontrivial.  If
    ``store_flat_tree=False``, tree edges are removed from explicit connection
    storage metadata but retained in the graph for topology/evaluation.
    """
    cfg = config or CompilerConfig()
    if not geometry.edges:
        geometry.metadata["gauge_fixed"] = True
        geometry.metadata["tree_edges"] = []
        geometry.metadata["residual_edges"] = []
        return geometry
    n = geometry.n_sites
    edges = list((int(i), int(j)) for i, j in geometry.edges)
    edge_arr = np.asarray(edges, dtype=np.int64)
    omega = geometry.omega if geometry.omega is not None else np.zeros((n, n), dtype=float)
    phases = edge_phases_from_omega(np.asarray(omega), edges, phase_from=cfg.phase_from)
    if cfg.use_numba and NUMBA_AVAILABLE:
        gauge, residual, tree_mask = tree_gauge_residual_numba(edge_arr, phases, n)
    else:
        gauge, residual, tree_mask = _tree_gauge_residual_python(edge_arr, phases, n)
    residual_edges: List[Tuple[int, int]] = []
    tree_edges: List[Tuple[int, int]] = []
    for k, e in enumerate(edges):
        if int(tree_mask[k]) == 1:
            tree_edges.append(e)
        else:
            residual_edges.append(e)
            if abs(float(residual[k])) > cfg.curvature_tol:
                geometry.curvature_defects.append(Defect(
                    "residual_holonomy",
                    charge=float(residual[k]),
                    support=e,
                    strength=abs(float(residual[k])),
                    metadata={"edge": e, "phase": float(phases[k]), "residual_phase": float(residual[k]), "gauge": cfg.gauge_group},
                ))
    geometry.metadata.update({
        "gauge_fixed": True,
        "gauge_group": cfg.gauge_group,
        "node_gauge": np.asarray(gauge).tolist(),
        "tree_edges": tree_edges,
        "residual_edges": residual_edges,
        "n_tree_edges": len(tree_edges),
        "n_residual_edges": len(residual_edges),
        "stored_flat_tree": bool(cfg.store_flat_tree),
    })
    geometry.compress_curvature(tol=cfg.curvature_tol, max_defects=cfg.max_curvature_defects)
    return geometry


def _tree_gauge_residual_python(edges: np.ndarray, phases: np.ndarray, n_nodes: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    adj: List[List[Tuple[int, int, float]]] = [[] for _ in range(n_nodes)]
    for e, (i, j) in enumerate(edges):
        adj[int(i)].append((int(j), e, 1.0))
        adj[int(j)].append((int(i), e, -1.0))
    gauge = np.zeros(n_nodes, dtype=float)
    seen = np.zeros(n_nodes, dtype=bool)
    tree = np.zeros(len(edges), dtype=np.uint8)
    from collections import deque
    for root in range(n_nodes):
        if seen[root]:
            continue
        seen[root] = True
        q = deque([root])
        while q:
            u = q.popleft()
            for v, e, sgn in adj[u]:
                if not seen[v]:
                    gauge[v] = gauge[u] + phases[e] * sgn
                    seen[v] = True
                    tree[e] = 1
                    q.append(v)
    residual = np.zeros(len(edges), dtype=float)
    for e, (i, j) in enumerate(edges):
        residual[e] = (phases[e] + gauge[int(i)] - gauge[int(j)] + np.pi) % (2 * np.pi) - np.pi
    return gauge, residual, tree


# ---------------------------------------------------------------------------
# Compiler and report
# ---------------------------------------------------------------------------


class DefectGeometryCompiler:
    """Advanced compiler for curvature-sparse state geometry."""

    def __init__(self, config: Optional[CompilerConfig] = None, **kwargs: Any):
        self.config = config or CompilerConfig(**{k: v for k, v in kwargs.items() if k in CompilerConfig.__annotations__})
        self.geometry: Optional[CurvatureSparseStateGeometry] = None
        self.report_: Optional[GeometryCompileReport] = None

    def fit(self, source: Any, **kwargs: Any) -> CurvatureSparseStateGeometry:
        """Fit from Omega, samples/features, MPS, dense tensor, Ising spec, or geometry."""
        cfg = self.config
        if isinstance(source, CurvatureSparseStateGeometry):
            geom = source
        elif isinstance(source, SimpleMPS):
            geom = geometry_from_mps(source, config=cfg)
        elif isinstance(source, Mapping) and str(source.get("type", "")).lower() in {"ising", "boltzmann", "phase_ising"}:
            geom = geometry_from_ising(source, config=cfg)
        elif isinstance(source, np.ndarray) and source.ndim >= 3:
            geom = geometry_from_mps(dense_to_mps(np.asarray(source), max_bond=cfg.max_bond, tol=cfg.svd_tol), config=cfg)
            geom.metadata["source"] = "dense_tensor"
        else:
            if isinstance(source, np.ndarray) and source.ndim == 2 and source.shape[0] == source.shape[1] and np.allclose(source, source.T.conj(), atol=1e-10):
                omega = np.asarray(source)
            else:
                omega = extract_omega(source, use_numba=cfg.use_numba)
            geom = self.fit_from_omega(omega, _return=True)
        if not (isinstance(source, np.ndarray) and source.ndim == 2 and source.shape[0] == source.shape[1]) and not isinstance(source, (CurvatureSparseStateGeometry, SimpleMPS)):
            # Generator/MPS paths still benefit from report/gauge metadata but do not need topology refit.
            if cfg.gauge_fix and geom.omega is not None and geom.edges:
                gauge_fix_geometry(geom, config=cfg)
            self.geometry = geom
            self.report_ = complexity_report(geom, config=cfg)
            return geom
        if isinstance(source, CurvatureSparseStateGeometry):
            if cfg.gauge_fix and geom.edges:
                gauge_fix_geometry(geom, config=cfg)
            self.geometry = geom
            self.report_ = complexity_report(geom, config=cfg)
            return geom
        self.geometry = geom
        self.report_ = complexity_report(geom, config=cfg)
        return geom

    def fit_from_omega(self, omega: np.ndarray, *, _return: bool = False) -> CurvatureSparseStateGeometry:
        cfg = self.config
        O = np.asarray(omega)
        edges, topo_meta = learn_adaptive_topology(O, config=cfg)
        conns = build_binary_connections(np.asarray(O, dtype=float), edges, use_numba=cfg.use_numba)
        defects = detect_curvature_loops(np.asarray(O, dtype=float), edges, tol=cfg.curvature_tol, max_defects=cfg.max_curvature_defects, use_numba=cfg.use_numba)
        s = np.linalg.svd(np.asarray(O, dtype=float), compute_uv=False)
        omega_rank = int(np.sum(s > max(cfg.svd_tol, cfg.svd_tol * float(s[0] if len(s) else 1.0)))) if len(s) else 0
        ordering = learn_mps_ordering(np.abs(O), method="greedy")
        bonds = estimate_bond_dimensions(np.asarray(O, dtype=float), ordering, tol=cfg.svd_tol, max_bond=cfg.max_bond)
        geom = CurvatureSparseStateGeometry(
            graph_nodes=list(range(O.shape[0])),
            edges=edges,
            local_dims=[int(cfg.local_dim)] * O.shape[0],
            connections=conns,
            curvature_defects=defects,
            omega=np.asarray(O, dtype=float),
            eval_rule="holonomy_product" if cfg.eval_rule == "auto" else cfg.eval_rule,
            metadata={
                "source": "omega_compiler",
                "compiler": "DefectGeometryCompiler",
                "omega_rank": omega_rank,
                "ordering": ordering,
                "bond_dims_est": bonds,
                "max_bond": int(max(bonds)) if bonds else 1,
                "max_dense_states": cfg.max_dense_states,
                **topo_meta,
            },
        )
        if cfg.gauge_fix:
            gauge_fix_geometry(geom, config=cfg)
        self.geometry = geom
        self.report_ = complexity_report(geom, config=cfg)
        return geom

    def fit_samples(self, samples: np.ndarray) -> CurvatureSparseStateGeometry:
        return self.fit(samples)

    def fit_generator(self, generator: Mapping[str, Any]) -> CurvatureSparseStateGeometry:
        return self.fit(generator)

    def fit_circuit(self, circuit: Mapping[str, Any]) -> CurvatureSparseStateGeometry:
        """Prototype compact circuit/generator storage.

        This stores the circuit as generator-compressed boundary data.  It does
        not simulate arbitrary quantum circuits; it honestly reports that eval
        cost depends on the supplied evaluator if any.
        """
        n = int(circuit.get("n", circuit.get("n_qubits", 0)))
        if n <= 0:
            raise ValueError("circuit spec requires n or n_qubits")
        gates = list(circuit.get("gates", []))
        geom = CurvatureSparseStateGeometry(
            graph_nodes=list(range(n)),
            edges=[],
            local_dims=[2] * n,
            boundary_data={"circuit": dict(circuit), "gates": gates, "evaluator": circuit.get("evaluator")},
            eval_rule="generator",
            metadata={"source": "circuit_generator", "generator_storage": len(gates), "note": "amplitude eval requires a user-supplied evaluator"},
        )
        self.geometry = geom
        self.report_ = complexity_report(geom, config=self.config)
        return geom

    def compile(self, target: Optional[str] = None, **kwargs: Any) -> Any:
        if self.geometry is None:
            raise RuntimeError("call fit() first")
        return compile_geometry(self.geometry, target=target or self.config.compile_target, config=self.config, **kwargs)

    def report(self) -> GeometryCompileReport:
        if self.geometry is None:
            raise RuntimeError("call fit() first")
        self.report_ = complexity_report(self.geometry, config=self.config)
        return self.report_

    def to_defect_state(self) -> DefectState:
        if self.geometry is None:
            raise RuntimeError("call fit() first")
        state = self.geometry.to_defect_state()
        state.observables["compiler_report"] = self.report().as_dict()
        return state


def compile_geometry(geometry: CurvatureSparseStateGeometry, *, target: str = "auto", config: Optional[CompilerConfig] = None, **kwargs: Any) -> Any:
    cfg = config or CompilerConfig()
    t = str(target or "auto").lower()
    if t == "auto":
        t = "mps" if geometry.n_sites <= 24 else "geometry"
    if t in {"geometry", "defect", "state_geometry"}:
        return geometry
    if t in {"mps", "tt", "tensor_train"}:
        return geometry.to_mps(max_bond=kwargs.get("max_bond", cfg.max_bond), tol=kwargs.get("tol", cfg.svd_tol), max_states=kwargs.get("max_states", cfg.max_dense_states))
    if t in {"state", "defect_state"}:
        return geometry.to_defect_state()
    raise ValueError(f"unknown compile target {target!r}")


def complexity_report(geometry: CurvatureSparseStateGeometry, *, config: Optional[CompilerConfig] = None) -> GeometryCompileReport:
    cfg = config or CompilerConfig()
    c = geometry.estimate_complexity()
    n_tree = int(geometry.metadata.get("n_tree_edges", max(0, geometry.n_sites - 1 if geometry.edges else 0)))
    n_res = int(geometry.metadata.get("n_residual_edges", max(0, len(geometry.edges) - n_tree)))
    curvature_sparsity = float(len(geometry.curvature_defects)) / float(max(1, len(geometry.edges)))
    sample_cost = c.eval_cost_estimate * int(geometry.n_sites) if geometry.n_sites <= 30 else None
    return GeometryCompileReport(
        storage_complexity=int(c.storage_parameters),
        eval_complexity=int(c.eval_cost_estimate),
        sample_complexity=sample_cost,
        learn_complexity_estimate=c.learn_cost_estimate or int(max(1, geometry.n_sites) ** 2),
        n_nodes=geometry.n_sites,
        n_edges=len(geometry.edges),
        n_tree_edges=n_tree,
        n_residual_edges=n_res,
        n_curvature_defects=len(geometry.curvature_defects),
        curvature_sparsity=curvature_sparsity,
        omega_rank=c.omega_rank,
        effective_bond_rank=c.max_bond,
        compression_ratio_vs_dense=c.compression_ratio_vs_dense,
        mode=normalize_mode(cfg.mode),
        backend=str(geometry.metadata.get("source", "geometry")),
        notes={"gauge_fixed": bool(geometry.metadata.get("gauge_fixed", False)), "eval_rule": geometry.eval_rule, "topology_method": geometry.metadata.get("topology_method")},
    )


# Convenience monkey-patched methods for v0.6 class users.  This keeps the public
# object ergonomic without forcing a large rewrite of state_geometry.py.

def _geom_fix_gauge(self: CurvatureSparseStateGeometry, **kwargs: Any) -> CurvatureSparseStateGeometry:
    cfg = CompilerConfig(**{k: v for k, v in kwargs.items() if k in CompilerConfig.__annotations__}) if kwargs else CompilerConfig()
    return gauge_fix_geometry(self, config=cfg)


def _geom_complexity_report(self: CurvatureSparseStateGeometry, **kwargs: Any) -> Dict[str, Any]:
    cfg = CompilerConfig(**{k: v for k, v in kwargs.items() if k in CompilerConfig.__annotations__}) if kwargs else CompilerConfig()
    return complexity_report(self, config=cfg).as_dict()


def _geom_compile_to_mps(self: CurvatureSparseStateGeometry, **kwargs: Any) -> SimpleMPS:
    return self.to_mps(**kwargs)


if not hasattr(CurvatureSparseStateGeometry, "fix_gauge"):
    setattr(CurvatureSparseStateGeometry, "fix_gauge", _geom_fix_gauge)
if not hasattr(CurvatureSparseStateGeometry, "complexity_report"):
    setattr(CurvatureSparseStateGeometry, "complexity_report", _geom_complexity_report)
if not hasattr(CurvatureSparseStateGeometry, "compile_to_mps"):
    setattr(CurvatureSparseStateGeometry, "compile_to_mps", _geom_compile_to_mps)


__all__ = [
    "CompilerConfig",
    "GeometryCompileReport",
    "DefectGeometryCompiler",
    "learn_adaptive_topology",
    "build_binary_connections",
    "detect_curvature_loops",
    "edge_phases_from_omega",
    "gauge_fix_geometry",
    "compile_geometry",
    "complexity_report",
]
