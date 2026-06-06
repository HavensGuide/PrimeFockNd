"""
graph_holonomy.py
=================

Graph holonomy backend for the Causal Defect Algebra Engine.

Selectable modes
----------------
legacy / exact
    DFS fundamental cycle basis, matching the v0.1 behavior.
efficient / local
    Spanning-tree chords with early cycle cap; avoids extra dedupe work.
auto
    Efficient mode for larger graphs, exact/legacy mode for small graphs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .core import CausalDefectEngine, Defect, DefectState, normalize_mode, wrap_angle
from .numba_kernels import graph_degree_imbalance_numba, cycle_holonomy_padded_numba


@dataclass
class GraphHolonomyConfig:
    cycle_basis_max: int = 5000
    defect_threshold: float = 0.25
    directed: bool = False
    mode: str = "auto"       # legacy | exact | efficient | local | auto
    max_exact_edges: int = 20000
    dense_phase_max_nodes: int = 2500  # enables Numba bulk holonomy when safe


class GraphHolonomyEngine(CausalDefectEngine):
    """Connection/holonomy defect engine for graphs."""

    def __init__(self, config: Optional[GraphHolonomyConfig] = None, **kwargs):
        config = config or GraphHolonomyConfig()
        kwargs.setdefault("mode", config.mode)
        super().__init__(**kwargs)
        self.config = config
        self.n_nodes: int = 0
        self.edges: np.ndarray = np.zeros((0, 2), dtype=np.int64)
        self.weights: np.ndarray = np.zeros(0, dtype=np.float64)
        self.phase: Dict[Tuple[int, int], float] = {}
        self.edge_phases: np.ndarray = np.zeros(0, dtype=np.float64)
        self._phase_dense: Optional[np.ndarray] = None
        self.adj: List[List[int]] = []
        self._resolved_mode = "legacy"

    def fit(self, events: Any, **kwargs) -> "GraphHolonomyEngine":
        phases = kwargs.get("phases", None)
        if isinstance(events, dict):
            edges = np.asarray(events["edges"], dtype=np.int64)
            weights = np.asarray(events.get("weights", np.ones(len(edges))), dtype=np.float64)
            n = int(events.get("n_nodes", edges.max() + 1 if len(edges) else 0))
            if phases is None:
                phases = events.get("phases", None)
        else:
            arr = np.asarray(events)
            if arr.ndim == 2 and arr.shape[0] == arr.shape[1]:
                src, dst = np.nonzero(arr)
                edges = np.column_stack([src, dst]).astype(np.int64)
                weights = arr[src, dst].astype(np.float64)
                n = arr.shape[0]
            else:
                edges = np.asarray(events, dtype=np.int64)
                weights = np.ones(len(edges), dtype=np.float64)
                n = int(kwargs.get("n_nodes", edges.max() + 1 if len(edges) else 0))

        self.n_nodes = n
        self.edges = edges
        self.weights = weights
        self._resolved_mode = self._choose_mode(len(edges))
        self._build_adjacency()
        self._install_phases(phases)

        self.state = DefectState(
            events={"n_nodes": n, "edges": edges, "weights": weights},
            metadata={
                "backend": "graph_holonomy",
                "numba": self.use_numba,
                "mode": self._resolved_mode,
                "n_nodes": n,
                "n_edges": len(edges),
            },
            valid=True,
        )
        return self

    def _choose_mode(self, e: int) -> str:
        mode = normalize_mode(self.config.mode or self.runtime.mode)
        if mode == "auto":
            return "exact" if e <= self.config.max_exact_edges else "efficient"
        if mode == "legacy":
            return "exact"
        if mode == "local":
            return "efficient"
        return mode

    def _build_adjacency(self) -> None:
        self.adj = [[] for _ in range(self.n_nodes)]
        seen = set()
        for i, j in self.edges:
            i, j = int(i), int(j)
            if i < 0 or j < 0 or i >= self.n_nodes or j >= self.n_nodes:
                continue
            self.adj[i].append(j)
            if not self.config.directed:
                self.adj[j].append(i)
                seen.add(tuple(sorted((i, j))))
            else:
                seen.add((i, j))

    def _install_phases(self, phases: Any) -> None:
        self.phase = {}
        if phases is None:
            scale = np.std(self.weights) + 1e-12
            for (i, j), w in zip(self.edges, self.weights):
                ph = float(w / scale)
                self.phase[(int(i), int(j))] = ph
                if not self.config.directed:
                    self.phase[(int(j), int(i))] = -ph
        elif isinstance(phases, dict):
            for key, val in phases.items():
                i, j = key
                self.phase[(int(i), int(j))] = float(val)
                if not self.config.directed and (int(j), int(i)) not in self.phase:
                    self.phase[(int(j), int(i))] = -float(val)
        else:
            ph_arr = np.asarray(phases, dtype=np.float64)
            if ph_arr.shape[0] != len(self.edges):
                raise ValueError("phases array must have length equal to number of edges")
            for (i, j), ph in zip(self.edges, ph_arr):
                self.phase[(int(i), int(j))] = float(ph)
                if not self.config.directed:
                    self.phase[(int(j), int(i))] = -float(ph)
        self.edge_phases = np.asarray([self.phase.get((int(i), int(j)), 0.0) for i, j in self.edges], dtype=np.float64)
        self._phase_dense = None
        if self.n_nodes <= self.config.dense_phase_max_nodes:
            dense = np.zeros((self.n_nodes, self.n_nodes), dtype=np.float64)
            for (i, j), ph in self.phase.items():
                if 0 <= i < self.n_nodes and 0 <= j < self.n_nodes:
                    dense[i, j] = ph
            self._phase_dense = dense

    def estimate_generators(self, **kwargs) -> Dict[str, Any]:
        if self.use_numba and len(self.edges):
            degree, imbalance = graph_degree_imbalance_numba(
                self.edges.astype(np.int64), self.weights.astype(np.float64),
                self.edge_phases.astype(np.float64), int(self.n_nodes), bool(self.config.directed)
            )
        else:
            degree = np.zeros(self.n_nodes, dtype=np.float64)
            imbalance = np.zeros(self.n_nodes, dtype=np.float64)
            for (i, j), w in zip(self.edges, self.weights):
                i, j = int(i), int(j)
                degree[i] += abs(w)
                imbalance[i] += self.phase.get((i, j), 0.0)
                if not self.config.directed:
                    degree[j] += abs(w)
                    imbalance[j] += self.phase.get((j, i), 0.0)
        gens = {"degree": degree, "phase_imbalance": imbalance}
        self.state.generators = gens
        self.state.observables.update({"degree": degree, "phase_imbalance": imbalance})
        return gens

    def build_connection(self, **kwargs) -> Dict[str, Any]:
        conn = {"edge_phase": dict(self.phase), "adjacency": self.adj}
        self.state.connection = conn
        return conn

    def detect_defects(self, **kwargs) -> List[Defect]:
        if self.state.connection is None:
            self.build_connection(**kwargs)
        cycles = self.cycle_basis(max_cycles=kwargs.get("max_cycles", self.config.cycle_basis_max))
        defects: List[Defect] = []
        if self.use_numba and self._phase_dense is not None and cycles:
            max_len = max(len(c) for c in cycles)
            cyc_arr = -np.ones((len(cycles), max_len), dtype=np.int64)
            lens = np.zeros(len(cycles), dtype=np.int64)
            for r, cyc in enumerate(cycles):
                lens[r] = len(cyc)
                cyc_arr[r, :len(cyc)] = np.asarray(cyc, dtype=np.int64)
            hol = cycle_holonomy_padded_numba(cyc_arr, lens, self._phase_dense)
            iterable = zip(cycles, hol)
        else:
            iterable = ((cyc, self.holonomy(cyc + [cyc[0]])) for cyc in cycles)
        for cyc, h in iterable:
            charge = float(h) / (2.0 * np.pi)
            if abs(charge) > self.config.defect_threshold:
                defects.append(
                    Defect(
                        kind="graph_holonomy",
                        charge=float(charge),
                        support=tuple(cyc),
                        strength=abs(float(h)),
                        metadata={"holonomy": float(h), "cycle_length": len(cyc), "numba_bulk": bool(self.use_numba and self._phase_dense is not None)},
                    )
                )
        self.state.defects = defects
        self.state.observables["n_cycles_checked"] = len(cycles)
        self.state.observables["n_defects"] = len(defects)
        return defects

    def holonomy(self, path: Sequence[int]) -> float:
        if len(path) < 2:
            return 0.0
        total = 0.0
        for a, b in zip(path[:-1], path[1:]):
            total += self.phase.get((int(a), int(b)), 0.0)
        return float(wrap_angle(total))

    def cycle_basis(self, max_cycles: int = 5000) -> List[List[int]]:
        if self._resolved_mode in {"efficient", "local"}:
            return self._cycle_basis_chords(max_cycles=max_cycles, dedupe=False)
        return self._cycle_basis_chords(max_cycles=max_cycles, dedupe=True)

    def _cycle_basis_chords(self, max_cycles: int, dedupe: bool) -> List[List[int]]:
        """Fundamental cycle basis from a spanning forest plus non-tree chords."""
        n = self.n_nodes
        parent = [-1] * n
        depth = [0] * n
        visited = [False] * n
        tree_edges = set()

        # Build spanning forest first. This avoids the v0.1 DFS timing bug
        # where a back-edge could be inspected before parent/depth stabilized.
        for root in range(n):
            if visited[root]:
                continue
            visited[root] = True
            stack = [root]
            while stack:
                u = stack.pop()
                for v in self.adj[u]:
                    if not visited[v]:
                        visited[v] = True
                        parent[v] = u
                        depth[v] = depth[u] + 1
                        ekey = tuple(sorted((u, v))) if not self.config.directed else (u, v)
                        tree_edges.add(ekey)
                        stack.append(v)

        cycles: List[List[int]] = []
        emitted = set()
        for u, v in self.edges:
            u, v = int(u), int(v)
            ekey = tuple(sorted((u, v))) if not self.config.directed else (u, v)
            if ekey in tree_edges or ekey in emitted:
                continue
            emitted.add(ekey)
            cyc = self._cycle_from_edge(u, v, parent, depth)
            if len(cyc) >= 3:
                cycles.append(cyc)
                if len(cycles) >= max_cycles:
                    return self._dedupe_cycles(cycles) if dedupe else cycles
        return self._dedupe_cycles(cycles)[:max_cycles] if dedupe else cycles[:max_cycles]

    def _dedupe_cycles(self, cycles: List[List[int]]) -> List[List[int]]:
        unique: Dict[Tuple[int, ...], List[int]] = {}
        for c in cycles:
            key = tuple(sorted(set(c)))
            unique.setdefault(key, c)
        return list(unique.values())

    def _cycle_from_edge(self, u: int, v: int, parent: List[int], depth: List[int]) -> List[int]:
        path_u = []
        path_v = []
        a, b = u, v
        seen_guard = 0
        while a != b and a != -1 and b != -1 and seen_guard < self.n_nodes + 5:
            seen_guard += 1
            if depth[a] >= depth[b]:
                path_u.append(a)
                a = parent[a]
            else:
                path_v.append(b)
                b = parent[b]
        if a == -1 or b == -1 or seen_guard >= self.n_nodes + 5:
            return []
        lca = a
        path_u.append(lca)
        return path_u + path_v[::-1]
