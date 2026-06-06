"""
optim.py
========

Optimization/search routines over defect states.

The module provides small, dependency-light tools that work across backends:
neutral-cycle search for parity matrices, low-action path search, holonomy path
ranking, and a simple defect annealer.  These are meant as baselines for the
DefectVM rather than as final state-of-the-art solvers.
"""

from __future__ import annotations

from dataclasses import dataclass
from heapq import heappop, heappush
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .core import Defect, DefectState, EPS

try:
    from numba import njit  # type: ignore
    NUMBA_AVAILABLE = True
except Exception:  # pragma: no cover
    NUMBA_AVAILABLE = False

    def njit(*args, **kwargs):  # type: ignore
        def deco(fn):
            return fn
        if args and callable(args[0]):
            return args[0]
        return deco


@dataclass
class OptimConfig:
    """Shared optimization configuration."""

    max_iter: int = 1000
    temperature: float = 1.0
    cooling: float = 0.995
    random_state: Optional[int] = None
    use_numba: bool = True


@dataclass
class SearchResult:
    """Generic search result."""

    success: bool
    value: Any = None
    score: float = np.inf
    path: Optional[List[int]] = None
    metadata: Dict[str, Any] = None


@njit(cache=True)
def _gf2_rref_numba(A: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    M = A.copy() & 1
    n_rows, n_cols = M.shape
    pivots = []
    r = 0
    for c in range(n_cols):
        piv = -1
        for i in range(r, n_rows):
            if M[i, c] == 1:
                piv = i
                break
        if piv < 0:
            continue
        if piv != r:
            tmp = M[r].copy()
            M[r] = M[piv]
            M[piv] = tmp
        for i in range(n_rows):
            if i != r and M[i, c] == 1:
                M[i] ^= M[r]
        pivots.append(c)
        r += 1
        if r == n_rows:
            break
    return M, np.asarray(pivots, dtype=np.int64)


def gf2_nullspace(A: np.ndarray, *, use_numba: bool = True) -> np.ndarray:
    """Return a basis for the right nullspace of binary matrix A over GF(2)."""
    M = (np.asarray(A, dtype=np.uint8) & 1)
    n_rows, n_cols = M.shape
    if use_numba and NUMBA_AVAILABLE:
        R, pivots = _gf2_rref_numba(M)
    else:
        R = M.copy()
        pivots_list: List[int] = []
        r = 0
        for c in range(n_cols):
            piv = None
            for i in range(r, n_rows):
                if R[i, c]:
                    piv = i
                    break
            if piv is None:
                continue
            if piv != r:
                R[[r, piv]] = R[[piv, r]]
            for i in range(n_rows):
                if i != r and R[i, c]:
                    R[i] ^= R[r]
            pivots_list.append(c)
            r += 1
            if r == n_rows:
                break
        pivots = np.asarray(pivots_list, dtype=int)
    pivot_set = set(int(p) for p in pivots)
    free_cols = [c for c in range(n_cols) if c not in pivot_set]
    if not free_cols:
        return np.zeros((0, n_cols), dtype=np.uint8)
    basis = []
    for f in free_cols:
        x = np.zeros(n_cols, dtype=np.uint8)
        x[f] = 1
        for row, p in enumerate(pivots):
            if R[row, f]:
                x[int(p)] = 1
        basis.append(x)
    return np.vstack(basis)


def neutral_cycle_search(parity_matrix: np.ndarray, *, left: bool = True, use_numba: bool = True) -> np.ndarray:
    """Find parity-neutral dependencies.

    If rows are relations and columns are defect channels, use ``left=True`` to
    return row-combination vectors ``c`` satisfying ``c^T M = 0``.
    """
    M = np.asarray(parity_matrix, dtype=np.uint8) & 1
    return gf2_nullspace(M.T if left else M, use_numba=use_numba)


def low_action_path(
    adjacency: Sequence[Sequence[Tuple[int, float]]],
    source: int,
    target: int,
    *,
    max_cost: float = np.inf,
) -> SearchResult:
    """Dijkstra search over a weighted adjacency list."""
    n = len(adjacency)
    dist = np.full(n, np.inf)
    prev = np.full(n, -1, dtype=int)
    dist[source] = 0.0
    heap = [(0.0, int(source))]
    while heap:
        d, u = heappop(heap)
        if d != dist[u]:
            continue
        if u == target:
            break
        if d > max_cost:
            break
        for v, w in adjacency[u]:
            nd = d + float(w)
            if nd < dist[v]:
                dist[v] = nd
                prev[v] = u
                heappush(heap, (nd, int(v)))
    if not np.isfinite(dist[target]):
        return SearchResult(False, score=np.inf, metadata={"reason": "unreachable"})
    path = []
    cur = int(target)
    while cur >= 0:
        path.append(cur)
        cur = int(prev[cur])
    path.reverse()
    return SearchResult(True, value=path, score=float(dist[target]), path=path, metadata={})


def rank_holonomy_paths(
    holonomy_fn: Callable[[Sequence[int]], Any],
    paths: Sequence[Sequence[int]],
    *,
    target: float = 0.0,
) -> List[SearchResult]:
    """Rank candidate paths by closeness of scalar holonomy to target."""
    results: List[SearchResult] = []
    for p in paths:
        h = holonomy_fn(p)
        hv = float(np.asarray(h).ravel()[0])
        score = abs(hv - target)
        results.append(SearchResult(True, value=h, score=score, path=list(p), metadata={"holonomy": hv}))
    results.sort(key=lambda r: r.score)
    return results


class DefectAnnealer:
    """Simple simulated annealer over user-defined defect states.

    The state can be any Python object. ``proposal_fn`` creates a nearby state;
    ``energy_fn`` evaluates it.  This is intentionally generic and can be used
    for relation subset search, graph path search, ordering search, etc.
    """

    def __init__(self, energy_fn: Callable[[Any], float], proposal_fn: Callable[[Any, np.random.Generator], Any], config: Optional[OptimConfig] = None):
        self.energy_fn = energy_fn
        self.proposal_fn = proposal_fn
        self.config = config or OptimConfig()
        self.rng = np.random.default_rng(self.config.random_state)

    def run(self, x0: Any) -> SearchResult:
        x = x0
        e = float(self.energy_fn(x))
        best_x, best_e = x, e
        T = float(self.config.temperature)
        history = []
        for it in range(int(self.config.max_iter)):
            xp = self.proposal_fn(x, self.rng)
            ep = float(self.energy_fn(xp))
            accept = ep <= e or self.rng.random() < np.exp(-(ep - e) / (T + EPS))
            if accept:
                x, e = xp, ep
                if e < best_e:
                    best_x, best_e = x, e
            history.append(best_e)
            T *= self.config.cooling
        return SearchResult(True, value=best_x, score=best_e, metadata={"history": np.asarray(history)})


def optimize_mps_ordering(omega: np.ndarray, *, random_state: Optional[int] = None, max_iter: int = 500) -> SearchResult:
    """Heuristic ordering optimizer: minimizes weighted edge span of Omega."""
    O = np.abs(np.asarray(omega, dtype=float))
    n = O.shape[0]
    rng = np.random.default_rng(random_state)

    def energy(order: np.ndarray) -> float:
        pos = np.empty(n, dtype=int)
        pos[order] = np.arange(n)
        # weighted linear arrangement
        ii, jj = np.triu_indices(n, 1)
        return float(np.sum(O[ii, jj] * np.abs(pos[ii] - pos[jj])))

    def proposal(order: np.ndarray, r: np.random.Generator) -> np.ndarray:
        out = order.copy()
        i, j = r.choice(n, size=2, replace=False)
        out[i], out[j] = out[j], out[i]
        return out

    annealer = DefectAnnealer(energy, proposal, OptimConfig(max_iter=max_iter, random_state=random_state))
    return annealer.run(np.arange(n, dtype=int))


def register_vm_ops(vm: Any) -> Any:
    """Register optimization helpers on a DefectVM instance."""
    vm.register("neutral_cycles", lambda parity_matrix=None, **kw: neutral_cycle_search(parity_matrix if parity_matrix is not None else vm.state.observables.get("parity_matrix"), **kw))
    vm.register("low_action_path", low_action_path)
    vm.register("rank_holonomy_paths", lambda paths, target=0.0: rank_holonomy_paths(vm.op_holonomy, paths, target=target))
    vm.register("optimize_mps_ordering", lambda omega=None, **kw: optimize_mps_ordering(vm.state.omega if omega is None else omega, **kw))
    return vm


__all__ = [
    "OptimConfig",
    "SearchResult",
    "gf2_nullspace",
    "neutral_cycle_search",
    "low_action_path",
    "rank_holonomy_paths",
    "DefectAnnealer",
    "optimize_mps_ordering",
    "register_vm_ops",
]
