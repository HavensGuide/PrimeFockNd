
"""
numba_kernels.py
================

Central optional Numba kernels for Causal Defect Engine v0.4.

Every public wrapper has a pure NumPy/Python fallback through the backend code;
these kernels are imported opportunistically by modules that can benefit from
compiled inner loops.  The goal is to accelerate hot loops without making Numba
semantically required.
"""

from __future__ import annotations

import math
from typing import Tuple

import numpy as np

from .core import njit, prange, NUMBA_AVAILABLE, EPS


@njit(cache=True, parallel=True)
def local_gaussian_score_knn_numba(X: np.ndarray, nn: np.ndarray, bw: float) -> np.ndarray:
    """Neighbor-local Gaussian KDE score: grad log rho at fitted points.

    score_i = - sum_j softmax(-||x_i-x_j||^2/2bw^2) (x_i-x_j) / bw^2
    """
    n, d = X.shape
    k = nn.shape[1]
    out = np.zeros((n, d), dtype=np.float64)
    inv_bw2 = 1.0 / (bw * bw + EPS)
    for i in prange(n):
        # stable softmax over neighbors
        maxlog = -1.0e300
        for q in range(k):
            j = int(nn[i, q])
            s = 0.0
            for c in range(d):
                z = X[i, c] - X[j, c]
                s += z * z
            lk = -0.5 * s * inv_bw2
            if lk > maxlog:
                maxlog = lk
        denom = 0.0
        # first accumulate denominator
        for q in range(k):
            j = int(nn[i, q])
            s = 0.0
            for c in range(d):
                z = X[i, c] - X[j, c]
                s += z * z
            w = math.exp(-0.5 * s * inv_bw2 - maxlog)
            denom += w
        if denom <= EPS:
            continue
        for q in range(k):
            j = int(nn[i, q])
            s = 0.0
            for c in range(d):
                z = X[i, c] - X[j, c]
                s += z * z
            w = math.exp(-0.5 * s * inv_bw2 - maxlog) / denom
            for c in range(d):
                out[i, c] -= w * (X[i, c] - X[j, c]) * inv_bw2
    return out


@njit(cache=True, parallel=True)
def local_time_drift_knn_numba(X: np.ndarray, t: np.ndarray, nn: np.ndarray) -> np.ndarray:
    """Local regression drift dX/dt over each point's neighbors."""
    n, d = X.shape
    k = nn.shape[1]
    out = np.zeros((n, d), dtype=np.float64)
    for i in prange(n):
        tm = 0.0
        for q in range(k):
            tm += t[int(nn[i, q])]
        tm /= max(k, 1)
        var = 0.0
        xm = np.zeros(d, dtype=np.float64)
        for q in range(k):
            j = int(nn[i, q])
            tc = t[j] - tm
            var += tc * tc
            for c in range(d):
                xm[c] += X[j, c]
        for c in range(d):
            xm[c] /= max(k, 1)
        if var <= EPS:
            continue
        for q in range(k):
            j = int(nn[i, q])
            tc = t[j] - tm
            for c in range(d):
                out[i, c] += tc * (X[j, c] - xm[c])
        for c in range(d):
            out[i, c] /= var
    return out


@njit(cache=True, parallel=True)
def rank2_frame_numba(u: np.ndarray, v: np.ndarray, alpha: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Leading frame of G = u u^T + alpha v v^T using its 2D span.

    Avoids per-point d x d eigendecomposition. Returns theta, strength,
    frame.  The frame sign is gauge-fixed by making its largest component
    positive.
    """
    n, d = u.shape
    frame = np.zeros((n, d), dtype=np.float64)
    strength = np.zeros(n, dtype=np.float64)
    theta = np.zeros(n, dtype=np.float64)
    sa = math.sqrt(max(alpha, 0.0))
    for i in prange(n):
        a = 0.0
        b = 0.0
        c = 0.0
        for q in range(d):
            a += u[i, q] * u[i, q]
            b += u[i, q] * (sa * v[i, q])
            c += (sa * v[i, q]) * (sa * v[i, q])
        tr = a + c
        disc = (a - c) * (a - c) + 4.0 * b * b
        lam = 0.5 * (tr + math.sqrt(max(disc, 0.0)))
        strength[i] = max(lam, 0.0)
        # eigenvector y of [[a,b],[b,c]] for largest eigenvalue
        y0 = b
        y1 = lam - a
        yn = math.sqrt(y0 * y0 + y1 * y1)
        if yn <= EPS:
            if a >= c:
                y0, y1 = 1.0, 0.0
            else:
                y0, y1 = 0.0, 1.0
        else:
            y0 /= yn
            y1 /= yn
        norm = 0.0
        for q in range(d):
            val = y0 * u[i, q] + y1 * sa * v[i, q]
            frame[i, q] = val
            norm += val * val
        norm = math.sqrt(norm)
        if norm > EPS:
            for q in range(d):
                frame[i, q] /= norm
        # sign gauge fix
        pivot = 0
        mx = 0.0
        for q in range(d):
            aa = abs(frame[i, q])
            if aa > mx:
                mx = aa
                pivot = q
        if frame[i, pivot] < 0.0:
            for q in range(d):
                frame[i, q] = -frame[i, q]
        if d >= 2:
            theta[i] = math.atan2(frame[i, 1], frame[i, 0])
    return theta, strength, frame


@njit(cache=True, parallel=True)
def graph_degree_imbalance_numba(edges: np.ndarray, weights: np.ndarray, phases: np.ndarray, n_nodes: int, directed: bool) -> Tuple[np.ndarray, np.ndarray]:
    """Degree and phase imbalance for edge arrays.

    Note: parallel atomic adds are not universally available in Numba CPU, so
    this uses a serial-safe implementation despite the parallel signature.
    For very large graphs, backends can chunk by edge externally.
    """
    degree = np.zeros(n_nodes, dtype=np.float64)
    imbalance = np.zeros(n_nodes, dtype=np.float64)
    for e in range(edges.shape[0]):
        i = int(edges[e, 0])
        j = int(edges[e, 1])
        if i < 0 or j < 0 or i >= n_nodes or j >= n_nodes:
            continue
        w = abs(weights[e])
        ph = phases[e]
        degree[i] += w
        imbalance[i] += ph
        if not directed:
            degree[j] += w
            imbalance[j] -= ph
    return degree, imbalance


@njit(cache=True)
def wrap_angle_numba(x: float) -> float:
    return (x + math.pi) % (2.0 * math.pi) - math.pi


@njit(cache=True, parallel=True)
def cycle_holonomy_padded_numba(cycles: np.ndarray, lengths: np.ndarray, phase_dense: np.ndarray) -> np.ndarray:
    """Bulk holonomy for padded closed/open cycles using dense phase lookup."""
    ncyc = cycles.shape[0]
    out = np.zeros(ncyc, dtype=np.float64)
    for r in prange(ncyc):
        L = int(lengths[r])
        total = 0.0
        if L >= 2:
            for k in range(L - 1):
                a = int(cycles[r, k])
                b = int(cycles[r, k + 1])
                total += phase_dense[a, b]
            # close cycle if not explicitly closed
            if int(cycles[r, 0]) != int(cycles[r, L - 1]):
                total += phase_dense[int(cycles[r, L - 1]), int(cycles[r, 0])]
        out[r] = wrap_angle_numba(total)
    return out


@njit(cache=True, parallel=True)
def trial_factor_candidates_numba(N: int, center: int, radius: int, base: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Trial factor Q(a)=a^2-N over factor base for a in [center-radius, center+radius].

    Returns a, Q, exponent matrix, remainder, sign column.
    """
    n_cand = 2 * radius + 1
    m = len(base)
    a_arr = np.empty(n_cand, dtype=np.int64)
    q_arr = np.empty(n_cand, dtype=np.int64)
    exps = np.zeros((n_cand, m), dtype=np.int64)
    rems = np.empty(n_cand, dtype=np.int64)
    signs = np.zeros(n_cand, dtype=np.int64)
    for idx in prange(n_cand):
        a = center - radius + idx
        Q = a * a - N
        a_arr[idx] = a
        q_arr[idx] = Q
        rem = Q
        if rem < 0:
            signs[idx] = 1
            rem = -rem
        for j in range(m):
            p = int(base[j])
            while rem % p == 0:
                exps[idx, j] += 1
                rem //= p
        rems[idx] = rem
    return a_arr, q_arr, exps, rems, signs


@njit(cache=True, parallel=True)
def pack_f2_rows_numba(M: np.ndarray) -> np.ndarray:
    n, m = M.shape
    words = (m + 63) // 64
    out = np.zeros((n, words), dtype=np.uint64)
    for i in prange(n):
        for j in range(m):
            if M[i, j] & 1:
                out[i, j // 64] |= np.uint64(1) << np.uint64(j % 64)
    return out


@njit(cache=True)
def f2_rref_numba(A: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int]:
    """Reduced row echelon form of uint8 matrix over F2.

    Returns transformed matrix, pivot columns array, rank.
    """
    M = (A.copy() & np.uint8(1))
    n_rows, n_cols = M.shape
    max_rank = min(n_rows, n_cols)
    pivots = -np.ones(max_rank, dtype=np.int64)
    r = 0
    for c in range(n_cols):
        pivot = -1
        for i in range(r, n_rows):
            if M[i, c] & 1:
                pivot = i
                break
        if pivot == -1:
            continue
        if pivot != r:
            tmp = M[r].copy()
            M[r] = M[pivot]
            M[pivot] = tmp
        for i in range(n_rows):
            if i != r and (M[i, c] & 1):
                for j in range(c, n_cols):
                    M[i, j] ^= M[r, j]
        pivots[r] = c
        r += 1
        if r == n_rows:
            break
    return M, pivots, r


@njit(cache=True, parallel=True)
def corr_from_features_numba(X: np.ndarray) -> np.ndarray:
    """Correlation matrix for samples x features."""
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
        stds[j] = math.sqrt(s / max(n - 1, 1)) + EPS
    C = np.empty((d, d), dtype=np.float64)
    for a in prange(d):
        for b in range(d):
            s = 0.0
            for i in range(n):
                s += (X[i, a] - means[a]) * (X[i, b] - means[b])
            C[a, b] = s / max(n - 1, 1) / (stds[a] * stds[b])
    return C

# ---------------------------------------------------------------------------
# v0.7 state-geometry compiler kernels
# ---------------------------------------------------------------------------

@njit(cache=True)
def threshold_edges_numba(W: np.ndarray, threshold: float, max_edges: int) -> np.ndarray:
    """Return undirected edges (i,j) with abs(W[i,j]) >= threshold, capped.

    The output is padded/truncated to max_edges rows; invalid rows are (-1,-1).
    This avoids Python-loop overhead for dense Omega threshold scans.
    """
    n = W.shape[0]
    out = -np.ones((max_edges, 2), dtype=np.int64)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if abs(W[i, j]) >= threshold:
                if c < max_edges:
                    out[c, 0] = i
                    out[c, 1] = j
                c += 1
    return out


@njit(cache=True, parallel=True)
def edge_connection_binary_numba(W: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Binary gluing matrices [[1+w,1-w],[1-w,1+w]] for edge weights."""
    m = edges.shape[0]
    conns = np.zeros((m, 2, 2), dtype=np.float64)
    for e in prange(m):
        i = int(edges[e, 0])
        j = int(edges[e, 1])
        if i < 0 or j < 0:
            continue
        w = abs(W[i, j])
        conns[e, 0, 0] = 1.0 + w
        conns[e, 1, 1] = 1.0 + w
        conns[e, 0, 1] = 1.0 - w
        conns[e, 1, 0] = 1.0 - w
    return conns


@njit(cache=True)
def triangle_curvature_candidates_numba(W: np.ndarray, edges: np.ndarray, tol: float, max_defects: int) -> tuple:
    """Triangle curvature proxy over an undirected edge list.

    Returns padded loops (max_defects,3) and charges.  This is intended for
    moderate sparse graphs; Python code sorts/trims the returned nonzero rows.
    """
    n = W.shape[0]
    adj = np.zeros((n, n), dtype=np.uint8)
    for e in range(edges.shape[0]):
        i = int(edges[e, 0])
        j = int(edges[e, 1])
        if i >= 0 and j >= 0:
            adj[i, j] = 1
            adj[j, i] = 1
    loops = -np.ones((max_defects, 3), dtype=np.int64)
    charges = np.zeros(max_defects, dtype=np.float64)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if adj[i, j] == 0:
                continue
            for k in range(j + 1, n):
                if adj[i, k] and adj[j, k]:
                    curv = W[i, j] * W[j, k] * W[k, i]
                    if abs(curv) > tol:
                        if c < max_defects:
                            loops[c, 0] = i
                            loops[c, 1] = j
                            loops[c, 2] = k
                            charges[c] = curv
                        c += 1
    return loops, charges


@njit(cache=True)
def tree_gauge_residual_numba(edges: np.ndarray, phases: np.ndarray, n_nodes: int) -> tuple:
    """Gauge-fix edge U(1) phases along a spanning forest.

    A BFS-like pass sets node gauges so tree-edge residual phases are zero.
    Returns node_gauge, residual_phases, tree_mask.
    """
    m = edges.shape[0]
    adj_count = np.zeros(n_nodes, dtype=np.int64)
    for e in range(m):
        i = int(edges[e, 0]); j = int(edges[e, 1])
        if i >= 0 and j >= 0:
            adj_count[i] += 1; adj_count[j] += 1
    offs = np.zeros(n_nodes + 1, dtype=np.int64)
    for i in range(n_nodes):
        offs[i + 1] = offs[i] + adj_count[i]
    ptr = offs.copy()
    nbr = np.empty(offs[-1], dtype=np.int64)
    eid = np.empty(offs[-1], dtype=np.int64)
    sgn = np.empty(offs[-1], dtype=np.float64)
    for e in range(m):
        i = int(edges[e, 0]); j = int(edges[e, 1])
        if i < 0 or j < 0:
            continue
        p = ptr[i]; nbr[p] = j; eid[p] = e; sgn[p] = 1.0; ptr[i] += 1
        p = ptr[j]; nbr[p] = i; eid[p] = e; sgn[p] = -1.0; ptr[j] += 1
    gauge = np.zeros(n_nodes, dtype=np.float64)
    seen = np.zeros(n_nodes, dtype=np.uint8)
    tree_mask = np.zeros(m, dtype=np.uint8)
    queue = np.empty(n_nodes, dtype=np.int64)
    for root in range(n_nodes):
        if seen[root]:
            continue
        head = 0; tail = 0
        queue[tail] = root; tail += 1; seen[root] = 1; gauge[root] = 0.0
        while head < tail:
            u = queue[head]; head += 1
            for p in range(offs[u], offs[u + 1]):
                v = nbr[p]; e = eid[p]
                if seen[v] == 0:
                    # Residual edge phase phi + g_u - g_v = 0 for orientation i->j.
                    phi = phases[e] * sgn[p]
                    gauge[v] = gauge[u] + phi
                    seen[v] = 1
                    tree_mask[e] = 1
                    queue[tail] = v; tail += 1
    residual = np.zeros(m, dtype=np.float64)
    for e in range(m):
        i = int(edges[e, 0]); j = int(edges[e, 1])
        val = phases[e] + gauge[i] - gauge[j]
        residual[e] = (val + math.pi) % (2.0 * math.pi) - math.pi
    return gauge, residual, tree_mask


@njit(cache=True, parallel=True)
def local_parity_energy_numba(configs: np.ndarray, edges: np.ndarray, weights: np.ndarray, fields: np.ndarray) -> np.ndarray:
    """Batch Ising/log-linear energy helper for boundary configurations."""
    ncfg, n = configs.shape
    m = edges.shape[0]
    out = np.zeros(ncfg, dtype=np.float64)
    for r in prange(ncfg):
        e = 0.0
        for i in range(n):
            s = 1.0 if configs[r, i] > 0 else -1.0
            e -= fields[i] * s
        for q in range(m):
            i = int(edges[q, 0]); j = int(edges[q, 1])
            si = 1.0 if configs[r, i] > 0 else -1.0
            sj = 1.0 if configs[r, j] > 0 else -1.0
            e -= weights[q] * si * sj
        out[r] = e
    return out
