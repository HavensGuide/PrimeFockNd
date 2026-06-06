#!/usr/bin/env python3
"""
Prime-Fock / Shor simulation benchmark.

What this script does
---------------------
1. Builds the exact modular-exponentiation Shor state for small semiprimes:

       |psi> = 1/sqrt(Q) sum_x |x>|a^x mod N>

   without building a gate-level circuit. This is the post-modexp state that
   matters for period structure.

2. Computes quantum-simulation diagnostics:
   - exact period r
   - success/factorization check
   - Schmidt ranks and entanglement entropy across the exponent|residue cut
   - optional binary-qubit cut entropies for tiny systems

3. Computes arithmetic / Prime-Fock-ish orbit compression diagnostics:
   - residue orbit support
   - valuation features of residues and residue differences
   - gcd-defect hits gcd(a^d - 1, N)
   - random sparse period recovery sample complexity
   - compressed orbit graph spectral diagnostics
   - value-vs-order compression diagnostics: valuation rank, cyclic successor rank,
     and block-Hankel rank of the Prime-Fock valuation sequence

4. Compares against simple baselines:
   - classical sequential orbit enumeration
   - random sample gcd-defect recovery

This is deliberately honest: it is a small-N benchmark to see whether the
arithmetic/orbit compression story has teeth. It is not a scalable Shor claim.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import random
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, Optional

import numpy as np


# ------------------------- number theory helpers -------------------------


def gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return abs(a)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def factor_trial(n: int) -> tuple[int, int] | None:
    for d in range(2, int(math.isqrt(n)) + 1):
        if n % d == 0:
            return d, n // d
    return None


def semiprimes(max_n: int, min_n: int = 15) -> list[int]:
    out = []
    for n in range(min_n, max_n + 1):
        fac = factor_trial(n)
        if fac is None:
            continue
        p, q = fac
        if p * q == n and is_prime(p) and is_prime(q):
            out.append(n)
    return out


def multiplicative_order(a: int, n: int, max_steps: Optional[int] = None) -> Optional[int]:
    if gcd(a, n) != 1:
        return None
    if max_steps is None:
        max_steps = n + 1
    x = 1
    for r in range(1, max_steps + 1):
        x = (x * a) % n
        if x == 1:
            return r
    return None


def valid_bases(n: int, max_bases: int = 8) -> list[int]:
    bases = []
    for a in range(2, n):
        if gcd(a, n) != 1:
            continue
        r = multiplicative_order(a, n)
        if r is None or r <= 1:
            continue
        bases.append(a)
        if len(bases) >= max_bases:
            break
    return bases


def shor_factor_from_period(a: int, n: int, r: Optional[int]) -> tuple[bool, int | None, int | None]:
    if r is None or r % 2 == 1:
        return False, None, None
    y = pow(a, r // 2, n)
    if y == n - 1 or y == 1:
        return False, None, None
    p = gcd(y - 1, n)
    q = gcd(y + 1, n)
    ok = 1 < p < n and 1 < q < n and p * q == n
    return ok, p if 1 < p < n else None, q if 1 < q < n else None


def next_power_of_two_at_least(n: int) -> int:
    return 1 << ((n - 1).bit_length())


# ------------------------- exact Shor state diagnostics -------------------------


def shor_state_matrix(n: int, a: int, q: Optional[int] = None) -> tuple[np.ndarray, np.ndarray]:
    """Return exponent|residue amplitude matrix M[x, y]."""
    if q is None:
        # Standard phase estimation often uses Q ~ power of two with N^2 <= Q < 2N^2.
        q = next_power_of_two_at_least(n * n)
    residues = np.fromiter((pow(a, x, n) for x in range(q)), dtype=np.int64, count=q)
    m = np.zeros((q, n), dtype=np.complex128)
    amp = 1.0 / math.sqrt(q)
    m[np.arange(q), residues] = amp
    return m, residues


def schmidt_entropy_from_matrix(m: np.ndarray, tol: float = 1e-12) -> tuple[int, float, float, np.ndarray]:
    """Schmidt rank, von Neumann entropy base2, participation rank, singular values."""
    s = np.linalg.svd(m, compute_uv=False)
    p = s * s
    p = p[p > tol]
    rank = int(p.size)
    entropy = float(-np.sum(p * np.log2(p))) if p.size else 0.0
    participation = float(1.0 / np.sum(p * p)) if p.size else 0.0
    return rank, entropy, participation, s


def binary_cut_entropies(n: int, a: int, q: Optional[int] = None, max_total_qubits: int = 18) -> list[float] | None:
    """Exact qubit cut entropies for the full padded binary state. Tiny systems only."""
    if q is None:
        q = next_power_of_two_at_least(n * n)
    qbits_x = int(math.log2(q))
    y_dim = next_power_of_two_at_least(n)
    qbits_y = int(math.log2(y_dim))
    total = qbits_x + qbits_y
    if total > max_total_qubits:
        return None
    state = np.zeros(q * y_dim, dtype=np.complex128)
    amp = 1.0 / math.sqrt(q)
    for x in range(q):
        y = pow(a, x, n)
        idx = (x << qbits_y) | y
        state[idx] = amp
    ent = []
    for cut in range(1, total):
        left = 1 << cut
        right = 1 << (total - cut)
        mat = state.reshape(left, right)
        _, entropy, _, _ = schmidt_entropy_from_matrix(mat)
        ent.append(entropy)
    return ent


# ------------------------- Prime-Fock-ish orbit diagnostics -------------------------


def first_primes(k: int) -> list[int]:
    primes = []
    x = 2
    while len(primes) < k:
        if is_prime(x):
            primes.append(x)
        x += 1
    return primes


def valuation_vector(x: int, primes: list[int]) -> np.ndarray:
    """Prime valuation vector for |x|, with zero vector for 0/1."""
    x = abs(int(x))
    v = np.zeros(len(primes), dtype=np.float64)
    if x <= 1:
        return v
    for i, p in enumerate(primes):
        while x % p == 0:
            v[i] += 1.0
            x //= p
        if x == 1:
            break
    return v


def orbit_residues(n: int, a: int, length: int) -> np.ndarray:
    return np.fromiter((pow(a, x, n) for x in range(length)), dtype=np.int64, count=length)


def valuation_orbit_features(n: int, a: int, r: int, n_primes: int = 16) -> dict[str, float]:
    primes = first_primes(n_primes)
    residues = orbit_residues(n, a, r)
    vals = np.stack([valuation_vector(int(y), primes) for y in residues], axis=0)
    diffs = []
    for i in range(r):
        d = int((residues[(i + 1) % r] - residues[i]) % n)
        diffs.append(valuation_vector(d, primes))
    dvals = np.stack(diffs, axis=0)
    active_modes = float(np.sum(vals.sum(axis=0) > 0))
    active_diff_modes = float(np.sum(dvals.sum(axis=0) > 0))
    support = int(len(set(int(x) for x in residues)))
    mode_mass = vals.sum(axis=0)
    total = float(mode_mass.sum())
    if total > 0:
        p = mode_mass / total
        entropy = float(-np.sum(p[p > 0] * np.log2(p[p > 0])))
    else:
        entropy = 0.0
    return {
        "orbit_support": float(support),
        "valuation_active_modes": active_modes,
        "valuation_active_diff_modes": active_diff_modes,
        "valuation_mass": float(vals.sum()),
        "valuation_diff_mass": float(dvals.sum()),
        "valuation_mode_entropy": entropy,
        "valuation_support_over_period": float(support / max(r, 1)),
        "valuation_active_modes_over_logN": float(active_modes / max(math.log2(n), 1.0)),
    }


def gcd_defect_period_candidates(n: int, a: int, max_d: int) -> list[tuple[int, int]]:
    """Return d,g where gcd(a^d-1,N) is a nontrivial factor."""
    out = []
    for d in range(1, max_d + 1):
        g = gcd(pow(a, d, n) - 1, n)
        if 1 < g < n:
            out.append((d, g))
    return out


def sequential_period_search(n: int, a: int, max_steps: Optional[int] = None) -> tuple[Optional[int], int]:
    if max_steps is None:
        max_steps = n + 1
    x = 1
    for steps in range(1, max_steps + 1):
        x = (x * a) % n
        if x == 1:
            return steps, steps
    return None, max_steps


def random_sparse_gcd_recovery(
    n: int,
    a: int,
    r: int,
    *,
    trials: int = 200,
    max_samples: int = 64,
    seed: int = 0,
) -> dict[str, float]:
    """Try random d samples; success if gcd(a^d - 1,N) gives nontrivial factor or d hits multiple of r."""
    rng = random.Random(seed)
    successes = []
    samples_needed = []
    for _ in range(trials):
        seen = set()
        ok = False
        for k in range(1, max_samples + 1):
            # Bias d logarithmically and uniformly mixed, otherwise hitting r exactly is too rare.
            if rng.random() < 0.5:
                d = rng.randint(1, max(2 * r, 2))
            else:
                d = rng.randint(1, max(n, 2 * r, 2))
            if d in seen:
                continue
            seen.add(d)
            g = gcd(pow(a, d, n) - 1, n)
            if (1 < g < n) or (d % r == 0):
                ok = True
                samples_needed.append(k)
                break
        successes.append(ok)
    if samples_needed:
        mean_samples = float(np.mean(samples_needed))
        med_samples = float(np.median(samples_needed))
    else:
        mean_samples = float("nan")
        med_samples = float("nan")
    return {
        "sparse_gcd_success_rate": float(np.mean(successes)),
        "sparse_gcd_mean_samples_if_success": mean_samples,
        "sparse_gcd_median_samples_if_success": med_samples,
    }


def orbit_graph_spectrum(n: int, a: int, r: int) -> dict[str, float]:
    """Simple residue orbit graph with edges by modular successor and gcd-defect weights."""
    residues = orbit_residues(n, a, r)
    idx = {int(y): i for i, y in enumerate(sorted(set(int(z) for z in residues)))}
    m = len(idx)
    if m <= 1:
        return {"orbit_graph_nodes": float(m), "orbit_graph_rank": 0.0, "orbit_graph_spectral_entropy": 0.0}
    W = np.zeros((m, m), dtype=np.float64)
    for i in range(r):
        u = idx[int(residues[i])]
        v = idx[int(residues[(i + 1) % r])]
        diff = abs(int(residues[(i + 1) % r]) - int(residues[i]))
        weight = math.log1p(gcd(diff, n)) + 1.0
        W[u, v] += weight
        W[v, u] += weight
    vals = np.linalg.eigvalsh(W)
    absvals = np.abs(vals)
    rank = float(np.sum(absvals > 1e-10))
    total = float(absvals.sum())
    if total > 0:
        p = absvals / total
        entropy = float(-np.sum(p[p > 0] * np.log2(p[p > 0])))
    else:
        entropy = 0.0
    return {"orbit_graph_nodes": float(m), "orbit_graph_rank": rank, "orbit_graph_spectral_entropy": entropy}





def normalized_spectral_entropy_from_values(values: np.ndarray, tol: float = 1e-12) -> tuple[float, float]:
    """Return raw spectral entropy and entropy/log2(rank) normalization."""
    vals = np.asarray(values, dtype=np.float64)
    vals = np.abs(vals[np.abs(vals) > tol])
    if vals.size == 0:
        return 0.0, 0.0
    p = vals / float(vals.sum())
    ent = float(-np.sum(p * np.log2(p)))
    norm = float(ent / math.log2(vals.size)) if vals.size > 1 else 0.0
    return ent, norm


def valuation_sequence_matrix(n: int, a: int, r: int, n_primes: int = 16) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return residues, residue valuation matrix V[t,mode], and successor-diff valuation matrix DV[t,mode]."""
    primes = first_primes(n_primes)
    residues = orbit_residues(n, a, r)
    V = np.stack([valuation_vector(int(y), primes) for y in residues], axis=0)
    diffs = []
    for i in range(r):
        d = int((residues[(i + 1) % r] - residues[i]) % n)
        diffs.append(valuation_vector(d, primes))
    DV = np.stack(diffs, axis=0)
    return residues, V, DV


def matrix_rank_svd(X: np.ndarray, tol: float = 1e-10) -> tuple[float, float, float, float]:
    """Rank, stable rank, spectral entropy, normalized spectral entropy of a matrix."""
    X = np.asarray(X, dtype=np.float64)
    if X.size == 0 or min(X.shape) == 0:
        return 0.0, 0.0, 0.0, 0.0
    s = np.linalg.svd(X, compute_uv=False)
    rank = float(np.sum(s > tol))
    if np.sum(s * s) > 0 and np.max(s) > 0:
        stable_rank = float(np.sum(s * s) / (np.max(s) ** 2))
    else:
        stable_rank = 0.0
    ent, ent_norm = normalized_spectral_entropy_from_values(s, tol=tol)
    return rank, stable_rank, ent, ent_norm


def block_hankel_rank(X: np.ndarray, window: Optional[int] = None, tol: float = 1e-10) -> tuple[float, float, float, float, float]:
    """Block-Hankel rank for a cyclic vector sequence X[t,d]. Tracks order/sequence complexity in valuation space."""
    X = np.asarray(X, dtype=np.float64)
    T, D = X.shape
    if T <= 2 or D == 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0
    if window is None:
        window = max(2, min(8, T // 2))
    window = int(max(2, min(window, T)))
    # Cyclic block-Hankel: row t = [X_t, X_{t+1}, ..., X_{t+w-1}]
    H = np.zeros((T, window * D), dtype=np.float64)
    for t in range(T):
        blocks = [X[(t + j) % T] for j in range(window)]
        H[t] = np.concatenate(blocks)
    rank, stable_rank, ent, ent_norm = matrix_rank_svd(H, tol=tol)
    return rank, stable_rank, ent, ent_norm, float(window)


def cycle_order_complexity(n: int, a: int, r: int, n_primes: int = 16, hankel_window: Optional[int] = None) -> dict[str, float]:
    """
    Separate value compression from order/cycle compression.

    Value compression:
      - rank/entropy of residue valuation matrix V[t,mode]
      - active modes and valuation rank ratios vs period

    Order compression:
      - directed successor adjacency/permutation rank for the orbit cycle
      - block-Hankel rank of the valuation sequence, which tests whether the
        time ordering is low-rank in Prime-Fock coordinates

    If valuation_rank << r but cycle/hankel ranks ~ r, values compress but the
    orbit order does not. If both are << r, the cycle itself is compressing.
    """
    residues, V, DV = valuation_sequence_matrix(n, a, r, n_primes=n_primes)
    support = len(set(int(z) for z in residues))

    v_rank, v_stable, v_ent, v_ent_norm = matrix_rank_svd(V)
    dv_rank, dv_stable, dv_ent, dv_ent_norm = matrix_rank_svd(DV)
    h_rank, h_stable, h_ent, h_ent_norm, h_window = block_hankel_rank(V, window=hankel_window)
    dh_rank, dh_stable, dh_ent, dh_ent_norm, dh_window = block_hankel_rank(DV, window=hankel_window)

    # Directed successor matrix in residue-state basis. For a true period-r orbit, this is a permutation cycle.
    idx = {int(y): i for i, y in enumerate(residues)}
    m = len(idx)
    P = np.zeros((m, m), dtype=np.float64)
    for i in range(r):
        u = idx[int(residues[i])]
        v = idx[int(residues[(i + 1) % r])]
        P[u, v] = 1.0
    p_rank, p_stable, p_ent, p_ent_norm = matrix_rank_svd(P)

    # Undirected unweighted cycle adjacency, separate from gcd-weighted orbit_graph_spectrum.
    A = np.maximum(P, P.T)
    eig = np.linalg.eigvalsh(A) if m else np.array([])
    a_rank = float(np.sum(np.abs(eig) > 1e-10))
    a_ent, a_ent_norm = normalized_spectral_entropy_from_values(eig)

    denom_r = float(max(r, 1))
    denom_support = float(max(support, 1))
    return {
        "valuation_matrix_rank": v_rank,
        "valuation_matrix_stable_rank": v_stable,
        "valuation_matrix_spectral_entropy": v_ent,
        "valuation_matrix_spectral_entropy_norm": v_ent_norm,
        "valuation_rank_over_period": float(v_rank / denom_r),
        "valuation_rank_over_support": float(v_rank / denom_support),
        "diff_valuation_matrix_rank": dv_rank,
        "diff_valuation_matrix_stable_rank": dv_stable,
        "diff_valuation_rank_over_period": float(dv_rank / denom_r),
        "valuation_hankel_window": h_window,
        "valuation_hankel_rank": h_rank,
        "valuation_hankel_stable_rank": h_stable,
        "valuation_hankel_spectral_entropy": h_ent,
        "valuation_hankel_spectral_entropy_norm": h_ent_norm,
        "valuation_hankel_rank_over_period": float(h_rank / denom_r),
        "diff_valuation_hankel_rank": dh_rank,
        "diff_valuation_hankel_stable_rank": dh_stable,
        "diff_valuation_hankel_rank_over_period": float(dh_rank / denom_r),
        "successor_permutation_rank": p_rank,
        "successor_permutation_stable_rank": p_stable,
        "successor_permutation_spectral_entropy": p_ent,
        "successor_permutation_spectral_entropy_norm": p_ent_norm,
        "successor_rank_over_period": float(p_rank / denom_r),
        "cycle_adjacency_rank": a_rank,
        "cycle_adjacency_spectral_entropy": a_ent,
        "cycle_adjacency_spectral_entropy_norm": a_ent_norm,
        "cycle_adjacency_rank_over_period": float(a_rank / denom_r),
        "value_compresses_but_order_not": float((v_rank / denom_r < 0.5) and (h_rank / denom_r > 0.75)),
        "order_compresses_in_valuation_basis": float(h_rank / denom_r < 0.5),
    }


@dataclass
class ShorBenchRow:
    N: int
    p: int | None
    q_factor: int | None
    a: int
    Q: int
    n_bits: int
    qbits_x: int
    qbits_y: int
    total_qubits_padded: int
    period_r: int | None
    shor_factor_success: bool
    factor_p_from_period: int | None
    factor_q_from_period: int | None
    exact_state_seconds: float
    exp_residue_schmidt_rank: int
    exp_residue_entropy_bits: float
    exp_residue_participation_rank: float
    max_binary_cut_entropy_bits: float | None
    mean_binary_cut_entropy_bits: float | None
    binary_cut_count: int
    sequential_steps_to_period: int | None
    gcd_defect_min_d: int | None
    gcd_defect_min_factor: int | None
    orbit_support: float
    valuation_active_modes: float
    valuation_active_diff_modes: float
    valuation_mass: float
    valuation_diff_mass: float
    valuation_mode_entropy: float
    valuation_support_over_period: float
    valuation_active_modes_over_logN: float
    sparse_gcd_success_rate: float
    sparse_gcd_mean_samples_if_success: float
    sparse_gcd_median_samples_if_success: float
    orbit_graph_nodes: float
    orbit_graph_rank: float
    orbit_graph_spectral_entropy: float
    orbit_graph_rank_over_period: float
    valuation_matrix_rank: float
    valuation_matrix_stable_rank: float
    valuation_matrix_spectral_entropy: float
    valuation_matrix_spectral_entropy_norm: float
    valuation_rank_over_period: float
    valuation_rank_over_support: float
    diff_valuation_matrix_rank: float
    diff_valuation_matrix_stable_rank: float
    diff_valuation_rank_over_period: float
    valuation_hankel_window: float
    valuation_hankel_rank: float
    valuation_hankel_stable_rank: float
    valuation_hankel_spectral_entropy: float
    valuation_hankel_spectral_entropy_norm: float
    valuation_hankel_rank_over_period: float
    diff_valuation_hankel_rank: float
    diff_valuation_hankel_stable_rank: float
    diff_valuation_hankel_rank_over_period: float
    successor_permutation_rank: float
    successor_permutation_stable_rank: float
    successor_permutation_spectral_entropy: float
    successor_permutation_spectral_entropy_norm: float
    successor_rank_over_period: float
    cycle_adjacency_rank: float
    cycle_adjacency_spectral_entropy: float
    cycle_adjacency_spectral_entropy_norm: float
    cycle_adjacency_rank_over_period: float
    value_compresses_but_order_not: float
    order_compresses_in_valuation_basis: float


def benchmark_case(n: int, a: int, args: argparse.Namespace) -> ShorBenchRow:
    fac = factor_trial(n)
    p, qfac = fac if fac else (None, None)
    qdim = next_power_of_two_at_least(n * n) if args.Q is None else int(args.Q)
    r = multiplicative_order(a, n)
    t0 = time.time()
    M, residues = shor_state_matrix(n, a, qdim)
    sch_rank, ent, part_rank, _ = schmidt_entropy_from_matrix(M)
    exact_seconds = time.time() - t0
    qbits_x = int(math.log2(qdim))
    ydim = next_power_of_two_at_least(n)
    qbits_y = int(math.log2(ydim))
    total_q = qbits_x + qbits_y
    cuts = binary_cut_entropies(n, a, qdim, max_total_qubits=args.max_binary_qubits)
    if cuts is None:
        max_cut = None
        mean_cut = None
        cut_count = 0
    else:
        max_cut = float(max(cuts)) if cuts else 0.0
        mean_cut = float(np.mean(cuts)) if cuts else 0.0
        cut_count = len(cuts)
    ok, fp, fq = shor_factor_from_period(a, n, r)
    seq_r, seq_steps = sequential_period_search(n, a)
    gcd_cands = gcd_defect_period_candidates(n, a, max_d=min(n, max(2 * (r or 1), 1)))
    if gcd_cands:
        gd, gf = gcd_cands[0]
    else:
        gd, gf = None, None
    if r is not None:
        val_feats = valuation_orbit_features(n, a, r, n_primes=args.n_primes)
        sparse = random_sparse_gcd_recovery(n, a, r, trials=args.sparse_trials, max_samples=args.sparse_max_samples, seed=args.seed + n + a)
        graph = orbit_graph_spectrum(n, a, r)
        graph["orbit_graph_rank_over_period"] = float(graph.get("orbit_graph_rank", float("nan")) / max(r, 1))
        order = cycle_order_complexity(n, a, r, n_primes=args.n_primes, hankel_window=args.hankel_window)
    else:
        val_feats = {
            "orbit_support": float("nan"),
            "valuation_active_modes": float("nan"),
            "valuation_active_diff_modes": float("nan"),
            "valuation_mass": float("nan"),
            "valuation_diff_mass": float("nan"),
            "valuation_mode_entropy": float("nan"),
            "valuation_support_over_period": float("nan"),
            "valuation_active_modes_over_logN": float("nan"),
        }
        sparse = {
            "sparse_gcd_success_rate": float("nan"),
            "sparse_gcd_mean_samples_if_success": float("nan"),
            "sparse_gcd_median_samples_if_success": float("nan"),
        }
        graph = {"orbit_graph_nodes": float("nan"), "orbit_graph_rank": float("nan"), "orbit_graph_spectral_entropy": float("nan"), "orbit_graph_rank_over_period": float("nan")}
        order = {
            "valuation_matrix_rank": float("nan"),
            "valuation_matrix_stable_rank": float("nan"),
            "valuation_matrix_spectral_entropy": float("nan"),
            "valuation_matrix_spectral_entropy_norm": float("nan"),
            "valuation_rank_over_period": float("nan"),
            "valuation_rank_over_support": float("nan"),
            "diff_valuation_matrix_rank": float("nan"),
            "diff_valuation_matrix_stable_rank": float("nan"),
            "diff_valuation_rank_over_period": float("nan"),
            "valuation_hankel_window": float("nan"),
            "valuation_hankel_rank": float("nan"),
            "valuation_hankel_stable_rank": float("nan"),
            "valuation_hankel_spectral_entropy": float("nan"),
            "valuation_hankel_spectral_entropy_norm": float("nan"),
            "valuation_hankel_rank_over_period": float("nan"),
            "diff_valuation_hankel_rank": float("nan"),
            "diff_valuation_hankel_stable_rank": float("nan"),
            "diff_valuation_hankel_rank_over_period": float("nan"),
            "successor_permutation_rank": float("nan"),
            "successor_permutation_stable_rank": float("nan"),
            "successor_permutation_spectral_entropy": float("nan"),
            "successor_permutation_spectral_entropy_norm": float("nan"),
            "successor_rank_over_period": float("nan"),
            "cycle_adjacency_rank": float("nan"),
            "cycle_adjacency_spectral_entropy": float("nan"),
            "cycle_adjacency_spectral_entropy_norm": float("nan"),
            "cycle_adjacency_rank_over_period": float("nan"),
            "value_compresses_but_order_not": float("nan"),
            "order_compresses_in_valuation_basis": float("nan"),
        }
    return ShorBenchRow(
        N=n,
        p=p,
        q_factor=qfac,
        a=a,
        Q=qdim,
        n_bits=int(math.ceil(math.log2(n))),
        qbits_x=qbits_x,
        qbits_y=qbits_y,
        total_qubits_padded=total_q,
        period_r=r,
        shor_factor_success=ok,
        factor_p_from_period=fp,
        factor_q_from_period=fq,
        exact_state_seconds=exact_seconds,
        exp_residue_schmidt_rank=sch_rank,
        exp_residue_entropy_bits=ent,
        exp_residue_participation_rank=part_rank,
        max_binary_cut_entropy_bits=max_cut,
        mean_binary_cut_entropy_bits=mean_cut,
        binary_cut_count=cut_count,
        sequential_steps_to_period=seq_steps if seq_r is not None else None,
        gcd_defect_min_d=gd,
        gcd_defect_min_factor=gf,
        **val_feats,
        **sparse,
        **graph,
        **order,
    )


def parse_list_int(s: str) -> list[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--Ns", type=str, default="15,21,33,35,39,55,65,77,91,143,187,221,247,299,323,437,667,899")
    ap.add_argument("--max-N", type=int, default=None, help="Use semiprimes up to max-N instead of --Ns.")
    ap.add_argument("--max-cases", type=int, default=80)
    ap.add_argument("--max-bases", type=int, default=4)
    ap.add_argument("--Q", type=int, default=None, help="Override exponent register dimension. Default next pow2 >= N^2.")
    ap.add_argument("--max-binary-qubits", type=int, default=18)
    ap.add_argument("--n-primes", type=int, default=16)
    ap.add_argument("--hankel-window", type=int, default=None, help="Block-Hankel window for valuation sequence order-complexity. Default min(8, r//2).")
    ap.add_argument("--sparse-trials", type=int, default=200)
    ap.add_argument("--sparse-max-samples", type=int, default=64)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=str, default="/mnt/data/shor_prime_fock_benchmark.csv")
    args = ap.parse_args()

    if args.max_N is not None:
        Ns = semiprimes(args.max_N)
    else:
        Ns = parse_list_int(args.Ns)
    cases: list[tuple[int, int]] = []
    for n in Ns:
        if factor_trial(n) is None:
            print(f"[SKIP] N={n} is not composite/semiprime by trial division")
            continue
        for a in valid_bases(n, max_bases=args.max_bases):
            cases.append((n, a))
            if len(cases) >= args.max_cases:
                break
        if len(cases) >= args.max_cases:
            break

    print(f"[CASES] {len(cases)} cases")
    rows = []
    for i, (n, a) in enumerate(cases, 1):
        row = benchmark_case(n, a, args)
        rows.append(row)
        print(
            f"[{i:03d}/{len(cases)}] N={n} a={a} r={row.period_r} "
            f"success={row.shor_factor_success} S_exp|res={row.exp_residue_entropy_bits:.3f} "
            f"rank={row.exp_residue_schmidt_rank} pf_modes={row.valuation_active_modes:.0f} "
            f"v_rank={row.valuation_matrix_rank:.0f} h_rank={row.valuation_hankel_rank:.0f} "
            f"succ_rank={row.successor_permutation_rank:.0f} sparse_success={row.sparse_gcd_success_rate:.2f}"
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()) if rows else [])
        if rows:
            writer.writeheader()
            for row in rows:
                writer.writerow(asdict(row))

    if rows:
        arr_r = np.array([row.period_r or np.nan for row in rows], dtype=float)
        arr_ent = np.array([row.exp_residue_entropy_bits for row in rows], dtype=float)
        arr_rank = np.array([row.exp_residue_schmidt_rank for row in rows], dtype=float)
        arr_success = np.array([row.sparse_gcd_success_rate for row in rows], dtype=float)
        arr_vrank_ratio = np.array([row.valuation_rank_over_period for row in rows], dtype=float)
        arr_hrank_ratio = np.array([row.valuation_hankel_rank_over_period for row in rows], dtype=float)
        arr_succ_ratio = np.array([row.successor_rank_over_period for row in rows], dtype=float)
        arr_val_comp_order_not = np.array([row.value_compresses_but_order_not for row in rows], dtype=float)
        arr_order_comp = np.array([row.order_compresses_in_valuation_basis for row in rows], dtype=float)
        shor_success = np.mean([row.shor_factor_success for row in rows])
        print("\n[SUMMARY]")
        print(f"rows: {len(rows)}")
        print(f"period r mean/median/max: {np.nanmean(arr_r):.2f} / {np.nanmedian(arr_r):.2f} / {np.nanmax(arr_r):.0f}")
        print(f"exp|res entropy mean/max bits: {np.mean(arr_ent):.3f} / {np.max(arr_ent):.3f}")
        print(f"exp|res Schmidt rank mean/max: {np.mean(arr_rank):.2f} / {np.max(arr_rank):.0f}")
        print(f"Shor factor success from exact r: {shor_success:.3f}")
        print(f"sparse gcd recovery success mean: {np.nanmean(arr_success):.3f}")
        print(f"valuation rank / period mean: {np.nanmean(arr_vrank_ratio):.3f}")
        print(f"valuation Hankel rank / period mean: {np.nanmean(arr_hrank_ratio):.3f}")
        print(f"successor rank / period mean: {np.nanmean(arr_succ_ratio):.3f}")
        print(f"value-compresses-but-order-not fraction: {np.nanmean(arr_val_comp_order_not):.3f}")
        print(f"order-compresses-in-valuation-basis fraction: {np.nanmean(arr_order_comp):.3f}")
        print(f"[SAVE] {out}")
    else:
        print("[NO ROWS]")


if __name__ == "__main__":
    main()
