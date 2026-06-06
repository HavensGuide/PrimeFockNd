#!/usr/bin/env python3
"""
Prime-Fock / Shor v4: wavefunction + QFT interference benchmark.

This benchmark tests the hypothesis that any Prime-Fock advantage for Shor-like
period finding should appear in the interference / wave dynamics, not merely in
static residue valuation compression.

For each semiprime N and base a:
  1. Compute the true period r of a mod N.
  2. Build the exact modular-exponentiation wavefunction implicitly:
         |psi> = 1/sqrt(Q) sum_x |x>|a^x mod N>
  3. Apply the QFT over the exponent register analytically/implicitly.
  4. Measure the exact QFT peak distribution and period recovery.
  5. Build Prime-Fock valuation wave signals over the same exponent lattice:
         V_x = valuation(a^x mod N)
     and compute the QFT energy of this compressed arithmetic signal.
  6. Compare whether Prime-Fock interference peaks align with exact Shor peaks,
     recover the same period, and concentrate signal more than shuffled controls.

This is a diagnostic benchmark, not a scalable factoring claim.
"""
from __future__ import annotations

import argparse
import csv
import math
import random
import time
from dataclasses import dataclass, asdict, fields
from fractions import Fraction
from pathlib import Path
from typing import Optional

import numpy as np


# ------------------------- number theory -------------------------

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
    y = 1
    for r in range(1, max_steps + 1):
        y = (y * a) % n
        if y == 1:
            return r
    return None


def valid_bases(n: int, max_bases: int = 8) -> list[int]:
    bases = []
    for a in range(2, n):
        if gcd(a, n) != 1:
            continue
        r = multiplicative_order(a, n)
        if r is not None and r > 1:
            bases.append(a)
        if len(bases) >= max_bases:
            break
    return bases


def shor_factor_from_period(a: int, n: int, r: Optional[int]) -> tuple[bool, int | None, int | None]:
    if r is None or r % 2 == 1:
        return False, None, None
    y = pow(a, r // 2, n)
    if y in (1, n - 1):
        return False, None, None
    p = gcd(y - 1, n)
    q = gcd(y + 1, n)
    ok = 1 < p < n and 1 < q < n and p * q == n
    return ok, p if 1 < p < n else None, q if 1 < q < n else None


def next_power_of_two_at_least(n: int) -> int:
    return 1 << ((n - 1).bit_length())


def choose_Q(n: int, q_mode: str, q_bits: Optional[int], max_q: int) -> int:
    if q_bits is not None:
        q = 1 << int(q_bits)
    elif q_mode == "n2":
        q = next_power_of_two_at_least(n * n)
    elif q_mode == "2n2":
        q = next_power_of_two_at_least(2 * n * n)
    elif q_mode == "n":
        q = next_power_of_two_at_least(n)
    else:
        raise ValueError(f"unknown q_mode {q_mode}")
    return min(q, max_q)


# ------------------------- Prime-Fock valuations -------------------------

def first_primes(k: int) -> list[int]:
    out = []
    x = 2
    while len(out) < k:
        if is_prime(x):
            out.append(x)
        x += 1
    return out


def valuation_vector(x: int, primes: list[int]) -> np.ndarray:
    x = abs(int(x))
    v = np.zeros(len(primes), dtype=float)
    if x <= 1:
        return v
    for i, p in enumerate(primes):
        while x % p == 0:
            v[i] += 1.0
            x //= p
        if x == 1:
            break
    return v


def orbit_residues(n: int, a: int, r: int) -> np.ndarray:
    return np.fromiter((pow(a, t, n) for t in range(r)), dtype=np.int64, count=r)


def valuation_sequence(n: int, a: int, r: int, n_primes: int) -> tuple[np.ndarray, np.ndarray, list[int]]:
    primes = first_primes(n_primes)
    residues = orbit_residues(n, a, r)
    V = np.stack([valuation_vector(int(y), primes) for y in residues], axis=0)
    return residues, V, primes


# ------------------------- wave / interference utilities -------------------------

def normalize_prob(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    p = np.maximum(p, 0.0)
    s = float(p.sum())
    if s <= 0:
        return np.ones_like(p) / max(len(p), 1)
    return p / s


def shor_qft_prob_by_period(Q: int, r: int) -> np.ndarray:
    """Exact QFT probability over exponent frequencies k, without dense QxN state.

    Since a^x mod N has period r and distinct residues for t=0..r-1, the residue
    register splits the sum into r congruence classes x == t mod r.  The QFT
    probability is:
        P(k) = (1/Q^2) sum_t |sum_{x<Q, x=t mod r} exp(2pi i k x/Q)|^2.
    """
    k = np.arange(Q, dtype=float)
    prob = np.zeros(Q, dtype=float)
    for t in range(r):
        count = (Q - 1 - t) // r + 1 if t < Q else 0
        if count <= 0:
            continue
        # Sum_m exp(2pi i k(t + m r)/Q).
        ratio_phase = 2.0 * np.pi * k * r / Q
        numer = np.sin(0.5 * count * ratio_phase)
        denom = np.sin(0.5 * ratio_phase)
        # Limit when denom -> 0 is count.
        mag = np.empty_like(k)
        mask = np.abs(denom) < 1e-12
        mag[mask] = count
        mag[~mask] = np.abs(numer[~mask] / denom[~mask])
        prob += mag * mag
    prob /= float(Q * Q)
    return normalize_prob(prob)


def repeat_periodic_rows(X: np.ndarray, Q: int) -> np.ndarray:
    idx = np.arange(Q) % X.shape[0]
    return X[idx]


def valuation_qft_energy(V: np.ndarray, Q: int, *, center: bool = True, whiten: bool = False) -> np.ndarray:
    """QFT energy of repeated Prime-Fock valuation sequence V[x mod r]."""
    X = repeat_periodic_rows(V, Q).astype(float)
    if center:
        X = X - X.mean(axis=0, keepdims=True)
    active = np.sum(np.abs(X), axis=0) > 1e-12
    X = X[:, active]
    if X.shape[1] == 0:
        return np.ones(Q) / Q
    if whiten:
        sd = X.std(axis=0, keepdims=True)
        X = X / np.where(sd > 1e-12, sd, 1.0)
    F = np.fft.fft(X, axis=0)
    e = np.sum(np.abs(F) ** 2, axis=1)
    return normalize_prob(e)


def shuffled_valuation_energy(V: np.ndarray, Q: int, seed: int, *, trials: int = 8) -> np.ndarray:
    rng = np.random.default_rng(seed)
    acc = np.zeros(Q, dtype=float)
    for _ in range(max(1, trials)):
        Vs = V.copy()
        rng.shuffle(Vs, axis=0)
        acc += valuation_qft_energy(Vs, Q, center=True, whiten=False)
    return normalize_prob(acc)


def entropy_bits(p: np.ndarray) -> float:
    p = normalize_prob(p)
    nz = p[p > 1e-15]
    return float(-np.sum(nz * np.log2(nz)))


def topk_indices(p: np.ndarray, k: int) -> np.ndarray:
    k = int(min(max(k, 1), len(p)))
    idx = np.argpartition(-p, kth=k-1)[:k]
    return idx[np.argsort(-p[idx])]


def ideal_peak_indices(Q: int, r: int) -> np.ndarray:
    return np.unique(np.array([int(round(j * Q / r)) % Q for j in range(r)], dtype=int))


def mass_near_indices(p: np.ndarray, centers: np.ndarray, radius: int) -> float:
    Q = len(p)
    mask = np.zeros(Q, dtype=bool)
    for c in centers:
        for d in range(-radius, radius + 1):
            mask[(int(c) + d) % Q] = True
    return float(np.sum(p[mask]))


def circular_peak_distance(k: int, centers: np.ndarray, Q: int) -> int:
    diffs = np.abs(((centers - int(k) + Q // 2) % Q) - Q // 2)
    return int(np.min(diffs)) if centers.size else Q


def distribution_corr(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if np.std(a) <= 1e-15 or np.std(b) <= 1e-15:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def period_from_frequency(k: int, Q: int, N: int, max_den: Optional[int] = None) -> int | None:
    if k <= 0 or k >= Q:
        return None
    if max_den is None:
        max_den = N
    frac = Fraction(int(k), int(Q)).limit_denominator(max_den)
    den = int(frac.denominator)
    if den <= 1:
        return None
    return den


def recover_period_from_distribution(p: np.ndarray, N: int, a: int, top_k: int) -> tuple[int | None, bool, int | None, int | None, int | None]:
    Q = len(p)
    candidates = []
    for k in topk_indices(p, top_k):
        cand = period_from_frequency(int(k), Q, N, max_den=N)
        if cand is None:
            continue
        # Try small multiples of denominator; continued fractions often return r/g.
        for m in range(1, min(8, N // max(cand, 1)) + 1):
            rr = cand * m
            if rr > 1 and pow(a, rr, N) == 1:
                ok, fp, fq = shor_factor_from_period(a, N, rr)
                return rr, ok, fp, fq, int(k)
        candidates.append(cand)
    if candidates:
        return candidates[0], False, None, None, None
    return None, False, None, None, None


def qft_peak_metrics(p: np.ndarray, r: int, *, top_k: int) -> dict[str, float]:
    Q = len(p)
    peaks = ideal_peak_indices(Q, r)
    top = topk_indices(p, top_k)
    peak_set = set(int(x) for x in peaks.tolist())
    top_hit_exact = sum(1 for x in top if int(x) in peak_set) / max(len(top), 1)
    dists = [circular_peak_distance(int(x), peaks, Q) for x in top]
    radius = max(1, int(round(Q / max(8 * r, 1))))
    return {
        "entropy_bits": entropy_bits(p),
        "entropy_norm": entropy_bits(p) / math.log2(Q) if Q > 1 else 0.0,
        "max_prob": float(np.max(p)),
        "top1_index": float(top[0]),
        "top1_prob": float(p[top[0]]),
        "topk_exact_peak_hit_rate": float(top_hit_exact),
        "topk_mean_distance_to_ideal_peak": float(np.mean(dists)) if dists else float("nan"),
        "ideal_peak_mass_exact": float(np.sum(p[peaks])),
        "ideal_peak_mass_radius": mass_near_indices(p, peaks, radius),
        "ideal_peak_radius": float(radius),
    }


def valuation_overlap_graph(V: np.ndarray) -> np.ndarray:
    """Simple arithmetic overlap graph using min valuation overlap."""
    V = np.asarray(V, dtype=float)
    r = V.shape[0]
    W = np.zeros((r, r), dtype=float)
    for i in range(r):
        W[i] = np.minimum(V[i][None, :], V).sum(axis=1)
    np.fill_diagonal(W, 0.0)
    return W


def phase_current_metrics(V: np.ndarray, Q: int, r: int, k: int) -> dict[str, float]:
    """Graph/wave current on valuation-overlap graph at Fourier frequency k."""
    W = valuation_overlap_graph(V)
    if np.sum(W) <= 1e-15:
        return {
            "phase_current_abs_mean": 0.0,
            "phase_current_net_abs": 0.0,
            "phase_current_cycle_alignment": 0.0,
            "phase_current_graph_energy": 0.0,
        }
    t = np.arange(r, dtype=float)
    theta = 2.0 * np.pi * float(k) * t / float(Q)
    dtheta = theta[None, :] - theta[:, None]
    J = W * np.sin(dtheta)
    net = J.sum(axis=1)
    # Cycle edge current t->t+1; positive coherent flow means wave phase agrees with modular successor.
    cycle_curr = []
    for i in range(r):
        j = (i + 1) % r
        cycle_curr.append(W[i, j] * math.sin(theta[j] - theta[i]))
    cycle_curr = np.array(cycle_curr, dtype=float)
    return {
        "phase_current_abs_mean": float(np.mean(np.abs(J))),
        "phase_current_net_abs": float(np.sum(np.abs(net)) / (np.sum(W) + 1e-12)),
        "phase_current_cycle_alignment": float(np.mean(cycle_curr)),
        "phase_current_graph_energy": float(np.mean(J * J)),
    }


def compare_pf_to_exact(exact: np.ndarray, pf: np.ndarray, r: int, top_k: int) -> dict[str, float]:
    Q = len(exact)
    exact_top = set(int(x) for x in topk_indices(exact, top_k).tolist())
    pf_top = set(int(x) for x in topk_indices(pf, top_k).tolist())
    inter = len(exact_top & pf_top)
    union = len(exact_top | pf_top)
    peaks = ideal_peak_indices(Q, r)
    return {
        "pf_exact_distribution_corr": distribution_corr(exact, pf),
        "pf_exact_topk_jaccard": float(inter / max(union, 1)),
        "pf_topk_ideal_peak_mass_exact": float(sum(1 for x in pf_top if x in set(peaks.tolist())) / max(len(pf_top), 1)),
        "pf_energy_on_exact_topk_mass": float(np.sum(pf[list(exact_top)])) if exact_top else 0.0,
        "exact_energy_on_pf_topk_mass": float(np.sum(exact[list(pf_top)])) if pf_top else 0.0,
    }


# ------------------------- rows / benchmark -------------------------

@dataclass
class WaveBenchRow:
    N: int
    p: int | None
    q_factor: int | None
    a: int
    period_r: int
    n_bits: int
    Q: int
    q_over_n2: float
    shor_factor_success_from_true_r: bool
    factor_p_from_true_r: int | None
    factor_q_from_true_r: int | None
    elapsed_seconds: float

    valuation_active_modes: float
    valuation_rank: float
    valuation_rank_over_period: float

    exact_qft_entropy_bits: float
    exact_qft_entropy_norm: float
    exact_qft_max_prob: float
    exact_qft_top1_index: float
    exact_qft_top1_prob: float
    exact_qft_topk_exact_peak_hit_rate: float
    exact_qft_topk_mean_distance_to_ideal_peak: float
    exact_qft_ideal_peak_mass_exact: float
    exact_qft_ideal_peak_mass_radius: float
    exact_period_guess: int | None
    exact_period_guess_success: bool
    exact_factor_success_from_qft: bool
    exact_qft_witness_k: int | None

    pf_qft_entropy_bits: float
    pf_qft_entropy_norm: float
    pf_qft_max_prob: float
    pf_qft_top1_index: float
    pf_qft_top1_prob: float
    pf_qft_topk_exact_peak_hit_rate: float
    pf_qft_topk_mean_distance_to_ideal_peak: float
    pf_qft_ideal_peak_mass_exact: float
    pf_qft_ideal_peak_mass_radius: float
    pf_period_guess: int | None
    pf_period_guess_success: bool
    pf_factor_success_from_qft: bool
    pf_qft_witness_k: int | None

    pf_exact_distribution_corr: float
    pf_exact_topk_jaccard: float
    pf_topk_ideal_peak_mass_exact: float
    pf_energy_on_exact_topk_mass: float
    exact_energy_on_pf_topk_mass: float

    shuffled_pf_distribution_corr: float
    shuffled_pf_ideal_peak_mass_radius: float
    pf_vs_shuffled_peak_mass_gain: float
    pf_vs_shuffled_corr_gain: float

    phase_current_abs_mean_at_exact_top1: float
    phase_current_net_abs_at_exact_top1: float
    phase_current_cycle_alignment_at_exact_top1: float
    phase_current_graph_energy_at_exact_top1: float
    phase_current_abs_mean_at_pf_top1: float
    phase_current_net_abs_at_pf_top1: float
    phase_current_cycle_alignment_at_pf_top1: float
    phase_current_graph_energy_at_pf_top1: float

    interference_alignment_win: float
    pf_period_recovery_win: float
    pf_peak_concentration_win: float
    value_compressed_but_interference_not: float


def matrix_rank(X: np.ndarray, tol: float = 1e-10) -> float:
    if X.size == 0:
        return 0.0
    s = np.linalg.svd(np.asarray(X, dtype=float), compute_uv=False)
    return float(np.sum(s > tol))


def benchmark_case(
    N: int,
    a: int,
    *,
    n_primes: int,
    q_mode: str,
    q_bits: Optional[int],
    max_q: int,
    top_k: int,
    shuffle_trials: int,
    seed: int,
) -> WaveBenchRow | None:
    t0 = time.time()
    r = multiplicative_order(a, N)
    if r is None or r <= 1:
        return None
    Q = choose_Q(N, q_mode, q_bits, max_q)
    if Q < 2 * r:
        # Too little Fourier resolution for meaningful diagnostics.
        return None
    fac = factor_trial(N)
    p, q = fac if fac is not None else (None, None)
    shor_ok, fp, fq = shor_factor_from_period(a, N, r)

    residues, V, primes = valuation_sequence(N, a, r, n_primes=n_primes)
    active_modes = float(np.sum(np.sum(np.abs(V), axis=0) > 0))
    v_rank = matrix_rank(V)

    exact_p = shor_qft_prob_by_period(Q, r)
    pf_p = valuation_qft_energy(V, Q, center=True, whiten=False)
    shuf_p = shuffled_valuation_energy(V, Q, seed=seed + 1009 * N + 9176 * a, trials=shuffle_trials)

    exact_m = qft_peak_metrics(exact_p, r, top_k=top_k)
    pf_m = qft_peak_metrics(pf_p, r, top_k=top_k)
    cmp_m = compare_pf_to_exact(exact_p, pf_p, r, top_k=top_k)
    shuf_m = qft_peak_metrics(shuf_p, r, top_k=top_k)

    exact_guess, exact_guess_factor_ok, efp, efq, exact_k = recover_period_from_distribution(exact_p, N, a, top_k=top_k)
    pf_guess, pf_guess_factor_ok, pfp, pfq, pf_k = recover_period_from_distribution(pf_p, N, a, top_k=top_k)

    exact_top1 = int(exact_m["top1_index"])
    pf_top1 = int(pf_m["top1_index"])
    cur_exact = phase_current_metrics(V, Q, r, exact_top1)
    cur_pf = phase_current_metrics(V, Q, r, pf_top1)

    pf_peak_gain = pf_m["ideal_peak_mass_radius"] - shuf_m["ideal_peak_mass_radius"]
    corr_shuf = distribution_corr(exact_p, shuf_p)
    corr_gain = cmp_m["pf_exact_distribution_corr"] - corr_shuf if not math.isnan(corr_shuf) else float("nan")

    interference_alignment_win = float(
        (cmp_m["pf_exact_distribution_corr"] > 0.25 if not math.isnan(cmp_m["pf_exact_distribution_corr"]) else False)
        and cmp_m["pf_exact_topk_jaccard"] > 0.10
    )
    pf_period_recovery_win = float(pf_guess == r or (pf_guess is not None and r % pf_guess == 0 and pow(a, pf_guess, N) == 1))
    pf_peak_concentration_win = float(pf_peak_gain > 0.05 and pf_m["ideal_peak_mass_radius"] > shuf_m["ideal_peak_mass_radius"] * 1.25)
    value_compressed_but_interference_not = float((v_rank / max(r, 1) < 0.5) and not (interference_alignment_win or pf_period_recovery_win))

    return WaveBenchRow(
        N=N,
        p=p,
        q_factor=q,
        a=a,
        period_r=r,
        n_bits=N.bit_length(),
        Q=Q,
        q_over_n2=float(Q / (N * N)),
        shor_factor_success_from_true_r=bool(shor_ok),
        factor_p_from_true_r=fp,
        factor_q_from_true_r=fq,
        elapsed_seconds=time.time() - t0,
        valuation_active_modes=active_modes,
        valuation_rank=v_rank,
        valuation_rank_over_period=v_rank / max(r, 1),
        exact_qft_entropy_bits=exact_m["entropy_bits"],
        exact_qft_entropy_norm=exact_m["entropy_norm"],
        exact_qft_max_prob=exact_m["max_prob"],
        exact_qft_top1_index=exact_m["top1_index"],
        exact_qft_top1_prob=exact_m["top1_prob"],
        exact_qft_topk_exact_peak_hit_rate=exact_m["topk_exact_peak_hit_rate"],
        exact_qft_topk_mean_distance_to_ideal_peak=exact_m["topk_mean_distance_to_ideal_peak"],
        exact_qft_ideal_peak_mass_exact=exact_m["ideal_peak_mass_exact"],
        exact_qft_ideal_peak_mass_radius=exact_m["ideal_peak_mass_radius"],
        exact_period_guess=exact_guess,
        exact_period_guess_success=bool(exact_guess == r or (exact_guess is not None and r % exact_guess == 0 and pow(a, exact_guess, N) == 1)),
        exact_factor_success_from_qft=bool(exact_guess_factor_ok),
        exact_qft_witness_k=exact_k,
        pf_qft_entropy_bits=pf_m["entropy_bits"],
        pf_qft_entropy_norm=pf_m["entropy_norm"],
        pf_qft_max_prob=pf_m["max_prob"],
        pf_qft_top1_index=pf_m["top1_index"],
        pf_qft_top1_prob=pf_m["top1_prob"],
        pf_qft_topk_exact_peak_hit_rate=pf_m["topk_exact_peak_hit_rate"],
        pf_qft_topk_mean_distance_to_ideal_peak=pf_m["topk_mean_distance_to_ideal_peak"],
        pf_qft_ideal_peak_mass_exact=pf_m["ideal_peak_mass_exact"],
        pf_qft_ideal_peak_mass_radius=pf_m["ideal_peak_mass_radius"],
        pf_period_guess=pf_guess,
        pf_period_guess_success=bool(pf_guess == r or (pf_guess is not None and r % pf_guess == 0 and pow(a, pf_guess, N) == 1)),
        pf_factor_success_from_qft=bool(pf_guess_factor_ok),
        pf_qft_witness_k=pf_k,
        pf_exact_distribution_corr=cmp_m["pf_exact_distribution_corr"],
        pf_exact_topk_jaccard=cmp_m["pf_exact_topk_jaccard"],
        pf_topk_ideal_peak_mass_exact=cmp_m["pf_topk_ideal_peak_mass_exact"],
        pf_energy_on_exact_topk_mass=cmp_m["pf_energy_on_exact_topk_mass"],
        exact_energy_on_pf_topk_mass=cmp_m["exact_energy_on_pf_topk_mass"],
        shuffled_pf_distribution_corr=corr_shuf,
        shuffled_pf_ideal_peak_mass_radius=shuf_m["ideal_peak_mass_radius"],
        pf_vs_shuffled_peak_mass_gain=pf_peak_gain,
        pf_vs_shuffled_corr_gain=corr_gain,
        phase_current_abs_mean_at_exact_top1=cur_exact["phase_current_abs_mean"],
        phase_current_net_abs_at_exact_top1=cur_exact["phase_current_net_abs"],
        phase_current_cycle_alignment_at_exact_top1=cur_exact["phase_current_cycle_alignment"],
        phase_current_graph_energy_at_exact_top1=cur_exact["phase_current_graph_energy"],
        phase_current_abs_mean_at_pf_top1=cur_pf["phase_current_abs_mean"],
        phase_current_net_abs_at_pf_top1=cur_pf["phase_current_net_abs"],
        phase_current_cycle_alignment_at_pf_top1=cur_pf["phase_current_cycle_alignment"],
        phase_current_graph_energy_at_pf_top1=cur_pf["phase_current_graph_energy"],
        interference_alignment_win=interference_alignment_win,
        pf_period_recovery_win=pf_period_recovery_win,
        pf_peak_concentration_win=pf_peak_concentration_win,
        value_compressed_but_interference_not=value_compressed_but_interference_not,
    )


def build_cases(args: argparse.Namespace) -> list[tuple[int, int]]:
    if args.Ns:
        Ns = [int(x.strip()) for x in args.Ns.split(",") if x.strip()]
    else:
        Ns = semiprimes(args.max_N, min_n=args.min_N)
    cases = []
    for N in Ns:
        for a in valid_bases(N, max_bases=args.max_bases):
            cases.append((N, a))
    if args.max_cases is not None:
        cases = cases[: args.max_cases]
    return cases


def write_rows(rows: list[WaveBenchRow], out: str | Path) -> None:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    names = [f.name for f in fields(WaveBenchRow)]
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=names)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def print_summary(rows: list[WaveBenchRow]) -> None:
    if not rows:
        print("[SUMMARY] no rows")
        return
    def arr(name: str) -> np.ndarray:
        return np.array([getattr(r, name) for r in rows], dtype=float)
    print("\n[SUMMARY]")
    print(f"rows: {len(rows)}")
    print(f"period r mean/median/max: {np.mean(arr('period_r')):.3f} / {np.median(arr('period_r')):.3f} / {np.max(arr('period_r')):.0f}")
    print(f"Q mean/median/max: {np.mean(arr('Q')):.1f} / {np.median(arr('Q')):.1f} / {np.max(arr('Q')):.0f}")
    print(f"valuation rank / period mean: {np.mean(arr('valuation_rank_over_period')):.3f}")
    print(f"exact QFT period recovery rate: {np.mean(arr('exact_period_guess_success')):.3f}")
    print(f"PF QFT period recovery rate: {np.mean(arr('pf_period_guess_success')):.3f}")
    print(f"PF/exact distribution corr mean: {np.nanmean(arr('pf_exact_distribution_corr')):.3f}")
    print(f"PF/exact top-k Jaccard mean: {np.mean(arr('pf_exact_topk_jaccard')):.3f}")
    print(f"PF ideal peak mass radius mean: {np.mean(arr('pf_qft_ideal_peak_mass_radius')):.3f}")
    print(f"shuffled ideal peak mass radius mean: {np.mean(arr('shuffled_pf_ideal_peak_mass_radius')):.3f}")
    print(f"PF vs shuffled peak-mass gain mean: {np.mean(arr('pf_vs_shuffled_peak_mass_gain')):.3f}")
    print(f"interference alignment win fraction: {np.mean(arr('interference_alignment_win')):.3f}")
    print(f"PF peak concentration win fraction: {np.mean(arr('pf_peak_concentration_win')):.3f}")
    print(f"value compressed but interference not fraction: {np.mean(arr('value_compressed_but_interference_not')):.3f}")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Prime-Fock Shor v4 wave/QFT interference benchmark")
    ap.add_argument("--Ns", type=str, default=None, help="Comma-separated explicit N values")
    ap.add_argument("--min-N", type=int, default=15)
    ap.add_argument("--max-N", type=int, default=300)
    ap.add_argument("--max-bases", type=int, default=4)
    ap.add_argument("--max-cases", type=int, default=None)
    ap.add_argument("--n-primes", type=int, default=64)
    ap.add_argument("--q-mode", type=str, default="n2", choices=["n", "n2", "2n2"])
    ap.add_argument("--q-bits", type=int, default=None, help="Override Q as 2^q_bits")
    ap.add_argument("--max-q", type=int, default=262144, help="Cap Fourier register size for laptop safety")
    ap.add_argument("--top-k", type=int, default=16)
    ap.add_argument("--shuffle-trials", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=str, default="shor_prime_fock_wave_interference_v4.csv")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    cases = build_cases(args)
    print(f"[CASES] {len(cases)} cases")
    rows: list[WaveBenchRow] = []
    for i, (N, a) in enumerate(cases, 1):
        try:
            row = benchmark_case(
                N,
                a,
                n_primes=args.n_primes,
                q_mode=args.q_mode,
                q_bits=args.q_bits,
                max_q=args.max_q,
                top_k=args.top_k,
                shuffle_trials=args.shuffle_trials,
                seed=args.seed,
            )
        except Exception as e:
            print(f"[{i:03d}/{len(cases)}] N={N} a={a} ERROR: {e}")
            continue
        if row is None:
            print(f"[{i:03d}/{len(cases)}] N={N} a={a} skipped")
            continue
        rows.append(row)
        print(
            f"[{i:03d}/{len(cases)}] N={N} a={a} r={row.period_r} Q={row.Q} "
            f"exact_rec={int(row.exact_period_guess_success)} pf_rec={int(row.pf_period_guess_success)} "
            f"pf_corr={row.pf_exact_distribution_corr:.3f} topJ={row.pf_exact_topk_jaccard:.3f} "
            f"pf_peak_gain={row.pf_vs_shuffled_peak_mass_gain:.3f} "
            f"v_rank/r={row.valuation_rank_over_period:.3f}"
        )
    write_rows(rows, args.out)
    print(f"[SAVE] {args.out}")
    print_summary(rows)


if __name__ == "__main__":
    main()
