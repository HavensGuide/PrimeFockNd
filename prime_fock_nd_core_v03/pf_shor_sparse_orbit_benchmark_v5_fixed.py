#!/usr/bin/env python3
"""
Prime-Fock / Shor v5: sparse / partial-orbit interference benchmark.

This benchmark asks the stricter question left open by v4:

    Can Prime-Fock valuation/interference features recover the modular period
    from only m << r sampled orbit points, without constructing the whole orbit?

For each semiprime N and base a, we sample sparse exponent points x in [0, Q),
compute y = a^x mod N, build several feature signals, and use a nonuniform
Fourier / candidate-period spectral score to infer the period.

Signals compared:
  - prime_fock: valuation(a^x mod N)
  - prime_fock_whitened: whitened valuation features
  - shuffled_prime_fock: shuffled orbit valuation control (requires full small orbit for control)
  - onehot_residue: one-hot residue label features
  - random_periodic_features: random vector assigned to orbit phase x mod r

Controls/diagnostics:
  - sparse gcd-defect factor recovery from sampled residues
  - candidate-period score of true r vs best guess
  - sample complexity by m/r

This is a diagnostic benchmark, not a scalable factoring claim.  The current
implementation is intentionally transparent and laptop-friendly for N <= O(1000).
"""
from __future__ import annotations

import argparse
import csv
import math
import random
import time
from dataclasses import dataclass, asdict, fields
from pathlib import Path
from typing import Optional, Iterable

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
    out: list[int] = []
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
    bases: list[int] = []
    for a in range(2, n):
        if gcd(a, n) != 1:
            continue
        r = multiplicative_order(a, n)
        if r is not None and r > 1:
            bases.append(a)
        if len(bases) >= max_bases:
            break
    return bases


def divisors(n: int) -> list[int]:
    out = []
    for d in range(1, int(math.isqrt(n)) + 1):
        if n % d == 0:
            out.append(d)
            if d * d != n:
                out.append(n // d)
    return sorted(out)


def reduce_period_candidate(a: int, n: int, c: int) -> Optional[int]:
    """If c is a valid period multiple, reduce to the minimal divisor period."""
    if c <= 0:
        return None
    if pow(a, c, n) != 1:
        return None
    for d in divisors(c):
        if d > 0 and pow(a, d, n) == 1:
            return d
    return c


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
    out: list[int] = []
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


def matrix_rank(X: np.ndarray, tol: float = 1e-10) -> int:
    if X.size == 0:
        return 0
    try:
        return int(np.linalg.matrix_rank(X, tol=tol))
    except Exception:
        return 0


# ------------------------- sparse feature signals -------------------------

def center_and_scale(X: np.ndarray, *, whiten: bool = False) -> np.ndarray:
    """Center/optionally whiten a sampled feature matrix while preserving row count.

    Important: in sparse sampling, some samples can be all-zero in a valuation
    basis, especially tiny cases like residues equal to 1.  The previous version
    returned zeros((n_features, 1)) when all columns were inactive, which changed
    the number of rows and caused xs/X broadcast errors such as:

        phase[:, None] has shape (m, 1) but X has shape (n_features, 1).

    The row axis must always remain the sampled exponent axis of length m.
    """
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    n_rows = X.shape[0]
    if n_rows == 0:
        return np.zeros((0, 1), dtype=float)
    X = X - X.mean(axis=0, keepdims=True)
    active = np.sum(np.abs(X), axis=0) > 1e-12
    X = X[:, active]
    if X.shape[1] == 0:
        return np.zeros((n_rows, 1), dtype=float)
    if whiten:
        sd = X.std(axis=0, keepdims=True)
        X = X / np.where(sd > 1e-12, sd, 1.0)
    return X


def onehot_features(values: np.ndarray, n: int) -> np.ndarray:
    X = np.zeros((len(values), n), dtype=float)
    for i, y in enumerate(values):
        X[i, int(y) % n] = 1.0
    return X


def random_periodic_features(phases: np.ndarray, r: int, dim: int, rng: np.random.Generator) -> np.ndarray:
    table = rng.normal(size=(r, dim))
    return table[np.asarray(phases, dtype=int) % r]


def sparse_gcd_factor_success(n: int, residues: np.ndarray) -> tuple[bool, int | None, int | None, int]:
    """Try pairwise gcd(y_i - y_j, N). Returns success and pairs checked."""
    m = len(residues)
    checked = 0
    for i in range(m):
        for j in range(i + 1, m):
            checked += 1
            d = gcd(int(residues[i]) - int(residues[j]), n)
            if 1 < d < n:
                return True, d, n // d if n % d == 0 else None, checked
    return False, None, None, checked


# ------------------------- sparse spectral period recovery -------------------------

def coprime_harmonics(c: int, max_harmonics: int, rng: np.random.Generator) -> list[int]:
    js = [j for j in range(1, c) if gcd(j, c) == 1]
    if len(js) <= max_harmonics:
        return js
    # deterministic-ish spread plus random tail for large c
    idx = np.linspace(0, len(js) - 1, max_harmonics, dtype=int)
    return [js[int(i)] for i in idx]


def nonuniform_energy_at_k(xs: np.ndarray, X: np.ndarray, Q: int, k: int) -> float:
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    if len(xs) != X.shape[0]:
        raise ValueError(f"row mismatch: len(xs)={len(xs)} but X.shape={X.shape}")
    phase = np.exp(-2j * np.pi * (float(k) / float(Q)) * xs.astype(float))
    amp = phase[:, None] * X
    z = amp.sum(axis=0)
    return float(np.real(np.vdot(z, z))) / max(1, X.shape[0] * X.shape[1])


def candidate_period_score(
    xs: np.ndarray,
    Xraw: np.ndarray,
    Q: int,
    c: int,
    *,
    max_harmonics: int,
    rng: np.random.Generator,
    whiten: bool = False,
    score_mode: str = "mean",
) -> float:
    X = center_and_scale(Xraw, whiten=whiten)
    if X.shape[0] < 2 or X.shape[1] < 1 or c <= 1:
        return 0.0
    js = coprime_harmonics(c, max_harmonics=max_harmonics, rng=rng)
    if not js:
        return 0.0
    vals: list[float] = []
    seen: set[int] = set()
    for j in js:
        k = int(round(j * Q / c)) % Q
        if k == 0 or k in seen:
            continue
        seen.add(k)
        vals.append(nonuniform_energy_at_k(xs, X, Q, k))
    if not vals:
        return 0.0
    arr = np.asarray(vals, dtype=float)
    if score_mode == "max":
        return float(arr.max())
    if score_mode == "top3mean":
        return float(np.sort(arr)[-min(3, len(arr)):].mean())
    return float(arr.mean())


def recover_period_sparse_spectral(
    xs: np.ndarray,
    Xraw: np.ndarray,
    n: int,
    a: int,
    true_r: int,
    Q: int,
    *,
    max_candidate_period: int,
    max_harmonics: int,
    rng: np.random.Generator,
    whiten: bool = False,
    score_mode: str = "mean",
    top_candidates: int = 8,
) -> dict:
    candidates = list(range(2, max(3, max_candidate_period) + 1))
    scores = []
    for c in candidates:
        s = candidate_period_score(xs, Xraw, Q, c, max_harmonics=max_harmonics, rng=rng, whiten=whiten, score_mode=score_mode)
        scores.append(s)
    scores_np = np.asarray(scores, dtype=float)
    order = np.argsort(-scores_np)
    best_raw = int(candidates[int(order[0])]) if len(order) else -1
    best_score = float(scores_np[int(order[0])]) if len(order) else 0.0
    true_score = float(scores_np[candidates.index(true_r)]) if true_r in candidates else 0.0

    guessed_reduced = None
    guessed_raw = best_raw
    # Prefer top scoring candidate that verifies as a period or period multiple.
    for idx in order[:max(1, top_candidates)]:
        c = int(candidates[int(idx)])
        red = reduce_period_candidate(a, n, c)
        if red is not None:
            guessed_raw = c
            guessed_reduced = red
            best_score = float(scores_np[int(idx)])
            break
    success = int(guessed_reduced == true_r)
    rank_true = int(np.where(np.asarray(candidates)[order] == true_r)[0][0] + 1) if true_r in candidates else -1
    # Gap positive means true period scored near/above chosen best raw candidate.
    score_gap_true_vs_best = true_score - best_score
    return {
        "guessed_raw_period": guessed_raw,
        "guessed_reduced_period": guessed_reduced if guessed_reduced is not None else -1,
        "period_success": success,
        "true_period_score": true_score,
        "best_period_score": best_score,
        "score_gap_true_minus_best": score_gap_true_vs_best,
        "true_period_score_rank": rank_true,
    }


# ------------------------- rows / cases -------------------------

@dataclass
class SparseOrbitRow:
    N: int
    a: int
    period_r: int
    Q: int
    sample_m: int
    sample_m_over_r: float
    trial: int
    signal_name: str
    feature_dim: int
    feature_rank: int
    observed_distinct_residues: int
    observed_phase_coverage: float
    guessed_raw_period: int
    guessed_reduced_period: int
    period_success: int
    true_period_score: float
    best_period_score: float
    score_gap_true_minus_best: float
    true_period_score_rank: int
    sparse_gcd_factor_success: int
    sparse_gcd_pairs_checked: int
    valuation_rank_over_period: float
    elapsed_sec: float


def make_cases(args) -> list[tuple[int, int, int]]:
    cases: list[tuple[int, int, int]] = []
    if args.Ns:
        Ns = [int(x.strip()) for x in args.Ns.split(",") if x.strip()]
    else:
        Ns = semiprimes(args.max_N, min_n=args.min_N)
    for n in Ns:
        bases = valid_bases(n, max_bases=args.max_bases)
        for a in bases:
            r = multiplicative_order(a, n)
            if r is None or r <= 1:
                continue
            cases.append((n, a, r))
            if args.max_cases and len(cases) >= args.max_cases:
                return cases
    return cases


def parse_int_list(s: str) -> list[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def run_case(args, n: int, a: int, r: int, case_idx: int) -> list[SparseOrbitRow]:
    rows: list[SparseOrbitRow] = []
    Q = choose_Q(n, args.q_mode, args.q_bits, args.max_q)
    n_primes = args.n_primes or max(8, min(128, int(math.ceil(math.sqrt(max(n, 2)))) + 8))
    full_residues, V_full, primes = valuation_sequence(n, a, r, n_primes=n_primes)
    v_rank = matrix_rank(V_full)
    valuation_rank_over_period = float(v_rank / max(1, r))

    sample_sizes = parse_int_list(args.sample_sizes)
    base_rng = np.random.default_rng(args.seed + 1000003 * case_idx + 7919 * a + n)

    for m0 in sample_sizes:
        m = int(max(2, min(m0, Q)))
        for trial in range(args.trials):
            t0 = time.time()
            rng = np.random.default_rng(int(base_rng.integers(0, 2**31 - 1)) + trial)
            replace = m > Q
            xs = np.sort(rng.choice(Q, size=m, replace=replace).astype(np.int64))
            phases = xs % r
            residues = np.array([pow(a, int(x), n) for x in xs], dtype=np.int64)
            observed_distinct = int(len(set(int(y) for y in residues)))
            observed_phase_coverage = float(len(set(int(p) for p in phases)) / max(1, r))
            gcd_ok, _, _, gcd_pairs = sparse_gcd_factor_success(n, residues)

            # Build signal feature matrices.
            pf_X = np.stack([valuation_vector(int(y), primes) for y in residues], axis=0)

            # Shuffled control: same valuation vectors assigned to orbit phase after a random permutation.
            perm = np.arange(r)
            rng.shuffle(perm)
            shuffled_full = V_full[perm]
            shuffled_X = shuffled_full[phases]

            onehot_X = onehot_features(residues, n)
            random_X = random_periodic_features(phases, r, dim=min(args.random_dim, max(2, pf_X.shape[1])), rng=rng)

            signal_map = {
                "prime_fock": (pf_X, False),
                "prime_fock_whitened": (pf_X, True),
                "shuffled_prime_fock": (shuffled_X, False),
                "onehot_residue": (onehot_X, False),
                "random_periodic_features": (random_X, False),
            }
            if args.signals != "all":
                keep = set(x.strip() for x in args.signals.split(",") if x.strip())
                signal_map = {k: v for k, v in signal_map.items() if k in keep}

            for name, (X, whiten) in signal_map.items():
                rec = recover_period_sparse_spectral(
                    xs,
                    X,
                    n,
                    a,
                    r,
                    Q,
                    max_candidate_period=min(args.max_candidate_period, n + 1),
                    max_harmonics=args.max_harmonics,
                    rng=rng,
                    whiten=whiten,
                    score_mode=args.score_mode,
                    top_candidates=args.top_candidates,
                )
                rows.append(SparseOrbitRow(
                    N=n,
                    a=a,
                    period_r=r,
                    Q=Q,
                    sample_m=m,
                    sample_m_over_r=float(m / max(1, r)),
                    trial=trial,
                    signal_name=name,
                    feature_dim=int(X.shape[1] if X.ndim == 2 else 1),
                    feature_rank=matrix_rank(center_and_scale(X, whiten=whiten)),
                    observed_distinct_residues=observed_distinct,
                    observed_phase_coverage=observed_phase_coverage,
                    guessed_raw_period=int(rec["guessed_raw_period"]),
                    guessed_reduced_period=int(rec["guessed_reduced_period"]),
                    period_success=int(rec["period_success"]),
                    true_period_score=float(rec["true_period_score"]),
                    best_period_score=float(rec["best_period_score"]),
                    score_gap_true_minus_best=float(rec["score_gap_true_minus_best"]),
                    true_period_score_rank=int(rec["true_period_score_rank"]),
                    sparse_gcd_factor_success=int(gcd_ok),
                    sparse_gcd_pairs_checked=int(gcd_pairs),
                    valuation_rank_over_period=valuation_rank_over_period,
                    elapsed_sec=float(time.time() - t0),
                ))
    return rows


def write_csv(path: str | Path, rows: list[SparseOrbitRow]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = [f.name for f in fields(SparseOrbitRow)]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for row in rows:
            w.writerow(asdict(row))


def summarize(rows: list[SparseOrbitRow]) -> None:
    if not rows:
        print("[SUMMARY] no rows")
        return
    print("\n[SUMMARY]")
    print(f"rows: {len(rows)}")
    # Basic per-signal means.
    sigs = sorted(set(r.signal_name for r in rows))
    for sig in sigs:
        rr = [r for r in rows if r.signal_name == sig]
        succ = np.mean([r.period_success for r in rr])
        rank = np.mean([r.true_period_score_rank for r in rr if r.true_period_score_rank > 0])
        gap = np.mean([r.score_gap_true_minus_best for r in rr])
        print(f"{sig:26s} success={succ:.3f} mean_true_rank={rank:.2f} mean_gap={gap:.4g} rows={len(rr)}")
    # Per sample size for the main signals.
    print("\n[BY SAMPLE SIZE]")
    for m in sorted(set(r.sample_m for r in rows)):
        chunk = [r for r in rows if r.sample_m == m]
        bits = []
        for sig in sigs:
            rr = [r for r in chunk if r.signal_name == sig]
            if rr:
                bits.append(f"{sig}={np.mean([x.period_success for x in rr]):.2f}")
        print(f"m={m:4d}: " + " | ".join(bits))
    # Compare Prime-Fock to controls.
    pf = [r.period_success for r in rows if r.signal_name == "prime_fock"]
    sh = [r.period_success for r in rows if r.signal_name == "shuffled_prime_fock"]
    oh = [r.period_success for r in rows if r.signal_name == "onehot_residue"]
    gd = [r.sparse_gcd_factor_success for r in rows if r.signal_name == "prime_fock"]
    if pf:
        print("\n[KEY]")
        print(f"Prime-Fock sparse period success: {np.mean(pf):.3f}")
        if sh: print(f"Shuffled PF sparse period success: {np.mean(sh):.3f}")
        if oh: print(f"One-hot residue sparse period success: {np.mean(oh):.3f}")
        if gd: print(f"Sparse gcd factor success: {np.mean(gd):.3f}")


def main() -> None:
    ap = argparse.ArgumentParser(description="Prime-Fock Shor sparse/partial orbit interference benchmark v5")
    ap.add_argument("--Ns", type=str, default="", help="Comma-separated N values. If empty, generate semiprimes.")
    ap.add_argument("--min-N", type=int, default=15)
    ap.add_argument("--max-N", type=int, default=300)
    ap.add_argument("--max-bases", type=int, default=4)
    ap.add_argument("--max-cases", type=int, default=80)
    ap.add_argument("--q-mode", type=str, default="n2", choices=["n", "n2", "2n2"])
    ap.add_argument("--q-bits", type=int, default=None)
    ap.add_argument("--max-q", type=int, default=65536)
    ap.add_argument("--n-primes", type=int, default=None)
    ap.add_argument("--sample-sizes", type=str, default="4,8,16,32,64")
    ap.add_argument("--trials", type=int, default=20)
    ap.add_argument("--max-candidate-period", type=int, default=512)
    ap.add_argument("--max-harmonics", type=int, default=24)
    ap.add_argument("--score-mode", type=str, default="top3mean", choices=["mean", "max", "top3mean"])
    ap.add_argument("--top-candidates", type=int, default=12)
    ap.add_argument("--signals", type=str, default="all", help="all or comma-list of signal names")
    ap.add_argument("--random-dim", type=int, default=16)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=str, default="shor_pf_sparse_orbit_v5.csv")
    args = ap.parse_args()

    cases = make_cases(args)
    print(f"[CASES] {len(cases)} cases")
    all_rows: list[SparseOrbitRow] = []
    for i, (n, a, r) in enumerate(cases, 1):
        t0 = time.time()
        rows = run_case(args, n, a, r, i)
        all_rows.extend(rows)
        # compact status for prime_fock only
        pf_rows = [x for x in rows if x.signal_name == "prime_fock"]
        pf_succ = np.mean([x.period_success for x in pf_rows]) if pf_rows else float("nan")
        gcd_succ = np.mean([x.sparse_gcd_factor_success for x in pf_rows]) if pf_rows else float("nan")
        print(f"[{i:03d}/{len(cases)}] N={n} a={a} r={r} rows={len(rows)} pf_success={pf_succ:.2f} gcd_factor={gcd_succ:.2f} elapsed={time.time()-t0:.2f}s")

    write_csv(args.out, all_rows)
    print(f"[SAVE] {args.out}")
    summarize(all_rows)


if __name__ == "__main__":
    main()
