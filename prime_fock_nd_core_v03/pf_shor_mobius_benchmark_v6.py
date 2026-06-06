#!/usr/bin/env python3
"""
Prime-Fock / Shor v6: Möbius / divisor-lattice period benchmark.

This benchmark tests the "native arithmetic transform" direction directly.
It does NOT try to replace the additive Fourier/QFT over the exponent register.
Instead it asks:

  Can divisor-zeta / Möbius features of the residue orbit y_x = a^x mod N
  recover the modular period using cheap lag/closure tests, and do they compress
  the orbit better than one-hot residues or simple prime valuations?

For each semiprime N and base a:
  1. Compute one exact modular orbit of length r.
  2. Build several feature representations of each residue y_x:
       - onehot_residue           exact residue label control
       - prime_valuation          prime exponent vector of y_x
       - squarefree_prime_support binary prime support vector
       - zeta_divisor             zeta transform: 1[d | y_x]
       - mobius_residue_delta     Möbius inversion of zeta, equals one-hot over divisor atoms
       - mobius_squarefree        μ(d) 1[d | y_x] for squarefree d
       - mobius_truncated_k       same, restricted to small squarefree divisors
  3. Score candidate lags c by circular self-similarity X_x · X_{x+c}.
     This is a no-FFT/no-QFT period probe. A representation wins if its first
     strong closure lag is r, not a false alias < r.
  4. Track rank/compression, dimension, alias count, and runtime.

Interpretation:
  If Möbius/zeta features recover r with dimension/rank << r and fewer false
  aliases than valuation features, the divisor-lattice transform is useful as an
  arithmetic front-end. If only one-hot/zeta-full wins, it is exact but not a
  compression win. If all compressed divisor features alias, the additive-cycle
  period is not captured by divisor structure alone.
"""
from __future__ import annotations

import argparse
import csv
import math
import time
from dataclasses import dataclass, asdict, fields
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


def valid_bases(n: int, max_bases: int) -> list[int]:
    out = []
    for a in range(2, n):
        if gcd(a, n) != 1:
            continue
        r = multiplicative_order(a, n)
        if r is not None and r > 1:
            out.append(a)
        if len(out) >= max_bases:
            break
    return out


def primes_upto(n: int) -> list[int]:
    sieve = np.ones(n + 1, dtype=bool)
    if n >= 0:
        sieve[:2] = False
    for p in range(2, int(math.isqrt(n)) + 1):
        if sieve[p]:
            sieve[p*p:n+1:p] = False
    return [i for i in range(2, n + 1) if sieve[i]]


def mobius_sieve(n: int) -> np.ndarray:
    """Return μ(0..n)."""
    mu = np.ones(n + 1, dtype=np.int8)
    is_p = np.ones(n + 1, dtype=bool)
    is_p[:2] = False
    primes = []
    mu[0] = 0
    for i in range(2, n + 1):
        if is_p[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            if i * p > n:
                break
            is_p[i * p] = False
            if i % p == 0:
                mu[i * p] = 0
                break
            else:
                mu[i * p] = -mu[i]
    return mu


def orbit_residues(a: int, n: int, r: int) -> np.ndarray:
    ys = np.empty(r, dtype=np.int64)
    y = 1
    ys[0] = 1
    for x in range(1, r):
        y = (y * a) % n
        ys[x] = y
    return ys


# ------------------------- features -------------------------

def center_scale(X: np.ndarray, normalize_rows: bool = True) -> np.ndarray:
    X = np.asarray(X, dtype=np.float64)
    if X.ndim == 1:
        X = X[:, None]
    if X.shape[1] == 0:
        X = np.zeros((X.shape[0], 1), dtype=np.float64)
    X = X - X.mean(axis=0, keepdims=True)
    std = X.std(axis=0, keepdims=True)
    keep = (std.ravel() > 1e-12)
    if keep.any():
        X = X[:, keep] / std[:, keep]
    else:
        X = np.zeros((X.shape[0], 1), dtype=np.float64)
    if normalize_rows:
        nr = np.linalg.norm(X, axis=1, keepdims=True)
        X = X / np.maximum(nr, 1e-12)
    return X


def onehot_residue_features(ys: np.ndarray, n: int) -> np.ndarray:
    X = np.zeros((len(ys), n), dtype=np.float64)
    X[np.arange(len(ys)), ys] = 1.0
    return X


def valuation_features(ys: np.ndarray, primes: list[int]) -> np.ndarray:
    X = np.zeros((len(ys), len(primes)), dtype=np.float64)
    for i, y0 in enumerate(ys):
        y = int(y0)
        for j, p in enumerate(primes):
            while y > 0 and y % p == 0:
                X[i, j] += 1.0
                y //= p
            if y == 1:
                break
    return X


def support_features(ys: np.ndarray, primes: list[int]) -> np.ndarray:
    X = valuation_features(ys, primes)
    return (X > 0).astype(np.float64)


def zeta_divisor_features(ys: np.ndarray, divisors: np.ndarray) -> np.ndarray:
    # Z_y(d)=1[d|y].  This is the divisor-zeta transform of residue delta_y.
    X = np.zeros((len(ys), len(divisors)), dtype=np.float64)
    for i, y in enumerate(ys):
        y = int(y)
        if y == 0:
            continue
        X[i, :] = ((y % divisors) == 0).astype(np.float64)
    return X


def mobius_invert_zeta_rows(Z: np.ndarray, divisors: np.ndarray, mu: np.ndarray) -> np.ndarray:
    """Invert Z_y(d)=sum_{m: d|m} delta_y(m) on the divisor poset.

    F_y(d) = sum_{e: d|e} μ(e/d) Z_y(e).  If all divisors 1..N are present,
    this gives a delta on d=y.  With a restricted divisor set, it is a projected
    Möbius atomization.
    """
    D = len(divisors)
    M = np.zeros((D, D), dtype=np.float64)
    div_list = divisors.tolist()
    idx = {int(d): j for j, d in enumerate(div_list)}
    for i, d in enumerate(div_list):
        for j, e in enumerate(div_list):
            if e % d == 0:
                q = e // d
                if q < len(mu):
                    M[j, i] = float(mu[q])  # row Z_e contributes to F_d
    return Z @ M


def mobius_squarefree_features(ys: np.ndarray, squarefree_ds: np.ndarray, mu: np.ndarray) -> np.ndarray:
    # χ_d(y) = μ(d) 1[d|y], squarefree divisor characters.
    X = zeta_divisor_features(ys, squarefree_ds)
    signs = np.array([mu[int(d)] for d in squarefree_ds], dtype=np.float64)
    return X * signs[None, :]


# ------------------------- period scoring -------------------------

def matrix_rank_stats(X: np.ndarray) -> tuple[int, float, float]:
    Xc = np.asarray(X, dtype=np.float64)
    if Xc.ndim == 1:
        Xc = Xc[:, None]
    if min(Xc.shape) == 0:
        return 0, 0.0, 0.0
    try:
        s = np.linalg.svd(Xc, compute_uv=False)
    except np.linalg.LinAlgError:
        return 0, 0.0, 0.0
    tol = max(Xc.shape) * np.finfo(float).eps * (s[0] if len(s) else 1.0)
    rank = int(np.sum(s > tol))
    stable = float((np.sum(s * s) ** 2) / max(np.sum(s ** 4), 1e-30)) if len(s) else 0.0
    p = (s * s) / max(float(np.sum(s * s)), 1e-30)
    ent = float(-np.sum(p[p > 0] * np.log2(p[p > 0]))) if len(p) else 0.0
    return rank, stable, ent


def lag_similarity_scores(X: np.ndarray, max_lag: int) -> np.ndarray:
    """Circular mean cosine/dot self-similarity at lags 1..max_lag."""
    Xn = center_scale(X, normalize_rows=True)
    r = Xn.shape[0]
    max_lag = min(max_lag, max(1, r))
    scores = np.zeros(max_lag + 1, dtype=np.float64)
    for c in range(1, max_lag + 1):
        scores[c] = float(np.mean(np.sum(Xn * np.roll(Xn, -c, axis=0), axis=1)))
    return scores


def recover_from_lag_scores(scores: np.ndarray, true_r: int, tol: float = 0.995) -> tuple[int, float, int, float, float]:
    """Return best lag, score, first strong closure, true score, false alias max."""
    if len(scores) <= 1:
        return -1, float('nan'), -1, float('nan'), float('nan')
    sc = scores.copy()
    sc[0] = -np.inf
    best = int(np.nanargmax(sc))
    best_score = float(sc[best])
    true_score = float(sc[true_r]) if true_r < len(sc) else float('nan')
    strong = np.where(sc >= tol)[0]
    first_strong = int(strong[0]) if len(strong) else best
    false_alias_max = float(np.nanmax(sc[1:true_r])) if true_r > 1 else float('-inf')
    return best, best_score, first_strong, true_score, false_alias_max


@dataclass
class Row:
    case_id: int
    N: int
    a: int
    period_r: int
    feature: str
    feature_dim: int
    feature_rank: int
    feature_stable_rank: float
    feature_spectral_entropy: float
    rank_over_period: float
    dim_over_period: float
    best_lag: int
    best_lag_score: float
    first_strong_lag: int
    true_period_score: float
    false_alias_max_score: float
    recovered_exact_best: int
    recovered_exact_first_strong: int
    false_alias_before_r: int
    alias_margin_true_minus_false: float
    build_seconds: float
    score_seconds: float


def build_feature_blocks(ys: np.ndarray, n: int, args) -> dict[str, np.ndarray]:
    primes = primes_upto(n)
    mu = mobius_sieve(n)
    all_ds = np.arange(1, n + 1, dtype=np.int64)
    squarefree_ds = np.array([d for d in range(2, n + 1) if mu[d] != 0], dtype=np.int64)
    small_squarefree = squarefree_ds[: min(len(squarefree_ds), args.mobius_truncate)]
    small_divisors = all_ds[: min(len(all_ds), args.zeta_truncate)]

    blocks: dict[str, np.ndarray] = {}
    blocks["onehot_residue"] = onehot_residue_features(ys, n)
    blocks["prime_valuation"] = valuation_features(ys, primes)
    blocks["squarefree_prime_support"] = support_features(ys, primes)

    Z_small = zeta_divisor_features(ys, small_divisors)
    blocks[f"zeta_divisor_trunc{len(small_divisors)}"] = Z_small
    blocks[f"mobius_delta_trunc{len(small_divisors)}"] = mobius_invert_zeta_rows(Z_small, small_divisors, mu)

    # Full divisor/zeta can be high-dimensional but for N<=1000 still cheap.
    if args.include_full_divisor:
        Z_full = zeta_divisor_features(ys, all_ds)
        blocks["zeta_divisor_full"] = Z_full
        blocks["mobius_delta_full"] = mobius_invert_zeta_rows(Z_full, all_ds, mu)

    blocks[f"mobius_squarefree_trunc{len(small_squarefree)}"] = mobius_squarefree_features(ys, small_squarefree, mu)
    if args.include_full_squarefree:
        blocks["mobius_squarefree_full"] = mobius_squarefree_features(ys, squarefree_ds, mu)
    return blocks


def run_case(args, n: int, a: int, r: int, case_id: int) -> list[Row]:
    ys = orbit_residues(a, n, r)
    rows: list[Row] = []
    t0 = time.time()
    blocks = build_feature_blocks(ys, n, args)
    build_total = time.time() - t0
    for name, X in blocks.items():
        tb = time.time()
        rank, stable, ent = matrix_rank_stats(center_scale(X, normalize_rows=False))
        scores = lag_similarity_scores(X, max_lag=min(args.max_lag, r))
        best, best_score, first, true_score, false_alias = recover_from_lag_scores(scores, r, tol=args.closure_tol)
        score_seconds = time.time() - tb
        rows.append(Row(
            case_id=case_id,
            N=n,
            a=a,
            period_r=r,
            feature=name,
            feature_dim=int(X.shape[1]),
            feature_rank=rank,
            feature_stable_rank=stable,
            feature_spectral_entropy=ent,
            rank_over_period=float(rank / max(r, 1)),
            dim_over_period=float(X.shape[1] / max(r, 1)),
            best_lag=best,
            best_lag_score=best_score,
            first_strong_lag=first,
            true_period_score=true_score,
            false_alias_max_score=false_alias,
            recovered_exact_best=int(best == r),
            recovered_exact_first_strong=int(first == r),
            false_alias_before_r=int(false_alias >= args.closure_tol),
            alias_margin_true_minus_false=float(true_score - false_alias),
            build_seconds=build_total / max(len(blocks), 1),
            score_seconds=score_seconds,
        ))
    return rows


def make_cases(args) -> list[tuple[int, int, int]]:
    cases: list[tuple[int, int, int]] = []
    if args.Ns:
        Ns = [int(x.strip()) for x in args.Ns.split(',') if x.strip()]
    else:
        Ns = semiprimes(args.max_N, args.min_N)
    for n in Ns:
        for a in valid_bases(n, args.max_bases):
            r = multiplicative_order(a, n)
            if r is None or r <= 1:
                continue
            if args.min_period and r < args.min_period:
                continue
            if args.max_period and r > args.max_period:
                continue
            cases.append((n, a, r))
            if args.max_cases and len(cases) >= args.max_cases:
                return cases
    return cases


def write_csv(path: str | Path, rows: list[Row]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = [f.name for f in fields(Row)]
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for row in rows:
            w.writerow(asdict(row))


def print_summary(rows: list[Row]) -> None:
    if not rows:
        print("[SUMMARY] no rows")
        return
    by: dict[str, list[Row]] = {}
    for row in rows:
        by.setdefault(row.feature, []).append(row)
    print("\n[SUMMARY BY FEATURE]")
    print("feature                              n  exact_best exact_first alias  rank/r   dim/r  margin  score_s")
    for feat, rr in sorted(by.items()):
        n = len(rr)
        exact = np.mean([x.recovered_exact_best for x in rr])
        first = np.mean([x.recovered_exact_first_strong for x in rr])
        alias = np.mean([x.false_alias_before_r for x in rr])
        rankr = np.mean([x.rank_over_period for x in rr])
        dimr = np.mean([x.dim_over_period for x in rr])
        margin = np.mean([x.alias_margin_true_minus_false for x in rr])
        sec = np.mean([x.score_seconds for x in rr])
        print(f"{feat:34s} {n:3d} {exact:10.3f} {first:11.3f} {alias:5.3f} {rankr:7.3f} {dimr:7.2f} {margin:7.3f} {sec:7.4f}")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Möbius / divisor-lattice Shor period benchmark")
    ap.add_argument("--Ns", type=str, default="", help="Comma-separated N values. Overrides --max-N.")
    ap.add_argument("--min-N", type=int, default=15)
    ap.add_argument("--max-N", type=int, default=300)
    ap.add_argument("--max-bases", type=int, default=3)
    ap.add_argument("--max-cases", type=int, default=80)
    ap.add_argument("--min-period", type=int, default=0)
    ap.add_argument("--max-period", type=int, default=0)
    ap.add_argument("--max-lag", type=int, default=4096)
    ap.add_argument("--closure-tol", type=float, default=0.995)
    ap.add_argument("--zeta-truncate", type=int, default=128, help="Use divisors 1..K for truncated zeta/Mobius delta.")
    ap.add_argument("--mobius-truncate", type=int, default=128, help="Use first K squarefree d for truncated Mobius squarefree features.")
    ap.add_argument("--include-full-divisor", action="store_true", help="Also include full zeta/Mobius over all d<=N.")
    ap.add_argument("--include-full-squarefree", action="store_true")
    ap.add_argument("--out", type=str, default="shor_pf_mobius_v6.csv")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    cases = make_cases(args)
    print(f"[CASES] {len(cases)} cases")
    all_rows: list[Row] = []
    t_all = time.time()
    for i, (n, a, r) in enumerate(cases, 1):
        t0 = time.time()
        rows = run_case(args, n, a, r, i)
        all_rows.extend(rows)
        # compact row for key features
        d = {row.feature: row for row in rows}
        pv = d.get("prime_valuation")
        ms = next((row for row in rows if row.feature.startswith("mobius_squarefree_trunc")), None)
        zd = next((row for row in rows if row.feature.startswith("zeta_divisor_trunc")), None)
        print(
            f"[{i:03d}/{len(cases)}] N={n} a={a} r={r} "
            f"prime_exact={pv.recovered_exact_best if pv else -1} "
            f"mobius_exact={ms.recovered_exact_best if ms else -1} "
            f"zeta_exact={zd.recovered_exact_best if zd else -1} "
            f"rank/r={ms.rank_over_period if ms else float('nan'):.3f} "
            f"elapsed={time.time()-t0:.3f}s"
        )
    write_csv(args.out, all_rows)
    print(f"[SAVE] {args.out}")
    print_summary(all_rows)
    print(f"[TOTAL SECONDS] {time.time()-t_all:.2f}")


if __name__ == "__main__":
    main()
