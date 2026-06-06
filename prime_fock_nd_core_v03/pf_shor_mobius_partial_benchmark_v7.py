#!/usr/bin/env python3
"""
Prime-Fock / Shor v7: partial-prefix Möbius / divisor-lattice closure benchmark.

v6 showed that full-orbit circular lag closure is extremely cheap, but too easy:
all sufficiently injective feature maps recover the period r once the full orbit is known.

v7 asks a stricter, more algorithmically relevant question:

  If we observe only a prefix of the modular orbit y_x = a^x mod N, how much
  orbit exposure is needed before cheap divisor/Prime-Fock closure scoring can
  recover the true period r?

Important detail:
  Pure prefix lag closure cannot identify a lag c unless the prefix length M
  exceeds c by at least a small overlap. Therefore v7 measures exposure as

      exposure_ratio = M / r
      true_overlap   = M - r
      overlap_ratio  = (M - r) / r

  A ratio below 1.0 cannot score the true period at all. This benchmark defaults
  to ratios just above 1.0, e.g. 1.05, 1.10, 1.25, 1.50, 2.00.

For each case and prefix length M:
  1. Compute y_0, ..., y_{M-1}.
  2. Build feature maps:
       - onehot_residue
       - prime_valuation
       - squarefree_prime_support
       - zeta_divisor_truncK
       - mobius_delta_truncK
       - mobius_squarefree_truncK
  3. For candidate lag c, score only observed pairs:

       score(c) = mean_{0 <= x < M-c} <F(y_x), F(y_{x+c})>

     with a minimum required overlap.
  4. Compare which feature map recovers r at the lowest M/r, with fewest false
     aliases c < r and the best true-minus-false margin.

This is still classical and still uses exact modular powers; it is a sample-
complexity / representation benchmark, not a claim of Shor-like speedup.
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
    if n < 2:
        return []
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(math.isqrt(n)) + 1):
        if sieve[p]:
            sieve[p*p:n+1:p] = False
    return [i for i in range(2, n + 1) if sieve[i]]


def mobius_sieve(n: int) -> np.ndarray:
    mu = np.ones(n + 1, dtype=np.int8)
    is_p = np.ones(n + 1, dtype=bool)
    if n >= 0:
        mu[0] = 0
    if n >= 1:
        is_p[:2] = False
    primes: list[int] = []
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
            mu[i * p] = -mu[i]
    return mu


def orbit_prefix_residues(a: int, n: int, m: int) -> np.ndarray:
    ys = np.empty(m, dtype=np.int64)
    y = 1
    for x in range(m):
        if x == 0:
            y = 1
        elif x > 0:
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
        if y <= 1:
            continue
        for j, p in enumerate(primes):
            while y > 0 and y % p == 0:
                X[i, j] += 1.0
                y //= p
            if y == 1:
                break
    return X


def support_features(ys: np.ndarray, primes: list[int]) -> np.ndarray:
    return (valuation_features(ys, primes) > 0).astype(np.float64)


def zeta_divisor_features(ys: np.ndarray, divisors: np.ndarray) -> np.ndarray:
    X = np.zeros((len(ys), len(divisors)), dtype=np.float64)
    for i, y0 in enumerate(ys):
        y = int(y0)
        if y == 0:
            continue
        X[i, :] = ((y % divisors) == 0).astype(np.float64)
    return X


def mobius_invert_zeta_rows(Z: np.ndarray, divisors: np.ndarray, mu: np.ndarray) -> np.ndarray:
    D = len(divisors)
    M = np.zeros((D, D), dtype=np.float64)
    div_list = [int(d) for d in divisors]
    for i, d in enumerate(div_list):
        for j, e in enumerate(div_list):
            if e % d == 0:
                q = e // d
                if q < len(mu):
                    M[j, i] = float(mu[q])
    return Z @ M


def mobius_squarefree_features(ys: np.ndarray, squarefree_ds: np.ndarray, mu: np.ndarray) -> np.ndarray:
    X = zeta_divisor_features(ys, squarefree_ds)
    signs = np.array([mu[int(d)] for d in squarefree_ds], dtype=np.float64)
    return X * signs[None, :]


def build_feature_blocks(ys: np.ndarray, n: int, args) -> dict[str, np.ndarray]:
    primes = primes_upto(n)
    mu = mobius_sieve(n)
    all_ds = np.arange(1, n + 1, dtype=np.int64)
    squarefree_ds = np.array([d for d in range(2, n + 1) if mu[d] != 0], dtype=np.int64)
    small_divisors = all_ds[: min(len(all_ds), args.zeta_truncate)]
    small_squarefree = squarefree_ds[: min(len(squarefree_ds), args.mobius_truncate)]

    blocks: dict[str, np.ndarray] = {
        "onehot_residue": onehot_residue_features(ys, n),
        "prime_valuation": valuation_features(ys, primes),
        "squarefree_prime_support": support_features(ys, primes),
    }
    Z = zeta_divisor_features(ys, small_divisors)
    blocks[f"zeta_divisor_trunc{len(small_divisors)}"] = Z
    blocks[f"mobius_delta_trunc{len(small_divisors)}"] = mobius_invert_zeta_rows(Z, small_divisors, mu)
    blocks[f"mobius_squarefree_trunc{len(small_squarefree)}"] = mobius_squarefree_features(ys, small_squarefree, mu)
    if args.include_full_squarefree:
        blocks["mobius_squarefree_full"] = mobius_squarefree_features(ys, squarefree_ds, mu)
    if args.include_full_divisor:
        Zf = zeta_divisor_features(ys, all_ds)
        blocks["zeta_divisor_full"] = Zf
        blocks["mobius_delta_full"] = mobius_invert_zeta_rows(Zf, all_ds, mu)
    return blocks


# ------------------------- scoring -------------------------

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
    if len(s) == 0:
        return 0, 0.0, 0.0
    tol = max(Xc.shape) * np.finfo(float).eps * max(float(s[0]), 1.0)
    rank = int(np.sum(s > tol))
    stable = float((np.sum(s * s) ** 2) / max(np.sum(s ** 4), 1e-30))
    p = (s * s) / max(float(np.sum(s * s)), 1e-30)
    ent = float(-np.sum(p[p > 0] * np.log2(p[p > 0])))
    return rank, stable, ent


def prefix_lag_scores(X: np.ndarray, max_candidate_lag: int, min_overlap: int) -> tuple[np.ndarray, np.ndarray]:
    """Non-circular prefix lag scores for candidate lags c.

    score(c) = mean_{x < M-c} dot(X_x, X_{x+c}) if M-c >= min_overlap, else NaN.
    """
    Xn = center_scale(X, normalize_rows=True)
    M = Xn.shape[0]
    max_candidate_lag = min(max_candidate_lag, M - min_overlap)
    scores = np.full(max_candidate_lag + 1, np.nan, dtype=np.float64)
    overlaps = np.zeros(max_candidate_lag + 1, dtype=np.int64)
    for c in range(1, max_candidate_lag + 1):
        ov = M - c
        if ov < min_overlap:
            continue
        overlaps[c] = ov
        scores[c] = float(np.mean(np.sum(Xn[:ov] * Xn[c:c+ov], axis=1)))
    return scores, overlaps


def recover_prefix(scores: np.ndarray, true_r: int, tol: float) -> tuple[int, float, int, float, float, int]:
    if len(scores) <= 1 or np.all(~np.isfinite(scores[1:])):
        return -1, float("nan"), -1, float("nan"), float("nan"), 0
    sc = scores.copy()
    sc[0] = np.nan
    finite = np.isfinite(sc)
    candidate_count = int(np.sum(finite))
    best = int(np.nanargmax(sc))
    best_score = float(sc[best])
    strong = np.where(finite & (sc >= tol))[0]
    first_strong = int(strong[0]) if len(strong) else best
    true_score = float(sc[true_r]) if true_r < len(sc) and np.isfinite(sc[true_r]) else float("nan")
    if true_r > 1 and np.any(np.isfinite(sc[1:true_r])):
        false_alias = float(np.nanmax(sc[1:true_r]))
    else:
        false_alias = float("-inf")
    return best, best_score, first_strong, true_score, false_alias, candidate_count


@dataclass
class Row:
    case_id: int
    N: int
    a: int
    period_r: int
    exposure_ratio_target: float
    observed_len: int
    observed_over_period: float
    true_overlap: int
    true_overlap_over_period: float
    min_overlap: int
    feature: str
    feature_dim: int
    feature_rank: int
    feature_stable_rank: float
    feature_spectral_entropy: float
    rank_over_period: float
    dim_over_period: float
    candidate_count: int
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


def exposure_lengths(r: int, ratios: list[float], min_overlap: int) -> list[tuple[float, int]]:
    out = []
    seen = set()
    for ratio in ratios:
        m = int(math.ceil(r * ratio))
        # To score true r at all, need at least min_overlap points after r.
        if ratio >= 1.0:
            m = max(m, r + min_overlap)
        if m not in seen:
            out.append((ratio, m))
            seen.add(m)
    return out


def run_case(args, n: int, a: int, r: int, case_id: int) -> list[Row]:
    ratios = [float(x) for x in args.exposure_ratios.split(',') if x.strip()]
    rows: list[Row] = []
    for ratio, M in exposure_lengths(r, ratios, args.min_overlap):
        ys = orbit_prefix_residues(a, n, M)
        tb = time.time()
        blocks = build_feature_blocks(ys, n, args)
        build_total = time.time() - tb
        max_c = min(args.max_candidate_period, M - args.min_overlap)
        for name, X in blocks.items():
            ts = time.time()
            rank, stable, ent = matrix_rank_stats(center_scale(X, normalize_rows=False))
            scores, overlaps = prefix_lag_scores(X, max_c, args.min_overlap)
            best, best_score, first, true_score, false_alias, cand_count = recover_prefix(scores, r, args.closure_tol)
            score_seconds = time.time() - ts
            true_overlap = int(overlaps[r]) if r < len(overlaps) else 0
            rows.append(Row(
                case_id=case_id,
                N=n,
                a=a,
                period_r=r,
                exposure_ratio_target=ratio,
                observed_len=M,
                observed_over_period=float(M / max(r, 1)),
                true_overlap=true_overlap,
                true_overlap_over_period=float(true_overlap / max(r, 1)),
                min_overlap=args.min_overlap,
                feature=name,
                feature_dim=int(X.shape[1]),
                feature_rank=rank,
                feature_stable_rank=stable,
                feature_spectral_entropy=ent,
                rank_over_period=float(rank / max(r, 1)),
                dim_over_period=float(X.shape[1] / max(r, 1)),
                candidate_count=cand_count,
                best_lag=best,
                best_lag_score=best_score,
                first_strong_lag=first,
                true_period_score=true_score,
                false_alias_max_score=false_alias,
                recovered_exact_best=int(best == r),
                recovered_exact_first_strong=int(first == r),
                false_alias_before_r=int(false_alias >= args.closure_tol),
                alias_margin_true_minus_false=float(true_score - false_alias) if np.isfinite(true_score) else float("nan"),
                build_seconds=build_total / max(len(blocks), 1),
                score_seconds=score_seconds,
            ))
    return rows


def make_cases(args) -> list[tuple[int, int, int]]:
    if args.Ns:
        Ns = [int(x.strip()) for x in args.Ns.split(',') if x.strip()]
    else:
        Ns = semiprimes(args.max_N, args.min_N)
    cases: list[tuple[int, int, int]] = []
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
    print("\n[SUMMARY BY FEATURE x EXPOSURE]")
    print("feature                         obs/r     n  exact_best first_exact alias  rank/r  margin score_s")
    groups: dict[tuple[str, float], list[Row]] = {}
    for row in rows:
        groups.setdefault((row.feature, round(row.observed_over_period, 3)), []).append(row)
    # show common aggregate features first; still print all groups.
    for (feat, obsr), rr in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        n = len(rr)
        exact = np.mean([x.recovered_exact_best for x in rr])
        first = np.mean([x.recovered_exact_first_strong for x in rr])
        alias = np.mean([x.false_alias_before_r for x in rr])
        rankr = np.mean([x.rank_over_period for x in rr])
        margin = np.nanmean([x.alias_margin_true_minus_false for x in rr])
        sec = np.mean([x.score_seconds for x in rr])
        print(f"{feat:31s} {obsr:6.3f} {n:5d} {exact:10.3f} {first:11.3f} {alias:5.3f} {rankr:7.3f} {margin:7.3f} {sec:7.4f}")

    print("\n[MIN EXPOSURE FOR >= TARGET SUCCESS]")
    target = 0.90
    by_feat: dict[str, list[Row]] = {}
    for row in rows:
        by_feat.setdefault(row.feature, []).append(row)
    for feat, rr in sorted(by_feat.items()):
        ratios = sorted(set(round(x.observed_over_period, 3) for x in rr))
        hit = None
        for obsr in ratios:
            sub = [x for x in rr if round(x.observed_over_period, 3) == obsr]
            if sub and np.mean([x.recovered_exact_best for x in sub]) >= target:
                hit = obsr
                break
        print(f"{feat:31s} min_obs/r@{target:.2f}={hit}")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Partial-prefix Möbius / divisor-lattice Shor closure benchmark v7")
    ap.add_argument("--Ns", type=str, default="", help="Comma-separated N values. Overrides --max-N.")
    ap.add_argument("--min-N", type=int, default=15)
    ap.add_argument("--max-N", type=int, default=300)
    ap.add_argument("--max-bases", type=int, default=3)
    ap.add_argument("--max-cases", type=int, default=80)
    ap.add_argument("--min-period", type=int, default=0)
    ap.add_argument("--max-period", type=int, default=0)
    ap.add_argument("--exposure-ratios", type=str, default="1.05,1.10,1.25,1.50,2.00", help="Comma-separated observed_len / period targets. Prefix closure needs >1.")
    ap.add_argument("--min-overlap", type=int, default=4, help="Minimum pair overlap M-c required to score candidate lag c.")
    ap.add_argument("--max-candidate-period", type=int, default=4096)
    ap.add_argument("--closure-tol", type=float, default=0.995)
    ap.add_argument("--zeta-truncate", type=int, default=128)
    ap.add_argument("--mobius-truncate", type=int, default=128)
    ap.add_argument("--include-full-divisor", action="store_true")
    ap.add_argument("--include-full-squarefree", action="store_true")
    ap.add_argument("--out", type=str, default="shor_pf_mobius_partial_v7.csv")
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
        # Compact case summary for prime and mobius-squarefree at smallest/largest exposure.
        pv = [x for x in rows if x.feature == "prime_valuation"]
        ms = [x for x in rows if x.feature.startswith("mobius_squarefree_trunc")]
        def curve(xs: list[Row]) -> str:
            return ",".join(f"{x.observed_over_period:.2f}:{x.recovered_exact_best}" for x in sorted(xs, key=lambda z: z.observed_over_period)[:5])
        print(f"[{i:03d}/{len(cases)}] N={n} a={a} r={r} prime={curve(pv)} mobius_sq={curve(ms)} elapsed={time.time()-t0:.3f}s")
    write_csv(args.out, all_rows)
    print(f"[SAVE] {args.out}")
    print_summary(all_rows)
    print(f"[TOTAL SECONDS] {time.time()-t_all:.2f}")


if __name__ == "__main__":
    main()
