#!/usr/bin/env python3
"""
Prime-Fock / Shor v8: predictive closure before seeing the repeat.

v6/v7 showed that divisor/Prime-Fock closure scoring is extremely cheap, but once
we observe y_x and y_{x+r} pairs, the task is too easy. v8 removes those pairs.

Question:
  Given only a prefix y_0, ..., y_{L-1} with L < r, can a cheap transition model
  in a feature space F(y) infer the modular period r by rolling forward and
  detecting predicted closure?

For each case y_x = a^x mod N:
  1. Choose L = ceil(train_ratio * r), with train_ratio < 1 by default.
  2. Build feature maps:
       - onehot_residue
       - prime_valuation
       - squarefree_prime_support
       - zeta_divisor_truncK
       - mobius_delta_truncK
       - mobius_squarefree_truncK
  3. Fit ridge linear transition on the prefix:
       z_{t+1} ~= A z_t + b,  t=0,...,L-2
  4. Roll z_0 forward without using future residues.
  5. Score candidate periods h by cosine(z_hat_h, z_0).
  6. Compare predicted best h to the true multiplicative order r.

This is a classical representation/sample-complexity diagnostic. It is not a
claim of Shor-like speedup. A meaningful win would be exact period recovery with
L/r < 1, especially for compressed arithmetic features beating one-hot residue.
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


def orbit_residues(a: int, n: int, m: int) -> np.ndarray:
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


# ------------------------- preprocessing / operator -------------------------

class Scaler:
    def __init__(self, mean: np.ndarray, std: np.ndarray, keep: np.ndarray, normalize_rows: bool = True):
        self.mean = mean
        self.std = std
        self.keep = keep
        self.normalize_rows = normalize_rows

    @classmethod
    def fit(cls, X: np.ndarray, normalize_rows: bool = True) -> "Scaler":
        X = np.asarray(X, dtype=np.float64)
        if X.ndim == 1:
            X = X[:, None]
        if X.shape[1] == 0:
            X = np.zeros((X.shape[0], 1), dtype=np.float64)
        mean = X.mean(axis=0, keepdims=True)
        std = X.std(axis=0, keepdims=True)
        keep = (std.ravel() > 1e-12)
        if not keep.any():
            keep = np.zeros(X.shape[1], dtype=bool)
        return cls(mean, std, keep, normalize_rows=normalize_rows)

    def transform(self, X: np.ndarray) -> np.ndarray:
        X = np.asarray(X, dtype=np.float64)
        if X.ndim == 1:
            X = X[:, None]
        if X.shape[1] == 0:
            X = np.zeros((X.shape[0], 1), dtype=np.float64)
        if self.keep.any():
            Z = (X[:, self.keep] - self.mean[:, self.keep]) / self.std[:, self.keep]
        else:
            Z = np.zeros((X.shape[0], 1), dtype=np.float64)
        if self.normalize_rows:
            nr = np.linalg.norm(Z, axis=1, keepdims=True)
            Z = Z / np.maximum(nr, 1e-12)
        return Z


def lift_features(Z: np.ndarray, mode: str) -> np.ndarray:
    if mode == "linear":
        return Z
    if mode == "linear_square":
        return np.concatenate([Z, Z * Z], axis=1)
    if mode == "linear_pair_lite":
        # Add pair interactions only up to a capped dimension to keep laptop-safe.
        d = Z.shape[1]
        cap = min(d, 24)
        pairs = []
        for i in range(cap):
            for j in range(i, cap):
                pairs.append((Z[:, i] * Z[:, j])[:, None])
        if pairs:
            return np.concatenate([Z] + pairs, axis=1)
        return Z
    raise ValueError(f"unknown lift mode: {mode}")


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


def fit_ridge_transition(Z: np.ndarray, alpha: float, use_bias: bool = True) -> tuple[np.ndarray, float, float, float]:
    """Fit Z[t] -> Z[t+1]. Returns W, train_r2, train_rmse, spectral_radius."""
    if Z.shape[0] < 2:
        d = Z.shape[1]
        return np.eye(d), float("nan"), float("nan"), 1.0
    X = Z[:-1]
    Y = Z[1:]
    if use_bias:
        Xd = np.concatenate([X, np.ones((X.shape[0], 1))], axis=1)
    else:
        Xd = X
    G = Xd.T @ Xd
    G.flat[::G.shape[0]+1] += alpha
    try:
        W = np.linalg.solve(G, Xd.T @ Y)
    except np.linalg.LinAlgError:
        W = np.linalg.pinv(G) @ Xd.T @ Y
    Yhat = Xd @ W
    resid = Y - Yhat
    rmse = float(np.sqrt(np.mean(resid * resid)))
    sse = float(np.sum(resid * resid))
    sst = float(np.sum((Y - Y.mean(axis=0, keepdims=True)) ** 2))
    r2 = float("nan") if sst < 1e-18 else float(1.0 - sse / sst)
    # spectral radius of linear part only
    A = W[:-1, :] if use_bias else W
    try:
        eig = np.linalg.eigvals(A)
        sr = float(np.max(np.abs(eig))) if len(eig) else 0.0
    except Exception:
        sr = float("nan")
    return W, r2, rmse, sr


def one_step_oos_r2(Z: np.ndarray, W: np.ndarray, L: int, horizon_end: int, use_bias: bool = True) -> tuple[float, float]:
    if L >= horizon_end - 1:
        return float("nan"), float("nan")
    X = Z[L:horizon_end-1]
    Y = Z[L+1:horizon_end]
    if X.shape[0] == 0:
        return float("nan"), float("nan")
    if use_bias:
        Xd = np.concatenate([X, np.ones((X.shape[0], 1))], axis=1)
    else:
        Xd = X
    Yhat = Xd @ W
    resid = Y - Yhat
    rmse = float(np.sqrt(np.mean(resid * resid)))
    sse = float(np.sum(resid * resid))
    sst = float(np.sum((Y - Y.mean(axis=0, keepdims=True)) ** 2))
    return (float("nan") if sst < 1e-18 else float(1.0 - sse / sst)), rmse


def rollout(W: np.ndarray, z0: np.ndarray, steps: int, use_bias: bool = True, renorm: bool = True, clip_norm: float = 1e6) -> tuple[np.ndarray, float, float]:
    d = z0.shape[0]
    Zhat = np.zeros((steps + 1, d), dtype=np.float64)
    z = z0.astype(np.float64).copy()
    n0 = np.linalg.norm(z)
    if renorm:
        z = z / max(n0, 1e-12)
    Zhat[0] = z
    max_norm = float(np.linalg.norm(z))
    diverged = 0.0
    for t in range(1, steps + 1):
        if use_bias:
            xd = np.concatenate([z, np.ones(1)])
        else:
            xd = z
        z = xd @ W
        nz = float(np.linalg.norm(z))
        if not np.isfinite(nz) or nz > clip_norm:
            diverged = 1.0
            z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
            nz = float(np.linalg.norm(z))
            if nz > 0:
                z = z / nz
        if renorm:
            z = z / max(float(np.linalg.norm(z)), 1e-12)
        Zhat[t] = z
        max_norm = max(max_norm, float(np.linalg.norm(z)))
    return Zhat, max_norm, diverged


def closure_recovery(Zhat: np.ndarray, true_r: int, min_candidate: int, max_candidate: int) -> dict[str, float | int]:
    z0 = Zhat[0]
    z0n = z0 / max(float(np.linalg.norm(z0)), 1e-12)
    max_candidate = min(max_candidate, Zhat.shape[0] - 1)
    scores = np.full(max_candidate + 1, np.nan, dtype=np.float64)
    for h in range(max(1, min_candidate), max_candidate + 1):
        zh = Zhat[h]
        scores[h] = float(np.dot(zh, z0n) / max(float(np.linalg.norm(zh)), 1e-12))
    if np.all(~np.isfinite(scores[1:])):
        return dict(best_lag=-1, best_score=float("nan"), true_score=float("nan"), false_alias=float("nan"), margin=float("nan"), exact=0, within_1=0, within_5pct=0)
    best = int(np.nanargmax(scores))
    best_score = float(scores[best])
    true_score = float(scores[true_r]) if true_r < len(scores) and np.isfinite(scores[true_r]) else float("nan")
    lo = max(1, min_candidate)
    if true_r > lo and np.any(np.isfinite(scores[lo:true_r])):
        false_alias = float(np.nanmax(scores[lo:true_r]))
    else:
        false_alias = float("-inf")
    return dict(
        best_lag=best,
        best_score=best_score,
        true_score=true_score,
        false_alias=false_alias,
        margin=float(true_score - false_alias) if np.isfinite(true_score) else float("nan"),
        exact=int(best == true_r),
        within_1=int(abs(best - true_r) <= 1),
        within_5pct=int(abs(best - true_r) <= max(1, int(round(0.05 * true_r)))),
    )


@dataclass
class Row:
    case_id: int
    N: int
    a: int
    period_r: int
    train_ratio_target: float
    train_len: int
    train_over_period: float
    roll_steps: int
    feature: str
    lift_mode: str
    raw_feature_dim: int
    model_dim: int
    feature_rank_prefix: int
    feature_stable_rank_prefix: float
    feature_spectral_entropy_prefix: float
    rank_over_period: float
    dim_over_period: float
    ridge_alpha: float
    train_one_step_r2: float
    train_one_step_rmse: float
    oos_one_step_r2: float
    oos_one_step_rmse: float
    spectral_radius: float
    rollout_max_norm: float
    rollout_diverged: float
    best_lag: int
    best_lag_score: float
    true_period_score: float
    false_alias_max_score: float
    alias_margin_true_minus_false: float
    recovered_exact_best: int
    recovered_within_1: int
    recovered_within_5pct: int
    build_seconds: float
    fit_roll_seconds: float


def train_lengths(r: int, ratios: list[float], min_train: int) -> list[tuple[float, int]]:
    out = []
    seen = set()
    for ratio in ratios:
        L = max(min_train, int(math.ceil(r * ratio)))
        # enforce prefix before closure by default; cap to r-1
        L = min(L, r - 1)
        if L >= 2 and L not in seen:
            out.append((ratio, L))
            seen.add(L)
    return out


def run_case(args, n: int, a: int, r: int, case_id: int) -> list[Row]:
    ratios = [float(x) for x in args.train_ratios.split(',') if x.strip()]
    # Need full feature trajectory to evaluate OOS one-step and compare true period closure.
    roll_steps = min(args.max_roll_steps, max(r * args.roll_periods, r + 1))
    full_len = roll_steps + 1
    ys_full = orbit_residues(a, n, full_len)
    tb = time.time()
    blocks = build_feature_blocks(ys_full, n, args)
    build_total = time.time() - tb
    rows: list[Row] = []
    for ratio, L in train_lengths(r, ratios, args.min_train_len):
        for name, Xraw_full in blocks.items():
            t0 = time.time()
            X_prefix = Xraw_full[:L]
            scaler = Scaler.fit(X_prefix, normalize_rows=args.normalize_rows)
            Xs_full = scaler.transform(Xraw_full)
            Z_full = lift_features(Xs_full, args.lift_mode)
            # Re-normalize lifted states for stable cosine dynamics.
            if args.normalize_rows:
                nr = np.linalg.norm(Z_full, axis=1, keepdims=True)
                Z_full = Z_full / np.maximum(nr, 1e-12)
            Z_train = Z_full[:L]
            rank, stable, ent = matrix_rank_stats(Z_train)
            W, tr_r2, tr_rmse, sr = fit_ridge_transition(Z_train, args.ridge_alpha, use_bias=not args.no_bias)
            oos_r2, oos_rmse = one_step_oos_r2(Z_full, W, L, min(r, full_len), use_bias=not args.no_bias)
            Zhat, max_norm, diverged = rollout(W, Z_train[0], roll_steps, use_bias=not args.no_bias, renorm=args.renorm_rollout, clip_norm=args.clip_norm)
            rec = closure_recovery(Zhat, r, args.min_candidate_period, min(args.max_candidate_period, roll_steps))
            dt = time.time() - t0
            rows.append(Row(
                case_id=case_id,
                N=n,
                a=a,
                period_r=r,
                train_ratio_target=ratio,
                train_len=L,
                train_over_period=float(L / max(r, 1)),
                roll_steps=roll_steps,
                feature=name,
                lift_mode=args.lift_mode,
                raw_feature_dim=int(Xraw_full.shape[1]),
                model_dim=int(Z_train.shape[1]),
                feature_rank_prefix=rank,
                feature_stable_rank_prefix=stable,
                feature_spectral_entropy_prefix=ent,
                rank_over_period=float(rank / max(r, 1)),
                dim_over_period=float(Z_train.shape[1] / max(r, 1)),
                ridge_alpha=float(args.ridge_alpha),
                train_one_step_r2=tr_r2,
                train_one_step_rmse=tr_rmse,
                oos_one_step_r2=oos_r2,
                oos_one_step_rmse=oos_rmse,
                spectral_radius=sr,
                rollout_max_norm=max_norm,
                rollout_diverged=diverged,
                best_lag=int(rec["best_lag"]),
                best_lag_score=float(rec["best_score"]),
                true_period_score=float(rec["true_score"]),
                false_alias_max_score=float(rec["false_alias"]),
                alias_margin_true_minus_false=float(rec["margin"]),
                recovered_exact_best=int(rec["exact"]),
                recovered_within_1=int(rec["within_1"]),
                recovered_within_5pct=int(rec["within_5pct"]),
                build_seconds=build_total / max(len(blocks), 1),
                fit_roll_seconds=dt,
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
    print("\n[SUMMARY BY FEATURE x TRAIN]")
    print("feature                         L/r      n exact within5 oos_r2   rank/r  trueS margin fit_s")
    groups: dict[tuple[str, float], list[Row]] = {}
    for row in rows:
        groups.setdefault((row.feature, round(row.train_over_period, 3)), []).append(row)
    for (feat, lr), rr in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        n = len(rr)
        exact = np.mean([x.recovered_exact_best for x in rr])
        within5 = np.mean([x.recovered_within_5pct for x in rr])
        vals_oos = [x.oos_one_step_r2 for x in rr if np.isfinite(x.oos_one_step_r2)]
        oos = float(np.mean(vals_oos)) if vals_oos else float("nan")
        rankr = np.mean([x.rank_over_period for x in rr])
        trueS = np.nanmean([x.true_period_score for x in rr])
        margin = np.nanmean([x.alias_margin_true_minus_false for x in rr])
        sec = np.mean([x.fit_roll_seconds for x in rr])
        print(f"{feat:31s} {lr:6.3f} {n:5d} {exact:5.3f} {within5:7.3f} {oos:7.3f} {rankr:8.3f} {trueS:6.3f} {margin:6.3f} {sec:6.4f}")

    print("\n[MIN TRAIN EXPOSURE FOR >= TARGET EXACT]")
    target = 0.50
    by_feat: dict[str, list[Row]] = {}
    for row in rows:
        by_feat.setdefault(row.feature, []).append(row)
    for feat, rr in sorted(by_feat.items()):
        ratios = sorted(set(round(x.train_over_period, 3) for x in rr))
        hit = None
        for lr in ratios:
            sub = [x for x in rr if round(x.train_over_period, 3) == lr]
            if sub and np.mean([x.recovered_exact_best for x in sub]) >= target:
                hit = lr
                break
        best_exact = max(np.mean([x.recovered_exact_best for x in rr if round(x.train_over_period,3)==lr]) for lr in ratios)
        best_w5 = max(np.mean([x.recovered_within_5pct for x in rr if round(x.train_over_period,3)==lr]) for lr in ratios)
        print(f"{feat:31s} min_L/r@{target:.2f}={hit} best_exact={best_exact:.3f} best_within5={best_w5:.3f}")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Predictive Prime-Fock/Mobius Shor closure benchmark v8")
    ap.add_argument("--Ns", type=str, default="", help="Comma-separated N values. Overrides --max-N.")
    ap.add_argument("--min-N", type=int, default=15)
    ap.add_argument("--max-N", type=int, default=300)
    ap.add_argument("--max-bases", type=int, default=3)
    ap.add_argument("--max-cases", type=int, default=80)
    ap.add_argument("--min-period", type=int, default=4)
    ap.add_argument("--max-period", type=int, default=0)
    ap.add_argument("--train-ratios", type=str, default="0.20,0.33,0.50,0.67,0.80,0.90", help="Comma-separated L/r training prefix ratios; capped below 1.")
    ap.add_argument("--min-train-len", type=int, default=4)
    ap.add_argument("--roll-periods", type=int, default=2, help="Roll at least this many true periods for evaluation.")
    ap.add_argument("--max-roll-steps", type=int, default=4096)
    ap.add_argument("--min-candidate-period", type=int, default=2)
    ap.add_argument("--max-candidate-period", type=int, default=4096)
    ap.add_argument("--ridge-alpha", type=float, default=1e-6)
    ap.add_argument("--no-bias", action="store_true")
    ap.add_argument("--normalize-rows", action="store_true", default=True)
    ap.add_argument("--no-normalize-rows", dest="normalize_rows", action="store_false")
    ap.add_argument("--renorm-rollout", action="store_true", default=True)
    ap.add_argument("--no-renorm-rollout", dest="renorm_rollout", action="store_false")
    ap.add_argument("--clip-norm", type=float, default=1e6)
    ap.add_argument("--lift-mode", type=str, default="linear", choices=["linear", "linear_square", "linear_pair_lite"])
    ap.add_argument("--zeta-truncate", type=int, default=128)
    ap.add_argument("--mobius-truncate", type=int, default=128)
    ap.add_argument("--include-full-divisor", action="store_true")
    ap.add_argument("--include-full-squarefree", action="store_true")
    ap.add_argument("--out", type=str, default="shor_pf_predictive_closure_v8.csv")
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
        # compact per-case trace: show prime and mobius_squarefree exact successes by train ratio
        def trace_for(feat: str) -> str:
            sub = [x for x in rows if x.feature.startswith(feat)]
            toks = []
            for x in sub:
                if x.recovered_exact_best:
                    toks.append(f"{x.train_over_period:.2f}:1")
            return ",".join(toks[:6]) if toks else "none"
        print(f"[{i:03d}/{len(cases)}] N={n} a={a} r={r} prime={trace_for('prime_valuation')} mobius_sq={trace_for('mobius_squarefree')} elapsed={time.time()-t0:.3f}s")
    write_csv(args.out, all_rows)
    print(f"[SAVE] {args.out}")
    print_summary(all_rows)
    print(f"[TOTAL SECONDS] {time.time()-t_all:.2f}")


if __name__ == "__main__":
    main()
