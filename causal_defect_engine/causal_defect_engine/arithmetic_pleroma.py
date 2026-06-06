"""
arithmetic_pleroma.py
=====================

Arithmetic/Pleroma backend for the Causal Defect Algebra Engine.

Selectable modes
----------------
legacy / exact
    Original dense uint8 F2 matrix and dense Gaussian elimination.
packed / efficient
    Bit-packed parity rows for lower memory plus dense fallback dependency
    extraction for small/medium problems.
large_prime
    Relation collector also keeps one-large-prime near-smooth relations and
    merges pairs sharing the same cofactor into full relations.

This module remains a correctness/prototype backend, not a replacement for a
production QS/NFS implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import math
import numpy as np

from .core import CausalDefectEngine, CompositionRule, Defect, DefectState, njit, normalize_mode
from .numba_kernels import trial_factor_candidates_numba, pack_f2_rows_numba, f2_rref_numba


@njit(cache=True)
def _trial_factor_over_base(n: int, base: np.ndarray) -> Tuple[np.ndarray, int]:
    rem = n
    exps = np.zeros(len(base), dtype=np.int64)
    if rem < 0:
        rem = -rem
    for i in range(len(base)):
        p = int(base[i])
        while rem % p == 0:
            exps[i] += 1
            rem //= p
    return exps, rem


@njit(cache=True)
def _mobius_numba(n: int) -> int:
    if n == 1:
        return 1
    m = n
    cnt = 0
    p = 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            cnt += 1
            if m % p == 0:
                return 0
            while m % p == 0:
                m //= p
        p += 1 if p == 2 else 2
    if m > 1:
        cnt += 1
    return -1 if cnt % 2 else 1


@njit(cache=True)
def _von_mangoldt_numba(n: int) -> float:
    if n <= 1:
        return 0.0
    m = n
    found = 0
    found_p = 0
    p = 2
    while p * p <= m:
        if m % p == 0:
            found += 1
            found_p = p
            while m % p == 0:
                m //= p
            if m != 1:
                return 0.0
        p += 1 if p == 2 else 2
    if found == 0:
        return math.log(n)
    return math.log(found_p)


def primes_up_to(B: int) -> np.ndarray:
    if B < 2:
        return np.zeros(0, dtype=np.int64)
    sieve = np.ones(B + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(B ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p : B + 1 : p] = False
    return np.flatnonzero(sieve).astype(np.int64)


def legendre_symbol(a: int, p: int) -> int:
    a %= p
    if a == 0:
        return 0
    val = pow(a, (p - 1) // 2, p)
    return -1 if val == p - 1 else val


class DirichletConvolution(CompositionRule):
    """Composition law for finite-support arithmetic functions."""

    def identity(self) -> Dict[int, complex]:
        return {1: 1.0}

    def compose(self, a: Dict[int, complex], b: Dict[int, complex]) -> Dict[int, complex]:
        out: Dict[int, complex] = {}
        for n, av in a.items():
            for m, bv in b.items():
                out[n * m] = out.get(n * m, 0.0) + av * bv
        return {k: v for k, v in out.items() if abs(v) > 0}

    def invert(self, a: Dict[int, complex]) -> Dict[int, complex]:
        if a.get(1, 0) == 0:
            raise ValueError("Dirichlet inverse requires f(1) != 0")
        max_n = max(a.keys())
        g: Dict[int, complex] = {1: 1.0 / a[1]}
        for n in range(2, max_n + 1):
            s = 0.0
            for d, fd in a.items():
                if d == 1:
                    continue
                if n % d == 0:
                    s += fd * g.get(n // d, 0.0)
            g[n] = -s / a[1]
        return {k: v for k, v in g.items() if abs(v) > 0}


@dataclass
class ArithmeticConfig:
    smooth_bound: int = 1000
    include_negative_one: bool = True
    quadratic_residue_filter: bool = True
    relation_margin: int = 10
    mode: str = "auto"           # legacy | dense | packed | efficient | large_prime | auto
    relation_mode: str = "smooth" # smooth | one_large_prime
    large_prime_factor: int = 20


class ArithmeticPleromaEngine(CausalDefectEngine):
    """Defect algebra backend for integers and factoring-style relation search."""

    def __init__(self, config: Optional[ArithmeticConfig] = None, **kwargs):
        config = config or ArithmeticConfig()
        kwargs.setdefault("mode", config.mode)
        super().__init__(**kwargs)
        self.config = config
        self.N: Optional[int] = None
        self.values: Optional[np.ndarray] = None
        self.factor_base: Optional[np.ndarray] = None
        self.sign_column: bool = self.config.include_negative_one
        self._resolved_mode: str = "dense"

    def fit(self, events: Any, **kwargs) -> "ArithmeticPleromaEngine":
        if isinstance(events, (int, np.integer)):
            self.N = int(events)
            self.values = None
        else:
            vals = np.asarray(list(events), dtype=np.int64)
            if np.any(vals == 0):
                raise ValueError("Arithmetic events must be nonzero integers")
            self.values = vals
            self.N = kwargs.get("N")

        B = int(kwargs.get("smooth_bound", self.config.smooth_bound))
        primes = primes_up_to(B)
        if self.N is not None and self.config.quadratic_residue_filter and self.N > 2:
            primes = np.asarray([p for p in primes if p == 2 or legendre_symbol(self.N, int(p)) >= 0], dtype=np.int64)
        self.factor_base = primes
        self._resolved_mode = self._choose_mode(len(primes))
        self.state = DefectState(
            events=self.values if self.values is not None else self.N,
            fusion_rules=DirichletConvolution(),
            metadata={
                "backend": "arithmetic_pleroma",
                "numba": self.use_numba,
                "mode": self._resolved_mode,
                "N": self.N,
                "smooth_bound": B,
                "factor_base_size": int(len(primes)),
                "sign_column": self.sign_column,
            },
            valid=True,
        )
        return self

    def _choose_mode(self, m: int) -> str:
        mode = normalize_mode(self.config.mode or self.runtime.mode)
        if mode == "auto":
            return "packed" if m >= 128 else "dense"
        if mode == "efficient":
            return "packed"
        if mode == "legacy":
            return "dense"
        return mode

    def estimate_generators(self, **kwargs) -> Dict[str, Any]:
        if self.factor_base is None:
            raise RuntimeError("Call fit first")
        if self.values is None:
            self.state.generators = {"factor_base": self.factor_base}
            return self.state.generators

        exps = []
        rems = []
        signs = []
        saturation = []
        for n in self.values:
            signs.append(1 if n < 0 else 0)
            e, r = _trial_factor_over_base(int(abs(n)), self.factor_base)
            exps.append(e)
            rems.append(int(r))
            resolved = float(np.dot(e, np.log(self.factor_base + 1e-300))) if len(self.factor_base) else 0.0
            saturation.append(resolved / (math.log(abs(int(n))) + 1e-12))
        E = np.vstack(exps) if exps else np.zeros((0, len(self.factor_base)), dtype=np.int64)
        parity = E % 2
        if self.sign_column:
            parity = np.column_stack([np.asarray(signs, dtype=np.int64), parity])
        gens = {
            "factor_base": self.factor_base,
            "exponents": E,
            "parity": parity.astype(np.uint8),
            "parity_packed": self.pack_f2_rows(parity.astype(np.uint8)),
            "remainders": np.asarray(rems, dtype=np.int64),
            "smooth_mask": np.asarray(rems, dtype=np.int64) == 1,
            "defect_saturation": np.asarray(saturation, dtype=np.float64),
            "lambda": np.asarray([_von_mangoldt_numba(int(abs(n))) for n in self.values]),
            "mobius": np.asarray([_mobius_numba(int(abs(n))) for n in self.values]),
        }
        self.state.generators = gens
        self.state.observables.update(
            {
                "smooth_count": int(gens["smooth_mask"].sum()),
                "smooth_rate": float(gens["smooth_mask"].mean()) if len(gens["smooth_mask"]) else 0.0,
                "mean_defect_saturation": float(np.mean(gens["defect_saturation"])) if len(gens["defect_saturation"]) else 0.0,
            }
        )
        return gens

    def build_connection(self, **kwargs) -> Dict[str, Any]:
        if self.state.generators is None:
            self.estimate_generators(**kwargs)
        parity = self.state.generators.get("parity")
        conn = {"boundary_matrix_f2": parity, "boundary_packed": self.state.generators.get("parity_packed")}
        self.state.connection = conn
        return conn

    def detect_defects(self, **kwargs) -> List[Defect]:
        if self.state.generators is None:
            self.estimate_generators(**kwargs)
        defects: List[Defect] = []
        if self.values is not None:
            rems = self.state.generators["remainders"]
            for i, r in enumerate(rems):
                if r != 1:
                    defects.append(Defect(kind="unresolved_prime_defect", charge=1.0, support=(i,), metadata={"cofactor": int(r)}))
        self.state.defects = defects
        return defects

    # ------------------------------------------------------------------
    # F2 representation helpers
    # ------------------------------------------------------------------

    @staticmethod
    def pack_f2_rows(M: np.ndarray) -> np.ndarray:
        M = np.asarray(M, dtype=np.uint8) & 1
        if M.ndim != 2:
            return np.zeros((0, 0), dtype=np.uint64)
        # Prefer central Numba kernel when available; pure Python fallback is kept
        # for environments where Numba import/compilation is disabled.
        try:
            return pack_f2_rows_numba(M)
        except Exception:
            n, m = M.shape
            words = (m + 63) // 64
            out = np.zeros((n, words), dtype=np.uint64)
            for i in range(n):
                for j in range(m):
                    if M[i, j]:
                        out[i, j // 64] |= np.uint64(1) << np.uint64(j % 64)
            return out

    @staticmethod
    def unpack_f2_rows(P: np.ndarray, n_cols: int) -> np.ndarray:
        P = np.asarray(P, dtype=np.uint64)
        out = np.zeros((P.shape[0], n_cols), dtype=np.uint8)
        for i in range(P.shape[0]):
            for j in range(n_cols):
                out[i, j] = np.uint8((P[i, j // 64] >> np.uint64(j % 64)) & np.uint64(1))
        return out

    # ------------------------------------------------------------------
    # Factoring relation helpers
    # ------------------------------------------------------------------

    def collect_quadratic_relations(
        self,
        radius: int,
        center: Optional[int] = None,
        relation_mode: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Collect B-smooth or one-large-prime Q(a)=a^2-N relations."""
        if self.N is None:
            raise ValueError("collect_quadratic_relations requires fit(N)")
        if self.factor_base is None:
            raise RuntimeError("Call fit first")
        mode = normalize_mode(relation_mode or self.config.relation_mode)
        if mode == "large_prime":
            mode = "one_large_prime"
        N = int(self.N)
        root = math.isqrt(N)
        if root * root < N:
            root += 1
        if center is not None:
            root = int(center)

        smooth_rows: List[Dict[str, Any]] = []
        partial_by_cofactor: Dict[int, List[Dict[str, Any]]] = {}
        B_lp = int(self.config.smooth_bound * self.config.large_prime_factor)

        if self.use_numba:
            a_arr, q_arr, exp_mat, rem_arr, sign_arr = trial_factor_candidates_numba(
                int(N), int(root), int(radius), self.factor_base.astype(np.int64)
            )
            logs = np.log(self.factor_base + 1e-300) if len(self.factor_base) else np.zeros(0)
            for idx in range(len(a_arr)):
                a = int(a_arr[idx])
                Q = int(q_arr[idx])
                e = exp_mat[idx].copy()
                rem = int(rem_arr[idx])
                sign = int(sign_arr[idx])
                resolved = float(np.dot(e, logs)) if len(logs) else 0.0
                sat = resolved / (math.log(abs(Q) + 1) + 1e-12)
                if rem == 1:
                    smooth_rows.append(self._relation_row(a, Q, e, sign, rem, sat, "smooth"))
                elif mode == "one_large_prime" and 1 < rem <= B_lp:
                    row = self._relation_row(a, Q, e, sign, int(rem), sat, "partial")
                    partial_by_cofactor.setdefault(int(rem), []).append(row)
        else:
            for t in range(-int(radius), int(radius) + 1):
                a = root + t
                Q = a * a - N
                sign = 1 if Q < 0 else 0
                e, rem = _trial_factor_over_base(abs(Q), self.factor_base)
                resolved = float(np.dot(e, np.log(self.factor_base + 1e-300))) if len(self.factor_base) else 0.0
                sat = resolved / (math.log(abs(Q) + 1) + 1e-12)
                if rem == 1:
                    smooth_rows.append(self._relation_row(a, Q, e, sign, rem, sat, "smooth"))
                elif mode == "one_large_prime" and 1 < rem <= B_lp:
                    row = self._relation_row(a, Q, e, sign, int(rem), sat, "partial")
                    partial_by_cofactor.setdefault(int(rem), []).append(row)

        # Merge pairs sharing the same large cofactor. The large-prime column
        # cancels in parity, and the exponent vector gains +2 in that cofactor,
        # which is irrelevant for the factor-base parity but needed for y only
        # if the large cofactor is included. We store merged rows with a note;
        # y-building ignores the cancelled large prime because it contributes an
        # even exponent to both Q product sides.
        merged_rows: List[Dict[str, Any]] = []
        if mode == "one_large_prime":
            for cof, rows in partial_by_cofactor.items():
                for i in range(0, len(rows) - 1, 2):
                    r1, r2 = rows[i], rows[i + 1]
                    e = r1["exponents"] + r2["exponents"]
                    # Include one exponent of cofactor in metadata for exact y if desired.
                    Qprod_sign = -1 if (r1["Q"] < 0) ^ (r2["Q"] < 0) else 1
                    parity = e % 2
                    sign = 1 if Qprod_sign < 0 else 0
                    if self.sign_column:
                        parity_full = np.concatenate([[sign], parity])
                    else:
                        parity_full = parity
                    merged_rows.append({
                        "a": None,
                        "a_list": [r1["a"], r2["a"]],
                        "Q": r1["Q"] * r2["Q"],
                        "Q_list": [r1["Q"], r2["Q"]],
                        "exponents": e,
                        "parity": parity_full.astype(np.uint8),
                        "cofactor": 1,
                        "large_prime_pair": int(cof),
                        "defect_saturation": min(1.0, 0.5 * (r1["defect_saturation"] + r2["defect_saturation"])),
                        "kind": "merged_large_prime",
                    })

        rows = smooth_rows + merged_rows
        width = len(self.factor_base) + int(self.sign_column)
        matrix = np.vstack([r["parity"] for r in rows]) if rows else np.zeros((0, width), dtype=np.uint8)
        relations = {
            "rows": rows,
            "matrix": matrix,
            "matrix_packed": self.pack_f2_rows(matrix),
            "mode": mode,
            "n_smooth": len(smooth_rows),
            "n_merged_large_prime": len(merged_rows),
            "n_partials_unmerged": sum(len(v) % 2 for v in partial_by_cofactor.values()),
        }
        self.state.observables["n_relations"] = len(rows)
        self.state.observables["n_smooth_relations"] = len(smooth_rows)
        self.state.observables["n_merged_large_prime"] = len(merged_rows)
        return relations

    def _relation_row(self, a: int, Q: int, e: np.ndarray, sign: int, rem: int, sat: float, kind: str) -> Dict[str, Any]:
        parity = e % 2
        if self.sign_column:
            parity = np.concatenate([[sign], parity])
        return {"a": int(a), "Q": int(Q), "exponents": e, "parity": parity.astype(np.uint8), "cofactor": int(rem), "defect_saturation": float(sat), "kind": kind}

    def find_f2_dependencies(self, M: np.ndarray) -> List[np.ndarray]:
        """Return basis vectors c such that c^T M = 0 over F2.

        The implementation intentionally remains dense for correctness and
        portability. Packed rows reduce storage in relation objects; future
        production work should swap this for block Lanczos/Wiedemann.
        """
        A = np.asarray(M, dtype=np.uint8) % 2
        if A.ndim != 2 or A.shape[0] == 0:
            return []
        B = A.T.copy()
        n_rows, n_cols = B.shape
        if self.use_numba:
            aug, piv_arr, rank = f2_rref_numba(B.astype(np.uint8))
            pivots = [int(piv_arr[i]) for i in range(int(rank)) if int(piv_arr[i]) >= 0]
        else:
            aug = B.copy()
            pivots: List[int] = []
            r = 0
            for c in range(n_cols):
                pivot = None
                for i in range(r, n_rows):
                    if aug[i, c]:
                        pivot = i
                        break
                if pivot is None:
                    continue
                if pivot != r:
                    aug[[r, pivot]] = aug[[pivot, r]]
                for i in range(n_rows):
                    if i != r and aug[i, c]:
                        aug[i] ^= aug[r]
                pivots.append(c)
                r += 1
                if r == n_rows:
                    break
        free = [c for c in range(n_cols) if c not in pivots]
        deps = []
        for f in free:
            x = np.zeros(n_cols, dtype=np.uint8)
            x[f] = 1
            for row_idx in range(len(pivots) - 1, -1, -1):
                pc = pivots[row_idx]
                s = 0
                for c in free:
                    if aug[row_idx, c] and x[c]:
                        s ^= 1
                x[pc] = s
            deps.append(x)
        return deps

    def try_factor_from_relations(self, relations: Dict[str, Any]) -> Optional[int]:
        """Try extracting a nontrivial factor from collected relations."""
        if self.N is None:
            raise ValueError("Requires fit(N)")
        rows = relations["rows"]
        M = relations["matrix"]
        deps = self.find_f2_dependencies(M)
        for dep in deps:
            idx = np.flatnonzero(dep)
            if len(idx) == 0:
                continue
            x = 1
            total_exp = np.zeros(len(self.factor_base), dtype=np.int64)
            large_square_factor = 1
            for j in idx:
                row = rows[int(j)]
                if row.get("a") is not None:
                    x = (x * int(row["a"])) % self.N
                else:
                    for aa in row.get("a_list", []):
                        x = (x * int(aa)) % self.N
                total_exp += row["exponents"]
                if row.get("large_prime_pair"):
                    large_square_factor = (large_square_factor * int(row["large_prime_pair"])) % self.N
            y = large_square_factor
            for p, e in zip(self.factor_base, total_exp):
                if e:
                    y = (y * pow(int(p), int(e // 2), self.N)) % self.N
            d = math.gcd((x - y) % self.N, self.N)
            if 1 < d < self.N:
                return d
            d = math.gcd((x + y) % self.N, self.N)
            if 1 < d < self.N:
                return d
        return None
