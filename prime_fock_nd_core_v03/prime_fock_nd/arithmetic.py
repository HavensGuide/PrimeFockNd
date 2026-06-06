from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np


ArithmeticKernelName = Literal[
    "log_gcd",
    "log_lcm",
    "arithmetic_distance",
    "jaccard_gcd_lcm",
    "overlap",
    "coprime_indicator",
]


@dataclass(frozen=True)
class ArithmeticKernelResult:
    matrix: np.ndarray
    kernel: str


class ValuationArithmetic:
    """
    Arithmetic operations over prime-valuation vectors.

    If an item has exponent vector v, it represents the implicit integer:

        n(v) = prod_k p_k ** v_k

    Literal construction of n(v) is intentionally avoided because it may be
    enormous. All operations happen in valuation space.

    Core identities:
        multiplication: v(nm) = v(n) + v(m)
        gcd:            v(gcd(n,m)) = min(v(n), v(m))
        lcm:            v(lcm(n,m)) = max(v(n), v(m))
        quotient:       v(n/m) = v(n) - v(m), when divisible
    """

    @staticmethod
    def multiply(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        return np.asarray(a, dtype=float) + np.asarray(b, dtype=float)

    @staticmethod
    def gcd(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        return np.minimum(np.asarray(a, dtype=float), np.asarray(b, dtype=float))

    @staticmethod
    def lcm(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        return np.maximum(np.asarray(a, dtype=float), np.asarray(b, dtype=float))

    @staticmethod
    def quotient(a: np.ndarray, b: np.ndarray, *, clip_negative: bool = False) -> np.ndarray:
        out = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
        if clip_negative:
            out = np.maximum(out, 0.0)
        return out

    @staticmethod
    def divides(a: np.ndarray, b: np.ndarray) -> bool:
        """
        Return True if implicit integer n(a) divides n(b).

        This holds when every valuation of a is <= corresponding valuation of b.
        """
        return bool(np.all(np.asarray(a, dtype=float) <= np.asarray(b, dtype=float)))

    @staticmethod
    def coprime(a: np.ndarray, b: np.ndarray, *, atol: float = 1e-12) -> bool:
        return bool(np.all(np.minimum(np.asarray(a, dtype=float), np.asarray(b, dtype=float)) <= atol))

    @staticmethod
    def log_integer(a: np.ndarray, log_primes: np.ndarray) -> np.ndarray:
        return np.asarray(a, dtype=float) @ np.asarray(log_primes, dtype=float)

    @staticmethod
    def log_gcd(a: np.ndarray, b: np.ndarray, log_primes: np.ndarray) -> np.ndarray:
        return np.minimum(np.asarray(a, dtype=float), np.asarray(b, dtype=float)) @ np.asarray(log_primes, dtype=float)

    @staticmethod
    def log_lcm(a: np.ndarray, b: np.ndarray, log_primes: np.ndarray) -> np.ndarray:
        return np.maximum(np.asarray(a, dtype=float), np.asarray(b, dtype=float)) @ np.asarray(log_primes, dtype=float)

    @staticmethod
    def arithmetic_distance(a: np.ndarray, b: np.ndarray, log_primes: np.ndarray) -> np.ndarray:
        """
        log(lcm/gcd) = sum_k |a_k - b_k| log p_k.
        """
        return np.abs(np.asarray(a, dtype=float) - np.asarray(b, dtype=float)) @ np.asarray(log_primes, dtype=float)

    @staticmethod
    def jaccard_gcd_lcm(a: np.ndarray, b: np.ndarray, log_primes: np.ndarray, eps: float = 1e-12) -> np.ndarray:
        lg = ValuationArithmetic.log_gcd(a, b, log_primes)
        ll = ValuationArithmetic.log_lcm(a, b, log_primes)
        return lg / (ll + eps)


def pairwise_valuation_kernel(
    valuations: np.ndarray,
    log_primes: np.ndarray,
    *,
    kernel: ArithmeticKernelName = "log_gcd",
    eps: float = 1e-12,
) -> ArithmeticKernelResult:
    """
    Compute a pairwise arithmetic kernel matrix from valuation vectors.

    Parameters
    ----------
    valuations:
        Shape (I, P), nonnegative exponent vectors.
    log_primes:
        Shape (P,), log prime weights defining the implicit integer scale.
    kernel:
        Arithmetic kernel to compute.

    Returns
    -------
    ArithmeticKernelResult with matrix shape (I, I).
    """
    v = np.asarray(valuations, dtype=float)
    logs = np.asarray(log_primes, dtype=float)

    if v.ndim != 2:
        raise ValueError("valuations must have shape (I, P).")
    if logs.shape != (v.shape[1],):
        raise ValueError(f"log_primes must have shape {(v.shape[1],)}, got {logs.shape}.")
    if not np.all(np.isfinite(v)) or not np.all(np.isfinite(logs)):
        raise ValueError("valuations and log_primes must be finite.")
    if np.any(v < 0):
        raise ValueError("valuations must be nonnegative.")

    if kernel in {"log_gcd", "overlap"}:
        shared = np.minimum(v[:, None, :], v[None, :, :])
        mat = np.sum(shared * logs.reshape(1, 1, -1), axis=2)
    elif kernel == "log_lcm":
        union = np.maximum(v[:, None, :], v[None, :, :])
        mat = np.sum(union * logs.reshape(1, 1, -1), axis=2)
    elif kernel == "arithmetic_distance":
        diff = np.abs(v[:, None, :] - v[None, :, :])
        mat = np.sum(diff * logs.reshape(1, 1, -1), axis=2)
    elif kernel == "jaccard_gcd_lcm":
        shared = np.minimum(v[:, None, :], v[None, :, :])
        union = np.maximum(v[:, None, :], v[None, :, :])
        lg = np.sum(shared * logs.reshape(1, 1, -1), axis=2)
        ll = np.sum(union * logs.reshape(1, 1, -1), axis=2)
        mat = lg / (ll + eps)
    elif kernel == "coprime_indicator":
        shared_mass = np.sum(np.minimum(v[:, None, :], v[None, :, :]), axis=2)
        mat = (shared_mass <= eps).astype(float)
    else:
        raise ValueError(f"Unknown arithmetic kernel: {kernel!r}")

    np.fill_diagonal(mat, 0.0)
    return ArithmeticKernelResult(matrix=mat, kernel=kernel)
