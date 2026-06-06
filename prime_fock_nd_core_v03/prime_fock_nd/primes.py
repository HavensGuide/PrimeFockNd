from __future__ import annotations

import math
from math import isqrt
from typing import List


def first_primes(n: int) -> List[int]:
    """
    Return the first n prime numbers.

    Uses a standard nth-prime upper-bound estimate to avoid repeated
    dynamic reallocation for wide feature spaces, with a fallback growth
    loop for safety.
    """
    if n <= 0:
        raise ValueError("n must be positive.")

    limit = _nth_prime_upper_bound(n)

    while True:
        primes = _sieve(limit)
        if len(primes) >= n:
            return primes[:n]
        # Extremely conservative fallback; should rarely trigger.
        limit = int(limit * 1.5) + 16


def _nth_prime_upper_bound(n: int) -> int:
    """
    Upper bound for the nth prime.

    For n >= 6, p_n < n(log n + log log n) is a classical usable bound.
    Small n are handled explicitly.
    """
    small = [2, 3, 5, 7, 11]
    if n <= len(small):
        return small[n - 1] + 1

    nf = float(n)
    return int(math.ceil(nf * (math.log(nf) + math.log(math.log(nf))))) + 8


def _sieve(limit: int) -> List[int]:
    if limit < 2:
        return []

    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"

    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            step = p
            sieve[start:limit + 1:step] = b"\x00" * (((limit - start) // step) + 1)

    return [i for i, flag in enumerate(sieve) if flag]
