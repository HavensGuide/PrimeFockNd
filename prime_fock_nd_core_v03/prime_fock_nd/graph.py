from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

from .arithmetic import ArithmeticKernelName, pairwise_valuation_kernel
from .encoder import EncodedND


@dataclass(frozen=True)
class SharedSupportGraph:
    """
    Dense arithmetic graph for each slice.

    adjacency:
        Shape (S, I, I), symmetric with zero diagonal.

    By default this is a log-gcd graph over implicit prime-product states:

        E_ij = log gcd(n_i, n_j)
             = sum_p min(v_ip, v_jp) log p

    where n_i = prod_p p ** v_ip is never materialized.
    """
    adjacency: np.ndarray
    active_modes: np.ndarray
    slices: list[tuple]
    items: list
    kernel: str = "log_gcd"


@dataclass(frozen=True)
class EdgeListGraph:
    """
    Sparse arithmetic graph represented as an edge list.

    edges columns:
        slice_idx, source_idx, target_idx, weight

    Metadata maps indices back to user-facing slice/item labels.
    """
    edges: pd.DataFrame
    n_slices: int
    n_items: int
    slices: list[tuple]
    items: list
    kernel: str = "log_gcd"

    @property
    def empty(self) -> bool:
        return self.edges.empty


class ArithmeticGraphBuilder:
    """
    Dense arithmetic kernel graph builder.

    Supported kernels:
        - log_gcd: log gcd of implicit prime-product states
        - log_lcm: log lcm
        - arithmetic_distance: log(lcm/gcd)
        - jaccard_gcd_lcm: log-gcd / log-lcm
        - coprime_indicator: 1 if no shared prime support
        - overlap: backward-compatible alias for log_gcd

    For large item counts, use SparseArithmeticGraphBuilder for log_gcd/overlap.
    """

    def __init__(
        self,
        mode_weights: Optional[np.ndarray] = None,
        *,
        kernel: ArithmeticKernelName = "log_gcd",
    ):
        self.mode_weights = mode_weights
        self.kernel = kernel

    def build(self, encoded: EncodedND) -> SharedSupportGraph:
        x = np.asarray(encoded.exponents, dtype=float)
        if x.ndim != 3:
            raise ValueError("encoded.exponents must have shape (S, I, P).")
        if not np.all(np.isfinite(x)):
            raise ValueError("encoded.exponents contains non-finite values.")
        if np.any(x < 0):
            raise ValueError("Fock/valuation exponents must be nonnegative.")

        log_primes = _resolve_mode_weights(encoded, self.mode_weights)

        S, I, P = x.shape
        adjacency = np.zeros((S, I, I), dtype=float)

        for s in range(S):
            adjacency[s] = pairwise_valuation_kernel(
                x[s],
                log_primes,
                kernel=self.kernel,
            ).matrix

        return SharedSupportGraph(
            adjacency=adjacency,
            active_modes=(x > 0).astype(float),
            slices=encoded.slices,
            items=encoded.items,
            kernel=self.kernel,
        )


class SharedSupportGraphBuilder(ArithmeticGraphBuilder):
    """
    Backward-compatible alias for ArithmeticGraphBuilder.

    Defaults to kernel='log_gcd', which is the arithmetic interpretation of
    the original weighted overlap graph.
    """
    pass


class SparseArithmeticGraphBuilder:
    """
    Sparse edge-list builder for large item universes.

    Currently supports sparse construction for:
        - log_gcd
        - overlap

    It uses an inverted-index over active modes. For each slice and mode,
    only items active in that mode are paired. This avoids allocating
    O(S*I*I) dense matrices.

    Warning:
        If a mode is active for nearly every item, that mode still induces
        O(k^2) pairs. Use max_items_per_mode, min_exponent, or more selective
        specs to control graph explosion.
    """

    def __init__(
        self,
        mode_weights: Optional[np.ndarray] = None,
        *,
        kernel: ArithmeticKernelName = "log_gcd",
        min_weight: float = 0.0,
        min_exponent: float = 0.0,
        max_items_per_mode: Optional[int] = None,
    ):
        if kernel not in {"log_gcd", "overlap"}:
            raise ValueError(
                "SparseArithmeticGraphBuilder currently supports only "
                "kernel='log_gcd' or kernel='overlap'. Use dense ArithmeticGraphBuilder "
                "for lcm/distance/jaccard/coprime kernels."
            )

        self.mode_weights = mode_weights
        self.kernel = kernel
        self.min_weight = float(min_weight)
        self.min_exponent = float(min_exponent)
        self.max_items_per_mode = max_items_per_mode

        if self.min_weight < 0:
            raise ValueError("min_weight must be nonnegative.")
        if self.min_exponent < 0:
            raise ValueError("min_exponent must be nonnegative.")
        if self.max_items_per_mode is not None and self.max_items_per_mode <= 1:
            raise ValueError("max_items_per_mode must be > 1 when provided.")

    def build(self, encoded: EncodedND) -> EdgeListGraph:
        x = np.asarray(encoded.exponents, dtype=float)
        if x.ndim != 3:
            raise ValueError("encoded.exponents must have shape (S, I, P).")
        if not np.all(np.isfinite(x)):
            raise ValueError("encoded.exponents contains non-finite values.")
        if np.any(x < 0):
            raise ValueError("Fock/valuation exponents must be nonnegative.")

        log_primes = _resolve_mode_weights(encoded, self.mode_weights)
        S, I, P = x.shape

        rows: list[tuple[int, int, int, float]] = []

        for s in range(S):
            accum: dict[tuple[int, int], float] = defaultdict(float)
            xs = x[s]

            for p in range(P):
                active = np.flatnonzero(xs[:, p] > self.min_exponent)
                k = active.size
                if k < 2:
                    continue

                if self.max_items_per_mode is not None and k > self.max_items_per_mode:
                    # Keep strongest exponents for this broad mode.
                    order = np.argsort(xs[active, p])[-self.max_items_per_mode:]
                    active = active[order]
                    k = active.size

                vals = xs[active, p]
                w = log_primes[p]

                for a_pos in range(k - 1):
                    i = int(active[a_pos])
                    vi = vals[a_pos]
                    for b_pos in range(a_pos + 1, k):
                        j = int(active[b_pos])
                        vj = vals[b_pos]
                        contribution = min(vi, vj) * w
                        if contribution > 0:
                            accum[(i, j)] += float(contribution)

            for (i, j), weight in accum.items():
                if weight > self.min_weight:
                    rows.append((s, i, j, weight))

        edges = pd.DataFrame(rows, columns=["slice_idx", "source_idx", "target_idx", "weight"])
        if not edges.empty:
            edges = edges.sort_values(["slice_idx", "source_idx", "target_idx"]).reset_index(drop=True)

        return EdgeListGraph(
            edges=edges,
            n_slices=S,
            n_items=I,
            slices=encoded.slices,
            items=encoded.items,
            kernel=self.kernel,
        )


class SparseSharedSupportGraphBuilder(SparseArithmeticGraphBuilder):
    """
    Backward-compatible alias for SparseArithmeticGraphBuilder.
    """
    pass


def _resolve_mode_weights(encoded: EncodedND, mode_weights: Optional[np.ndarray]) -> np.ndarray:
    if mode_weights is None:
        return np.log(encoded.primes.astype(float))

    weights = np.asarray(mode_weights, dtype=float)
    if weights.shape != (encoded.exponents.shape[2],):
        raise ValueError(f"mode_weights must have shape {(encoded.exponents.shape[2],)}, got {weights.shape}.")
    if not np.all(np.isfinite(weights)):
        raise ValueError("mode_weights contains non-finite values.")
    return weights
