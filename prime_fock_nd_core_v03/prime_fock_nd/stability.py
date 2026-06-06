from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .config import StabilityConfig
from .graph import SharedSupportGraph
from .validation import ensure_square_matrix, ensure_symmetric


@dataclass(frozen=True)
class StabilityResult:
    stable_adjacency: np.ndarray
    spectral_radius: np.ndarray
    local_stability: np.ndarray
    slices: list[tuple]
    items: list


class StabilityFilter:
    """
    General cone-style stability filter.

    Given H_int(s) and optional background H0(s), computes
    C = H0^-1/2 H_int H0^-1/2.

    If no background is passed, identity is used.
    """

    def __init__(self, config: StabilityConfig | None = None):
        self.config = config or StabilityConfig()

    def fit_transform(
        self,
        graph: SharedSupportGraph,
        background_metric: Optional[np.ndarray] = None,
    ) -> StabilityResult:
        adjacency = np.asarray(graph.adjacency, dtype=float)
        if adjacency.ndim != 3:
            raise ValueError("graph.adjacency must have shape (S, I, I).")
        S, I, _ = adjacency.shape

        if background_metric is None:
            H0_all = np.broadcast_to(np.eye(I), (S, I, I)).copy()
        else:
            H0_all = np.asarray(background_metric, dtype=float)
            if H0_all.shape == (I, I):
                H0_all = np.broadcast_to(H0_all, (S, I, I)).copy()
            elif H0_all.shape != (S, I, I):
                raise ValueError(f"background_metric must have shape {(I, I)} or {(S, I, I)}, got {H0_all.shape}.")

        stable = np.zeros_like(adjacency)
        local_stability = np.zeros_like(adjacency)
        spectral_radius = np.zeros(S, dtype=float)

        for s in range(S):
            H_int = ensure_square_matrix(adjacency[s], I, "H_int")
            H0 = ensure_square_matrix(H0_all[s], I, "H0")
            ensure_symmetric(H_int, "H_int")
            ensure_symmetric(H0, "H0")

            H0 = H0 + self.config.ridge * np.eye(I)
            inv_half = _inverse_symmetric_sqrt(H0)
            C = inv_half @ H_int @ inv_half
            C = (C + C.T) / 2.0

            eigvals = np.linalg.eigvalsh(C)
            spectral_radius[s] = float(np.max(np.abs(eigvals))) if eigvals.size else 0.0

            diag = np.sqrt(np.outer(np.diag(H0), np.diag(H0)))
            diag = np.maximum(diag, self.config.ridge)
            local = np.abs(H_int) / diag
            local_stability[s] = local

            keep = (local <= self.config.local_bond_threshold).astype(float)
            keep *= (1.0 - np.eye(I))
            stable[s] = H_int * keep

            if spectral_radius[s] > self.config.max_spectral_radius:
                shrink = self.config.max_spectral_radius / (spectral_radius[s] + self.config.ridge)
                stable[s] *= shrink

        return StabilityResult(
            stable_adjacency=stable,
            spectral_radius=spectral_radius,
            local_stability=local_stability,
            slices=graph.slices,
            items=graph.items,
        )


def _inverse_symmetric_sqrt(matrix: np.ndarray) -> np.ndarray:
    eigvals, eigvecs = np.linalg.eigh((matrix + matrix.T) / 2.0)
    eigvals = np.maximum(eigvals, 1e-12)
    return eigvecs @ np.diag(1.0 / np.sqrt(eigvals)) @ eigvecs.T
