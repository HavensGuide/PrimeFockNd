from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

from .config import StabilityConfig
from .graph import EdgeListGraph


@dataclass(frozen=True)
class SparseStabilityResult:
    stable_edges: pd.DataFrame
    spectral_radius: np.ndarray
    n_slices: int
    n_items: int
    slices: list[tuple]
    items: list


class SparseStabilityFilter:
    """
    Edge-list stability filter for large graphs.

    This avoids dense adjacency allocation. It supports:
    - diagonal background capacity via background_diag
    - local bond filtering
    - approximate spectral-radius shrinkage by power iteration

    For exact dense generalized eigenvalue behavior, use StabilityFilter.
    """

    def __init__(self, config: StabilityConfig | None = None, power_iter: int = 30):
        self.config = config or StabilityConfig()
        self.power_iter = int(power_iter)
        if self.power_iter <= 0:
            raise ValueError("power_iter must be positive.")

    def fit_transform(
        self,
        graph: EdgeListGraph,
        background_diag: Optional[np.ndarray] = None,
    ) -> SparseStabilityResult:
        if background_diag is None:
            diag = np.ones((graph.n_slices, graph.n_items), dtype=float)
        else:
            diag = np.asarray(background_diag, dtype=float)
            if diag.shape == (graph.n_items,):
                diag = np.broadcast_to(diag, (graph.n_slices, graph.n_items)).copy()
            elif diag.shape != (graph.n_slices, graph.n_items):
                raise ValueError(
                    f"background_diag must have shape {(graph.n_items,)} or "
                    f"{(graph.n_slices, graph.n_items)}, got {diag.shape}."
                )
            if not np.all(np.isfinite(diag)):
                raise ValueError("background_diag contains non-finite values.")
            if np.any(diag <= 0):
                raise ValueError("background_diag must be strictly positive.")

        if graph.edges.empty:
            return SparseStabilityResult(
                stable_edges=graph.edges.copy(),
                spectral_radius=np.zeros(graph.n_slices),
                n_slices=graph.n_slices,
                n_items=graph.n_items,
                slices=graph.slices,
                items=graph.items,
            )

        edges = graph.edges.copy()
        spectral = np.zeros(graph.n_slices, dtype=float)
        kept_parts = []

        for s in range(graph.n_slices):
            part = edges[edges["slice_idx"] == s].copy()
            if part.empty:
                continue

            i = part["source_idx"].to_numpy(dtype=int)
            j = part["target_idx"].to_numpy(dtype=int)
            w = part["weight"].to_numpy(dtype=float)

            cap = np.sqrt(diag[s, i] * diag[s, j])
            cap = np.maximum(cap, self.config.ridge)
            local = np.abs(w) / cap
            part["local_stability"] = local

            keep = local <= self.config.local_bond_threshold
            part = part.loc[keep].copy()

            if part.empty:
                continue

            i = part["source_idx"].to_numpy(dtype=int)
            j = part["target_idx"].to_numpy(dtype=int)
            w = part["weight"].to_numpy(dtype=float)

            normalized_w = w / np.sqrt(np.maximum(diag[s, i] * diag[s, j], self.config.ridge))
            rho = _edge_power_spectral_radius(graph.n_items, i, j, normalized_w, self.power_iter)
            spectral[s] = rho

            if rho > self.config.max_spectral_radius:
                shrink = self.config.max_spectral_radius / (rho + self.config.ridge)
                part["weight"] = part["weight"] * shrink

            kept_parts.append(part)

        if kept_parts:
            stable_edges = pd.concat(kept_parts, ignore_index=True)
            stable_edges = stable_edges.sort_values(["slice_idx", "source_idx", "target_idx"]).reset_index(drop=True)
        else:
            stable_edges = edges.iloc[0:0].copy()
            if "local_stability" not in stable_edges.columns:
                stable_edges["local_stability"] = []

        return SparseStabilityResult(
            stable_edges=stable_edges,
            spectral_radius=spectral,
            n_slices=graph.n_slices,
            n_items=graph.n_items,
            slices=graph.slices,
            items=graph.items,
        )


def _edge_power_spectral_radius(n: int, i: np.ndarray, j: np.ndarray, w: np.ndarray, iters: int) -> float:
    if w.size == 0:
        return 0.0

    x = np.ones(n, dtype=float) / np.sqrt(max(n, 1))
    last_norm = 0.0

    for _ in range(iters):
        y = np.zeros(n, dtype=float)
        np.add.at(y, i, w * x[j])
        np.add.at(y, j, w * x[i])

        norm = float(np.linalg.norm(y))
        if norm <= 1e-12:
            return 0.0
        x = y / norm
        last_norm = norm

    return float(last_norm)
