from __future__ import annotations

import numpy as np
import pandas as pd

from .config import FeatureEmitterConfig
from .encoder import EncodedND
from .graph import EdgeListGraph
from .sparse_stability import SparseStabilityResult


class SparseFeatureEmitter:
    """
    Emits generic features from sparse edge-list graphs without dense adjacency.
    """

    def __init__(self, config: FeatureEmitterConfig | None = None):
        self.config = config or FeatureEmitterConfig()

    def transform(
        self,
        encoded: EncodedND,
        graph: EdgeListGraph,
        stability: SparseStabilityResult,
    ) -> pd.DataFrame:
        x = np.asarray(encoded.exponents, dtype=float)
        S, I, _ = x.shape

        if S != graph.n_slices or I != graph.n_items:
            raise ValueError("Encoded tensor and graph shapes are inconsistent.")

        raw = _degrees_from_edges(graph.edges, S, I)
        stable = _degrees_from_edges(stability.stable_edges, S, I)

        rows = []
        prev_x = None
        prev_raw_degree = None
        prev_stable_degree = None

        for s in range(S):
            raw_degree = raw["degree"][s]
            stable_degree = stable["degree"][s]
            raw_edge_count = raw["edge_count"][s]
            stable_edge_count = stable["edge_count"][s]

            total_stable = max(float(stable_degree.sum()), self.config.eps)
            stable_centrality = stable_degree / total_stable
            unstable_mass = np.maximum(raw_degree - stable_degree, 0.0)

            active_mode_count = (x[s] > 0).sum(axis=1)
            fock_mass = x[s].sum(axis=1)

            if prev_x is None:
                defect_count = np.zeros(I)
                defect_mass = np.zeros(I)
                graph_delta = np.zeros(I)
                stable_graph_delta = np.zeros(I)
            else:
                created = np.maximum((x[s] > 0).astype(float) - (prev_x > 0).astype(float), 0.0)
                defect_count = created.sum(axis=1)
                defect_mass = np.maximum(x[s] - prev_x, 0.0).sum(axis=1)
                graph_delta = np.abs(raw_degree - prev_raw_degree)
                stable_graph_delta = np.abs(stable_degree - prev_stable_degree)

            for item_idx, item in enumerate(encoded.items):
                rows.append(
                    {
                        "slice": encoded.slices[s],
                        "item": item,
                        "fock_active_modes": float(active_mode_count[item_idx]),
                        "fock_mass": float(fock_mass[item_idx]),
                        "raw_degree": float(raw_degree[item_idx]),
                        "stable_degree": float(stable_degree[item_idx]),
                        "raw_edge_count": float(raw_edge_count[item_idx]),
                        "stable_edge_count": float(stable_edge_count[item_idx]),
                        "stable_centrality": float(stable_centrality[item_idx]),
                        "unstable_mass": float(unstable_mass[item_idx]),
                        "defect_count": float(defect_count[item_idx]),
                        "defect_mass": float(defect_mass[item_idx]),
                        "graph_delta": float(graph_delta[item_idx]),
                        "stable_graph_delta": float(stable_graph_delta[item_idx]),
                        "spectral_radius": float(stability.spectral_radius[s]),
                        "graph_kernel": getattr(graph, "kernel", "log_gcd"),
                    }
                )

            prev_x = x[s].copy()
            prev_raw_degree = raw_degree.copy()
            prev_stable_degree = stable_degree.copy()

        return pd.DataFrame(rows)


def _degrees_from_edges(edges: pd.DataFrame, n_slices: int, n_items: int) -> dict[str, np.ndarray]:
    degree = np.zeros((n_slices, n_items), dtype=float)
    edge_count = np.zeros((n_slices, n_items), dtype=float)

    if edges.empty:
        return {"degree": degree, "edge_count": edge_count}

    for row in edges.itertuples(index=False):
        s = int(row.slice_idx)
        i = int(row.source_idx)
        j = int(row.target_idx)
        w = float(row.weight)
        degree[s, i] += w
        degree[s, j] += w
        edge_count[s, i] += 1
        edge_count[s, j] += 1

    return {"degree": degree, "edge_count": edge_count}
