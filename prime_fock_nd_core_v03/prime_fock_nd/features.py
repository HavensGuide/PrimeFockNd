from __future__ import annotations

import numpy as np
import pandas as pd

from .config import FeatureEmitterConfig
from .encoder import EncodedND
from .graph import SharedSupportGraph
from .stability import StabilityResult


class GenericFeatureEmitter:
    """
    Emits model-ready generic features from encoded Fock geometry.

    This is output/task agnostic. It creates features; it does not train labels.
    """

    def __init__(self, config: FeatureEmitterConfig | None = None):
        self.config = config or FeatureEmitterConfig()

    def transform(
        self,
        encoded: EncodedND,
        graph: SharedSupportGraph,
        stability: StabilityResult,
    ) -> pd.DataFrame:
        x = np.asarray(encoded.exponents, dtype=float)
        A = np.asarray(graph.adjacency, dtype=float)
        Sg = np.asarray(stability.stable_adjacency, dtype=float)

        if x.shape[:2] != A.shape[:2] or A.shape != Sg.shape:
            raise ValueError("Encoded tensor, graph, and stability result shapes are inconsistent.")

        S, I, P = x.shape
        rows = []

        prev_x = None
        prev_A = None
        prev_Sg = None

        for s in range(S):
            raw_degree = A[s].sum(axis=1)
            stable_degree = Sg[s].sum(axis=1)
            raw_edge_count = (A[s] > 0).sum(axis=1)
            stable_edge_count = (Sg[s] > 0).sum(axis=1)

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
                graph_delta = np.abs(A[s] - prev_A).sum(axis=1)
                stable_graph_delta = np.abs(Sg[s] - prev_Sg).sum(axis=1)

            for i, item in enumerate(encoded.items):
                row = {
                    "slice": encoded.slices[s],
                    "item": item,
                    "fock_active_modes": float(active_mode_count[i]),
                    "fock_mass": float(fock_mass[i]),
                    "raw_degree": float(raw_degree[i]),
                    "stable_degree": float(stable_degree[i]),
                    "raw_edge_count": float(raw_edge_count[i]),
                    "stable_edge_count": float(stable_edge_count[i]),
                    "stable_centrality": float(stable_centrality[i]),
                    "unstable_mass": float(unstable_mass[i]),
                    "defect_count": float(defect_count[i]),
                    "defect_mass": float(defect_mass[i]),
                    "graph_delta": float(graph_delta[i]),
                    "stable_graph_delta": float(stable_graph_delta[i]),
                    "spectral_radius": float(stability.spectral_radius[s]),
                        "graph_kernel": getattr(graph, "kernel", "log_gcd"),
                }
                rows.append(row)

            prev_x = x[s].copy()
            prev_A = A[s].copy()
            prev_Sg = Sg[s].copy()

        return pd.DataFrame(rows)
