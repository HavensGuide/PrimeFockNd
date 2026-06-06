import numpy as np
import pandas as pd

from prime_fock_nd import (
    NumericBinningSpec,
    CategoricalSpec,
    PrimeFockNDConfig,
    PrimeFockNDEncoder,
    SharedSupportGraphBuilder,
    SparseSharedSupportGraphBuilder,
    SparseStabilityFilter,
    SparseFeatureEmitter,
    StabilityConfig,
)


def _fixture():
    df = pd.DataFrame({
        "time": [1, 1, 1, 2, 2, 2],
        "entity": ["A", "B", "C", "A", "B", "C"],
        "x": [0.1, 0.2, 2.0, 0.15, 1.8, 2.2],
        "kind": ["red", "red", "blue", "red", "blue", "blue"],
    })
    cfg = PrimeFockNDConfig(
        index_cols=["time", "entity"],
        item_col="entity",
        slice_cols=["time"],
        specs=[
            NumericBinningSpec("x", bins=[1.0]),
            CategoricalSpec("kind"),
        ],
    )
    return PrimeFockNDEncoder(cfg).fit_transform(df)


def test_sparse_builder_matches_dense_edges_for_fixture():
    encoded = _fixture()
    dense = SharedSupportGraphBuilder().build(encoded)
    sparse = SparseSharedSupportGraphBuilder().build(encoded)

    reconstructed = np.zeros_like(dense.adjacency)
    for row in sparse.edges.itertuples(index=False):
        s = int(row.slice_idx)
        i = int(row.source_idx)
        j = int(row.target_idx)
        reconstructed[s, i, j] = row.weight
        reconstructed[s, j, i] = row.weight

    assert np.allclose(dense.adjacency, reconstructed)


def test_sparse_stability_and_features():
    encoded = _fixture()
    sparse = SparseSharedSupportGraphBuilder().build(encoded)
    stable = SparseStabilityFilter(StabilityConfig(local_bond_threshold=10.0)).fit_transform(sparse)
    features = SparseFeatureEmitter().transform(encoded, sparse, stable)

    assert len(features) == 6
    assert {"slice", "item", "stable_degree", "defect_count", "spectral_radius"}.issubset(features.columns)
    assert stable.stable_edges.shape[0] <= sparse.edges.shape[0]
