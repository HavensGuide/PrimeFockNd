import numpy as np
import pandas as pd

from prime_fock_nd import (
    NumericBinningSpec,
    CategoricalSpec,
    PrimeFockNDConfig,
    PrimeFockNDEncoder,
    SharedSupportGraphBuilder,
    StabilityFilter,
    StabilityConfig,
    GenericFeatureEmitter,
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
    encoded = PrimeFockNDEncoder(cfg).fit_transform(df)
    graph = SharedSupportGraphBuilder().build(encoded)
    return encoded, graph


def test_graph_shared_support():
    encoded, graph = _fixture()
    assert graph.adjacency.shape == (2, 3, 3)
    assert graph.adjacency[0, 0, 1] > 0
    assert graph.adjacency[0, 0, 2] == 0
    assert np.allclose(graph.adjacency, np.swapaxes(graph.adjacency, 1, 2))


def test_stability_and_features():
    encoded, graph = _fixture()
    stable = StabilityFilter(StabilityConfig(local_bond_threshold=10.0)).fit_transform(graph)
    features = GenericFeatureEmitter().transform(encoded, graph, stable)

    assert stable.stable_adjacency.shape == graph.adjacency.shape
    assert stable.spectral_radius.shape == (2,)
    assert len(features) == 6
    assert {"slice", "item", "stable_degree", "defect_count", "spectral_radius"}.issubset(features.columns)
    assert features.loc[features["slice"] == (1,), "defect_count"].sum() == 0
    assert features.loc[features["slice"] == (2,), "defect_count"].sum() > 0
