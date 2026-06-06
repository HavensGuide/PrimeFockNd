import numpy as np
import pandas as pd

from prime_fock_nd import (
    NumericBinningSpec,
    PrimeFockNDConfig,
    PrimeFockNDEncoder,
    ArithmeticGraphBuilder,
    SparseArithmeticGraphBuilder,
    ValuationArithmetic,
)


def _encoded():
    df = pd.DataFrame({
        "time": [1, 1, 1],
        "entity": ["A", "B", "C"],
        "x": [0.1, 0.2, 2.0],
        "y": [5.0, 9.0, 5.0],
    })
    cfg = PrimeFockNDConfig(
        index_cols=["time", "entity"],
        item_col="entity",
        slice_cols=["time"],
        specs=[
            NumericBinningSpec("x", bins=[1.0]),
            NumericBinningSpec("y", bins=[7.0]),
        ],
    )
    return PrimeFockNDEncoder(cfg).fit_transform(df)


def test_dense_graph_is_log_gcd():
    encoded = _encoded()
    graph = ArithmeticGraphBuilder(kernel="log_gcd").build(encoded)

    logs = np.log(encoded.primes.astype(float))
    a = encoded.exponents[0, 0]
    b = encoded.exponents[0, 1]
    expected = ValuationArithmetic.log_gcd(a, b, logs)

    assert graph.kernel == "log_gcd"
    assert np.isclose(graph.adjacency[0, 0, 1], expected)


def test_distance_kernel_differs_from_gcd():
    encoded = _encoded()
    gcd_graph = ArithmeticGraphBuilder(kernel="log_gcd").build(encoded)
    dist_graph = ArithmeticGraphBuilder(kernel="arithmetic_distance").build(encoded)

    assert gcd_graph.adjacency.shape == dist_graph.adjacency.shape
    assert not np.allclose(gcd_graph.adjacency, dist_graph.adjacency)


def test_sparse_log_gcd_matches_dense_log_gcd():
    encoded = _encoded()
    dense = ArithmeticGraphBuilder(kernel="log_gcd").build(encoded)
    sparse = SparseArithmeticGraphBuilder(kernel="log_gcd").build(encoded)

    reconstructed = np.zeros_like(dense.adjacency)
    for row in sparse.edges.itertuples(index=False):
        s = int(row.slice_idx)
        i = int(row.source_idx)
        j = int(row.target_idx)
        reconstructed[s, i, j] = row.weight
        reconstructed[s, j, i] = row.weight

    assert np.allclose(dense.adjacency, reconstructed)
