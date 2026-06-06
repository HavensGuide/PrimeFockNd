# Prime Fock ND Core v0.3

A general-purpose, n-dimensional prime-valuation arithmetic and stability-filtering core.

This package is **domain agnostic**. It is not tied to markets, materials, text, biology, images, or any specific data source. It provides the reusable core for mapping arbitrary records into sparse arithmetic Fock modes, constructing shared-support graphs, filtering unstable relationships, and emitting generic downstream features.

No real data and no synthetic data generator are included.

## What it does

Given a dataframe with arbitrary index dimensions, for example:

```text
time, entity
experiment, sample
patient, visit, tissue
document, section, paragraph
x, y, z, channel
asset, date
```

it maps input columns into prime/Fock exponent vectors:

```math
v(x) = (v_{p_1}, v_{p_2}, \dots, v_{p_K})
```

and computes shared-support structure:

```math
E_{ij}=\sum_p \min(v_{ip},v_{jp}) \log p
```

Then it optionally applies a cone-style stability filter:

```math
C = H_0^{-1/2} H_{int} H_0^{-1/2}
```

and outputs generic features such as:

- raw degree
- stable degree
- stable centrality
- unstable mass
- defect creation count
- defect mass
- graph delta
- stable graph delta
- spectral radius
- active mode count
- Fock mass

## Install

```bash
pip install -e .
pip install -e ".[test]"
```

## Run tests

```bash
pytest -q
```

## Minimal example

```python
import pandas as pd
from prime_fock_nd import (
    NumericBinningSpec,
    CategoricalSpec,
    PrimeFockNDConfig,
    PrimeFockNDEncoder,
    SharedSupportGraphBuilder,
    StabilityFilter,
    GenericFeatureEmitter,
)

df = pd.DataFrame({
    "time": [1, 1, 1, 2, 2, 2],
    "entity": ["A", "B", "C", "A", "B", "C"],
    "x": [0.1, 0.2, 2.0, 0.15, 1.8, 2.2],
    "kind": ["red", "red", "blue", "red", "blue", "blue"],
})

config = PrimeFockNDConfig(
    index_cols=["time", "entity"],
    item_col="entity",
    slice_cols=["time"],
    specs=[
        NumericBinningSpec("x", bins=[0.0, 1.0, 2.0]),
        CategoricalSpec("kind"),
    ],
)

encoded = PrimeFockNDEncoder(config).fit_transform(df)
graph = SharedSupportGraphBuilder().build(encoded)
stable = StabilityFilter().fit_transform(graph)
features = GenericFeatureEmitter().transform(encoded, graph, stable)

print(features)
```

## Design

The core separates five concerns:

1. **Spec layer**: declare how arbitrary columns become modes.
2. **Encoder layer**: dataframe to `(slice, item, mode)` Fock tensor.
3. **Graph layer**: shared arithmetic support graph.
4. **Stability layer**: risk/background-normalized filtering.
5. **Emitter layer**: generic model-ready output features.

This lets you later attach any data adapter or task-specific model without rewriting the arithmetic geometry.


## Scaling path

For modest item counts, the dense builder is simplest:

```python
graph = SharedSupportGraphBuilder().build(encoded)
```

For large item counts, use the sparse edge-list builder:

```python
from prime_fock_nd import (
    SparseSharedSupportGraphBuilder,
    SparseStabilityFilter,
    SparseFeatureEmitter,
)

graph = SparseSharedSupportGraphBuilder(
    min_weight=0.0,
    max_items_per_mode=500,  # optional guard against broad-mode explosion
).build(encoded)

stable = SparseStabilityFilter().fit_transform(graph)
features = SparseFeatureEmitter().transform(encoded, graph, stable)
```

Dense complexity:

```text
memory = O(S * I^2)
```

Sparse edge-list complexity:

```text
memory = O(E)
```

where `E` is the number of retained shared-support edges.

The sparse path still cannot defeat a fully dense mode where every item shares the same active prime. Use selective feature specs, `max_items_per_mode`, or mode-frequency pruning in those cases.


## v0.3 arithmetic upgrade

Earlier versions used the graph formula:

```math
E_{ij}=\sum_k \min(v_{ik},v_{jk})\log p_k
```

This is now explicitly treated as the **log-gcd kernel** over implicit prime-product states:

```math
n_i=\prod_k p_k^{v_{ik}}
```

so:

```math
E_{ij}=\log\gcd(n_i,n_j)
```

The package now exposes valuation arithmetic directly:

```python
from prime_fock_nd import ValuationArithmetic

ValuationArithmetic.multiply(a, b)            # a + b
ValuationArithmetic.gcd(a, b)                 # min(a, b)
ValuationArithmetic.lcm(a, b)                 # max(a, b)
ValuationArithmetic.quotient(a, b)            # a - b
ValuationArithmetic.arithmetic_distance(a,b,log_primes)  # log(lcm/gcd)
```

Dense graph kernels:

```python
from prime_fock_nd import ArithmeticGraphBuilder

graph = ArithmeticGraphBuilder(kernel="log_gcd").build(encoded)
graph = ArithmeticGraphBuilder(kernel="log_lcm").build(encoded)
graph = ArithmeticGraphBuilder(kernel="arithmetic_distance").build(encoded)
graph = ArithmeticGraphBuilder(kernel="jaccard_gcd_lcm").build(encoded)
graph = ArithmeticGraphBuilder(kernel="coprime_indicator").build(encoded)
```

Sparse graph construction currently supports `log_gcd` / `overlap`, because those can be built efficiently through inverted active-mode indices.
