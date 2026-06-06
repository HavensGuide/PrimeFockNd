# Causal Defect Engine v0.2

Extensible computational substrate for defect-algebra experiments.

## Backends

- `core.py` — shared abstractions, runtime modes, Numba kernels, kNN helpers.
- `pointcloud_spin.py` — point-cloud spin / Madelung flow frame backend.
- `arithmetic_pleroma.py` — arithmetic defect algebra and relation search backend.
- `graph_holonomy.py` — graph connection / holonomy backend.

## Execution modes

Each backend accepts mode controls. The legacy behavior is preserved for testing.

```python
PointCloudSpinEngine(PointCloudSpinConfig(mode="legacy"))    # all-pairs exact score
PointCloudSpinEngine(PointCloudSpinConfig(mode="local"))     # kNN-local score, lazy metric
ArithmeticPleromaEngine(ArithmeticConfig(mode="dense"))      # dense F2 matrix
ArithmeticPleromaEngine(ArithmeticConfig(mode="packed"))     # bit-packed matrix storage
GraphHolonomyEngine(GraphHolonomyConfig(mode="exact"))       # v0.1-style DFS cycles
GraphHolonomyEngine(GraphHolonomyConfig(mode="efficient"))   # chord cycles, capped early
```

Global runtime knobs are also accepted:

```python
engine = PointCloudSpinEngine(mode="local", neighbor_backend="auto", use_numba=True)
```

## Important naming fix

The point-cloud tensor is now called `G_flow`, not a Fisher metric:

\[
G_{flow} = u u^T + \alpha v v^T + \lambda I.
\]

It is a rank-2-plus-floor Madelung/defect frame tensor. The score is recoverable by

\[
\nabla \log \rho = u / D_{eff}.
\]

`PointCloudSpinEngine.recover_score()` returns this. `recover_log_density()` integrates the score over the kNN graph and recovers log density up to an additive constant.

## Quick smoke test

```python
import numpy as np
from causal_defect_engine import PointCloudSpinEngine, PointCloudSpinConfig

X = np.random.default_rng(0).normal(size=(300, 4))
eng = PointCloudSpinEngine(PointCloudSpinConfig(mode="local", n_pca_components=2))
state = eng.run(X)
print(state.summary())
print(eng.recover_score().shape)
```

## Dependencies

Required: `numpy`, `numba`.
Optional but recommended: `scipy` for cKDTree nearest-neighbors, Delaunay triangulation, and sparse density recovery.

## v0.3 additions

### DefectVM

`DefectVM` is a small runtime for programs over a fitted `DefectState` or engine.
It supports built-in operations such as `summary`, `observe`, `compose`,
`project`, `holonomy`, `defects`, and context references via `$name`.

```python
from causal_defect_engine import DefectVM

vm = DefectVM(engine)
out = vm.run([
    {"op": "summary", "save_as": "s"},
    {"op": "defect_table", "args": {"min_strength": 0.2}, "save_as": "defects"},
])
print(vm.context["s"])
```

### TensorDefectEngine

`TensorDefectEngine` bridges defect states to MPS-like compression.  It extracts
or estimates an `Omega`/connected-correlation matrix, learns a chain ordering,
estimates bond dimensions from cross-cut ranks, and can convert small dense
tensors to a dependency-free `SimpleMPS` via TT-SVD.

```python
from causal_defect_engine import TensorDefectEngine, dense_to_mps

engine = TensorDefectEngine(mode="efficient")
engine.fit(samples_x_features).estimate_generators()
engine.build_connection()
engine.detect_defects()
print(engine.tensor_state.summary())

mps = dense_to_mps(tensor, max_bond=16, tol=1e-8)
print(mps.bond_dims, mps.entanglement_profile())
```


## v0.4 performance update

Added `causal_defect_engine.numba_kernels` and wired it into the main hot loops:

- point-cloud local Gaussian score over kNN neighborhoods,
- point-cloud local time-drift regression,
- rank-2 `G_flow = uu^T + alpha vv^T + lambda I` frame extraction without dense eigensolves,
- arithmetic quadratic-relation trial division over a factor base,
- packed F2 row packing and optional Numba RREF for parity dependencies,
- graph degree/imbalance accumulation and dense-phase bulk cycle holonomy when the graph is small enough,
- shared feature-correlation/Omega extraction for tensor defect work.

All accelerated paths keep Python/NumPy fallbacks and can be disabled with `use_numba=False`.


## v0.5 additions

Added the remaining computational-medium modules:

- `memory.py` — Pleroma associative trace memory with write/recall/compose/forget/compress operations.
- `coarse_grain.py` — defect merging, Omega compression, quotienting by labels, and VM registration helpers.
- `optim.py` — GF(2) neutral-cycle search, low-action path search, holonomy path ranking, annealing, and MPS ordering optimization.
- `benchmarks.py` and `benchmarks/run_smoke.py` — lightweight proof-of-work benchmarks.

These modules preserve the existing backend style: explicit data structures, optional Numba acceleration where useful, and dependency-light fallbacks.


## v0.6 addition: Curvature-sparse state geometry

Added `state_geometry.py`, a higher-level representation for compressed
random-looking but generatively structured exponential states.

The central object is `CurvatureSparseStateGeometry`:

\[
\mathfrak D_\Psi=(\Gamma,\mathcal F,\mathcal A,\Omega,\Pi,\mathcal B,\mathcal R,\operatorname{Eval})
\]

It stores an adaptive gluing graph, local fibers, edge connections, sparse
curvature defects, optional MPS data, boundary/generator data, and explicit
evaluation semantics.  It separates storage complexity from evaluation and
learning complexity.

```python
import numpy as np
from causal_defect_engine import CurvatureSparseStateEngine, StateGeometryConfig

# From Omega / connected correlations
Omega = np.eye(8)
Omega[0, 1] = Omega[1, 0] = 0.8
Omega[1, 2] = Omega[2, 1] = 0.7
eng = CurvatureSparseStateEngine(StateGeometryConfig(connection_threshold=0.2))
eng.fit(Omega).estimate_generators()
print(eng.estimate_defect_rank())
print(eng.eval_amplitude([0,1,0,1,0,1,0,1]))

# From a compact Ising/Boltzmann generator
spec = {"type": "ising", "n": 4, "J": {(0,1): 0.7, (1,2): -0.2}, "beta": 0.5}
eng.fit(spec)
print(eng.geometry.summary())
```

Supported construction paths:

- `geometry_from_omega(Omega)` — adaptive graph from connected correlations.
- `geometry_from_mps(mps)` — wraps a `SimpleMPS`; bonds become gluing channels.
- `geometry_from_ising(spec)` — compact generator storage for Boltzmann/phase states.
- dense tensor input — TT-SVD to `SimpleMPS`, then geometry wrapper.

The key rank estimate is exposed through:

```python
geom.estimate_defect_rank(epsilon=1e-8)
```

which reports rough `r_def_store`, `r_def_curvature`, Omega rank, MPS bond rank,
and storage ratio versus the dense state table when enumeration is meaningful.

## v0.7 advanced compiler update

Added `geometry_compiler.py`, which upgrades curvature-sparse state geometry into
an adaptive compiler/runtime layer.

New features:

- `DefectGeometryCompiler` for fitting from Omega, samples/features, MPS, dense tensors, Ising specs, or compact circuit/generator specs.
- Adaptive topology discovery with modes: `threshold`, `mst`, `mst_chords`, `knn`, `legacy`, and `auto`.
- U(1)/Z2-proxy gauge fixing over a spanning forest: tree transport is flattened and only nontrivial residual edges/curvature defects are stored.
- Explicit complexity report separating storage, evaluation, sampling, and learning estimates.
- Compile targets: `geometry`, `defect_state`, and small-state `mps`/`tt` via TT-SVD.
- New Numba kernels for dense threshold edge scans, binary edge-connection construction, triangle curvature candidates, tree-gauge residuals, and batch Ising energies.

Example:

```python
import numpy as np
from causal_defect_engine import DefectGeometryCompiler, CompilerConfig

Omega = np.eye(12)
Omega[0,1] = Omega[1,0] = 0.9
Omega[1,2] = Omega[2,1] = 0.8
Omega[2,0] = Omega[0,2] = -0.7   # signed loop creates residual holonomy proxy

compiler = DefectGeometryCompiler(CompilerConfig(
    mode="efficient",
    topology_method="mst_chords",
    connection_threshold=0.15,
    curvature_tol=0.02,
    gauge_fix=True,
    use_numba=True,
))

geom = compiler.fit(Omega)
print(compiler.report().as_dict())
print(geom.complexity_report())
print(geom.eval_amplitude([0,1,0,1,0,1,0,1,0,1,0,1]))
```

The key compression move is now explicit:

\[
\text{learn connection} \rightarrow \text{gauge-fix flat tree} \rightarrow \text{store residual curvature only.}
\]

## v0.8: Evolution layer

`evolution.py` adds a unified dynamics mechanism for compressed defect geometries:

```python
from causal_defect_engine import (
    DefectState, Defect, TransportEvolution, CurvatureRelaxation,
    DefectDiffusion, EvolutionPipeline, EvolutionConfig
)

state = DefectState(
    generators=np.random.randn(5, 2),
    connection={"edges": [(0,1), (1,2), (2,3), (3,4)], "weights": [1,1,1,1]},
    defects=[Defect(kind="charge", charge=1.0, support=(0,))],
    omega=np.eye(5),
    valid=True,
)

pipe = EvolutionPipeline([
    TransportEvolution(damping=0.01),
    DefectDiffusion(diffusion=0.1),
    CurvatureRelaxation(eta=0.05),
], EvolutionConfig(steps=10, dt=0.1))

state2, report = pipe.evolve(state)
print(report.as_dict())
```

The conceptual update is:

```text
D_{t+dt} = R_epsilon ∘ G ∘ Phi_dt(D_t)
```

where `Phi_dt` is a raw evolution rule, optional gauge fixing removes flat transport, and optional recompression/coarse-graining keeps the state compact.

Available rules:

- `TransportEvolution`: message passing / connection transport over graph-like connections.
- `CurvatureRelaxation`: relaxes curvature/connection inconsistency and weak defects.
- `DefectDiffusion`: diffuses defect charge over graph supports.
- `ReactionDiffusionEvolution`: diffusion plus a stable local reaction term.
- `HamiltonianEvolution`: small/medium dense state-vector evolution and imaginary-time evolution.
- `EvolutionPipeline`: composes multiple rules with optional gauge-fix / recompress steps.

Numba-enhanced kernels accelerate edge diffusion, feature transport, dense Euler steps, and curvature-energy diagnostics when Numba is available.
