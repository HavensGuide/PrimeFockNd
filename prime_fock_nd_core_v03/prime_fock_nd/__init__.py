"""
Prime-Fock signal-extraction and dynamics package.

Pipeline:
    specs        column -> Fock-mode mappings
    encoder      dataframe -> (slices, items, prime-modes) valuation tensor
    arithmetic   valuation-space algebra + pairwise arithmetic kernels
    graph        valuation tensor -> relatedness graph (dense or sparse)
    stability    cone-style stability filter (dense or sparse)
    features     graph + stability -> model-ready feature rows
    dynamics     density (diffusion) and amplitude (wave/interference) models
"""

from .specs import (
    BaseModeSpec,
    NumericBinningSpec,
    CategoricalSpec,
    MultiHotSpec,
    PreEncodedSpec,
)
from .config import PrimeFockNDConfig, StabilityConfig, FeatureEmitterConfig
from .primes import first_primes
from .encoder import EncodedND, PrimeFockNDEncoder

# Arithmetic valuation algebra + kernels
from .arithmetic import (
    ValuationArithmetic,
    ArithmeticKernelResult,
    ArithmeticKernelName,
    pairwise_valuation_kernel,
)

# Graph builders — dense and sparse
from .graph import (
    SharedSupportGraph,
    EdgeListGraph,
    ArithmeticGraphBuilder,
    SharedSupportGraphBuilder,
    SparseArithmeticGraphBuilder,
    SparseSharedSupportGraphBuilder,
)

# Stability filters — dense and sparse
from .stability import StabilityResult, StabilityFilter
from .sparse_stability import SparseStabilityResult, SparseStabilityFilter

# Feature emitters — dense and sparse
from .features import GenericFeatureEmitter
from .sparse_features import SparseFeatureEmitter


# Dynamics — diffusion, configurable wave/open-flow dynamics, and features
from .dynamics_energy_flux import (
    GraphDiffusionDynamics,
    GraphWaveDynamics,
    DynamicsFeatureEmitter,
    PreparedOperator,
    GeneratorName,
)

__all__ = [
    # specs
    "BaseModeSpec",
    "NumericBinningSpec",
    "CategoricalSpec",
    "MultiHotSpec",
    "PreEncodedSpec",
    # config
    "PrimeFockNDConfig",
    "StabilityConfig",
    "FeatureEmitterConfig",
    # primes / encoder
    "first_primes",
    "EncodedND",
    "PrimeFockNDEncoder",
    # arithmetic
    "ValuationArithmetic",
    "ArithmeticKernelResult",
    "ArithmeticKernelName",
    "pairwise_valuation_kernel",
    # graph
    "SharedSupportGraph",
    "EdgeListGraph",
    "ArithmeticGraphBuilder",
    "SharedSupportGraphBuilder",
    "SparseArithmeticGraphBuilder",
    "SparseSharedSupportGraphBuilder",
    # stability
    "StabilityResult",
    "StabilityFilter",
    "SparseStabilityResult",
    "SparseStabilityFilter",
    # features
    "GenericFeatureEmitter",
    "SparseFeatureEmitter",
    # dynamics
    "GraphDiffusionDynamics",
    "GraphWaveDynamics",
    "DynamicsFeatureEmitter",
    "PreparedOperator",
    "GeneratorName",
]
