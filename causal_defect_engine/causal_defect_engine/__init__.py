"""Causal Defect Algebra Engine package."""

from .core import CausalDefectEngine, DefectState, Defect, CompositionRule, RuntimeConfig, NUMBA_AVAILABLE
from . import numba_kernels
from .pointcloud_spin import PointCloudSpinEngine, PointCloudSpinConfig
from .arithmetic_pleroma import ArithmeticPleromaEngine, ArithmeticConfig, DirichletConvolution, primes_up_to
from .graph_holonomy import GraphHolonomyEngine, GraphHolonomyConfig

from .vm import DefectVM, VMTraceEntry

from .memory import MemoryConfig, MemoryTrace, RecallResult, PleromaMemoryEngine
from .coarse_grain import (
    CoarseGrainConfig, CoarseGrainReport, merge_nearby_defects,
    compress_omega, quotient_by_labels, coarse_grain_state,
)
from .optim import (
    OptimConfig, SearchResult, gf2_nullspace, neutral_cycle_search,
    low_action_path, rank_holonomy_paths, DefectAnnealer, optimize_mps_ordering,
)


from .evolution import (
    EvolutionConfig, EvolutionReport, EvolutionRule, TransportEvolution,
    CurvatureRelaxation, DefectDiffusion, ReactionDiffusionEvolution,
    HamiltonianEvolution, EvolutionPipeline, graph_laplacian_from_edges,
    defect_charges_to_vector, vector_to_defects, evolve_state,
    transport_along_connection, relax_curvature, diffuse_defects,
    evolve_and_recompress,
)

from .tensor_defect import (
    TensorDefectEngine, TensorDefectConfig, TensorDefectState, SimpleMPS,
    dense_to_mps, defect_guided_mps_from_state, estimate_bond_dimensions,
    learn_mps_ordering, connected_correlator,
)

from .state_geometry import (
    StateGeometryConfig, DefectComplexity, CurvatureSparseStateGeometry,
    CurvatureSparseStateEngine, geometry_from_omega, geometry_from_mps,
    geometry_from_ising, fit_curvature_sparse_geometry,
)

from .geometry_compiler import (
    CompilerConfig, GeometryCompileReport, DefectGeometryCompiler,
    learn_adaptive_topology, build_binary_connections, detect_curvature_loops,
    edge_phases_from_omega, gauge_fix_geometry, compile_geometry, complexity_report,
)

__all__ = [
    "CausalDefectEngine",
    "DefectState",
    "Defect",
    "CompositionRule",
    "RuntimeConfig",
    "NUMBA_AVAILABLE",
    "numba_kernels",
    "PointCloudSpinEngine",
    "PointCloudSpinConfig",
    "ArithmeticPleromaEngine",
    "ArithmeticConfig",
    "DirichletConvolution",
    "primes_up_to",
    "GraphHolonomyEngine",
    "GraphHolonomyConfig",
    "DefectVM",
    "VMTraceEntry",
    "TensorDefectEngine",
    "TensorDefectConfig",
    "TensorDefectState",
    "SimpleMPS",
    "dense_to_mps",
    "defect_guided_mps_from_state",
    "estimate_bond_dimensions",
    "learn_mps_ordering",
    "connected_correlator",
    "MemoryConfig",
    "MemoryTrace",
    "RecallResult",
    "PleromaMemoryEngine",
    "CoarseGrainConfig",
    "CoarseGrainReport",
    "merge_nearby_defects",
    "compress_omega",
    "quotient_by_labels",
    "coarse_grain_state",
    "OptimConfig",
    "SearchResult",
    "gf2_nullspace",
    "neutral_cycle_search",
    "low_action_path",
    "rank_holonomy_paths",
    "DefectAnnealer",
    "optimize_mps_ordering",
    "EvolutionConfig",
    "EvolutionReport",
    "EvolutionRule",
    "TransportEvolution",
    "CurvatureRelaxation",
    "DefectDiffusion",
    "ReactionDiffusionEvolution",
    "HamiltonianEvolution",
    "EvolutionPipeline",
    "graph_laplacian_from_edges",
    "defect_charges_to_vector",
    "vector_to_defects",
    "evolve_state",
    "transport_along_connection",
    "relax_curvature",
    "diffuse_defects",
    "evolve_and_recompress",
    "StateGeometryConfig",
    "DefectComplexity",
    "CurvatureSparseStateGeometry",
    "CurvatureSparseStateEngine",
    "geometry_from_omega",
    "geometry_from_mps",
    "geometry_from_ising",
    "fit_curvature_sparse_geometry",
    "CompilerConfig",
    "GeometryCompileReport",
    "DefectGeometryCompiler",
    "learn_adaptive_topology",
    "build_binary_connections",
    "detect_curvature_loops",
    "edge_phases_from_omega",
    "gauge_fix_geometry",
    "compile_geometry",
    "complexity_report",
]
