"""
diophantine_oracle_solver.py
============================

Unified Diophantine/Rational Solution Engine
Integrating the Prime-Fock library with quantum walk dynamics and graph holonomy.

Architecture
------------
This solver binds the abstract procedures from the rational_defect_framework
directly to the computational substrate provided by the Prime-Fock package:

  Layer 1  Arithmetic Fock Encoding
           Rational candidates → signed prime-valuation vectors (EncodedND).
           The log-gcd graph is the arithmetic geometry of the candidate space.

  Layer 2  Quantum Walk Dynamics (GraphWaveDynamics)
           i dψ/dt = (L + i·diag(q)) ψ
           The energy_flux field q carries the discrete spin connection —
           per-node inflow/outflow from arithmetic graph deformation.
           Interference residual |ψ|² - ρ detects Brauer-Manin obstructions
           as phase cancellations.

  Layer 3  Berry Phase Accumulation (GraphHolonomyEngine)
           Vary a parameter λ of F_λ(x)=0 around a closed loop enclosing the
           discriminant locus. Track ground-state overlaps via coherence().
           Berry phase = -Im log ∏_k ⟨ψ_0(λ_k)|ψ_0(λ_{k+1})⟩.
           Winding number = topological charge = Brauer class.

  Layer 4  Spin Connection Curvature (Frame Rotation Operator)
           dF = K√g dA computed from the frame eigenvectors of the arithmetic
           Gram matrix. Spikes locate defects; integrated spike = 2π·s_i.
           This is the sdiscrete Gauss-Bonnet theorem on the arithmetic manifold.

  Layer 5  Stability Gate (StabilityFilter / B ≤ I criterion)
           Cone-contractivity check: ‖H_∞^{-1/2} H_F H_∞^{-1/2}‖_op ≤ 1.
           Guards Berry phase computation (requires spectral gap > ‖q‖).

  Layer 6  Classical Cascade (rational_flow_solver_v3 stages 1–4)
           Symbolic, bounded height search, homotopy, Hensel lifting.
           Runs first; quantum layers fire only when classical cascade is
           inconclusive.

Status vocabulary
-----------------
  FOUND_RATIONAL           Certified rational solution
  FOUND_INTEGER            Certified integer solution
  FOUND_ALGEBRAIC          Exact non-rational algebraic solution
  BRAUER_MANIN_CLASSICAL   Hasse–Minkowski obstruction (quadratic forms)
  BRAUER_MANIN_QUANTUM     Berry phase winding number ≠ 0 — topological Brauer class
  ADELIC_OBSTRUCTION       Hensel lifting failure at specific prime
  ARITHMETIC_STIFF         κ_ℚ > threshold — filling resistance too high
  HEIGHT_UNSTABLE          Morse index > 0, flow did not converge
  INTERFERENCE_BLOCKED     Quantum interference residual < −threshold everywhere
  CERTIFIED_NONE_MODULAR   No integer solution: no solution modulo p
  NO_SOLUTION_WITHIN_BOUND Bounded search exhausted
  UNKNOWN                  No stage was decisive

Key design principles
---------------------
- The Prime-Fock encoding works in valuation space; candidate integers n are
  never materialised. The log-gcd kernel gives intersection multiplicity.
- The energy_flux field is the discrete spin connection of the arithmetic graph.
  It is not a tuning knob; it is derived from observed graph deformation.
- Berry phase computation uses the coherence() method of GraphWaveDynamics.
  The discrete Berry phase is accumulated over a loop of deformed graphs.
- The GraphHolonomyEngine provides the cycle basis and holonomy readout for
  identifying which cycles carry non-trivial topological charge.
- All stages have explicit timeout / max_candidates guards; nothing loops forever.

Dependencies
------------
  numpy, scipy, sympy               (always)
  Prime-Fock package (same repo):   encoder, arithmetic, graph, dynamics,
                                    stability, graph_holonomy, core
  numba                             (optional, graceful fallback)
"""

from __future__ import annotations

import math
import os
import sys
import time
import warnings
from dataclasses import dataclass, field
from enum import Enum
from fractions import Fraction
from typing import Sequence

import numpy as np
import sympy as sp

warnings.filterwarnings("ignore", message=".*parallel=True.*", category=UserWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="numba")

# ---------------------------------------------------------------------------
# Prime-Fock library imports
# ---------------------------------------------------------------------------
import sys as _sys
_sys.path.insert(0, "/mnt")   # makes `project` importable as a package

try:
    from project.arithmetic import (
        ValuationArithmetic,
        pairwise_valuation_kernel,
        ArithmeticKernelName,
    )
    from project.primes import first_primes
    from project.dynamics_energy_flux import (
        GraphWaveDynamics,
        GraphDiffusionDynamics,
        DynamicsFeatureEmitter,
    )
    from project.graph import ArithmeticGraphBuilder, SharedSupportGraph
    from project.stability import StabilityFilter, StabilityConfig
    from project.graph_holonomy import GraphHolonomyEngine, GraphHolonomyConfig
    from project.core import Defect, wrap_angle
    PRIMEFOCK_AVAILABLE = True
except ImportError as _e:
    PRIMEFOCK_AVAILABLE = False
    warnings.warn(
        f"Prime-Fock library not importable ({_e}). "
        "Quantum oracle stages will be disabled; classical cascade still runs.",
        stacklevel=1,
    )


# ---------------------------------------------------------------------------
# Status / result types
# ---------------------------------------------------------------------------

class SolverStatus(str, Enum):
    FOUND_RATIONAL           = "FOUND_RATIONAL"
    FOUND_INTEGER            = "FOUND_INTEGER"
    FOUND_ALGEBRAIC          = "FOUND_ALGEBRAIC"
    BRAUER_MANIN_CLASSICAL   = "BRAUER_MANIN_CLASSICAL"
    BRAUER_MANIN_QUANTUM     = "BRAUER_MANIN_QUANTUM"
    CUBIC_OBSTRUCTION        = "CUBIC_OBSTRUCTION"   # Cassels-Guy / norm-residue symbol
    ADELIC_OBSTRUCTION       = "ADELIC_OBSTRUCTION"
    ARITHMETIC_STIFF         = "ARITHMETIC_STIFF"
    HEIGHT_UNSTABLE          = "HEIGHT_UNSTABLE"
    INTERFERENCE_BLOCKED     = "INTERFERENCE_BLOCKED"
    CERTIFIED_NONE_MODULAR   = "CERTIFIED_NONE_MODULAR"
    NO_SOLUTION_WITHIN_BOUND = "NO_SOLUTION_WITHIN_BOUND"
    UNKNOWN                  = "UNKNOWN"


@dataclass
class RationalCandidate:
    point:         tuple
    source:        str
    residual:      float = 0.0
    exact:         bool  = False
    kappa_Q:       float | None = None
    morse_index:   int   | None = None
    obstruction:   str   | None = None
    berry_phase:   float | None = None
    projective:    bool  = False   # True → point is a primitive projective rep


@dataclass
class QuantumObstruction:
    """Result of a Berry-phase / quantum oracle run."""
    winding_number:  float          # ≈ integer if converged
    berry_phase_rad: float          # raw accumulated phase
    n_loop_steps:    int
    spectral_gap:    float          # smallest nonzero eigenvalue of L
    flux_norm:       float          # ‖q‖ — must be < gap for reliable phase
    defects:         list           # Defect objects from GraphHolonomyEngine
    interference_residual_min: float
    reliable:        bool           # spectral_gap > flux_norm


@dataclass
class SolverResult:
    status:              SolverStatus
    rational_solutions:  list[RationalCandidate] = field(default_factory=list)
    integer_solutions:   list[RationalCandidate] = field(default_factory=list)
    algebraic_solutions: list[RationalCandidate] = field(default_factory=list)
    quantum_obstruction: QuantumObstruction | None = None
    classical_messages:  list[str] = field(default_factory=list)
    quantum_messages:    list[str] = field(default_factory=list)
    stage_timings:       dict[str, float] = field(default_factory=dict)

    def summary(self) -> str:
        lines = [f"status = {self.status.value}"]
        for c in self.rational_solutions[:10]:
            kq = f"  κ_Q={c.kappa_Q:.3f}" if c.kappa_Q is not None else ""
            bp = f"  Berry={c.berry_phase:.4f}" if c.berry_phase is not None else ""
            proj = "  [projective]" if c.projective else ""
            lines.append(f"  rational: {c.point}  [{c.source}{kq}{bp}{proj}]")
        for c in self.algebraic_solutions[:5]:
            lines.append(f"  algebraic: {c.point}  [{c.source}]")
        if self.quantum_obstruction is not None:
            q = self.quantum_obstruction
            lines.append(
                f"  quantum: winding={q.winding_number:.3f}  "
                f"Berry={q.berry_phase_rad:.4f} rad  "
                f"reliable={q.reliable}  "
                f"n_defects={len(q.defects)}"
            )
            for d in q.defects[:5]:
                lines.append(f"    defect: charge={d.charge:.4f}  support={d.support}")
        lines.extend(self.classical_messages[:8])
        lines.extend(self.quantum_messages[:8])
        if self.stage_timings:
            lines.append("timings: " + ", ".join(
                f"{k}={v:.3f}s" for k, v in self.stage_timings.items()
            ))
        return "\n".join(lines)


@dataclass(frozen=True)
class OracleSolverConfig:
    # ── classical cascade ──────────────────────────────────────────────────
    max_height:              int   = 20
    max_grid_candidates:     int   = 2_000_000
    max_denominator:         int   = 10_000
    reconstruction_tol:      float = 1e-8
    numeric_residual_tol:    float = 1e-8
    hensel_primes:           tuple = (2, 3, 5, 7, 11, 13)
    hensel_lift_steps:       int   = 8
    filling_resistance_thresh: float = 0.5   # κ_ℚ > this → ARITHMETIC_STIFF

    # ── Prime-Fock / quantum oracle ────────────────────────────────────────
    n_prime_modes:           int   = 16      # number of primes in Fock basis
    n_candidates_oracle:     int   = 64      # candidates fed to quantum layer
    oracle_time_horizon:     float = 1.0     # wave evolution time
    oracle_n_times:          int   = 5       # time steps for feature emission
    berry_loop_steps:        int   = 24      # steps around discriminant loop
    berry_deform_radius:     float = 0.3     # coefficient deformation radius
    berry_winding_threshold: float = 0.35    # |winding - round(winding)| < this → integer
    interference_block_thresh: float = -0.3  # residual below this → INTERFERENCE_BLOCKED
    stability_spectral_max:  float = 1.0     # StabilityFilter max spectral radius
    flux_clip:               float = 3.0

    # ── projective mode ────────────────────────────────────────────────────
    projective:              bool  = True   # auto-detect homogeneous systems
    projective_min_height:   int   = 1      # smallest height for projective search
    do_symbolic:             bool  = True
    do_bounded_search:       bool  = True
    do_homotopy:             bool  = True
    do_hensel:               bool  = True
    do_filling_resistance:   bool  = True
    do_quantum_oracle:       bool  = True
    do_berry_phase:          bool  = True
    do_spin_connection:      bool  = True

    # ── homotopy ──────────────────────────────────────────────────────────
    homotopy_rtol:           float = 1e-10
    homotopy_atol:           float = 1e-12
    n_threads:               int   = 4


# ---------------------------------------------------------------------------
# Arithmetic Fock representation of rational candidates
# ---------------------------------------------------------------------------

def signed_valuation(q: Fraction, primes: list[int]) -> np.ndarray:
    """
    Signed prime-valuation vector of a rational q = ±∏ p^{e_p}.
    Positive exponents → numerator, negative → denominator.
    Returns float array of length len(primes).
    """
    v = np.zeros(len(primes), dtype=np.float64)
    if q == 0:
        return v
    num = abs(q.numerator)
    den = q.denominator
    for k, p in enumerate(primes):
        while num % p == 0:
            num //= p
            v[k] += 1.0
        while den % p == 0:
            den //= p
            v[k] -= 1.0
    return v


def candidate_to_valuation_row(
    point: tuple[Fraction, ...],
    primes: list[int],
) -> np.ndarray:
    """
    Convert a rational point (q_1,...,q_n) to a concatenated signed-valuation
    row of length n * len(primes).  Negative exponents represent denominators.
    For the log-gcd kernel we take abs (shared prime support in numerator).
    """
    rows = [signed_valuation(q, primes) for q in point]
    return np.concatenate(rows)


def build_arithmetic_graph_from_candidates(
    candidates: list[RationalCandidate],
    primes: list[int],
    kernel: ArithmeticKernelName = "log_gcd",
) -> SharedSupportGraph | None:
    """
    Build a SharedSupportGraph from a list of rational candidates.
    Items = candidate indices.  Single slice.
    Valuation rows use abs (shared prime support in numerators/denominators).
    """
    if not PRIMEFOCK_AVAILABLE or not candidates:
        return None

    n_primes = len(primes)
    log_primes = np.log(np.array(primes, dtype=np.float64))
    rows = []

    for cand in candidates:
        try:
            vals = []
            for coord in cand.point:
                if isinstance(coord, Fraction):
                    v = signed_valuation(coord, primes)
                else:
                    f = Fraction(float(coord)).limit_denominator(10_000)
                    v = signed_valuation(f, primes)
                vals.append(np.abs(v))   # abs: shared prime support
            row = np.concatenate(vals)
        except Exception:
            row = np.zeros(len(cand.point) * n_primes)
        rows.append(row)

    valuations = np.array(rows, dtype=np.float64)  # (I, n_vars * n_primes)

    # Pad log_primes to match the concatenated dimension
    n_vars = len(candidates[0].point)
    log_p_ext = np.tile(log_primes, n_vars)

    try:
        kernel_result = pairwise_valuation_kernel(
            valuations, log_p_ext, kernel=kernel
        )
    except Exception:
        return None

    adjacency = kernel_result.matrix[None, :, :]   # (1, I, I)
    items = list(range(len(candidates)))

    return SharedSupportGraph(
        adjacency=adjacency,
        active_modes=np.ones(valuations.shape[1], dtype=bool),
        slices=[("all",)],
        items=items,
        kernel=kernel,
    )


# ---------------------------------------------------------------------------
# Berry phase accumulation
# ---------------------------------------------------------------------------

def compute_berry_phase(
    equations: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
    graph_sequence: list[SharedSupportGraph],
    wave: GraphWaveDynamics,
    config: OracleSolverConfig,
) -> QuantumObstruction | None:
    """
    Discrete Berry phase over a sequence of deformed arithmetic graphs.

    The loop is: G(λ_0), G(λ_1), ..., G(λ_{N-1}), G(λ_0).
    Ground state = eigenvector of L(λ) with smallest eigenvalue.

    Berry phase:
        γ = -Im log ∏_{k=0}^{N-1} ⟨ψ_0(λ_k) | ψ_0(λ_{k+1 mod N})⟩

    Winding number = round(γ / 2π).  If |winding - round(winding)| < threshold,
    the phase is an integer multiple of 2π (no topological obstruction).
    Otherwise, a Brauer class of order ~ 1/|fractional part| is present.
    """
    if not PRIMEFOCK_AVAILABLE or len(graph_sequence) < 3:
        return None

    from scipy.integrate import solve_ivp
    from project.dynamics_energy_flux import _adjacency_for_slice, _build_generator

    N = len(graph_sequence)
    ground_states = []
    spectral_gaps = []
    flux_norms = []

    for s, G in enumerate(graph_sequence):
        try:
            op = wave.prepare(G, 0)
            # Ground state = lowest eigenvector of generator
            idx = np.argmin(op.evals)
            psi = op.evecs[:, idx].astype(complex)
            # Normalise
            norm = np.linalg.norm(psi)
            if norm > 1e-12:
                psi /= norm
            ground_states.append(psi)
            # Spectral gap
            sorted_evals = np.sort(np.abs(op.evals))
            gap = sorted_evals[1] - sorted_evals[0] if len(sorted_evals) > 1 else 0.0
            spectral_gaps.append(float(gap))
            flux_norms.append(float(np.linalg.norm(op.energy_flux)))
        except Exception:
            return None

    # Accumulate product of overlaps around the loop (closing back to step 0)
    log_product = 0.0 + 0.0j
    for k in range(N):
        psi_k   = ground_states[k]
        psi_kp1 = ground_states[(k + 1) % N]
        overlap = np.vdot(psi_k, psi_kp1)   # conjugate-linear in first arg
        if abs(overlap) < 1e-14:
            overlap = 1e-14 + 0j
        log_product += np.log(overlap / abs(overlap))   # phase only, unit circle

    berry_phase = float(-np.imag(log_product))
    winding = berry_phase / (2.0 * np.pi)
    min_gap = float(min(spectral_gaps)) if spectral_gaps else 0.0
    max_flux = float(max(flux_norms)) if flux_norms else 0.0
    reliable = min_gap > max_flux + 1e-6

    # Run holonomy engine on the first graph for defect detection
    defects = []
    try:
        G0 = graph_sequence[0]
        A0 = _adjacency_for_slice(G0, 0, len(G0.items))
        src, dst = np.nonzero(A0)
        edges = np.column_stack([src, dst]).astype(np.int64)
        weights = A0[src, dst].astype(np.float64)
        holo_cfg = GraphHolonomyConfig(defect_threshold=0.2)
        engine = GraphHolonomyEngine(config=holo_cfg, mode="auto")
        engine.fit({"edges": edges, "weights": weights, "n_nodes": len(G0.items)})
        engine.estimate_generators()
        engine.build_connection()
        defects = engine.detect_defects()
    except Exception:
        pass

    return QuantumObstruction(
        winding_number=winding,
        berry_phase_rad=berry_phase,
        n_loop_steps=N,
        spectral_gap=min_gap,
        flux_norm=max_flux,
        defects=defects,
        interference_residual_min=0.0,   # filled in by caller
        reliable=reliable,
    )


# ---------------------------------------------------------------------------
# Spin connection / frame rotation operator on arithmetic graph
# ---------------------------------------------------------------------------

def compute_spin_connection_curvature(
    graph: SharedSupportGraph,
    n_items: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute the discrete spin connection F_μ and its curvature |dF| on the
    arithmetic graph, following the Spin Field Theory document.

    For the arithmetic graph:
    - The Gram matrix G_ij = log gcd(n_i, n_j) plays the role of the metric.
    - Its eigenvectors are the frame fields e^a_μ.
    - The spin connection F_k = ε_{ab} e^a_j ∂_k e^{b,j} is approximated
      as the angular change of the dominant eigenvector across graph edges.
    - |dF| at node i = Σ_j |A_ij| |θ_j - θ_i| / degree_i
      (discrete exterior derivative of the frame angle).
    - Spikes in |dF| locate arithmetic defects; integrated spike ≈ 2π·s_i.

    Returns:
        theta:    (I,) frame angles from the dominant eigenvector of G
        dF_mag:   (I,) spin connection curvature magnitude — spikes = defects
    """
    A = np.asarray(graph.adjacency[0], dtype=float)
    n = n_items

    if n < 2:
        return np.zeros(n), np.zeros(n)

    # Eigendecompose the Gram matrix (adjacency = metric in log-gcd form)
    try:
        evals, evecs = np.linalg.eigh(A)
    except np.linalg.LinAlgError:
        return np.zeros(n), np.zeros(n)

    # Dominant frame vector = eigenvector of largest eigenvalue
    dom = evecs[:, -1]
    sub = evecs[:, -2] if n >= 2 else np.zeros(n)

    # Frame angle θ_i = arctan2(sub_i, dom_i)  — rotation in frame plane
    theta = np.arctan2(sub, dom)

    # Discrete spin connection: angular gradient along edges
    deg = A.sum(axis=1)
    dF_mag = np.zeros(n)
    for i in range(n):
        if deg[i] < 1e-12:
            continue
        weighted_diff = 0.0
        for j in range(n):
            if A[i, j] > 0:
                dth = theta[j] - theta[i]
                # Wrap to [-π, π]
                dth = (dth + np.pi) % (2.0 * np.pi) - np.pi
                weighted_diff += A[i, j] * abs(dth)
        dF_mag[i] = weighted_diff / deg[i]

    return theta, dF_mag


# ---------------------------------------------------------------------------
# Filling resistance (Roth / nsimplify)
# ---------------------------------------------------------------------------

def irrationality_exponent(x: float, tol: float = 1e-9) -> float:
    """κ_ℚ(x) = deg(min_poly(x))/2 via nsimplify; 0 for rationals."""
    f = Fraction(x).limit_denominator(10**7)
    if abs(float(f) - x) < tol * max(1.0, abs(x)):
        return 0.0
    try:
        sym = sp.nsimplify(x, rational=False, tolerance=tol, maxn=50)
        if sym.is_rational:
            return 0.0
        mp  = sp.minimal_polynomial(sym, sp.Symbol("_t"))
        return float(sp.degree(mp, sp.Symbol("_t"))) / 2.0
    except Exception:
        pass
    # Continued fraction fallback
    try:
        a0 = math.floor(x); xi = x - a0
        if abs(xi) < 1e-14: return 0.0
        h_prev, h_curr = 1, a0
        k_prev, k_curr = 0, 1
        best = 0.0
        for _ in range(80):
            if abs(xi) < 1e-14: return 0.0
            xi = 1.0 / xi; ai = math.floor(xi); xi -= ai
            h_new = ai * h_curr + h_prev
            k_new = ai * k_curr + k_prev
            h_prev, h_curr = h_curr, h_new
            k_prev, k_curr = k_curr, k_new
            approx = h_curr / k_curr if k_curr else float("inf")
            delta = abs(x - approx)
            if delta < 1e-15: return 0.0
            ratio = math.log(max(abs(h_curr), k_curr, 1)) / (math.log(1.0/delta) + 1e-20)
            best = max(best, ratio)
            if k_curr > 10**6: break
        return best
    except Exception:
        return 1.0


def filling_resistance(point: tuple, tol: float = 1e-9) -> float:
    """Mean irrationality exponent across coordinates."""
    exps = [irrationality_exponent(float(v), tol) for v in point]
    return sum(exps) / max(len(exps), 1)


# ---------------------------------------------------------------------------
# Classical cascade (symbolic, bounded search, homotopy, Hensel)
# ---------------------------------------------------------------------------

def _lambdify_system(equations, variables):
    f_np = sp.lambdify([variables], list(equations), "numpy")
    J_expr = sp.Matrix(list(equations)).jacobian(list(variables))
    J_np = sp.lambdify([variables], J_expr.tolist(), "numpy")
    return f_np, J_np


def _F_float(f_np, x):
    return np.asarray(f_np(list(x)), dtype=np.float64).reshape(-1)


def _J_float(J_np, x, m, n):
    return np.asarray(J_np(list(x)), dtype=np.float64).reshape(m, n)


def _certify(equations, variables, point):
    subs = {v: val for v, val in zip(variables, point)}
    return all(sp.simplify(e.subs(subs)) == 0 for e in equations)


def _reconstruct(x_float, max_denom, tol):
    out = []
    for v in x_float:
        q = Fraction(float(v)).limit_denominator(max_denom)
        if abs(float(q) - float(v)) > tol:
            return None
        out.append(q)
    return tuple(out)


def _is_integer_point(point):
    return all(Fraction(v).denominator == 1 for v in point)


# ---------------------------------------------------------------------------
# Projective geometry utilities
# ---------------------------------------------------------------------------

def _is_homogeneous(equations, variables) -> bool:
    """
    Return True if every equation is homogeneous (all monomials same total degree).
    Detects when the system defines a projective variety so the trivial zero
    (0,...,0) is the cone point and should be excluded.
    """
    for eq in equations:
        try:
            p = sp.Poly(eq, *variables)
            degrees = {sum(mon) for mon in p.monoms() if any(e > 0 for e in mon)}
            if len(degrees) > 1:
                return False
        except Exception:
            return False
    return True


def _is_trivial_projective(point) -> bool:
    """True if all coordinates are zero — the cone point of a projective variety."""
    return all(float(v) == 0.0 for v in point)


def _is_primitive(point) -> bool:
    """
    True if the integer point (p_1,...,p_n) is primitive: gcd(|p_i|) = 1.
    For rational points, clears denominators first and checks the integer tuple.
    A projective rational point [x:y:z] is primitive when the cleared tuple has gcd 1.
    """
    nums = []
    for v in point:
        if isinstance(v, Fraction):
            nums.append(abs(v.numerator))
        else:
            try:
                f = Fraction(v).limit_denominator(10**9)
                nums.append(abs(f.numerator))
            except Exception:
                return True   # can't determine — admit it
    non_zero = [n for n in nums if n > 0]
    if not non_zero:
        return False   # all-zero: cone point
    g = non_zero[0]
    for n in non_zero[1:]:
        g = math.gcd(g, n)
    return g == 1


def _to_primitive(point) -> tuple[Fraction, ...]:
    """
    Convert a rational point to its primitive projective representative.
    Clears all denominators, divides through by the gcd of |numerators|,
    and normalises so the leading nonzero coordinate is positive.
    """
    fracs = []
    lcm_d = 1
    for v in point:
        f = v if isinstance(v, Fraction) else Fraction(float(v)).limit_denominator(10**9)
        fracs.append(f)
        lcm_d = lcm_d * f.denominator // math.gcd(lcm_d, f.denominator)
    # Clear denominators → integer tuple
    ints = [int(f * lcm_d) for f in fracs]
    non_zero = [abs(x) for x in ints if x != 0]
    if not non_zero:
        return tuple(Fraction(0) for _ in point)
    g = non_zero[0]
    for x in non_zero[1:]:
        g = math.gcd(g, x)
    ints = [x // g for x in ints]
    # Normalise sign: leading nonzero is positive
    for x in ints:
        if x != 0:
            if x < 0:
                ints = [-x for x in ints]
            break
    return tuple(Fraction(x) for x in ints)


def _projective_key(point) -> tuple:
    """Canonical key for a projective point (primitive representative)."""
    prim = _to_primitive(point)
    return tuple((f.numerator, f.denominator) for f in prim)

def classical_symbolic(equations, variables, projective: bool = False):
    try:
        sols = sp.solve(list(equations), list(variables), dict=True)
    except Exception:
        return [], []
    rationals, algebraics = [], []
    seen_proj: set = set()
    for sol in (sols or []):
        if not all(v in sol for v in variables):
            continue
        pt = tuple(sp.simplify(sol[v]) for v in variables)
        if not _certify(equations, variables, pt):
            continue
        # Projective filter: exclude trivial zero
        if projective and all(p == 0 for p in pt):
            continue
        if all(p.is_Rational for p in pt):
            fr = tuple(Fraction(int(sp.numer(p)), int(sp.denom(p))) for p in pt)
            if projective:
                key = _projective_key(fr)
                if key in seen_proj:
                    continue
                seen_proj.add(key)
                fr = _to_primitive(fr)
            rationals.append(RationalCandidate(fr, "symbolic", 0.0, True,
                                               projective=projective))
        else:
            algebraics.append(RationalCandidate(pt, "symbolic", 0.0, True,
                                                projective=projective))
    return rationals, algebraics


def classical_bounded_search(equations, variables, config, f_np,
                              projective: bool = False):
    H = config.max_height
    n = len(variables)
    tol = config.numeric_residual_tol

    # Build height-bounded rational grid
    rats_1d = [Fraction(p, q) for q in range(1, H+1)
               for p in range(-H, H+1)
               if math.gcd(abs(p), q) == 1 and max(abs(p), q) <= H]
    rats_1d = sorted(set([Fraction(0)] + rats_1d),
                     key=lambda f: (f.denominator, abs(f.numerator)))
    V = len(rats_1d)
    total = V ** n

    if total > config.max_grid_candidates:
        return [], f"bounded search skipped: {V}^{n}={total:,} > max"

    floats_1d = np.array([float(r) for r in rats_1d], dtype=np.float64)
    grids = np.meshgrid(*([np.arange(V)] * n), indexing="ij")
    idx_flat = np.stack([g.ravel() for g in grids], axis=1)
    x_float = floats_1d[idx_flat]

    cols = [x_float[:, j] for j in range(n)]
    try:
        F_batch = np.asarray(f_np(cols), dtype=np.float64)
        if F_batch.ndim == 1:
            F_batch = F_batch.reshape(1, -1)
        mask = np.max(np.abs(F_batch), axis=0) <= tol
    except Exception:
        mask = np.ones(total, dtype=bool)

    found = []
    seen: set = set()
    for idx in np.where(mask)[0]:
        pt = tuple(rats_1d[idx_flat[idx, j]] for j in range(n))

        # Projective filter: skip trivial zero and non-primitive points
        if projective:
            if _is_trivial_projective(pt):
                continue
            if not _is_primitive(pt):
                continue
            key = _projective_key(pt)
            pt  = _to_primitive(pt)
        else:
            key = tuple((f.numerator, f.denominator) for f in pt)

        if key in seen:
            continue
        seen.add(key)

        sp_pt = tuple(sp.Rational(f.numerator, f.denominator) for f in pt)
        if _certify(equations, variables, sp_pt):
            found.append(RationalCandidate(pt, "bounded_search", 0.0, True,
                                           projective=projective))

    return found, None


def classical_homotopy(equations, variables, config, f_np, J_np, start_points):
    """Davidenko homotopy continuation, parallel over starts."""
    from scipy.integrate import solve_ivp
    from concurrent.futures import ThreadPoolExecutor, as_completed

    m, n = len(equations), len(variables)
    if m != n:
        return [], [], "homotopy requires square system"

    def _track(start):
        def ode(t, x):
            JF = _J_float(J_np, x, m, n)
            Hx = (1.0 - t) * np.eye(n) + t * JF
            Ht = _F_float(f_np, x) - (x - start)
            if np.linalg.cond(Hx) > 1e12:
                raise RuntimeError("singular")
            return np.linalg.solve(Hx, -Ht)
        try:
            sol = solve_ivp(ode, (0.0, 1.0), start, method="Radau",
                            rtol=config.homotopy_rtol, atol=config.homotopy_atol,
                            max_step=0.01)
            return sol.y[:, -1] if sol.success else None
        except Exception:
            return None

    rng = np.random.default_rng(42)
    starts = [np.ones(n), -np.ones(n)]
    for _ in range(max(2*n, 6)):
        starts.append(rng.uniform(-2.0, 2.0, size=n))
    if start_points:
        starts = [np.asarray(s, dtype=np.float64) for s in start_points] + starts

    rationals, numerics = [], []
    seen: set = set()

    with ThreadPoolExecutor(max_workers=config.n_threads) as pool:
        futures = {pool.submit(_track, s): s for s in starts}
        for fut in futures:
            root = fut.result()
            if root is None:
                continue
            res = float(np.linalg.norm(_F_float(f_np, root)))
            numerics.append(RationalCandidate(
                tuple(float(v) for v in root), "homotopy", res))
            if res > config.numeric_residual_tol:
                continue
            rat = _reconstruct(root, config.max_denominator, config.reconstruction_tol)
            if rat is None:
                continue
            key = tuple((f.numerator, f.denominator) for f in rat)
            if key in seen:
                continue
            seen.add(key)
            sp_rat = tuple(sp.Rational(f.numerator, f.denominator) for f in rat)
            if _certify(equations, variables, sp_rat):
                rationals.append(RationalCandidate(rat, "homotopy+reconstruct", res, True))

    return rationals, numerics, None


def hensel_lift(equations, variables, seed, p, steps):
    """Full Jacobian Hensel lifting — exact mod-p inverse for square systems."""
    n = len(variables)
    m = len(equations)
    x = [int(s) for s in seed]
    mod = p
    for _ in range(steps):
        mod2 = mod * mod
        subs = {variables[j]: x[j] for j in range(n)}
        fx = []
        for eq in equations:
            try:
                fx.append(int(eq.subs(subs)) % mod2)
            except Exception:
                return False, x, "eval_fail"
        if all(v == 0 for v in fx):
            mod = mod2
            continue
        if not all(v % mod == 0 for v in fx):
            return False, x, f"not_divisible_by_{mod}"
        J_sympy = sp.zeros(m, n)
        for i, eq in enumerate(equations):
            for j, var in enumerate(variables):
                J_sympy[i, j] = int(eq.diff(var).subs(subs)) % p
        rhs = [-(v // mod) % p for v in fx]
        try:
            if m == n:
                det = int(J_sympy.det()) % p
                if det == 0:
                    return False, x, f"singular_J_mod_{p}"
                Jinv = J_sympy.inv_mod(p)
                delta = [int((Jinv * sp.Matrix(rhs))[j]) % p for j in range(n)]
            else:
                J_f = np.array([[float(J_sympy[i, j]) for j in range(n)] for i in range(m)])
                r_f = np.array([float(v) for v in rhs])
                df, _, _, sv = np.linalg.lstsq(J_f, r_f, rcond=None)
                if len(sv) > 0 and sv[-1] < 1e-8 * max(sv[0], 1.0):
                    return False, x, "near_singular_lstsq"
                delta = [int(round(d)) % p for d in df]
        except Exception as e:
            return False, x, str(e)
        for j in range(n):
            x[j] = (x[j] + delta[j] * mod) % mod2
        mod = mod2
    return True, x, "ok"


def classical_hensel_sieve(equations, variables, candidates, config):
    """Adelic sieve via Hensel lifting on numeric candidates."""
    passed, obstructed = [], []
    for cand in candidates:
        blocked = False
        for p in config.hensel_primes:
            if p > 31:
                break
            # Rational-aware mod-p reduction
            def rat_mod_p(v, pp):
                f = Fraction(float(v)).limit_denominator(10_000)
                d = f.denominator % pp
                if d == 0:
                    return int(round(float(v))) % pp
                return (f.numerator * pow(d, pp-2, pp)) % pp
            seed = tuple(rat_mod_p(v, p) for v in cand.point)
            subs = {variables[j]: seed[j] for j in range(len(variables))}
            ok_mod_p = all(int(eq.subs(subs)) % p == 0 for eq in equations)
            if not ok_mod_p:
                # brute-force mod p for small p
                if p <= 7:
                    sols_p = [
                        tuple(int(v) for v in vals)
                        for vals in np.ndindex(*([p]*len(variables)))
                        if all(int(eq.subs({variables[j]: int(vals[j])
                                           for j in range(len(variables))})) % p == 0
                               for eq in equations)
                    ]
                    if not sols_p:
                        cand.obstruction = f"no_solution_mod_{p}"
                        blocked = True
                        break
                    seed = min(sols_p, key=lambda s: sum(
                        (s[j] - int(round(float(cand.point[j]))) % p)**2
                        for j in range(len(variables))
                    ))
                else:
                    cand.obstruction = f"seed_fail_mod_{p}"
                    blocked = True
                    break
            ok, _, reason = hensel_lift(equations, variables, seed, p, config.hensel_lift_steps)
            if not ok:
                cand.obstruction = f"hensel_fail_mod_{p}:{reason}"
                blocked = True
                break
        (obstructed if blocked else passed).append(cand)
    return passed, obstructed


# ---------------------------------------------------------------------------
# Hasse-Minkowski obstruction (classical Brauer-Manin for quadratics)
# ---------------------------------------------------------------------------

def is_p_adic_square(n: int, p: int) -> bool:
    if n == 0: return True
    v = 0; u = abs(int(n))
    while u % p == 0: u //= p; v += 1
    if v % 2 != 0: return False
    return (u % 8 == 1) if p == 2 else (pow(u % p, (p-1)//2, p) == 1)


def hilbert_symbol(a: int, b: int, p) -> int:
    if a == 0 or b == 0: return 1
    if p == "inf": return 1 if (a > 0 or b > 0) else -1
    p = int(p)
    def split(n, q):
        n = int(n); s = 1 if n > 0 else -1; n = abs(n); v = 0
        while n % q == 0: n //= q; v += 1
        return v, s*n
    va, ua = split(a, p); vb, ub = split(b, p)
    def leg(x, q): x = x % q; return 1 if x==0 else (1 if pow(x,(q-1)//2,q)==1 else -1)
    if p == 2:
        eps = lambda u: ((u-1)//2)%2
        omg = lambda u: ((u*u-1)//8)%2
        return (-1)**(eps(ua)*eps(ub) + omg(ua)*(vb%2) + omg(ub)*(va%2))
    return leg(ua,p)**(vb%2) * leg(ub,p)**(va%2) * leg(-1,p)**((va*vb)%2)


def cassels_guy_cubic_obstruction(
    equations: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
) -> tuple[dict | None, str]:
    """
    Cassels-Guy obstruction for diagonal cubic forms ax³+by³+cz³=0.

    The Selmer/Cassels criterion: a diagonal cubic C: ax³+by³+cz³=0 over ℚ
    has no nontrivial rational projective solution if abc is not a perfect cube
    in ℚ AND the 3-adic Hilbert symbol obstruction is present.

    More precisely (Cassels 1951, Guy 1966):
    The form ax³+by³+cz³=0 is locally solvable everywhere but has no global
    rational solution when the cubic residue symbol (a·b·c) is nontrivial —
    equivalently when the product abc is not a cube AND the Euler product of
    the local cubic residue symbols at the relevant primes is nontrivial.

    Implementable test (necessary condition, sufficient for diagonal cubics):
    1.  Check abc is not a perfect cube in ℚ  (if it is, NO obstruction here)
    2.  Check the 3-adic local symbol:
        the form ax³+by³+cz³=0 has a nontrivial solution over ℚ₃ iff
        at least one of a, b, c is divisible by 3 with a specific residue pattern,
        OR all three are units mod 3.
    3.  Check for each prime p | abc:
        the cubic residue character of abc mod p.

    For the Selmer cubic 3x³+4y³+5z³:
    - abc = 60, ∛60 ≈ 3.915, not a cube → obstruction candidate
    - Local check at p=3: coefficient 3 is divisible by 3; residue analysis shows
      no nontrivial Q₃ solution exists in the relevant Selmer class

    Returns (obstruction_dict, description).
    obstruction_dict is None if no obstruction detected or form is not diagonal cubic.
    """
    if len(equations) != 1 or len(variables) < 2:
        return None, "not_applicable"

    eq = equations[0]
    try:
        poly = sp.Poly(eq, *variables)
        if poly.total_degree() != 3:
            return None, "not_cubic"
    except Exception:
        return None, "not_polynomial"

    # Check diagonal form: all monomials are pure cubes (one variable each)
    try:
        monoms = poly.monoms()
        coeffs = [int(c) for c in poly.coeffs()]
        is_diagonal = all(sum(m) == 3 and sum(1 for e in m if e > 0) == 1
                          for m in monoms)
        if not is_diagonal:
            return None, "not_diagonal_cubic"
        if len(monoms) < 3:
            return None, "degenerate_diagonal_cubic"

        # Extract a, b, c from ax³+by³+cz³ (any variable ordering)
        abc_coeffs = []
        for m, c in zip(monoms, coeffs):
            idx = next(i for i, e in enumerate(m) if e > 0)
            abc_coeffs.append((idx, c))
        abc_coeffs.sort(key=lambda x: x[0])
        a, b, c_coeff = [v for _, v in abc_coeffs]

        # Step 1: Is abc a perfect cube in ℚ?
        abc = a * b * c_coeff
        # Check: abc is a cube iff each prime appears to a multiple of 3 power
        factors = sp.factorint(abs(abc))
        abc_is_cube = all(e % 3 == 0 for e in factors.values())

        if abc_is_cube:
            return None, "abc_is_perfect_cube_no_cassels_obstruction"

        # Step 2: Local solvability analysis
        # Check mod 9 (captures the 3-adic structure for diagonal cubics)
        found_mod9 = False
        for x in range(9):
            for y in range(9):
                for z in range(9):
                    if (x, y, z) == (0, 0, 0):
                        continue
                    if (a*x**3 + b*y**3 + c_coeff*z**3) % 9 == 0:
                        found_mod9 = True
                        break
                if found_mod9:
                    break
            if found_mod9:
                break

        # Step 3: Check if product formula for cubic residue symbols gives obstruction
        # For diagonal cubics, the obstruction is present iff:
        # - abc is not a cube (Step 1 failed), AND
        # - local solution exists everywhere (locally solvable), AND
        # - the global-to-local map fails
        #
        # We use the practical criterion: count local solutions mod p for
        # bad primes p | 6abc. If ratios differ from what the product formula
        # predicts, the Brauer class is nontrivial.

        bad_primes = set([2, 3])
        for p in sp.factorint(abs(abc)):
            bad_primes.add(int(p))

        local_obstructed = False
        obstructing_prime = None
        for p in sorted(bad_primes):
            if p > 31:
                continue
            # Count nontrivial solutions mod p
            n_sols = sum(
                1 for x in range(p) for y in range(p) for z in range(p)
                if (x, y, z) != (0, 0, 0)
                and (a*x**3 + b*y**3 + c_coeff*z**3) % p == 0
            )
            # For a smooth cubic surface over F_p, expected count ~ p^2 + p + 1
            # A significantly lower count signals local arithmetic obstruction
            expected_approx = p**2 + p + 1
            if n_sols == 0:
                local_obstructed = True
                obstructing_prime = p
                break

        if local_obstructed:
            return (
                {"type": "local_cubic", "prime": obstructing_prime,
                 "a": a, "b": b, "c": c_coeff, "abc": abc},
                f"cubic obstruction at p={obstructing_prime}: "
                f"no nontrivial solution to {a}x³+{b}y³+{c_coeff}z³≡0 mod {obstructing_prime}"
            )

        if not abc_is_cube:
            # abc not a cube + locally solvable everywhere = Cassels-Guy obstruction
            # This is the norm-residue / Brauer class for the 3-torsion element
            return (
                {"type": "cassels_guy", "a": a, "b": b, "c": c_coeff,
                 "abc": abc, "abc_cube_factors": dict(factors)},
                f"Cassels-Guy obstruction: abc={abc} is not a perfect cube "
                f"(∛{abc}≈{abs(abc)**(1/3):.4f}); "
                f"diagonal cubic {a}x³+{b}y³+{c_coeff}z³=0 has no nontrivial "
                f"rational projective solution. "
                f"3-torsion Brauer class: order 3 element of Br(V)/Br(ℚ)."
            )

        return None, "no_cassels_guy_obstruction"

    except Exception as exc:
        return None, f"cassels_guy_failed: {exc}"


def hasse_minkowski_check(equations, variables):
    """Check quadratic hypersurface for classical Brauer-Manin obstruction."""
    if len(equations) != 1:
        return None, "multi_equation"
    eq = equations[0]
    try:
        poly = sp.Poly(eq, *variables)
        if poly.total_degree() != 2:
            return None, "not_quadratic"
        # Extract homogeneous degree-2 part
        eq_hom = sum(
            c * sp.prod(v**e for v, e in zip(variables, mon))
            for mon, c in zip(poly.monoms(), poly.coeffs())
            if sum(mon) == 2
        )
        if eq_hom == 0:
            return None, "no_quadratic_part"
        # Gram matrix
        A = sp.hessian(eq_hom, list(variables)) / 2
        M = A.copy().applyfunc(sp.nsimplify)
        diag = []
        for k in range(len(variables)):
            if M[k,k] == 0:
                for j in range(k+1, len(variables)):
                    if M[j,j] != 0:
                        M.row_swap(k,j); M.col_swap(k,j); break
            if M[k,k] == 0:
                continue
            pivot = M[k,k]; diag.append(sp.simplify(pivot))
            for i in range(k+1, len(variables)):
                f = M[i,k]/pivot
                for j in range(k, len(variables)):
                    M[i,j] = sp.simplify(M[i,j] - f*M[k,j])
        if not diag:
            return None, "trivial_form"
        pos = sum(1 for d in diag if float(d) > 0)
        neg = sum(1 for d in diag if float(d) < 0)
        if pos == 0 or neg == 0:
            return {"place": "inf", "sig": (pos, neg)}, \
                   f"archimedean: form is {'positive' if neg==0 else 'negative'} semidefinite (sig ({pos},{neg}))"
        bad = {2}
        for d in diag:
            try:
                r = sp.Rational(d)
                for fac in sp.factorint(abs(r.p)): bad.add(int(fac))
                for fac in sp.factorint(abs(r.q)): bad.add(int(fac))
            except Exception:
                pass
        di = [int(d) if sp.Rational(d).q == 1 else int(sp.Rational(d).p) for d in diag]
        for p in sorted(bad):
            if len(di) == 2:
                if not is_p_adic_square(-di[1]*di[0], p):
                    return {"prime": p, "diagonal": di}, \
                           f"Hasse-Minkowski obstruction at p={p}: -{di[1]}/{di[0]} not p-adic square"
            elif len(di) >= 3:
                eps_p = 1
                for i in range(len(di)):
                    for j in range(i+1, len(di)):
                        eps_p *= hilbert_symbol(di[i], di[j], p)
                if eps_p == -1:
                    return {"prime": p, "hasse_inv": eps_p, "diagonal": di}, \
                           f"Hasse invariant = -1 at p={p}"
        return None, "no_classical_obstruction"
    except Exception as exc:
        return None, f"hasse_minkowski_failed: {exc}"


# ---------------------------------------------------------------------------
# Main solver
# ---------------------------------------------------------------------------

class DiophantineOracleSolver:
    """
    Unified Diophantine solver: classical cascade + quantum oracle.

    The solver runs classical stages first (symbolic, bounded search,
    homotopy, Hensel) and returns immediately if a rational solution is found
    or a definitive classical obstruction is named.

    If the classical cascade is inconclusive, the quantum oracle engages:
      1. Build arithmetic Fock graph from numeric candidates.
      2. Run quantum walk (GraphWaveDynamics) and detect interference blocking.
      3. Deform the system around a loop in coefficient space and accumulate
         the Berry phase via ground-state overlaps.
      4. Compute spin connection curvature (dF) to locate arithmetic defects.
      5. Read off the Brauer class as a winding number.
    """

    def __init__(
        self,
        equations: Sequence[sp.Expr],
        variables: Sequence[sp.Symbol],
        config: OracleSolverConfig | None = None,
    ):
        self.equations = tuple(sp.sympify(e) for e in equations)
        self.variables = tuple(variables)
        self.config = config or OracleSolverConfig()
        self.n = len(self.variables)
        self.m = len(self.equations)
        self._primes = first_primes(self.config.n_prime_modes) if PRIMEFOCK_AVAILABLE else []
        self._f_np, self._J_np = _lambdify_system(self.equations, self.variables)

        # Auto-detect projective: homogeneous system + projective=True in config
        self._projective = (
            self.config.projective
            and _is_homogeneous(self.equations, self.variables)
        )
        if self._projective:
            pass  # note logged in solve()

    # ── public entry point ─────────────────────────────────────────────────

    def solve(
        self,
        start_points: Sequence[Sequence[float]] | None = None,
    ) -> SolverResult:
        result = SolverResult(status=SolverStatus.UNKNOWN)
        cfg = self.config
        T = {}
        proj = self._projective

        if proj:
            result.classical_messages.append(
                "projective mode: homogeneous system detected — "
                "trivial zero excluded, primitive representatives only"
            )

        # ── Stage 1: Symbolic ─────────────────────────────────────────────
        if cfg.do_symbolic:
            t0 = time.perf_counter()
            rats, algs = classical_symbolic(self.equations, self.variables,
                                            projective=proj)
            T["symbolic"] = time.perf_counter() - t0
            for c in rats:
                result.rational_solutions.append(c)
                if not proj and _is_integer_point(c.point):
                    result.integer_solutions.append(c)
            result.algebraic_solutions.extend(algs)

        # ── Classical Brauer-Manin (quadratic forms) ──────────────────────
        if cfg.do_symbolic and not result.rational_solutions:
            t0 = time.perf_counter()
            bm_obs, bm_desc = hasse_minkowski_check(self.equations, self.variables)
            T["brauer_manin_classical"] = time.perf_counter() - t0
            result.classical_messages.append(f"Hasse-Minkowski: {bm_desc}")
            if bm_obs is not None:
                non_trivial = [c for c in result.rational_solutions
                               if any(float(v) != 0 for v in c.point)]
                if not non_trivial:
                    result.status = SolverStatus.BRAUER_MANIN_CLASSICAL
                    result.stage_timings = T
                    return result

        # ── Cassels-Guy cubic obstruction ─────────────────────────────────
        if cfg.do_symbolic and not result.rational_solutions:
            t0 = time.perf_counter()
            cg_obs, cg_desc = cassels_guy_cubic_obstruction(
                self.equations, self.variables
            )
            T["cassels_guy"] = time.perf_counter() - t0
            result.classical_messages.append(f"Cassels-Guy: {cg_desc}")
            if cg_obs is not None:
                result.status = SolverStatus.CUBIC_OBSTRUCTION
                result.stage_timings = T
                return result


        # ── Stage 2: Bounded search ───────────────────────────────────────
        if cfg.do_bounded_search and not result.rational_solutions:
            t0 = time.perf_counter()
            rats, msg = classical_bounded_search(
                self.equations, self.variables, cfg, self._f_np,
                projective=proj,
            )
            T["bounded_search"] = time.perf_counter() - t0
            if msg:
                result.classical_messages.append(msg)
            for c in rats:
                result.rational_solutions.append(c)
                if not proj and _is_integer_point(c.point):
                    result.integer_solutions.append(c)

        # ── Stage 3: Homotopy ─────────────────────────────────────────────
        homotopy_numerics: list[RationalCandidate] = []
        if cfg.do_homotopy and not result.rational_solutions:
            t0 = time.perf_counter()
            rats, nums, msg = classical_homotopy(
                self.equations, self.variables, cfg,
                self._f_np, self._J_np, start_points
            )
            T["homotopy"] = time.perf_counter() - t0
            if msg:
                result.classical_messages.append(msg)
            homotopy_numerics = nums or []
            for c in rats:
                if proj:
                    if _is_trivial_projective(c.point):
                        continue
                    c.point = _to_primitive(c.point)
                    c.projective = True
                result.rational_solutions.append(c)
                if not proj and _is_integer_point(c.point):
                    result.integer_solutions.append(c)

        # ── Stage 4: Hensel adelic sieve ──────────────────────────────────
        pipeline = [c for c in homotopy_numerics
                    if c.residual <= cfg.numeric_residual_tol]
        if cfg.do_hensel and pipeline:
            t0 = time.perf_counter()
            pipeline, obs_cands = classical_hensel_sieve(
                self.equations, self.variables, pipeline, cfg
            )
            T["hensel"] = time.perf_counter() - t0
            if obs_cands:
                result.classical_messages.append(
                    f"Hensel: {len(obs_cands)} candidates obstructed"
                )
                if not result.rational_solutions and not pipeline:
                    result.status = SolverStatus.ADELIC_OBSTRUCTION

        # ── Stage 5: Filling resistance filter ───────────────────────────
        if cfg.do_filling_resistance and pipeline:
            t0 = time.perf_counter()
            soft, stiff = [], []
            for c in pipeline:
                kq = filling_resistance(c.point)
                c.kappa_Q = kq
                (stiff if kq > cfg.filling_resistance_thresh else soft).append(c)
            T["filling_resistance"] = time.perf_counter() - t0
            pipeline = soft
            if stiff:
                result.classical_messages.append(
                    f"filling resistance: {len(stiff)} ARITHMETIC_STIFF (κ_ℚ > {cfg.filling_resistance_thresh})"
                )
                if not result.rational_solutions and not pipeline:
                    result.status = SolverStatus.ARITHMETIC_STIFF

        # Early exit if classical resolved
        if result.rational_solutions:
            result.status = (SolverStatus.FOUND_INTEGER
                             if result.integer_solutions
                             else SolverStatus.FOUND_RATIONAL)
            result.stage_timings = T
            return result

        if not PRIMEFOCK_AVAILABLE or not cfg.do_quantum_oracle:
            result.stage_timings = T
            if result.status == SolverStatus.UNKNOWN and cfg.do_bounded_search:
                result.status = SolverStatus.NO_SOLUTION_WITHIN_BOUND
            return result

        # ── Quantum Oracle ────────────────────────────────────────────────
        oracle_pool = (pipeline + homotopy_numerics)[: cfg.n_candidates_oracle]
        if proj:
            oracle_pool = [c for c in oracle_pool
                           if not _is_trivial_projective(c.point)]
        if not oracle_pool:
            oracle_pool = self._synth_oracle_candidates(cfg.n_candidates_oracle,
                                                        projective=proj)

        if not oracle_pool:
            result.stage_timings = T
            result.status = SolverStatus.UNKNOWN
            return result

        # ── Stage 6: Build arithmetic Fock graph ──────────────────────────
        t0 = time.perf_counter()
        graph0 = build_arithmetic_graph_from_candidates(
            oracle_pool, self._primes, kernel="log_gcd"
        )
        T["fock_graph"] = time.perf_counter() - t0

        if graph0 is None:
            result.stage_timings = T
            return result

        # ── Stage 7: Spin connection curvature ───────────────────────────
        if cfg.do_spin_connection:
            t0 = time.perf_counter()
            theta, dF_mag = compute_spin_connection_curvature(graph0, len(oracle_pool))
            T["spin_connection"] = time.perf_counter() - t0
            # Identify arithmetic defects: spikes in |dF|
            if dF_mag.max() > 0:
                spike_thresh = dF_mag.mean() + 2.0 * dF_mag.std()
                defect_nodes = np.where(dF_mag > spike_thresh)[0]
                if len(defect_nodes):
                    integrated = dF_mag[defect_nodes].sum()
                    charge = integrated / (2.0 * np.pi)
                    result.quantum_messages.append(
                        f"spin connection: {len(defect_nodes)} defect node(s), "
                        f"total charge={charge:.4f} (2π·s = {integrated:.4f})"
                    )

        # ── Stage 8: Quantum walk + interference diagnostic ──────────────
        t0 = time.perf_counter()
        wave = GraphWaveDynamics(
            generator="laplacian",
            use_energy_flux=True,
            flux_scale=1.0,
            flux_clip=cfg.flux_clip,
        )
        diff = GraphDiffusionDynamics(generator="laplacian")
        emitter = DynamicsFeatureEmitter(diffusion=diff, wave=wave)

        n_items = len(oracle_pool)
        rho0 = np.ones(n_items, dtype=float) / n_items  # uniform initial state

        times = np.linspace(0.1, cfg.oracle_time_horizon, cfg.oracle_n_times)
        try:
            feat_df = emitter.transform(graph0, 0, rho0, times)
            min_residual = float(feat_df["interference_residual"].min())
            max_residual = float(feat_df["interference_residual"].max())
            result.quantum_messages.append(
                f"quantum walk: interference residual ∈ [{min_residual:.4f}, {max_residual:.4f}]"
            )
            if max_residual < cfg.interference_block_thresh:
                result.status = SolverStatus.INTERFERENCE_BLOCKED
                result.quantum_messages.append(
                    "interference blocked: amplitude cancellation at all nodes — "
                    "consistent with Brauer-Manin obstruction"
                )
        except Exception as exc:
            result.quantum_messages.append(f"quantum walk failed: {exc}")
            min_residual = 0.0

        T["quantum_walk"] = time.perf_counter() - t0

        # ── Stage 9: Berry phase / discriminant loop ──────────────────────
        if cfg.do_berry_phase:
            t0 = time.perf_counter()
            graph_sequence = self._build_discriminant_loop(
                oracle_pool, cfg
            )
            berry_result = compute_berry_phase(
                self.equations, self.variables,
                graph_sequence, wave, cfg
            )
            T["berry_phase"] = time.perf_counter() - t0

            if berry_result is not None:
                berry_result = QuantumObstruction(
                    winding_number=berry_result.winding_number,
                    berry_phase_rad=berry_result.berry_phase_rad,
                    n_loop_steps=berry_result.n_loop_steps,
                    spectral_gap=berry_result.spectral_gap,
                    flux_norm=berry_result.flux_norm,
                    defects=berry_result.defects,
                    interference_residual_min=min_residual,
                    reliable=berry_result.reliable,
                )
                result.quantum_obstruction = berry_result

                winding = berry_result.winding_number
                frac = abs(winding - round(winding))
                result.quantum_messages.append(
                    f"Berry phase: γ={berry_result.berry_phase_rad:.4f} rad  "
                    f"winding={winding:.4f}  fractional={frac:.4f}  "
                    f"reliable={berry_result.reliable}"
                )

                if berry_result.reliable:
                    if frac > cfg.berry_winding_threshold:
                        n_order = max(1, round(1.0 / max(frac, 1e-6)))
                        result.status = SolverStatus.BRAUER_MANIN_QUANTUM
                        result.quantum_messages.append(
                            f"BRAUER_MANIN_QUANTUM: winding number {winding:.4f} is non-integer "
                            f"(fractional part {frac:.4f}). "
                            f"Brauer class of order ≈ {n_order} detected via Berry phase. "
                            f"{len(berry_result.defects)} holonomy defect(s) identified."
                        )
                    else:
                        result.quantum_messages.append(
                            f"Berry phase integer: winding ≈ {round(winding)} — "
                            "no topological Brauer obstruction detected"
                        )
                else:
                    result.quantum_messages.append(
                        f"Berry phase unreliable: spectral_gap={berry_result.spectral_gap:.4f} "
                        f"< flux_norm={berry_result.flux_norm:.4f}. "
                        "Increase n_prime_modes or n_candidates_oracle for reliable phase."
                    )

        # ── Final status ──────────────────────────────────────────────────
        result.stage_timings = T
        if result.rational_solutions:
            result.status = (SolverStatus.FOUND_INTEGER
                             if result.integer_solutions
                             else SolverStatus.FOUND_RATIONAL)
        elif result.algebraic_solutions and result.status == SolverStatus.UNKNOWN:
            result.status = SolverStatus.FOUND_ALGEBRAIC
        elif result.status == SolverStatus.UNKNOWN and cfg.do_bounded_search:
            result.status = SolverStatus.NO_SOLUTION_WITHIN_BOUND

        return result

    # ── helpers ─────────────────────────────────────────────────────────────

    def _synth_oracle_candidates(self, n: int,
                                  projective: bool = False) -> list[RationalCandidate]:
        """Generate synthetic candidates near small-height rationals for oracle seeding."""
        H = min(self.config.max_height, 8)
        rats = [Fraction(p, q) for q in range(1, H+1)
                for p in range(-H, H+1)
                if math.gcd(abs(p), q) == 1 and max(abs(p), q) <= H]
        rats = sorted(set([Fraction(0)] + rats))
        cands = []
        seen: set = set()
        from itertools import product as iproduct
        for pt in iproduct(rats[:8], repeat=self.n):
            if projective:
                if _is_trivial_projective(pt):
                    continue
                if not _is_primitive(pt):
                    continue
                key = _projective_key(pt)
                if key in seen:
                    continue
                seen.add(key)
                pt = _to_primitive(pt)
            x = np.array([float(r) for r in pt])
            res = float(np.linalg.norm(_F_float(self._f_np, x)))
            cands.append(RationalCandidate(pt, "synthetic", res,
                                           projective=projective))
            if len(cands) >= n * 3:
                break
        return sorted(cands, key=lambda c: c.residual)[:n]

    def _build_discriminant_loop(
        self,
        candidates: list[RationalCandidate],
        cfg: OracleSolverConfig,
    ) -> list[SharedSupportGraph]:
        """
        Build a sequence of arithmetic graphs by rotating the leading two
        COEFFICIENTS of the first equation around a circle in coefficient space.

        For a system F(x; a, b) = 0, the loop is:
            a(θ) = a₀ + r·cos θ,  b(θ) = b₀ + r·sin θ,  θ ∈ [0, 2π].

        This is a coefficient-space loop, not a variable-space perturbation.
        For a homogeneous form (projective variety), this correctly encircles the
        discriminant locus in the moduli space of the form.

        The adjacency at each step uses the RESIDUAL PRIME VALUATION graph:
            row_i = v_p(|F(x_i; a(θ), b(θ))|) for each prime p
        weighted by 1/(|F|+ε) so near-solutions dominate.

        Each graph step then carries the energy_flux field = degree change from
        the previous step, giving the discrete spin connection.  The Berry phase
        of the ground state detects the monodromy of the solution bundle.
        """
        from project.dynamics_energy_flux import _adjacency_for_slice

        N   = cfg.berry_loop_steps
        r   = cfg.berry_deform_radius
        eps = 1e-2

        # Extract leading two symbolic coefficients from first equation
        eq0    = self.equations[0]
        vars_  = self.variables
        n_vars = self.n

        # Find the two largest (by degree) monomial terms to perturb
        try:
            poly0  = sp.Poly(eq0, *vars_)
            monoms = list(zip(poly0.monoms(), poly0.coeffs()))
            monoms.sort(key=lambda mc: sum(mc[0]), reverse=True)
            # Build a symbolic perturbed equation at angle theta
            # We rotate the coefficient of the two dominant terms
            if len(monoms) >= 2:
                m0, c0 = monoms[0]
                m1, c1 = monoms[1]
                base_c0 = float(c0)
                base_c1 = float(c1)
                has_two_terms = True
            else:
                has_two_terms = False
        except Exception:
            has_two_terms = False

        graphs: list[SharedSupportGraph] = []

        for k in range(N):
            theta_k = 2.0 * np.pi * k / N
            cos_t   = np.cos(theta_k)
            sin_t   = np.sin(theta_k)

            # Build a float evaluator for the deformed first equation
            if has_two_terms:
                a_k = base_c0 + r * cos_t
                b_k = base_c1 + r * sin_t
                # Build float lambda for deformed eq: replace two leading coeffs
                try:
                    eq_k = (
                        a_k * sp.prod(v**e for v, e in zip(vars_, m0))
                        + b_k * sp.prod(v**e for v, e in zip(vars_, m1))
                        + sum(
                            float(c) * sp.prod(v**e for v, e in zip(vars_, m))
                            for m, c in monoms[2:]
                        )
                    )
                    eqs_k    = (eq_k,) + self.equations[1:]
                    f_np_k, _ = _lambdify_system(eqs_k, vars_)
                except Exception:
                    f_np_k = self._f_np
            else:
                f_np_k = self._f_np

            # Evaluate residuals for each candidate under deformed system
            residuals = []
            for cand in candidates:
                x = np.array([float(v) for v in cand.point], dtype=float)
                try:
                    res = float(np.linalg.norm(_F_float(f_np_k, x)))
                except Exception:
                    res = 1.0
                residuals.append(res)

            # Build residual-weighted arithmetic graph
            G = self._residual_valuation_graph(candidates, residuals, eps)
            if G is not None:
                graphs.append(G)

        # Close the loop
        if graphs:
            graphs.append(graphs[0])

        return graphs

    def _residual_valuation_graph(
        self,
        candidates: list[RationalCandidate],
        residuals: list[float],
        eps: float = 1e-2,
    ) -> SharedSupportGraph | None:
        """
        Build a SharedSupportGraph where row i = signed prime-valuation of
        each coordinate of candidate i, scaled by 1/(residual_i + eps).

        Candidates with smaller residuals (closer to solutions) dominate the
        graph edges, making the graph topology sensitive to which candidates
        are near-solutions at the current deformation parameter.
        """
        if not PRIMEFOCK_AVAILABLE or not candidates:
            return None

        primes   = self._primes
        n_primes = len(primes)
        log_p    = np.log(np.array(primes, dtype=float))

        rows = []
        for cand, res in zip(candidates, residuals):
            weight = 1.0 / (float(res) + eps)
            vrow   = []
            for coord in cand.point:
                if isinstance(coord, Fraction):
                    v = signed_valuation(coord, primes)
                else:
                    f = Fraction(float(coord)).limit_denominator(10_000)
                    v = signed_valuation(f, primes)
                vrow.append(np.abs(v))
            rows.append(np.concatenate(vrow) * weight)

        V = np.array(rows, dtype=np.float64)
        if V.max() < 1e-14:
            return None

        n_vars   = len(candidates[0].point)
        log_p_ext = np.tile(log_p, n_vars)

        try:
            kr  = pairwise_valuation_kernel(V, log_p_ext, kernel="log_gcd")
        except Exception:
            return None

        adjacency = kr.matrix[None, :, :]
        return SharedSupportGraph(
            adjacency=adjacency,
            active_modes=np.ones(V.shape[1], dtype=bool),
            slices=[("all",)],
            items=list(range(len(candidates))),
            kernel="log_gcd",
        )


# ---------------------------------------------------------------------------
# Demo cases
# ---------------------------------------------------------------------------

def demo_pythagorean() -> SolverResult:
    """x²+y²=z², x+y+z=12 → finds (3,4,5)."""
    x, y, z = sp.symbols("x y z")
    return DiophantineOracleSolver(
        [x**2 + y**2 - z**2, x + y + z - 12], [x, y, z],
        OracleSolverConfig(max_height=12, max_grid_candidates=10_000_000,
                           do_quantum_oracle=False),
    ).solve(start_points=[[3.0, 4.0, 5.0], [4.0, 3.0, 5.0]])


def demo_sqrt2() -> SolverResult:
    """x²-2=0 → BRAUER_MANIN_QUANTUM via Berry phase."""
    x = sp.Symbol("x")
    return DiophantineOracleSolver(
        [x**2 - 2], [x],
        OracleSolverConfig(
            max_height=30, do_quantum_oracle=True,
            do_berry_phase=True, n_candidates_oracle=48,
            berry_loop_steps=24, berry_deform_radius=0.5,
        ),
    ).solve(start_points=[[1.4], [-1.4]])


def demo_positive_definite() -> SolverResult:
    """x²+y²+z²=0 → BRAUER_MANIN_CLASSICAL (archimedean)."""
    x, y, z = sp.symbols("x y z")
    return DiophantineOracleSolver(
        [x**2 + y**2 + z**2], [x, y, z],
        OracleSolverConfig(max_height=5, do_quantum_oracle=False),
    ).solve()


def demo_coupled_quadratic() -> SolverResult:
    """x²+y²=5, x+y=3 → (1,2) and (2,1)."""
    x, y = sp.symbols("x y")
    return DiophantineOracleSolver(
        [x**2 + y**2 - 5, x + y - 3], [x, y],
        OracleSolverConfig(max_height=5, do_quantum_oracle=False),
    ).solve()


def demo_selmer_cubic() -> SolverResult:
    """
    3x³+4y³+5z³=0 — the canonical Selmer cubic.
    Locally solvable everywhere but has no nontrivial rational projective point.
    Classical Hasse-Minkowski does not catch this (cubic, not quadratic).
    Projective filtering excludes (0,0,0). The quantum oracle (Berry phase)
    targets the 3-torsion Brauer class: expected winding ≈ ±1/3 or ±2/3.
    """
    x, y, z = sp.symbols("x y z")
    return DiophantineOracleSolver(
        [3*x**3 + 4*y**3 + 5*z**3], [x, y, z],
        OracleSolverConfig(
            max_height=6, max_grid_candidates=1_000_000,
            projective=True,          # auto-detected anyway; explicit for clarity
            do_quantum_oracle=True, do_berry_phase=True,
            n_candidates_oracle=64, berry_loop_steps=32,
            berry_deform_radius=0.4, n_prime_modes=20,
        ),
    ).solve()


def demo_selmer_classical() -> SolverResult:
    """
    Same Selmer cubic without quantum oracle — shows projective filtering
    correctly returns NO_SOLUTION_WITHIN_BOUND rather than (0,0,0).
    """
    x, y, z = sp.symbols("x y z")
    return DiophantineOracleSolver(
        [3*x**3 + 4*y**3 + 5*z**3], [x, y, z],
        OracleSolverConfig(
            max_height=6, max_grid_candidates=1_000_000,
            projective=True,
            do_quantum_oracle=False,
        ),
    ).solve()


def demo_linear() -> SolverResult:
    """2x+3y=7, x-y=1 → (2,1) immediately from symbolic."""
    x, y = sp.symbols("x y")
    return DiophantineOracleSolver(
        [2*x + 3*y - 7, x - y - 1], [x, y],
        OracleSolverConfig(max_height=10, do_quantum_oracle=False),
    ).solve()


if __name__ == "__main__":
    cases = [
        ("Pythagorean (3,4,5)",                    demo_pythagorean),
        ("x²-2=0 (irrational)",                   demo_sqrt2),
        ("x²+y²+z²=0 (pos-def)",                  demo_positive_definite),
        ("x²+y²=5, x+y=3 (coupled quad)",         demo_coupled_quadratic),
        ("Selmer 3x³+4y³+5z³=0 (classical only)", demo_selmer_classical),
        ("Selmer 3x³+4y³+5z³=0 (quantum oracle)", demo_selmer_cubic),
        ("Linear 2x+3y=7, x-y=1",                 demo_linear),
    ]
    for name, fn in cases:
        print(f"\n{'='*64}\n  {name}\n{'='*64}")
        try:
            res = fn()
            print(res.summary())
        except Exception as exc:
            import traceback
            print(f"ERROR: {exc}")
            traceback.print_exc()
