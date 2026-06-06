"""
dynamics.py — Density and amplitude dynamics on the arithmetic relatedness graph
================================================================================

Turns the signal-extraction layer into a *model*: given a relatedness graph
(dense SharedSupportGraph or sparse EdgeListGraph), evolve a state forward in
time on the graph and read out predictions / features.

Three layers, following the pleroma/CFG symplectic-lift picture:

  STAGE 1 — GraphDiffusionDynamics  (real, dissipative, the honest default)
      dρ/dt = -L ρ        (graph heat flow; L = D - A is the Laplacian)
      ρ(t) = V exp(-Λ t) Vᵀ ρ(0)
      Density spreads and relaxes toward the stationary distribution.
      Irreversible, mass-conserving (for the combinatorial Laplacian),
      arrow-of-time correct for an open system. NO interference is possible:
      a nonnegative density can only spread and overlap, never cancel.

  STAGE 2 — GraphWaveDynamics  (complex amplitude, OPTIONALLY damped)
      Unitary (gamma=0):   i dψ/dt = L ψ                (quantum walk)
      Damped  (gamma>0):   i dψ/dt = (L - i·gamma·D) ψ   (dissipative)
      ψ = √ρ·e^{iS} is the Madelung/symplectic lift. Complex amplitudes let
      two assets' spreading influence INTERFERE — reinforce at some nodes,
      cancel at others. This is the one genuinely new phenomenon over diffusion.
      *** This layer is a PHASE-SENSITIVE FEATURE GENERATOR, not a claim that
      market states evolve unitarily. Use the damped variant (gamma>0) for
      dissipative systems so amplitudes decay rather than ring forever. ***

  FEATURES — DynamicsFeatureEmitter
      Emits model-ready features from both layers. The headline feature is the
      INTERFERENCE RESIDUAL  |ψ(t)|² - ρ(t): same graph, same initial state,
      same horizon, one propagation with interference and one without. Whether
      that residual predicts future residual returns / vol expansion /
      correlation breaks is THE empirical test of whether the wave layer earns
      its place.

Both real and complex layers use the SAME generator L, so any predictive
difference is attributable to the lift (phase/interference), not to a different
model — exactly the A/B test the lift demands.

CAVEAT (markets): unitary evolution is reversible and norm-conserving, a poor
prior for dissipative, arrow-of-time markets. Stage 1 (diffusion) is the
defensible model; Stage 2 is exploratory and its value must be demonstrated
empirically. The damped-wave option exists precisely so the wave layer is not
forced to be lossless on a lossy system.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional

import numpy as np


GeneratorName = Literal["laplacian", "adjacency", "normalized_laplacian"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _normalize_density(rho: np.ndarray, *, clip: bool = True, eps: float = 1e-12) -> np.ndarray:
    """
    Clip tiny negative numerical artifacts and renormalize to a probability
    vector. The combinatorial-Laplacian heat kernel preserves nonnegativity in
    exact arithmetic, but eigendecomposition can introduce small negatives;
    this cleans them. Optional so raw spectral output remains available for
    debugging.
    """
    rho = np.asarray(rho, dtype=float)
    if clip:
        rho = np.maximum(rho, 0.0)
    total = rho.sum()
    if total <= eps:
        return np.ones_like(rho) / max(len(rho), 1)
    return rho / total


def _adjacency_for_slice(graph, s: int, n_items: int) -> np.ndarray:
    """
    Return the dense (I, I) symmetric adjacency for slice s, from either a
    dense SharedSupportGraph (.adjacency shape (S,I,I)) or a sparse
    EdgeListGraph (.edges DataFrame).
    """
    if hasattr(graph, "adjacency"):
        A = np.asarray(graph.adjacency, dtype=float)[s]
        return 0.5 * (A + A.T)
    if hasattr(graph, "edges"):
        A = np.zeros((n_items, n_items), dtype=float)
        part = graph.edges[graph.edges["slice_idx"] == s]
        for row in part.itertuples(index=False):
            i, j, w = int(row.source_idx), int(row.target_idx), float(row.weight)
            A[i, j] += w
            A[j, i] += w
        return A
    raise TypeError("graph must expose either .adjacency or .edges")


def _degrees(A: np.ndarray) -> np.ndarray:
    return A.sum(axis=1)


def _build_generator(A: np.ndarray, generator: GeneratorName) -> np.ndarray:
    """
    Build the real symmetric generator from the adjacency.

    laplacian             L = D - A.  Satisfies 1ᵀL = 0, so diffusion conserves
                          total mass in ordinary coordinates. This is the only
                          generator for which 'probability-conserving density
                          diffusion' is literally true.

    adjacency             H = A.  Tight-binding generator (standard for quantum
                          walks). Not mass-conserving as a diffusion generator;
                          intended for the wave layer.

    normalized_laplacian  L_sym = I - D^{-1/2} A D^{-1/2}.  GEOMETRY-NORMALIZED,
                          NOT ordinary density diffusion: 1ᵀL_sym ≠ 0 in general,
                          so exp(-L_sym t) does NOT conserve mass in ordinary
                          coordinates. Use only when you explicitly want the
                          degree-normalized geometry, and renormalize outputs.
    """
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    if generator == "adjacency":
        return A
    deg = _degrees(A)
    if generator == "laplacian":
        return np.diag(deg) - A
    if generator == "normalized_laplacian":
        d_inv_sqrt = np.where(deg > 0, 1.0 / np.sqrt(deg), 0.0)
        L = np.diag(deg) - A
        return (d_inv_sqrt[:, None] * L) * d_inv_sqrt[None, :]
    raise ValueError(f"unknown generator {generator!r}")


# ---------------------------------------------------------------------------
# Shared prepared-operator container
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class PreparedOperator:
    """
    Eigendecomposition of a generator for one slice, reusable across times.

    evals, evecs : spectrum of the real symmetric generator H.
    degrees      : node degrees of the slice adjacency (used by the damped
                   wave generator as the per-node dissipation scale).
    generator    : which generator was used (for honest downstream labeling).
    """
    evals: np.ndarray
    evecs: np.ndarray
    items: list
    degrees: np.ndarray
    generator: GeneratorName


def graph_n_items(graph) -> int:
    if hasattr(graph, "n_items"):
        return int(graph.n_items)
    if hasattr(graph, "adjacency"):
        return int(np.asarray(graph.adjacency).shape[1])
    return len(graph.items)


def _prepare(graph, slice_idx: int, generator: GeneratorName) -> PreparedOperator:
    A = _adjacency_for_slice(graph, slice_idx, graph_n_items(graph))
    H = _build_generator(A, generator)
    evals, evecs = np.linalg.eigh(H)
    return PreparedOperator(
        evals=evals, evecs=evecs, items=list(graph.items),
        degrees=_degrees(0.5 * (A + A.T)), generator=generator,
    )


# ---------------------------------------------------------------------------
# STAGE 1 — Diffusion (real density)
# ---------------------------------------------------------------------------
@dataclass
class GraphDiffusionDynamics:
    """
    Heat flow on the relatedness graph:  dρ/dt = -L ρ.

    Eigendecompose the generator once, then evolve any initial density to any
    time(s) in closed form. With generator='laplacian' the density stays
    nonnegative and mass-conserving (up to numerical noise, which
    normalize_output cleans) and relaxes toward the stationary distribution.

    For generator='normalized_laplacian' the flow is geometry-normalized and
    NOT mass-conserving in ordinary coordinates; normalize_output then rescales
    to a probability vector, which is the only sensible readout in that case.
    """
    generator: GeneratorName = "laplacian"
    normalize_output: bool = True
    clip_negative: bool = True

    def __post_init__(self):
        if self.generator == "adjacency":
            import warnings
            warnings.warn(
                "GraphDiffusionDynamics with generator='adjacency' is not a "
                "mass-conserving diffusion. Use 'laplacian' for density "
                "diffusion, or treat this only as an exploratory smoother.",
                stacklevel=2,
            )

    def prepare(self, graph, slice_idx: int) -> PreparedOperator:
        return _prepare(graph, slice_idx, self.generator)

    def _raw_evolve(self, op: PreparedOperator, rho0: np.ndarray, t: float) -> np.ndarray:
        rho0 = np.asarray(rho0, dtype=float)
        coeff = op.evecs.T @ rho0
        return op.evecs @ (np.exp(-op.evals * t) * coeff)

    def evolve(self, op: PreparedOperator, rho0: np.ndarray, t: float) -> np.ndarray:
        rho = self._raw_evolve(op, rho0, t)
        if self.normalize_output:
            rho = _normalize_density(rho, clip=self.clip_negative)
        return rho

    def trajectory(self, op: PreparedOperator, rho0: np.ndarray,
                   times: np.ndarray) -> np.ndarray:
        """Return (T, I) array of densities at each requested time."""
        out = np.empty((len(times), len(op.items)))
        for k, t in enumerate(times):
            out[k] = self.evolve(op, rho0, float(t))
        return out


# ---------------------------------------------------------------------------
# STAGE 2 — Wave (complex amplitude), optionally damped
# ---------------------------------------------------------------------------
@dataclass
class GraphWaveDynamics:
    """
    Continuous-time quantum walk on the relatedness graph.

        gamma = 0 :  i dψ/dt = L ψ                  (UNITARY — reversible)
        gamma > 0 :  i dψ/dt = (L - i·gamma·D) ψ     (DAMPED — dissipative)

    The damped variant adds a non-Hermitian, per-node dissipation proportional
    to degree (D = diag(degrees)), so amplitudes DECAY rather than ring forever.
    This makes the wave layer suitable as a phase-sensitive feature generator on
    dissipative systems (e.g. markets), addressing the objection that pure
    unitary evolution is the wrong prior for lossy, arrow-of-time data.

    *** This is a FEATURE GENERATOR, not a claim that market states evolve
    unitarily (or even by a damped wave). It produces phase-sensitive graph
    propagation features whose predictive value must be shown empirically. ***
    """
    generator: GeneratorName = "laplacian"
    gamma: float = 0.0   # dissipation rate; 0 = unitary, >0 = damped

    def __post_init__(self):
        if self.gamma < 0:
            raise ValueError("gamma (dissipation rate) must be >= 0.")

    def prepare(self, graph, slice_idx: int) -> PreparedOperator:
        return _prepare(graph, slice_idx, self.generator)

    # ----- core propagator -----
    def _U(self, op: PreparedOperator, t: float) -> np.ndarray:
        """Return the (I, I) propagator exp(-i H' t)."""
        if self.gamma == 0.0:
            # Unitary: reuse the real eigendecomposition.
            phase = np.exp(-1j * op.evals * t)
            return (op.evecs * phase[None, :]) @ op.evecs.T
        # Damped: H' = H_real - i·gamma·diag(degrees). Build and exponentiate.
        H_real = (op.evecs * op.evals[None, :]) @ op.evecs.T
        Hp = H_real.astype(complex)
        Hp[np.diag_indices_from(Hp)] -= 1j * self.gamma * op.degrees
        w, V = np.linalg.eig(-1j * t * Hp)
        return V @ np.diag(np.exp(w)) @ np.linalg.inv(V)

    def evolve(self, op: PreparedOperator, psi0: np.ndarray, t: float) -> np.ndarray:
        """Evolve a complex amplitude to time t (unitary if gamma=0)."""
        psi0 = np.asarray(psi0, dtype=complex)
        return self._U(op, t) @ psi0

    def propagator_column(self, op: PreparedOperator, i: int, t: float) -> np.ndarray:
        """Amplitude field from a unit impulse on asset i: ψ_i(t) = U(t) e_i."""
        return self._U(op, t)[:, i].copy()

    # ----- interference observables (no diffusion analogue) -----
    def interference_field(self, op: PreparedOperator, i: int, j: int,
                           t: float) -> np.ndarray:
        """
        Node-resolved interference between assets i and j at time t:
            χ_ij(node) = 2·Re( conj(ψ_i[node]) · ψ_j[node] )
        χ > 0  : the two assets' influence REINFORCES (constructive)
        χ < 0  : it CANCELS (destructive)
        For the unitary case Σ_node χ = 0 exactly (interference redistributes,
        never creates/destroys probability). Under damping the totals decay but
        the sign-changing structure — the interference signature — remains.
        """
        U = self._U(op, t)
        psi_i = U[:, i]
        psi_j = U[:, j]
        return 2.0 * np.real(np.conj(psi_i) * psi_j)

    def interference_strength(self, op: PreparedOperator, i: int, j: int,
                              t: float) -> float:
        """Total constructive interference mass between i and j at time t."""
        chi = self.interference_field(op, i, j, t)
        return float(np.maximum(chi, 0.0).sum())

    def interference_matrix(self, op: PreparedOperator, t: float) -> np.ndarray:
        """(I, I) symmetric matrix of pairwise interference strengths at time t."""
        n = len(op.items)
        U = self._U(op, t)
        M = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                chi = 2.0 * np.real(np.conj(U[:, i]) * U[:, j])
                s = float(np.maximum(chi, 0.0).sum())
                M[i, j] = s
                M[j, i] = s
        return M

    def probability(self, op: PreparedOperator, psi0: np.ndarray, t: float) -> np.ndarray:
        """|ψ(t)|² — the wave-propagated density (decays if gamma>0)."""
        psi = self.evolve(op, psi0, t)
        return np.abs(psi) ** 2

    def coherence(self, op: PreparedOperator, i: int, j: int, t: float) -> complex:
        """
        Complex coherence ⟨ψ_i(t)|ψ_j(t)⟩ summed over nodes. For unitary impulse
        starts this is exactly δ_ij, which is WHY interference must be read
        node-by-node via interference_field. Sanity-check only, not a signal.
        """
        U = self._U(op, t)
        return complex(np.vdot(U[:, i], U[:, j]))


# ---------------------------------------------------------------------------
# FEATURES — diffusion vs wave, headline = interference residual
# ---------------------------------------------------------------------------
@dataclass
class DynamicsFeatureEmitter:
    """
    Emits model-ready features by propagating an initial density through both
    the diffusion and wave layers over a set of horizons.

    The headline feature is the INTERFERENCE RESIDUAL  |ψ(t)|² - ρ(t): same
    graph, same initial state, same horizon, one propagation with interference
    and one without. If this residual carries predictive information over plain
    diffusion, the wave layer earns its place; if not, diffusion is the model
    and interference is a diagnostic. That comparison is the whole point.

    Per (slice, item, t) features:
        diffusion_density        ρ(t)_item
        wave_density             |ψ(t)|²_item
        interference_residual    |ψ(t)|² - ρ(t)         <-- headline
    Per (slice, t) scalars (broadcast to items):
        diffusion_entropy        Shannon entropy of ρ(t)
        diffusion_concentration  Σ ρ²  (inverse participation)
        wave_norm                Σ|ψ(t)|²  (1 if unitary, <1 if damped)
        total_interference       Σ_item |interference_residual|
    """
    diffusion: Optional[GraphDiffusionDynamics] = None
    wave: Optional[GraphWaveDynamics] = None

    def __post_init__(self):
        if self.diffusion is None:
            self.diffusion = GraphDiffusionDynamics(generator="laplacian")
        if self.wave is None:
            self.wave = GraphWaveDynamics(generator="laplacian", gamma=0.0)

    def transform(self, graph, slice_idx: int, rho0: np.ndarray, times):
        import pandas as pd
        rho0 = np.asarray(rho0, dtype=float)
        psi0 = np.sqrt(np.maximum(rho0, 0.0)).astype(complex)

        op_d = self.diffusion.prepare(graph, slice_idx)
        op_w = self.wave.prepare(graph, slice_idx)
        items = list(op_d.items)
        slice_label = graph.slices[slice_idx] if hasattr(graph, "slices") else slice_idx
        rows = []

        for t in np.atleast_1d(times).astype(float):
            rho = self.diffusion.evolve(op_d, rho0, t)
            wave_density = self.wave.probability(op_w, psi0, t)
            residual = wave_density - rho

            p = np.clip(rho, 1e-12, None)
            entropy = float(-(p * np.log(p)).sum())
            concentration = float((rho ** 2).sum())
            wave_norm = float(wave_density.sum())
            total_interf = float(np.abs(residual).sum())

            for k, item in enumerate(items):
                rows.append({
                    "slice": slice_label,
                    "item": item,
                    "t": float(t),
                    "diffusion_density": float(rho[k]),
                    "wave_density": float(wave_density[k]),
                    "interference_residual": float(residual[k]),
                    "diffusion_entropy": entropy,
                    "diffusion_concentration": concentration,
                    "wave_norm": wave_norm,
                    "total_interference": total_interf,
                    "wave_gamma": self.wave.gamma,
                    "generator": op_d.generator,
                })
        return pd.DataFrame(rows)
