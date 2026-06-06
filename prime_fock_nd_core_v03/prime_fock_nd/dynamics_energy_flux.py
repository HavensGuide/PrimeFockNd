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

  STAGE 2 — GraphWaveDynamics  (complex amplitude, data-driven open flow)
      Closed:          i dψ/dt = L ψ
      Open/data-flow:  i dψ/dt = (L + i·diag(q)) ψ
      ψ = √ρ·e^{iS} is the Madelung/symplectic lift. Complex amplitudes let
      two assets' spreading influence INTERFERE — reinforce at some nodes,
      cancel at others. q is not a tunable damping constant: it is a per-node
      inflow/outflow field inferred from the observed slice-to-slice graph
      change, so energy can enter or leave where the data geometry says it did.

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
prior for arrow-of-time markets. Stage 1 (diffusion) is the defensible model;
Stage 2 is exploratory and its value must be demonstrated empirically. The
open-flow wave option prevents the wave layer from being forced to be lossless,
without introducing an exogenous damping knob.
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
    degrees      : node degrees of the slice adjacency.
    energy_flux  : per-node observed inflow/outflow field. Positive values are
                   source/inflow; negative values are sink/outflow.
    generator    : which generator was used (for honest downstream labeling).
    """
    evals: np.ndarray
    evecs: np.ndarray
    items: list
    degrees: np.ndarray
    energy_flux: np.ndarray
    generator: GeneratorName


def graph_n_items(graph) -> int:
    if hasattr(graph, "n_items"):
        return int(graph.n_items)
    if hasattr(graph, "adjacency"):
        return int(np.asarray(graph.adjacency).shape[1])
    return len(graph.items)



def _slice_degree_flux(graph, slice_idx: int, n_items: int, *, eps: float = 1e-12) -> np.ndarray:
    """
    Data-derived open-system source/sink field q for one slice.

    q_i > 0 means node i gained graph energy/capacity versus the previous
    slice; q_i < 0 means it lost graph energy/capacity. The normalization keeps
    the field dimensionless and comparable across slices without introducing a
    tuned damping coefficient.
    """
    if slice_idx <= 0:
        return np.zeros(n_items, dtype=float)

    A_now = _adjacency_for_slice(graph, slice_idx, n_items)
    A_prev = _adjacency_for_slice(graph, slice_idx - 1, n_items)
    d_now = _degrees(A_now)
    d_prev = _degrees(A_prev)
    raw = d_now - d_prev
    scale = np.mean(np.abs(d_prev)) + eps
    return raw / scale


def _resolve_energy_flux(
    graph,
    slice_idx: int,
    n_items: int,
    energy_flux: Optional[np.ndarray] = None,
    *,
    flux_scale: float = 1.0,
    flux_clip: Optional[float] = None,
) -> np.ndarray:
    """Resolve the per-node open-flow field used by the wave layer."""
    if energy_flux is None:
        q = _slice_degree_flux(graph, slice_idx, n_items)
    else:
        arr = np.asarray(energy_flux, dtype=float)
        if arr.shape == (n_items,):
            q = arr.copy()
        elif arr.ndim == 2 and arr.shape[1] == n_items:
            q = arr[slice_idx].copy()
        else:
            raise ValueError(
                f"energy_flux must have shape {(n_items,)} or (S, {n_items}), got {arr.shape}."
            )

    if not np.all(np.isfinite(q)):
        raise ValueError("energy_flux contains non-finite values.")
    q = float(flux_scale) * q
    if flux_clip is not None:
        c = float(flux_clip)
        if c <= 0:
            raise ValueError("flux_clip must be positive when provided.")
        q = np.clip(q, -c, c)
    return q

def _prepare(graph, slice_idx: int, generator: GeneratorName, *, energy_flux: Optional[np.ndarray] = None, flux_scale: float = 1.0, flux_clip: Optional[float] = None) -> PreparedOperator:
    n_items = graph_n_items(graph)
    A = _adjacency_for_slice(graph, slice_idx, n_items)
    H = _build_generator(A, generator)
    evals, evecs = np.linalg.eigh(H)
    degrees = _degrees(0.5 * (A + A.T))
    q = _resolve_energy_flux(
        graph, slice_idx, n_items,
        energy_flux=energy_flux,
        flux_scale=flux_scale,
        flux_clip=flux_clip,
    )
    return PreparedOperator(
        evals=evals, evecs=evecs, items=list(graph.items),
        degrees=degrees, energy_flux=q, generator=generator,
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
    Continuous-time amplitude walk on the relatedness graph.

        closed:          i dψ/dt = H ψ
        open/data-flow:  i dψ/dt = (H + i·diag(q)) ψ

    q is a per-node source/sink field. By default it is inferred from the
    observed change in node degree from the previous slice:
        q_i = (degree_i[s] - degree_i[s-1]) / mean(abs(degree[s-1]))

    This removes the old scalar damping knob. Energy can flow in or out by node,
    and the flow itself is emitted downstream as a feature. If you want a fully
    closed wave A/B test, pass use_energy_flux=False.
    """
    generator: GeneratorName = "laplacian"
    use_energy_flux: bool = True
    energy_flux: Optional[np.ndarray] = None
    flux_scale: float = 1.0
    flux_clip: Optional[float] = 3.0

    def __post_init__(self):
        if self.flux_scale <= 0:
            raise ValueError("flux_scale must be positive.")
        if self.flux_clip is not None and self.flux_clip <= 0:
            raise ValueError("flux_clip must be positive when provided.")

    def prepare(self, graph, slice_idx: int) -> PreparedOperator:
        if not self.use_energy_flux:
            n = graph_n_items(graph)
            q = np.zeros(n, dtype=float)
        else:
            q = self.energy_flux
        return _prepare(
            graph,
            slice_idx,
            self.generator,
            energy_flux=q,
            flux_scale=self.flux_scale,
            flux_clip=self.flux_clip,
        )

    # ----- core propagator -----
    def _U(self, op: PreparedOperator, t: float) -> np.ndarray:
        """Return the (I, I) propagator exp(-i H' t)."""
        H_real = (op.evecs * op.evals[None, :]) @ op.evecs.T
        q = np.asarray(op.energy_flux, dtype=float)

        if (not self.use_energy_flux) or np.allclose(q, 0.0):
            phase = np.exp(-1j * op.evals * t)
            return (op.evecs * phase[None, :]) @ op.evecs.T

        # Open system: H' = H + i·diag(q). Since exp(-i H' t) contains exp(q t),
        # positive q is observed inflow/source and negative q is outflow/sink.
        Hp = H_real.astype(complex)
        Hp[np.diag_indices_from(Hp)] += 1j * q
        w, V = np.linalg.eig(-1j * t * Hp)
        return V @ np.diag(np.exp(w)) @ np.linalg.inv(V)

    def evolve(self, op: PreparedOperator, psi0: np.ndarray, t: float) -> np.ndarray:
        """Evolve a complex amplitude to time t."""
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
        """|ψ(t)|² — the wave-propagated density with observed open-flow."""
        psi = self.evolve(op, psi0, t)
        return np.abs(psi) ** 2

    def coherence(self, op: PreparedOperator, i: int, j: int, t: float) -> complex:
        """Complex coherence diagnostic."""
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
        wave_norm                Σ|ψ(t)|²  (changes when open flow is active)
        total_interference       Σ_item |interference_residual|
        energy_flux              node-level observed source/sink field
        net_energy_flux          Σ_i q_i
    """
    diffusion: Optional[GraphDiffusionDynamics] = None
    wave: Optional[GraphWaveDynamics] = None

    def __post_init__(self):
        if self.diffusion is None:
            self.diffusion = GraphDiffusionDynamics(generator="laplacian")
        if self.wave is None:
            self.wave = GraphWaveDynamics(generator="laplacian", use_energy_flux=True)

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
                    "energy_flux": float(op_w.energy_flux[k]),
                    "net_energy_flux": float(op_w.energy_flux.sum()),
                    "source_energy_flux": float(np.maximum(op_w.energy_flux, 0.0).sum()),
                    "sink_energy_flux": float(np.maximum(-op_w.energy_flux, 0.0).sum()),
                    "generator": op_d.generator,
                })
        return pd.DataFrame(rows)
