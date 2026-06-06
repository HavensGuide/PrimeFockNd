"""
coarse_grain.py
================

Renormalization/compression utilities for DefectState objects.

The functions here operate on the common DefectState container and therefore
work across point-cloud, arithmetic, graph, tensor, and memory backends.  The
module deliberately keeps the transformations explicit: merge nearby defects,
compress Omega, quotient events by labels, and produce a new coarse state.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .core import Defect, DefectState, EPS


@dataclass
class CoarseGrainConfig:
    """Configuration for coarse-graining defect states."""

    merge_radius: float = 0.1
    charge_tolerance: Optional[float] = None
    omega_rank: Optional[int] = None
    min_strength: float = 0.0
    preserve_observables: bool = True


@dataclass
class CoarseGrainReport:
    """Diagnostics returned by coarse-graining helpers."""

    n_events_in: Optional[int]
    n_events_out: Optional[int]
    n_defects_in: int
    n_defects_out: int
    omega_rank: Optional[int] = None
    retained_energy: Optional[float] = None
    metadata: Dict[str, Any] = None


def _n_events(events: Any) -> Optional[int]:
    try:
        return len(events)
    except Exception:
        return None


def merge_nearby_defects(
    defects: Sequence[Defect],
    *,
    radius: float = 0.1,
    charge_tolerance: Optional[float] = None,
    min_strength: float = 0.0,
) -> List[Defect]:
    """Merge defects whose locations are within ``radius`` and compatible charge."""
    ds = [d for d in defects if d.strength >= min_strength]
    if not ds:
        return []
    locs = []
    loc_valid = []
    for d in ds:
        if d.location is None:
            locs.append(np.zeros(1))
            loc_valid.append(False)
        else:
            locs.append(np.asarray(d.location, dtype=float).ravel())
            loc_valid.append(True)
    max_dim = max(len(x) for x in locs)
    X = np.zeros((len(locs), max_dim), dtype=float)
    for i, x in enumerate(locs):
        X[i, : len(x)] = x
    used = np.zeros(len(ds), dtype=bool)
    out: List[Defect] = []
    for i, d in enumerate(ds):
        if used[i]:
            continue
        group = [i]
        used[i] = True
        if loc_valid[i]:
            for j in range(i + 1, len(ds)):
                if used[j] or not loc_valid[j]:
                    continue
                if charge_tolerance is not None and abs(ds[j].charge - d.charge) > charge_tolerance:
                    continue
                if np.linalg.norm(X[j] - X[i]) <= radius:
                    group.append(j)
                    used[j] = True
        charges = np.array([ds[g].charge for g in group], dtype=float)
        strengths = np.array([ds[g].strength for g in group], dtype=float)
        weights = strengths / (strengths.sum() + EPS)
        loc = None
        if any(loc_valid[g] for g in group):
            loc = (X[group] * weights[:, None]).sum(axis=0)
        support: List[int] = []
        for g in group:
            if ds[g].support is not None:
                support.extend(ds[g].support)
        out.append(
            Defect(
                kind=ds[i].kind,
                charge=float((charges * weights).sum()),
                location=loc,
                support=tuple(sorted(set(support))) if support else None,
                strength=float(strengths.sum()),
                metadata={
                    "merged": len(group),
                    "members": [ds[g].as_dict() for g in group],
                },
            )
        )
    return out


def compress_omega(omega: Any, rank: Optional[int] = None, energy: Optional[float] = None) -> Tuple[Any, Dict[str, Any]]:
    """Low-rank SVD compression of a dense Omega/correlation matrix."""
    if omega is None:
        return None, {"omega_rank": None, "retained_energy": None}
    O = np.asarray(omega, dtype=float)
    if O.ndim != 2:
        return omega, {"omega_rank": None, "retained_energy": None, "reason": "omega_not_matrix"}
    U, s, Vt = np.linalg.svd(O, full_matrices=False)
    if rank is None:
        if energy is None:
            rank = len(s)
        else:
            cumulative = np.cumsum(s ** 2) / (np.sum(s ** 2) + EPS)
            rank = int(np.searchsorted(cumulative, energy) + 1)
    r = max(1, min(int(rank), len(s)))
    Oc = (U[:, :r] * s[:r]) @ Vt[:r]
    retained = float(np.sum(s[:r] ** 2) / (np.sum(s ** 2) + EPS))
    return Oc, {"omega_rank": r, "retained_energy": retained, "singular_values": s}


def quotient_by_labels(state: DefectState, labels: Sequence[Any]) -> DefectState:
    """Aggregate event-level generators/omega by user-supplied labels."""
    labels_arr = np.asarray(labels)
    if _n_events(state.events) is not None and len(labels_arr) != len(state.events):
        raise ValueError("labels length must match number of events")
    uniq, inv = np.unique(labels_arr, return_inverse=True)
    gen = state.generators
    gen_out = None
    if gen is not None:
        G = np.asarray(gen, dtype=float)
        if G.ndim >= 2 and len(G) == len(labels_arr):
            gen_out = np.zeros((len(uniq),) + G.shape[1:], dtype=float)
            counts = np.bincount(inv).astype(float)
            for i in range(len(labels_arr)):
                gen_out[inv[i]] += G[i]
            gen_out /= counts.reshape((-1,) + (1,) * (G.ndim - 1))
    omega_out = None
    if state.omega is not None:
        O = np.asarray(state.omega, dtype=float)
        if O.shape == (len(labels_arr), len(labels_arr)):
            omega_out = np.zeros((len(uniq), len(uniq)), dtype=float)
            counts = np.zeros((len(uniq), len(uniq)), dtype=float)
            for i in range(len(labels_arr)):
                for j in range(len(labels_arr)):
                    omega_out[inv[i], inv[j]] += O[i, j]
                    counts[inv[i], inv[j]] += 1.0
            omega_out /= counts + EPS
    return DefectState(
        events=list(uniq),
        generators=gen_out if gen_out is not None else state.generators,
        connection=state.connection,
        fusion_rules=state.fusion_rules,
        projectors=state.projectors,
        omega=omega_out if omega_out is not None else state.omega,
        defects=list(state.defects),
        observables=dict(state.observables),
        metadata={**state.metadata, "coarse_grain": "quotient_by_labels", "n_labels": len(uniq)},
        valid=state.valid,
    )


def coarse_grain_state(
    state: DefectState,
    config: Optional[CoarseGrainConfig] = None,
    *,
    labels: Optional[Sequence[Any]] = None,
) -> Tuple[DefectState, CoarseGrainReport]:
    """Coarse-grain a DefectState by optional quotient, defect merge, and Omega compression."""
    cfg = config or CoarseGrainConfig()
    src = quotient_by_labels(state, labels) if labels is not None else state
    defects_out = merge_nearby_defects(
        src.defects,
        radius=cfg.merge_radius,
        charge_tolerance=cfg.charge_tolerance,
        min_strength=cfg.min_strength,
    )
    omega_out, info = compress_omega(src.omega, rank=cfg.omega_rank)
    observables = dict(src.observables) if cfg.preserve_observables else {}
    observables.update({
        "n_defects_coarse": len(defects_out),
        "omega_rank_coarse": info.get("omega_rank"),
        "omega_retained_energy": info.get("retained_energy"),
    })
    out = DefectState(
        events=src.events,
        generators=src.generators,
        connection=src.connection,
        fusion_rules=src.fusion_rules,
        projectors=src.projectors,
        omega=omega_out,
        defects=defects_out,
        observables=observables,
        metadata={**src.metadata, "coarse_grain": "state", **{k: v for k, v in info.items() if k != "singular_values"}},
        valid=src.valid,
    )
    report = CoarseGrainReport(
        n_events_in=_n_events(state.events),
        n_events_out=_n_events(out.events),
        n_defects_in=len(state.defects),
        n_defects_out=len(defects_out),
        omega_rank=info.get("omega_rank"),
        retained_energy=info.get("retained_energy"),
        metadata={"singular_values": info.get("singular_values")},
    )
    return out, report


def register_vm_ops(vm: Any) -> Any:
    """Register coarse-graining helpers on a DefectVM instance."""

    def _coarse_state(**kwargs: Any):
        cfg = kwargs.pop("config", None)
        if cfg is None:
            cfg = CoarseGrainConfig(**{k: kwargs.pop(k) for k in list(kwargs.keys()) if k in CoarseGrainConfig.__annotations__})
        return coarse_grain_state(vm.state, cfg, **kwargs)

    vm.register("coarse_grain", _coarse_state)
    vm.register("merge_defects", lambda **kw: merge_nearby_defects(vm.state.defects, **kw))
    vm.register("compress_omega", lambda **kw: compress_omega(vm.state.omega, **kw))
    return vm


__all__ = [
    "CoarseGrainConfig",
    "CoarseGrainReport",
    "merge_nearby_defects",
    "compress_omega",
    "quotient_by_labels",
    "coarse_grain_state",
    "register_vm_ops",
]
