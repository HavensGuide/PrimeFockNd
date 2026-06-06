from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .specs import BaseModeSpec


@dataclass(frozen=True)
class PrimeFockNDConfig:
    """
    Configuration for n-dimensional encoding.

    index_cols:
        All columns defining an observation. Usually slice_cols + [item_col].
    item_col:
        Column identifying nodes/items compared within each slice.
    slice_cols:
        Columns defining independent slices. If empty, all rows form one slice.
    specs:
        Column-to-mode mapping specs.
    """
    index_cols: Sequence[str]
    item_col: str
    slice_cols: Sequence[str]
    specs: Sequence[BaseModeSpec]
    dtype: str = "float64"

    def __post_init__(self) -> None:
        if not self.item_col:
            raise ValueError("item_col must be non-empty.")
        if not self.index_cols:
            raise ValueError("index_cols must be non-empty.")
        if self.item_col not in self.index_cols:
            raise ValueError("item_col must appear in index_cols.")
        for c in self.slice_cols:
            if c not in self.index_cols:
                raise ValueError(f"slice_col {c!r} must appear in index_cols.")
        if len(self.specs) == 0:
            raise ValueError("At least one mode spec is required.")
        if self.dtype not in {"float64", "float32"}:
            raise ValueError("dtype must be 'float64' or 'float32'.")


@dataclass(frozen=True)
class StabilityConfig:
    max_spectral_radius: float = 1.0
    ridge: float = 1e-6
    local_bond_threshold: float = 1.0

    def __post_init__(self) -> None:
        if self.max_spectral_radius <= 0:
            raise ValueError("max_spectral_radius must be positive.")
        if self.ridge <= 0:
            raise ValueError("ridge must be positive.")
        if self.local_bond_threshold <= 0:
            raise ValueError("local_bond_threshold must be positive.")


@dataclass(frozen=True)
class FeatureEmitterConfig:
    eps: float = 1e-12
    include_pairwise_summary: bool = True
