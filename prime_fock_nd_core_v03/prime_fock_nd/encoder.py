from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Tuple

import numpy as np
import pandas as pd

from .config import PrimeFockNDConfig
from .primes import first_primes
from .specs import BaseModeSpec
from .validation import require_columns, sorted_unique_tuples


@dataclass(frozen=True)
class EncodedND:
    """
    N-dimensional encoded Fock tensor.

    exponents:
        Shape (S, I, P), where S=slices, I=items, P=prime modes.
    slices:
        Slice keys length S. Each key is a tuple.
    items:
        Item labels length I.
    primes:
        Prime values length P.
    mode_metadata:
        Dict mapping global mode index to metadata.
    config:
        Fitted PrimeFockNDConfig.
    """
    exponents: np.ndarray
    slices: list[tuple]
    items: list[Any]
    primes: np.ndarray
    mode_metadata: Dict[int, dict]
    config: PrimeFockNDConfig

    @property
    def shape(self) -> Tuple[int, int, int]:
        return self.exponents.shape


class PrimeFockNDEncoder:
    """
    General dataframe-to-Fock encoder.

    This supports arbitrary dimensions by collapsing all non-item index
    dimensions into slice keys. Within each slice, items become graph nodes.
    """

    def __init__(self, config: PrimeFockNDConfig):
        self.config = config
        self.fitted_specs: list[BaseModeSpec] | None = None
        self.mode_slices: dict[str, slice] = {}
        self.mode_metadata_: Dict[int, dict] = {}

    def fit(self, df: pd.DataFrame) -> "PrimeFockNDEncoder":
        self._validate_columns(df)
        self.fitted_specs = [spec.fit(df) for spec in self.config.specs]
        self._build_mode_map()
        return self

    def fit_transform(self, df: pd.DataFrame) -> EncodedND:
        return self.fit(df).transform(df)

    def transform(self, df: pd.DataFrame) -> EncodedND:
        if self.fitted_specs is None:
            raise RuntimeError("PrimeFockNDEncoder must be fitted before transform().")
        self._validate_columns(df)

        slice_cols = list(self.config.slice_cols)
        item_col = self.config.item_col
        index_cols = list(self.config.index_cols)

        slices = sorted_unique_tuples(df, slice_cols)
        items = sorted(pd.Index(df[item_col]).unique().tolist())

        slice_to_i = {s: i for i, s in enumerate(slices)}
        item_to_i = {item: i for i, item in enumerate(items)}

        p_count = sum(spec.width() for spec in self.fitted_specs)
        primes = np.asarray(first_primes(p_count), dtype=np.int64)
        dtype = np.float64 if self.config.dtype == "float64" else np.float32
        exponents = np.zeros((len(slices), len(items), p_count), dtype=dtype)

        dupes = df.duplicated(index_cols)
        if dupes.any():
            example = df.loc[dupes, index_cols].head(3).to_dict("records")
            raise ValueError(f"Duplicate index rows are not allowed. Examples: {example}")

        for _, row in df.iterrows():
            s_key = tuple(row[c] for c in slice_cols) if slice_cols else ("__all__",)
            s_idx = slice_to_i[s_key]
            i_idx = item_to_i[row[item_col]]

            for spec in self.fitted_specs:
                sl = self.mode_slices[spec.name]
                local = spec.encode_row(row)
                for local_idx, value in local.items():
                    exponents[s_idx, i_idx, sl.start + local_idx] = value

        return EncodedND(
            exponents=exponents,
            slices=slices,
            items=items,
            primes=primes,
            mode_metadata=dict(self.mode_metadata_),
            config=self.config,
        )

    def _validate_columns(self, df: pd.DataFrame) -> None:
        cols = list(self.config.index_cols)
        for spec in self.config.specs:
            cols.extend(spec.required_columns())
        require_columns(df, sorted(set(cols)))

    def _build_mode_map(self) -> None:
        assert self.fitted_specs is not None
        self.mode_slices = {}
        self.mode_metadata_ = {}

        start = 0
        for spec in self.fitted_specs:
            width = spec.width()
            self.mode_slices[spec.name] = slice(start, start + width)
            for local_idx in range(width):
                idx = start + local_idx
                meta = spec.mode_metadata(local_idx)
                meta["global_mode"] = idx
                self.mode_metadata_[idx] = meta
            start += width
