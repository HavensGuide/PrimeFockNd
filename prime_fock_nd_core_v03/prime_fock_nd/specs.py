from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, Sequence

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class BaseModeSpec:
    """
    Base class for a column-to-Fock-mode mapping.
    """
    name: str
    missing_policy: Literal["zero", "raise"] = "zero"
    exponent_scale: float = 1.0

    def required_columns(self) -> list[str]:
        return [self.name]

    def fit(self, df: pd.DataFrame) -> "BaseModeSpec":
        return self

    def width(self) -> int:
        raise NotImplementedError

    def encode_row(self, row: pd.Series) -> dict[int, float]:
        raise NotImplementedError

    def mode_metadata(self, local_idx: int) -> dict[str, Any]:
        return {"source": self.name, "local_idx": local_idx, "spec_type": type(self).__name__}

    def _handle_missing(self, value: Any) -> bool:
        if pd.isna(value):
            if self.missing_policy == "raise":
                raise ValueError(f"Missing value for {self.name!r}.")
            return True
        return False

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("Spec name must be non-empty.")
        if self.exponent_scale <= 0:
            raise ValueError("exponent_scale must be positive.")
        if self.missing_policy not in {"zero", "raise"}:
            raise ValueError("missing_policy must be 'zero' or 'raise'.")


@dataclass(frozen=True)
class NumericBinningSpec(BaseModeSpec):
    """
    Maps a numeric column into one active mode using bin edges.
    """
    bins: Sequence[float] = field(default_factory=list)

    def __post_init__(self) -> None:
        super().__post_init__()
        if any(self.bins[i] >= self.bins[i + 1] for i in range(len(self.bins) - 1)):
            raise ValueError(f"Bins for {self.name!r} must be strictly increasing.")

    def width(self) -> int:
        return len(self.bins) + 1

    def encode_row(self, row: pd.Series) -> dict[int, float]:
        value = row[self.name]
        if self._handle_missing(value):
            return {}
        idx = int(np.digitize(float(value), np.asarray(self.bins, dtype=float), right=False))
        return {idx: float(self.exponent_scale)}

    def mode_metadata(self, local_idx: int) -> dict[str, Any]:
        return {
            "source": self.name,
            "local_idx": local_idx,
            "spec_type": type(self).__name__,
            "bins": list(self.bins),
            "exponent_scale": self.exponent_scale,
        }


@dataclass(frozen=True)
class CategoricalSpec(BaseModeSpec):
    """
    Maps a categorical column into one active mode per category.

    categories=None means categories are inferred during fit_transform.
    """
    categories: Sequence[Any] | None = None
    unknown_policy: Literal["zero", "raise"] = "zero"

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.unknown_policy not in {"zero", "raise"}:
            raise ValueError("unknown_policy must be 'zero' or 'raise'.")
        if self.categories is not None and len(set(self.categories)) != len(self.categories):
            raise ValueError("categories must be unique.")

    def fit(self, df: pd.DataFrame) -> "CategoricalSpec":
        if self.categories is not None:
            return self
        values = [x for x in pd.Index(df[self.name].dropna()).unique().tolist()]
        values = sorted(values)
        return CategoricalSpec(
            name=self.name,
            missing_policy=self.missing_policy,
            exponent_scale=self.exponent_scale,
            categories=values,
            unknown_policy=self.unknown_policy,
        )

    def width(self) -> int:
        if self.categories is None:
            raise ValueError("CategoricalSpec must be fitted before width() is known.")
        return len(self.categories)

    def encode_row(self, row: pd.Series) -> dict[int, float]:
        value = row[self.name]
        if self._handle_missing(value):
            return {}
        assert self.categories is not None
        try:
            idx = list(self.categories).index(value)
        except ValueError:
            if self.unknown_policy == "raise":
                raise ValueError(f"Unknown category {value!r} for {self.name!r}.")
            return {}
        return {idx: float(self.exponent_scale)}

    def mode_metadata(self, local_idx: int) -> dict[str, Any]:
        assert self.categories is not None
        return {
            "source": self.name,
            "local_idx": local_idx,
            "spec_type": type(self).__name__,
            "category": self.categories[local_idx],
            "exponent_scale": self.exponent_scale,
        }


@dataclass(frozen=True)
class MultiHotSpec(BaseModeSpec):
    """
    Maps a column containing an iterable of tags into multiple active modes.
    """
    vocabulary: Sequence[Any] | None = None
    separator: str | None = None
    unknown_policy: Literal["zero", "raise"] = "zero"

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.unknown_policy not in {"zero", "raise"}:
            raise ValueError("unknown_policy must be 'zero' or 'raise'.")
        if self.vocabulary is not None and len(set(self.vocabulary)) != len(self.vocabulary):
            raise ValueError("vocabulary must be unique.")

    def fit(self, df: pd.DataFrame) -> "MultiHotSpec":
        if self.vocabulary is not None:
            return self
        vocab = set()
        for value in df[self.name].dropna():
            for token in self._tokens(value):
                vocab.add(token)
        return MultiHotSpec(
            name=self.name,
            missing_policy=self.missing_policy,
            exponent_scale=self.exponent_scale,
            vocabulary=sorted(vocab),
            separator=self.separator,
            unknown_policy=self.unknown_policy,
        )

    def width(self) -> int:
        if self.vocabulary is None:
            raise ValueError("MultiHotSpec must be fitted before width() is known.")
        return len(self.vocabulary)

    def encode_row(self, row: pd.Series) -> dict[int, float]:
        value = row[self.name]
        if self._handle_missing(value):
            return {}
        assert self.vocabulary is not None
        vocab = list(self.vocabulary)
        out = {}
        for token in self._tokens(value):
            try:
                idx = vocab.index(token)
            except ValueError:
                if self.unknown_policy == "raise":
                    raise ValueError(f"Unknown token {token!r} for {self.name!r}.")
                continue
            out[idx] = float(self.exponent_scale)
        return out

    def _tokens(self, value: Any) -> list[Any]:
        if self.separator is not None and isinstance(value, str):
            return [x.strip() for x in value.split(self.separator) if x.strip()]
        if isinstance(value, str):
            return [value]
        try:
            return list(value)
        except TypeError:
            return [value]

    def mode_metadata(self, local_idx: int) -> dict[str, Any]:
        assert self.vocabulary is not None
        return {
            "source": self.name,
            "local_idx": local_idx,
            "spec_type": type(self).__name__,
            "token": self.vocabulary[local_idx],
            "exponent_scale": self.exponent_scale,
        }


@dataclass(frozen=True)
class PreEncodedSpec(BaseModeSpec):
    """
    Maps a set of already numeric exponent columns directly into Fock modes.

    Example:
        PreEncodedSpec(name="latent", columns=["z1", "z2", "z3"])

    The `name` is used as a namespace; required columns are `columns`.
    """
    columns: Sequence[str] = field(default_factory=list)
    clip_negative: bool = True

    def required_columns(self) -> list[str]:
        return list(self.columns)

    def width(self) -> int:
        return len(self.columns)

    def encode_row(self, row: pd.Series) -> dict[int, float]:
        out = {}
        for i, col in enumerate(self.columns):
            value = row[col]
            if pd.isna(value):
                if self.missing_policy == "raise":
                    raise ValueError(f"Missing value for {col!r}.")
                continue
            val = float(value) * self.exponent_scale
            if self.clip_negative:
                val = max(0.0, val)
            if val != 0.0:
                out[i] = val
        return out

    def mode_metadata(self, local_idx: int) -> dict[str, Any]:
        return {
            "source": self.name,
            "column": self.columns[local_idx],
            "local_idx": local_idx,
            "spec_type": type(self).__name__,
            "exponent_scale": self.exponent_scale,
        }

    def __post_init__(self) -> None:
        super().__post_init__()
        if not self.columns:
            raise ValueError("PreEncodedSpec.columns must be non-empty.")
        if len(set(self.columns)) != len(self.columns):
            raise ValueError("PreEncodedSpec.columns must be unique.")
