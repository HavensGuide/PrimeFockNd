from __future__ import annotations

import numpy as np
import pandas as pd


def require_columns(df: pd.DataFrame, columns: list[str]) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def ensure_square_matrix(matrix: np.ndarray, n: int, name: str) -> np.ndarray:
    arr = np.asarray(matrix, dtype=float)
    if arr.shape != (n, n):
        raise ValueError(f"{name} must have shape {(n, n)}, got {arr.shape}.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values.")
    return arr


def ensure_symmetric(matrix: np.ndarray, name: str, atol: float = 1e-8) -> None:
    if not np.allclose(matrix, matrix.T, atol=atol):
        raise ValueError(f"{name} must be symmetric.")


def normalize_key(value):
    if isinstance(value, tuple):
        return value
    return (value,)


def sorted_unique_tuples(df: pd.DataFrame, cols: list[str]) -> list[tuple]:
    if not cols:
        return [("__all__",)]
    return sorted([tuple(x) for x in df[cols].drop_duplicates().to_numpy().tolist()])
