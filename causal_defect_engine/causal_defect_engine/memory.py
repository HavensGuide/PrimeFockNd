"""
memory.py
=========

Pleroma Memory Algebra utilities.

The memory module treats memory as a persistent algebraic trace rather than as a
list of stored copies.  It is deliberately backend-agnostic: events can be raw
vectors, symbols, defects, or outputs from another CausalDefectEngine.

Core idea
---------
A memory state is a compressed trace

    M_{t+1} = decay * M_t + write_weight * e_t + assoc_weight * compose(M_t, e_t)

where composition defaults to a low-rank outer-association update.  Recall is a
cue-to-trace projection.  Compression is SVD-based and optional.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np

from .core import CausalDefectEngine, Defect, DefectState, EPS, normalize_mode

try:
    from numba import njit, prange  # type: ignore
    NUMBA_AVAILABLE = True
except Exception:  # pragma: no cover
    NUMBA_AVAILABLE = False

    def njit(*args, **kwargs):  # type: ignore
        def deco(fn):
            return fn
        if args and callable(args[0]):
            return args[0]
        return deco

    prange = range  # type: ignore


@njit(cache=True, fastmath=True)
def _cosine_scores_numba(memory: np.ndarray, cue: np.ndarray) -> np.ndarray:
    n, d = memory.shape
    out = np.zeros(n, dtype=np.float64)
    cue_norm = 0.0
    for j in range(d):
        cue_norm += cue[j] * cue[j]
    cue_norm = np.sqrt(cue_norm) + 1e-12
    for i in prange(n):
        dot = 0.0
        mn = 0.0
        for j in range(d):
            dot += memory[i, j] * cue[j]
            mn += memory[i, j] * memory[i, j]
        out[i] = dot / ((np.sqrt(mn) + 1e-12) * cue_norm)
    return out


@dataclass
class MemoryConfig:
    """Configuration for :class:`PleromaMemoryEngine`."""

    dim: Optional[int] = None
    decay: float = 0.98
    write_weight: float = 1.0
    assoc_weight: float = 0.15
    max_items: Optional[int] = None
    compression_rank: Optional[int] = None
    normalize_writes: bool = True
    mode: str = "auto"  # legacy/exact keeps append-only; efficient compresses/prunes
    use_numba: bool = True


@dataclass
class MemoryTrace:
    """A single memory event and its embedding/metadata."""

    key: str
    vector: np.ndarray
    weight: float = 1.0
    timestamp: Optional[float] = None
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class RecallResult:
    """Cue recall result."""

    indices: np.ndarray
    keys: List[str]
    scores: np.ndarray
    vectors: np.ndarray
    metadata: List[Dict[str, Any]]


class PleromaMemoryEngine(CausalDefectEngine):
    """Associative trace memory over defect/event embeddings.

    The engine stores a weighted memory matrix ``M`` and an association matrix
    ``A``.  ``M`` is the compressed trace bank. ``A`` is a low-rank-like
    association operator accumulated from cue/value outer products.

    This is intentionally not a neural network.  It is a small algebraic memory
    substrate with explicit write, compose, recall, forget, project, and
    compress operations.
    """

    def __init__(self, config: Optional[MemoryConfig] = None, **kwargs: Any):
        self.config = config or MemoryConfig()
        super().__init__(
            mode=self.config.mode,
            use_numba=self.config.use_numba,
            **kwargs,
        )
        self.traces: List[MemoryTrace] = []
        self.memory = np.zeros((0, self.config.dim or 0), dtype=np.float64)
        self.weights = np.zeros(0, dtype=np.float64)
        self.assoc = np.zeros((self.config.dim or 0, self.config.dim or 0), dtype=np.float64)

    # ------------------------------------------------------------------
    # Fitting / backend API
    # ------------------------------------------------------------------

    def fit(self, events: Any = None, **kwargs: Any) -> "PleromaMemoryEngine":
        if events is not None:
            for i, ev in enumerate(events):
                self.write(ev, key=kwargs.get("key_prefix", "e") + str(i))
        self._refresh_state()
        return self

    def estimate_generators(self, **kwargs: Any) -> np.ndarray:
        return self.memory

    def build_connection(self, **kwargs: Any) -> np.ndarray:
        return self.assoc

    def detect_defects(self, threshold: float = 0.85, **kwargs: Any) -> List[Defect]:
        """Detect memory collisions/attractors as high-similarity trace clusters."""
        defects: List[Defect] = []
        if len(self.memory) < 2:
            self.state.defects = defects
            return defects
        X = self._row_normalize(self.memory)
        S = X @ X.T
        np.fill_diagonal(S, 0.0)
        used = np.zeros(len(X), dtype=bool)
        for i in range(len(X)):
            if used[i]:
                continue
            nbrs = np.where(S[i] >= threshold)[0]
            if len(nbrs) == 0:
                continue
            group = np.unique(np.r_[i, nbrs])
            used[group] = True
            charge = float(np.mean(S[np.ix_(group, group)])) if len(group) > 1 else 1.0
            loc = self.memory[group].mean(axis=0)
            defects.append(
                Defect(
                    kind="memory_attractor",
                    charge=charge,
                    location=loc,
                    support=tuple(int(x) for x in group),
                    strength=float(len(group)),
                    metadata={"keys": [self.traces[j].key for j in group]},
                )
            )
        self.state.defects = defects
        return defects

    # ------------------------------------------------------------------
    # Memory operations
    # ------------------------------------------------------------------

    def vectorize(self, event: Any) -> np.ndarray:
        """Convert an event into a numeric vector.

        Accepted inputs:
        - numpy/list vector
        - mapping with a ``vector`` key
        - :class:`Defect` with a location or metadata vector
        - scalar/symbol, hashed into a deterministic random-like vector
        """
        if isinstance(event, MemoryTrace):
            return np.asarray(event.vector, dtype=np.float64)
        if isinstance(event, Defect):
            if event.location is not None:
                return np.asarray(event.location, dtype=np.float64).ravel()
            if "vector" in event.metadata:
                return np.asarray(event.metadata["vector"], dtype=np.float64).ravel()
            vals = [event.charge, event.strength]
            return np.asarray(vals, dtype=np.float64)
        if isinstance(event, Mapping) and "vector" in event:
            return np.asarray(event["vector"], dtype=np.float64).ravel()
        arr = np.asarray(event, dtype=object)
        if arr.dtype != object:
            return arr.astype(np.float64).ravel()
        # deterministic symbol embedding
        dim = self.config.dim or 16
        seed = abs(hash(str(event))) % (2**32)
        rng = np.random.default_rng(seed)
        return rng.standard_normal(dim)

    def _ensure_dim(self, d: int) -> None:
        if self.config.dim is None:
            self.config.dim = int(d)
        if int(d) != int(self.config.dim):
            raise ValueError(f"event vector has dim {d}, expected {self.config.dim}")
        if self.assoc.shape != (d, d):
            self.assoc = np.zeros((d, d), dtype=np.float64)
        if self.memory.shape[1] == 0:
            self.memory = np.zeros((0, d), dtype=np.float64)

    @staticmethod
    def _row_normalize(X: np.ndarray) -> np.ndarray:
        return X / (np.linalg.norm(X, axis=1, keepdims=True) + EPS)

    def write(
        self,
        event: Any,
        *,
        key: Optional[str] = None,
        weight: Optional[float] = None,
        timestamp: Optional[float] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> MemoryTrace:
        """Write an event into memory and update the association operator."""
        x = self.vectorize(event).astype(np.float64)
        self._ensure_dim(len(x))
        if self.config.normalize_writes:
            x = x / (np.linalg.norm(x) + EPS)
        w = float(self.config.write_weight if weight is None else weight)

        # decay old traces
        if len(self.weights):
            self.weights *= self.config.decay
            self.memory *= self.config.decay
            # algebraic association trace: old memory composes with new event
            mean_trace = self.memory.mean(axis=0) if len(self.memory) else np.zeros_like(x)
            self.assoc = self.config.decay * self.assoc + self.config.assoc_weight * np.outer(mean_trace, x)
        else:
            self.assoc = self.config.assoc_weight * np.outer(x, x)

        trace = MemoryTrace(
            key=key or f"m{len(self.traces)}",
            vector=x.copy(),
            weight=w,
            timestamp=timestamp,
            metadata=dict(metadata or {}),
        )
        self.traces.append(trace)
        self.memory = np.vstack([self.memory, w * x])
        self.weights = np.r_[self.weights, w]

        self._maybe_prune_or_compress()
        self._refresh_state()
        return trace

    def compose(self, a: Any, b: Any) -> np.ndarray:  # type: ignore[override]
        """Compose two cues/traces via the learned association operator."""
        va = self.vectorize(a)
        vb = self.vectorize(b)
        self._ensure_dim(len(va))
        return va + vb + self.config.assoc_weight * (va @ self.assoc + self.assoc @ vb)

    def recall(self, cue: Any, *, k: int = 5, use_assoc: bool = True) -> RecallResult:
        """Recall the top-k traces associated with ``cue``."""
        if len(self.memory) == 0:
            return RecallResult(np.zeros(0, dtype=int), [], np.zeros(0), np.zeros((0, self.config.dim or 0)), [])
        q = self.vectorize(cue).astype(np.float64)
        self._ensure_dim(len(q))
        if self.config.normalize_writes:
            q = q / (np.linalg.norm(q) + EPS)
        if use_assoc:
            q = q + q @ self.assoc
        if self.use_numba and NUMBA_AVAILABLE:
            scores = _cosine_scores_numba(self.memory.astype(np.float64), q.astype(np.float64))
        else:
            scores = self._row_normalize(self.memory) @ (q / (np.linalg.norm(q) + EPS))
        kk = min(int(k), len(scores))
        idx = np.argpartition(-scores, kk - 1)[:kk]
        idx = idx[np.argsort(-scores[idx])]
        return RecallResult(
            indices=idx,
            keys=[self.traces[int(i)].key for i in idx],
            scores=scores[idx],
            vectors=self.memory[idx],
            metadata=[self.traces[int(i)].metadata for i in idx],
        )

    def forget(self, *, threshold: float = 1e-4, keys: Optional[Sequence[str]] = None) -> int:
        """Remove traces by weight threshold or explicit keys."""
        if keys is not None:
            keyset = set(keys)
            keep = np.array([tr.key not in keyset for tr in self.traces], dtype=bool)
        else:
            keep = self.weights >= float(threshold)
        removed = int((~keep).sum())
        self.traces = [tr for tr, ok in zip(self.traces, keep) if ok]
        self.memory = self.memory[keep]
        self.weights = self.weights[keep]
        self._refresh_state()
        return removed

    def compress(self, rank: Optional[int] = None) -> np.ndarray:
        """SVD-compress the memory matrix in-place and return singular values."""
        if len(self.memory) == 0:
            return np.zeros(0)
        r = int(rank or self.config.compression_rank or min(self.memory.shape))
        U, s, Vt = np.linalg.svd(self.memory, full_matrices=False)
        r = max(1, min(r, len(s)))
        self.memory = (U[:, :r] * s[:r]) @ Vt[:r]
        self.state.metadata["memory_compression_rank"] = r
        self.state.observables["memory_singular_values"] = s
        return s

    def project(self, sector: Union[Sequence[int], np.ndarray], x: Any = None) -> np.ndarray:  # type: ignore[override]
        """Project memory or vector onto a coordinate sector."""
        idx = np.asarray(sector, dtype=int)
        if x is None:
            return self.memory[:, idx]
        v = self.vectorize(x)
        out = np.zeros_like(v)
        out[idx] = v[idx]
        return out

    def _maybe_prune_or_compress(self) -> None:
        mode = normalize_mode(self.config.mode)
        if self.config.max_items is not None and len(self.traces) > self.config.max_items:
            # keep highest-weight recent-ish traces
            order = np.argsort(-self.weights)[: self.config.max_items]
            order = np.sort(order)
            self.traces = [self.traces[int(i)] for i in order]
            self.memory = self.memory[order]
            self.weights = self.weights[order]
        if mode in {"efficient", "local", "auto"} and self.config.compression_rank is not None:
            if len(self.memory) > self.config.compression_rank:
                self.compress(self.config.compression_rank)

    def _refresh_state(self) -> None:
        self.state = DefectState(
            events=[tr.key for tr in self.traces],
            generators=self.memory,
            connection=self.assoc,
            omega=(self.memory.T @ self.memory if len(self.memory) else None),
            defects=self.state.defects if self.state.valid else [],
            observables={
                "n_memory_items": len(self.traces),
                "memory_dim": self.config.dim,
                "memory_norm": float(np.linalg.norm(self.memory)) if len(self.memory) else 0.0,
                "assoc_norm": float(np.linalg.norm(self.assoc)) if self.assoc.size else 0.0,
            },
            metadata={"backend": "memory", "mode": self.config.mode, "numba": self.use_numba},
            valid=True,
        )


__all__ = [
    "MemoryConfig",
    "MemoryTrace",
    "RecallResult",
    "PleromaMemoryEngine",
]
