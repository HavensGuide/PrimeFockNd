"""
pointcloud_spin.py
==================

Point-cloud spin/connection backend for the Causal Defect Algebra Engine.

Selectable modes
----------------
legacy / exact
    Old all-pairs Gaussian score, full dense metric storage, all-pairs kNN.
local / efficient
    Neighbor-local score/drift, lazy rank-2 tensor storage, kNN output only.
auto
    Uses local/efficient mode for larger clouds and exact mode for small clouds.

The tensor called ``G_flow`` here is intentionally *not* called a Fisher metric.
It is a Madelung/defect frame tensor:

    G_flow = u u^T + alpha v v^T + lambda I

where u = D_eff * grad log rho and v is an empirical current/drift.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .core import (
    EPS,
    CausalDefectEngine,
    Defect,
    DefectState,
    normalize_mode,
    pairwise_sqdist,
    knn_indices,
    winding_loop_numba,
)
from .numba_kernels import (
    local_gaussian_score_knn_numba,
    local_time_drift_knn_numba,
    rank2_frame_numba,
)

try:
    from scipy.spatial import Delaunay  # type: ignore
    SCIPY_AVAILABLE = True
except Exception:  # pragma: no cover
    Delaunay = None
    SCIPY_AVAILABLE = False

try:
    from scipy.sparse import coo_matrix, vstack  # type: ignore
    from scipy.sparse.linalg import lsqr  # type: ignore
    SCIPY_SPARSE_AVAILABLE = True
except Exception:  # pragma: no cover
    coo_matrix = None
    vstack = None
    lsqr = None
    SCIPY_SPARSE_AVAILABLE = False


@dataclass
class PointCloudSpinConfig:
    k_neighbors: int = 15
    n_pca_components: Optional[int] = None
    diffusion_coeff: float = 0.5
    alpha_v: float = 1.0
    lambda_reg: float = 1e-4
    score_bandwidth: Optional[float] = None
    normalize_uv: bool = True
    defect_threshold: float = 0.25
    use_delaunay: bool = True
    mode: str = "auto"          # legacy | exact | local | efficient | auto
    score_mode: str = "auto"    # exact | local | auto
    metric_mode: str = "lazy"   # dense | lazy
    store_dense_metric: bool = False
    max_exact_points: int = 2500


class PointCloudSpinEngine(CausalDefectEngine):
    """Direct point-cloud estimator of a spin/connection defect field."""

    def __init__(self, config: Optional[PointCloudSpinConfig] = None, **kwargs):
        config = config or PointCloudSpinConfig()
        kwargs.setdefault("mode", config.mode)
        super().__init__(**kwargs)
        self.config = config
        self.X_raw: Optional[np.ndarray] = None
        self.t: Optional[np.ndarray] = None
        self.X_work: Optional[np.ndarray] = None
        self.pca_mean_: Optional[np.ndarray] = None
        self.pca_components_: Optional[np.ndarray] = None
        self.score_bandwidth_: Optional[float] = None
        self.nn_: Optional[np.ndarray] = None
        self.nn_d2_: Optional[np.ndarray] = None
        self._resolved_mode: str = "legacy"

    # ------------------------------------------------------------------
    # Fitting / projection
    # ------------------------------------------------------------------

    def fit(self, events: Any, **kwargs) -> "PointCloudSpinEngine":
        if isinstance(events, dict):
            X = np.asarray(events["X"], dtype=np.float64)
            t = None if events.get("t") is None else np.asarray(events["t"], dtype=np.float64)
        elif isinstance(events, tuple) and len(events) == 2:
            X = np.asarray(events[0], dtype=np.float64)
            t = np.asarray(events[1], dtype=np.float64)
        else:
            X = np.asarray(events, dtype=np.float64)
            t = None

        if X.ndim != 2:
            raise ValueError("PointCloudSpinEngine expects X with shape (T, N)")
        if len(X) < 4:
            raise ValueError("Need at least 4 points")
        if t is not None and len(t) != len(X):
            raise ValueError("t must have same length as X")

        self.X_raw = X
        self.t = t
        self.X_work = self._project_fit_transform(X)
        self._resolved_mode = self._choose_mode(len(self.X_work))

        k = min(self.config.k_neighbors, len(self.X_work) - 1)
        self.nn_, self.nn_d2_ = knn_indices(
            self.X_work,
            k,
            backend=self.runtime.neighbor_backend,
            use_numba=self.use_numba,
        )
        self.score_bandwidth_ = self._resolve_bandwidth(self.X_work)

        self.state = DefectState(
            events=self.X_raw,
            metadata={
                "backend": "pointcloud_spin",
                "numba": self.use_numba,
                "mode": self._resolved_mode,
                "score_mode": self._score_mode(),
                "metric_mode": self._metric_mode(),
                "n_points": int(X.shape[0]),
                "n_dim_raw": int(X.shape[1]),
                "n_dim_work": int(self.X_work.shape[1]),
                "score_bandwidth": float(self.score_bandwidth_),
                "has_time": self.t is not None,
                "uses_knn_only": self._score_mode() == "local",
            },
            valid=True,
        )
        return self

    def _choose_mode(self, n: int) -> str:
        m = normalize_mode(self.config.mode or self.runtime.mode)
        if m == "auto":
            return "exact" if n <= self.config.max_exact_points else "local"
        if m == "efficient":
            return "local"
        if m == "legacy":
            return "legacy"
        return m

    def _score_mode(self) -> str:
        sm = normalize_mode(self.config.score_mode)
        if sm == "auto":
            if self._resolved_mode in {"legacy", "exact"}:
                return "exact"
            return "local"
        if sm == "efficient":
            return "local"
        if sm == "legacy":
            return "exact"
        return sm

    def _metric_mode(self) -> str:
        mm = normalize_mode(self.config.metric_mode)
        if self.config.store_dense_metric:
            return "dense"
        if mm in {"legacy", "exact"}:
            return "dense"
        return "lazy"

    def _project_fit_transform(self, X: np.ndarray) -> np.ndarray:
        d = self.config.n_pca_components
        if d is None or d >= X.shape[1]:
            self.pca_mean_ = None
            self.pca_components_ = None
            return X.copy()
        if d < 2:
            raise ValueError("n_pca_components must be at least 2 for spin defects")
        mu = X.mean(axis=0)
        Xc = X - mu
        _, _, Vt = np.linalg.svd(Xc, full_matrices=False)
        self.pca_mean_ = mu
        self.pca_components_ = Vt[:d]
        return Xc @ self.pca_components_.T

    def _project(self, X: np.ndarray) -> np.ndarray:
        if self.pca_components_ is None:
            return X
        return (X - self.pca_mean_) @ self.pca_components_.T

    def _unproject_vector(self, V: np.ndarray) -> np.ndarray:
        if self.pca_components_ is None:
            return V
        return V @ self.pca_components_

    def _resolve_bandwidth(self, X: np.ndarray) -> float:
        if self.config.score_bandwidth is not None:
            return float(self.config.score_bandwidth)
        if self.nn_d2_ is not None and self.nn_d2_.size:
            nn = np.sqrt(np.maximum(self.nn_d2_[:, 0], 0.0))
        else:
            d2 = pairwise_sqdist(X, self.use_numba)
            np.fill_diagonal(d2, np.inf)
            nn = np.sqrt(np.min(d2, axis=1))
        bw = float(np.median(nn) * max(1.0, X.shape[0] ** (-1.0 / (X.shape[1] + 4)) * 4.0))
        return max(bw, 1e-6)

    # ------------------------------------------------------------------
    # Generator estimation
    # ------------------------------------------------------------------

    def estimate_generators(self, **kwargs) -> Dict[str, np.ndarray]:
        if self.X_work is None:
            raise RuntimeError("Call fit first")
        X = self.X_work
        score_mode = self._score_mode()
        if score_mode == "exact":
            score = self._kernel_score_exact(X, X, self.score_bandwidth_)
        elif score_mode == "local":
            score = self._kernel_score_local(X, self.score_bandwidth_)
        else:
            raise ValueError(f"Unknown score_mode: {score_mode}")

        u_work = self.config.diffusion_coeff * score
        v_work = self._estimate_drift_work(X)

        u_metric = u_work.copy()
        v_metric = v_work.copy()
        if self.config.normalize_uv:
            u_scale = np.sqrt(np.mean(np.sum(u_metric * u_metric, axis=1))) + EPS
            v_scale = np.sqrt(np.mean(np.sum(v_metric * v_metric, axis=1))) + EPS
            u_metric /= u_scale
            v_metric /= v_scale
        else:
            u_scale = v_scale = 1.0

        theta, strength, frame = self._theta_from_rank2(u_metric, v_metric)
        u_raw = self._unproject_vector(u_work)
        v_raw = self._unproject_vector(v_work)

        metric_mode = self._metric_mode()
        g_obj: Any
        if metric_mode == "dense":
            g_obj = self._metric_dense(u_metric, v_metric)
        else:
            g_obj = {
                "kind": "rank2_plus_floor",
                "u_metric": u_metric,
                "v_metric": v_metric,
                "alpha_v": self.config.alpha_v,
                "lambda_reg": self.config.lambda_reg,
            }

        gens = {
            "X_work": X,
            "u_work": u_work,
            "v_work": v_work,
            "u": u_raw,
            "v": v_raw,
            "G_flow": g_obj,
            "theta": theta,
            "frame": frame,
            "strength": strength,
            "u_scale": np.asarray(u_scale),
            "v_scale": np.asarray(v_scale),
            "score": score,
            "score_recoverable_from_u": True,
        }
        self.state.generators = gens
        self.state.observables.update(
            {
                "u_mag": np.linalg.norm(u_raw, axis=1),
                "v_mag": np.linalg.norm(v_raw, axis=1),
                "theta": theta,
                "frame_strength": strength,
                "density_recoverable_up_to_constant": True,
            }
        )
        return gens

    def _kernel_score_exact(self, query: np.ndarray, X: np.ndarray, bw: float) -> np.ndarray:
        diff = query[:, None, :] - X[None, :, :]
        dist2 = np.einsum("ijk,ijk->ij", diff, diff)
        logk = -0.5 * dist2 / (bw * bw)
        logk -= logk.max(axis=1, keepdims=True)
        k = np.exp(logk)
        denom = k.sum(axis=1, keepdims=True) + EPS
        return -np.einsum("ij,ijd->id", k / denom, diff) / (bw * bw)

    def _kernel_score_local(self, X: np.ndarray, bw: float) -> np.ndarray:
        if self.nn_ is None:
            raise RuntimeError("Neighbors missing; call fit first")
        if self.use_numba:
            return local_gaussian_score_knn_numba(np.asarray(X, dtype=np.float64), self.nn_.astype(np.int64), float(bw))
        n, d = X.shape
        score = np.zeros((n, d), dtype=np.float64)
        for i in range(n):
            idx = self.nn_[i]
            diff = X[i][None, :] - X[idx]
            dist2 = np.sum(diff * diff, axis=1)
            logk = -0.5 * dist2 / (bw * bw)
            logk -= float(np.max(logk))
            w = np.exp(logk)
            w /= w.sum() + EPS
            score[i] = -(w[:, None] * diff).sum(axis=0) / (bw * bw)
        return score

    def _estimate_drift_work(self, X: np.ndarray) -> np.ndarray:
        if self.t is None:
            return np.zeros_like(X)
        if self.nn_ is None:
            raise RuntimeError("Neighbors missing; call fit first")
        if self.use_numba:
            return local_time_drift_knn_numba(np.asarray(X, dtype=np.float64), np.asarray(self.t, dtype=np.float64), self.nn_.astype(np.int64))
        n, d = X.shape
        v = np.zeros_like(X)
        for i in range(n):
            idx = self.nn_[i]
            tt = self.t[idx]
            XX = X[idx]
            tm = tt.mean()
            tc = tt - tm
            var = float(np.dot(tc, tc)) + EPS
            Xm = XX.mean(axis=0)
            v[i] = (tc[:, None] * (XX - Xm)).sum(axis=0) / var
        return v

    def _metric_dense(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        n, d = u.shape
        eye = np.eye(d)[None, :, :]
        return u[:, :, None] * u[:, None, :] + self.config.alpha_v * v[:, :, None] * v[:, None, :] + self.config.lambda_reg * eye

    def metric_matvec(self, i: int, z: np.ndarray) -> np.ndarray:
        """Lazy G_flow[i] @ z without forming dense metric."""
        if self.state.generators is None:
            raise RuntimeError("Run estimate_generators first")
        G = self.state.generators["G_flow"]
        if isinstance(G, np.ndarray):
            return G[int(i)] @ z
        u = G["u_metric"][int(i)]
        v = G["v_metric"][int(i)]
        return u * float(np.dot(u, z)) + G["alpha_v"] * v * float(np.dot(v, z)) + G["lambda_reg"] * z

    def _theta_from_rank2(self, u: np.ndarray, v: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.use_numba:
            return rank2_frame_numba(np.asarray(u, dtype=np.float64), np.asarray(v, dtype=np.float64), float(self.config.alpha_v))
        n, d = u.shape
        frame = np.zeros((n, d), dtype=np.float64)
        strength = np.zeros(n, dtype=np.float64)
        for i in range(n):
            G = np.outer(u[i], u[i]) + self.config.alpha_v * np.outer(v[i], v[i])
            vals, vecs = np.linalg.eigh(G)
            frame[i] = vecs[:, -1]
            strength[i] = max(vals[-1], 0.0)
            pivot = int(np.argmax(np.abs(frame[i])))
            if frame[i, pivot] < 0:
                frame[i] *= -1.0
        theta = np.arctan2(frame[:, 1], frame[:, 0]) if d >= 2 else np.zeros(n)
        return theta, strength, frame

    # ------------------------------------------------------------------
    # Score / density recovery
    # ------------------------------------------------------------------

    def recover_score(self, ambient: bool = True) -> np.ndarray:
        """Return grad log rho from stored u and D_eff."""
        if self.state.generators is None:
            raise RuntimeError("Run estimate_generators first")
        key = "u" if ambient else "u_work"
        return self.state.generators[key] / (self.config.diffusion_coeff + EPS)

    def recover_log_density(self, k_neighbors: Optional[int] = None, ambient: bool = False) -> np.ndarray:
        """Recover log density on the fitted points up to an additive constant.

        Uses edge integration least squares:
            ell_j - ell_i ~= score_i dot (x_j - x_i)
        This is more direct than a graph-divergence Poisson approximation.
        """
        if self.state.generators is None or self.X_work is None:
            raise RuntimeError("Run estimate_generators first")
        if not SCIPY_SPARSE_AVAILABLE:
            raise ImportError("scipy is required for recover_log_density()")
        X = self.X_work if not ambient else self.X_raw
        score = self.recover_score(ambient=ambient)
        k = k_neighbors or self.config.k_neighbors
        nn, _ = knn_indices(X, k, backend=self.runtime.neighbor_backend, use_numba=self.use_numba)
        rows = []
        cols = []
        vals = []
        b = []
        r = 0
        for i in range(len(X)):
            for j in nn[i]:
                j = int(j)
                rows.extend([r, r])
                cols.extend([j, i])
                vals.extend([1.0, -1.0])
                b.append(float(np.dot(score[i], X[j] - X[i])))
                r += 1
        A = coo_matrix((vals, (rows, cols)), shape=(r, len(X))).tocsr()
        # Anchor mean density log to zero: sum ell_i = 0.
        anchor = coo_matrix((np.ones(len(X)), (np.zeros(len(X)), np.arange(len(X)))), shape=(1, len(X))).tocsr()
        A = vstack([A, 1e-3 * anchor])
        rhs = np.asarray(b + [0.0], dtype=np.float64)
        ell = lsqr(A, rhs, atol=1e-8, btol=1e-8, iter_lim=1000)[0]
        ell -= np.max(ell)  # numerically stable density has max 1 before normalization
        return ell

    # ------------------------------------------------------------------
    # Connection / defects
    # ------------------------------------------------------------------

    def build_connection(self, **kwargs) -> Dict[str, Any]:
        if self.state.generators is None:
            self.estimate_generators(**kwargs)
        X = self.state.generators["X_work"]
        theta = self.state.generators["theta"]
        if self.nn_ is None:
            k = min(self.config.k_neighbors, len(X) - 1)
            self.nn_, self.nn_d2_ = knn_indices(X, k, backend=self.runtime.neighbor_backend, use_numba=self.use_numba)
        edges = []
        transports = []
        for i in range(len(X)):
            for j in self.nn_[i]:
                j = int(j)
                if i < j:
                    delta = (theta[j] - theta[i] + np.pi) % (2 * np.pi) - np.pi
                    edges.append((i, j))
                    transports.append(delta)
        conn = {"edges": np.asarray(edges, dtype=np.int64), "u1_delta": np.asarray(transports), "nn": self.nn_}
        self.state.connection = conn
        return conn

    def detect_defects(self, **kwargs) -> List[Defect]:
        if self.state.connection is None:
            self.build_connection(**kwargs)
        X = self.state.generators["X_work"]
        theta = self.state.generators["theta"]
        strength = self.state.generators["strength"]
        chart = self._chart_2d(X)
        tris = self._triangulate(chart, self.state.connection["nn"])
        theta2 = 2.0 * theta
        defects: List[Defect] = []
        for tri in tris:
            winding2 = float(winding_loop_numba(theta2, np.asarray(tri, dtype=np.int64)))
            charge = 0.5 * winding2
            if abs(charge) > self.config.defect_threshold:
                loc = self.X_raw[list(tri)].mean(axis=0)
                defects.append(
                    Defect(
                        kind="spin_winding",
                        charge=float(np.round(charge * 2) / 2),
                        location=loc,
                        support=tuple(int(x) for x in tri),
                        strength=float(np.min(strength[list(tri)])),
                        metadata={"raw_charge": charge, "winding_2theta": winding2},
                    )
                )
        self.state.defects = self._merge_defects(defects)
        self.state.observables["n_defects"] = len(self.state.defects)
        self.state.observables["triangles_checked"] = int(len(tris))
        return self.state.defects

    def _chart_2d(self, X: np.ndarray) -> np.ndarray:
        if X.shape[1] == 1:
            return np.column_stack([X[:, 0], np.zeros(len(X))])
        if X.shape[1] == 2:
            return X
        Xc = X - X.mean(axis=0)
        _, _, Vt = np.linalg.svd(Xc, full_matrices=False)
        return Xc @ Vt[:2].T

    def _triangulate(self, chart: np.ndarray, nn: np.ndarray) -> np.ndarray:
        if self.config.use_delaunay and SCIPY_AVAILABLE and len(chart) >= 4:
            try:
                return np.asarray(Delaunay(chart).simplices, dtype=np.int64)
            except Exception:
                pass
        tris = set()
        n, k = nn.shape
        neigh = [set(int(x) for x in nn[i]) for i in range(n)]
        for i in range(n):
            nbrs = list(neigh[i])
            for a in range(len(nbrs)):
                for b in range(a + 1, len(nbrs)):
                    j, l = nbrs[a], nbrs[b]
                    if j in neigh[l] or l in neigh[j]:
                        tris.add(tuple(sorted((i, j, l))))
        return np.asarray(sorted(tris), dtype=np.int64) if tris else np.zeros((0, 3), dtype=np.int64)

    def _merge_defects(self, defects: List[Defect]) -> List[Defect]:
        if not defects:
            return []
        locs = np.asarray([d.location for d in defects], dtype=np.float64)
        if len(locs) > 1:
            nn, d2 = knn_indices(locs, 1, backend=self.runtime.neighbor_backend, use_numba=self.use_numba)
            radius = float(np.sqrt(np.median(d2))) if d2.size else 0.0
            full = pairwise_sqdist(locs, self.use_numba) if len(locs) <= 5000 else None
        else:
            radius = 0.0
            full = np.zeros((1, 1))
        radius = max(radius, EPS)
        used = np.zeros(len(defects), dtype=bool)
        merged: List[Defect] = []
        for i, d in enumerate(defects):
            if used[i]:
                continue
            group = [i]
            used[i] = True
            if full is not None:
                candidates = np.where(full[i] <= radius * radius)[0]
            else:
                # local fallback: merge only immediate nearest clusters at large scale
                candidates = [int(x) for x in nn[i]] if len(locs) > 1 else []
            for j in candidates:
                if j != i and not used[j]:
                    group.append(int(j))
                    used[j] = True
            charge = sum(defects[g].charge for g in group)
            loc = np.mean([defects[g].location for g in group], axis=0)
            merged.append(
                Defect(
                    kind="spin_winding_cluster",
                    charge=float(np.round(charge * 2) / 2),
                    location=loc,
                    support=tuple(group),
                    strength=float(min(defects[g].strength for g in group)),
                    metadata={"n_triangles": len(group)},
                )
            )
        return merged

    def holonomy(self, path: Sequence[int]) -> float:
        if self.state.generators is None:
            raise RuntimeError("Run estimate_generators first")
        theta = self.state.generators["theta"]
        total = 0.0
        for a, b in zip(path[:-1], path[1:]):
            total += (theta[int(b)] - theta[int(a)] + np.pi) % (2 * np.pi) - np.pi
        return float(total)
