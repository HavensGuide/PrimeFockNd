import numpy as np
from causal_defect_engine import DefectGeometryCompiler, CompilerConfig, NUMBA_AVAILABLE

rng = np.random.default_rng(0)
Omega = np.eye(10)
for i in range(9):
    Omega[i, i+1] = Omega[i+1, i] = 0.6 + 0.02*i
Omega[0, 2] = Omega[2, 0] = -0.7
Omega[3, 6] = Omega[6, 3] = 0.5

compiler = DefectGeometryCompiler(CompilerConfig(
    mode="efficient",
    topology_method="mst_chords",
    connection_threshold=0.1,
    curvature_tol=0.01,
    gauge_fix=True,
    use_numba=True,
    max_curvature_defects=64,
))
geom = compiler.fit(Omega)
rep = compiler.report().as_dict()
print("numba", NUMBA_AVAILABLE)
print("summary", geom.summary())
print("report", rep)
assert rep["n_edges"] >= 9
assert rep["n_tree_edges"] >= 9
assert "gauge_fixed" in geom.metadata
print("amp", geom.eval_amplitude([0,1,0,1,0,1,0,1,0,1]))
