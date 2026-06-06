import numpy as np
from causal_defect_engine import (
    DefectState, Defect, TransportEvolution, CurvatureRelaxation,
    DefectDiffusion, ReactionDiffusionEvolution, HamiltonianEvolution,
    EvolutionPipeline, EvolutionConfig, evolve_and_recompress,
)


def main():
    rng = np.random.default_rng(0)
    state = DefectState(
        generators=rng.normal(size=(6, 3)),
        connection={"edges": [(0,1), (1,2), (2,3), (3,4), (4,5)], "weights": np.ones(5)},
        defects=[Defect(kind="charge", charge=1.0, support=(0,)), Defect(kind="charge", charge=-0.5, support=(5,))],
        omega=np.eye(6),
        valid=True,
    )
    pipe = EvolutionPipeline(
        [TransportEvolution(damping=0.01), DefectDiffusion(diffusion=0.1), CurvatureRelaxation(eta=0.05)],
        EvolutionConfig(steps=3, dt=0.1, copy=True),
    )
    state2, report = pipe.evolve(state)
    assert state2.generators.shape == (6, 3)
    assert len(report.energies) == 4

    rd = ReactionDiffusionEvolution(EvolutionConfig(steps=2, dt=0.1), diffusion=0.1, growth=0.2)
    state3, report2 = rd.evolve(state)
    assert "defect_charges" in state3.observables

    psi = np.array([1.0, 0.0], dtype=complex)
    H = np.array([[0, 1], [1, 0]], dtype=float)
    qs = DefectState(generators=psi, observables={"state_vector": psi}, valid=True)
    he = HamiltonianEvolution(EvolutionConfig(steps=1, dt=0.1), H=H, method="exact")
    qs2, _ = he.evolve(qs)
    assert np.isclose(np.linalg.norm(qs2.observables["state_vector"]), 1.0)
    print("v0.8 evolution smoke passed", report.as_dict())

if __name__ == "__main__":
    main()
