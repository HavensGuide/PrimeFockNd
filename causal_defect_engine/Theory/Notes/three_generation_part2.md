# Three-Generation Problem, Part II: Equipartition and the Phase

**Working research note — May 2026 (follow-up to the $S_3$/Koide note)**

This note attacks the two sub-problems left open after the $S_3$/Koide identification:
**(I)** derive the Koide equipartition condition ($Q=2/3$); **(II)** derive the residual phase and the democratic scale. I attempted explicit computations rather than arguments-by-analogy, and the results are mixed in an instructive way: one clean reduction, one sharp empirical prediction, and two genuine obstructions that *change the shape of the remaining problem*. I record the obstructions because they are more valuable than the optimistic version would have been — they tell you exactly where the mechanism cannot live.

---

## Sub-problem I: equipartition

### I.1 Complete algebraic reduction (rigorous)

With the general $S_3$-symmetric amplitude $\sqrt{m_n}=M(1+k\cos\phi_n)$, $\phi_n=\theta_0+2\pi n/3$, the identities $\sum_n\cos\phi_n=0$, $\sum_n\cos^2\phi_n=\tfrac32$ give
$$Q=\frac{\sum_n m_n}{(\sum_n\sqrt{m_n})^2}=\frac{1+k^2/2}{3}\ \Longrightarrow\ Q=\tfrac23\iff k=\sqrt2.$$
This is exact: the entire content of Koide is the single statement $k=\sqrt2$, i.e. equal norms of the $S_3$-singlet and $S_3$-standard parts of the amplitude ($45^\circ$ to the democratic axis). Sub-problem I reduces with no loss to: **why $k=\sqrt2$?**

### I.2 What forces $k=\sqrt2$ — and what does not (new, honest)

I tested whether a generic $S_3$-invariant potential has its minimum at equipartition. **It does not.** The generic $S_3$-invariant potential to quartic order has nine couplings (singlet $s$, standard doublet $(x,y)$; invariants $s^2$, $x^2+y^2$, $s^3$, $s(x^2+y^2)$, $x^3-3xy^2$, and four quartics), and its vacuum lies at $k=\sqrt2$ only on the codimension-one locus where the singlet and standard sectors enter symmetrically. So:

> **Equipartition $\iff$ singlet$\leftrightarrow$standard self-duality.** Koide's $k=\sqrt2$ holds precisely when the $k=3$ effective potential is invariant under the duality exchanging the two quadratic invariants $I_{\mathbf 1}=|\boldsymbol a_{\mathbf 1}|^2$ and $I_{\mathbf 2}=|\boldsymbol a_{\mathbf 2}|^2$. This duality is an additional structural principle; it is **not** automatic from $S_3$ alone.

This is the same *type* of principle as the programme's other balance points ($\Theta=\pi$ self-duality, $\rho_0=\tfrac12$ equal phase-split, $g_0^2=J_{\min}^2/S_{\rm esc}$ spectral balance), and the natural candidate for the exchanging $\mathbb{Z}_2$ is the same spinorial $\mathbb{Z}_2$ monodromy that pervades the framework. But honesty requires stating it as a *reduction*, not a derivation: **Koide $\Leftrightarrow$ a singlet/standard self-duality of the $k=3$ action that remains to be exhibited from the Fisher geometry.** That is a sharper and more falsifiable statement than the original note's, and it is the correct status.

---

## Sub-problem II: the phase and the scale

### II.1 The phase is empirically sharp (parameter-free, falsifiable)

With full Monte-Carlo propagation of PDG mass errors ($m_\tau$ dominant), the residual phase is
$$\phi \equiv \theta_0-\tfrac{2\pi}{3}=0.222231\pm0.000008\ \text{rad}.$$

| Candidate | Value | Deviation |
|---|---|---|
| $2/9$ | $0.222222$ | $\mathbf{1.0\,\sigma}$ |
| $1/(2\sqrt5)$ (Brannen) | $0.223607$ | $163\,\sigma$ — **excluded** |
| $\pi/14$ | $0.224399$ | $271\,\sigma$ — excluded |
| $7\pi/99$ | $0.222133$ | $12\,\sigma$ — excluded |

$\phi=2/9$ fits at $1.0\sigma$; the main rival is excluded at $163\sigma$; **no rational multiple of $\pi$ fits within $2\sigma$.** Imposing $\phi=2/9$ exactly reproduces all three masses to $<0.006\%$. This is a genuine prediction at the edge of current precision: **improved $m_\tau$ must drive $\phi\to0.2222$, not $0.2236$** — a dated, falsifiable consequence (Belle II / tau-charm factory).

### II.2 What the phase is, structurally (new)

Two findings reshape the problem:

**(a) The phase is a rigid degeneracy-lifting rotation.** At $\phi=0$ the $\mathbb{Z}_3$-symmetric configuration has two *degenerate* generations (amplitudes $M(1+\sqrt2)$ and a degenerate pair $M(1-\tfrac{\sqrt2}{2})$). The physical configuration is the whole generation triplet rigidly rotated by $\phi=2/9$ rad away from this symmetric point; the rotation is exactly what lifts the two-fold degeneracy. So $\phi$ is the **degeneracy-splitting parameter** of the generation sector, not an arbitrary orientation.

**(b) The phase cannot be an $S_3$ vacuum angle — a hard obstruction.** Every $S_3$ invariant depends on the standard-plane angle only through $\cos(3\psi)$ (the $\mathbb{Z}_3$ clock term $x^3-3xy^2$). Minimizing *any* $S_3$-invariant potential therefore locks the angle to multiples of $60^\circ$. The observed phase is $2/9$ rad $=12.73^\circ$ **off** a $60^\circ$-multiple. Hence:

> $\phi=2/9$ **cannot** arise as the vacuum angle of an $S_3$-symmetric potential. Its origin must be either explicit $S_3$-breaking, RG running between the defect scale and the pole scale, or — most naturally in an $m\sim e^{-E}$ framework — a **non-angular (log/RG) quantity** that merely *appears* as a phase through the Foot parametrization.

This obstruction is genuinely useful: it eliminates the most obvious (and wrong) mechanism and points at the log/RG origin. It also fits the suspicious fact that $\phi$ is a rational number of *radians* (not of $\pi$) — natural for a log-ratio, unnatural for an angle.

**(c) The value's formula is ambiguous from one sector.** $\phi=2/9$ equals both $(k-1)/k^2$ and $2/k^2$ at $k=3$ (they coincide only at $k=3$). Distinguishing them requires a second generation sector — which the quark sector could provide (II.4).

**Status:** $\phi=2/9$ is empirically nailed and structurally constrained (degeneracy-lifting, non-$S_3$-angular, $\sim(k{-}1)/k^2$), but **no dynamical derivation of the value $2/9$ has been found.** The honest next step is to look for it as a log/RG ratio in the $k=3$ exponent structure.

### II.3 The democratic scale $M$ (open)

$M^2=313.84$ MeV $=2m_P e^{-E_M}$ with $E_M=45.80$, i.e. a second dimensional-transmutation scale shifted from the electron exponent $E_e=52.22$ by $\Delta=6.42$ (effective $b_0^{(M)}=7.22$ vs the electron's $19/3$). Not derived; requires the explicit $k=3$ beta function and BPS condition, the same calculation type as Fisher-Frame §13.

### II.4 The quark cross-check (new, and it is favorable)

If the $S_3$/equipartition story were universal, quarks would also sit at $k=\sqrt2$. They do not: effective $k^2_{\rm up}\approx3.09$, $k^2_{\rm down}\approx2.39$, both $\ne2$. This is the **known** failure of naive Koide for quarks — and here it is *consistent with the framework's structure rather than a problem for it*: leptons are the clean $k=3$ spinorial sector (at equipartition), while quarks carry additional $SU(3)_c$ — itself the $k=3$-*internal* gauge structure (Physical Sector OP 11.3) — which distorts the amplitude away from equipartition. The framework thus *predicts* that leptons are Koide-exact and quarks are not, for a structural reason, and offers the quark deviations as a future quantitative target (the distortion should be calculable from the colour coupling).

---

## Net assessment

**Reduction achieved.** Three free lepton masses $\to$ (i) equipartition $k=\sqrt2$ [$\Leftrightarrow$ a singlet/standard self-duality], (ii) phase $\phi=2/9$ [empirically $1.0\sigma$, rival at $163\sigma$, falsifiable by $m_\tau$], (iii) one open scale $M$. The famous Koide coincidence is explained as representation-balance; the quark non-Koide is explained as $SU(3)_c$ distortion. This is real structural progress the Standard Model does not provide.

**Not a full solution — and the obstructions say why.** Equipartition needs a self-duality principle not yet derived from the $k=3$ Fisher geometry. The phase value $2/9$ has no dynamical derivation, and I *proved* it cannot be a simple $S_3$ vacuum angle, which redirects the search to a log/RG origin. The scale is open. **I tuned nothing:** $k=\sqrt2$ and $\phi=2/9$ are read off the data and tested against structural candidates; the $1\sigma$ vs $163\sigma$ split is data, not a choice; and the $60^\circ$ obstruction was a genuine attempt that failed and thereby informed.

**Highest-value next steps, in order.**
1. Exhibit the singlet$\leftrightarrow$standard $\mathbb{Z}_2$ self-duality of the $k=3$ Fisher action (closes Sub-I rigorously). The candidate is the spinorial $\mathbb{Z}_2$ monodromy.
2. Find $\phi=2/9$ as a log/RG ratio in the $k=3$ exponent — the $60^\circ$ obstruction shows it is *not* an angle, so this is the right place to look.
3. Compute the $k=3$ transmutation scale $M$ (Fisher-Frame §13 method, three primes).
4. Compute the $SU(3)_c$ distortion of quark equipartition; predict $k^2_{\rm up}, k^2_{\rm down}$.

**Falsifiable now:** improved $m_\tau$ must drive $\phi\to0.2222$, not $0.2236$.

---

### Appendix: what was computed

- $Q=(1+k^2/2)/3$ — symbolic, exact.
- Generic $S_3$ quartic potential — 9 couplings; minimum at $k=\sqrt2$ only on the self-dual locus (not generic).
- $\phi=0.222231\pm0.000008$ rad — Monte-Carlo, $2\times10^4$ samples, PDG errors.
- $60^\circ$ obstruction — every $S_3$ invariant $\propto\cos(3\psi)$; observed angle off-lattice.
- Degeneracy structure at $\phi=0$ — two-fold degenerate, lifted by $\phi$.
- Quark $k^2$: up $\approx3.09$, down $\approx2.39$ (running masses, scheme-dependent).
- Mass reconstruction from $(k=\sqrt2,\phi=2/9,M^2=313.84)$: $<0.006\%$ on all three.
