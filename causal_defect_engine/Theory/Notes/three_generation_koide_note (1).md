# The Three-Generation Problem in the Fisher-Frame Programme: an $S_3$/Koide Route

**Working research note — May 2026**

Companion to *Fisher-Frame Defects and Spinorial Matter* v0.8 (Open Problem 2) and *Arithmetic Vacuum: Physical Sector* v3.1 (Open Problem 11.3).

---

## 0. Summary

The Fisher-Frame paper derives the electron mass as $m_e = 2m_P\,e^{-32\pi^3/19}$ but leaves the muon and tau unexplained: the naive next hemisphere eigenvalues ($\lambda \in \{6,12\}$) give $m_\mu/m_e \sim \sqrt{3},\ \sqrt{6}$, i.e. $\sim 1.7,\ 2.4$, against the empirical $206.77$ and $3477.2$. The naive eigenvalue ladder fails by two orders of magnitude.

This note argues that the correct structure is **not** a one-particle eigenvalue ladder but a **three-generation $S_3$ multiplet**, and that the empirical content of the lepton spectrum is almost entirely captured by **Koide's relation**

$$Q \;=\; \frac{m_e + m_\mu + m_\tau}{\bigl(\sqrt{m_e}+\sqrt{m_\mu}+\sqrt{m_\tau}\bigr)^2} \;=\; 0.666661 \;\approx\; \frac{2}{3},$$

which holds to within $1.5\times10^{-5}$ of the exact value $2/3$ and is currently unexplained in the Standard Model. The central observation of this note is:

> **Koide's $Q = 2/3$ is the statement that, in the $S_3$ representation decomposition of the amplitude vector $\boldsymbol{a} = (\sqrt{m_e},\sqrt{m_\mu},\sqrt{m_\tau})$, the $S_3$-singlet (democratic) component and the $S_3$-standard (2-dimensional) component have exactly equal magnitude.**

The Fisher-Frame programme has a natural home for this: the three-prime ($k=3$) sector of the Arithmetic Vacuum has abstract automorphism group $\mathrm{Aut}_{\mathrm{abs}} = S_3$ (Physical Sector OP 11.3), and the Fisher-frame identification of mass with squared amplitude, $m \sim |\psi|^2$, makes $\sqrt{m}$ — not $m$ — the natural carrier of the $S_3$ representation. The three charged leptons would then be the three components of an $S_3$ amplitude triplet, and Koide becomes a statement about the singlet/standard balance of that triplet.

This reduces the three-generation problem to two sharp sub-problems, both stated precisely in §4.

---

## 1. Reframing: generations live in the exponent, and the observable is the amplitude

The mass law is exponential, $m = 2m_P e^{-E}$ with $E_e = 32\pi^3/19 = 52.2211$. Writing each generation as a shift,
$$E_\mu = E_e - \ln(m_\mu/m_e),\qquad E_\tau = E_e - \ln(m_\tau/m_e),$$
the empirical shifts are $\delta_\mu = \ln(m_\mu/m_e) = 5.3316$ and $\delta_\tau = \ln(m_\tau/m_e) = 8.1540$. These are $O(10)$ shifts on an exponent of size $\sim 52$, i.e. the generations are *modest perturbations of the same dimensional-transmutation scale*, not wildly different eigenvalues. This already tells us the right object is a near-degenerate multiplet split by a sub-leading mechanism, not three unrelated eigenvalues.

The Fisher-frame defect carries an amplitude $\psi$; its energy (mass) scales as $m \sim |\psi|^2$. Hence the quantity transforming linearly under any internal symmetry acting on amplitudes is $\sqrt{m} \sim |\psi|$. This is the single most important modelling choice in this note, and it is *forced* by the amplitude-squared reading of mass that the gravity paper already uses ($\rho = |\psi|^2$ in the Madelung/symplectic-lift picture). The Koide combination is built from $\sqrt{m}$ precisely because $\sqrt{m}$ is the amplitude.

---

## 2. Koide as $S_3$ representation balance

Let $\boldsymbol{a} = (\sqrt{m_e},\sqrt{m_\mu},\sqrt{m_\tau})$. Decompose under the permutation action of $S_3$ on $\mathbb{R}^3$:
$$\mathbb{R}^3 \;=\; \underbrace{\mathbf{1}}_{\text{trivial (democratic)}} \;\oplus\; \underbrace{\mathbf{2}}_{\text{standard rep}}.$$
The projections are
$$\boldsymbol{a}_{\mathbf 1} = M(1,1,1),\quad M = \tfrac13\textstyle\sum_n\sqrt{m_n},\qquad \boldsymbol{a}_{\mathbf 2} = \boldsymbol{a} - \boldsymbol{a}_{\mathbf 1}.$$
A short computation (verified numerically) gives, with empirical pole masses,
$$|\boldsymbol{a}_{\mathbf 1}| = M\sqrt{3} = 30.6843,\qquad |\boldsymbol{a}_{\mathbf 2}| = 30.6840,$$
so that
$$\frac{|\boldsymbol{a}_{\mathbf 2}|}{|\boldsymbol{a}_{\mathbf 1}|} = 0.999991.$$

**Proposition (Koide = equipartition).** *$Q = 2/3$ if and only if $|\boldsymbol{a}_{\mathbf 1}| = |\boldsymbol{a}_{\mathbf 2}|$, equivalently the angle between $\boldsymbol{a}$ and the democratic axis $(1,1,1)$ is exactly $45^\circ$.*

*Proof.* $\cos^2\theta = |\boldsymbol a_{\mathbf 1}|^2/|\boldsymbol a|^2 = (\sum\sqrt{m_n})^2/(3\sum m_n) = 1/(3Q)$. So $Q = 1/(3\cos^2\theta)$, and $Q = 2/3 \Leftrightarrow \cos^2\theta = 1/2 \Leftrightarrow \theta = 45^\circ \Leftrightarrow |\boldsymbol a_{\mathbf 1}| = |\boldsymbol a_{\mathbf 2}|$. $\square$

The empirical angle is $\theta = 44.9997^\circ$. The lepton amplitude vector bisects the angle between the fully democratic configuration ($\theta = 0$) and the fully traceless/standard configuration ($\theta = 90^\circ$). In representation-theoretic language, **the $S_3$-symmetric and $S_3$-breaking parts of the lepton amplitude are in exact balance.** This is a clean, scale-free, parameter-free statement — exactly the kind of statement the Fisher-Frame programme is built to produce.

---

## 3. What is left after Koide: one phase

Koide fixes one of the two invariants of $\boldsymbol a$ (the singlet/standard balance). The orientation of $\boldsymbol a_{\mathbf 2}$ within the 2-dimensional standard rep is one residual angle. Using Foot's geometric parametrization,
$$\sqrt{m_n} = M\bigl(1 + \sqrt{2}\,\cos(\theta_0 + \tfrac{2\pi n}{3})\bigr),\qquad n = 0,1,2,$$
the equipartition (Koide) condition is the coefficient $\sqrt2$ in front of the cosine — it is *built in*. A single phase $\theta_0$ then reproduces all three masses. A least-squares fit gives
$$\theta_0 = 132.733^\circ = 2.31662\ \text{rad} = 0.73741\,\pi,$$
reproducing $(m_e, m_\mu, m_\tau)$ to $\sim 10^{-3}$ relative accuracy (electron $0.5108$ vs $0.5110$; muon $105.657$ vs $105.658$; tau $1776.88$ vs $1776.86$).

So the entire charged-lepton spectrum is, empirically:
1. an overall scale $M^2 = 313.84$ MeV (the democratic / generation scale);
2. the Koide equipartition condition (the $\sqrt 2$ coefficient), and
3. one phase $\theta_0 = 132.73^\circ$.

The framework must supply all three. Item 2 is the structural heart and is addressed in §2. Items 1 and 3 are the open targets of §4.

---

## 4. The two sub-problems the framework must now solve

**Sub-problem I — Derive Koide equipartition from $S_3$.** Show that the Fisher-frame three-generation amplitude triplet, transforming in the $\mathbf 1 \oplus \mathbf 2$ of the $k=3$ sector's $S_3$, is dynamically forced to the equipartition point $|\boldsymbol a_{\mathbf 1}| = |\boldsymbol a_{\mathbf 2}|$. Candidate mechanisms, in increasing order of ambition:
- *(a) Maximal democratic-symmetry breaking.* The vacuum sits where the $S_3$-breaking order parameter saturates against the $S_3$-symmetric condensate — a balance point enforced by a BPS-like or extremization condition analogous to the $J_{\min}^2/S_{\mathrm{esc}}$ balance that fixes $g_0^2$ in the one-generation sector. The recurring motif in the programme is that physical points are *balance/escape* conditions; equipartition is exactly such a condition.
- *(b) Self-consistent eigenvalue (Sumino-type).* Koide is the statement that the mass matrix is $\propto (\,\mathbb 1 + \text{rank-1 democratic}\,)^2$ with a specific coefficient; a Fisher-frame mass operator whose democratic and traceless parts have equal Hilbert–Schmidt norm would yield it. This needs the explicit $S_3$ mass operator on the $k=3$ sector.

**Sub-problem II — Derive the phase $\theta_0$ and the scale $M$.** The phase $\theta_0 = 132.73^\circ$ and the democratic scale $M^2 = 313.8$ MeV are the two remaining numbers. Note $\theta_0$ is close to but not obviously equal to simple framework angles ($\gamma_{\mathrm{total}}=\pi$, $\Theta=\pi$, Berry $=\pi/2$); a derivation must come from the $S_3$ sector's own phase structure, not the one-generation angles. The scale $M$ should connect to the same UV stiffness $M_*$ that fixes the electron, via a $k=3$-sector dimensional-transmutation analogous to §13 of the Fisher-Frame paper.

---

## 5. Honest caveats and falsification handles

1. **Pole vs running masses.** Koide $Q = 2/3$ holds to $10^{-5}$ for *pole* masses. If the framework predicts a Koide relation at the *defect/RGI* scale rather than the pole scale, the running between scales must preserve the relation to that precision. Whether QED+QCD running preserves $Q = 2/3$ from a high scale down to pole masses is a known subtlety (Sumino introduced family gauge bosons precisely to protect it); the framework needs its own answer. **This is the sharpest internal consistency test.**

2. **Quarks do not obviously obey Koide.** With naive running masses, up-type $Q \approx 0.85$ and down-type $Q \approx 0.73$, neither close to $2/3$. If the $S_3$/amplitude story is fundamental it must either (i) be lepton-specific for a structural reason (e.g. leptons are the "clean" $k=3$ spinorial sector while quarks carry additional $SU(3)_c$ = $k=3$-internal structure that distorts the relation), or (ii) hold for quarks only at a specific common scale. Either way this is a genuine open problem, not to be papered over. The lepton result stands on its own as the cleanest case.

3. **The amplitude ($\sqrt m$) choice is a modelling input.** The whole construction rests on $S_3$ acting on $\sqrt m$, justified by $m \sim |\psi|^2$. This is natural in the Fisher-frame/symplectic-lift picture but should be derived, not assumed, from the $k=3$ sector's representation content.

4. **This note explains a *constraint*, not yet the spectrum.** Koide-as-equipartition reduces three masses to (scale, equipartition, phase) = effectively one structural condition plus two numbers. That is real progress — it converts "three unexplained masses" into "one symmetry condition (Koide) + two numbers (scale, phase)" — but the framework has not yet *derived* any of the three. The deliverable of this note is the **identification of the right structure** (an $S_3$ amplitude triplet at the equipartition point), which is the necessary precondition for a derivation.

---

## 6. Why this is the right track (and what would confirm it)

The naive eigenvalue ladder was the wrong structure: it treated generations as excitations of *one* defect. The $S_3$ reading treats them as the *three components of one multiplet*, which is both more economical and aligned with the arithmetic-vacuum origin ($k=3$ sector, $\mathrm{Aut}_{\mathrm{abs}}=S_3$). Three independent reasons to prefer it:

- It explains the otherwise-mysterious **Koide coincidence** ($Q=2/3$ to $10^{-5}$) as a representation-balance condition rather than an accident — something the Standard Model has never done.
- It uses a structure (**$S_3$ on the $k=3$ sector**) that the programme *already independently posits* for the colour/three-ness sector (Physical Sector OP 11.3), rather than introducing new machinery.
- It is consistent with the **amplitude-squared** reading of mass already used in the gravity sector.

**Confirmation would look like:** a derivation of the equipartition condition from a balance/extremization principle on the $k=3$ Fisher-frame sector (Sub-problem I), followed by a derivation of $\theta_0$ from the $S_3$ phase structure (Sub-problem II). If the same construction that gives equipartition *also* outputs $\theta_0 = 132.73^\circ$ with no free parameter, the three-generation problem is solved and the result is at least as strong as the electron-mass and $G$ results — arguably stronger, because Koide is a famously unexplained precision relation.

**Disconfirmation would look like:** the $k=3$ sector failing to produce equipartition, or producing it only with a tuned coefficient; or the pole-vs-RGI running destroying $Q=2/3$ at the required precision with no structural protection mechanism.

---

## Appendix: numerical summary

| Quantity | Value | Role |
|---|---|---|
| $m_e, m_\mu, m_\tau$ (MeV) | $0.5110,\ 105.658,\ 1776.86$ | inputs (pole) |
| $Q$ | $0.666661$ | Koide; $=2/3$ to $1.5\times10^{-5}$ |
| $\theta$ (amplitude vs democratic axis) | $44.9997^\circ$ | $=45^\circ$: equipartition |
| $|\boldsymbol a_{\mathbf 2}|/|\boldsymbol a_{\mathbf 1}|$ | $0.999991$ | $=1$: singlet/standard balance |
| $\theta_0$ (Foot phase) | $132.733^\circ = 0.7374\,\pi$ | residual phase (open) |
| $M^2$ (democratic scale) | $313.84$ MeV | generation scale (open) |
| $\delta_\mu, \delta_\tau$ (exponent shifts) | $5.3316,\ 8.1540$ | generations as exponent shifts |

All values verified numerically against PDG pole masses.
