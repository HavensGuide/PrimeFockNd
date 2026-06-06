# Rational Defect Resolution and Arithmetic Filling Resistance:
## Abstract Procedures for Arbitrary Rational/Diophantine Solving

**Research notes — companion to *rational_flow_solver.py* and the Stable Defect Resolution framework**
*Draft, June 2026*

---

## Abstract

The current hybrid rational solver — symbolic exact pass, bounded height search, homotopy continuation plus rational reconstruction — is an honest but incomplete attack framework. Its incompleteness has a specific shape: it uses the topology of the solution variety only as a *search* instrument, never as a *certification* or *obstruction* instrument. The Stable Defect Resolution (SDR) framework and the Causal Filling Modulus (CFM) procedure together suggest how to repair this. The organizing move in both is the same: topology classifies a raw space of candidates; a *filling cost* is assigned to topological sectors; filters — spectral, cohomological, anomaly-theoretic — reduce the raw space to an admissible stable spectrum. This note develops six abstract procedures that carry that logic into the arithmetic setting. Topology of the real solution variety, adelic local-global structure, a rational filling resistance modulus, a height-gradient stability flow, and an obstruction-filtered adelic projection together constitute a cascade which, for systems where the Hasse principle holds and global structure is accessible, approaches arbitrariness. For systems where it fails, the framework produces an explicit named obstruction — a Brauer class, a Tamagawa defect — rather than silence.

---

## 1. The Gap: Why Topology Is Only the First Sieve

The SDR paper opens with a precise statement of the problem it solves: given a vacuum manifold $\mathcal{M}$, the raw topological spectrum $\mathscr{D}_\mathrm{raw}(\mathcal{M}) = \bigcup_{q \geq 1} \pi_{q-1}(\mathcal{M})$ is too large. A homotopy class may admit no finite-energy representative, may carry an anomaly, may have a negative fluctuation mode, or may decompose into lower-charge sectors. The physically real spectrum is the *stable, finite, anomaly-free, spectrally non-negative* part.

The same gap appears in Diophantine geometry. Given a system of equations $F(x) = 0$ over $\mathbb{Q}$, the naive search space — all height-bounded rationals, all homotopy tracks — is too large and unstructured. A homotopy path finds numeric candidates without knowing whether the topology of the real solution variety $V(\mathbb{R})$ admits rational points at all. A bounded search exhausts height without knowing whether arithmetic curvature forces rational points into inaccessible regions. The solver currently has no procedure for asking: *what is the admissible stable rational spectrum of this system?*

The procedures below build that filter cascade. The guiding chain is:

$$
\boxed{F = 0 \;\longrightarrow\; V_\mathbb{C} \;\longrightarrow\; V_\mathbb{R} \;\longrightarrow\; \text{adelic sections} \;\longrightarrow\; \text{height-stable fixed points} \;\longrightarrow\; V(\mathbb{Q})}
$$

Each arrow is a derived procedure, not a search. The rational points are the *indecomposable, obstruction-free, height-stable defects* of the arithmetic vacuum defined by $F$.

---

## 2. The Arithmetic Vacuum Datum

**Definition 2.1** (Arithmetic vacuum datum). An *arithmetic vacuum datum* for a system $F : \mathbb{A}^n \to \mathbb{A}^m$ of polynomial equations is the tuple

$$\mathfrak{V}_\mathrm{arith} = (V_\mathbb{C},\; V_\mathbb{R},\; \{V(\mathbb{Q}_p)\}_{p \leq \infty},\; \pi_\mathrm{top},\; \kappa_\mathbb{Q},\; \mathscr{D}_\mathrm{obs})$$

where $V_\mathbb{C}$ is the complex variety (the complete solution locus, accessible by Bézout and numerical algebraic geometry), $V_\mathbb{R}$ is its real slice, $\{V(\mathbb{Q}_p)\}$ is the family of local solution spaces at each prime $p$ and at the archimedean place $p = \infty$, $\pi_\mathrm{top}$ is the homotopy/cohomology functor applied to $V_\mathbb{R}$, $\kappa_\mathbb{Q}$ is the rational filling resistance modulus defined in Section 4, and $\mathscr{D}_\mathrm{obs}$ is the sheaf of local obstruction data.

The rational points $V(\mathbb{Q})$ are the *admissible stable defects* of this datum: those classes $\alpha \in \pi_0(V_\mathbb{R})$ which survive every local filter, carry finite height, and have zero arithmetic filling resistance. A raw homotopy candidate that fails any filter is either rejected or assigned an explicit obstruction class.

**Remark 2.2.** This exactly parallels the SDR vacuum datum $\mathfrak{V} = (X, \mathcal{M}, G, H, \mathcal{A}, D, S)$. The base space $X$ corresponds to $\mathrm{Spec}(\mathbb{Z})$, the vacuum manifold $\mathcal{M}$ to $V_\mathbb{C}$, the symmetry group $G$ to the Galois group $\mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q})$, and the spectral operator $D$ to the height function $h$ and its Hessian.

---

## 3. Procedure A — Rational Defect Spectrum (Adelic Sieve)

The SDR paper's first filter is finite energy: a homotopy class must have a representative of finite completed action. In the arithmetic setting, the analogue of finite completed action is *local solvability at every place*.

**Definition 3.1** (Rational defect class). A *rational defect class* is a connected component $\alpha \in \pi_0(V_\mathbb{R})$ such that:

1. **Local solvability.** For every prime $p$, $V(\mathbb{Q}_p) \cap \overline{\alpha} \neq \emptyset$, where $\overline{\alpha}$ denotes the $p$-adic closure of the component.
2. **Archimedean solvability.** $\alpha \neq \emptyset$ as a subset of $V(\mathbb{R})$.

These two conditions together constitute the *adelic solvability filter*.

**Definition 3.2** (Raw rational defect spectrum). The *raw rational defect spectrum* is

$$\mathscr{D}_\mathrm{raw}(\mathfrak{V}_\mathrm{arith}) = \pi_0(V_\mathbb{R}).$$

**Definition 3.3** (Adelic defect spectrum). The *adelic defect spectrum* is

$$\mathscr{D}_\mathrm{adel}(\mathfrak{V}_\mathrm{arith}) = \left\{ \alpha \in \pi_0(V_\mathbb{R}) : V(\mathbb{Q}_p) \cap \bar\alpha \neq \emptyset \;\forall p \leq \infty \right\}.$$

**Procedure A.** Given $F$, enumerate $\pi_0(V_\mathbb{R})$ by topological methods (interval arithmetic, Morse theory applied to $\|F\|^2$). For each component, check local solvability at small primes by Hensel lifting from mod-$p$ solutions. Discard components failing any local test. Report each discarded component with the prime $p$ at which solvability fails — this is the *topological charge vanishing* analogue.

**Remark 3.4** (Hasse principle). For many classical systems — quadrics, smooth cubic curves with a rational point, varieties satisfying the conditions of the Hasse-Minkowski theorem — adelic solvability is *sufficient* for rational solvability. In these cases Procedure A is already decisive. The procedures below address the gap when it is not.

---

## 4. Procedure B — Rational Filling Resistance Modulus

The causal filling modulus $\kappa_\mathrm{caus}(\Gamma) = S_\mathrm{fill}(\Gamma)/\tau_\Gamma^2$ assigns a stiffness to each causal cycle by minimizing a Fisher information cost over all admissible fillings. Low modulus means the region is geometrically soft — filling is cheap. High modulus signals curvature concentration.

The arithmetic analogue assigns a *height cost* to filling a loop in $V_\mathbb{R}$ with rational data.

**Definition 4.1** (Rational filling resistance). Let $\Gamma \subset V_\mathbb{R}$ be a closed loop and let $\tau_\Gamma = \max_{x,y \in \Gamma} d_\mathbb{R}(x,y)$ be its real diameter. The *rational filling resistance* of $\Gamma$ is

$$\kappa_\mathbb{Q}(\Gamma) = \inf_{\substack{\Sigma : \partial\Sigma = \Gamma \\ \Sigma \cap V_\mathbb{Q} \text{ dense in } \Sigma}} \frac{h_\mathrm{Weil}(\Sigma)}{\tau_\Gamma^2}$$

where $h_\mathrm{Weil}(\Sigma) = \inf_{x \in \Sigma \cap V(\mathbb{Q})} h(x)$ is the minimal Weil height of any rational point on the filling surface $\Sigma$, and the infimum runs over all surfaces bounded by $\Gamma$ whose rational points are dense.

**Properties.** Three properties follow immediately, in precise parallel to the causal filling modulus.

$(\mathrm{i})$ $\kappa_\mathbb{Q}(\Gamma) \geq 0$ since both numerator and denominator are non-negative.

$(\mathrm{ii})$ $\kappa_\mathbb{Q}(\Gamma) = 0$ if and only if $\Gamma$ bounds a rational surface containing arbitrarily low-height points — the arithmetic geometry is *flat* inside $\Gamma$.

$(\mathrm{iii})$ Normalization by $\tau_\Gamma^2$ removes the trivial scale dependence and makes $\kappa_\mathbb{Q}$ sensitive to the *shape* of the arithmetic curvature rather than the accidental size of the loop.

**Interpretation.** A region of $V_\mathbb{R}$ with low rational filling resistance is arithmetically soft: rational points exist and are accessible at small height. A region with high resistance is arithmetically stiff: rational points, if they exist, are forced to large height or may be absent. The resistance modulus replaces the ad hoc `max_height` parameter in the current solver with a geometrically motivated bound: search is concentrated in low-resistance components; high-resistance components are flagged as *arithmetically hard* sectors requiring either a Diophantine certificate of sparseness or a deeper obstruction analysis.

**Definition 4.2** (Arithmetic stiffness spectrum). The *arithmetic stiffness spectrum* of $\mathfrak{V}_\mathrm{arith}$ is the assignment

$$\alpha \mapsto \kappa_\mathbb{Q}(\alpha) = \inf_{\Gamma \subset \alpha} \kappa_\mathbb{Q}(\Gamma)$$

taken over representative loops generating $H_1(\alpha; \mathbb{Z})$.

Components with $\kappa_\mathbb{Q}(\alpha) = 0$ are *arithmetically flat* (abundant rational points, low height). Components with $\kappa_\mathbb{Q}(\alpha) = \infty$ are *arithmetic insulators* — this is the arithmetic analogue of infinite filling modulus indicating a non-physical defect sector.

---

## 5. Procedure C — Brauer–Manin Obstruction Sheaf

When the Hasse principle fails — a system is locally solvable everywhere but has no rational point — the defect is not topological but cohomological. The SDR paper captures the analogous failure through the positivity-sheaf completion analysis: local positivity cones may fail to glue into a global positive cone, and the gluing failure is an obstruction class $\Omega_\alpha \in H^1(X, \mathrm{Aut}(\mathcal{P}_\alpha))$.

The arithmetic analogue is the Brauer–Manin obstruction.

**Definition 5.1** (Positivity sheaf of local solution cones). For each open $U \subset \mathrm{Spec}(\mathbb{Z})$, define the *local solution cone*

$$\mathcal{S}_\alpha(U) = \left\{ x \in V(\mathbb{Q}_p) : p \in U,\; x \in \bar\alpha \right\}.$$

The assignment $U \mapsto \mathcal{S}_\alpha(U)$ forms a sheaf of local rational candidates over the arithmetic base $\mathrm{Spec}(\mathbb{Z})$.

**Definition 5.2** (Adelic admissibility projection). Let $\mathbb{A}_\mathbb{Q} = \prod_p' \mathbb{Q}_p \times \mathbb{R}$ be the adele ring. The *adelic admissibility projection* is

$$P_\mathcal{A} : \prod_{p \leq \infty} V(\mathbb{Q}_p) \twoheadrightarrow \mathscr{D}_\mathrm{adel},$$

the natural map from the full product of local solution sets to the adelic defect spectrum.

**Definition 5.3** (Brauer obstruction class). For a component $\alpha \in \mathscr{D}_\mathrm{adel}$, the *Brauer obstruction class* is

$$\Omega_\alpha \in \mathrm{Br}(V) = H^2_\mathrm{ét}(V, \mathbb{G}_m)$$

computed from the Azumaya algebra obstruction to extending local sections to a global rational point. Vanishing $\Omega_\alpha = 0$ is necessary for $\alpha$ to contain a rational point; it is equivalent to

$$V(\mathbb{A}_\mathbb{Q})^{\mathrm{Br}} \neq \emptyset \implies V(\mathbb{Q}) \neq \emptyset$$

in the cases where the Brauer–Manin obstruction is the only one.

**Procedure C.** For each component surviving Procedure A, compute finite-dimensional approximations to $\mathrm{Br}(V)/\mathrm{Br}(\mathbb{Q})$ using the étale cohomology of the variety. For many curves and surfaces of interest this is accessible via SageMath's algebraic geometry modules. Report $\Omega_\alpha$ explicitly when non-zero. This is the *anomaly cancellation* filter in the SDR language: a component with $\Omega_\alpha \neq 0$ is anomalous and removed from the admissible spectrum.

**Remark 5.4.** This is the cohomological analogue of SDR §5.1: the gluing obstruction $\Omega_\alpha \in H^1(X, \mathrm{Aut}(\mathcal{P}_\alpha))$ for positivity sheaves is replaced by $\Omega_\alpha \in \mathrm{Br}(V)$ for rational section sheaves. In both cases, vanishing of the obstruction is necessary but not sufficient for the existence of a global object.

---

## 6. Procedure D — Height-Flow Stability and the Rational Morse Index

The SDR paper's most powerful procedure is Stable Defect Resolution itself: instead of checking $\mathrm{Spec}(\mathcal{H}_\alpha) \geq 0$ directly, a path $\mathcal{H}_\alpha(t) = \mathcal{H}_\mathrm{arch} - t\mathcal{H}_\mathrm{def}$ is deformed from a known-positive operator to the target Hessian, and the Morse index $\mu(\alpha) = \#\{\lambda < 0 : \lambda \in \mathrm{Spec}(\mathcal{H}_\alpha)\}$ is tracked. Instability is born only at zero crossings. A stable sector has $\mu(\alpha) = 0$ and no zero crossings on $t \in [0,1]$.

The arithmetic analogue introduces a *height-gradient deformation* of the equation system.

**Definition 6.1** (Height flow). Let $h : V_\mathbb{R} \to \mathbb{R}_{\geq 0}$ be the real extension of the logarithmic Weil height. The *height flow* on the real variety is the one-parameter family of deformed systems

$$F_t(x) = F(x) + t \cdot \nabla_x h(x), \qquad t \in [0, \infty).$$

At $t = 0$ this is the original system. As $t \to \infty$, the height gradient dominates and pushes solutions toward height-minimizing configurations.

**Definition 6.2** (Rational fixed point). A point $x^* \in V_\mathbb{R}$ is a *rational fixed point of the height flow* if $x^* \in V(\mathbb{Q})$ and $x^*$ is a local minimum of $h|_{V_\mathbb{R}}$.

Rational points are precisely the points where $F(x^*) = 0$ and $h$ has a local minimum on the variety — they are not moved by the height flow because they already occupy the height-minimal configuration in their neighborhood.

**Definition 6.3** (Rational Morse index). The *rational Morse index* of a component $\alpha$ is

$$\mu_\mathbb{Q}(\alpha) = \#\left\{ \lambda < 0 : \lambda \in \mathrm{Spec}\left(\nabla^2 h|_{V_\mathbb{R}} \big|_\alpha\right) \right\}$$

counted with multiplicity after projecting out directions tangent to $V_\mathbb{R}$. A component with $\mu_\mathbb{Q}(\alpha) = 0$ is *height-stable*: the height function restricted to $\alpha$ has no negative curvature directions, and rational points in $\alpha$ are genuine minima rather than saddle-point artifacts.

**Procedure D.** For each component surviving Procedures A–C, integrate the height flow $F_t$ from $t = 0$ to $t = T$ for moderate $T$. Track zero crossings of the height stability determinant

$$D(t) = \det_\mathrm{reg}\left(I - t \cdot B_\alpha\right), \qquad B_\alpha = \mathcal{H}_h^{-1/2} \mathcal{H}_F \mathcal{H}_h^{-1/2}$$

where $\mathcal{H}_h$ is the Hessian of $h$ on $V_\mathbb{R}$ and $\mathcal{H}_F$ is the Hessian of $\|F\|^2$. A component with no zeros of $D(t)$ on $[0,1]$ is *height-stable* and is admitted to the rational defect spectrum. Height-unstable components are flagged as possessing non-minimal arithmetic configurations — their homotopy candidates are saddle points in height space and are unlikely to be rational.

**Remark 6.4.** The no-crossing theorem in the SDR heat-flow model (Proposition 6.1 of SDR) states that if all zeros of $H_{t_0}$ are real and simple, the real-zero configuration is forward-stable under heat time. The arithmetic analogue: if all height-flow fixed points in $\alpha$ at $t = 0$ are non-degenerate rational minima, then the rational configuration is forward-stable under height deformation and no spurious rational candidates are born.

---

## 7. Procedure E — Obstruction-Filtered Rational Analysis

SDR §8 defines obstruction-filtered stable defect analysis: a kernel $K$ may fail positivity on the raw test space while succeeding on the *admissible* defect subspace. The physically relevant condition is $P_\mathcal{A} K P_\mathcal{A} \geq 0$, not $K \geq 0$ globally. The diagnostic procedure is: find a negative mode, determine its admissibility sector, project to the admissible subspace, test whether the obstruction survives.

This structure applies directly to rational candidates from homotopy continuation.

**Definition 7.1** (Admissible adelic subspace). For a component $\alpha \in \mathscr{D}_\mathrm{adel}$, the *admissible adelic subspace* is

$$\mathcal{A}_\alpha = \left\{ (x_p)_{p \leq \infty} \in \prod_{p \leq \infty} V(\mathbb{Q}_p) : x_p \in \bar\alpha,\; \sum_p \log|F(x_p)|_p \leq \varepsilon \right\},$$

the set of adelic points near $\alpha$ satisfying a uniform residual bound in each local norm.

**Definition 7.2** (Obstruction-filtered rational analysis). Given a homotopy candidate $\hat x \in V_\mathbb{R}$ with small real residual $\|F(\hat x)\|$:

1. Compute the local reduction $\hat x \bmod p$ at each small prime.
2. Attempt Hensel lifting to obtain $\hat x_p \in V(\mathbb{Q}_p)$.
3. Form the adelic candidate $(\hat x_p)_{p \leq \infty} \in \mathcal{A}_\alpha$.
4. Project onto $P_{\mathcal{A}_\alpha}$ — the projection onto those adelic points admissible to $\alpha$.
5. Compute the residual of the projected candidate and test whether it survives the adelic admissibility filter.

If the projection kills the residual — the obstruction was non-physical — the candidate is admitted. If the obstruction survives projection, it is a genuine local-global defect, and the candidate is recorded as an *obstruction-type negative diagnostic* with the prime(s) at which the Hensel lift fails.

**Remark 7.3.** This is the arithmetic analogue of SDR §9's $\mathrm{PF}_5$ diagnostic for the de Bruijn–Newman kernel. There, the raw negative minor of the Toeplitz matrix survived the even/cosine admissibility projection, establishing that total positivity is not available as a global certificate. Here, a homotopy candidate whose Hensel lift fails at some prime $p$ and whose failure survives the adelic projection constitutes a *certified non-rational* negative diagnostic — it names the obstruction rather than reporting only `UNKNOWN`.

---

## 8. Procedure F — Arithmetic Filling Recovery

The SDR paper develops several *recovery* procedures — prime bridge recovery, gamma/archimedean recovery, scattering recovery — each of which reconstructs a global operator from its local components by exploiting the structure of the completed filling form $\mathcal{H}_\mathrm{fill}^\zeta = \mathcal{H}_\infty - \mathcal{H}_\mathrm{prime}$. The cone-contractivity criterion then says the global positivity question reduces to $\mathcal{B}_\zeta = \mathcal{H}_\infty^{-1/2} \mathcal{H}_\mathrm{prime} \mathcal{H}_\infty^{-1/2} \leq I$.

The arithmetic analogue for rational solving is *height recovery from local data*: reconstruct the global height function from its local $p$-adic and archimedean components, and use the cone-contractivity condition to certify whether rational points exist.

**Definition 8.1** (Completed height filling form). Let

$$h_\mathrm{fill}(x) = h_\infty(x) - \sum_p h_p(x)$$

where $h_\infty(x) = \log \max(1, |x|_\infty)$ is the archimedean height contribution and $h_p(x) = -\log |x|_p$ is the $p$-adic contribution for each finite prime. By the product formula, $h_\mathrm{fill}(x) = h_\mathrm{Weil}(x)$ for $x \in \mathbb{Q}^n$.

The height Hessian splits as

$$\mathcal{H}_h = \mathcal{H}_\infty - \mathcal{H}_p^\mathrm{fin}$$

in exact parallel to $\mathcal{H}_\mathrm{fill}^\zeta = \mathcal{H}_\infty - \mathcal{H}_\mathrm{prime}$.

**Procedure F.** For a system $F$ with suspected rational solution in component $\alpha$:

1. Compute the archimedean height Hessian $\mathcal{H}_\infty$ on $V_\mathbb{R}$ near $\alpha$.
2. For each prime $p$ below a bound, compute the local height Hessian $\mathcal{H}_p$ on $V(\mathbb{Q}_p) \cap \bar\alpha$.
3. Assemble the finite-defect transfer operator $B = \mathcal{H}_\infty^{-1/2} \mathcal{H}_p^\mathrm{fin} \mathcal{H}_\infty^{-1/2}$.
4. Compute the regularized filling determinant $D_\mathrm{arith}(t) = \det_\mathrm{reg}(I - tB)$.
5. If $D_\mathrm{arith}(t) > 0$ for all $t \in [0,1]$, the height transfer is *contractive in the archimedean completion norm*: rational points in $\alpha$ are arithmetically stable and the component is admitted.

**Theorem 8.2** (Arithmetic cone-contractivity). Under the assumption $\mathcal{H}_\infty > 0$ on $\mathcal{A}_\alpha$,

$$\mathcal{H}_\infty - \mathcal{H}_p^\mathrm{fin} \geq 0 \iff B \leq I.$$

Equivalently, finite local height transfer must be a *contraction in the archimedean completion norm*. A component satisfying this is height-stable; a component violating it has a local prime at which rational points are pushed to higher height than the archimedean completion permits — a concrete arithmetic obstruction.

---

## 9. The Completed Solver Architecture

The six procedures form a cascade exactly parallel to the SDR filter sequence. The combined solver architecture is:

```
Input: system F(x) = 0 over Q

Stage 1  [Symbolic/Gröbner]          → exact rational solutions, algebraic non-rationals
Stage 2  [Bounded height search]      → rational candidates up to height H
Stage 3  [Homotopy continuation]      → numeric candidates + rational reconstruction
─────────────────────────────────────────────────────────────────────────
Stage 4A [Topology of V_R]            → π₀(V_R), raw rational defect spectrum
Stage 4B [Adelic sieve]               → D_adel, filter by local solvability at all primes
Stage 5  [Rational filling resistance]→ κ_Q(α) per component, concentrate search
Stage 6  [Brauer–Manin obstruction]   → Ω_α ∈ Br(V), filter anomalous components
Stage 7  [Height flow stability]      → μ_Q(α), Morse index filter, D(t) no-crossing
Stage 8  [Obstruction-filtered proj.] → Hensel lift + adelic projection, name survivors
Stage 9  [Arithmetic filling recovery]→ B ≤ I check, height contraction certificate
─────────────────────────────────────────────────────────────────────────
Output:  Admissible rational defect spectrum with explicit obstruction data
```

The new status vocabulary extends the current `SolverStatus` enum:

| Status | Meaning |
|---|---|
| `FOUND_RATIONAL` | Certified rational solution |
| `FOUND_INTEGER` | Certified integer solution |
| `FOUND_ALGEBRAIC_NONRATIONAL` | Exact non-rational algebraic solution |
| `ADELIC_OBSTRUCTION` | Locally solvable everywhere, rational point absent; Brauer class named |
| `HEIGHT_UNSTABLE` | Component found but height Morse index > 0; candidates are saddle points |
| `ARITHMETIC_STIFF` | $\kappa_\mathbb{Q}(\alpha) > \kappa_\mathrm{thresh}$; rational points forced to large height or absent |
| `NO_SOLUTION_WITHIN_BOUND` | Bounded search exhausted; no rational solution up to height $H$ |
| `CERTIFIED_NONE_MODULAR` | Integer obstruction: no solution modulo some prime $p$ |
| `UNKNOWN` | No procedure was decisive |

---

## 10. What New Abstract Objects Are Required

The procedures above require three genuinely new mathematical objects that are not present in the current solver.

**The rational filling resistance modulus** $\kappa_\mathbb{Q}(\Gamma)$ is the central new invariant. It plays the role that $\kappa_\mathrm{caus}(\Gamma) = S_\mathrm{fill}/\tau_\Gamma^2$ plays in the causal filling programme and that the hemisphere escape modulus $\kappa(D) = \mathcal{S}_\mathrm{esc}(D)/\lambda_D^{n+1}$ plays in the Fisher-frame spinorial defect programme. All three assign a stiffness or cost to a topological structure by a minimax variational principle, and all three produce coupling constants or height bounds as inverse stiffnesses. The rational filling resistance requires a theory of *height-cost minimization over constrained surfaces in algebraic varieties*, which sits at the intersection of Arakelov geometry and variational Diophantine approximation. This is the core invention needed.

**The arithmetic Morse index** $\mu_\mathbb{Q}(\alpha)$ is new. Ordinary Morse theory applied to smooth manifolds counts negative eigenvalues of a Hessian along a path. The arithmetic version must be applied to the height function restricted to a real algebraic variety, and it must be computable from the Jacobian of $F$ and the gradient of $h$ without requiring a smooth ambient manifold. A correct definition requires the notion of a *height-critical point on a singular variety*, which is accessible via the theory of Lagrange multipliers on algebraic sets but has not been systematically developed for Diophantine applications.

**The height filling determinant** $D_\mathrm{arith}(t) = \det_\mathrm{reg}(I - tB)$ with $B = \mathcal{H}_\infty^{-1/2} \mathcal{H}_p^\mathrm{fin} \mathcal{H}_\infty^{-1/2}$ is the arithmetic regularized stability determinant. Its zero crossings correspond to primes at which local height transfer ceases to be contractive — i.e., where local arithmetic geometry conflicts with global rational structure. Computing this determinant requires a regularization scheme for the infinite product over primes, analogous to the zeta-function regularization used for spectral determinants in the SDR pleroma context. The natural candidate is a height-zeta regularization:

$$D_\mathrm{arith}(t) = \exp\left( -\frac{d}{ds}\bigg|_{s=0} \zeta_{B,t}(s) \right), \qquad \zeta_{B,t}(s) = \sum_{\lambda \in \mathrm{Spec}(I - tB)} \lambda^{-s}.$$

This connects the arithmetic filling recovery procedure to the spectral theory of the completed height form, and in particular to the question of whether the adelic height zeta function has analytic continuation properties that force the determinant to be non-vanishing on $[0,1]$ — the arithmetic analogue of the no-double-zero theorem for the de Bruijn–Newman heat flow.

---

## 11. Relation to the Broader Programme

The three filling-cost objects across the programme share a common abstract structure:

| Domain | Object | Topology | Filling Cost | Stiffness | Coupling |
|---|---|---|---|---|---|
| Gravity | $\kappa_\mathrm{caus}(\Gamma)$ | Causal cycles $\Gamma$ | Fisher info $\int \mathcal{I}[\rho, K]$ | $\kappa_\mathrm{grav} = G^{-1}$ | Newton's $G$ |
| Spinorial matter | $\kappa(D)$ | Hemisphere defect $D$ | Escape action $\mathcal{S}_\mathrm{esc}$ | BPS condition | $g_0^2 = 3/4\pi$ |
| Arithmetic | $\kappa_\mathbb{Q}(\Gamma)$ | Loops in $V_\mathbb{R}$ | Weil height infimum | $\kappa_\mathbb{Q}^{-1}$ | Height bound, $V(\mathbb{Q})$ |

In each case: topology classifies the raw object; a filling variational problem assigns a cost; the cost determines the admissible stable spectrum and the effective coupling. The rational defect framework is the third member of this family.

The Brauer–Manin obstruction sheaf (Procedure C) is the arithmetic analogue of the positivity-sheaf obstruction $\Omega_\alpha \in H^1(X, \mathrm{Aut}(\mathcal{P}_\alpha))$. Both are cohomological obstructions to gluing local data into global objects — one for positive quadratic forms, one for rational sections.

The height filling determinant (Procedure F) is the arithmetic analogue of the completed pleroma Hessian $\mathcal{H}_\mathrm{fill}^\zeta = \mathcal{H}_\infty - \mathcal{H}_\mathrm{prime}$ and its cone-contractivity condition $\mathcal{B}_\zeta \leq I$. The no-zero-crossing of $D_\mathrm{arith}(t)$ on $[0,1]$ is the arithmetic analogue of the no-double-zero theorem for the de Bruijn–Newman heat flow — and both may ultimately require the same global mechanism to establish: a self-adjoint realization of the completed operator in a Hilbert space that forces spectral non-negativity.

---

## 12. Open Obligations

The framework as stated leaves several obligations open.

**Computable rational filling resistance.** The modulus $\kappa_\mathbb{Q}(\Gamma)$ as defined involves an infimum over all rational-dense surfaces bounded by $\Gamma$. A computable approximation is needed. The natural approach is via Arakelov height pairing on arithmetic surfaces: for each finite triangulation of $\Sigma$, the height is a finite sum of local terms, and the infimum over triangulations is a combinatorial optimization over the space of rational paths. This is computable in finite cases but a general algorithm is open.

**Hensel lifting radius.** Procedure E requires Hensel lifting of $\hat x \bmod p$ to $\hat x_p \in V(\mathbb{Q}_p)$. For singular fibers of $F \bmod p$, the standard Hensel lemma does not apply. A desingularization step — or a generalized Newton–Puiseux expansion over $\mathbb{Q}_p$ — is required. This is the arithmetic analogue of the completed representative construction in SDR §4 (construct $\mathcal{H}_\alpha$ at a *completed* representative, not an arbitrary one).

**Height zeta regularization.** The regularized height filling determinant $D_\mathrm{arith}(t)$ requires a theory of spectral determinants for the operator $B = \mathcal{H}_\infty^{-1/2}\mathcal{H}_p^\mathrm{fin}\mathcal{H}_\infty^{-1/2}$ acting on the $L^2$ space of the adelic height form. This connects to the analytic theory of height zeta functions and, in particular, to whether the no-crossing property $D_\mathrm{arith}(t) \neq 0$ for $t \in [0,1]$ is forced by global arithmetic structure — the deepest open question in this programme.

**Absolute height bound from resistance.** The filling procedure should ultimately produce a *derived* height bound $H_\alpha = \kappa_\mathbb{Q}(\alpha)^{-1/2}$ for each component, replacing the user-supplied `max_height` parameter with a geometrically motivated one. This requires the full Arakelov-theoretic development of $\kappa_\mathbb{Q}$ and is the counterpart to the derivation of $G = 2\ell_0^2$ from the causal kernel in the filling modulus programme: a coupling constant produced as an inverse stiffness rather than inserted by hand.

---

*This note is a companion to* rational_flow_solver.py *and the Stable Defect Resolution preprint (Kirby, May 2026). The procedures described are program-level proposals in the spirit of the SDR executive summary §1.1. None constitutes a proof of decidability for arbitrary Diophantine systems; the framework is designed to be honest about what it does and does not decide, in exact parallel to the SDR paper's careful distinction between equivalent reformulations and genuine proof mechanisms.*
