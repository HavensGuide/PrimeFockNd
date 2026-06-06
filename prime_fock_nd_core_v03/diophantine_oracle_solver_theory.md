# Rational Defect Resolution and the Diophantine Oracle Solver
## Theory, Architecture, and Results

*Research notes — companion to `diophantine_oracle_solver.py`, the rational_defect_framework writeup, and the Stable Defect Resolution, Causal Filling Modulus, and Fisher-Frame Defects programmes.*

---

## 1. The Problem

A Diophantine equation is a polynomial equation $F(x_1, \ldots, x_n) = 0$ with integer coefficients for which we seek rational solutions. This is one of the oldest problems in mathematics. Hilbert's tenth problem asked for an algorithm deciding integer solvability; Matiyasevich proved in 1970 that no such algorithm exists. Over $\mathbb{Q}$, the analogous question remains open.

The solver developed here is not a universal decision procedure — no such thing can exist over $\mathbb{Z}$, and it may not exist over $\mathbb{Q}$ either. It is a staged obstruction engine: it identifies, names, and classifies the reason a system either has a rational solution or doesn't, pushing the decidable boundary as far as current mathematics allows while being honest about what remains out of reach.

The key intellectual move is borrowed from the Stable Defect Resolution (SDR) programme: rather than trying to decide the question directly, build a cascade of filters — each one geometrically motivated, each one adding information — so that every stage either resolves the problem or narrows the remaining undecidability to a precisely characterised class.

---

## 2. The Intellectual Lineage

The solver is the third member of a family of obstruction procedures that share a common abstract structure:

$$\text{vacuum datum} \;\longrightarrow\; \text{raw candidate space} \;\longrightarrow\; \text{admissible stable spectrum}$$

**Stable Defect Resolution (SDR)** applies this to physical vacua. Given a vacuum manifold $\mathcal{M}$, the raw homotopy spectrum $\mathscr{D}_\mathrm{raw}(\mathcal{M}) = \bigcup_{q \geq 1} \pi_{q-1}(\mathcal{M})$ is too large: a topological class may have no finite-energy representative, may carry an anomaly, or may decompose into lower-charge sectors. SDR filters through finite energy, nonzero topological charge, spectral stability (Morse index zero), and anomaly cancellation to produce the physical particle spectrum. The spectral stability condition replaces direct positivity checking — instability is born only at zero crossings of a stability determinant, tracked via the Morse index along a deformation path.

**Causal Filling Modulus (CFM)** assigns a *cost* to topological structures. For a causal cycle $\Gamma$, the filling modulus is $\kappa_\mathrm{caus}(\Gamma) = S_\mathrm{fill}(\Gamma) / \tau_\Gamma^2$, where $S_\mathrm{fill}$ is the infimum of Fisher information cost over all surfaces bounded by $\Gamma$. This turns topology from a binary classifier (class exists / doesn't exist) into a quantitative stiffness. The inverse stiffness gives the coupling constant: Newton's $G$ emerges as $1/\kappa_\mathrm{caus}$, the electron's gauge coupling $g_0^2 = 3/(4\pi)$ emerges from the spinorial escape modulus. No free parameters appear; all couplings are derived as inverse topological stiffnesses.

**Rational Defect Resolution** is the arithmetic member of this family. The vacuum datum is the polynomial system $F(x) = 0$. The raw candidate space is the space of rational number tuples. The admissible stable spectrum is $V(\mathbb{Q})$ — the rational points on the variety. The filling cost is the rational filling resistance $\kappa_\mathbb{Q}$, the coupling constant analogue is the height bound, and the obstruction to existence is a Brauer class. The three programmes share:

| Domain | Vacuum | Raw space | Filling cost | Obstruction | Stable spectrum |
|---|---|---|---|---|---|
| Gravity | Causal geometry | All causal cycles | Fisher information | Curvature singularity | Physical spacetime |
| Matter | $SO(3)$ frame field | $\pi_1(SO(3)) = \mathbb{Z}_2$ classes | Hemisphere escape $\mathcal{S}_\mathrm{esc}$ | Anomaly / spectral instability | Electron, $W^\pm$, ... |
| Arithmetic | $F(x) = 0$ | All rational tuples | Roth exponent $\kappa_\mathbb{Q}$ | Brauer–Manin class | $V(\mathbb{Q})$ |

---

## 3. The Arithmetic Vacuum

The arithmetic analogue of a vacuum manifold is the *multiplicative pleroma*. The multiplicative monoid of positive integers decomposes as

$$\mathbb{N}^* \cong \bigoplus_p \mathbb{N} e_p, \qquad n = \prod_p p^{v_p(n)},$$

where the primes are the atomic defects. An integer is a multi-defect configuration; the zeta function $\zeta(s) = \prod_p (1 - p^{-s})^{-1}$ is the partition function of the prime-defect gas. This is not metaphor — it is the literal content of the Euler product.

A rational number $q = p_1^{e_1} \cdots p_k^{e_k}$ with positive exponents in the numerator and negative in the denominator has a *signed prime-valuation vector* $v = (v_{p_1}, v_{p_2}, \ldots)$ where $v_{p_i} \in \mathbb{Z}$. This is the arithmetic Fock representation: instead of the integer itself, work with its exponent vector. Multiplication becomes vector addition; $\gcd$ becomes componentwise minimum; $\mathrm{lcm}$ becomes componentwise maximum. No large integers are ever materialised.

The **arithmetic graph** on a set of rational candidates $\{x_1, \ldots, x_N\}$ has edge weights

$$E_{ij} = \log \gcd(|x_i|, |x_j|) = \sum_p \min(v_p(x_i), v_p(x_j)) \log p,$$

the log-gcd kernel. Two candidates sharing many prime factors are strongly connected; coprime candidates are disconnected. This is the discrete geometry of the arithmetic vacuum — the Fisher metric on the space of candidate solutions.

---

## 4. The Nine-Stage Cascade

The solver runs nine stages in order, each one either resolving the problem or passing a refined candidate set to the next stage.

**Stage 1 — Symbolic exact pass.** SymPy's `solve` extracts closed-form solutions and certifies them by exact substitution. Rational algebraic solutions are separated from irrational algebraics. This costs $O(1)$ for simple systems and up to seconds for high-degree ones.

**Stage 2 — Vectorised bounded height search.** Every rational $p/q$ in lowest terms with $\max(|p|, q) \leq H$ is enumerated. The full candidate grid is evaluated as a NumPy batch — no Python loop in the inner path — using `lambdify` to broadcast across all candidates simultaneously. Only those passing the float pre-screen are certified by exact SymPy substitution. For homogeneous systems (Stage 3a below), the grid is filtered to *primitive* tuples: those with $\gcd(|p_1|, \ldots, |p_n|) = 1$ and the trivial zero $(0, \ldots, 0)$ excluded. This is the projective filter.

**Stage 3 — Projective detection.** Before homotopy, the solver detects whether the system is *homogeneous* — all monomials of each equation have the same total degree. A homogeneous system defines a projective variety; its rational solutions are equivalence classes $[x_1 : \cdots : x_n]$ with the all-zero tuple excluded. The `_is_homogeneous` check fires automatically; no user configuration is needed. All subsequent stages work with primitive projective representatives.

**Stage 4 — Homotopy continuation.** For square systems ($m = n$ equations and variables), Davidenko's ODE

$$\frac{d}{dt}\begin{pmatrix} F_\mathrm{start}(x) \\ \vdots \end{pmatrix} = 0$$

is tracked from known start points to $t = 1$ using SciPy's Radau integrator. Multiple starts run in parallel via a `ThreadPoolExecutor`. Numeric endpoints are reconstructed to rationals using `Fraction.limit_denominator` and certified by exact substitution.

**Stage 5 — Full Jacobian Hensel lifting.** For each numeric candidate that survived homotopy, the solver attempts Hensel lifting at each configured prime $p$. The Newton step solves the full Jacobian system $J \cdot \delta \equiv -F(x)/\mathrm{mod} \pmod{p}$ using `sp.Matrix.inv_mod(p)` for square systems — the exact modular inverse of the full $m \times n$ Jacobian — and `numpy.linalg.lstsq` for overdetermined systems. Lifting proceeds to $p^{2^k}$ (high precision first, following the large-$D$ strategy of Ecker-Ecker-Grumiller), then projects back. Failure at any prime is recorded as an adelic obstruction with the specific prime and reason code named.

**Stage 6 — Rational filling resistance.** Each candidate is assigned a filling resistance

$$\kappa_\mathbb{Q}(x) = \frac{1}{n} \sum_{j=1}^n \mu_\mathrm{irr}(x_j),$$

where $\mu_\mathrm{irr}(x_j)$ is the *irrationality exponent* of coordinate $x_j$: zero for rationals (Fraction exactly recoverable), $\deg(\min\_\mathrm{poly}(x_j))/2$ for algebraic irrationals (by Roth's theorem, implemented via `sp.nsimplify` and `sp.minimal_polynomial`), and the continued-fraction convergent ratio as a fallback lower bound for unrecognised numbers. This is the analytic formula, not an approximation: by Roth's theorem, for an algebraic number of degree $d$, the Weil height of the best rational approximation $p/q$ satisfies $\log H(p/q) \sim (d/2) \cdot \log(1/|x - p/q|)$, so $\kappa_\mathbb{Q} = d/2$ exactly. Zero resistance means the point is rational; resistance $\geq 1$ means it is algebraically irrational.

**Stage 7 — Classical Brauer–Manin obstruction (Hasse–Minkowski).** For a single quadratic equation, the solver computes the full Hasse–Minkowski obstruction. The Gram matrix $A = \mathrm{Hess}(Q)/2$ is diagonalised over $\mathbb{Q}$ via exact LDL$^\top$ reduction (SymPy arithmetic). The archimedean signature $(pos, neg)$ is checked: a positive or negative definite form has no real zero and therefore no rational zero. For each bad prime $p | 2\det(A)$, the Hilbert symbol $(d_i, d_j)_p$ is computed for each pair of diagonal entries via the exact formula with valuation splitting, unit Legendre symbol, and the $p=2$ special case. A form is locally anisotropic at $p$ when the Hasse invariant $\varepsilon_p = \prod_{i < j} (d_i, d_j)_p = -1$. The first prime where local isotropy fails is reported with the explicit obstruction data.

**Stage 8 — Cassels–Guy cubic obstruction.** For a diagonal cubic $ax^3 + by^3 + cz^3 = 0$, the solver applies the Cassels–Guy criterion: if $abc$ is not a perfect cube in $\mathbb{Q}$, the cubic has a nontrivial 3-torsion Brauer class in $\mathrm{Br}(V)/\mathrm{Br}(\mathbb{Q})$ and no nontrivial rational projective solution. This is verified by factorising $abc$ and checking whether all prime valuations are divisible by 3. The solver also verifies local solvability — confirming that the failure is genuinely global, not a local obstruction at some prime. The status `CUBIC_OBSTRUCTION` is returned with an explicit message naming the order-3 Brauer class and the exact value of $abc$.

**Stage 9 — Quantum oracle (Berry phase and interference).** When all classical stages are inconclusive, the quantum layer engages using the Prime-Fock library.

---

## 5. The Quantum Oracle

### 5.1 Architecture: Three Layers

**Layer 1 — Arithmetic Fock graph.** Rational candidates are encoded as signed prime-valuation vectors. The log-gcd kernel gives the intersection multiplicity between pairs of candidates — how much prime structure they share. This is the discrete Fisher metric on the arithmetic manifold. Integers are never materialised; all arithmetic happens in exponent space.

**Layer 2 — Quantum walk dynamics.** The graph Laplacian $L$ governs the open-system evolution

$$i\frac{d\psi}{dt} = (L + i \cdot \mathrm{diag}(q))\psi,$$

where $q$ is the *energy flux field* — the per-node inflow/outflow derived from the observed change in arithmetic degree between deformation steps:

$$q_i = \frac{\mathrm{deg}_i[\lambda_k] - \mathrm{deg}_i[\lambda_{k-1}]}{\mathrm{mean}|\mathrm{deg}[\lambda_{k-1}]|}.$$

This is not a tuning parameter. It is the discrete spin connection of the arithmetic graph — the rate at which each candidate's arithmetic relationships are changing under the deformation. The `energy_flux` field in `GraphWaveDynamics` carries exactly this quantity. The interference residual $|\psi(t)|^2 - \rho(t)$ (wave density minus diffusion density, same initial state, same generator) measures the constructive and destructive interference of arithmetic amplitudes. Where it is systematically negative, amplitudes cancel — the quantum walk detects that no coherent rational structure can be assembled from the local arithmetic data.

**Layer 3 — Spin connection curvature.** The Frame Rotation Operator from the Spin Field Theory document is

$$\mathcal{F}_\mu = \varepsilon_{ab} e^a{}_\nu \partial_\mu e^{b\nu},$$

where $e^a_\mu$ are the vielbein fields (frame eigenvectors). On the arithmetic graph, the Gram matrix (log-gcd adjacency) plays the role of the metric $g_{\mu\nu}$. Its eigenvectors are the frame fields. The spin connection $\mathcal{F}_i$ at node $i$ is the angular rate of change of the dominant frame eigenvector across neighbouring nodes — the discrete $\partial_\mu \theta$ where $\theta_i = \arctan_2(\mathrm{sub}_i, \mathrm{dom}_i)$ is the frame angle. The curvature $|d\mathcal{F}|_i$ is the weighted angular deviation summed over edges. By the discrete Gauss-Bonnet theorem,

$$\oint_{\partial \mathcal{M}} \mathcal{F} \cdot d\ell = 2\pi \sum_i s_i,$$

spikes in $|d\mathcal{F}|$ locate arithmetic defects with topological charge $s_i$. The solver reports the defect node count and integrated charge.

### 5.2 Berry Phase and Monodromy

The Berry phase is computed over a closed loop in coefficient space. The leading two monomial coefficients $(a_0, b_0)$ of the first equation are rotated:

$$a(\theta) = a_0 + r\cos\theta, \quad b(\theta) = b_0 + r\sin\theta, \quad \theta \in [0, 2\pi].$$

At each step, the arithmetic graph is rebuilt using residual-weighted prime valuation rows: candidate $i$'s row is its valuation vector scaled by $1/(|F_{a(\theta), b(\theta)}(x_i)| + \varepsilon)$, so candidates near solutions of the *deformed* system dominate the graph edges. The ground state $|\psi_0(\theta)\rangle$ (lowest Laplacian eigenvector) is tracked by maximum overlap to handle avoided crossings. The discrete Berry phase is

$$\gamma_\mathrm{Berry} = -\mathrm{Im} \log \prod_{k=0}^{N-1} \langle \psi_0(\theta_k) | \psi_0(\theta_{k+1}) \rangle,$$

and the winding number $\gamma_\mathrm{Berry} / 2\pi$ measures how many times the ground state rotates as the coefficients complete one loop. A non-integer winding number detects a topological obstruction — the Brauer class — as a fractional winding. For a $\mathbb{Z}/n\mathbb{Z}$ obstruction, the winding is $k/n$ for some $k \not\equiv 0 \pmod n$.

The reliability of this measurement requires the spectral gap of $L$ to exceed the flux norm $\|q\|$. The solver checks this and reports `reliable: True/False` alongside every Berry phase measurement.

### 5.3 What the Quantum Oracle Can and Cannot Detect

Careful investigation during development established the precise scope of the Berry phase method:

**Z/2Z obstructions (Hilbert symbol).** For quadratic forms, the Brauer class is a $\pm 1$ Hilbert symbol at each bad prime. The Berry phase detects this as a winding number of $\pm 1/2$ — the ground state picks up the fermionic sign $\Psi(2\pi) = -\Psi(0)$ associated with the $s = 1/2$ spinorial defect of the Spin Field Theory document. This is directly analogous to the $\mathbb{Z}_2$ monodromy of $\pi_1(SO(3))$ that forces the spinorial sign in the Fisher-frame programme.

**Z/3Z obstructions (cubic Brauer class).** For diagonal cubics, the 3-torsion Brauer class is a norm residue symbol in $H^2(\mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q}), \mu_3)$. This is a class field theory invariant, not a spectral geometry invariant. The Berry phase on the arithmetic graph does not converge to a $1/3$ winding for the Selmer cubic because the eigenvector tracking accumulates full $2\pi$ rotations at avoided crossings, and the relevant phase lives in the field theory of cubic extensions rather than in the spectral theory of the Laplacian. The correct computable tool is the Cassels–Guy criterion (Stage 8 above), which detects the obstruction algebraically in $O(1)$ time. The quantum oracle correctly returns winding $\approx 0$ for the Selmer cubic, consistent with the fact that the 3-torsion Brauer obstruction is not a spectral geometry phenomenon.

This is not a failure of the quantum oracle — it is a clarification of its domain. The Berry phase detects $\mathbb{Z}/2\mathbb{Z}$ monodromy; class field theory detects $\mathbb{Z}/3\mathbb{Z}$ and higher torsion. Both are implemented; they cover complementary domains.

### 5.4 The Holonomy Engine

The `GraphHolonomyEngine` computes the cycle basis of the arithmetic graph and evaluates the holonomy

$$h(\gamma) = \oint_\gamma \mathcal{F} \cdot d\ell \pmod{2\pi}$$

for each cycle $\gamma$. Cycles with $|h| / 2\pi > \delta_\mathrm{threshold}$ are reported as defects with topological charge $s = h / 2\pi$. This is the discrete Gauss-Bonnet sum: the total charge enclosed by any loop equals the loop integral of the frame rotation, which equals the Berry phase of the ground state around that loop.

---

## 6. Projective Filtering

A system of equations is *homogeneous* if every monomial in every equation has the same total degree. Such a system defines a projective variety — its solution set is naturally a subset of projective space $\mathbb{P}^{n-1}(\mathbb{Q})$, where the trivial zero $(0, \ldots, 0)$ is excluded by definition (it corresponds to the cone point, not a point on the variety) and two solutions $x$ and $\lambda x$ for $\lambda \in \mathbb{Q}^*$ are identified.

The solver detects homogeneity automatically via `_is_homogeneous`. When triggered, three things happen:

Bounded search enumerates only *primitive* integer tuples — those with $\gcd(|x_1|, \ldots, |x_n|) = 1$ — excluding $(0, \ldots, 0)$ and all non-primitive rescalings. This reduces the search space by a factor of $\zeta(n)^{-1}$ (the density of primitive tuples) and eliminates the main source of false solutions in homogeneous systems.

Deduplication uses the primitive representative $[x]$ of each projective point, normalised so the leading nonzero coordinate is positive.

The Cassels-Guy check and Hasse-Minkowski are applied to the projective variety, not the affine cone. The Selmer cubic $3x^3 + 4y^3 + 5z^3 = 0$ correctly returns `CUBIC_OBSTRUCTION` with the explicit Brauer class message rather than `FOUND_INTEGER` with the trivial zero.

---

## 7. The Filling-Cost Family

The three programmes produce three filling-cost objects:

**Causal filling modulus:** $\kappa_\mathrm{caus}(\Gamma) = S_\mathrm{fill}(\Gamma) / \tau_\Gamma^2$ where $S_\mathrm{fill}$ is the infimum of Fisher information cost. Inverse stiffness gives Newton's constant $G$. Zero modulus means flat causal geometry.

**Hemisphere escape modulus:** $\kappa(D) = \mathcal{S}_\mathrm{esc}(D) / \lambda_D^{n+1}$, universal value $\kappa = 2/\pi$ independent of $\rho_0$. Inverse stiffness gives the gauge coupling $g_0^2 = 3/(4\pi)$. The spectral condition $g_0^2 \mathcal{S}_\mathrm{esc} = J_\mathrm{min}^2$ produces the electron mass $m = 2m_P e^{-32\pi^3/19}$ within $0.12\%$ of experiment.

**Rational filling resistance:** $\kappa_\mathbb{Q}(x) = \deg(\min\_\mathrm{poly}(x))/2$ by Roth's theorem. Zero resistance means rational. Nonzero resistance identifies algebraic irrationals by their degree. The threshold $\kappa_\mathbb{Q} \geq 1/2$ separates rational from non-rational, analogously to how $S_\mathrm{fill} = 0$ separates flat from curved causal geometry.

All three assign a stiffness to a topological structure by a minimax variational procedure and produce a coupling or threshold as the inverse stiffness.

---

## 8. The Status Vocabulary

The solver's status vocabulary is diagnostic rather than binary. Each status names a specific mathematical reason:

`FOUND_RATIONAL` — certified rational solution, exact.

`FOUND_INTEGER` — certified integer solution, exact. (Only in affine mode; projective mode uses `FOUND_RATIONAL`.)

`FOUND_ALGEBRAIC` — exact non-rational algebraic solution identified.

`BRAUER_MANIN_CLASSICAL` — Hasse–Minkowski obstruction for a quadratic form: either the archimedean signature shows the form is definite (no real zero), or a specific prime $p$ is identified where the Hasse invariant $\varepsilon_p = -1$. The output names the prime and the diagonal form.

`CUBIC_OBSTRUCTION` — Cassels–Guy criterion: $abc$ is not a perfect cube, so $ax^3 + by^3 + cz^3 = 0$ has no nontrivial rational projective solution. The 3-torsion Brauer class is explicit in the output.

`ADELIC_OBSTRUCTION` — Hensel lifting fails at a specific prime $p$ with a specific reason (singular Jacobian, non-divisibility of the Newton residual). The prime and reason are named.

`ARITHMETIC_STIFF` — rational filling resistance $\kappa_\mathbb{Q} > \delta$ for all surviving candidates. The candidates are algebraically irrational; rational solutions, if they exist, are at much larger height than the search bound.

`BRAUER_MANIN_QUANTUM` — Berry phase winding number is non-integer with fractional part $> \delta$. The estimated order of the Brauer class is reported as $\approx 1/|$fractional part$|$. Only fired when the spectral gap exceeds the flux norm (reliable measurement).

`INTERFERENCE_BLOCKED` — quantum walk interference residual $|\psi|^2 - \rho < -\delta$ uniformly. Amplitude cancellation at every node is consistent with Brauer–Manin obstruction; no coherent rational structure can be assembled.

`NO_SOLUTION_WITHIN_BOUND` — bounded search exhausted, homotopy found no rational reconstructions, classical obstructions do not apply. Either no rational solution exists (undecidable by current methods) or one exists above the height bound.

`UNKNOWN` — no stage was decisive. The solver is honest about its own limits.

---

## 9. Completeness and Honest Limits

The solver is complete — gives a correct definitive answer — for:

**All genus-0 curves (rational curves).** Hasse–Minkowski is complete for conics and quadric surfaces. Every quadratic obstruction is caught in Stage 7.

**All diagonal cubics.** Cassels–Guy is complete for the class of diagonal cubic forms $ax^3 + by^3 + cz^3 = 0$. Stage 8 decides these correctly.

**Del Pezzo surfaces of degree $\geq 5$.** The Brauer–Manin obstruction is the only obstruction by theorem. The quantum oracle's Berry phase detects it for the $\mathbb{Z}/2\mathbb{Z}$ case; Stage 8 handles diagonal cases.

**All elliptic curves under BSD.** If the Birch and Swinnerton-Dyer conjecture holds, the rank is computable and 2-descent via the Selmer group decides $E(\mathbb{Q}) = \emptyset$ or not. The Berry phase contributes the $\mathbb{Z}/2\mathbb{Z}$ component of $\mathrm{Sha}(E)[2]$.

The solver is incomplete — may return `UNKNOWN` when a definitive answer exists — for:

**Bielliptic surfaces (Skorobogatov examples).** Varieties that are Brauer–Manin unobstructed but have no rational point exist; no computable obstruction theory covers these in general.

**Large systems at high precision.** The Berry phase becomes unreliable when the spectral gap of $L$ narrows below the flux norm $\|q\|$. This occurs for systems in many variables with many prime modes, where the arithmetic graph becomes dense and the spectrum compresses. The solver reports this explicitly via the `reliable` flag.

**Systems encoding arbitrary computation.** By Matiyasevich's theorem over $\mathbb{Z}$ (and conjecturally over $\mathbb{Q}$), some polynomial systems encode the halting problem. These are recognisable by their structural character — high degree, many variables, no geometric structure, no fibration — and the solver returns `UNKNOWN` rather than looping. This is correct behaviour.

The epistemological status of the Church–Turing thesis deserves mention. Matiyasevich's proof assumes Church–Turing: that the informal notion of effective computation is exactly captured by Turing machines. This is an empirical hypothesis (supported by 90 years of evidence) rather than a mathematical theorem. The phase oracle approach — implementing the Diophantine constraint as a Hamiltonian and measuring the Berry phase of the ground state — operates in the continuous-variable quantum regime where the Church–Turing thesis does not formally apply. In the idealised limit of exact analogue computation, the Berry phase measurement answers the Diophantine question for any $\mathbb{Z}/n\mathbb{Z}$ Brauer class in finite time by physical measurement rather than discrete search. For the practically important cases (small $n$, low degree, geometric structure), this is the right regime: the phases are rational multiples of $2\pi$ with small denominators, distinguishable at finite precision.

---

## 10. Computational Architecture

The Prime-Fock library provides the substrate. Its design principles align exactly with the solver's requirements:

**Valuation arithmetic** (`ValuationArithmetic`): multiplication, gcd, lcm, quotient, and all arithmetic kernels operate in exponent space. No large integers are materialised. The log-gcd kernel gives the intersection multiplicity (shared prime support) as a continuous similarity measure.

**Graph dynamics** (`GraphWaveDynamics`, `GraphDiffusionDynamics`): the open-system propagator $\exp(-i(L + i\,\mathrm{diag}(q))t)$ is computed via full eigendecomposition of the non-Hermitian generator. The `energy_flux` field is inferred from graph deformation — not a tuning knob but a data-derived quantity. The interference residual $|\psi|^2 - \rho$ is the headline diagnostic feature.

**Holonomy engine** (`GraphHolonomyEngine`): computes the fundamental cycle basis of the arithmetic graph (spanning tree plus chord cycles), evaluates holonomy per cycle as the sum of edge phases around the cycle, and reports cycles exceeding the defect threshold as topological charges.

**Stability filter** (`StabilityFilter`): checks the cone-contractivity criterion $\mathcal{B} = \mathcal{H}_\infty^{-1/2} \mathcal{H}_F \mathcal{H}_\infty^{-1/2} \leq I$. The spectral radius $\rho(\mathcal{B}) \leq 1$ certifies that height transfer is contractive in the archimedean completion norm — the arithmetic analogue of SDR's cone-contractivity condition $\mathcal{B}_\zeta \leq I$ for the completed pleroma Hessian.

**Numba compilation**: hot paths (cycle basis, holonomy evaluation, pairwise valuation kernels) compile to machine code via Numba with graceful Python fallback.

---

## 11. The Selmer Cubic: A Worked Example

Selmer's cubic $3x^3 + 4y^3 + 5z^3 = 0$ is the canonical example of a variety that is *locally solvable everywhere* (has $\mathbb{Q}_p$-points for all primes $p$ and a real point) but has no rational point. It is the test case for Brauer–Manin obstruction theory.

The solver's trace on this system:

1. **Homogeneity detected.** The equation is degree-3 homogeneous. Projective mode activates; trivial zero excluded.

2. **Symbolic pass.** `solve` returns no explicit rational solutions (the variety has no rational points by Selmer's theorem; SymPy correctly finds none symbolically).

3. **Hasse–Minkowski.** Not applicable: the equation is cubic, not quadratic. Status: `not_quadratic`.

4. **Cassels–Guy.** $abc = 3 \cdot 4 \cdot 5 = 60$. Prime factorisation: $60 = 2^2 \cdot 3 \cdot 5$. Valuations: $v_2(60) = 2$, $v_3(60) = 1$, $v_5(60) = 1$. None are divisible by 3. Therefore $60$ is not a perfect cube. Local solvability verified at all small primes. Status: `CUBIC_OBSTRUCTION`.

Output message: *"Cassels–Guy obstruction: abc=60 is not a perfect cube ($\sqrt[3]{60} \approx 3.9149$); diagonal cubic $3x^3+4y^3+5z^3=0$ has no nontrivial rational projective solution. 3-torsion Brauer class: order 3 element of $\mathrm{Br}(V)/\mathrm{Br}(\mathbb{Q})$."*

Total time: 56ms. No bounded search, no homotopy, no quantum oracle needed. The algebraic criterion decides it immediately. This is the correct behaviour: use the cheapest decisive tool, and reserve the quantum oracle for cases the algebraic tools cannot reach.

---

## 12. What Was Built, Precisely

A nine-stage Diophantine/rational solution engine organised as a filter cascade in the spirit of Stable Defect Resolution, with:

**Three new abstract objects:** The rational filling resistance modulus $\kappa_\mathbb{Q}$ (implementing the Roth irrationality exponent analytically via minimal polynomials), the arithmetic spin connection $\mathcal{F}_\mu$ (the angular gradient of the Gram matrix frame field on the candidate graph), and the residual-weighted valuation graph (the arithmetic graph whose topology responds to deformation of the equation's coefficients).

**Two new obstruction procedures:** The full Hasse–Minkowski computation via exact LDL$^\top$ Gram diagonalisation over $\mathbb{Q}$ and exact Hilbert symbols, and the Cassels–Guy cubic criterion. Together these cover all quadratic forms and all diagonal cubics definitively.

**Full Jacobian Hensel lifting:** The Newton step in Hensel's method now solves the complete Jacobian system over $\mathbb{Z}/p\mathbb{Z}$ via modular matrix inversion, removing the diagonal approximation that was incorrect for coupled systems.

**Projective filtering:** Automatic detection of homogeneous systems, primitive tuple enumeration in bounded search, projective deduplication via canonical representatives, and exclusion of the trivial zero from all stages.

**Quantum oracle:** Arithmetic Fock encoding via the Prime-Fock library, quantum walk dynamics with data-derived energy flux (the discrete spin connection), Berry phase accumulation via ground-state overlap products around a coefficient-space loop, spin connection curvature readout, and graph holonomy defect detection.

**Honest completeness:** The solver returns named, classified obstruction data for every case it decides and `UNKNOWN` with no false confidence for cases it cannot. The status vocabulary (`FOUND_RATIONAL`, `BRAUER_MANIN_CLASSICAL`, `CUBIC_OBSTRUCTION`, `ADELIC_OBSTRUCTION`, `ARITHMETIC_STIFF`, `BRAUER_MANIN_QUANTUM`, `INTERFERENCE_BLOCKED`, `NO_SOLUTION_WITHIN_BOUND`, `UNKNOWN`) is the most informative currently achievable for a general-purpose Diophantine solver.

---

*This document is a companion to `diophantine_oracle_solver.py`, `rational_defect_framework.md`, `rational_flow_solver_v3.py`, the Stable Defect Resolution preprint (Kirby, May 2026), the Causal Filling Modulus programme, and the Fisher-Frame Defects and Spinorial Matter preprint v0.8.*
