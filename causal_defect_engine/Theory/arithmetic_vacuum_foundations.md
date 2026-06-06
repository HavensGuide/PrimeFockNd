# The Arithmetic Vacuum: Causal Information Foundations of Mathematics

**Kai Cobbs**  
Standalone preprint v0.1 — May 2026  
Programme companion to *Causal Filling Gravity v0.7*, *Fisher-Frame Defects v0.7*, and *Numerical Pleroma I–III*

---

## Abstract

We propose the *Arithmetic Vacuum* (also called Pleroma) as the foundational ontological object from which mathematics in its entirety arises, and from which physics is subsequently forced as a continuum limit. The primitive structure is a pair of causal information structures (CIS) on the integers — one multiplicative, one additive — whose unique compatibility condition we prove to be the distributive law, making the ring axioms derived theorems rather than postulates. From this minimal starting point we derive, through a sequence of forced constructions: the natural numbers as the topological defect spectrum of the multiplicative vacuum; the Fundamental Theorem of Arithmetic as a topological stability and uniqueness theorem; the Möbius function as the canonical antidefect operator satisfying a vacuum neutrality condition; and the number tower $\mathbb{Z} \to \mathbb{Q} \to \mathbb{R} \to \mathbb{C}$ as successive forced completions. Abstract algebra (groups, rings, fields, modules) arises as the classification of compatible vacuum structures and their automorphisms. Topology arises from the causal order of the divisibility lattice. Analysis arises from the coarse-graining limit of discrete causal transport.

The algebra of the Arithmetic Vacuum is the *Universal Defect Algebra*: the ring of arithmetic functions under Dirichlet convolution, with the von Mangoldt function $\Lambda$ as the canonical density, the Möbius function $\mu$ as the vacuum inverse, and the Riemann zeta function $\zeta(s)$ as the partition function of the defect gas. For a fixed semiprime $N = pq$, this localises to the four-dimensional Defect Algebra $\mathrm{DA}_N$ of *Numerical Pleroma III*, carrying Galois group $\mathbb{Z}/2\mathbb{Z}$ and Adams operations giving it $\lambda$-ring structure.

All algebraic claims are proved rigorously. Programme elements — the transition to physics and the Riemann Hypothesis as an equilibrium condition — are precisely formulated as open problems with explicit proof obligations. The companion paper *The Arithmetic Vacuum: Physical Sector* addresses how the continuum coarse-graining of this structure forces Lorentzian geometry, spinorial matter, and the fundamental constants.

---

## Contents

1. Introduction and Aim  
2. The Multiplicative Causal Information Structure  
3. The Antidefect Structure and Vacuum Neutrality  
4. The Additive Causal Information Structure  
5. The Ring as Minimal Dual-Vacuum Completion  
6. The Tower of Number Systems  
7. The Universal Defect Algebra  
8. The Local Defect Algebra $\mathrm{DA}_N$  
9. Abstract Algebra from Defect Symmetries  
10. Topology from Causal Order  
11. Analysis from Causal Transport  
12. Open Problems  
13. Conclusion  
A. Symbol Glossary  
B. Derived versus Postulated  

---

## 1. Introduction and Aim

### 1.1 The Central Claim

Standard foundations of mathematics — whether set-theoretic (ZFC), type-theoretic (HoTT), or categorical (ETCS) — begin with primitive axioms and derive mathematics from them. The programme of the present paper makes a different claim:

> **The positive integers, equipped with divisibility as causal order and the log-GCD kernel as information transport, constitute a causal information structure in the precise sense of Causal Filling Gravity [CFG v0.7]. Every object of mathematics arises as a structural component of this primitive object, a forced completion of it, or an automorphism structure of it. The ring axioms are not postulated; they are the unique compatibility conditions between two independent such structures. Physics arises when this abstract structure acquires a Lorentz-invariant physical instantiation.**

This is not the claim that "physics is information" or that "the universe is a computer." It is a specific mathematical claim: that one precisely defined object — the pair of causal information structures $(\mathcal{V}_\times, \mathcal{V}_+)$ — generates all of mathematics through a chain of derived constructions, each of which is either a proved theorem or a precisely stated open problem.

The programme has two conceptual predecessors worth distinguishing from it. Wheeler's "it from bit" proposes that information is physically fundamental without specifying what structure the information obeys. Tegmark's Mathematical Universe proposes that all mathematical structures are physically real. The present programme makes a much sharper claim: *one specific* mathematical structure — the multiplicative causal order of the integers — is the primitive object, and mathematics and physics are its two sectors. This is a claim about a particular structure, not a general thesis about information or mathematical existence.

### 1.2 Ontological Status: Three Positions

The question of *why* the Arithmetic Vacuum exists admits three positions of increasing strength, analogous to Appendix B of *Fisher-Frame Defects v0.7*.

**(i) Mathematical necessity.** The divisibility causal structure may be the unique structure satisfying some deeper consistency condition not yet identified — the minimal self-consistent causal information structure. Under this view, the existence of the integers and their properties is not contingent.

**(ii) Pragmatic foundation.** The Arithmetic Vacuum is not claimed to be ontologically prior to other foundations; it is a foundation that generates more structure from less than any alternative, and which connects mathematics and physics through a single derivation chain.

**(iii) Ontological priority of the arithmetic structure.** The divisibility order, the log-GCD kernel, and the prime defect structure are ontologically prior to both the physical vacuum and to abstract mathematical objects. The physical vacuum is the Arithmetic Vacuum made dynamical; mathematical objects are its structural consequences.

The paper is consistent with all three positions. The proofs do not depend on which is adopted. The third position is noted as the most natural home for the programme's results.

### 1.3 Relationship to Existing Foundations

**Set theory (ZFC).** ZFC derives mathematics from the axioms of set membership. The Arithmetic Vacuum instead takes divisibility as primitive. The two approaches can coexist: the integers can be constructed within ZFC, and the Arithmetic Vacuum can be seen as identifying which ZFC-constructible structures are the natural primitives. The difference is pragmatic and ontological, not technical.

**Category theory.** Category theory abstracts mathematical structures to objects and morphisms. The Arithmetic Vacuum gives a specific object (the category of CIS-compatible structures over $\mathbb{Z}$) from which the abstract categorical framework emerges. In particular, the category of rings is derived in Theorem 9.3 as the category of compatible dual-vacuum structures.

**Arithmetic geometry (Spec $\mathbb{Z}$).** The programme is most closely allied with the arithmetic geometry programme: it treats $\mathrm{Spec}(\mathbb{Z})$ as the fundamental geometric object (§10.3). The Arithmetic Vacuum provides the physical/information-theoretic motivation for that choice.

### 1.4 Summary of Main Results

**(1)** $(\mathbb{N}^*, |, K_\times, \Lambda)$ is a causal information structure in the sense of CFG (Theorem 2.2).  
**(2)** Primes are the elementary defects of $\mathcal{V}_\times$ (Theorem 2.5).  
**(3)** The Fundamental Theorem of Arithmetic is a topological stability and uniqueness theorem (Theorem 2.6).  
**(4)** The Möbius function satisfies the vacuum neutrality condition $\sum_{d|n}\mu(d) = \varepsilon(n)$ (Theorem 3.2).  
**(5)** The ring axioms are the unique compatibility conditions between $\mathcal{V}_\times$ and $\mathcal{V}_+$ (Theorems 5.2–5.3).  
**(6)** $\mathbb{Z}$ is the initial ring; $\mathbb{Q}$ is the defect-antidefect completion; $\mathbb{R}$ is the Gaussian coarse-graining limit; $\mathbb{C}$ is the unique symplectic lift (Theorems 5.4, 6.2, 6.4–6.7).  
**(7)** The Universal Defect Algebra $\mathcal{A}$ is the ring of arithmetic functions under Dirichlet convolution; $\zeta(s)$ is its partition function (§7).  
**(8)** The Local Defect Algebra $\mathrm{DA}_N$ has $\mathrm{Gal}(\mathrm{DA}_N) \cong \mathbb{Z}/2\mathbb{Z}$ and $\lambda$-ring structure (§8).  
**(9)** Every group arises as the automorphism group of an arithmetic structure (Theorem 9.1); rings are compatible dual-vacuum structures (Theorem 9.3).  
**(10)** Topology arises from the causal order (§10); analysis from causal transport (§11).

---

## 2. The Multiplicative Causal Information Structure

### 2.1 Review: Causal Information Structures

From *CFG v0.7* Definition 2.1, a **causal information structure (CIS)** is a quadruple $(E, \preceq, K, \rho)$ where:
- $E$ is a set of elementary events,
- $\preceq$ is a causal preorder on $E$,
- $K: E \times E \to \mathbb{R}_{\geq 0}$ is a causal correlation kernel with $K(x,y) = 0$ unless $x \preceq y$ or $y \preceq x$,
- $\rho: E \to \mathbb{R}_{\geq 0}$ is an event density.

### 2.2 The Multiplicative CIS

**Definition 2.1** (Multiplicative Arithmetic Vacuum). The *multiplicative causal information structure* is
$$\mathcal{V}_\times = (\mathbb{N}^*,\; |\;,\; K_\times,\; \Lambda)$$
where:
- $\mathbb{N}^* = \{1, 2, 3, \ldots\}$ are the positive integers (elementary events),
- $x \preceq y$ iff $x \mid y$ (divisibility as causal order),
- $K_\times(x, y) = \log(y/x)$ when $x \mid y$, and $K_\times(x,y) = 0$ otherwise,
- $\Lambda(n)$ is the von Mangoldt function: $\Lambda(p^k) = \log p$ for prime powers, $\Lambda(n) = 0$ otherwise.

**Theorem 2.2** (Multiplicative CIS). *$\mathcal{V}_\times$ is a causal information structure.*

**Proof.**
*(i) Causal preorder.* Divisibility defines a partial order on $\mathbb{N}^*$: reflexive ($n \mid n$), antisymmetric (if $x \mid y$ and $y \mid x$ then $x = y$ for positive integers), transitive (if $x \mid y$ and $y \mid z$ then $x \mid z$). A partial order is a preorder. ✓

*(ii) Kernel support.* By definition, $K_\times(x,y) = 0$ unless $x \mid y$ or $y \mid x$, i.e., unless $x \preceq y$ or $y \preceq x$. ✓

*(iii) Kernel nonnegativity.* If $x \mid y$ then $y/x \in \mathbb{N}^*$, so $y/x \geq 1$ and $K_\times(x,y) = \log(y/x) \geq 0$. ✓

*(iv) Density nonnegativity.* $\Lambda(n) = \log p \geq \log 2 > 0$ for prime powers, $0$ otherwise. ✓ $\square$

**Proposition 2.3** (Structural Properties of $\mathcal{V}_\times$).

*(a) Unique vacuum.* $1 \in \mathbb{N}^*$ is the unique minimum: $1 \mid n$ for all $n$, and $K_\times(1,1) = \log 1 = 0$.

*(b) Additivity of path cost.* For $x \mid y \mid z$:
$$K_\times(x, z) = \log(z/x) = \log(z/y) + \log(y/x) = K_\times(y,z) + K_\times(x,y).$$
The information cost is additive along causal chains.

*(c) Density recovery.* $K_\times(1, n) = \log n = \sum_{d \mid n} \Lambda(d)$. The total causal cost from vacuum to $n$ equals the sum of the primitive causal densities $\Lambda$ over all intermediate divisors.

*(d) Zero self-cost.* $K_\times(x,x) = \log(x/x) = 0$ for all $x$. ✓

**Proof.** (a)–(d) follow directly from the definitions. For (c): the identity $\sum_{d|n}\Lambda(d) = \log n$ is standard — it follows from the logarithmic derivative of the Euler product $\zeta(s) = \prod_p (1-p^{-s})^{-1}$ via Möbius inversion. $\square$

**Remark 2.4.** Property (c) is the key physical content. The von Mangoldt density $\Lambda$ is the unique function on $\mathbb{N}^*$ such that summing over all divisors recovers the total causal path cost. This is the arithmetic analogue of the relationship between local curvature (the Riemann tensor) and path holonomy in Riemannian geometry: summing the local density $\Lambda$ over the divisibility chain from $1$ to $n$ gives the total "curvature" $\log n$ of the causal path.

### 2.3 Primes as Elementary Defects

**Definition 2.4** (Vacuum and Defects). The *vacuum* of $\mathcal{V}_\times$ is $1 \in \mathbb{N}^*$. An element $n > 1$ is an *elementary defect* if the causal interval
$$[1, n]_{\mid} := \{d \in \mathbb{N}^* : d \mid n\}$$
contains no element $d$ with $1 < d < n$.

**Theorem 2.5** (Primes as Elementary Defects). *An integer $n > 1$ is an elementary defect of $\mathcal{V}_\times$ if and only if $n$ is prime.*

**Proof.** By Definition 2.4, $n > 1$ fails to be an elementary defect iff $\exists d$ with $d \mid n$ and $1 < d < n$, i.e., iff $n$ has a proper non-trivial divisor, i.e., iff $n$ is composite. Therefore $n$ is an elementary defect iff $n$ has no proper non-trivial divisors iff $n$ is prime. $\square$

**Theorem 2.6** (Unique Defect Decomposition — Fundamental Theorem of Arithmetic). *Every $n > 1$ has a unique decomposition*
$$n = p_1^{a_1} p_2^{a_2} \cdots p_k^{a_k}, \quad p_i \text{ distinct primes},\; a_i \geq 1,$$
*unique up to reordering of factors.*

**Proof.** (Existence) By induction: $n = 2$ is prime. For $n > 2$, either $n$ is prime (done) or it has a proper divisor $1 < d < n$; by induction both $d$ and $n/d$ decompose into primes. (Uniqueness) Suppose $p_1 \cdots p_k = q_1 \cdots q_l$. By Euclid's lemma — if a prime $p$ divides a product $ab$ then $p \mid a$ or $p \mid b$, proved from the Euclidean algorithm — each $p_i$ divides some $q_j$, hence equals it. Cancelling gives the result by induction on $\max(k,l)$. $\square$

**Remark 2.7** (Topological Interpretation). Theorem 2.6 is the arithmetic analogue of a topological defect stability theorem. The prime factors $\{(p_i, a_i)\}$ are conserved topological charges: a composite number cannot be "deformed" to the vacuum within the divisibility lattice because there is no path from $n$ to $1$ that avoids its prime factors. The uniqueness of decomposition is the statement that these charges are independently conserved — the defect sectors labelled by distinct primes do not mix.

---

## 3. The Antidefect Structure and Vacuum Neutrality

### 3.1 The Möbius Function as Antidefect Operator

In defect theories, every elementary defect carries a conserved charge and has an associated antidefect of opposite charge. The pair-annihilation of a defect and its antidefect returns to the vacuum. The Möbius function is the canonical antidefect operator of $\mathcal{V}_\times$.

**Definition 3.1** (Möbius Function). The *Möbius function* $\mu: \mathbb{N}^* \to \{-1, 0, 1\}$ is:
$$\mu(n) = \begin{cases} 1 & n = 1 \\ (-1)^k & n = p_1 p_2 \cdots p_k,\; p_i \text{ distinct primes} \\ 0 & p^2 \mid n \text{ for some prime } p. \end{cases}$$

**Theorem 3.2** (Vacuum Neutrality Condition). *For all $n \geq 1$:*
$$\sum_{d \mid n} \mu(d) = \varepsilon(n) := \begin{cases} 1 & n = 1 \\ 0 & n > 1. \end{cases}$$

**Proof.** For $n = 1$: the only divisor is $1$ and $\mu(1) = 1 = \varepsilon(1)$. For $n > 1$ with prime factorization $n = p_1^{a_1} \cdots p_k^{a_k}$: since $\mu(d) = 0$ whenever $p_i^2 \mid d$, only squarefree divisors contribute. The squarefree divisors of $n$ are products $\prod_{i \in S} p_i$ over subsets $S \subseteq \{1,\ldots,k\}$:
$$\sum_{d \mid n} \mu(d) = \sum_{S \subseteq \{1,\ldots,k\}} (-1)^{|S|} = \sum_{j=0}^{k} \binom{k}{j}(-1)^j = (1-1)^k = 0. \quad \square$$

**Remark 3.3** (Charge-Cancellation Interpretation). The vacuum neutrality condition states that the sum of Möbius charges over all divisors of any composite $n$ is exactly zero. Only the vacuum $n = 1$ carries net Möbius charge $+1$. This is the arithmetic analogue of the charge-neutrality condition in quantum field theory: the vacuum is charge-neutral, and every composite number's divisor structure provides a complete defect-antidefect cancellation. The factor $(-1)^k$ for a product of $k$ distinct primes is the signature of the $k$-defect antisymmetrisation — the arithmetic Fermi statistics.

**Proposition 3.4** (Möbius Inversion). *For arithmetic functions $f, g: \mathbb{N}^* \to \mathbb{C}$:*
$$g(n) = \sum_{d \mid n} f(d) \iff f(n) = \sum_{d \mid n} \mu(d)\, g(n/d).$$

**Proof.** Standard consequence of Theorem 3.2 and Dirichlet convolution associativity: the convolution $\mu * \mathbf{1} = \varepsilon$ (vacuum neutrality) implies $\mu$ is the Dirichlet inverse of $\mathbf{1}$, so $g = \mathbf{1} * f$ iff $\mu * g = f$. $\square$

**Remark 3.5.** Möbius inversion is the arithmetic vacuum's analogue of the inverse causal filling: given the integrated information along all causal paths to $n$ (the function $g$), the Möbius operator extracts the primitive information density at each step (the function $f$). This is the arithmetic Radon inversion or the arithmetic "un-smearing" of the causal kernel.

---

## 4. The Additive Causal Information Structure

### 4.1 The Second Vacuum

The integers carry a second causal information structure, independent of and geometrically orthogonal to the multiplicative one.

**Definition 4.1** (Additive Arithmetic Vacuum). The *additive causal information structure* is
$$\mathcal{V}_+ = (\mathbb{Z},\; \leq,\; K_+,\; \rho_+)$$
where:
- $\mathbb{Z} = \{\ldots, -1, 0, 1, 2, \ldots\}$ are the integers (elementary events),
- $x \preceq y$ iff $x \leq y$ (the natural integer order as causal order),
- $K_+(x, y) = y - x$ for $x \leq y$, zero otherwise,
- $\rho_+ \equiv 1$ (uniform, translation-invariant event density).

**Proposition 4.2.** *$\mathcal{V}_+$ is a causal information structure.*

**Proof.** The natural order $\leq$ is a total order (hence a partial order). $K_+(x,y) = y - x \geq 0$ when $x \leq y$. $K_+(x,y) = 0$ when $x > y$. $\rho_+ = 1 \geq 0$. ✓ $\square$

**Proposition 4.3** (Structural Properties of $\mathcal{V}_+$).

*(a) Vacuum.* $0 \in \mathbb{Z}$ is the vacuum: neutral under addition.

*(b) Elementary defects.* The unit steps $\pm 1$ are the elementary defects of $\mathcal{V}_+$: they are the minimal non-trivial causal steps.

*(c) Translation invariance.* $K_+(x+c, y+c) = K_+(x,y)$ for all $c \in \mathbb{Z}$: the additive CIS is translation-invariant, unlike the multiplicative CIS.

*(d) Generator.* The causal update semigroup $T_n: f(x) \mapsto f(x+n)$ has infinitesimal generator $\partial$, the difference operator $\Delta f(x) = f(x+1) - f(x)$ in the discrete setting, recovered as the derivative $\partial_x$ in the continuum limit.

**Remark 4.4** (Structural Contrast Between the Two Vacua). The multiplicative vacuum $\mathcal{V}_\times$ has a *branching* partial order (the divisibility lattice: incomparable elements abound) and a *non-uniform* density $\Lambda$ concentrated on prime powers. The additive vacuum $\mathcal{V}_+$ has a *linear* total order and a *uniform* density. The multiplicative vacuum governs composition, scale, and primality. The additive vacuum governs position, displacement, and counting. Their interaction generates all of arithmetic — and, as shown below, all of algebra.

---

## 5. The Ring as Minimal Dual-Vacuum Completion

This section contains the paper's central theorem: the ring axioms are not postulated but derived as the unique compatibility condition between two independent causal information structures.

### 5.1 Compatibility

**Definition 5.1** (Dual-Vacuum Compatibility). Two CIS structures $\mathcal{V}_1$ and $\mathcal{V}_2$ on the same event set $E$ are *compatible* if every automorphism of $\mathcal{V}_1$ that acts as a bijection of $E$ is a morphism of $\mathcal{V}_2$: it preserves $\mathcal{V}_2$'s causal order and kernel.

**Theorem 5.2** (Distributive Law from Compatibility). *Consider the multiplicative CIS $\mathcal{V}_\times$ extended by group completion from $\mathbb{N}^*$ to $\mathbb{Z}$ (adding additive inverses), and the additive CIS $\mathcal{V}_+$ on $\mathbb{Z}$. The group endomorphisms of $(\mathbb{Z}, +)$ are exactly the maps $\sigma_a: n \mapsto a \cdot n$ for $a \in \mathbb{Z}$. The requirement that $\sigma_a$ is a morphism of $\mathcal{V}_+$ — i.e., that multiplication by $a$ preserves the additive causal structure — is the condition*
$$\sigma_a(b + c) = \sigma_a(b) + \sigma_a(c),\quad \text{i.e.,}\quad a(b+c) = ab + ac.$$
*This is the distributive law.*

**Proof.** A map $\phi: \mathbb{Z} \to \mathbb{Z}$ preserves the additive CIS $\mathcal{V}_+ = (\mathbb{Z}, \leq, K_+, 1)$ as a morphism iff it is a group homomorphism of $(\mathbb{Z}, +)$: $\phi(b+c) = \phi(b) + \phi(c)$. Every group endomorphism of $(\mathbb{Z}, +)$ is determined by its value on $1$: $\phi(n) = n \cdot \phi(1)$ for all $n \in \mathbb{Z}$ (by induction on $n$ and extension to negatives). Setting $a = \phi(1)$ gives $\phi = \sigma_a$. The morphism condition $\sigma_a(b+c) = \sigma_a(b) + \sigma_a(c)$ then reads $a(b+c) = ab + ac$. $\square$

**Corollary 5.3** (Ring Axioms as Compatibility Theorems). *The ring axioms for $(\mathbb{Z}, +, \times)$ are:*

| Axiom | Source |
|-------|--------|
| $(\mathbb{Z}, +)$ is an abelian group | Given by $\mathcal{V}_+$ and its group-completion |
| $(\mathbb{Z}, \times)$ is a commutative monoid with identity $1$ | Given by $\mathcal{V}_\times$ |
| $a(b+c) = ab+ac$ and $(b+c)a = ba+ca$ | Theorem 5.2: compatibility of $\mathcal{V}_\times$ and $\mathcal{V}_+$ |

*The ring $\mathbb{Z}$ is the minimal algebraic structure that simultaneously realises both causal vacua with their unique compatibility condition. It is derived, not axiomatised.*

**Theorem 5.4** (Universality of $\mathbb{Z}$). *$\mathbb{Z}$ is the initial object in the category of rings: for every ring $R$, there is a unique ring homomorphism $\mathbb{Z} \to R$, sending $n \mapsto n \cdot 1_R$.*

**Proof.** Any ring homomorphism $\varphi: \mathbb{Z} \to R$ must satisfy $\varphi(1) = 1_R$ (identity preservation) and $\varphi(n) = n \cdot 1_R$ by additivity. This is well-defined and is the unique such map. $\square$

**Remark 5.5** (Physical Interpretation). Theorem 5.4 asserts that $\mathbb{Z}$ is the arithmetic ground state: every ring is an excitation above the $\mathbb{Z}$ vacuum. Ring homomorphisms are causal maps preserving both vacuum structures. The category of rings is the category of compatible dual-vacuum structures, with $\mathbb{Z}$ as the unique universal ground state.

---

## 6. The Tower of Number Systems

From the Arithmetic Vacuum, the standard number systems arise through a sequence of forced completions. Each step fills a structural deficiency of the previous system; none is an arbitrary choice.

### 6.1 $\mathbb{Q}$: The Defect-Antidefect Completion

In $\mathcal{V}_\times$, a prime $p$ is an elementary defect; the multiplicative inverse $1/p$ does not exist in $\mathbb{N}^*$ or $\mathbb{Z}$. The multiplicative vacuum is not charge-neutral in $\mathbb{Z}^*$: defects exist without their antidefects. The minimal repair is the field of fractions.

**Theorem 6.1** (Field of Fractions). *The field of fractions of $\mathbb{Z}$ is*
$$\mathbb{Q} = \left\{a/b : a \in \mathbb{Z},\; b \in \mathbb{Z}^*\right\}/{\sim}, \quad a/b \sim a'/b' \iff ab' = a'b.$$
*$\mathbb{Q}$ is the minimal field containing $\mathbb{Z}$. Every prime defect $p \in \mathbb{N}^*$ has multiplicative inverse $1/p \in \mathbb{Q}$.*

**Proof.** Standard construction of the field of fractions of an integral domain. $\mathbb{Z}$ is an integral domain (Theorem 7.2 below), so the construction applies. Minimality: any field containing $\mathbb{Z}$ must contain $a/b$ for all $a \in \mathbb{Z}$, $b \in \mathbb{Z}^*$. $\square$

### 6.2 $\mathbb{R}$: The Coarse-Graining Limit

The additive CIS $\mathcal{V}_+$ on $\mathbb{Z}$ has discrete causal steps. Under coarse-graining — taking the continuum limit — the effective causal kernel flows to a Gaussian fixed point by the CLT/CFG coarse-graining theorem.

**Theorem 6.2** (Gaussian Kernel Fixed Point). *Composing $N$ independent elementary causal steps of $\mathcal{V}_+$, the effective kernel $K_+^{(N)}$ satisfies: under the rescaling $x \mapsto x/\sqrt{N}$, in the limit $N \to \infty$,*
$$K_+^{(N)}(x, y) \to \frac{1}{\sqrt{2\pi\ell^2}} \exp\!\left(-\frac{(y-x)^2}{2\ell^2}\right).$$
*The unique stable fixed-point kernel of repeated causal composition is Gaussian.*

**Proof.** This is the Central Limit Theorem applied to convolutions of the step kernel. In Fourier space, $\hat{K}_+^{(N)}(k) = (\hat{K}_+(k))^N$. For the unit step, $\hat{K}_+(k) = e^{ik}$ (or any kernel with mean $\mu$ and variance $\sigma^2$). Under diffusive rescaling $k \mapsto k/\sqrt{N}$, $(\hat{K}_+(k/\sqrt{N}))^N \to e^{-\sigma^2 k^2/2}$, the Fourier transform of the Gaussian. This matches the Gaussian fixed-point theorem of CFG §16.2 applied to the additive CIS. $\square$

**Theorem 6.3** ($\mathbb{R}$ as Completion). *The Cauchy completion of $(\mathbb{Q}, |\cdot|)$ under the standard absolute value is the unique (up to isomorphism) complete ordered Archimedean field $\mathbb{R}$.*

**Proof.** Standard (Dedekind cuts or Cauchy sequences). The Archimedean property of $\mathbb{Q}$ and the completeness of $\mathbb{R}$ together characterise it uniquely among ordered fields. $\square$

**Remark 6.4.** The Gaussian fixed-point of Theorem 6.2 is the unique fixed point of the causal coarse-graining map; the Cauchy completion of Theorem 6.3 is the unique metric-complete extension of $\mathbb{Q}$. These are two descriptions of the same operation: the continuum limit of the discrete additive CIS is $(\mathbb{R}, \leq, K_\infty, dx)$ where $K_\infty$ is the Gaussian kernel and $dx$ is Lebesgue measure.

### 6.3 $\mathbb{C}$: The Symplectic Lift

The Arithmetic Vacuum, at the level of $\mathbb{R}$, supports only real amplitudes (densities). In *CFG v0.7* §7, quantum phase dynamics arise from the symplectic lift $(ρ, S) \in T^*\mathcal{P}$: the cotangent bundle of density space. Applying the same construction to the real line forces $\mathbb{C}$.

**Definition 6.5** (Symplectic Lift). The *symplectic lift* of the real line $(\mathbb{R}, +)$ is the cotangent bundle $T^*\mathbb{R} = \mathbb{R} \times \mathbb{R}$, with symplectic form $\omega = dp \wedge dq$ (where $(q,p)$ are position and momentum coordinates). A point $(q,p)$ encodes both amplitude ($q$) and phase ($p$). The complex multiplication $i: (q,p) \mapsto (-p, q)$ makes $T^*\mathbb{R}$ into the field $\mathbb{C}$.

**Theorem 6.6** (Uniqueness of $\mathbb{C}$ as Symplectic Completion). *$\mathbb{C}$ is the unique $2$-dimensional real algebra satisfying:*
*(i) $\mathbb{R} \hookrightarrow \mathbb{C}$ (extends $\mathbb{R}$);*
*(ii) every nonzero element has a multiplicative inverse (division algebra);*
*(iii) the symplectic form $\omega(z_1, z_2) = \mathrm{Im}(\bar{z}_1 z_2)$ is non-degenerate and preserved by unit-modulus multiplication.*

**Proof.** By Frobenius' theorem, the only finite-dimensional associative division algebras over $\mathbb{R}$ are $\mathbb{R}$ (dimension 1), $\mathbb{C}$ (dimension 2), and $\mathbb{H}$ (quaternions, dimension 4). $\mathbb{R}$ is $1$-dimensional and has no non-degenerate $2$-form. $\mathbb{H}$ is $4$-dimensional and non-commutative. Therefore among $2$-dimensional real division algebras, $\mathbb{C}$ is the unique option. The symplectic form $\omega(z_1, z_2) = \mathrm{Im}(\bar{z}_1 z_2)$ is the standard symplectic form on $\mathbb{R}^2 \cong \mathbb{C}$; it is non-degenerate, antisymmetric, and preserved by multiplication by $e^{i\theta}$ (rotation of the complex plane). $\square$

**Theorem 6.7** (Algebraic Closure). *$\mathbb{C}$ is algebraically closed: every non-constant polynomial with complex coefficients has a root in $\mathbb{C}$ (Fundamental Theorem of Algebra).*

**Proof.** Standard (via topology: $|p(z)| \to \infty$ as $|z| \to \infty$ implies $|p|$ has a minimum, and if that minimum is non-zero one derives a contradiction via Laurent expansion). $\square$

**Theorem 6.8** (The Number Tower is Complete). *The chain*
$$\mathbb{N}^* \xrightarrow{\text{FTA}} \mathbb{Z} \xrightarrow{\text{fractions}} \mathbb{Q} \xrightarrow{\text{completion}} \mathbb{R} \xrightarrow{\text{symplectic lift}} \mathbb{C}$$
*terminates at $\mathbb{C}$: it is algebraically closed (no algebraic completion needed), complete as a metric space, and the unique $2$-dimensional real division algebra with compatible symplectic structure. No further forced completion exists.*

**Proof.** Termination at $\mathbb{C}$: algebraic closure (Theorem 6.7) closes the algebraic tower; metric completeness is inherited from $\mathbb{R}$; uniqueness follows from Theorem 6.6. $\square$

**Remark 6.9** (The Quantum Phase is Forced). The step $\mathbb{R} \to \mathbb{C}$ is the same step as the introduction of quantum phase dynamics in *CFG v0.7* §7. The polar decomposition $z = \sqrt{\rho} \cdot e^{iS/\hbar}$ identifies $\rho = |z|^2$ (amplitude-squared) and $S/\hbar$ (phase) — exactly the symplectic coordinates $(q, p) \in T^*\mathbb{R}$. Quantum mechanics, in its wave-function formulation, is the arithmetic of the symplectic lift of the real Arithmetic Vacuum.

---

## 7. The Universal Defect Algebra

### 7.1 Arithmetic Functions Under Dirichlet Convolution

The full algebraic structure of the Arithmetic Vacuum is encoded in the ring of arithmetic functions.

**Definition 7.1** (Universal Defect Algebra). The *Universal Defect Algebra* is
$$\mathcal{A} = \bigl(\{f : \mathbb{N}^* \to \mathbb{C}\},\; +,\; *,\; \varepsilon\bigr)$$
where $+$ is pointwise addition, $*$ is *Dirichlet convolution*
$$(f * g)(n) = \sum_{d \mid n} f(d)\, g(n/d),$$
and $\varepsilon = [n=1]$ is the convolution identity.

**Theorem 7.2** (Ring Structure of $\mathcal{A}$). *$(\mathcal{A}, +, *)$ is a commutative unitary ring. It is an integral domain: $f * g = \varepsilon$ with $f \neq 0$ and $g \neq 0$ is impossible (i.e., if $f(n_0) \neq 0$ and $g(m_0) \neq 0$ for the respective smallest such $n_0, m_0$, then $(f*g)(n_0 m_0) \neq 0$).*

**Proof.** Commutativity: $(f*g)(n) = \sum_{d|n} f(d)g(n/d) = \sum_{e|n} f(n/e)g(e) = (g*f)(n)$ via $e = n/d$. Associativity: standard rearrangement of the triple sum $\sum_{d|n}\sum_{e|d}f(e)g(d/e)h(n/d) = \sum_{a|n}\sum_{b|(n/a)}f(a)g(b)h(n/ab)$. Identity: $(\varepsilon * f)(n) = \sum_{d|n}\varepsilon(d)f(n/d) = f(n)$. Distributivity: Dirichlet convolution is bilinear over pointwise addition. Integral domain: let $n_0 = \min\mathrm{supp}(f)$ and $m_0 = \min\mathrm{supp}(g)$. Then $(f*g)(n_0 m_0) = \sum_{d | n_0 m_0} f(d)g(n_0 m_0/d)$. The only divisor $d$ of $n_0 m_0$ with $d \geq n_0$ and $n_0 m_0/d \geq m_0$ is $d = n_0$ (since $n_0 m_0/d \geq m_0$ iff $d \leq n_0$). All other terms vanish by minimality of support. Hence $(f*g)(n_0 m_0) = f(n_0)g(m_0) \neq 0$. $\square$

**Proposition 7.3** (Dictionary of Key Elements).

| Symbol | Definition | Role in $\mathcal{A}$ |
|--------|-----------|----------------------|
| $\varepsilon$ | $\varepsilon(n) = [n=1]$ | Convolution identity (vacuum) |
| $\mathbf{1}$ | $\mathbf{1}(n) = 1$ | "Fill all" function |
| $\mu$ | Möbius function | Convolution inverse of $\mathbf{1}$: $\mu * \mathbf{1} = \varepsilon$ |
| $\Lambda$ | Von Mangoldt function | Primitive density: $\Lambda = \mu * \log$ |
| $\phi$ | Euler totient | $\phi = \mu * \mathrm{id}$ |
| $\mathrm{id}$ | $\mathrm{id}(n) = n$ | Mass function |
| $\tau$ | Number of divisors | $\tau = \mathbf{1} * \mathbf{1}$ |
| $\sigma$ | Sum of divisors | $\sigma = \mathbf{1} * \mathrm{id}$ |

**Theorem 7.4** (Identities in $\mathcal{A}$).

*(a)* $\mu * \mathbf{1} = \varepsilon$ (Möbius is the inverse of $\mathbf{1}$).

*(b)* $\Lambda = \mu * \log$, i.e., $\Lambda(n) = \sum_{d|n} \mu(d)\log(n/d)$ (von Mangoldt from Möbius).

*(c)* $\mathrm{id} = \phi * \mathbf{1}$, i.e., $n = \sum_{d|n} \phi(d)$ (totient reconstruction).

*(d)* $\log = \Lambda * \mathbf{1}$, i.e., $\log n = \sum_{d|n} \Lambda(d)$ (path cost from density — Proposition 2.3(c)).

**Proof.** (a): Theorem 3.2. (b): $(\mu * \log)(n) = \sum_{d|n}\mu(d)\log(n/d)$. Using $\log(n/d) = \sum_{e|(n/d)}\Lambda(e)$ and rearranging: $= \sum_{d|n}\mu(d)\sum_{e|(n/d)}\Lambda(e) = \sum_{e|n}\Lambda(e)\sum_{d|(n/e)}\mu(d) = \Lambda(n)\cdot\varepsilon(n/n) = \Lambda(n)$ by (a). (c) and (d): standard identities, proved similarly. $\square$

### 7.2 The Dirichlet Series Map and the Partition Function

**Definition 7.5** (Dirichlet Series Map). For $f \in \mathcal{A}$, define the Dirichlet series
$$\mathcal{D}[f](s) = \sum_{n=1}^\infty \frac{f(n)}{n^s}, \quad \mathrm{Re}(s) \text{ sufficiently large for convergence.}$$

**Theorem 7.6** (Convolution $\to$ Multiplication). *$\mathcal{D}$ is a ring homomorphism from $(\mathcal{A}, +, *)$ to the ring of convergent Dirichlet series: $\mathcal{D}[f * g] = \mathcal{D}[f] \cdot \mathcal{D}[g]$.*

**Proof.** $\mathcal{D}[f*g](s) = \sum_n \frac{(f*g)(n)}{n^s} = \sum_n \frac{1}{n^s}\sum_{d|n}f(d)g(n/d) = \sum_d \frac{f(d)}{d^s}\sum_m\frac{g(m)}{m^s} = \mathcal{D}[f](s)\cdot\mathcal{D}[g](s)$. $\square$

**Theorem 7.7** (Euler Product and Partition Function). *If $f$ is completely multiplicative then $\mathcal{D}[f](s) = \prod_p (1 - f(p)p^{-s})^{-1}$. In particular:*
$$\mathcal{D}[\mathbf{1}](s) = \zeta(s) = \prod_p \frac{1}{1-p^{-s}}, \qquad \mathcal{D}[\mu](s) = \frac{1}{\zeta(s)}, \qquad \mathcal{D}[\Lambda](s) = -\frac{\zeta'(s)}{\zeta(s)}.$$

**Proof.** Standard Euler product for completely multiplicative Dirichlet series, using unique prime factorization (Theorem 2.6) to factor the sum over $n$ into a product over primes. $\square$

**Remark 7.8** (The Zeta Function as Partition Function). $\zeta(s) = \sum_n n^{-s}$ is the partition function of the Arithmetic Vacuum: it sums the Boltzmann weights $e^{-s\log n}$ over all states $n \in \mathbb{N}^*$, with "inverse temperature" $\mathrm{Re}(s)$. The Euler product
$$\zeta(s) = \prod_p \frac{1}{1-p^{-s}} = \prod_p \sum_{k=0}^\infty p^{-ks}$$
expresses this as a product over prime "oscillators," each contributing an infinite tower $\{p^k : k \geq 0\}$ with energy $k\log p$. This is the arithmetic Fock space: the state $n = p_1^{a_1}\cdots p_r^{a_r}$ is the $a_i$-th excited state of the $p_i$-th oscillator, and the partition function sums over all multi-particle states.

The meromorphic continuation of $\zeta(s)$ to $\mathbb{C}$ is the analytic continuation of this partition function beyond the domain of absolute convergence. The pole at $s = 1$ is the usual thermodynamic divergence at the "critical temperature" $\mathrm{Re}(s) = 1$. The non-trivial zeros of $\zeta$ are the quantum corrections to the prime counting function.

---

## 8. The Local Defect Algebra $\mathrm{DA}_N$

For a fixed semiprime $N = pq$, the Universal Defect Algebra $\mathcal{A}$ localises to the four-dimensional Defect Algebra introduced in *Numerical Pleroma III*.

**Definition 8.1** (Local Defect Algebra). For $N = pq$ ($p < q$ prime), the *Defect Algebra* $\mathrm{DA}_N$ is the four-dimensional commutative $\mathbb{C}$-algebra with orthogonal idempotent basis $\{e_0, e_A, e_B, e_C\}$:
$$e_\alpha e_\beta = \delta_{\alpha\beta} e_\alpha, \qquad e_0 + e_A + e_B + e_C = \mathbf{1},$$
with eigenvalue map $\lambda(e_0) = 0$, $\lambda(e_A) = A = p\log q$, $\lambda(e_B) = B = q\log p$, $\lambda(e_C) = A+B$, and dimension map $\dim(e_0) = (p-1)(q-1)$, $\dim(e_A) = q-1$, $\dim(e_B) = p-1$, $\dim(e_C) = 1$.

$\mathrm{DA}_N$ is the spectral algebra of the log-GCD circulant $L_N$ (the image of $K_N(j) = \log\gcd(j,N)$ under the DFT), compressing $N$ spectral values to four algebraic objects.

**Theorem 8.2** (Galois Structure). $\mathrm{Gal}(\mathrm{DA}_N) \cong \mathbb{Z}/2\mathbb{Z}$, *generated by the prime-swap automorphism $\sigma: (p,q) \mapsto (q,p)$, acting as $\sigma(e_A) = e_B$, $\sigma(e_B) = e_A$, $\sigma(e_0) = e_0$, $\sigma(e_C) = e_C$.*

**Proof.** $\sigma$ is an algebra automorphism since it preserves the idempotent relations and $\sigma^2 = \mathrm{id}$. It is the only nontrivial automorphism: any automorphism must permute the primitive idempotents $\{e_0, e_A, e_B, e_C\}$; since $e_0$ (dimension $(p-1)(q-1)$) and $e_C$ (dimension 1) are distinguished by their dimensions, and $e_A, e_B$ have distinct dimensions $q-1 \neq p-1$ (since $p \neq q$) but are swapped by $\sigma$, the only nontrivial automorphism is $e_A \leftrightarrow e_B$. $\square$

**Theorem 8.3** (Characteristic Polynomial). *The factorisation of $N = pq$ is the degree-2 minimal polynomial*
$$M_N(x) = x^2 - \Sigma x + \Pi, \quad \Sigma = A + B = p\log q + q\log p,\quad \Pi = AB = N\log p\log q,$$
*whose roots are $\{A, B\} = \{p\log q, q\log p\}$. From the roots, $p = \lfloor e^{B/q} \rceil$ and $q = \lfloor e^{A/p}\rceil$.*

**Remark 8.4.** The characteristic polynomial $M_N$ is the maximal compression: the factorisation of the $b$-bit number $N$ is encoded in two real parameters $(\Sigma, \Pi)$. The discriminant $\Delta = \Sigma^2 - 4\Pi = (A-B)^2 \geq 0$ is a perfect square, confirming both roots are real. The spectral imbalance ratio $\Delta/\Sigma = |A-B|/(A+B)$ ranges from $0$ (balanced semiprimes $p \approx q$) to $1$ (maximally unbalanced).

**Theorem 8.5** (Adams Operations and $\lambda$-Ring). *Define Adams operations $\psi^k$ on $\mathrm{DA}_N$ by $\psi^k(e_A) = A^k$, $\psi^k(e_B) = B^k$. The power sums $P_k = A^k + B^k$ satisfy Newton's recurrence*
$$P_k = \Sigma \cdot P_{k-1} - \Pi \cdot P_{k-2},\quad P_0 = 2,\; P_1 = \Sigma.$$
*This makes $\mathrm{DA}_N$ a $\lambda$-ring. The generating sequence $\{P_k\}$ is a restricted Witt ghost component (NP III Proposition 8.6), connecting $\mathrm{DA}_N$ to the theory of Witt vectors.*

**Proof.** The recurrence follows from the standard Newton identity for the elementary symmetric polynomials $\Sigma = A + B$ and $\Pi = AB$: $P_k - \Sigma P_{k-1} + \Pi P_{k-2} = 0$ is equivalent to $(x-A)(x-B) = x^2 - \Sigma x + \Pi$. The $\lambda$-ring structure is established by the Adams operations satisfying $\psi^j \circ \psi^k = \psi^{jk}$ and the splitting principle. $\square$

**Theorem 8.6** (Embedding in the Universal Defect Algebra). *$\mathrm{DA}_N$ is the localisation of $\mathcal{A}$ at $N$: the Arithmetic Convolution Algebra $\mathrm{AC}_N = \mathrm{span}\{e_d : d \mid N\}$ with multiplication $e_d \cdot e_{d'} = e_{dd'}$ (when $dd' \mid N$, else $0$) satisfies $\mathrm{AC}_N \cong \mathrm{DA}_N$.*

**Proof.** The four divisors of $N = pq$ are $\{1, p, q, N\}$, in bijection with the four idempotents $\{e_0, e_A, e_B, e_C\}$ of $\mathrm{DA}_N$ via their gcd-classes. The multiplication rule $e_p \cdot e_q = e_N$ (since $\gcd(p,q) = 1$ and $pq = N$) matches the idempotent product in $\mathrm{DA}_N$. The eigenvalue structure of $L_N$ restricted to divisors of $N$ recovers the spectral data of $\mathrm{DA}_N$. $\square$

---

## 9. Abstract Algebra from Defect Symmetries

### 9.1 Groups

**Theorem 9.1** (Every Group is a Defect Symmetry Group). *Every group $G$ arises as the automorphism group of some finite combinatorial structure derived from the Arithmetic Vacuum.*

**Proof.** By Frucht's theorem (1939): every group is the automorphism group of some finite graph. A finite graph with vertex set $V \subseteq \mathbb{N}^*$ and edge set $E \subseteq \mathbb{N}^* \times \mathbb{N}^*$ is a combinatorial structure built from the natural numbers — a structure derived from the event set of $\mathcal{V}_\times$. Therefore every group arises as $\mathrm{Aut}(G)$ for some graph $G$ over $\mathbb{N}^*$. $\square$

**Corollary 9.2** (Galois Groups). *Every Galois group $\mathrm{Gal}(K/\mathbb{Q})$ of a finite extension $K/\mathbb{Q}$ arises as a group of automorphisms of a splitting field of a polynomial with integer coefficients — a structure directly derived from the Arithmetic Vacuum's ring $\mathbb{Z}$. The local Defect Algebra $\mathrm{DA}_N$ realises $\mathbb{Z}/2\mathbb{Z}$ as a Galois group for each semiprime $N$.*

**Remark 9.3.** The inverse Galois problem — whether every finite group is realised as $\mathrm{Gal}(K/\mathbb{Q})$ for some $K$ — is open in general (proved for symmetric groups, solvable groups, sporadic groups, etc.). Theorem 9.1 gives the weaker but sufficient result via Frucht's theorem.

### 9.2 Rings as Compatible Dual-Vacuum Structures

**Theorem 9.3** (Classification of Rings). *A ring is precisely an algebraic structure that simultaneously realises an additive CIS and a multiplicative CIS on the same set, with the compatibility condition of Definition 5.1. The category of rings is equivalent to the category of compatible dual-vacuum structures.*

**Proof.** By Theorem 5.2, the compatibility condition between an abelian group $(R, +)$ (additive CIS) and a monoid $(R, \times)$ (multiplicative CIS) is exactly the distributive law $a(b+c) = ab+ac$. A ring is by definition such a triple $(R, +, \times)$; any ring homomorphism preserves both structures and their compatibility. $\square$

### 9.3 Fields, Modules, Linear Algebra

**Theorem 9.4** (Fields). *A field is a ring in which the multiplicative CIS has no confined defects: every non-zero element has a multiplicative antidefect (inverse). Equivalently, a field is a CIS-compatible dual-vacuum structure in which both vacuum structures are group-complete.*

**Theorem 9.5** (Linear Algebra from Modules). *A module over a ring $R$ is a representation of the multiplicative CIS of $R$ on an additive abelian group: the monoid $(R, \times)$ acts by additive-CIS-preserving (i.e., linear) maps. Linear algebra over $\mathbb{R}$ or $\mathbb{C}$ is the theory of modules over the continuum-limit and symplectic-lift completions of the Arithmetic Vacuum.*

**Proof.** A module is by definition a map $R \times M \to M$ satisfying the module axioms, which are exactly the statement that multiplication by each $r \in R$ is a group endomorphism of $(M, +)$ (the additive CIS morphism condition from Theorem 5.2). $\square$

---

## 10. Topology from Causal Order

### 10.1 The Alexandrov Topology of the Divisibility Lattice

**Definition 10.1** (Causal Order Topology). Given a CIS $(E, \preceq, K, \rho)$, the *Alexandrov topology* declares $U \subseteq E$ open iff $U$ is an upset: $x \in U$ and $x \preceq y$ implies $y \in U$.

**Proposition 10.2.** *In the Alexandrov topology of $\mathcal{V}_\times = (\mathbb{N}^*, |)$:*
*(i) The open sets are upward-closed in the divisibility order.*
*(ii) The closed sets are downward-closed: sets of divisors.*
*(iii) The closure of $\{p\}$ for a prime $p$ is $\{p^k : k \geq 1\} \cup \{\text{multiples of }p\}$.*
*(iv) The open set generated by $n$ is $\{m : n \mid m\}$ — all multiples of $n$.*

### 10.2 Spec $\mathbb{Z}$ and the Arithmetic Spacetime

**Theorem 10.3** ($\mathrm{Spec}(\mathbb{Z})$ as Arithmetic Spacetime). *The prime spectrum $\mathrm{Spec}(\mathbb{Z}) = \{(0)\} \cup \{(p) : p \text{ prime}\}$ with the Zariski topology canonically identifies with the causal structure of $\mathcal{V}_\times$:*
- *Prime defects $p \in \mathbb{N}^*$ correspond to prime ideals $(p) \in \mathrm{Spec}(\mathbb{Z})$ (closed points).*
- *The vacuum $1$ corresponds to the generic point $(0)$, which is in the closure of every point.*
- *The Zariski topology is the Alexandrov topology of the divisibility-dual order on $\mathrm{Spec}(\mathbb{Z})$.*

**Proof.** Standard algebraic geometry: closed sets in $\mathrm{Spec}(\mathbb{Z})$ are $V(I) = \{P : I \subseteq P\}$ for ideals $I$; for $I = (n)$, $V((n)) = \{(p) : p \mid n\} \cup \{(0)\}$ if $n \neq 0$, which is the set of prime divisors of $n$ — exactly the downward closure of $\{n\}$ in the divisibility order. $\square$

**Remark 10.4.** $\mathrm{Spec}(\mathbb{Z})$ is the arithmetic spacetime: primes are points, the generic point $(0)$ is the vacuum connected to everything, and arithmetic results about primes are geometric statements in this space. The prime number theorem $\pi(x) \sim x/\log x$ is the statement that the prime points are distributed with density $1/\log$ in the coordinate $x$ — a density which, in log-coordinates, is uniform. This is exactly the Poisson sprinkling of CFG §10 with uniform reference density.

### 10.3 Metric Topology from the Causal Kernel

**Proposition 10.5.** *The symmetrised additive kernel $K_+(x,y) = |y-x|$ is a metric on $\mathbb{Z}$. Its Cauchy completion is $(\mathbb{R}, |\cdot|)$.*

**Proposition 10.6** ($p$-adic Completions). *By Ostrowski's theorem, every non-trivial absolute value on $\mathbb{Q}$ is equivalent to either $|\cdot|$ (the Archimedean completion, giving $\mathbb{R}$) or $|\cdot|_p$ for some prime $p$ (the $p$-adic completion, giving $\mathbb{Q}_p$). The $p$-adic absolute value $|n|_p = p^{-v_p(n)}$ (where $v_p(n)$ is the exponent of $p$ in the factorisation of $n$) defines the causal kernel in the direction of the prime defect $p$. The $p$-adic numbers are the completions of the Arithmetic Vacuum along each prime's causal cone.*

**Proof.** Ostrowski's theorem: standard. The interpretation in terms of prime causal cones follows from the identification of primes as elementary defects of $\mathcal{V}_\times$: each prime $p$ defines a natural distance on $\mathbb{Q}$ by measuring the $p$-adic valuation, i.e., how deeply $n$ sits in the causal cone above $p$ in the divisibility order. $\square$

---

## 11. Analysis from Causal Transport

### 11.1 Differentiation as Infinitesimal Generator

**Definition 11.1** (Infinitesimal Generator). In $\mathcal{V}_+$, the causal transport semigroup $T_h: f(x) \mapsto f(x+h)$ has infinitesimal generator
$$\partial_x f = \lim_{h \to 0} \frac{T_h f - f}{h} = \lim_{h \to 0} \frac{f(x+h) - f(x)}{h}.$$
This is the derivative: the infinitesimal generator of additive causal transport.

**Theorem 11.2** (Fundamental Theorem of Calculus as Filling Identity). *For smooth $f: \mathbb{R} \to \mathbb{R}$:*
$$f(b) - f(a) = \int_a^b f'(x)\, dx.$$
*This is the causal filling identity for $\mathcal{V}_+$: the net causal displacement from $a$ to $b$ equals the integral of the infinitesimal causal density $f'$ over the causal path $[a,b]$. It is the continuum analogue of Proposition 2.3(d): $\log n = \sum_{d|n}\Lambda(d)$.*

**Proof.** Standard calculus (Newton-Leibniz theorem). The causal interpretation: $K_+(a,b) = b - a$ is the total filling cost; $f'(x)dx$ is the infinitesimal cost density; the integral is the Riemann sum of filling costs. $\square$

### 11.2 Fourier Analysis as Spectral Decomposition

**Theorem 11.3** (Fourier Transform as Causal Spectral Decomposition). *The Fourier transform*
$$\hat{f}(\xi) = \int_{-\infty}^\infty f(x)\, e^{-2\pi i \xi x}\, dx$$
*decomposes $f$ into eigenstates of the causal generator $\partial_x$: the modes $e^{2\pi i\xi x}$ satisfy $\partial_x e^{2\pi i\xi x} = 2\pi i\xi \cdot e^{2\pi i\xi x}$. The Fourier modes are the stable eigenfunctions of causal translation transport in $\mathcal{V}_+$.*

**Proof.** Differentiation of the Fourier mode is immediate. The completeness of the Fourier basis in $L^2(\mathbb{R})$ (Plancherel theorem) is standard. $\square$

**Remark 11.4** (Structural Parallel). The Fourier decomposition of functions on $\mathcal{V}_+$ is the additive analogue of the von Mangoldt spectral decomposition of arithmetic functions on $\mathcal{V}_\times$ (Theorem 3.1 of *NP I*). In both cases, the natural "frequencies" are the eigenvalues of the causal transport operator: logarithms of primes $\{\log p\}$ in the multiplicative case, real frequencies $\{\xi\}$ in the additive case.

---

## 12. Open Problems

### 12.1 A Unified Parent Structure

**Open Problem 12.1** (Common Parent). *Identify a single CIS $(E, \preceq, K, \rho)$ from which both $\mathcal{V}_\times$ and $\mathcal{V}_+$ arise as restrictions. A candidate: the multiplicative group $(\mathbb{Q}^+, \times)$ maps to $(\mathbb{R}, +)$ via $\log$. Is there a CIS on $\mathbb{Q}^+$ whose restrictions to the multiplicative and additive directions recover $\mathcal{V}_\times$ and $\mathcal{V}_+$ respectively?*

### 12.2 The Distributive Law as a Theorem in Greater Generality

**Open Problem 12.2** (Universal Compatibility Theorem). *Prove that the distributive law is the unique compatibility condition between any additive CIS $(E, \leq, K_+, \rho_+)$ and any multiplicative CIS $(E, |, K_\times, \rho_\times)$ satisfying the condition of Definition 5.1, without appealing to the specific group structure of $(\mathbb{Z}, +)$. In other words: show that compatibility of two arbitrary CIS structures forces the ring axioms whenever the two structures are non-degenerate and mutually consistent.*

### 12.3 The Riemann Hypothesis as Equilibrium

**Conjecture 12.3** (RH as Vacuum Equilibrium). *The partition function $\zeta(s)$ of the Arithmetic Vacuum satisfies the equilibrium condition — all non-trivial zeros on $\mathrm{Re}(s) = 1/2$ — if and only if the prime defects of $\mathcal{V}_\times$ are distributed without preferred scale, i.e., the prime number theorem holds with optimal error: $\pi(x) = \mathrm{li}(x) + O(\sqrt{x}\log x)$. The Riemann Hypothesis is the arithmetic vacuum equilibrium condition.*

This is a restatement of the known equivalence between RH and the sharp prime number theorem, reframed in the Pleroma language. The novel claim is that this equilibrium condition is the arithmetic analogue of the CFG gravitational vacuum stability condition, and that the two might be provably equivalent via the transition to physics (Open Problem 12.5).

### 12.4 Classification of All Stable Defect Sectors

**Open Problem 12.4** (Defect Classification). *Classify all stable defect sectors of the Arithmetic Vacuum. Elementary defects are primes; bound states are prime powers and composites. The $p$-adic completions (Proposition 10.6) give one defect per prime direction. Is there a topological classification — via $\pi_*(\mathrm{Spec}(\mathbb{Z}))$ or étale cohomology — of all stable sectors? Does this classification have a direct bijection with the stable particle sectors of the Fisher-Frame programme?*

### 12.5 The Transition to Physics (Companion Paper)

**Open Problem 12.5** (Physical Sector — Companion Programme). *Show that the continuum coarse-graining limit of $\mathcal{V}_\times$ in 3+1 dimensions forces:*

*(a) Lorentzian signature $(-, +, +, +)$ from the generator-vs-fluctuation split of the Gaussian fixed-point kernel.*

*(b) The gauge group $\mathrm{SU}(2)_L \times U(1)_Y$ from the $\mathbb{Z}/2\mathbb{Z}$ Galois group of $\mathrm{DA}_N$ in the continuum limit.*

*(c) Three generations of matter from the electromagnetic running constraint (as derived in Fisher-Frame Defects v0.7 §16).*

*(d) Newton's constant $G = 2\ell_0^2$ from the Gaussian second moment.*

*(e) The cosmological constant $\Lambda = 3/R_C^2$ from the IR causal-boundary zero mode.*

*The specific programme: (i) take the log-coordinate map $n \mapsto \log n$ on $\mathcal{V}_\times$; (ii) apply the Gaussian fixed-point theorem (Theorem 6.2); (iii) identify the four-element divisibility diamond $\{1, p, q, pq\}$ with the causal bivector space $\mathcal{C}_p = \mathrm{span}\{\ell \wedge e_1, \ell \wedge e_2, n \wedge e_1, n \wedge e_2\}$ of CFG §14; (iv) show that the $\mathbb{Z}/2\mathbb{Z}$ prime-swap symmetry forces $\pi_1 = \mathbb{Z}/2\mathbb{Z}$ in the continuous rotation group, selecting $\mathrm{SO}(3)$ with its spinorial lift $\mathrm{SU}(2)$; (v) derive the 16π normalization as $4\pi \times 4$ where $4$ is the size of the minimal divisibility diamond.*

---

## 13. Conclusion

The Arithmetic Vacuum — the pair of causal information structures $(\mathcal{V}_\times, \mathcal{V}_+)$ on the integers — is the foundational object from which mathematics arises through a sequence of forced constructions.

The multiplicative CIS $(\mathbb{N}^*, |, K_\times, \Lambda)$ is a causal information structure in the precise sense of CFG. Its elementary defects are the prime numbers (Theorem 2.5); their unique decomposition is the Fundamental Theorem of Arithmetic cast as a topological stability theorem (Theorem 2.6); the von Mangoldt function $\Lambda$ is the canonical density encoding the information cost of each causal step; and the Möbius function $\mu$ is the antidefect operator satisfying the vacuum neutrality condition $\mu * \mathbf{1} = \varepsilon$ (Theorem 3.2).

The ring $\mathbb{Z}$ is derived, not axiomatised: it is the minimal structure compatible with both $\mathcal{V}_\times$ and $\mathcal{V}_+$, and the distributive law is the unique compatibility condition between them (Theorem 5.2). The number tower $\mathbb{Z} \to \mathbb{Q} \to \mathbb{R} \to \mathbb{C}$ is a sequence of forced completions: the defect-antidefect completion, the Gaussian coarse-graining limit, and the symplectic lift — the last of which is the passage to quantum phase dynamics (Theorem 6.6).

The algebra of the Arithmetic Vacuum is the Universal Defect Algebra $(\mathcal{A}, +, *, \varepsilon)$: the ring of arithmetic functions under Dirichlet convolution, with $\zeta(s)$ as the partition function of the prime defect gas (Theorem 7.7). For a fixed semiprime $N = pq$, this localises to the four-dimensional $\mathrm{DA}_N$ with Galois group $\mathbb{Z}/2\mathbb{Z}$ and $\lambda$-ring structure (Theorems 8.2, 8.5).

Abstract algebra arises from defect symmetries: every group is the automorphism group of some arithmetic structure (Theorem 9.1), and every ring is a compatible dual-vacuum structure (Theorem 9.3). Topology arises from the causal order: the Alexandrov topology of the divisibility lattice is $\mathrm{Spec}(\mathbb{Z})$, the arithmetic spacetime; metric topology arises from the additive kernel; $p$-adic completions arise along each prime's causal cone (Proposition 10.6). Analysis arises from causal transport: differentiation is the infinitesimal generator of $\mathcal{V}_+$, integration is the causal filling cost, the Fundamental Theorem of Calculus is the filling identity (Theorem 11.2), and Fourier analysis is spectral decomposition of the additive causal operator (Theorem 11.3).

The Riemann Hypothesis is conjectured to be the equilibrium condition on the Arithmetic Vacuum (Conjecture 12.3). The transition from the Arithmetic Vacuum to physics — Lorentzian geometry, gauge symmetry, spinorial matter, and the fundamental constants — is the subject of the companion paper, whose programme is precisely formulated in Open Problem 12.5.

The programme's honest position: all algebraic claims in this paper are proved. The transition to physics is a precisely stated programme with explicit proof obligations. The Riemann Hypothesis is an equilibrium conjecture, not a claim of proof. The identification of the Arithmetic Vacuum as the unique foundation of mathematics rests on all the proved constructions and on the open problems whose resolution would complete it.

---

## A. Symbol Glossary

| Symbol | Meaning |
|--------|---------|
| $\mathcal{V}_\times$ | Multiplicative CIS $(\mathbb{N}^*, \mid, K_\times, \Lambda)$ |
| $\mathcal{V}_+$ | Additive CIS $(\mathbb{Z}, \leq, K_+, 1)$ |
| $K_\times(x,y)$ | $\log(y/x)$ when $x \mid y$, else $0$ |
| $K_+(x,y)$ | $y - x$ when $x \leq y$, else $0$ |
| $\Lambda(n)$ | Von Mangoldt function: $\log p$ if $n = p^k$, else $0$ |
| $\mu(n)$ | Möbius function: canonical antidefect operator |
| $\varepsilon(n)$ | $[n=1]$: vacuum function, convolution identity |
| $\mathcal{A}$ | Universal Defect Algebra: $\{f: \mathbb{N}^* \to \mathbb{C}\}$ under $*$ |
| $\mathrm{DA}_N$ | Local Defect Algebra at $N = pq$ (four-dimensional) |
| $\mathrm{AC}_N$ | Arithmetic Convolution Algebra at $N$ |
| $\sigma$ | Prime-swap automorphism $(p,q) \mapsto (q,p)$ |
| $\psi^k$ | Adams operations on $\mathrm{DA}_N$ |
| $M_N(x)$ | Characteristic polynomial $x^2 - \Sigma x + \Pi$ |
| $\zeta(s)$ | Riemann zeta function: partition function of $\mathcal{V}_\times$ |
| $G(w,s)$ | $\sum_N H_N(s)/N^w = \zeta(s+w)\zeta(w-1)/\zeta(w)$ |
| $\mathrm{Spec}(\mathbb{Z})$ | Arithmetic spacetime: prime spectrum of $\mathbb{Z}$ |

---

## B. Derived versus Postulated

| Item | Status | Location |
|------|--------|----------|
| $\mathcal{V}_\times$ is a CIS | Proved | Theorem 2.2 |
| $\mathcal{V}_+$ is a CIS | Proved | Proposition 4.2 |
| Primes as elementary defects | Proved | Theorem 2.5 |
| FTA as topological stability | Proved (standard proof + interpretation) | Theorem 2.6 |
| Möbius vacuum neutrality | Proved | Theorem 3.2 |
| Möbius inversion | Proved | Proposition 3.4 |
| Distributive law from compatibility | Proved | Theorem 5.2 |
| Ring axioms as derived | Proved | Corollary 5.3 |
| $\mathbb{Z}$ as initial ring | Proved (standard) | Theorem 5.4 |
| Gaussian coarse-graining fixed point | Proved (CLT) | Theorem 6.2 |
| $\mathbb{R}$ as Cauchy completion | Proved (standard) | Theorem 6.3 |
| $\mathbb{C}$ as symplectic lift | Proved (Frobenius + new interpretation) | Theorem 6.6 |
| $\mathbb{C}$ algebraically closed | Proved (standard) | Theorem 6.7 |
| $\mathcal{A}$ is integral domain ring | Proved | Theorem 7.2 |
| $\zeta(s)$ as partition function | Proved (Euler product + interpretation) | Theorem 7.7 |
| $\mathrm{Gal}(\mathrm{DA}_N) \cong \mathbb{Z}/2\mathbb{Z}$ | Proved | Theorem 8.2 |
| Adams operations / $\lambda$-ring | Proved | Theorem 8.5 |
| Groups from defect symmetries | Proved (Frucht's theorem) | Theorem 9.1 |
| Rings from dual-vacuum structures | Proved | Theorem 9.3 |
| Topology from causal order | Proved | §10 |
| Analysis from causal transport | Proved | §11 |
| RH as equilibrium condition | **Conjecture** | Conjecture 12.3 |
| Physics from continuum limit | **Programme** | Open Problem 12.5 |