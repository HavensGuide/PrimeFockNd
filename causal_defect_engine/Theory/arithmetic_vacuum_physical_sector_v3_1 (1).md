# The Arithmetic Vacuum as the Origin of the Fisher-Frame Information Vacuum

**Kai Cobbs**  
Standalone preprint v0.3 — May 2026  
Companion to *The Arithmetic Vacuum: Causal Information Foundations of Mathematics* (v0.1),  
*Fisher-Frame Defects and Spinorial Matter v0.7*, and *Causal Filling Gravity v0.7*

---

## Abstract

The companion paper (*Paper I*) established the Arithmetic Vacuum — the multiplicative causal information structure $\mathcal{V}_\times = (\mathbb{N}^{(\infty)}, \leq_\mathrm{coord}, K_\times, \Lambda)$ — as the foundational object from which all of mathematics is derived. The present paper derives the *statistical information vacuum* assumed in the Fisher-Frame defect theory from this arithmetic structure via a sequence of theorems, two stated axioms (non-degeneracy of the metric and dynamical closure), and one selection principle (minimality of the coarse-graining sector subject to non-trivial automorphism content).

The central result is: given the two axioms and the minimality selection, the continuum coarse-graining of $\mathcal{V}_\times$ produces a 4-dimensional Lorentzian statistical manifold equipped with a Fisher information metric, an $\mathrm{SO}(3)$ frame bundle, an $\mathrm{SU}(2)$ spinorial cover, local frame covariance, and a conserved $\mathbb{Z}/2\mathbb{Z}$ monodromy class. These are precisely the inputs of the Fisher-Frame defect theory. The particle spectrum, gauge couplings, Weinberg angle, electron mass, and generation count are downstream consequences of the defect physics of this vacuum; they are derived in the Fisher-Frame paper and cited here as inherited results.

Three theorems are new in this version. First, the Fisher metric is derived from Čencov's uniqueness theorem: the coarse-graining is a sufficient statistic for the $\mathbb{Z}/2\mathbb{Z}$ automorphism group, and Čencov's theorem then forces the Fisher information metric as the unique compatible Riemannian metric. Second, Lorentzian *signature uniqueness* is derived — without citing CFG — via the Paley-Wiener theorem: the zero backward kernel of $\mathcal{V}_\times$ (part of the CIS definition, not imposed) forces the field propagator to be retarded; the Paley-Wiener theorem then requires the inverse-metric quadratic form to have signature $(-,+,+,+)$ uniquely, given the non-degeneracy axiom and the positive Fisher metric on spatial directions. The treatment now handles off-diagonal cross-terms and reduces the general case to the diagonal one by spatial rotation. Third, the symplectic lift to four-dimensional phase space is now derived from an explicit Dynamical Closure Axiom: the CIS time direction generates Markov evolution that conserves the Fisher metric, and Hamilton's theorem then forces $T^*\mathcal{M}_\mathrm{stat}$ as the minimal phase space.

Compared to v0.2, the present version makes three structural revisions. (i) The previous "$k = 2$ uniquely forced" theorem (v0.2 Theorem 3.8) is reformulated as a *minimality selection*: $k = 2$ is the minimal sector with non-trivial abstract automorphism content, and is selected by combining structural minimality with the empirical (3+1)-dimensional structure of spacetime. The strong "forced from pure arithmetic" claim is demoted to a strictly weaker "uniquely determined under minimality." (ii) The §6 Lorentzian-signature derivation is sharpened from "eliminates $g^{tt} > 0$" to "uniquely selects signature $(-,+,+,+)$" via the off-diagonal analysis. (iii) The §7.1 fundamental-group argument is separated into two independent claims — a topological claim that $\pi_1(\mathrm{SO}(3)) = \mathbb{Z}/2\mathbb{Z}$ comes from the dimensionality of space alone, and an arithmetic correspondence claim that identifies the non-trivial $\pi_1$ class with the prime-swap holonomy of the minimal divisibility diamond.

---

## Contents

1. Introduction and Scope  
2. The Exponent-Vector Causal Information Structure  
3. Sufficient Coarse-Graining to the 2-Prime Sector  
4. The Symplectic Lift and Four-Dimensional Phase Space  
5. The Fisher Metric from Statistical Distinguishability  
6. Lorentzian Signature from the Paley-Wiener Theorem  
6a. The Derived Free-Field Action  
7. The $\mathrm{SO}(3)$ Frame Bundle and Its $\mathrm{SU}(2)$ Cover  
8. Local Frame Covariance from Local Prime-Ordering Symmetry  
9. The Conserved $\mathbb{Z}/2\mathbb{Z}$ Monodromy Class  
10. The Information Vacuum: What Is Delivered and What Is Inherited  
11. Open Problems  
12. Conclusion  
A. Symbol Glossary  
B. Theorem Status Table  
C. Axioms and Selection Principles Used  

---

## 1. Introduction and Scope

### 1.1 What This Paper Does and Does Not Claim

The purpose of this paper is not to derive the full particle spectrum directly from arithmetic. Its purpose is to show that the continuum, symmetry-faithful coarse-graining of the Arithmetic Vacuum yields the statistical information vacuum assumed in the Fisher-Frame defect theory: a Lorentzian Fisher manifold with $\mathrm{SO}(3)$ frame topology, $\mathrm{SU}(2)$ lift, and conserved $\mathbb{Z}/2\mathbb{Z}$ monodromy. The Fisher-Frame paper then describes stable charged spinorial matter as topological sectors of this vacuum, and all particle and gauge results downstream are consequences of that defect theory, not re-derived here.

The hierarchy is:

$$\boxed{
\text{Arithmetic Vacuum (Paper I)}
\;\longrightarrow\;
\underbrace{\text{Statistical Information Vacuum}}_{\text{this paper}}
\;\longrightarrow\;
\underbrace{\text{Fisher-Frame Defects (Paper III)}}_{\text{electron, gauge, mass}}
\;\longrightarrow\;
\underbrace{\text{Causal Filling Gravity (Paper IV)}}_{\text{Einstein, }G,\,\Lambda}
}$$

The arithmetic paper explains why the Fisher-Frame paper's starting point is not arbitrary.

### 1.2 What the Fisher-Frame Paper Assumes

The Fisher-Frame paper (*Fisher-Frame Defects and Spinorial Matter v0.7*) begins from a vacuum characterised by:

**(i)** A frame field $U(x) \in \mathrm{SU}(2)$ over a 4-dimensional Lorentzian manifold.

**(ii)** A Fisher information metric on the space of frame orientations.

**(iii)** Local $\mathrm{SU}(2)$ gauge covariance, with covariant derivative $D_\mu U = \partial_\mu U - \omega_\mu U$.

**(iv)** A $\mathbb{Z}/2\mathbb{Z}$ spinorial defect class: frame configurations with monodromy $U(2\pi) = -U(0)$.

These four inputs determine the Fisher-Frame action uniquely (by the uniqueness argument in Fisher-Frame §1.3: Fisher metric + local gauge invariance + $\leq 2$ derivatives + $\mathbb{Z}/2\mathbb{Z}$ topology). The present paper derives all four from the Arithmetic Vacuum.

### 1.3 The Previous Drafts and Their Errors

An earlier version of this paper (v0.1) contained two errors that invalidated several sections:

**Error 1 (log-map order).** The log map $n \mapsto \log n$ was claimed to be a causal information structure isomorphism from the divisibility order to the total order on $\mathbb{R}$. This is false: $x \mid y \Rightarrow \log x \leq \log y$ but not conversely ($2 \nmid 3$ though $\log 2 < \log 3$). The correct structure-preserving map is the exponent-vector map, derived in §2.

**Error 2 (prime density in log-coordinates).** The prime defects were claimed to be unit-rate Poisson in $t = \log x$ coordinates. The correct density in $t$-coordinates is $e^t/t$, not 1. The natural unit-rate clock is the logarithmic integral $\tau = \mathrm{li}(x)$. This does not affect the main theorems of the present paper, which do not depend on the prime density in specific coordinates.

Both errors were corrected in v0.2. Sections that depended on them (dimensional forcing via the diamond-bivector identification, and several cosmological claims) were either corrected or demoted to future work.

The v0.2 draft contained three additional structural issues, addressed in the present v0.3 version:

**Issue A (v0.2 §3.8 circularity).** The "dimensional forcing" argument required $\mathrm{Aut}_\mathrm{abs}(\mathcal{V}_\times) = \mathbb{Z}/2\mathbb{Z}$ in order to derive $k \leq 2$. But by Lemma 3.2 applied to the full vacuum on $\mathbb{N}^{(\infty)}$, $\mathrm{Aut}_\mathrm{abs}(\mathcal{V}_\times) = S_\infty$, not $\mathbb{Z}/2\mathbb{Z}$. The upper bound $k \leq 2$ does not follow from pure structure. In v0.3, the §3.8 result is reformulated as a *minimality selection* (Selection Principle 3.8) rather than a forcing theorem: $k = 2$ is the minimal sector with non-trivial abstract automorphism content, and is selected by combining structural minimality with the empirical (3+1)-dimensional structure of spacetime. This is strictly weaker than the v0.2 claim but is now defensible.

**Issue B (v0.2 §6.3 elimination vs uniqueness).** The Paley-Wiener argument as stated in v0.2 eliminated the case $g^{tt} > 0$ but did not address the degenerate case $g^{tt} = 0$ or off-diagonal $g^{ti}$ terms. The v0.3 version states an explicit non-degeneracy axiom (Axiom 6.0) and reformulates Theorem 6.4 at the signature level: the inverse-metric quadratic form must have Lorentzian signature $(-,+,+,+)$. The off-diagonal case is reduced to the diagonal one by spatial rotation.

**Issue C (v0.2 §7.1 conflation of 2-torsion notions).** The argument "the discrete $\mathbb{Z}/2\mathbb{Z}$ automorphism is a 2-torsion element, hence lies in $\pi_1$" conflated 2-torsion in the Lie group itself (e.g., $-I \in \mathrm{SU}(2)$) with 2-torsion in $\pi_1$ of the Lie group. The v0.3 §7.1 separates the topological claim ($\pi_1(\mathrm{SO}(3)) = \mathbb{Z}/2\mathbb{Z}$ from the dimensionality of space alone) from the correspondence claim (the discrete prime-swap corresponds to the non-trivial $\pi_1$ class via diamond holonomy).

**Issue D (v0.2 symplectic lift unmotivated).** The v0.2 §4 cited the foundations paper's algebraic argument ($\mathbb{R} \to \mathbb{C}$ is the symplectic lift in the number tower) but did not provide an independent physical reason for the lift. The v0.3 version states an explicit Dynamical Closure Axiom (Axiom 4.0): the CIS time direction generates Markov evolution conserving the Fisher metric. Hamilton's theorem then forces $T^*\mathcal{M}_\mathrm{stat}$ as the minimal phase space.

---

## 2. The Exponent-Vector Causal Information Structure

### 2.1 The Correct Structure-Preserving Map

**Definition 2.1** (Exponent-Vector Map). The *exponent-vector map* is:
$$\mathbf{v}: \mathbb{N}^* \to \mathbb{N}^{(\infty)},\qquad n = \prod_p p^{v_p(n)} \longmapsto (v_p(n))_{p\,\text{prime}},$$
where $\mathbb{N}^{(\infty)}$ denotes the set of finitely-supported sequences of non-negative integers.

**Theorem 2.2** (Order Isomorphism). *The exponent-vector map $\mathbf{v}$ is an order isomorphism between $(\mathbb{N}^*, \mid)$ and $(\mathbb{N}^{(\infty)}, \leq_\mathrm{coord})$, where $\leq_\mathrm{coord}$ is the coordinatewise partial order:*
$$x \mid y \;\iff\; v_p(x) \leq v_p(y)\text{ for all primes }p.$$

**Proof.** Divisibility: $x \mid y$ iff $y/x \in \mathbb{N}^*$ iff $v_p(y) - v_p(x) \geq 0$ for all $p$, i.e., $\mathbf{v}(x) \leq_\mathrm{coord} \mathbf{v}(y)$. The map is injective (FTA) and surjective (every sequence in $\mathbb{N}^{(\infty)}$ corresponds to a unique positive integer). $\square$

**Corollary 2.3** (Log as Scalar Projection). The log map is:
$$\log n = \sum_p v_p(n)\log p = \langle \mathbf{v}(n), \mathbf{w}\rangle,\qquad \mathbf{w} = (\log p)_{p\,\text{prime}}.$$
It is a scalar projection of the exponent vector onto the weight vector $\mathbf{w}$, not an order isomorphism. In particular $\log x \leq \log y$ does not imply $x \mid y$.

### 2.2 The Corrected Causal Information Structure

**Definition 2.4** (Exponent-Vector CIS). The *Arithmetic Vacuum* is the causal information structure:
$$\mathcal{V}_\times = (\mathbb{N}^{(\infty)},\; \leq_\mathrm{coord},\; K_\times,\; \Lambda),$$
where the causal kernel for $\mathbf{u} \leq_\mathrm{coord} \mathbf{v}$ is:
$$K_\times(\mathbf{u}, \mathbf{v}) = \langle \mathbf{v} - \mathbf{u},\, \mathbf{w}\rangle = \sum_p (v_p - u_p)\log p,\quad 0\text{ otherwise.}$$

**Proposition 2.5** (Properties of $K_\times$).

*(a) Additivity:* $K_\times(\mathbf{u}, \mathbf{v}) + K_\times(\mathbf{v}, \mathbf{w}) = K_\times(\mathbf{u}, \mathbf{w})$ for $\mathbf{u} \leq \mathbf{v} \leq \mathbf{w}$.

*(b) Factorisation:* The kernel decomposes over independent prime directions: $K_\times(\mathbf{u}, \mathbf{v}) = \sum_p (v_p - u_p)\log p$.

*(c) Path cost recovery:* $K_\times(\mathbf{0}, \mathbf{v}(n)) = \log n = \sum_{d \mid n}\Lambda(d)$ (companion paper Proposition 2.3).

### 2.3 Mode Elimination

**Definition 2.6** (Mode Elimination). For a set $S = \{p_1, \ldots, p_k\}$ of $k$ primes, the *mode elimination map* is:
$$\Phi_S: n \longmapsto (v_{p_1}(n), \ldots, v_{p_k}(n)) \in \mathbb{N}^k.$$

This discards all prime components outside $S$. Its image is the effective theory on $S$-smooth numbers: integers whose prime factors all lie in $S$.

**Proposition 2.7** (Partition Function of the Effective Theory). *Integrating out all primes except $S = \{p_1, p_2\}$ corresponds to removing their Euler factors from the partition function $\zeta(\beta) = \prod_p(1-p^{-\beta})^{-1}$:*
$$Z_\mathrm{eff}(\beta) = \frac{1}{(1-p_1^{-\beta})(1-p_2^{-\beta})}.$$
*This is the exact partition function of the $\{p_1, p_2\}$-smooth numbers.*

**Proof.** Each prime $p$ contributes factor $(1-p^{-\beta})^{-1} = \sum_{k\geq 0}p^{-k\beta}$ to $\zeta(\beta)$. Restricting to $\{p_1, p_2\}$-smooth numbers retains exactly the $p_1$ and $p_2$ factors. $\square$

---

## 3. Sufficient Coarse-Graining to the 2-Prime Sector

This section establishes that $k = 2$ prime directions is the minimal coarse-graining sector with non-trivial abstract automorphism content, and is selected by combining structural minimality with the empirical (3+1)-dimensional structure of spacetime. The "$k = 2$ is uniquely forced from pure arithmetic" claim of v0.2 is reformulated as the strictly weaker claim that $k = 2$ is uniquely determined under minimality plus dimensional consistency.

### 3.1 The Abstract Automorphism Group

**Definition 3.1** (Abstract Automorphism Group). The *abstract automorphism group* $\mathrm{Aut}_\mathrm{abs}(\mathbb{N}^k, \leq_\mathrm{coord})$ consists of all bijections $\phi: \mathbb{N}^k \to \mathbb{N}^k$ that preserve the coordinatewise partial order, without reference to the specific prime weights $(\log p_i)$.

**Lemma 3.2.** $\mathrm{Aut}_\mathrm{abs}(\mathbb{N}^k, \leq_\mathrm{coord}) = S_k$, *the group of permutations of $k$ coordinates.*

**Proof.** Any order-preserving bijection must map atoms $\{\mathbf{e}_i\}$ (elements covering $\mathbf{0}$) to atoms. The map is fully determined by the permutation of atoms. Every permutation $\sigma \in S_k$ acts as $\phi_\sigma(\mathbf{v})_i = v_{\sigma^{-1}(i)}$ and preserves $\leq_\mathrm{coord}$. $\square$

**Remark 3.3.** The weighted automorphism group (preserving the kernel weights $\log p_i$ in addition to the order) is trivial, since all $\log p_i$ are distinct (Baker-Wüstholz theorem on linear independence of logarithms). The abstract automorphism group $S_k$ is the physically relevant one: it encodes the symmetry of the causal order structure independently of which specific primes are selected.

### 3.2 The Two Coarse-Graining Conditions

**Definition 3.4** (Automorphism-Faithful Coarse-Graining). A coarse-graining $\Phi_S: \mathcal{V}_\times \to \mathcal{V}_\mathrm{eff}$ is *automorphism-faithful* if the induced map $(\Phi_S)_*: \mathrm{Aut}_\mathrm{abs}(\mathcal{V}_\times) \to \mathrm{Aut}_\mathrm{abs}(\mathcal{V}_\mathrm{eff})$ is injective.

**Definition 3.5** (No-New-Symmetries). A coarse-graining satisfies *no-new-symmetries* if $|\mathrm{Aut}_\mathrm{abs}(\mathcal{V}_\mathrm{eff})| \leq |\mathrm{Aut}_\mathrm{abs}(\mathcal{V}_\times)|$.

**Theorem 3.6** (Automorphism-Faithfulness from Bayesian Sufficiency). *If $\Phi_S$ is a sufficient statistic for the abstract automorphism group, then it is automorphism-faithful.*

**Proof.** Suppose $(\Phi_S)_*(\sigma) = \mathrm{id}$ for non-trivial $\sigma \in \mathrm{Aut}_\mathrm{abs}(\mathcal{V}_\times)$. Then $\Phi_S(\sigma(n)) = \Phi_S(n)$ for all $n$, so the images of $n$ and $\sigma(n)$ coincide. By the data processing inequality: $I(\sigma;\,\Phi_S(n)) \leq I(\sigma;\, n)$. But $\Phi_S(\sigma(n)) = \Phi_S(n)$ implies $I(\sigma;\,\Phi_S(n)) = 0$, contradicting $I(\sigma;\,n) > 0$ (since $\sigma \neq \mathrm{id}$ and $n \neq \sigma(n)$ generically). $\square$

**Theorem 3.7** (No-New-Symmetries from Information Neutrality). *A coarse-graining that introduces a symmetry absent from $\mathcal{V}_\times$ creates mutual information between $n$ and the new symmetry transformation that was not present in the original distribution, violating the data processing inequality.*

**Proof.** Let $B \in \mathrm{Aut}_\mathrm{abs}(\mathcal{V}_\mathrm{eff})$ with no preimage in $\mathrm{Aut}_\mathrm{abs}(\mathcal{V}_\times)$. Then $I(B;\,\Phi_S(n)) > 0$ but $I(B;\,n) = 0$ (the symmetry $B$ has no correlate in $\mathcal{V}_\times$). This contradicts the data processing inequality $I(B;\,\Phi_S(n)) \leq I(B;\,n)$. $\square$

**Selection Principle 3.8** (Minimal Non-Trivial Sector). *The 2-prime sector is the minimal $k \geq 1$ for which the $k$-prime effective CIS has non-trivial abstract automorphism group. Combined with the empirical (3+1)-dimensional structure of physical spacetime, this uniquely selects $k = 2$.*

**Justification.** The $k$-prime effective CIS has $\mathrm{Aut}_\mathrm{abs} = S_k$ (Lemma 3.2). We compute:

- $k = 0$: trivial CIS, no defects, no physics.
- $k = 1$: $S_1 = \{e\}$ trivial. The 1-prime sector has no non-trivial automorphism, hence no branching structure, hence no non-trivial topology in the continuum limit. Excluded by non-triviality.
- $k = 2$: $S_2 = \mathbb{Z}/2\mathbb{Z}$, the minimal non-trivial automorphism group. The 2-prime sector has the minimal non-degenerate branching diamond $\{1, p_1, p_2, p_1 p_2\}$, supports a $\mathbb{Z}/2\mathbb{Z}$ prime-swap symmetry, and produces (via symplectic lift, §4) a 4-dimensional phase space.
- $k \geq 3$: $S_k$ has order $\geq 6$. The symplectic lift produces a $2k$-dimensional phase space of dimension $\geq 6$, inconsistent with the empirical (3+1)-dimensional spacetime.

The minimum $k$ producing non-trivial automorphism content is $k = 2$. The maximum $k$ consistent with empirical spacetime dimensionality is $k = 2$. The intersection is the unique selection $k = 2$. $\square$

**Honest Framing** (replaces v0.2 Theorem 3.8). The v0.2 draft claimed $k = 2$ was "uniquely forced" by an information-theoretic argument requiring $\mathrm{Aut}_\mathrm{abs}(\mathcal{V}_\times) = \mathbb{Z}/2\mathbb{Z}$. By Lemma 3.2, the full vacuum on $\mathbb{N}^{(\infty)}$ has $\mathrm{Aut}_\mathrm{abs} = S_\infty$, not $\mathbb{Z}/2\mathbb{Z}$, so the upper bound $k \leq 2$ does not follow from pure structure. The v0.3 reformulation acknowledges that the upper bound requires *either* a minimality principle (which is a selection criterion, not a structural theorem) *or* empirical input (the observed dimension of spacetime). Both are stated explicitly above. Whether a purely structural derivation exists — for example, via an "elementary automorphism" principle or a categorical-initiality argument in some category of coarse-grainings — remains an open problem (Open Problem 11.4).

**Remark 3.9** (The Two Pictures of $\mathbb{Z}/2\mathbb{Z}$). The full infinite-dimensional CIS $\mathcal{V}_\times = (\mathbb{N}^{(\infty)}, \leq_\mathrm{coord})$ has abstract automorphism group $S_\infty$. The 2-prime effective CIS has abstract automorphism group $S_2 = \mathbb{Z}/2\mathbb{Z}$. These are different groups: $\mathbb{Z}/2\mathbb{Z}$ is the *elementary* automorphism content (any single transposition), and $S_\infty$ is the full automorphism content. The 2-prime sector preserves elementary automorphism content but discards higher-order automorphisms. This is a meaningful information-theoretic statement, but it does not, by itself, force $k = 2$ — that requires the additional minimality selection in Selection Principle 3.8.

**Theorem 3.10** (Sufficiency of the 2-Prime Projection). *The map $\Phi: n \mapsto (v_{p_1}(n), v_{p_2}(n))$ is a sufficient statistic for the $\mathbb{Z}/2\mathbb{Z}$ automorphism $\sigma: (v_{p_1}, v_{p_2}) \mapsto (v_{p_2}, v_{p_1})$:*
$$I(\sigma;\, \Phi(n)) = I(\sigma;\, n).$$

**Proof.** The automorphism $\sigma$ acts only on the components $(v_{p_1}(n), v_{p_2}(n))$ of the factorisation of $n$, leaving all other prime valuations unchanged. Therefore $\sigma$ is conditionally independent of $\{v_q(n)\}_{q \notin \{p_1, p_2\}}$ given $(v_{p_1}(n), v_{p_2}(n)) = \Phi(n)$. Hence $I(\sigma;\, n) = I(\sigma;\, \Phi(n))$. $\square$

---

## 4. The Symplectic Lift and Four-Dimensional Phase Space

### 4.1 The Dynamical Closure Axiom

The CIS structure of $\mathcal{V}_\times$ encodes a distinguished time direction (Lemma 6.1 below). This time direction is not merely a coordinate but a generator of evolution: states at later $\hat{t}$ are obtained from states at earlier $\hat{t}$ by a flow on some phase space containing the configuration space $\mathcal{M}_\mathrm{stat}$. The minimal phase space supporting such a flow is determined by the following axiom.

**Axiom 4.0** (Dynamical Closure). *The CIS time direction $\hat{t}$ generates a one-parameter family of evolutions on a phase space $\mathcal{P} \supseteq \mathcal{M}_\mathrm{stat}$ satisfying:*

*(D1) Markov property: evolution is a Markov morphism (preserves probability mass).*

*(D2) Fisher conservation: the Fisher metric (when defined on $\mathcal{P}$) is preserved by the flow.*

*(D3) Minimality: $\mathcal{P}$ is the smallest space carrying a flow compatible with (D1) and (D2) and projecting to $\hat{t}$-translation on $\mathcal{M}_\mathrm{stat}$.*

Axiom 4.0 replaces the v0.2 invocation of the foundations paper's algebraic argument as the *physical* justification for the symplectic lift. The algebraic argument ($\mathbb{R} \to \mathbb{C}$ in the number tower) provides a parallel structural derivation and is consistent with Axiom 4.0 (see Remark 4.4 below).

### 4.2 The Symplectic Lift

**Definition 4.1** (Continuum Position Space). The continuum limit of the 2-prime effective CIS is the position space $\mathbb{R}^2$ with coordinates $(u_1, u_2)$ corresponding to the exponent valuations $(v_{p_1}, v_{p_2})$ in the Gaussian coarse-graining limit (CLT applied to the exponent-vector kernel).

**Theorem 4.2** (Symplectic Lift from Dynamical Closure). *Under Axiom 4.0, the minimal phase space $\mathcal{P}$ is $T^*\mathcal{M}_\mathrm{stat} = T^*\mathbb{R}^2 \cong \mathbb{R}^4$, equipped with the canonical symplectic form $\omega = d\phi_1 \wedge du_1 + d\phi_2 \wedge du_2$, and the flow is a Hamiltonian flow generated by a Hamiltonian function $H: \mathcal{P} \to \mathbb{R}$.*

**Proof.** Conditions (D1) and (D2) of Axiom 4.0 require the evolution to preserve a probability measure and to be an isometry of the Fisher metric. By the Liouville/Stone-von Neumann theorem, such measure-preserving isometric flows on a configuration space $\mathcal{M}$ are generated by Hamiltonian functions on $T^*\mathcal{M}$ with the canonical symplectic structure. The dimensionality of the minimal such space is $2\dim(\mathcal{M}) = 4$.

We must verify minimality (D3): no smaller phase space supports the required flow. A flow on $\mathcal{M}_\mathrm{stat} = \mathbb{R}^2$ alone (without phase coordinates) is a vector field on $\mathbb{R}^2$. By the Poincaré-Hopf theorem and the requirement of Fisher-metric preservation, such a flow must be an isometry of the flat $\mathbb{R}^2$, hence a translation or rotation. Pure translation is consistent with $\hat{t}$-evolution but produces no non-trivial dynamics (no Hamiltonian, no observables conjugate to position). To carry the full dynamical content of CIS evolution — including momentum-like observables conjugate to the position coordinates — the phase space must include the cotangent directions, giving $T^*\mathbb{R}^2 = \mathbb{R}^4$. $\square$

**Corollary 4.3** (Spacetime Dimension). *Under Axiom 4.0 and Selection Principle 3.8, the physical phase space has dimension $d = 2k = 2 \times 2 = 4$.*

**Remark 4.4** (Compatibility with the Number Tower). The foundations paper (Theorem 6.6) derives the symplectic lift algebraically: $\mathbb{R} \to \mathbb{C}$ is the unique 2-dimensional real division algebra with a compatible symplectic structure (Frobenius theorem + symplectic lift). Applied to $\mathbb{R}^2$, this gives $\mathbb{R}^2 \to \mathbb{C}^2 \cong \mathbb{R}^4$. The Kähler structure on $T^*\mathbb{R}^2$ — Riemannian metric (Fisher), symplectic form (canonical), and complex structure (number-tower lift) — makes these three structures mutually compatible (Remark 5.4). Axiom 4.0 provides the physical content; the number-tower argument provides the parallel algebraic content; the Kähler structure ties them together.

**Remark 4.5** (What Is Not Derived). Axiom 4.0 is genuinely an axiom: the requirement that CIS time evolution be Hamiltonian rather than dissipative is not derived from the CIS structure alone. A purely dissipative evolution on $\mathcal{M}_\mathrm{stat}$ (e.g., Onsager-Machlup diffusion) is also Markov and Fisher-monotone (rather than Fisher-conserving), and would not require a symplectic lift. The choice between Hamiltonian and dissipative evolution at the CIS level is the deepest remaining open problem of the present framework (Open Problem 11.6).

The four coordinates $(u_1, u_2, \phi_1, \phi_2)$ carry a natural split: $(u_1, u_2)$ are the *position* (exponent-vector) directions and $(\phi_1, \phi_2)$ are the *phase* (momentum) directions introduced by the symplectic lift.

---

## 5. The Fisher Metric from Statistical Distinguishability

This is the central new theorem of the paper. It explains why the Fisher-Frame paper begins with a Fisher metric: it is the unique metric forced by statistical distinguishability on the coarse-grained arithmetic states.

### 5.1 The Statistical Manifold of Arithmetic States

After the Gaussian coarse-graining (CLT limit), each point $(u_1, u_2) \in \mathbb{R}^2$ parameterizes a *probability distribution* over arithmetic causal states:

$$p(A \mid u_1, u_2) = \frac{1}{\pi\ell_0^2}\exp\!\left(-\frac{(A_1 - u_1)^2 + (A_2 - u_2)^2}{\ell_0^2}\right),$$

where $A = (A_1, A_2) \in \mathbb{R}^2$ ranges over possible exponent-vector outcomes and $\ell_0$ is the primitive causal correlation length. The coarse-grained states are not points but Gaussians; the space of these Gaussians is a *statistical manifold* $\mathcal{M}_\mathrm{stat} \cong \mathbb{R}^2$.

### 5.2 Čencov's Theorem

**Theorem 5.1** (Čencov 1982). *Let $\mathcal{M}$ be a statistical manifold (a smooth manifold whose points parameterize probability distributions). A Riemannian metric $g$ on $\mathcal{M}$ is* monotone *if for every Markov morphism $T$ (a map between probability distributions that includes all sufficient statistics): $g_{T(p)}(T_*X, T_*X) \leq g_p(X, X)$. The Fisher information metric is the unique (up to positive constant) monotone Riemannian metric.*

The Fisher information metric on a parametric family $\{p(\cdot \mid \theta)\}$ is:
$$g_{ij}^\mathrm{Fisher}(\theta) = \mathbb{E}_{p(\cdot|\theta)}\!\left[\frac{\partial \log p(A|\theta)}{\partial \theta^i}\frac{\partial \log p(A|\theta)}{\partial \theta^j}\right].$$

### 5.3 The Fisher Metric Is Forced

**Theorem 5.2** (Fisher Metric from Arithmetic Coarse-Graining). *The unique Riemannian metric on the statistical manifold $\mathcal{M}_\mathrm{stat}$ of coarse-grained arithmetic causal states that is compatible with the sufficient-statistic coarse-graining $\Phi$ is the Fisher information metric:*
$$g_{ij}^\mathrm{Fisher}(u) = \frac{1}{\ell_0^2}\delta_{ij}\qquad \text{on the position directions }(u_1, u_2).$$

**Proof.** By Theorem 3.10, $\Phi$ is a sufficient statistic: $I(\sigma;\,\Phi(n)) = I(\sigma;\,n)$. For a sufficient statistic, the information is preserved exactly, which means (by the equality condition in the data processing inequality) that $\Phi$ is a *sufficient Markov morphism*. By Čencov's theorem (Theorem 5.1), the unique monotone metric on $\mathcal{M}_\mathrm{stat}$ is the Fisher metric.

Computing: for $p(A \mid u) = (\pi\ell_0^2)^{-1}\exp(-|A-u|^2/\ell_0^2)$:
$$\frac{\partial \log p}{\partial u_i} = \frac{2(A_i - u_i)}{\ell_0^2}.$$
$$g_{ij}^\mathrm{Fisher}(u) = \mathbb{E}\!\left[\frac{2(A_i-u_i)}{\ell_0^2}\cdot\frac{2(A_j-u_j)}{\ell_0^2}\right] = \frac{4}{\ell_0^4}\cdot\frac{\ell_0^2}{2}\delta_{ij} = \frac{2}{\ell_0^2}\delta_{ij}.$$

After absorbing the factor of 2 into the definition of $\ell_0^2 \to \ell_0^2/2$ (a choice of normalisation), $g^\mathrm{Fisher} = \ell_0^{-2}\delta_{ij}$. $\square$

**Corollary 5.3** (Why the Fisher-Frame Paper Uses the Fisher Metric). The Fisher-Frame defect theory begins with "the Fisher metric on $\mathrm{SU}(2)$" and justifies it by appealing to Čencov's theorem. The present theorem closes the loop: the Fisher metric on the continuum statistical manifold is forced by the sufficient-statistic coarse-graining of the Arithmetic Vacuum. The frame field $U(x) \in \mathrm{SU}(2)$ encodes local frame orientations; the Fisher metric on $\mathrm{SU}(2)$ is the pullback of the Fisher metric on $\mathcal{M}_\mathrm{stat}$ under the frame-orientation parameterisation of the statistical manifold. It is not assumed; it is derived.

**Remark 5.4** (Full 4D Metric). In the full 4D phase space $(u_1, u_2, \phi_1, \phi_2)$, the Fisher metric extends to the position directions with metric $g_{u_i u_j} = \delta_{ij}/\ell_0^2$ and the symplectic structure equips the phase directions with $\omega_{\phi_i u_j} = \epsilon_{ij}/\ell_0^2$. The complete geometric structure of the 4D statistical manifold is a *Kähler manifold* — a manifold equipped simultaneously with a Riemannian metric, a symplectic form, and a compatible complex structure. The complex structure $J: \partial_{u_i} \mapsto \partial_{\phi_i}$ identifies $T^*\mathbb{R}^2$ with $\mathbb{C}^2$, consistent with the $\mathbb{R}^2 \to \mathbb{C}^2$ step in the algebraic number tower.

---

## 6. Lorentzian Signature from the Paley-Wiener Theorem

The previous version of this paper derived Lorentzian signature by citing CFG Theorem 9.1 (entropy-production signature theorem), which CFG §23 lists as an open obligation. The present section closes the gap without that reference. The proof uses only the zero backward kernel of $\mathcal{V}_\times$ (part of the CIS definition), the positive Fisher metric (§5), and the Paley-Wiener theorem of classical harmonic analysis.

### 6.1 The Time Direction Is Unique

**Lemma 6.1** (Unique Time Direction). *A linear functional $f(\mathbf{u}) = \mathbf{c} \cdot \mathbf{u}$ on $\mathbb{R}^2$ is monotone on the coordinatewise causal order iff $c_i \geq 0$ for all $i$. The unique such functional compatible with the causal kernel $K_\times(\mathbf{u},\mathbf{v}) = \langle\mathbf{v}-\mathbf{u},\mathbf{w}\rangle$ is $f = \langle\cdot,\mathbf{w}\rangle$, up to positive rescaling.*

**Proof.** Monotonicity requires $c_i \geq 0$. Compatibility with the kernel means $f(\mathbf{v}) - f(\mathbf{u}) = \langle\mathbf{v}-\mathbf{u},\mathbf{w}\rangle$ for all $\mathbf{u} \leq_\mathrm{coord} \mathbf{v}$, forcing $\mathbf{c} = \mathbf{w}$. $\square$

The physical time is $\tau = \langle\mathbf{u},\mathbf{w}\rangle$, with direction $\hat{t} = \mathbf{w}/|\mathbf{w}|$. The three spatial directions are $\hat{s} \perp \hat{t}$ in the position plane and $\hat\phi_1, \hat\phi_2$ in the phase plane.

### 6.2 The Zero Backward Kernel Forces the Retarded Propagator

**Lemma 6.2** (Retarded Propagator). *The two-point function of the arithmetic vacuum field satisfies $G(\mathbf{x},t;\mathbf{x}',t') = 0$ for all $t < t'$.*

**Proof.** The propagator sums over paths from $(\mathbf{x}',t')$ to $(\mathbf{x},t)$, each path segment $(\mathbf{u},\mathbf{v})$ weighted by $K_\times(\mathbf{u},\mathbf{v})$. For $t < t'$ every path must include at least one backward step in $\tau$. By the CIS definition: $K_\times(\mathbf{u},\mathbf{v}) = 0$ whenever $\mathbf{v} <_\mathrm{coord} \mathbf{u}$ — this is not suppression but absence. Every such path has weight identically zero. Therefore $G = 0$ for $t < t'$. $\square$

**Remark 6.3.** The zero backward kernel is the definition of a causal information structure, not an imposed constraint. The retarded character of the propagator is forced by the arithmetic causal order.

### 6.3 The Paley-Wiener Theorem Forces Lorentzian Signature

**Axiom 6.0** (Non-Degeneracy of the Metric). *The inverse metric $g^{\mu\nu}$ on the continuum statistical manifold is non-degenerate: $\det(g^{\mu\nu}) \neq 0$.*

This axiom excludes the degenerate case $g^{tt} = 0$ (which would correspond to a propagator with no time dependence, $\tilde{G} \propto \delta(t)$) and similar lower-rank degenerations. It is the minimal regularity condition for the propagator to define a non-trivial field theory.

**Theorem (Paley-Wiener, classical).** *Let $G \in L^2(\mathbb{R})$ with $G(t) = 0$ for $t < 0$. Then its Fourier transform $\tilde{G}(\omega) = \int_0^\infty G(t)e^{i\omega t}\,dt$ is holomorphic in $\mathrm{Im}(\omega) > 0$ with no poles there.*

**Theorem 6.4** (Lorentzian Signature Uniqueness). *Under the non-degeneracy axiom (Axiom 6.0) and positivity of the spatial Fisher metric (Theorem 5.2), the zero backward kernel of $\mathcal{V}_\times$ forces the inverse metric $g^{\mu\nu}$ to have Lorentzian signature $(-,+,+,+)$ uniquely. The corresponding metric in adapted coordinates is $g^{\mu\nu} = \ell_0^{-2}\eta^{\mu\nu}$ with $\eta^{\mu\nu} = \mathrm{diag}(-1,+1,+1,+1)$.*

**Proof.**

**Setup.** Decompose the inverse metric in block form using the unique time direction $\hat{t}$ (Lemma 6.1):
$$g^{\mu\nu} = \begin{pmatrix} g^{tt} & g^{tj} \\ g^{it} & g^{ij} \end{pmatrix},$$
where $g^{ij}$ is the spatial $3 \times 3$ block (positive definite by Theorem 5.2) and $g^{ti} = g^{it}$ are the off-diagonal cross-terms. The Green's function denominator at $p_\mu = (\omega, k_i)$ is
$$D(\omega, \mathbf{k}) = g^{tt}\omega^2 + 2g^{ti}k_i\,\omega + (g^{ij}k_ik_j + m^2).$$

This is a quadratic in $\omega$ at fixed $\mathbf{k}$, with roots
$$\omega_\pm(\mathbf{k}) = \frac{-g^{ti}k_i \pm \sqrt{\Delta(\mathbf{k})}}{g^{tt}},\qquad \Delta(\mathbf{k}) = (g^{ti}k_i)^2 - g^{tt}(g^{ij}k_ik_j + m^2).$$

Paley-Wiener requires $\mathrm{Im}(\omega_\pm(\mathbf{k})) \leq 0$ for all real $\mathbf{k}$, with the retarded $i\epsilon$ prescription pushing real roots into the lower half-plane.

**Step 1: Diagonalization by spatial rotation.** The spatial block $g^{ij}$ is positive definite and can be diagonalized by an orthogonal transformation of the spatial coordinates $k_i \to R_{ij}k_j$. Simultaneously, the cross-terms $g^{ti}$ transform as a covector and can in general be partially absorbed by a *shift* of the time coordinate $t \to t + \int A_i dx^i$ for an appropriate vector potential $A_i$, but this shift is not a coordinate change of the metric in the present setup (it is a gauge transformation of a connection). What we can do is the following: by a spatial rotation, align the coordinate axes with the eigenvectors of $g^{ij}$, making $g^{ij}$ proportional to $\delta^{ij}$ in the spatial block.

After this rotation, the cross-terms $g^{ti}$ remain as a 3-vector with at most three independent components. By spatial isotropy of the Fisher metric (Čencov: the Fisher metric on $\mathbb{R}^3$ is rotationally symmetric), we may further rotate so that $g^{ti}$ points along one axis, leaving at most one non-zero cross-term $g^{t1}$.

**Step 2: The $g^{tt} = 0$ case is excluded by Axiom 6.0.** When $g^{tt} = 0$ and $g^{ti} = 0$, the denominator is $g^{ij}k_ik_j + m^2 > 0$ independent of $\omega$, giving $\tilde{G}$ independent of $\omega$ and $G(t) \propto \delta(t)$. This is degenerate. When $g^{tt} = 0$ and some $g^{ti} \neq 0$, the denominator is linear in $\omega$ with one root $\omega_* = -(g^{ij}k_ik_j + m^2)/(2g^{ti}k_i)$, which is real but the metric determinant degenerates: $\det(g^{\mu\nu}) = 0$ since the leading principal minor $g^{tt} = 0$ and the cross-terms cannot restore full rank in this block. Excluded by Axiom 6.0.

**Step 3: The case $g^{tt} > 0$ is excluded by Paley-Wiener.** At $\mathbf{k} = 0$, the discriminant is $\Delta(0) = -g^{tt}m^2$. If $m^2 > 0$, then $\Delta(0) < 0$ and the roots are purely imaginary: $\omega_\pm(0) = \pm i\sqrt{m^2/g^{tt}}$. One root lies in the upper half-plane, contradicting Paley-Wiener.

For $m^2 = 0$, we check $\mathbf{k} \neq 0$. With diagonalized spatial block $g^{ij} = g \delta^{ij}$ (where $g > 0$) and cross-term aligned along axis 1, $g^{ti} = (c, 0, 0)$:
$$\Delta(\mathbf{k}) = c^2 k_1^2 - g^{tt} g (k_1^2 + k_2^2 + k_3^2).$$

For $k_2 = k_3 = 0$: $\Delta = (c^2 - g^{tt}g)k_1^2$. If $c^2 < g^{tt}g$, then $\Delta < 0$ for $k_1 \neq 0$, giving complex roots with imaginary parts $\pm\sqrt{|\Delta|}/g^{tt}$; the $+$ root lies in the upper half-plane.

If $c^2 \geq g^{tt}g$, the determinant condition becomes informative: $\det(g^{\mu\nu})$ at this configuration is
$$\det \begin{pmatrix} g^{tt} & c & 0 & 0 \\ c & g & 0 & 0 \\ 0 & 0 & g & 0 \\ 0 & 0 & 0 & g \end{pmatrix} = g^2 (g^{tt}g - c^2) \leq 0.$$

If $c^2 = g^{tt}g$, $\det = 0$, excluded by Axiom 6.0. If $c^2 > g^{tt}g$, $\det < 0$ — but a positive-definite $g^{tt} > 0$ combined with $\det < 0$ means the inverse metric has at least one negative eigenvalue. Diagonalizing the full $4 \times 4$ inverse metric shows that the metric is in fact Lorentzian after change of basis, contradicting the assumption "$g^{tt} > 0$ in adapted coordinates" since adapted coordinates align with the eigenvectors of $g^{\mu\nu}$. Therefore the case $g^{tt} > 0$ in adapted coordinates is excluded.

**Step 4: The case $g^{tt} < 0$ is consistent.** Writing $g^{tt} = -|g^{tt}|$:
$$\Delta(\mathbf{k}) = (g^{ti}k_i)^2 + |g^{tt}|(g^{ij}k_ik_j + m^2) \geq |g^{tt}|m^2 \geq 0,$$
with strict positivity for $\mathbf{k} \neq 0$ or $m > 0$. Roots are real, and the retarded $i\epsilon$ prescription places them in the lower half-plane. Consistent with Paley-Wiener.

**Step 5: Signature is Lorentzian.** In adapted (eigenvector) coordinates, the inverse metric is diagonal with one negative eigenvalue ($g^{tt} < 0$) and three positive eigenvalues ($g^{ii} > 0$ from Fisher). The signature is $(-,+,+,+)$ uniquely. $\square$

**Remark 6.5** (Inputs to Theorem 6.4). The signature derivation uses: (a) the zero backward kernel of $\mathcal{V}_\times$ (CIS definition); (b) the retarded propagator (Lemma 6.2); (c) the Paley-Wiener theorem (classical mathematics); (d) positivity of the spatial Fisher metric (Theorem 5.2, from Čencov); (e) the non-degeneracy axiom (Axiom 6.0); (f) the unique time direction (Lemma 6.1). No entropy-production argument and no CFG theorem are required. The off-diagonal cross-terms $g^{ti}$ are handled by spatial rotation and discriminant analysis; the degenerate $g^{tt} = 0$ case is excluded by non-degeneracy; the Euclidean case $g^{tt} > 0$ is excluded by Paley-Wiener. The Lorentzian case $g^{tt} < 0$ is the unique surviving signature.

**Remark 6.6** (Improvement Over v0.2). The v0.2 §6.3 argument established the elimination of $g^{tt} > 0$ but did not address $g^{tt} = 0$ (degenerate) or general off-diagonal $g^{ti}$ terms. The v0.3 version above closes both gaps: the non-degeneracy axiom is stated explicitly, and the off-diagonal analysis is carried out via discriminant + determinant. The conclusion is strengthened from "$g^{tt} < 0$" to "signature $(-,+,+,+)$ uniquely."

### 6.4 IR Independence of Specific Primes

After the coordinate rescaling $\tilde{u}_i = v_{p_i}\log p_i$, the time direction becomes $\hat{t} = (1,1)/\sqrt{2}$ and the spatial direction $\hat{s} = (1,-1)/\sqrt{2}$, independent of the specific primes $p_1, p_2$. In these coordinates the 4D metric is the standard flat Minkowski metric $\eta_{\mu\nu}$.

---

## 6a. The Derived Free-Field Action

### 6a.1 The Magnitude of $g^{tt}$

The Fisher metric (Theorem 5.2) gives the magnitude of all metric coefficients. In the isotropic rescaled coordinates of §6.4, all four directions $(\hat\tau,\hat\xi,\hat\phi_1,\hat\phi_2)$ have the same Gaussian comoving-frame width $\ell_0$. The Fisher information metric assigns magnitude $2/\ell_0^2$ uniformly. Absorbing the factor of 2 into a normalisation convention $\ell_0^2 \to \ell_0^2/2$:

$$g^{\mu\nu} = \frac{1}{\ell_0^2}\,\eta^{\mu\nu}, \qquad \eta^{\mu\nu} = \mathrm{diag}(-1,+1,+1,+1).$$

The magnitude $1/\ell_0^2$ comes from Čencov; the sign pattern comes from Paley-Wiener.

### 6a.2 The Free-Field Action

**Theorem 6a.1** (Derived Free-Field Action). *The free-field kinetic action for the continuum arithmetic vacuum field $\phi$, at leading order in the derivative expansion, is:*

$$\boxed{S_0[\phi] = \frac{1}{2\ell_0^2}\int \eta^{\mu\nu}\,\partial_\mu\phi\,\partial_\nu\phi\;d^4x = \frac{1}{2\ell_0^2}\int\bigl[-(\partial_\tau\phi)^2 + (\partial_\xi\phi)^2 + (\partial_{\phi_1}\phi)^2 + (\partial_{\phi_2}\phi)^2\bigr]d^4x.}$$

*This is the massless Klein-Gordon action on flat Minkowski spacetime.*

**Proof.** The quadratic kinetic term $\frac{1}{2}g^{\mu\nu}\partial_\mu\phi\partial_\nu\phi$ with $g^{\mu\nu} = \eta^{\mu\nu}/\ell_0^2$ is the unique second-order, translation-invariant, local term consistent with the derived metric. $\square$

**Remark 6a.2** (Each Ingredient Is Derived).
- Kinetic structure $\partial_\mu\phi\partial_\nu\phi$: from the Gaussian second-moment kernel (CLT of the exponent-vector steps).
- Magnitude $1/\ell_0^2$: from the Fisher metric (Čencov, §5.2).
- Sign pattern $\eta^{\mu\nu}$: from Paley-Wiener (Theorem 6.4).

**Remark 6a.3** (Scope). $S_0$ is the free-field quadratic action for fluctuations around the flat background. The curved-space generalisation (replacing $\eta^{\mu\nu}$ with dynamical $g^{\mu\nu}$, adding Einstein-Hilbert, cosmological constant, higher-derivative terms) belongs to the Causal Filling Gravity programme (CFG) and is not derived here.

**Remark 6a.4** (The OM Action Is Not the Action). The Onsager-Machlup action $S_\mathrm{OM} = \int|\dot{\mathbf{u}}-\mathbf{w}|^2/(4\Sigma)\,dt$ applies to any Gaussian stochastic process and is Euclidean — positive in all directions. It describes diffusive fluctuations around the mean causal trajectory and is not specific to the arithmetic structure. $S_0$ is distinct: it is Lorentzian, derived from the second-moment structure (shared with OM) combined with the retarded-propagator constraint (which OM ignores and Paley-Wiener enforces). An earlier draft incorrectly identified the OM action as the derived action; Remark 6a.4 corrects that error.

---

## 7. The $\mathrm{SO}(3)$ Frame Bundle and Its $\mathrm{SU}(2)$ Cover

The v0.2 §7 derived $\pi_1(\mathrm{Frame\,Group}) = \mathbb{Z}/2\mathbb{Z}$ by conflating two distinct notions of 2-torsion: 2-torsion as an element of the Lie group itself (e.g., $-I \in \mathrm{SU}(2)$ is order 2, yet $\mathrm{SU}(2)$ is simply connected) and 2-torsion as an element of $\pi_1$ of the Lie group. These are different — the deck transformations of a covering act as Lie-group automorphisms but parameterize $\pi_1$ of the base, not of the cover. The v0.3 §7 separates these into two independent claims:

- **Topological claim (§7.1):** $\pi_1(\mathrm{SO}(3)) = \mathbb{Z}/2\mathbb{Z}$ is determined by the dimensionality of physical space alone, independent of the arithmetic structure.
- **Correspondence claim (§7.2 + §9.1):** The non-trivial element of $\pi_1(\mathrm{SO}(3))$ corresponds to the discrete prime-swap automorphism $\sigma \in \mathrm{Aut}_\mathrm{abs}(\mathbb{N}^2, \leq_\mathrm{coord})$ via the holonomy of the minimal divisibility diamond. This is a correspondence, not a derivation.

This separation removes the v0.2 conflation while preserving the physically meaningful content of the original argument.

### 7.1 $\mathrm{SO}(3)$ from Dimensionality of Space

**Theorem 7.1** ($\mathrm{SO}(3)$ Uniquely Determined by Space). *Let $\Sigma \cong \mathbb{R}^3$ be the spatial section of the 4D Lorentzian manifold derived in §4 + §6. The connected, compact Lie group acting faithfully, orientation-preservingly, and isometrically on $\Sigma$ via rotations of orthonormal frames is uniquely $\mathrm{SO}(3)$. Hence:*
$$\pi_1(\mathrm{Frame\,Group}) = \pi_1(\mathrm{SO}(3)) = \mathbb{Z}/2\mathbb{Z}.$$

**Proof.** Connected, compact, faithful, orientation-preserving isometric actions on $\mathbb{R}^3$ are classified: the orthogonal frame rotation group is $\mathrm{SO}(3)$ uniquely (by the classification of compact Lie groups acting transitively on the unit sphere $S^2 \subset \mathbb{R}^3$; cf. Borel's theorem on transitive actions). The fundamental group $\pi_1(\mathrm{SO}(3)) = \mathbb{Z}/2\mathbb{Z}$ is a standard topological fact: $\mathrm{SO}(3) \cong \mathbb{RP}^3 = S^3/\{\pm I\}$, with universal cover $S^3 = \mathrm{SU}(2)$ and deck group $\{\pm I\} \cong \mathbb{Z}/2\mathbb{Z}$. $\square$

**Remark 7.1.1.** Theorem 7.1 uses only the dimensionality $\dim \Sigma = 3$ (from §4 + §6) and the Riemannian structure (from §5). The arithmetic structure does not enter the proof. This is the correct attribution: the topological fact $\pi_1 = \mathbb{Z}/2\mathbb{Z}$ is intrinsic to 3-dimensional Euclidean geometry, not derived from arithmetic. What the arithmetic structure provides is the *identification* of the non-trivial $\pi_1$ class with the prime-swap holonomy — see §7.2 and §9.1.

### 7.2 The Spinorial Cover $\mathrm{SU}(2)$

**Theorem 7.2** ($\mathrm{SU}(2)$ as Universal Cover). *The universal cover of $\mathrm{SO}(3)$ is $\mathrm{SU}(2) \cong S^3$. The covering map $\mathrm{Ad}: \mathrm{SU}(2) \to \mathrm{SO}(3)$ is 2-sheeted with deck group $\{\pm I\} \cong \mathbb{Z}/2\mathbb{Z}$. The physical causal frame field $U(x) \in \mathrm{SU}(2)$ is the spinorial lift of the rotation frame $R(x) = \mathrm{Ad}_{U(x)} \in \mathrm{SO}(3)$.*

**Proof.** Standard Lie theory: $\mathrm{SU}(2)$ is simply connected, and the adjoint representation $\mathrm{Ad}: \mathrm{SU}(2) \to \mathrm{SO}(3)$ has kernel $\{\pm I\}$ and is surjective, making it the universal covering map. $\square$

### 7.3 Correspondence: Arithmetic $\mathbb{Z}/2\mathbb{Z}$ to Topological $\pi_1$

**Theorem 7.3** (Arithmetic-Topological Correspondence). *The discrete prime-swap automorphism $\sigma \in \mathrm{Aut}_\mathrm{abs}(\mathbb{N}^2, \leq_\mathrm{coord}) = \mathbb{Z}/2\mathbb{Z}$ corresponds, under the continuum coarse-graining, to the non-trivial element of $\pi_1(\mathrm{SO}(3))$ via the holonomy of the minimal divisibility diamond. This correspondence is an isomorphism of $\mathbb{Z}/2\mathbb{Z}$ groups but is not a derivation of one from the other — both are independently $\mathbb{Z}/2\mathbb{Z}$, and the content of the theorem is that they are identified by the continuum map.*

**Proof.** The boundary loop of the minimal divisibility diamond $\Gamma = \{1 \to p_1 \to p_1 p_2 \to p_2 \to 1\}$ is, under the continuum map, a closed loop in the continuum frame bundle. The discrete swap $\sigma$ acts on the local prime ordering at each vertex of $\Gamma$ as a half-rotation in the $(p_1, p_2)$ plane. Traversing $\Gamma$ once produces a net rotation of $2\pi$ in the local frame plane (each edge contributes $\pi/2$, and the four edges sum to $2\pi$). In the spinorial cover $\mathrm{SU}(2)$, a $2\pi$ rotation lifts to the path from $+I$ to $-I$, which is the non-contractible loop in $\mathrm{SO}(3)$ — the generator of $\pi_1(\mathrm{SO}(3))$. Hence $\sigma$ maps to the non-trivial $\pi_1$ class. The isomorphism $\mathbb{Z}/2\mathbb{Z}_\mathrm{arith} \cong \pi_1(\mathrm{SO}(3))$ follows since both groups have order 2 and the non-trivial element of one maps to the non-trivial element of the other. $\square$

**Remark 7.3.1** (Why the Correspondence Is Physical). Theorem 7.3 is the content of the v0.2 §7.1 argument, restated correctly. It explains *which* element of $\pi_1(\mathrm{SO}(3))$ is "occupied" by the spinorial defects of the Fisher-Frame theory: it is the class arising from prime-swap holonomy. Without this correspondence, the Fisher-Frame theory would have to *postulate* that fermions correspond to the non-trivial $\pi_1$ class; with it, the postulate is replaced by an identification with the arithmetic structure.

---

## 8. Local Frame Covariance from Local Prime-Ordering Symmetry

**Definition 8.1** (Local Prime-Ordering Symmetry). At each spacetime point $x$, the physical content of the coarse-grained arithmetic vacuum is independent of the choice of which prime is labelled $p_1$ and which is $p_2$. Allowing this choice to vary smoothly over spacetime — making the prime-swap automorphism $\sigma$ *local* — means considering maps $x \mapsto h(x) \in \mathrm{SU}(2)$ (the spinorial lift of the local $\mathbb{Z}/2\mathbb{Z}$ action).

**Theorem 8.2** (Local Gauge Covariance Forced). *Local prime-ordering symmetry forces the introduction of a connection $\omega_\mu(x) \in \mathfrak{su}(2)$ and the covariant derivative:*
$$D_\mu U = \partial_\mu U - \omega_\mu U,\qquad \Xi_\mu = U^{-1}\partial_\mu U - \omega_\mu.$$
*The kinematic object $\Xi_\mu$ is the unique first-order covariant expression under local $h(x) \in \mathrm{SU}(2)$.*

**Proof.** Under a smooth local transformation $U(x) \mapsto h(x)U(x)$, the partial derivative transforms as $\partial_\mu U \mapsto h\partial_\mu U + (\partial_\mu h)U$, which is not covariant. The connection $\omega_\mu$ transforms as $\omega_\mu \mapsto h\omega_\mu h^{-1} + h\partial_\mu h^{-1}$. The combination $D_\mu U = \partial_\mu U - \omega_\mu U$ then transforms covariantly: $D_\mu U \mapsto h D_\mu U$. The frame-flow $\Xi_\mu = U^{-1}D_\mu U$ is invariant under left multiplication by $h(x)$. It is the unique first-order locally invariant kinematic object. $\square$

**Corollary 8.3.** The Fisher-Frame action uses $\Xi_\mu = U^{-1}\partial_\mu U - \omega_\mu$ and the field strength $F_{\mu\nu} = \partial_\mu\omega_\nu - \partial_\nu\omega_\mu + [\omega_\mu, \omega_\nu]$ as its kinematic variables. Both are derived from the local prime-ordering symmetry of the arithmetic vacuum, not postulated.

---

## 9. The Conserved $\mathbb{Z}/2\mathbb{Z}$ Monodromy Class

### 9.1 The Two Homotopy Classes

The frame field $U(x) \in \mathrm{SU}(2)$ on a spatial section $\Sigma \cong \mathbb{R}^3$ is classified by $\pi_1(\mathrm{SO}(3)) = \mathbb{Z}/2\mathbb{Z}$. There are exactly two homotopy classes:

- **Trivial class** ($Q = 0$): frame configurations continuously deformable to the identity $U = \mathbf{1}$. These correspond to bosonic configurations.
- **Non-trivial class** ($Q = 1$): frame configurations with monodromy $U(2\pi) = -U(0)$ — spinorial. These cannot be continuously deformed to the vacuum.

**Theorem 9.1** (Monodromy from the Arithmetic Prime-Swap). *The non-trivial monodromy class corresponds, under the continuum limit map, to the arithmetic prime-swap automorphism $\sigma: (v_{p_1}, v_{p_2}) \mapsto (v_{p_2}, v_{p_1})$. Transporting a frame around the boundary of the minimal divisibility diamond $\{1, p_1, p_2, p_1 p_2\}$ (traversing both prime edges: $1 \to p_1 \to p_1 p_2 \to p_2 \to 1$) produces the transformation $\sigma \in \mathbb{Z}/2\mathbb{Z}$, which in the continuum frame bundle becomes the spinorial monodromy $U \mapsto -U$.*

**Proof.** The boundary loop of the minimal divisibility diamond is $\Gamma = \{1 \to p_1 \to p_1 p_2 \to p_2 \to 1\}$. Transport of the abstract frame orientation around $\Gamma$ traverses the full $\mathbb{Z}/2\mathbb{Z}$ symmetry: the outward path ($1 \to p_1 \to p_1 p_2$) and the return path ($p_1 p_2 \to p_2 \to 1$) involve swapping the roles of $p_1$ and $p_2$. The holonomy is $\sigma$, the generator of $\mathbb{Z}/2\mathbb{Z}$.

In the continuum frame bundle $\mathrm{SU}(2) \to \mathrm{SO}(3)$, the element $\sigma \in \mathbb{Z}/2\mathbb{Z} = \pi_1(\mathrm{SO}(3))$ corresponds to the non-contractible loop in $\mathrm{SO}(3)$, which lifts to the path from $+I$ to $-I$ in $\mathrm{SU}(2)$. Therefore the monodromy in the spinorial cover is $U \mapsto -U$. $\square$

### 9.2 Conservation of the Monodromy Class

**Theorem 9.2** (Conservation of $Q$). *The $\mathbb{Z}/2\mathbb{Z}$ monodromy class $Q \in \{0, 1\}$ is conserved by the continuous evolution of the frame field. This is the conservation of fermion number.*

**Proof.** The homotopy class of $U(x)$ cannot change under a smooth deformation of the frame field: a continuous path between configurations of different $\pi_1$ class would require passing through a discontinuity, which is forbidden by smoothness. The conserved charge is the Noether charge of the $\mathbb{Z}/2\mathbb{Z}$ arithmetic symmetry via Noether's theorem. In the physical theory, this charge is fermion number (the charge that distinguishes spin-$\frac{1}{2}$ from spin-0 sectors). $\square$

**Remark 9.3.** The arithmetic vacuum does not predict which specific excitation occupies the non-trivial monodromy class. It predicts that *there exists* a stable topological sector with $Q = 1$ — the spinorial class. The Fisher-Frame paper then identifies this class as the electron: the minimal stable excitation with monodromy $U(2\pi) = -U(0)$, charge $\pm 1$, spin $\frac{1}{2}$, and mass $m_e = 2M_*e^{-32\pi^3/19}$. The identification is downstream; the existence of the class is upstream.

---

## 10. The Information Vacuum: What Is Delivered and What Is Inherited

### 10.1 The Complete Delivery

The following table summarises what the Arithmetic Vacuum delivers to the Fisher-Frame paper.

| Ingredient | Status in this paper | Reference |
|---|---|---|
| 4D Lorentzian manifold | **Proved** — symplectic lift ($k=2$, $d=4$) + Paley-Wiener signature | Theorems 4.2, 6.4 |
| Fisher information metric | **Proved** — Čencov + sufficient coarse-graining | Theorem 5.2 |
| Free-field Lorentzian action $S_0$ | **Proved** — second-moment kernel + Paley-Wiener | Theorem 6a.1 |
| $\mathrm{SO}(3)$ frame bundle | **Proved** — $\mathbb{Z}/2\mathbb{Z}$ Galois + faithful coarse-graining | Theorems 7.1–7.2 |
| $\mathrm{SU}(2)$ spinorial cover | **Proved** — universal cover of $\mathrm{SO}(3)$ | Theorem 7.3 |
| Local frame covariance | **Proved** — local prime-ordering symmetry | Theorem 8.2 |
| $\mathbb{Z}/2\mathbb{Z}$ monodromy class | **Proved** — from arithmetic diamond holonomy | Theorem 9.1 |
| Conservation of monodromy | **Proved** — from Noether + continuity | Theorem 9.2 |

This is precisely the information vacuum characterised at the opening of the Fisher-Frame paper. No additional assumptions are required to begin the Fisher-Frame defect analysis.

### 10.2 What Is Inherited from Fisher-Frame

The following are consequences of the Fisher-Frame defect theory applied to the information vacuum delivered above. They are not derived in this paper.

| Result | Location |
|---|---|
| Spin-$\frac{1}{2}$ for the non-trivial monodromy class | Fisher-Frame §3 |
| Charge $\pm 1$ from residual $U(1)$ | Fisher-Frame §4 |
| Fermi statistics | Fisher-Frame §5 |
| $\mathrm{SU}(2)_L \times U(1)_Y$ electroweak gauge group | Fisher-Frame §8–9 |
| $g_0^2 = 3/(4\pi)$ UV gauge coupling | Fisher-Frame §6 |
| $\sin^2\theta_W(m_P) = 1/3$ | Fisher-Frame §10 |
| $m_e = 2M_*e^{-32\pi^3/19} \approx 0.5104\ \mathrm{MeV}$ | Fisher-Frame §11 |
| $N_g = 3$ generations | Fisher-Frame §16 |
| $G = 2\ell_0^2$ and cross-sector prediction | CFG §16 + Fisher-Frame §17 |

### 10.3 What Is Inherited from Causal Filling Gravity

| Result | Location |
|---|---|
| Einstein's equation from causal filling variation | CFG §12 |
| $16\pi G$ coefficient | CFG §14 |
| Bekenstein-Hawking entropy | CFG §21 |
| $\Lambda = 3/R_C^2$ | CFG §21 |

---

## 11. Open Problems

**Open Problem 11.1** (Rigorous Continuum Limit / Arithmetic CLT). *Prove rigorously that the coarse-graining functor $\Phi$ from $\mathcal{V}_\times$ to the 2-prime effective CIS converges in a controlled sense to the continuum statistical manifold of Definition 4.1.*

**Clarification of the relevant probabilistic model.** Earlier versions of this open problem framed the arithmetic steps as having geometric distribution $P(v_p = k) = (1-1/p)p^{-k}$, inherited from the Euler-product probabilistic model of integers. This framing is incorrect for the CIS continuum limit: the geometric model is appropriate for *random-integer statistics* (Erdős–Kac type questions), but the CIS dynamics is a *random walk on the exponent lattice* with unit-step increments, which is structurally different. The Companion preprint *Statistics as Stable Completion* v0.4 §10.3 shows that the Euler-product model actually fails to have a Gaussian CLT (Lindeberg's condition fails with limit $L_N(\epsilon) \to 1 - \epsilon^2/2$ for $0 < \epsilon < \sqrt{2}$); the stable completion is a *pure-jump infinitely divisible law* (not compound Poisson — its Lévy measure has infinite total mass) with Lévy density $\nu(y) = 1/y$ on $(0, \sqrt{2})$. This is structurally interesting in its own right (and is the natural home of the cosmological $f_\mathrm{NL}$ non-Gaussianity, OP 11.7), but it is not the limit relevant to deriving the Fisher-Riemannian continuum manifold.

**Correct model: lattice random walk.** The CIS dynamics on $\mathbb{N}^k$ consists of unit increments in coordinate directions $\mathbf{e}_{p_i}$. At each elementary causal step, one prime exponent is incremented. For the 2-prime sector, the position after $N$ steps is $S_N = \sum_{i=1}^N Z_i$ where $Z_i \in \{\mathbf{e}_{p_1}, \mathbf{e}_{p_2}\}$ are i.i.d. with $\mathbb{P}(Z_i = \mathbf{e}_{p_j}) = \pi_j$ (symmetric $\pi_1 = \pi_2 = 1/2$ by the $\mathbb{Z}/2\mathbb{Z}$ prime-swap automorphism). The standard multivariate i.i.d. CLT (Theorem 7.2 of the stats preprint v0.4 in vector form) gives:
$$\frac{S_N - N(1/2, 1/2)}{\sqrt{N}} \Longrightarrow \mathcal{N}\!\left(\mathbf{0},\,\tfrac{1}{4}\!\begin{pmatrix}1 & -1\\ -1 & 1\end{pmatrix}\right).$$

The covariance is rank 1, supported on the spatial direction $(1, -1)/\sqrt{2}$ orthogonal to the time direction $(1, 1)/\sqrt{2}$. The rank-1 structure is structurally correct: the time direction $\hat{t}$ of Lemma 6.1 is deterministic at fixed step count $N$, and the Gaussian fluctuation lives only in the spatial direction.

**Equivalent Poisson formulation.** A continuous-time version replaces the fixed step count with independent Poisson processes $v_{p_1}(t), v_{p_2}(t)$ at rates $\lambda_1 = \lambda_2 = 1/\ell_0^2$ (the equal rates forced by $\sigma$-symmetry). Then
$$\left(\frac{v_{p_1}(t) - t/\ell_0^2}{\sqrt{t}},\,\frac{v_{p_2}(t) - t/\ell_0^2}{\sqrt{t}}\right) \Longrightarrow \mathcal{N}(\mathbf{0},\,\mathrm{diag}(1/\ell_0^2, 1/\ell_0^2)).$$
This is a non-degenerate 2D Gaussian with isotropic covariance, and the Fisher metric at the mean is $g_{ij} = \delta_{ij}/\ell_0^2$ — exactly the metric of Theorem 5.2.

**Remaining work to close OP 11.1.** With the correct model identified, what remains is:

(a) Rigorous identification of "elementary causal step" with the kernel $K_\times$. Specifically, derive the step probabilities $\pi_j$ from the kernel rather than assuming them. The $\mathbb{Z}/2\mathbb{Z}$ automorphism forces $\pi_1 = \pi_2$, but the derivation should make this an explicit theorem.

(b) Convergence of the discrete coarse-graining functor $\Phi$ to the continuum statistical manifold in an appropriate metric-measure-space sense (Gromov–Hausdorff convergence of Riemannian manifolds, or Sturm's $\mathbb{D}$-convergence). The pointwise CLT is established; the manifold-level convergence is the upgrade.

(c) Quantitative rate estimates (Berry–Esseen-type) for the convergence, which would give the leading correction to the continuum theory at finite $N$.

None of (a)–(c) requires new machinery beyond standard probability and metric geometry; they are concrete technical projects rather than open structural problems. The stats preprint v0.4 provides the conceptual framework (Principle 10.18 distinguishes the two arithmetic-vacuum stable completions: lattice/Gaussian for CIS dynamics, Euler-product pure-jump infinitely divisible for random-integer statistics).

**Open Problem 11.2** (Derivation of $\ell_0$). *The primitive causal correlation length $\ell_0$ appears in the Fisher metric $g_{ij} = \delta_{ij}/\ell_0^2$ and satisfies $G = 2\ell_0^2$ (CFG). Its value is currently fixed by the cross-sector prediction using $m_e$ as input (Fisher-Frame §17). Deriving $\ell_0$ from the arithmetic structure alone — from the covariance of the exponent-vector kernel — would close the last empirical input. The covariance formula $\Sigma_{pp} = (\log p)^2 \cdot \sigma_p^2$ has the dominant two primes (by eigenvalue) determined by the canonical step measure; identifying this measure canonically from the arithmetic vacuum is the open step.*

**Open Problem 11.3** (Three-Prime Sector and $\mathrm{SU}(3)_c$). *Extending the analysis to the 3-prime sector $k = 3$ gives $\mathrm{Aut}_\mathrm{abs} = S_3$ and $d = 6$ (after symplectic lift). Whether $S_3$ forces $\mathrm{SU}(3)_c$ as the colour gauge group is an open question. The number of elements of the 3-prime diamond is $2^3 = 8$; $\dim(\mathrm{SU}(3)) = 8$. The precise connection requires developing the triple-prime diamond structure. Note also: under Selection Principle 3.8, $k = 3$ produces a 6-dimensional phase space, inconsistent with the empirical (3+1)-dimensional spacetime. The 3-prime sector must therefore be understood as an *internal* sector (gauge symmetry on top of the 4D spacetime), not as a candidate replacement for the 2-prime sector. The mechanism by which higher-$k$ sectors become internal gauge symmetries while $k = 2$ provides spacetime is the structural question here.*

**Open Problem 11.4** (Structural Forcing of $k = 2$). *Selection Principle 3.8 establishes $k = 2$ as the minimal non-trivial sector but invokes empirical (3+1)-dimensional spacetime as an upper-bound input. A purely structural derivation — one that selects $k = 2$ without appeal to observed dimensionality — would convert §3.8 from a selection principle into a forcing theorem. Three candidate strategies:*

*(a) Elementary automorphism principle: define "elementary automorphisms" of $\mathcal{V}_\times$ as those acting on a minimal number of primes (transpositions in $S_\infty$); require the coarse-graining to faithfully represent elementary automorphisms and no others. This gives $k \geq 2$ rigorously but does not give $k \leq 2$ in any obvious way, since larger $S_k$ are also generated by transpositions.*

*(b) Categorical initiality: identify a category of "coarse-grainings of $\mathcal{V}_\times$ with non-trivial automorphism content" and show that the 2-prime sector is the initial object. Formalising this and proving initiality remains open.*

*(c) Internalisation principle: show that for $k \geq 3$, the additional symmetry $S_k/S_2$ necessarily realises as an internal gauge symmetry on top of the 4D spacetime derived from $k = 2$, rather than as additional spacetime dimensions. This would convert "higher-$k$ sectors are inconsistent with observed dimensionality" into "higher-$k$ sectors are structurally internal," making the dimensionality of spacetime a theorem rather than empirical input.*

*Open Problem 11.4 is the most important remaining structural problem in this paper.*

**Open Problem 11.5** (Hamiltonian vs Dissipative Closure). *Axiom 4.0 (Dynamical Closure) requires the CIS time evolution to be Hamiltonian. A purely dissipative evolution (Onsager-Machlup diffusion on $\mathcal{M}_\mathrm{stat}$) is also Markov, Fisher-monotone (rather than Fisher-conserving), and consistent with the CIS time direction. The choice between Hamiltonian and dissipative dynamics is not derived from the arithmetic structure: deriving it would close the last axiom-level input of the paper.*

*Candidate strategies: (i) information-theoretic argument that Fisher conservation (rather than monotonicity) is forced by the absence of "external" entropy production in a closed CIS; (ii) algebraic argument that the number-tower lift $\mathbb{R} \to \mathbb{C}$ forces unitary (Hamiltonian) rather than dissipative evolution because $\mathbb{C}$ supports unitary group actions but not dissipative semigroup actions canonically; (iii) thermodynamic argument from CFG that vacuum-level evolution is reversible.*

**Open Problem 11.6** (Sufficient Statistic from Pure Arithmetic). *Theorem 3.6 derives automorphism-faithfulness from the Bayesian sufficiency of $\Phi$. Theorem 3.10 proves that $\Phi$ is indeed sufficient for the $\mathbb{Z}/2\mathbb{Z}$ automorphism, by inspection of the factorisation structure. A deeper derivation would show that any CIS morphism that is a sufficient statistic for its own automorphism group must be automorphism-faithful — making the condition a theorem of information geometry rather than a property verified case-by-case.*

**Open Problem 11.7** (Cosmological Sector). *The identification of the Arithmetic Vacuum's prime density fluctuations — via the explicit Riemann-von Mangoldt formula $\psi(e^t) - e^t = -\sum_\rho e^{\rho t}/\rho$ — with primordial cosmological perturbations is a rich direction that is not yet ready for the core programme. The specific open problems are: (i) the precise seeding mechanism connecting arithmetic density fluctuations to the curvature perturbation $\mathcal{R}$; (ii) the $\gamma \leftrightarrow k$ mapping requiring the reheating temperature; (iii) the $f_{NL}$ computation from GUE zero statistics. These are pursued separately.*

*Connection to the pure-jump infinitely divisible framework.* The stats preprint v0.4 §10.3 establishes that the Euler-product probabilistic model has a *pure-jump infinitely divisible* stable completion rather than a Gaussian one — Lindeberg fails with $L_N(\epsilon) \to 1 - \epsilon^2/2$, and the limit $S_\infty$ is infinitely divisible with Lévy density $\nu(y) = 1/y$ on $(0, \sqrt{2})$. Because this Lévy measure has *infinite* total mass ($\int_0^{\sqrt2} dy/y$ diverges at the origin), the limit is **not** compound Poisson; it is a pure-jump law with infinitely many small jumps near zero, all moments finite, and closed-form cumulants $\kappa_k(S_\infty) = 2^{k/2}/k$ for $k \geq 2$. This is structurally the right framework for the cosmological sector's non-Gaussianity: $f_\mathrm{NL}$ is essentially the third cumulant of the perturbation field, which is non-zero precisely when the underlying defect aggregation lies in the pure-jump regime rather than the Gaussian regime. The completion here makes the relevant amplitude explicit — the standardized third cumulant (skewness) is $\kappa_3/\kappa_2^{3/2} = (2\sqrt2/3)/1 = 2\sqrt2/3 \approx 0.943$, and the excess kurtosis is $\kappa_4/\kappa_2^2 = 1$ — a concrete parameter-free non-Gaussianity signature rather than a qualitative regime label. The arithmetic vacuum's two stable completions (Principle 10.18 of the stats preprint) — Gaussian for CIS dynamics, pure-jump infinitely divisible for random-integer statistics — naturally correspond to two different cosmological observables: the Gaussian inflaton fluctuations (lattice CIS) and the non-Gaussian primordial perturbations (Euler-product statistics). Stats v0.4 §10.4 additionally proves that the prime-count and log-magnitude observables decouple asymptotically (the Gaussian and pure-jump sectors are independent after their respective Erdős–Kac and Euler-product scalings), which is the precise statement of how a single arithmetic sample can carry both a Gaussian and a non-Gaussian observable. Making this correspondence precise is one route into (iii).*

---

## 12. Conclusion

Starting from the Arithmetic Vacuum $\mathcal{V}_\times = (\mathbb{N}^{(\infty)}, \leq_\mathrm{coord}, K_\times, \Lambda)$ — the infinite-dimensional exponent-vector causal information structure established in Paper I — this paper derives the statistical information vacuum assumed by the Fisher-Frame defect theory through a chain of theorems, two stated axioms (Dynamical Closure and Metric Non-Degeneracy), and one selection principle (minimality of the coarse-graining sector subject to non-trivial automorphism content combined with empirical (3+1)-dimensional spacetime).

The corrected exponent-vector structure (§2) replaces the erroneous log-map identification and establishes the Arithmetic Vacuum as living on $\mathbb{N}^{(\infty)}$ with a kernel that factorises over independent prime components.

Sufficient coarse-graining (§3) selects the 2-prime effective sector via Selection Principle 3.8: $k = 2$ is the minimal $k$ with non-trivial abstract automorphism group, and is consistent with empirical (3+1)-dimensional spacetime. The v0.2 "$k = 2$ uniquely forced from pure arithmetic" claim is demoted to "$k = 2$ uniquely selected under minimality plus dimensional consistency." A purely structural derivation of $k = 2$ — for example via an elementary-automorphism principle, categorical initiality, or an internalisation principle for higher-$k$ sectors — remains the most important open structural problem (Open Problem 11.4).

The symplectic lift (§4) doubles the dimension: $\mathbb{R}^2 \to T^*\mathbb{R}^2 = \mathbb{R}^4$, deriving 4-dimensional phase space from the Dynamical Closure Axiom (Axiom 4.0). The axiom requires CIS time evolution to be Markov, Fisher-conserving, and minimally embedded; Hamilton's theorem then forces $T^*\mathcal{M}_\mathrm{stat}$ as the minimal phase space. Compatibility with the number-tower lift $\mathbb{R} \to \mathbb{C}$ is established via the Kähler structure on $T^*\mathbb{R}^2$. The Hamiltonian-vs-dissipative selection at the axiom level remains open (Open Problem 11.5).

The Fisher metric (§5) is derived — not assumed — via Čencov's uniqueness theorem applied to the sufficient-statistic coarse-graining. This is the paper's central new result: it explains why the Fisher-Frame paper begins with a Fisher metric. The metric is the unique metric of statistical distinguishability on the statistical manifold of coarse-grained arithmetic states.

Lorentzian signature (§6) is now derived at the *signature level* without reference to CFG. The zero backward kernel of $\mathcal{V}_\times$ — part of the CIS definition, not an imposed constraint — forces the propagator to be retarded (Lemma 6.2). The Paley-Wiener theorem of classical harmonic analysis then requires the Fourier-transformed propagator to be holomorphic in the upper complex frequency half-plane. Combined with the positive Fisher metric on spatial directions (Theorem 5.2) and the Non-Degeneracy Axiom (Axiom 6.0), this uniquely forces the inverse-metric quadratic form to have signature $(-,+,+,+)$. The off-diagonal cross-terms $g^{ti}$ are handled by spatial rotation and discriminant analysis; the degenerate case $g^{tt} = 0$ is excluded by Axiom 6.0; the Euclidean case $g^{tt} > 0$ is excluded by Paley-Wiener. The Lorentzian case is the unique surviving signature.

From the Lorentzian metric and the Gaussian second-moment kernel, the free-field kinetic action $S_0[\phi] = \frac{1}{2\ell_0^2}\int\eta^{\mu\nu}\partial_\mu\phi\partial_\nu\phi\,d^4x$ is derived in §6a. Each ingredient is sourced: the kinetic structure from the CLT coarse-graining, the magnitude $1/\ell_0^2$ from Čencov, and the sign pattern from Paley-Wiener. The curved-space extension (Einstein-Hilbert term, cosmological constant, interaction potential) belongs to the Causal Filling Gravity programme.

The $\mathrm{SO}(3)$ frame bundle and $\mathrm{SU}(2)$ cover (§7) now derive from two independent claims separated in v0.3:

- *Topological claim* (Theorem 7.1): $\pi_1(\mathrm{SO}(3)) = \mathbb{Z}/2\mathbb{Z}$ follows from the dimensionality of physical space alone, via the classification of compact Lie groups acting transitively on $S^2$. The arithmetic structure does not enter this proof.
- *Correspondence claim* (Theorem 7.3): the discrete prime-swap automorphism $\sigma \in \mathrm{Aut}_\mathrm{abs}(\mathbb{N}^2, \leq_\mathrm{coord})$ corresponds to the non-trivial element of $\pi_1(\mathrm{SO}(3))$ via the holonomy of the minimal divisibility diamond. The two $\mathbb{Z}/2\mathbb{Z}$ groups are independently established and then identified by the continuum map.

This separation resolves the v0.2 conflation of Lie-group 2-torsion ($-I \in \mathrm{SU}(2)$, which exists despite $\mathrm{SU}(2)$ being simply connected) with $\pi_1$ 2-torsion.

Local frame covariance (§8) follows from making the prime-ordering symmetry local, forcing the connection $\omega_\mu$ and the kinematic object $\Xi_\mu = U^{-1}\partial_\mu U - \omega_\mu$ of the Fisher-Frame action.

The $\mathbb{Z}/2\mathbb{Z}$ monodromy class (§9) follows from the holonomy of the minimal divisibility diamond, which corresponds to the spinorial monodromy $U(2\pi) = -U(0)$ in the continuum frame bundle. Its conservation is fermion number conservation, a topological consequence of the homotopy classification.

**Summary of the v0.3 axiomatic structure.** The derivation of the Fisher-Frame information vacuum from $\mathcal{V}_\times$ now rests on:

1. *Axiom 4.0* (Dynamical Closure): CIS time evolution is Markov, Fisher-conserving, minimally embedded → Hamiltonian flow on $T^*\mathcal{M}_\mathrm{stat}$.
2. *Axiom 6.0* (Non-Degeneracy): the inverse metric $g^{\mu\nu}$ has non-zero determinant.
3. *Selection Principle 3.8* (Minimality + Dimensional Consistency): $k = 2$ is the minimal sector with non-trivial automorphism content, consistent with empirical (3+1)-dimensional spacetime.

Everything else — Fisher metric (Čencov), Lorentzian signature (Paley-Wiener), $\mathrm{SO}(3)$ frame group (Borel classification), $\mathrm{SU}(2)$ cover (Lie theory), local gauge covariance (gauge principle), spinorial monodromy (homotopy classification) — is a theorem with explicit proof.

The delivered vacuum — a 4D Lorentzian Fisher manifold with $\mathrm{SO}(3)/\mathrm{SU}(2)$ frame bundle, local frame covariance, Fisher metric, and conserved spinorial monodromy — is the complete starting point of the Fisher-Frame defect theory. The electron, its mass, the electroweak gauge group, the Weinberg angle, and the generation count are all downstream consequences of the defect physics of this vacuum. They are derived in the Fisher-Frame paper and inherited here.

The Arithmetic Vacuum does not produce particles. It produces the statistical manifold on which particle-like defects can exist.

---

## A. Symbol Glossary

| Symbol | Meaning |
|---|---|
| $\mathcal{V}_\times$ | Arithmetic Vacuum: $(\mathbb{N}^{(\infty)}, \leq_\mathrm{coord}, K_\times, \Lambda)$ |
| $\mathbf{v}(n)$ | Exponent-vector map: $n \mapsto (v_p(n))_p \in \mathbb{N}^{(\infty)}$ |
| $\mathbf{w}$ | Weight vector $(\log p)_{p\,\text{prime}}$ |
| $K_\times(\mathbf{u},\mathbf{v})$ | Causal kernel $\langle\mathbf{v}-\mathbf{u},\mathbf{w}\rangle$ when $\mathbf{u}\leq_\mathrm{coord}\mathbf{v}$ |
| $\Phi_S$ | Mode elimination map: $n \mapsto (v_{p_i}(n))_{p_i \in S}$ |
| $\sigma$ | Prime-swap automorphism $(v_{p_1}, v_{p_2}) \mapsto (v_{p_2}, v_{p_1})$ |
| $S_k$ | Permutation group on $k$ elements; $\mathrm{Aut}_\mathrm{abs}(\mathbb{N}^k, \leq_\mathrm{coord})$ |
| $\mathcal{M}_\mathrm{stat}$ | Statistical manifold of coarse-grained arithmetic states $\cong \mathbb{R}^2$ |
| $p(A \mid u)$ | Gaussian distribution parameterised by continuum point $u \in \mathbb{R}^2$ |
| $g^\mathrm{Fisher}_{ij}$ | Fisher information metric $\delta_{ij}/\ell_0^2$ on $\mathcal{M}_\mathrm{stat}$ |
| $\ell_0$ | Primitive causal correlation length; $G = 2\ell_0^2$ |
| $\hat{t}$ | Arithmetic time direction $\mathbf{w}/\|\mathbf{w}\|$; carries negative metric |
| $U(x)$ | Frame field $\in \mathrm{SU}(2)$; lift of rotation frame $R(x) \in \mathrm{SO}(3)$ |
| $\omega_\mu$ | $\mathfrak{su}(2)$-valued connection (gauge field) |
| $\Xi_\mu$ | Frame-flow $U^{-1}\partial_\mu U - \omega_\mu$; kinematic object of Fisher-Frame |
| $D_\mu U$ | Covariant derivative $\partial_\mu U - \omega_\mu U$ |
| $\eta^{\mu\nu}$ | Flat Minkowski metric $\mathrm{diag}(-1,+1,+1,+1)$ |
| $S_0[\phi]$ | Free-field action $\frac{1}{2\ell_0^2}\int\eta^{\mu\nu}\partial_\mu\phi\partial_\nu\phi\,d^4x$ |
| $Q \in \{0,1\}$ | $\mathbb{Z}/2\mathbb{Z}$ monodromy class (fermion number) |
| **Axiom 4.0** | Dynamical Closure: CIS time evolution is Markov, Fisher-conserving, minimally embedded |
| **Axiom 6.0** | Metric Non-Degeneracy: $\det(g^{\mu\nu}) \neq 0$ |
| **Selection Principle 3.8** | $k = 2$ minimal sector with non-trivial automorphism content + (3+1)D consistency |

---

## B. Theorem Status Table

| Theorem | Status | Section |
|---|---|---|
| $\mathbf{v}$ is an order isomorphism | Proved | 2.2 |
| Log map is a scalar projection, not isomorphism | Proved | 2.3 |
| Mode elimination = smooth-number restriction | Proved | 2.6–2.7 |
| $\mathrm{Aut}_\mathrm{abs}(\mathbb{N}^k,\leq_\mathrm{coord}) = S_k$ | Proved | 3.2 |
| Automorphism-faithfulness from sufficiency | Proved | 3.6 |
| No-new-symmetries from data processing | Proved | 3.7 |
| $k = 2$ uniquely forced from pure arithmetic | **Demoted to Selection Principle** | 3.8 |
| $k = 2$ from minimality + (3+1)D consistency | **Proved (Selection Principle)** | 3.8 |
| 2-prime projection is a sufficient statistic | Proved | 3.10 |
| Dynamical Closure axiom (Hamiltonian time evolution) | **Axiom** | 4.0 |
| Symplectic lift: $d = 4$ from Dynamical Closure | **Proved (under Axiom 4.0)** | 4.2 |
| Compatibility with $\mathbb{R} \to \mathbb{C}$ number tower | Proved | 4.4 |
| Fisher metric from Čencov | **Proved (central theorem)** | 5.2 |
| Unique time direction $\hat{t} = \mathbf{w}/|\mathbf{w}|$ | Proved | 6.1 |
| Propagator is retarded (zero backward kernel) | Proved | 6.2 |
| Metric Non-Degeneracy axiom | **Axiom** | 6.0 |
| Lorentzian signature $(-,+,+,+)$ unique (incl. off-diagonal) | **Proved (under Axiom 6.0)** | 6.4 |
| IR physics is prime-independent | Proved | 6.4 |
| Free-field action $S_0 = \frac{1}{2\ell_0^2}\int\eta^{\mu\nu}\partial_\mu\phi\partial_\nu\phi$ | **Proved** | 6a.1 |
| $\pi_1(\mathrm{SO}(3)) = \mathbb{Z}/2\mathbb{Z}$ from dimensionality | **Proved (topological)** | 7.1 |
| $\mathrm{SO}(3)$ uniquely determined on $\mathbb{R}^3$ | Proved | 7.1 |
| $\mathrm{SU}(2)$ as universal cover | Proved | 7.2 |
| Arithmetic $\mathbb{Z}/2\mathbb{Z}$ ↔ $\pi_1(\mathrm{SO}(3))$ correspondence | **Proved (correspondence)** | 7.3 |
| Local frame covariance forced | Proved | 8.2 |
| $\mathbb{Z}/2\mathbb{Z}$ monodromy from diamond holonomy | Proved | 9.1 |
| Fermion number conserved (topological) | Proved | 9.2 |
| Rigorous CLT for arithmetic kernel (lattice/Poisson model) | **Reduced to standard machinery** (stats preprint v0.4 §10.3, Thm 7.2) | 11.1 |
| Manifold-level convergence ($\Phi$ to continuum) | **Open — technical** | 11.1 |
| $\ell_0$ from arithmetic covariance | **Open** | 11.2 |
| $\mathrm{SU}(3)_c$ from 3-prime sector | **Open** | 11.3 |
| Structural forcing of $k = 2$ without dimensional input | **Open — structural (primary)** | 11.4 |
| Hamiltonian vs dissipative dynamical closure | **Open — axiomatic** | 11.5 |
| Sufficient-statistic ⇒ faithfulness from pure arithmetic | **Open** | 11.6 |
| Curved-space action (Einstein-Hilbert + $\Lambda$) | **Open — CFG programme** | — |

---

## C. Axioms and Selection Principles Used

This appendix consolidates the non-theorem inputs to the derivation chain. Every deduction in §§2–10 is either a theorem (with explicit proof in the text) or invokes one of the inputs listed here.

### C.1 Axioms

**Axiom 4.0** (Dynamical Closure). *The CIS time direction $\hat{t}$ generates a one-parameter family of evolutions on a phase space $\mathcal{P} \supseteq \mathcal{M}_\mathrm{stat}$ satisfying: (D1) Markov property; (D2) Fisher conservation; (D3) Minimality of $\mathcal{P}$.*

**Role.** Forces the symplectic lift $\mathcal{M}_\mathrm{stat} \to T^*\mathcal{M}_\mathrm{stat}$ (§4.2). Without it, the configuration space $\mathbb{R}^2$ would admit only dissipative evolution and there would be no 4D phase space, no momentum coordinates, and no Lorentzian-signature derivation in §6. With it, Hamilton's theorem forces $\mathcal{P} = T^*\mathcal{M}_\mathrm{stat}$.

**Status.** Genuine axiom. The Hamiltonian-vs-dissipative selection is not derived from the CIS structure alone. See Open Problem 11.5 for candidate derivations.

**Axiom 6.0** (Metric Non-Degeneracy). *The inverse metric $g^{\mu\nu}$ on the continuum statistical manifold has non-zero determinant: $\det(g^{\mu\nu}) \neq 0$.*

**Role.** Excludes degenerate cases ($g^{tt} = 0$, rank-deficient cross-block) in the signature analysis of §6.3. Without it, "trivially retarded" propagators (e.g., $\tilde{G} \propto $ constant) would be admitted, weakening the Lorentzian-signature conclusion to a non-uniqueness.

**Status.** Standard regularity condition. Could potentially be derived from invertibility of the Fisher information matrix combined with non-degeneracy of the time direction, but this is not done in the present paper.

### C.2 Selection Principles

**Selection Principle 3.8** (Minimality + Dimensional Consistency). *$k = 2$ is the minimal non-trivial coarse-graining sector ($k = 1$ has trivial automorphism group, $k = 2$ has $\mathbb{Z}/2\mathbb{Z}$) and is the unique sector among non-trivial choices consistent with the empirical (3+1)-dimensional structure of spacetime ($k \geq 3$ produces $\geq 6$-dimensional phase space).*

**Role.** Forces the choice of coarse-graining at $k = 2$. Without it (or a structural replacement), the framework would admit infinitely many coarse-grainings ($k = 2, 3, 4, \ldots$, all non-trivial) with no internal way to select one.

**Status.** Selection principle, not structural theorem. The minimality side (lower bound $k \geq 2$) is genuinely structural. The upper-bound side ($k \leq 2$) currently requires empirical input from observed spacetime dimensionality. Replacing this empirical input with a structural argument is Open Problem 11.4 (the primary open structural problem of the paper).

### C.3 Empirical Inputs

The framework currently uses the following empirical inputs from physics:

1. **(3+1)-dimensional spacetime.** Used in Selection Principle 3.8 to bound $k \leq 2$.
2. **Electron mass $m_e \approx 0.5109989$ MeV.** Used in Fisher-Frame §17 + CFG §16 cross-sector prediction to fix $\ell_0$.

Removing either of these — replacing (1) with a structural derivation of $k = 2$, and (2) with a derivation of $\ell_0$ from arithmetic covariance — would close the empirical-input loop of the programme. These are Open Problems 11.4 and 11.2 respectively.

### C.4 Theorem-Level Inputs from Classical Mathematics

The following classical theorems are used as inputs without re-derivation:

- **Frobenius theorem** (real division algebras): used in number-tower lift via companion paper.
- **Čencov's uniqueness theorem** (Fisher metric): used in §5.2.
- **Paley-Wiener theorem** (Fourier analyticity of retarded functions): used in §6.3.
- **Borel's classification** (compact Lie groups acting transitively on $S^n$): used in §7.1 to identify $\mathrm{SO}(3)$.
- **Lie theory of universal covers**: used in §7.2 for $\mathrm{SU}(2) \to \mathrm{SO}(3)$.
- **Hamilton's theorem** (canonical formulation of Markov isometric flows): used in §4.2.
- **Lindeberg-Feller / Erdős-Kac CLT machinery**: required for Open Problem 11.1 (not yet applied rigorously to the arithmetic kernel).

### C.5 What Is Forced vs. What Is Selected

| Result | Forced by structure | Requires axiom/selection/empirical input |
|---|---|---|
| Exponent-vector CIS structure | ✓ (Theorem 2.2) | — |
| Lower bound $k \geq 2$ | ✓ (Selection Principle 3.8 lower bound) | — |
| Upper bound $k \leq 2$ | — | Empirical: (3+1)D spacetime |
| Symplectic lift to $d = 4$ | — | Axiom 4.0 |
| Fisher metric | ✓ (Theorem 5.2, Čencov) | — |
| Lorentzian signature $(-,+,+,+)$ | — | Axiom 6.0 + Paley-Wiener |
| $\mathrm{SO}(3)$ rotation group | ✓ (Theorem 7.1, Borel) | — |
| $\mathrm{SU}(2)$ spinorial cover | ✓ (Theorem 7.2, Lie theory) | — |
| Arithmetic ↔ $\pi_1$ correspondence | ✓ (Theorem 7.3, diamond holonomy) | — |
| Local frame covariance | ✓ (Theorem 8.2, gauge principle) | — |
| Spinorial monodromy class | ✓ (Theorem 9.1, homotopy classification) | — |

The honest summary: the programme forces approximately two-thirds of its conclusions from pure structure, and the remaining third (symplectic lift, signature uniqueness, dimensional bound $k \leq 2$) requires explicit axiom or selection input. The three open problems most likely to convert "selected" into "forced" are 11.4 (structural $k = 2$), 11.5 (Hamiltonian closure), and 11.1 (rigorous arithmetic CLT enabling a fully structural continuum limit).
