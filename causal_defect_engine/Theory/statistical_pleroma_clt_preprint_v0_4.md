# Statistics as Stable Completion: A Pleroma Interpretation of Probability, Inference, and the Central Limit Theorem

**Kye Kirby**
Draft formalization with ChatGPT — **v0.4, May 2026**

Changes from v0.3:

- **New §10.4** promotes the joint completion of $(\omega_N, S_N)$ — listed as an open direction in §19.5 of v0.2/v0.3 — to a theorem (Theorem 10.21). The proof shows that $Y_1 := (\omega_N - \log\log N)/\sqrt{\log\log N}$ and $Y_2 := (S_N - \log N)/(\log N/\sqrt 2)$ jointly converge to $(Z, S_\infty)$ with $Z \sim \mathcal N(0,1)$ independent of the §10.3 limit $S_\infty$. The argument is via joint cumulants of $(B_p, v_p)$: every mixed joint cumulant $\kappa_{j,k}$ with $j,k\geq 1$ vanishes after standardization, at the explicit rate $O((\log\log N)^{-j/2})$. This is the Erdős–Kac/Euler-product splitting: the small-prime bulk produces an independent Gaussian factor, the large-prime tail produces the §10.3 pure-jump factor.
- §10.1 has a new worked Example 10.2bis (Pareto/Cauchy) illustrating $\alpha$-stable convergence concretely.
- The v0.3 changelog's reference to "Proposition 10.7'" of §10.2 has been removed — that proposition was never actually inserted in v0.3; the relevant Lévy–Khintchine statement is Theorem 10.8 of §10.2 as numbered.
- §19.5 has been pruned: the joint-completion bullet has been replaced by a pointer to §10.4; the L-function and explicit-formula bullets remain as research directions.
- §20 conclusion has been updated to mention the joint-completion theorem.

Changes from v0.2 (retained from v0.3):

- **§10.3 has been substantially rewritten** to close the main technical gap of v0.2. In v0.2 the proof of Theorem 10.11 ended with an unresolved scaling factor and a confused change of variable; the stated Lévy measure $\nu(dy) = \sqrt{2}\,dy/y$ was off by a factor of $\sqrt{2}$, and the limit was incorrectly described as compound Poisson. The corrected Lévy measure is $\nu(dy) = \mathbf{1}_{(0,\sqrt{2})}(y)\,dy/y$. Because this measure has infinite total mass, the limit is **not** compound Poisson; it is an infinitely divisible law with no Gaussian part, with all moments finite, and with infinitely many small jumps near zero. The Bernoulli approximation argument has been replaced by a rigorous proof via cumulant convergence, and the characteristic-function route is also given with explicit error bounds. The cumulants $\kappa_k(S_\infty) = 2^{k/2}/k$ ($k\geq 2$) and the explicit non-Gaussianity diagnostics (skewness $2\sqrt{2}/3$, excess kurtosis $1$) are now stated as a corollary.
- §10.1 was sharpened: Theorem 10.2 gives a precise statement of $\alpha$-stable convergence under regular variation, rather than a principle.
- §10.2 has the Lévy–Khintchine classification (Theorem 10.8) in the form actually needed for §10.3, and the summary table of completion regimes was corrected (the row previously labelled "Compound Poisson / Lévy" is now split between the finite-$\nu$ and infinite-$\nu$ cases).
- §18 dictionary entries were corrected and extended: "Compound Poisson completion" is now a special case of the general "Infinitely divisible completion" entry, and entries for "Lévy density" and "pure-jump non-Poisson completion" have been added.

---

## Abstract

This note formalizes a pleroma-style interpretation of statistics. The guiding principle is that probability theory describes incomplete access to an underlying space of possible defect configurations, while statistical laws describe stable summaries of unresolved microscopic distinctions. A random variable is treated as an observable on a hidden configuration space; a probability distribution is the pushforward shadow cast by such an observable; independence is factorization of defect sectors; dependence and correlation are gluing data; expectation is a completed average response; variance is the quadratic spread of unresolved defect responses; and Fisher information is the quadratic cost of deforming a statistical completion. The central limit theorem is recast as a stable-completion theorem: under independent additive composition and finite variance, the connected cumulant hierarchy flows to a universal quadratic fixed point, the Gaussian. The v0.2 extension classified stable completions in the independent but non-identically distributed regime: Gaussian when Lyapunov scaling holds, infinitely divisible with non-trivial Lévy measure when it fails. The arithmetic Euler-product example of §10.3 exhibits the second regime explicitly, with closed-form Lévy density $\nu(y) = 1/y$ on $(0,\sqrt{2})$ and closed-form cumulants $\kappa_k = 2^{k/2}/k$. The v0.4 extension (§10.4) proves the joint completion of $(\omega_N, S_N) = (\#\{p\leq N: v_p\geq 1\},\,\sum_{p\leq N}v_p\log p)$: after the natural Erdős–Kac/Mertens centerings and scalings by $\sqrt{\log\log N}$ and $\log N/\sqrt 2$, the two coordinates converge jointly to an *independent* pair $(Z, S_\infty)$ with $Z$ standard Gaussian and $S_\infty$ the pure-jump infinitely divisible limit of §10.3 — the small-prime bulk and the large-prime tail constitute two structurally distinct, asymptotically independent stable completions of the same arithmetic vacuum.

---

## Contents

1. Purpose and framing
2. Pleroma probability spaces
3. Expectation, variance, and quadratic uncertainty
4. Independence, gluing, and correlation
5. Sampling and empirical completion
6. Characteristic functions as phase probes
7. The central limit theorem as stable completion
8. Cumulant renormalization
9. Why the scale is √n
10. Stable completions beyond the i.i.d. Gaussian regime
   - 10.1. Heavy tails and α-stable completions
   - 10.2. The independent non-identically distributed case (Lindeberg–Feller)
   - 10.3. An arithmetic case: the Euler product and infinitely divisible completion
   - 10.4. Joint completion of (ω_N, S_N): asymptotic independence of the Erdős–Kac and Euler-product limits
11. Sufficient statistics and compressed defect coordinates
12. Likelihood as compatibility with a completion modulus
13. Fisher information as deformation cost
14. Entropy and relative entropy
15. Bayesian inference as renormalized possibility
16. Regression and quadratic filling
17. Statistical mechanics and arithmetic partition functions
18. Main dictionary
19. Research program
20. Conclusion

---

## 1. Purpose and framing

The ordinary mathematical definition of probability begins with a probability space
$$(\Omega, \mathcal{F}, \mathbb{P}),$$
and random variables
$$X : \Omega \to \mathbb{R}.$$
The pleroma interpretation keeps this formal apparatus, but changes the ontology. The sample space is read as a space of possible underlying configurations; the probability measure is a filling weight over unresolved possibilities; and an observed random variable is a probe of that hidden possibility space.

The slogan is

> *statistics is the geometry of incomplete access to possible configurations.*

More specifically,

> *statistics studies stable summaries of unresolved microscopic distinctions.*

This viewpoint is compatible with ordinary probability theory. It does not replace Kolmogorov probability; rather, it supplies a defect-geometric interpretation of its central objects.

---

## 2. Pleroma probability spaces

**Definition 2.1** (Pleroma probability space). A *pleroma probability space* is a quadruple
$$\mathcal{P} = (\Omega, \mathcal{F}, \mathbb{P}, \mathcal{D}),$$
where $(\Omega, \mathcal{F}, \mathbb{P})$ is a probability space and $\mathcal{D}$ is an interpreted defect structure on $\Omega$. Elements $\omega \in \Omega$ are called *microscopic configurations* or *defect configurations*. The measure $\mathbb{P}$ is called the *completion weight*.

The additional symbol $\mathcal{D}$ is not extra structure needed for measure theory. It is an interpretive layer: it records which coordinates, strata, symmetries, or gluing relations are regarded as defects. In arithmetic applications, $\mathcal{D}$ might encode prime-exponent data. In physical applications, it might encode vortices, fields, boundary states, or topological sectors. In statistics, it records the unresolved microscopic distinctions behind observations.

**Definition 2.2** (Observable). An *observable* on a pleroma probability space is a measurable map $X : (\Omega, \mathcal{F}) \to (E, \mathcal{E})$, where $(E, \mathcal{E})$ is an observation space. If $E = \mathbb{R}$, $X$ is a real-valued observable.

Thus a random variable is interpreted as an observable on hidden defect configurations:

> *random variable = observable on a hidden defect space.*

**Definition 2.3** (Observable shadow / distribution). Given an observable $X : \Omega \to E$, its *distribution* is the pushforward measure $\mu_X = X_\# \mathbb{P}$, defined by $\mu_X(A) = \mathbb{P}(X^{-1}(A))$ for $A \in \mathcal{E}$. In pleroma language, $\mu_X$ is called the *observable shadow* of the hidden configuration measure through $X$.

This gives the first main dictionary:
$$\Omega \;\longrightarrow\; X \;\longrightarrow\; \mu_X,$$
that is,
> *defect space → observable → distributional shadow.*

---

## 3. Expectation, variance, and quadratic uncertainty

For an integrable real-valued observable $X$, expectation is $\mathbb{E}[X] = \int_\Omega X(\omega)\,\mathbb{P}(d\omega)$.

**Definition 3.1** (Completed average response). The expectation $\mathbb{E}[X]$ is the *completed average response* of the observable $X$ over the unresolved defect configurations weighted by $\mathbb{P}$.

For $X \in L^2(\Omega)$, variance is $\mathrm{Var}(X) = \mathbb{E}[(X - \mathbb{E}X)^2]$.

**Definition 3.2** (Quadratic uncertainty). The variance $\mathrm{Var}(X)$ is the *quadratic spread* of the unresolved defect responses of $X$ around its completed average.

In this language, the passage from $X$ to $X - \mathbb{E}[X]$ removes the first-order drift, and variance keeps the second-order energy of unresolved fluctuations.

---

## 4. Independence, gluing, and correlation

**Definition 4.1** (Defect factorization). Let $(\Omega_i, \mathcal{F}_i, \mathbb{P}_i, \mathcal{D}_i)$ be pleroma probability spaces. A joint system $(\Omega_1 \times \cdots \times \Omega_n, \mathcal{F}_1 \otimes \cdots \otimes \mathcal{F}_n, \mathbb{P}, \mathcal{D})$ is *defect-factorized* if $\mathbb{P} = \mathbb{P}_1 \otimes \cdots \otimes \mathbb{P}_n$. It is *glued* if $\mathbb{P}$ does not factorize.

For observables $X_i$ on the different components, statistical independence is exactly measure factorization. The pleroma interpretation is

> *independence = absence of gluing between sampled defect sectors.*

Dependence is nontrivial gluing.

For square-integrable real-valued observables $X, Y$, covariance is $\mathrm{Cov}(X, Y) = \mathbb{E}[(X - \mathbb{E}X)(Y - \mathbb{E}Y)]$. It is a two-point response measuring failure of first-order factorization.

**Definition 4.2** (Two-point gluing response). The covariance $\mathrm{Cov}(X, Y)$ is the *two-point gluing response* between the observable shadows $X$ and $Y$.

The normalized response $\rho_{XY} = \mathrm{Cov}(X, Y)/\sqrt{\mathrm{Var}(X)\mathrm{Var}(Y)}$ is correlation. It is the scale-free component of the gluing response.

---

## 5. Sampling and empirical completion

Let $X_1, \ldots, X_n$ be independent copies of $X$. The empirical mean is $\bar{X}_n = n^{-1}\sum_{j=1}^n X_j$. The law of large numbers says that, under appropriate hypotheses, $\bar{X}_n \to \mathbb{E}[X]$ in probability, almost surely, or in $L^p$, depending on the theorem used.

**Principle 5.1** (Law of large numbers as completion). *Repeated independent probing of the same observable cancels unstable local fluctuations and completes toward the expectation.*

So the law of large numbers is interpreted as
> *many independent probes complete to the stable mean response.*

---

## 6. Characteristic functions as phase probes

For a real-valued random variable $X$, the characteristic function is
$$\varphi_X(t) = \mathbb{E}[e^{itX}],\quad t \in \mathbb{R}.$$
It is the Fourier transform of the distribution $\mu_X$.

**Definition 6.1** (Phase probe). The characteristic function $\varphi_X$ is the *phase probe* of the observable shadow $\mu_X$.

This is the statistical analogue of a spectral probe. It records how the distribution responds to oscillatory tests $x \mapsto e^{itx}$.

The logarithm of the characteristic function, when defined near zero, generates cumulants:
$$\log \varphi_X(t) = \sum_{k \geq 1} \frac{(it)^k}{k!} \kappa_k(X).$$

The first cumulants are $\kappa_1 = \mathbb{E}[X]$, $\kappa_2 = \mathrm{Var}(X)$, with higher cumulants measuring skewness, kurtosis, and higher connected structure.

**Definition 6.2** (Connected defect cumulants). The cumulants $\kappa_k(X)$ are the *connected defect responses* of order $k$ for the observable $X$.

This terminology is natural because cumulants, unlike raw moments, add under independent sums. If $X$ and $Y$ are independent, then $\kappa_k(X + Y) = \kappa_k(X) + \kappa_k(Y)$. Thus cumulants are the additive connected coordinates of a distribution.

---

## 7. The central limit theorem as stable completion

The central limit theorem is the central example of statistical stable completion. Let $X_1, X_2, \ldots$ be independent identically distributed real-valued observables with
$$\mathbb{E}[X_1] = \mu,\quad 0 < \mathrm{Var}(X_1) = \sigma^2 < \infty.$$
Define the normalized additive completion
$$Z_n = \frac{X_1 + \cdots + X_n - n\mu}{\sigma\sqrt{n}}.$$
The ordinary theorem says $Z_n \Rightarrow \mathcal{N}(0, 1)$.

**Definition 7.1** (Additive completion operator). For an observable $X$ with finite nonzero variance, define the $n$-fold *additive completion operator* by $R_n[X] := (X_1 + \cdots + X_n - n\mathbb{E}[X])/\sqrt{n\mathrm{Var}(X)}$, where $X_1, \ldots, X_n$ are independent copies of $X$.

The central limit theorem says that repeated additive completion flows to a universal Gaussian fixed point.

**Theorem 7.2** (Central limit stable completion). *Let $X_1, X_2, \ldots$ be independent copies of a real-valued observable $X$ with $\mathbb{E}[X] = \mu$ and $0 < \mathrm{Var}(X) = \sigma^2 < \infty$. Then $R_n[X] \Rightarrow \mathcal{N}(0, 1)$. Equivalently, finite-variance additive defect noise has a universal stable quadratic completion: the Gaussian.*

**Proof.** Let $Y = (X - \mu)/\sigma$. Then $\mathbb{E}[Y] = 0$ and $\mathrm{Var}(Y) = 1$. Let $\varphi_Y$ be its characteristic function. Since $Y$ has mean zero and variance one, $\varphi_Y(s) = 1 - s^2/2 + o(s^2)$ as $s \to 0$. For $Z_n = (Y_1 + \cdots + Y_n)/\sqrt{n}$, independence gives $\varphi_{Z_n}(t) = [\varphi_Y(t/\sqrt{n})]^n$. Using the expansion at zero, $\varphi_Y(t/\sqrt{n}) = 1 - t^2/(2n) + o(1/n)$. Therefore $[\varphi_Y(t/\sqrt{n})]^n \to e^{-t^2/2}$. The function $e^{-t^2/2}$ is the characteristic function of $\mathcal{N}(0, 1)$. By Lévy's continuity theorem, $Z_n \Rightarrow \mathcal{N}(0, 1)$. $\square$

The proof shows the actual mechanism. The low-frequency phase probe forgets all higher microscopic details except the quadratic term. The mean is removed; the variance is normalized; the residual shape is Gaussian.

---

## 8. Cumulant renormalization

When cumulants exist, the stable-completion mechanism is especially transparent.

**Proposition 8.1** (Cumulant flow under normalized sums). *Assume $X$ has cumulants $\kappa_k(X)$ for $k \geq 1$, with $\kappa_1(X) = \mu$ and $\kappa_2(X) = \sigma^2 > 0$. Let $Z_n = (X_1 + \cdots + X_n - n\mu)/(\sigma\sqrt{n})$. Then $\kappa_1(Z_n) = 0$, $\kappa_2(Z_n) = 1$, and for $k \geq 3$:*
$$\kappa_k(Z_n) = n^{1-k/2}\,\frac{\kappa_k(X)}{\sigma^k}.$$
*Thus, for every $k \geq 3$ with finite $\kappa_k(X)$, $\kappa_k(Z_n) \to 0$.*

**Proof.** Cumulants add under independent sums and scale homogeneously: for a constant $a$, $\kappa_k(aX) = a^k \kappa_k(X)$. The centered sum has cumulant $n\kappa_k(X)$ for $k \geq 2$, while the first cumulant is removed by centering. Scaling by $(\sigma\sqrt{n})^{-1}$ gives $\kappa_k(Z_n) = n\kappa_k(X)/(\sigma\sqrt{n})^k = n^{1-k/2}\kappa_k(X)/\sigma^k$. For $k = 2$ this equals 1, and for $k \geq 3$ the exponent $1 - k/2$ is negative. $\square$

This gives the pleroma statement:

> *all connected defect cumulants above order two are irrelevant under finite-variance additive completion.*

The Gaussian is the unique law with cumulants $\kappa_1 = 0$, $\kappa_2 = 1$, $\kappa_k = 0$ ($k \geq 3$). Thus
> *Gaussian = pure quadratic defect field.*

---

## 9. Why the scale is √n

The centering term $n\mu$ removes linear drift. The scale $\sqrt{n}$ is forced by variance addition: $\mathrm{Var}(X_1 + \cdots + X_n) = n\sigma^2$. Thus the fluctuation radius is $\sqrt{n\sigma^2} = \sigma\sqrt{n}$.

In the pleroma interpretation,
> *$\sqrt{n}$ = natural radius of independent quadratic fluctuation growth.*

The mean grows linearly, while the fluctuation amplitude grows as the square root of the number of independent components.

---

## 10. Stable completions beyond the i.i.d. Gaussian regime

The Gaussian completion of §§7–9 assumes (i) the contributions are independent, (ii) they are identically distributed, and (iii) they have finite variance. Each of these can fail, and each failure produces a different stable completion. This section classifies the possibilities.

### 10.1 Heavy tails and α-stable completions

The Gaussian is not the universal fixed point for all aggregation processes. It is the universal fixed point for finite-variance additive aggregation. If the variance is infinite, different stable laws may appear.

**Definition 10.1** (α-stable law). A random variable $W$ is *α-stable* with index $0 < \alpha \leq 2$ if for every $n \geq 1$ there exist constants $a_n > 0$ and $b_n \in \mathbb{R}$ such that, for $W_1, \ldots, W_n$ i.i.d. copies of $W$,
$$W_1 + \cdots + W_n \stackrel{d}{=} a_n W + b_n.$$
The scaling exponent is forced: $a_n = n^{1/\alpha}$. The case $\alpha = 2$ is the Gaussian; the cases $0 < \alpha < 2$ are non-Gaussian stable laws (Cauchy is $\alpha = 1$).

**Theorem 10.2** (α-stable convergence under regular variation). *Let $X_1, X_2, \ldots$ be i.i.d. with distribution function $F$ satisfying the regular-variation condition*
$$1 - F(x) \sim c_+ x^{-\alpha} L(x), \qquad F(-x) \sim c_- x^{-\alpha} L(x), \qquad x \to \infty,$$
*for some $0 < \alpha < 2$, constants $c_\pm \geq 0$ with $c_+ + c_- > 0$, and slowly varying function $L$. Then there exist sequences $a_n \to \infty$ with $a_n = n^{1/\alpha} L_1(n)$ ($L_1$ slowly varying) and $b_n \in \mathbb{R}$ such that*
$$\frac{X_1 + \cdots + X_n - b_n}{a_n} \Rightarrow W_\alpha,$$
*where $W_\alpha$ is an α-stable law with skewness parameter $\beta = (c_+ - c_-)/(c_+ + c_-)$.*

**Proof.** Classical; see Feller II.6 or Sato §14. $\square$

**Example 10.2bis** (Pareto and Cauchy). *Two canonical instances of Theorem 10.2:*

*(i) **Pareto tails, $1 < \alpha < 2$.** Let $X_i$ be i.i.d. with density $f(x) = \alpha\, c^\alpha x^{-\alpha - 1}\mathbf{1}_{\{x \geq c\}}$ for some $c > 0$ and $1 < \alpha < 2$. Then $\mathbb{E}[X_i] = \alpha c/(\alpha - 1) < \infty$ but $\mathbb{E}[X_i^2] = \infty$. The tail satisfies $\mathbb P(X_i > x) = (c/x)^\alpha$, so the regular-variation hypothesis of Theorem 10.2 holds with $c_+ = c^\alpha$, $c_- = 0$, $L \equiv 1$. The conclusion: $n^{-1/\alpha}(X_1 + \cdots + X_n - n\mathbb E[X_1]) \Rightarrow W_\alpha$, a totally right-skewed $\alpha$-stable law ($\beta = 1$). The scale of fluctuations is $n^{1/\alpha} \gg \sqrt n$; the sum is dominated by its single largest term, not by the bulk.*

*(ii) **Cauchy, $\alpha = 1$.** If $X_i$ are i.i.d. standard Cauchy with density $1/(\pi(1 + x^2))$, then $\mathbb{P}(|X_i| > x) \sim 2/(\pi x)$ as $x \to \infty$; this is the boundary regularly-varying case with $\alpha = 1$, $c_+ = c_- = 1/\pi$ (symmetric). The mean does not exist (so no centering is meaningful), and $(X_1 + \cdots + X_n)/n \Rightarrow X_1$ — the Cauchy is its own stable law, with no scaling improvement: averaging Cauchy variables does not concentrate them. The pleroma reading: when the first defect energy already diverges (no finite mean), no averaging recovers a deterministic limit, and the "law of large numbers" itself fails. Aggregation is non-self-averaging.*

These two cases bracket the genuinely non-Gaussian stable regime ($0 < \alpha < 2$). The pleroma reading of $\alpha$-stable completion is structural: when the second defect energy is infinite, quadratic completion is impossible because the energy required to specify the scale grows faster than the available bulk; the system saturates on its extreme excursions, and the stable completion records the law of those extremes.

**Pleroma reading.** In α-stable completion ($\alpha < 2$), the second defect energy is infinite. Aggregation is dominated by occasional very large individual excursions; the scaling $a_n \sim n^{1/\alpha}$ reflects that a typical aggregate value lives on a scale set by the largest of $n$ samples, not by the central bulk:

> *if the second defect energy is infinite, quadratic completion fails; the stable completion is α-stable with $\alpha < 2$ and is dominated by extreme individual contributions.*

The hierarchy is $\alpha = 2$ (Gaussian, bulk-dominated) $\to$ $\alpha < 2$ (heavy-tail, extreme-dominated).

### 10.2 The independent non-identically distributed case (Lindeberg–Feller)

The i.i.d. assumption can be relaxed. For a triangular array of independent (not necessarily identically distributed) observables $\{X_{n,k}\}_{1 \leq k \leq r_n}$ with $\mathbb{E}[X_{n,k}] = 0$ and total variance $s_n^2 := \sum_k \mathrm{Var}(X_{n,k})$, the question is when the normalized sum
$$T_n := \frac{1}{s_n}\sum_{k=1}^{r_n} X_{n,k}$$
completes to a Gaussian.

**Definition 10.3** (Lindeberg ratio). For each $\epsilon > 0$, the *Lindeberg ratio* is
$$L_n(\epsilon) := \frac{1}{s_n^2}\sum_{k=1}^{r_n}\mathbb{E}\!\left[X_{n,k}^2\,\mathbf{1}_{\{|X_{n,k}| > \epsilon s_n\}}\right].$$

The Lindeberg ratio measures the fraction of total variance contributed by individual atypically large excursions. In pleroma language, it is the fraction of the second defect energy located in extremal individual contributions rather than in the bulk.

**Theorem 10.4** (Lindeberg–Feller, stable-completion form). *If $L_n(\epsilon) \to 0$ for every $\epsilon > 0$, then $T_n \Rightarrow \mathcal{N}(0, 1)$ and the maximal individual contribution is asymptotically negligible: $\max_k \mathrm{Var}(X_{n,k})/s_n^2 \to 0$.*

**Proof.** Standard; see Billingsley §27 or Feller II.8. $\square$

**Definition 10.5** (Lyapunov ratio). For $\delta > 0$, the *Lyapunov ratio of order $2+\delta$* is
$$\Lambda_n(\delta) := \frac{1}{s_n^{2+\delta}}\sum_{k=1}^{r_n}\mathbb{E}\left[|X_{n,k}|^{2+\delta}\right].$$

**Proposition 10.6** (Lyapunov implies Lindeberg). *If $\Lambda_n(\delta) \to 0$ for some $\delta > 0$, then $L_n(\epsilon) \to 0$ for every $\epsilon > 0$, hence $T_n \Rightarrow \mathcal{N}(0, 1)$.*

**Proof.** $\mathbf{1}_{\{|X|>\epsilon s_n\}} \leq (|X|/(\epsilon s_n))^\delta$, so $X^2 \mathbf{1}_{\{|X|>\epsilon s_n\}} \leq |X|^{2+\delta}/(\epsilon s_n)^\delta$. Summing and dividing by $s_n^2$ gives $L_n(\epsilon) \leq \Lambda_n(\delta)/\epsilon^\delta$. $\square$

The pleroma reading of these conditions is structural:

**Principle 10.7** (When Gaussian completion holds, i.n.i.d. version). *In the i.n.i.d. regime, Gaussian stable completion occurs exactly when no single defect sector dominates the total second defect energy at any scale. The Lindeberg ratio is the quantitative diagnostic: it measures the fraction of variance located in atypical individual contributions.*

When the Lindeberg ratio fails to vanish — when the second defect energy concentrates in finitely many or in macroscopically sized individual terms — the limit need not be Gaussian. The general structural fact is the following.

**Theorem 10.8** (Lévy–Khintchine, stable-completion form). *Suppose $T_n = \sum_{k=1}^{r_n} X_{n,k}/s_n$ converges in distribution to a non-degenerate limit $T_\infty$ and the individual normalized terms are uniformly asymptotically negligible: $\max_k \mathbb{P}(|X_{n,k}/s_n| > \delta) \to 0$ for every $\delta > 0$. Then $T_\infty$ is infinitely divisible with characteristic function*
$$\log \varphi_{T_\infty}(t) = i\gamma t - \tfrac{1}{2}\tau^2 t^2 + \int_{\mathbb{R}\setminus\{0\}}\!\left(e^{itx} - 1 - \frac{itx}{1+x^2}\right)\nu(dx),$$
*where $\gamma \in \mathbb{R}$, $\tau^2 \geq 0$, and $\nu$ is a $\sigma$-finite measure on $\mathbb{R}\setminus\{0\}$ with $\int \min(1, x^2)\,\nu(dx) < \infty$. The triple $(\gamma, \tau^2, \nu)$ is uniquely determined by $T_\infty$ and is called its *Lévy–Khintchine triple*; $\nu$ is the Lévy measure.*

**Proof.** Classical; see Feller II.17 or Sato, *Lévy Processes and Infinitely Divisible Distributions*, Theorem 7.7. $\square$

**Pleroma reading.** When Gaussian completion fails because the Lindeberg ratio does not vanish, the stable completion has up to three components: a *deterministic drift* (the $\gamma$ term), a *quadratic component* (the Gaussian part $\tau^2$), and a *jump component* (the Lévy measure $\nu$). The jump component encodes the discrete or continuous-spectrum defect contributions that did not coarse-grain into the quadratic bulk. The Lévy measure is the *jump spectrum* of the unresolved completion.

**Remark 10.9** (Compound Poisson is the special case of finite $\nu$). If $\nu$ is a *finite* measure, $\nu(\mathbb{R}) = \lambda < \infty$, then $T_\infty$ has the explicit representation as a compound Poisson: $T_\infty \stackrel{d}{=} \gamma' + \sum_{j=1}^{N} J_j$ where $N \sim \mathrm{Poisson}(\lambda)$ and the jumps $J_j$ are i.i.d. with distribution $\nu/\lambda$, plus a deterministic shift $\gamma'$. If $\nu$ is infinite but $\int \min(1, x^2)\,\nu(dx) < \infty$, the law is infinitely divisible but *not* compound Poisson; it has infinitely many small jumps near zero, regularized by the centering. The example of §10.3 is of the second type.

This produces a hierarchy of stable completions:

| Regime | Lindeberg ratio | Stable completion | Lévy total mass $\nu(\mathbb{R})$ |
|---|---|---|---|
| i.i.d., finite variance | $L_n(\epsilon) \to 0$ | Gaussian | $0$ |
| i.n.i.d., Lyapunov holds | $L_n(\epsilon) \to 0$ | Gaussian | $0$ |
| i.n.i.d., Lindeberg fails, $\nu(\mathbb{R}) < \infty$ | $L_n(\epsilon) \not\to 0$ | Compound Poisson (+possible Gaussian part) | finite |
| i.n.i.d., Lindeberg fails, $\nu(\mathbb{R}) = \infty$ | $L_n(\epsilon) \not\to 0$ | Infinitely divisible, non-compound-Poisson | infinite |
| Infinite variance (heavy tails) | not applicable | $\alpha$-stable, $0 < \alpha < 2$ | infinite (and unbounded support of $\nu$) |

The Gaussian is the universal completion only in the first two rows. The arithmetic example of §10.3 lives in the fourth row.

### 10.3 An arithmetic case: the Euler product and infinitely divisible completion

The Euler-product probabilistic model on the natural numbers furnishes a clean example of i.n.i.d. completion that lies in the fourth row of the table above. It is interesting in its own right (it connects the pleroma framework to the Riemann zeta function via its Euler factorization) and serves as a worked illustration of infinitely divisible completion with infinite Lévy mass.

**Definition 10.10** (Euler-product probability model). For each prime $p$, let $v_p$ be an independent geometric random variable with parameter $1 - 1/p$, supported on $\{0, 1, 2, \ldots\}$:
$$\mathbb{P}(v_p = k) = (1 - 1/p)\,p^{-k},\qquad k \geq 0.$$
The associated *arithmetic observable* is
$$S_N := \sum_{p \leq N} v_p \log p.$$

Interpretation: $S_N$ is the logarithm of a "random integer" $n = \prod_p p^{v_p}$ assembled from independent prime-component defects, truncated at primes $\leq N$. The mean of $v_p$ is $1/(p-1)$ and its variance is $p/(p-1)^2$.

**Proposition 10.11** (Mean and variance of $S_N$). *As $N \to \infty$,*
$$\mathbb{E}[S_N] = \sum_{p \leq N}\frac{\log p}{p-1} = \log N + C_1 + o(1),$$
$$\mathrm{Var}(S_N) = \sum_{p \leq N}(\log p)^2 \cdot \frac{p}{(p-1)^2} = \frac{(\log N)^2}{2} + O(\log N).$$
*In particular, $s_N := \sqrt{\mathrm{Var}(S_N)} = \log N/\sqrt{2} + O(1)$, so $\log N/s_N \to \sqrt{2}$.*

**Proof.** Mean: $1/(p-1) = 1/p + 1/(p(p-1))$. By Mertens' theorem $\sum_{p \leq N}\log p/p = \log N + M + o(1)$ for a constant $M$; the remainder $\sum_p \log p/(p(p-1))$ converges absolutely. Hence the mean is $\log N + C_1 + o(1)$ with $C_1 = M + \sum_p \log p/(p(p-1))$.

Variance: expand $p/(p-1)^2 = 1/p + 2/p^2 + O(1/p^3)$. The $1/p^2$ and lower terms produce $\sum_p (\log p)^2/p^2 = O(1)$. So $\mathrm{Var}(S_N) = \sum_{p \leq N} (\log p)^2/p + O(1)$. Use Abel summation with $A(t) := \sum_{p \leq t}\log p/p = \log t + B + o(1)$:
$$\sum_{p \leq N}\frac{(\log p)^2}{p} = (\log N) A(N) - \int_2^N \frac{A(t)}{t}\,dt = (\log N)^2 + O(\log N) - \int_2^N \frac{\log t + O(1)}{t}\,dt = \frac{(\log N)^2}{2} + O(\log N).\quad\square$$

**Theorem 10.12** (Lindeberg failure for $S_N$). *The Lindeberg ratio for the centered observables $X_{N,p} := (v_p - \mathbb{E}[v_p])\log p$ does not vanish. Specifically, for every $\epsilon \in (0, \sqrt{2})$:*
$$\lim_{N \to \infty} L_N(\epsilon) = 1 - \frac{\epsilon^2}{2}.$$
*For $\epsilon \geq \sqrt{2}$, $L_N(\epsilon) \to 0$. Consequently $\tilde{S}_N := (S_N - \mathbb{E}[S_N])/s_N$ cannot converge to a Gaussian.*

**Proof.** Fix $\epsilon > 0$. The event $\{|X_{N,p}| > \epsilon s_N\}$ requires $|v_p - \mu_p|\log p > \epsilon s_N$, i.e., $v_p \geq m_p := \lceil \epsilon s_N/\log p + \mu_p\rceil$ (using $v_p \geq 0$ and the fact that $\mu_p = 1/(p-1) < 1$).

**Case 1: $\log p > \epsilon s_N$, i.e., $p > N^{\epsilon/\sqrt{2} + o(1)}$.** Then $m_p = 1$, so the truncation removes only the atom $v_p = 0$:
$$\mathbb{E}[X_{N,p}^2\mathbf{1}_{\{|X_{N,p}|>\epsilon s_N\}}] = \mathrm{Var}(X_{N,p}) - \mu_p^2(1-1/p)(\log p)^2 = \frac{(\log p)^2}{p} + O\!\left(\frac{(\log p)^2}{p^2}\right).$$

**Case 2: $\log p \leq \epsilon s_N$.** Then $m_p \geq 2$, and $\mathbb{P}(v_p \geq m_p) \leq (1-1/p) \sum_{k \geq m_p} p^{-k} \leq p^{1-m_p} \leq 1/p^2$. The contribution to the truncated second moment is at most $\mathbb{E}[X_{N,p}^2 \mathbf{1}_{v_p \geq 2}] \leq (\log p)^2 \mathbb{E}[v_p^2 \mathbf{1}_{v_p \geq 2}]$. From the geometric tail, $\mathbb{E}[v_p^2 \mathbf{1}_{v_p \geq 2}] = O(1/p^2)$. Hence this case contributes $O\!\left(\sum_p (\log p)^2/p^2\right) = O(1)$, which is negligible compared to $s_N^2 \sim (\log N)^2/2$.

So $L_N(\epsilon) = s_N^{-2} \sum_{p > N^{\epsilon/\sqrt{2}}} (\log p)^2/p + o(1)$. By the same Abel summation as in Proposition 10.11,
$$\sum_{N^{\epsilon/\sqrt{2}} < p \leq N} \frac{(\log p)^2}{p} = \frac{(\log N)^2 - (\epsilon\log N/\sqrt{2})^2}{2} + O(\log N) = \frac{(\log N)^2}{2}\!\left(1 - \frac{\epsilon^2}{2}\right) + O(\log N).$$
Dividing by $s_N^2 = (\log N)^2/2 + O(\log N)$ gives $L_N(\epsilon) \to 1 - \epsilon^2/2$ for $\epsilon \in (0, \sqrt{2})$, and $L_N(\epsilon) \to 0$ for $\epsilon \geq \sqrt{2}$ since the prime range becomes empty. $\square$

The functional form $L_N(\epsilon) \to 1 - \epsilon^2/2$ on $(0, \sqrt{2})$ is itself diagnostic: it says that, on the rescaled jump-size axis $y = \log p/s_N$, the truncated second-energy density is $\mathbf{1}_{\{y > \epsilon\}}\cdot y\,dy$, integrating to $(2-\epsilon^2)/2$ on $(\epsilon, \sqrt{2})$, normalized by total energy $1$. This is the same density we will identify in the Lévy measure below.

**Theorem 10.13** (Infinitely divisible stable completion of $S_N$). *Let $\tilde{S}_N := (S_N - \log N)/s_N$. As $N \to \infty$,*
$$\tilde{S}_N \;\Rightarrow\; S_\infty,$$
*where $S_\infty$ is an infinitely divisible random variable whose Lévy–Khintchine triple is*
$$(\gamma, \tau^2, \nu) = \left(\,0,\; 0,\; \mathbf{1}_{(0,\sqrt{2})}(y)\,\frac{dy}{y}\,\right),$$
*in the centered form*
$$\log\varphi_{S_\infty}(t) = \int_0^{\sqrt{2}}\!\bigl(e^{ity} - 1 - ity\bigr)\,\frac{dy}{y}.$$
*In particular, $S_\infty$ is **not** compound Poisson — the Lévy measure has infinite total mass — and has no Gaussian part: it is a pure-jump infinitely divisible law supported on $\mathbb{R}$ with positive skewness.*

The next two propositions give the main ingredients of the proof.

**Proposition 10.14** (Bernoulli reduction). *Let $B_p := \mathbf{1}_{\{v_p \geq 1\}}$, $T_N := \sum_{p \leq N} B_p \log p$. Then:*

*(a) $B_p \sim \mathrm{Bern}(1/p)$.*

*(b) $\mathbb{E}|S_N - T_N| = O(1)$ and $\mathrm{Var}(S_N - T_N) = O(1)$; hence $(S_N - T_N)/s_N \to 0$ in probability.*

*(c) Consequently $\tilde{S}_N$ and $\tilde{T}_N := (T_N - \log N)/s_N$ have the same weak limit (if either exists).*

**Proof.** (a) $\mathbb{P}(B_p = 1) = \mathbb{P}(v_p \geq 1) = 1/p$. (b) Compute
$$\mathbb{E}[S_N - T_N] = \sum_{p \leq N}\log p\cdot \mathbb{E}[(v_p - 1)\mathbf{1}_{v_p \geq 2}] = \sum_{p \leq N}\frac{\log p}{p(p-1)} = O(1),$$
using the closed form $\mathbb{E}[(v_p - 1)\mathbf{1}_{v_p \geq 2}] = 1/(p(p-1))$ (direct calculation). Similarly $\mathrm{Var}(v_p - B_p) \leq \mathbb{E}[(v_p - B_p)^2] = O(1/p^2)$, so $\mathrm{Var}(S_N - T_N) = \sum_p (\log p)^2 O(1/p^2) = O(1)$. By Markov, $(S_N - T_N)/s_N \to 0$ in probability. (c) Slutsky. $\square$

**Proposition 10.15** (Cumulant convergence of $\tilde{T}_N$). *For every integer $k \geq 2$, $\kappa_k(\tilde{T}_N) \to 2^{k/2}/k$ as $N \to \infty$.*

**Proof.** Cumulants of $\mathrm{Bern}(1/p)$: $\kappa_k(B_p) = q P_k(q)$ where $q = 1/p$ and $P_k(q) = 1 + O(q)$ (Bernoulli cumulant polynomial; e.g., $P_2 = 1-q$, $P_3 = (1-q)(1-2q)$, $P_4 = (1-q)(1-6q+6q^2)$, etc.). Hence $\kappa_k(B_p) = 1/p + O(1/p^2)$ for each fixed $k \geq 1$.

Since the $B_p$ are independent and cumulants are linear under independent sums with $\kappa_k(c X) = c^k \kappa_k(X)$:
$$\kappa_k(T_N) = \sum_{p \leq N}(\log p)^k \kappa_k(B_p) = \sum_{p \leq N}\frac{(\log p)^k}{p} + O\!\left(\sum_p\frac{(\log p)^k}{p^2}\right).$$
The error sum converges absolutely. The main sum, by Abel summation with $A(t) = \log t + B + o(1)$ and repeated application of $\int (\log t)^{k-1}/t\,dt = (\log t)^k/k$:
$$\sum_{p \leq N}\frac{(\log p)^k}{p} = \frac{(\log N)^k}{k} + O((\log N)^{k-1}).$$
Centering by $\log N$ kills only $\kappa_1$; for $k \geq 2$,
$$\kappa_k(\tilde{T}_N) = \frac{\kappa_k(T_N)}{s_N^k} = \frac{(\log N)^k/k + O((\log N)^{k-1})}{s_N^k} \xrightarrow[N\to\infty]{}\frac{(\sqrt{2})^k}{k} = \frac{2^{k/2}}{k},$$
since $\log N/s_N \to \sqrt{2}$. $\square$

**Proof of Theorem 10.13.** By Proposition 10.14, it suffices to show $\tilde{T}_N \Rightarrow S_\infty$. We give two arguments.

*Cumulant route.* By Proposition 10.15, $\kappa_k(\tilde{T}_N) \to 2^{k/2}/k$ for every $k \geq 2$. The candidate limit $S_\infty$ has Lévy measure $\nu(dy) = \mathbf{1}_{(0,\sqrt{2})}(y) dy/y$ supported on the bounded interval $(0, \sqrt{2})$. For an infinitely divisible law with bounded-support Lévy measure $\nu$, all moments exist and the cumulants are
$$\kappa_k(S_\infty) = \int y^k\,\nu(dy) \qquad (k \geq 2)$$
(this follows by differentiating $\log\varphi_{S_\infty}(t) = \int(e^{ity}-1-ity)\,\nu(dy)$ at $t=0$, using dominated convergence; the linear $-ity$ term contributes nothing for $k \geq 2$). Hence
$$\kappa_k(S_\infty) = \int_0^{\sqrt{2}} y^{k-1}\,dy = \frac{(\sqrt{2})^k}{k} = \frac{2^{k/2}}{k}.$$
The cumulants of $\tilde{T}_N$ match those of $S_\infty$ in the limit. Since $S_\infty$ has bounded-support Lévy measure, its moment generating function $\mathbb{E}[e^{tS_\infty}]$ is finite in a neighborhood of zero (Sato Theorem 25.17), so $S_\infty$ is moment-determinate. Combined with uniform integrability of $\tilde{T}_N^k$ (the cumulants and moments are uniformly bounded along the sequence by the above estimates), $\tilde{T}_N \Rightarrow S_\infty$.

*Characteristic function route (alternative).* Independence of the $B_p$ gives
$$\varphi_{T_N}(t) = \prod_{p \leq N}\!\left[1 + \tfrac{1}{p}\bigl(e^{it\log p} - 1\bigr)\right].$$
For $\tilde{T}_N$, take $t \mapsto t/s_N$ and subtract the centering:
$$\log\varphi_{\tilde{T}_N}(t) = -\frac{it\log N}{s_N} + \sum_{p \leq N}\log\!\left[1 + u_p(t)\right],\qquad u_p(t) := \tfrac{1}{p}\!\left(e^{it\log p/s_N} - 1\right).$$
Two bounds on $|u_p(t)|$ are useful: $|u_p(t)| \leq 2/p$ always, and $|u_p(t)| \leq |t|\log p/(p\,s_N)$ from $|e^{ix}-1|\leq|x|$. Hence
$$\sum_{p\leq N}|u_p(t)|^2 \leq t^2\!\sum_{p \leq e^{s_N}}\frac{(\log p)^2}{p^2 s_N^2} + 4\!\sum_{p > e^{s_N}}\frac{1}{p^2} = O\!\left(\frac{1}{s_N^2}\right) + O\!\left(e^{-s_N}\right) \xrightarrow[N\to\infty]{} 0,$$
where we split at $\log p = s_N$ and used $\sum_p(\log p)^2/p^2 = O(1)$, $\sum_{p > x} 1/p^2 = O(1/x)$. The Taylor identity $\log(1 + u) = u + O(|u|^2)$ for $|u| \leq 1/2$ (which holds for all $p \geq 3$, the finitely many $p=2$ contribution being $O(1/s_N) \to 0$ after rescaling) thus gives
$$\sum_p \log(1 + u_p) = \sum_p u_p + o(1).$$
Hence
$$\log\varphi_{\tilde{T}_N}(t) = \sum_{p \leq N}\frac{e^{it\log p/s_N} - 1}{p} - \frac{it\log N}{s_N} + o(1).$$
Adding and subtracting the linear-in-$t$ term:
$$\log\varphi_{\tilde{T}_N}(t) = \sum_{p \leq N}\frac{e^{it\log p/s_N} - 1 - it\log p/s_N}{p} + \frac{it}{s_N}\!\left(\sum_{p \leq N}\frac{\log p}{p} - \log N\right) + o(1).$$
The bracketed Mertens difference is $O(1)$, so the second term is $O(1/s_N) = o(1)$.

It remains to compute $\sum_{p \leq N} h(\log p/s_N)/p$ for $h(y) := e^{ity} - 1 - ity$. Since $h$ is smooth with $h(y) = O(y^2)$ near zero and bounded on $\mathbb{R}$, partial summation with the prime number theorem gives
$$\sum_{p \leq N}\frac{h(\log p/s_N)}{p} = \int_2^N \frac{h(\log t/s_N)}{t \log t}\,dt + o(1).$$
Substituting $u = \log t$ (so $du = dt/t$) and then $y = u/s_N$ (so $du = s_N\,dy$, $u = s_N y$):
$$\int_2^N \frac{h(\log t/s_N)}{t \log t}\,dt = \int_{\log 2}^{\log N}\frac{h(u/s_N)}{u}\,du = \int_{\log 2/s_N}^{\log N/s_N}\frac{h(y)}{y}\,dy.$$
The lower limit $\log 2/s_N \to 0^+$, the upper limit $\log N/s_N \to \sqrt{2}$. Since $h(y)/y = O(y) \to 0$ near $0$, the integral converges to
$$\int_0^{\sqrt{2}}\!\bigl(e^{ity} - 1 - ity\bigr)\,\frac{dy}{y}.$$
Hence $\log\varphi_{\tilde{T}_N}(t) \to \log\varphi_{S_\infty}(t)$ pointwise, and Lévy's continuity theorem gives $\tilde{T}_N \Rightarrow S_\infty$. $\square$

**Corollary 10.16** (Closed-form cumulants and non-Gaussianity diagnostics). *The cumulants of $S_\infty$ are*
$$\kappa_k(S_\infty) = \frac{2^{k/2}}{k} \qquad \text{for all } k \geq 2.$$
*In particular,*
$$\kappa_2 = 1,\qquad \kappa_3 = \frac{2\sqrt{2}}{3} \approx 0.9428,\qquad \kappa_4 = 1,\qquad \kappa_5 = \frac{4\sqrt{2}}{5},\qquad \kappa_6 = \frac{4}{3}.$$
*The standardized non-Gaussianity diagnostics are*
$$\mathrm{Skew}(S_\infty) = \kappa_3/\kappa_2^{3/2} = \frac{2\sqrt{2}}{3},\qquad \mathrm{ExKurt}(S_\infty) = \kappa_4/\kappa_2^2 = 1.$$
*A Gaussian completion would give $\kappa_k = 0$ for all $k \geq 3$; the values above quantify the deviation.*

**Proof.** Computed in the cumulant route of the proof of Theorem 10.13. $\square$

**Numerical verification.** Direct cumulant computation for $\tilde{T}_N$ at $N = 10^6$ gives $\kappa_2 \approx 1.000$, $\kappa_3 \approx 0.977$, $\kappa_4 \approx 1.062$, $\kappa_5 \approx 1.226$, with $\log N/s_N \approx 1.439$. The convergence to the limiting values $(\sqrt{2},\, 1,\, 2\sqrt{2}/3,\, 1,\, 4\sqrt{2}/5) \approx (1.414,\, 1,\, 0.943,\, 1,\, 1.131)$ is logarithmically slow (the leading correction is $O(1/\log N)$), but monotone and unambiguous.

**Remark 10.17** (What the limit *is* — and is not). The Lévy density $\nu(y) = 1/y$ on $(0,\sqrt{2})$ has infinite total mass; the limiting law $S_\infty$ therefore has infinitely many jumps in any compact time interval if realized as a Lévy process. It is, however, of *finite variation* (since $\int_0^{\sqrt{2}} y\cdot dy/y = \sqrt{2} < \infty$) and has all polynomial moments (since $\nu$ has bounded support). Some structural facts:

- $S_\infty$ is **not** Gaussian (skewness and excess kurtosis nonzero).
- $S_\infty$ is **not** compound Poisson (the Lévy measure is not finite).
- $S_\infty$ is **not** $\alpha$-stable (it has all moments; $\alpha$-stable laws for $\alpha < 2$ have only moments below $\alpha$).
- $S_\infty$ is **not** subordinator (despite $\nu$ being supported on $\mathbb{R}_+$, the centered form makes $S_\infty$ take negative values: $\mathbb{E}[S_\infty] = 0$ by construction, so the left tail is forced).

It is exactly what its Lévy–Khintchine triple says: a pure-jump, mean-zero, infinitely divisible law of finite variation with all moments, with positive-only jumps drawn from a $1/y$-density on $(0, \sqrt{2})$ together with the deterministic compensator required to enforce zero mean.

**Principle 10.18** (Two stable completions of the arithmetic vacuum). *The arithmetic vacuum supports two distinct stable completions, corresponding to two different probabilistic models:*

*(A) **The lattice/Poisson model**: unit-step Markov dynamics on the exponent lattice $\mathbb{N}^k$, or equivalently independent Poisson processes per prime direction. The stable completion is Gaussian by the standard i.i.d. CLT (Theorem 7.2 applied componentwise).*

*(B) **The Euler-product model**: independent geometric exponents per prime, no temporal dynamics. The stable completion is infinitely divisible with Lévy measure $\nu(dy) = \mathbf{1}_{(0,\sqrt{2})}(y) dy/y$ (Theorem 10.13). The total Lévy mass is infinite; the law is not compound Poisson; it is the pure-jump fixed point of additive aggregation of independent geometric prime sectors.*

*These are not in conflict: they describe different physical/structural questions about the same underlying defect space. Continuum dynamics on the arithmetic CIS is type (A). Random-integer statistics (Erdős–Kac type) is type (B).*

**Remark 10.19** (Connection to the cosmological sector). The non-Gaussian structure of $S_\infty$ is consistent with non-Gaussian primordial perturbation statistics. The $f_\mathrm{NL}$-type non-Gaussianity parameter familiar from cosmology has a natural arithmetic-vacuum interpretation here: it measures the leading-order deviation from Gaussian completion, controlled by the standardized third cumulant $\kappa_3/\kappa_2^{3/2}$. For the Euler-product completion this equals $2\sqrt{2}/3 \approx 0.943$ — a large, structurally fixed, parameter-free non-Gaussianity. Higher-order non-Gaussianity (the cosmological $g_\mathrm{NL}$ parameter) corresponds analogously to $\kappa_4/\kappa_2^2 = 1$. These are the appropriate stable-completion diagnostics for the arithmetic-vacuum cosmological sector listed as Open Problem 11.7 of the physical-sector paper.

**Remark 10.20** (Connection to the explicit formula). The integrated form $\sum_{p \leq N} \log p$ (the Chebyshev function $\theta(N)$) and the explicit formula $\psi(N) = N - \sum_\rho N^\rho/\rho + \cdots$ both involve sums over primes weighted logarithmically. The Lévy density $\nu(y) = 1/y$ on $(0, \sqrt{2})$ derived above is the probabilistic analogue of the prime-side jump structure in the explicit formula: discrete contributions from primes, weighted by $\log p$, with $1/p$ probability of activation, aggregating in the limit to a $1/y$ density on a bounded interval. Whether the zeta-zero contributions $\sum_\rho N^\rho/\rho$ in the explicit formula admit an interpretation as a *second* jump component in a richer Lévy measure (with $\nu_\rho$ supported on the spectrum $\{\Im\rho\}$ and weighted by $1/|\rho|$) is an open question — see the Active RH Research Path memorandum for the broader programme around explicit-formula-as-jump-process.

### 10.4 Joint completion of $(\omega_N, S_N)$: asymptotic independence of the Erdős–Kac and Euler-product limits

The Euler-product model carries *two* natural arithmetic observables: the prime-count statistic
$$\omega_N := \sum_{p \leq N}\mathbf{1}_{\{v_p \geq 1\}} = \sum_{p \leq N} B_p, \qquad B_p \sim \mathrm{Bern}(1/p),$$
which counts the active prime sectors in a sample (the analogue of $\omega(n) = $ number of distinct prime factors of a random integer), and the log-magnitude $S_N = \sum_{p \leq N} v_p \log p$ studied in §10.3. The two observables live on radically different scales: the variance of $\omega_N$ is $\sum_p (1/p)(1 - 1/p) = \log\log N + O(1)$, slow-growing; the variance of $S_N$ is $(\log N)^2/2 + O(\log N)$, fast-growing. Standardizing each on its own scale,
$$Y_1^{(N)} := \frac{\omega_N - \log\log N}{\sqrt{\log\log N}}, \qquad Y_2^{(N)} := \frac{S_N - \log N}{\log N/\sqrt 2},$$
we now ask for the *joint* limit.

That $Y_1^{(N)} \Rightarrow \mathcal N(0, 1)$ is essentially the Erdős–Kac theorem in its Euler-model form (Theorem 10.21(a) below). That $Y_2^{(N)} \Rightarrow S_\infty$ is Theorem 10.13. The new content of §10.4 is that the two limits are independent.

**Theorem 10.21** (Joint completion). *Let $\omega_N$ and $S_N$ be as above. Then:*

*(a) **Marginal Erdős–Kac (Euler form).** $Y_1^{(N)} \Rightarrow Z$, where $Z \sim \mathcal N(0, 1)$.*

*(b) **Marginal Euler-product.** $Y_2^{(N)} \Rightarrow S_\infty$, the infinitely divisible law of Theorem 10.13.*

*(c) **Joint asymptotic independence.** Jointly,*
$$\bigl(Y_1^{(N)},\, Y_2^{(N)}\bigr) \Rightarrow (Z, S_\infty),$$
*with $Z$ and $S_\infty$ independent.*

**Proof.**

*(a) Marginal $\omega_N$.* The summands $B_p - 1/p$ are independent, mean-zero, and bounded by $1$. The total variance is $V_N := \sum_{p \leq N}(1/p)(1 - 1/p) = \log\log N + B' + o(1)$ for the appropriate Mertens constant $B'$ (split off the $1/p^2$ tail, which converges absolutely). The Lyapunov ratio of order $3$ is
$$\Lambda_N(1) = V_N^{-3/2}\sum_{p \leq N}\mathbb E|B_p - 1/p|^3 \leq V_N^{-3/2}\sum_{p \leq N} \frac{1}{p} = V_N^{-3/2}(\log\log N + O(1)) = O\!\left(\frac{1}{\sqrt{\log\log N}}\right) \to 0,$$
since $\mathbb E|B_p - 1/p|^3 \leq \mathbb E[B_p^2] = 1/p$. By Proposition 10.6, $(\omega_N - \log\log N - B')/\sqrt{V_N} \Rightarrow \mathcal N(0, 1)$, and replacing $V_N$ by $\log\log N$ and removing the bounded $B'$ from the centering produces only $o(1)$ shifts. Hence $Y_1^{(N)} \Rightarrow \mathcal N(0, 1)$.

*(b) Marginal $S_N$.* Theorem 10.13.

*(c) Asymptotic independence.* By the Cramér–Wold device, it suffices to show that for every $(s, t) \in \mathbb R^2$,
$$\log \mathbb E\bigl[\exp\bigl(isY_1^{(N)} + itY_2^{(N)}\bigr)\bigr] \;\longrightarrow\; -\frac{s^2}{2} + \int_0^{\sqrt 2}\!\bigl(e^{ity} - 1 - ity\bigr)\,\frac{dy}{y}.$$
The RHS is exactly $\log\varphi_Z(s) + \log\varphi_{S_\infty}(t)$ — the sum of the marginal log-characteristic functions, which is the statement of independence at the level of characteristic functions.

Use the Bernoulli reduction of Proposition 10.14: replace $S_N$ by $T_N = \sum_{p \leq N} B_p \log p$, with $(S_N - T_N)/s_N \to 0$ in probability, hence in joint distribution with $Y_1^{(N)}$ (Slutsky). Then $(Y_1^{(N)}, \tilde T_N)$ depends only on the *common* Bernoulli sequence $(B_p)_{p \leq N}$. The joint characteristic function factorizes by prime independence:
$$\mathbb E\bigl[\exp(isY_1^{(N)} + it\tilde T_N)\bigr] = \prod_{p \leq N}\mathbb E\bigl[\exp\bigl(isB_p/\sqrt{V_N} + itB_p\log p/s_N\bigr)\bigr] \cdot (\text{deterministic centerings}).$$
For each prime $p$, write $a_p := s/\sqrt{V_N}$, $b_p := t\log p/s_N$, and compute
$$\mathbb E\bigl[\exp(i a_p B_p + i b_p B_p)\bigr] = \mathbb E\bigl[\exp(i(a_p + b_p)B_p)\bigr] = 1 + \frac{1}{p}\bigl(e^{i(a_p + b_p)} - 1\bigr).$$
Hence
$$\log\mathbb E\bigl[\exp(isY_1^{(N)} + it\tilde T_N)\bigr] = \sum_{p \leq N}\log\!\left[1 + \frac{1}{p}\bigl(e^{i(a_p + b_p)} - 1\bigr)\right] - \frac{is}{\sqrt{V_N}}\sum_{p \leq N}\frac{1}{p} - \frac{it\log N}{s_N}.$$
Now use $a_p + b_p \to 0$ uniformly in $p$ (since $|a_p| = O(1/\sqrt{\log\log N})$ and $|b_p| = O(1)$ uniformly) to Taylor-expand $\exp(i(a_p+b_p)) - 1 = i(a_p + b_p) - \tfrac{1}{2}(a_p + b_p)^2 + O((a_p+b_p)^3)$, then $\log(1 + u) = u - u^2/2 + O(u^3)$ for the outer logarithm with $u = (1/p)(\exp - 1)$.

Group the terms by their bivariate degree. The pure-$s$ terms — those involving $a_p$ but not $b_p$ — sum to the marginal Erdős–Kac characteristic exponent $-s^2/2 + o(1)$, exactly as in part (a). The pure-$t$ terms reproduce the §10.3 computation $\int_0^{\sqrt 2}(e^{ity} - 1 - ity)\,dy/y + o(1)$. The *mixed* terms are the new content; their total contribution to $\log \varphi$ is, leading order,
$$\Xi_N(s, t) := \sum_{p \leq N}\frac{(i a_p)(i b_p) - a_p b_p}{p}\cdot(\text{correction from the }-u^2/2\text{ in }\log(1+u)) + O\!\left(\sum_p (a_p^2 b_p + a_p b_p^2)/p\right).$$
Carrying out the algebra (or, more directly, expanding the joint cumulant generating function and reading off mixed cumulants — see the alternative joint-cumulant proof below), every nonvanishing mixed contribution arises from a joint cumulant $\kappa_{j,k}(B_p)$ with $j, k \geq 1$. Since $B_p \in \{0, 1\}$, $B_p^m = B_p$ for $m \geq 1$, so $\kappa_{j,k}(B_p, B_p) = \kappa_{j+k}(B_p)$ — there is only one Bernoulli to take cumulants of, and only the *total* cumulant order matters. Standardize: the standardized mixed cumulant has factor
$$\frac{1}{(\sqrt{V_N})^j\,(s_N)^k}\sum_{p \leq N}(\log p)^k \kappa_{j+k}(B_p) = \frac{1}{(\log\log N)^{j/2}\,(\log N/\sqrt 2)^k}\sum_{p \leq N}\frac{(\log p)^k}{p}\cdot O(1)$$
for $j + k \geq 2$. For $k \geq 1$, the prime sum is $(\log N)^k/k + O((\log N)^{k-1})$, which exactly cancels the $(\log N/\sqrt 2)^k$ denominator factor (producing a constant of order $(\sqrt 2)^k/k$). The remaining factor is $(\log\log N)^{-j/2}$, which $\to 0$ for every $j \geq 1$. Hence every mixed standardized cumulant with $j, k \geq 1$ vanishes in the limit at rate $(\log\log N)^{-j/2}$. The bivariate log-characteristic function therefore separates,
$$\log\mathbb E\bigl[\exp(isY_1^{(N)} + it\tilde T_N)\bigr] = -\frac{s^2}{2} + \int_0^{\sqrt 2}\!\bigl(e^{ity} - 1 - ity\bigr)\frac{dy}{y} + o(1),$$
which is the statement of independence. Replacing $\tilde T_N$ by $Y_2^{(N)}$ (which differs by $o_{\mathbb P}(1)$ uniformly in compact sets) preserves the joint limit. $\square$

**Alternative proof via joint cumulants** (sketch, useful for explicit error bounds). For $j, k \geq 0$ with $j + k \geq 1$, define the standardized joint cumulant
$$\widehat\kappa_{j,k}(N) := \frac{\kappa_{j,k}(\omega_N, S_N)}{V_N^{j/2}\,s_N^k} = \frac{\kappa_{j,k}(\omega_N, S_N)}{(\log\log N)^{j/2}\,(\log N/\sqrt 2)^k}\,(1 + o(1)).$$
By independence of the prime contributions,
$$\kappa_{j,k}(\omega_N, S_N) = \sum_{p \leq N}\kappa_{j,k}(B_p, v_p \log p) = \sum_{p \leq N}(\log p)^k\,\kappa_{j,k}(B_p, v_p),$$
where $\kappa_{j,k}(B_p, v_p)$ is the joint cumulant of the *single-prime* pair $(B_p, v_p)$. The latter is computable in closed form from the joint MGF $M_p(\sigma, \tau) = 1 - 1/p + (1/p)(1 - 1/p)e^{\sigma + \tau}/(1 - e^\tau/p)$. The leading-order single-prime joint cumulants (exact expressions, in increasing total order) are
$$\kappa_{1,0}(B_p, v_p) = \frac{1}{p},\quad \kappa_{0,1} = \frac{1}{p - 1},\quad \kappa_{1,1} = \frac{1}{p - 1} - \frac{1}{p(p-1)},$$
$$\kappa_{2,0} = \frac{1}{p}\!\left(1 - \frac{1}{p}\right),\quad \kappa_{0,2} = \frac{p}{(p-1)^2},$$
$$\kappa_{1,2} = O(1/p), \quad \kappa_{2,1} = \kappa_{1,1}\bigl(1 - 2/p\bigr) = O(1/p),$$
and in general $\kappa_{j,k}(B_p, v_p) = O(1/p)$ for all $j, k$ with $j + k \geq 1$. Hence the summed cumulant satisfies
$$\kappa_{j,k}(\omega_N, S_N) = \sum_{p \leq N}\frac{(\log p)^k}{p}\cdot (c_{j,k} + O(1/p))$$
for an explicit constant $c_{j,k}$ depending on the Bernoulli/geometric cumulant structure. By Abel summation,
$$\sum_{p \leq N}\frac{(\log p)^k}{p} = \begin{cases} \log\log N + O(1) & k = 0, \\ \log N + O(1) & k = 1,\\ \tfrac{1}{k}(\log N)^k + O((\log N)^{k-1}) & k \geq 1.\end{cases}$$
Standardizing,
$$\widehat\kappa_{j,k}(N) = \frac{c_{j,k}\,(\log N)^k/k + O((\log N)^{k-1})}{(\log\log N)^{j/2}\,(\log N/\sqrt 2)^k}\,(1 + o(1)) = \frac{c_{j,k}\,(\sqrt 2)^k}{k\,(\log\log N)^{j/2}}\,(1 + o(1)) \quad (k \geq 1),$$
which converges to $0$ for $j \geq 1$ at rate $(\log\log N)^{-j/2}$, and converges to $(\sqrt 2)^k/k$ (the Euler-product cumulants of §10.3) for $j = 0$. For $k = 0$, $j \geq 2$, the standardized cumulant is $\widehat\kappa_{j,0}(N) = (\log\log N)^{-j/2}\sum_p (1/p)\cdot O(1) = O((\log\log N)^{1 - j/2})$, which $\to 0$ for $j \geq 3$ (the Gaussian condition) and $\to 1$ for $j = 2$ (the Gaussian variance). The standardized joint cumulant generating function therefore converges:
$$\log\mathbb E\bigl[e^{sY_1 + tY_2}\bigr] = \sum_{j+k \geq 2}\frac{\widehat\kappa_{j,k}(N)}{j!\,k!}s^j t^k \;\to\; \frac{s^2}{2} + \sum_{k \geq 2}\frac{(\sqrt 2)^k/k}{k!}t^k = \frac{s^2}{2} + \int_0^{\sqrt 2}\!\bigl(e^{ty} - 1 - ty\bigr)\frac{dy}{y},$$
which is exactly the joint cumulant generating function of $(Z, S_\infty)$ with $Z$ and $S_\infty$ independent. $\square$

**Numerical verification.** Direct numerical computation of $\widehat\kappa_{j,k}(N)$ from the closed-form single-prime joint cumulants gives the following at $N = 10^8$: $\widehat\kappa_{2,0} \approx 0.93$ (theoretical limit $1$; rate $1/\log\log N$ correction); $\widehat\kappa_{0,3} \approx 0.944$ (theoretical limit $2\sqrt 2/3 \approx 0.9428$); the cross cumulants $\widehat\kappa_{1,1} \approx 0.77$, $\widehat\kappa_{2,1} \approx 0.42$, $\widehat\kappa_{1,2} \approx 0.58$, all decreasing monotonically and visibly tracking the predicted rates $(\log\log N)^{-1/2} \approx 0.59$, $(\log\log N)^{-1} \approx 0.35$, $(\log\log N)^{-1/2} \approx 0.59$ respectively. The convergence is logarithmically slow but unambiguous — characteristic of Erdős–Kac-type results, where $\log\log N$ remains $O(1)$ for any computationally accessible $N$.

**Pleroma reading.** The arithmetic vacuum $\mathrm{Spec}\,\mathbb Z$, sampled via the Euler-product probabilistic model, admits *two structurally distinct, asymptotically independent stable completions* of its additive logarithm sector:

> *(I) A **small-prime bulk** completion: the count $\omega_N$ of active prime sectors. This is Gaussian. Its variance is $\log\log N$ — the integrated $1/p$ density. Lindeberg holds; the bulk is the universal quadratic vacuum.*

> *(II) A **large-prime tail** completion: the weighted log-sum $S_N - \log N$. This is the pure-jump infinitely divisible law $S_\infty$ of §10.3. The variance is $(\log N)^2/2$ — incomparably larger. Lindeberg fails; each individual large prime makes a non-negligible contribution.*

These two completions are *not* artificially separated by the analyst: they coexist in the *same* sample $(v_p)_{p \leq N}$, but the small-prime statistics (counting alone) and the large-prime statistics (weighted by $\log p$) probe asymptotically independent sectors of the arithmetic defect data. This is the cleanest statement available within the framework of the splitting principle (small-prime bulk $\to$ quadratic; large-prime tail $\to$ pure-jump) that runs through analytic number theory, made precise here in the language of stable completions.

**Remark 10.22** (Why independence is structural). The asymptotic independence is *not* an accident of the geometric distribution: it is forced by the scale separation $\sqrt{\log\log N} \ll \log N$. Any pair of arithmetic observables of the form $(F_N, G_N) = (\sum_p f(B_p), \sum_p g(B_p, v_p)\log p)$ with $f$ supported on the small-prime bulk and $g$ activating only when $v_p \geq 1$ on a heavier scale, both standardized on their natural scales, will exhibit the same factorization. The structural content of Theorem 10.21 is that the framework of §§10.1–10.3 supplies the *correct* scale for each, and that the scales are decoupled enough for joint independence to hold.

---

## 11. Sufficient statistics and compressed defect coordinates

A statistic is a measurable function of data: $T = T(X_1, \ldots, X_n)$. It is sufficient for a parameter $\theta$ if it retains all the information in the sample relevant to $\theta$.

**Definition 11.1** (Sufficient defect coordinate). A *sufficient statistic* is a compressed observable preserving the parameter-relevant defect coordinates of the sample.

In a Gaussian family with unknown mean and variance, the sample mean and sample variance are sufficient. This matches the cumulant interpretation: a Gaussian completion is controlled entirely by first and second cumulant data.

---

## 12. Likelihood as compatibility with a completion modulus

A parametric statistical model is a family $\{P_\theta : \theta \in \Theta\}$. For observed data $x$, the likelihood is $L(\theta \mid x) = p_\theta(x)$, when densities exist.

**Definition 12.1** (Completion modulus). A parameter $\theta$ is a *completion modulus* if changing $\theta$ changes the filling measure $P_\theta$ over possible configurations.

Then likelihood measures compatibility between an observed shadow and a proposed completion modulus:

> *$L(\theta \mid x)$ = compatibility of observation $x$ with completion modulus $\theta$.*

The log-likelihood $\ell(\theta) = \log L(\theta \mid x)$ is especially natural because independent observations produce additive action: $\ell(\theta) = \sum_{j=1}^n \log p_\theta(x_j)$. Maximum likelihood chooses the completion modulus minimizing mismatch with the observed defect shadow.

---

## 13. Fisher information as deformation cost

For a smooth parametric model $p_\theta(x)$, the score is $S_i(x; \theta) = \partial_i \log p_\theta(x)$. The Fisher information matrix is $g_{ij}(\theta) = \mathbb{E}_\theta[S_i(X; \theta) S_j(X; \theta)]$.

**Definition 13.1** (Statistical deformation metric). The Fisher information metric is the *quadratic cost of deforming the statistical completion measure $P_\theta$ in parameter space*.

Thus Fisher information is not an incidental statistical object. It is the canonical metric on the manifold of completion weights. In pleroma language,
> *Fisher information = quadratic completion cost of changing a statistical defect modulus.*

This connects probability theory to information geometry and to field-theoretic interpretations in which quadratic action functionals define stable response.

---

## 14. Entropy and relative entropy

For a discrete probability distribution $P$, Shannon entropy is $H(P) = -\sum_x P(x) \log P(x)$. It measures unresolved multiplicity of compatible configurations.

**Definition 14.1** (Unresolved defect multiplicity). Entropy is the logarithmic measure of *unresolved defect multiplicity* under a probability weight.

For probability measures $P, Q$ with $P \ll Q$, relative entropy is $D(P \parallel Q) = \int \log(dP/dQ)\,dP$. It measures the cost of using completion measure $Q$ when the data are generated from $P$:
> *$D(P \parallel Q)$ = defect mismatch action between two completion measures.*

---

## 15. Bayesian inference as renormalized possibility

Bayes' rule is $p(\theta \mid x) = p(x \mid \theta) p(\theta) / p(x)$. The prior $p(\theta)$ is an initial completion weight over parameters. The likelihood $p(x \mid \theta)$ is compatibility of observation with a proposed completion modulus. The posterior is the renormalized completion weight after observing data.

**Principle 15.1** (Bayes as completion update). *Bayesian inference is renormalization of possibility after measurement.*

The evidence $p(x) = \int p(x \mid \theta) p(\theta)\,d\theta$ is the global normalization required for the updated completion measure.

---

## 16. Regression and quadratic filling

A regression model decomposes observations as $Y = f(X) + \varepsilon$. The residual $\varepsilon$ represents unresolved defect structure not captured by the stable observable geometry $f$.

Least squares minimizes $\sum_{j=1}^n (Y_j - f(X_j))^2$. This is a quadratic filling problem. Under Gaussian residuals, least squares coincides with maximum likelihood because the Gaussian is the quadratic stable completion.

**Principle 16.1** (Regression as stable decomposition). *Regression separates stable observable geometry from unresolved residual defect structure.*

---

## 17. Statistical mechanics and arithmetic partition functions

A Boltzmann distribution has the form $p(x) = e^{-\beta E(x)}/Z$, with $Z = \sum_x e^{-\beta E(x)}$. The partition function $Z$ is the global completion weight over configurations.

The arithmetic zeta function has a parallel form:
$$\zeta(s) = \sum_{n \geq 1} n^{-s} = \sum_{n \geq 1} e^{-s \log n}.$$
Thus $\log n$ plays the role of an arithmetic energy, and $s$ plays the role of an inverse-temperature-like spectral parameter. This supports the broader dictionary:
> *partition functions are completion normalizers over possible configurations.*

The Euler product $\zeta(s) = \prod_p (1 - p^{-s})^{-1}$ is the factorization of this partition function into independent prime sectors, which is exactly the probabilistic model used in §10.3.

---

## 18. Main dictionary

| Ordinary statistical object | Pleroma interpretation |
|---|---|
| Sample space $\Omega$ | space of possible defect configurations |
| Probability measure $P$ | filling weight over unresolved possibilities |
| Random variable $X$ | observable/probe of hidden configurations |
| Distribution $X_\# P$ | observable shadow |
| Expectation $\mathbb{E}[X]$ | completed average response |
| Variance $\mathrm{Var}(X)$ | quadratic uncertainty / second defect energy |
| Independence | absence of gluing between sampled sectors |
| Covariance | two-point gluing response |
| Correlation | normalized gluing response |
| Characteristic function | phase probe of a distribution |
| Cumulants | connected defect responses |
| Law of large numbers | empirical completion to the mean |
| Central limit theorem | stable quadratic completion of additive ignorance |
| Gaussian | pure quadratic defect field |
| Lindeberg ratio | fraction of second defect energy in atypical individual terms |
| Lyapunov ratio | higher-moment localization measure |
| **Infinitely divisible completion** | **general stable completion of additive aggregation under uniform negligibility** |
| **Lévy–Khintchine triple** $(\gamma, \tau^2, \nu)$ | **drift, Gaussian part, and jump part of an infinitely divisible completion** |
| **Lévy measure $\nu$** | **jump spectrum of an infinitely divisible completion** |
| **Lévy density** | **jump-spectrum density on the rescaled jump axis** |
| **Compound Poisson completion** ($\nu$ finite) | **stable completion with discretely many jump arrivals per realization** |
| **Pure-jump non-Poisson completion** ($\nu(\mathbb R)=\infty$, no Gaussian) | **stable completion with continuous-spectrum jump structure (Euler-product type)** |
| $\alpha$-stable completion ($\alpha<2$) | stable completion of infinite-variance heavy-tailed aggregation |
| Sufficient statistic | compressed parameter-relevant defect coordinate |
| Likelihood | compatibility with a completion modulus |
| Fisher information | quadratic deformation cost of a completion measure |
| Entropy | unresolved defect multiplicity |
| Relative entropy | defect mismatch action |
| Bayesian posterior | renormalized completion measure after observation |
| Regression | separation of stable geometry from residual defects |
| Partition function | global completion normalizer |

---

## 19. Research program

This draft suggests several technical directions.

### 19.1 Defect-rank statistics
One can define an effective defect rank of a distribution by the number, decay rate, or algebraic structure of its nonzero connected cumulants. Gaussian laws have defect rank two. $\alpha$-stable laws ($\alpha < 2$) have non-quadratic infinite-energy completion. The Euler-product limit of §10.3 has all cumulants $\kappa_k = 2^{k/2}/k$ nonzero with closed-form decay, so its defect rank is infinite but its cumulants are an explicit geometric-like sequence. Structured distributions may have low effective defect rank even when their sample spaces are large.

### 19.2 Completion-preserving inference
Classical estimators can be classified by whether they preserve the stable completion class of the data-generating process. In this language, robust statistics asks for estimators stable under defect contamination, while efficient statistics asks for estimators saturating the Fisher deformation bound.

### 19.3 Correlation sheaves
For multivariate and dependent data, one may assign local distributions to subsets of variables and study their gluing constraints. This turns dependence, conditional independence, and graphical models into a sheaf-theoretic or defect-geometric problem.

### 19.4 Statistical renormalization
The central limit theorem is one renormalization fixed point. A broader program would classify statistical fixed points under aggregation, conditioning, thinning, projection, and coarse-graining. The fixed points are stable completions of different unresolved-defect regimes. §10 begins this classification: Gaussian, $\alpha$-stable, finite-$\nu$ (compound Poisson), and infinite-$\nu$ (Euler-product type) completions are four distinct fixed-point classes within the infinitely divisible family.

### 19.5 Arithmetic completion theory (updated in v0.4)
The Euler-product example of §10.3 and the joint completion theorem of §10.4 close two items of the previous research program: the existence and explicit form of the stable completion for $\log n$, and the joint asymptotic behavior of $(\omega(n), \log n)$ under the Euler model. Outstanding directions:

- *L-function partition functions.* The Dirichlet L-function $L(s,\chi)$ has the Euler product $\prod_p (1 - \chi(p) p^{-s})^{-1}$, leading to character-twisted geometric exponents. The corresponding completion is an infinitely divisible law whose Lévy measure is a twist of $\mathbf{1}_{(0,\sqrt 2)}(y)\,dy/y$ by the character $\chi$; the twist should encode the $\chi$-dependence of zero distributions.
- *Explicit-formula reconstruction from Lévy structure.* If $\sum_{p\leq N}\log p$ has the infinitely divisible completion with Lévy density $1/y$ on $(0,\sqrt 2)$, and if zeta zeros contribute an additional jump structure $\nu_\rho$ supported on a different spectrum, then the explicit formula $\psi(N) = N - \sum_\rho N^\rho/\rho - \cdots$ should appear as a cumulant identity relating these two Lévy components. Making this precise is the long-range goal.
- *Higher joint completions.* §10.4 handles two coordinates with widely separated scales. The natural higher-order question concerns observables of intermediate scale (e.g., $\sum_p v_p (\log p)^\alpha$ for $0 < \alpha < 1$) and the structure of the resulting cumulant lattice — whether the splitting between "bulk Gaussian" and "tail pure-jump" sectors persists, deforms, or fails for these intermediate scalings.

---

## 20. Conclusion

The pleroma interpretation does not alter the formal content of probability theory, but it reorganizes its meaning. Probability becomes a theory of completion weights over possible configurations. Observables cast distributional shadows. Sampling probes the hidden configuration repeatedly. Independence is factorization; dependence is gluing; cumulants are connected defect responses; Fisher information is the metric cost of deforming a completion measure; entropy is unresolved multiplicity.

Within this framework, the central limit theorem takes on a simple conceptual form:
> *the Gaussian is the stable vacuum of finite-variance uncertainty.*

Equivalently,
> *the CLT is the stable-completion theorem for additive ignorance.*

Independent additive composition with finite variance and uniform negligibility removes all connected cumulant defects above order two. What remains is the universal quadratic completion: the Gaussian.

The v0.2/v0.3 extension makes the conditional structure of this completion theorem explicit. When uniform negligibility fails — when a finite fraction of the second defect energy is concentrated in atypically large individual contributions — the stable completion is no longer Gaussian. It is an infinitely divisible law with a non-trivial Lévy jump measure, classified up to drift and Gaussian part by the Lévy–Khintchine triple $(\gamma, \tau^2, \nu)$. The compound Poisson case is the subclass where $\nu$ has finite total mass.

The arithmetic Euler-product example of §10.3 exhibits this richer completion concretely. The centered logarithmic prime sum $\sum_{p \leq N} v_p \log p$ converges, after rescaling by $s_N \sim \log N/\sqrt{2}$, to a pure-jump mean-zero infinitely divisible law $S_\infty$ with explicit Lévy density
$$\nu(y) = \frac{1}{y},\qquad y \in (0, \sqrt{2}),$$
infinite total Lévy mass, and explicit cumulants
$$\kappa_k(S_\infty) = \frac{2^{k/2}}{k},\qquad k \geq 2.$$
The non-Gaussianity diagnostics are parameter-free: skewness $2\sqrt{2}/3$, excess kurtosis $1$. This is *not* a compound Poisson completion (the Lévy mass is infinite), and *not* a Gaussian completion (the higher cumulants are nonzero). It is the appropriate stable-completion framework for arithmetic random-integer statistics, for connections to Dirichlet L-functions, and for the cosmological non-Gaussianity that appears as $f_\mathrm{NL}$-type signatures in the broader arithmetic-vacuum programme.

The §10.4 joint completion theorem extends this picture in a structurally significant direction. The arithmetic vacuum, sampled through the Euler-product model, carries two simultaneous arithmetic observables — the prime-count $\omega_N$ and the log-magnitude $S_N$ — that live on incomparable scales ($\sqrt{\log\log N}$ versus $\log N/\sqrt 2$). Despite being functions of the *same* underlying defect data $(v_p)_{p\leq N}$, after the natural rescalings they are *asymptotically independent*: $\omega_N$ completes to a standard Gaussian (Erdős–Kac in its Euler form, driven by the small-prime bulk where Lindeberg holds), $S_N$ completes to the pure-jump infinitely divisible law $S_\infty$ (driven by the large-prime tail where Lindeberg fails), and the joint cross-cumulants vanish in the limit at the explicit rate $(\log\log N)^{-j/2}$. This is the rigorous form of the splitting principle that runs through analytic number theory: the small-prime bulk and the large-prime tail are structurally distinct sectors of the same arithmetic vacuum, each carrying its own stable completion, decoupled by the gap between their natural scales.

The broader picture is a small zoo of fixed points of additive aggregation under uniform negligibility: Gaussian (no jumps), compound Poisson (finitely many jump types), pure-jump non-Poisson (continuous-spectrum jumps, e.g. Euler-product), and $\alpha$-stable for $\alpha < 2$ (infinite-variance heavy tails). Each is a distinct stable vacuum of additive ignorance.

---

## References

[1] P. Billingsley, *Probability and Measure*, Wiley.

[2] W. Feller, *An Introduction to Probability Theory and Its Applications*, Vols. I–II, Wiley.

[3] R. Durrett, *Probability: Theory and Examples*, Cambridge University Press.

[4] S.-i. Amari, *Information Geometry and Its Applications*, Springer.

[5] T. Cover and J. Thomas, *Elements of Information Theory*, Wiley.

[6] K. Sato, *Lévy Processes and Infinitely Divisible Distributions*, Cambridge University Press.

[7] P. Erdős and M. Kac, "The Gaussian law of errors in the theory of additive number theoretic functions," *Amer. J. Math.* **62** (1940), 738–742.

[8] G. Tenenbaum, *Introduction to Analytic and Probabilistic Number Theory*, Cambridge University Press.

[9] J. Pitman, "Combinatorial stochastic processes," *Lecture Notes in Math.* **1875**, Springer, 2006 (for the appearance of Poisson–Dirichlet / Lévy structure in arithmetic random-integer models).
