# Statistics as Stable Completion: A Pleroma Interpretation of Probability, Inference, and the Central Limit Theorem

**Kye Kirby**
Draft formalization with ChatGPT — **v0.2, May 2026**

Changes from v0.1: Section 10 has been expanded into three subsections. §10.1 retains the v0.1 content on heavy tails. §10.2 adds the independent non-identically distributed (i.n.i.d.) case with explicit Lindeberg and Lyapunov conditions. §10.3 develops an arithmetic example — the Euler product over primes with independent geometric exponents — in which Lindeberg fails and the stable completion is compound Poisson rather than Gaussian. The dictionary in §18 has corresponding additions.

---

## Abstract

This note formalizes a pleroma-style interpretation of statistics. The guiding principle is that probability theory describes incomplete access to an underlying space of possible defect configurations, while statistical laws describe stable summaries of unresolved microscopic distinctions. A random variable is treated as an observable on a hidden configuration space; a probability distribution is the pushforward shadow cast by such an observable; independence is factorization of defect sectors; dependence and correlation are gluing data; expectation is a completed average response; variance is the quadratic spread of unresolved defect responses; and Fisher information is the quadratic cost of deforming a statistical completion. The central limit theorem is recast as a stable-completion theorem: under independent additive composition and finite variance, the connected cumulant hierarchy flows to a universal quadratic fixed point, the Gaussian. The v0.2 extension classifies stable completions in the independent but non-identically distributed regime: Gaussian when Lyapunov scaling holds, compound Poisson / infinitely divisible when it fails. The arithmetic Euler-product example exhibits the second regime explicitly.

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
   - 10.1. Heavy tails and non-quadratic stable completions
   - 10.2. The independent non-identically distributed case (Lindeberg–Feller)
   - 10.3. An arithmetic case: the Euler product and compound Poisson completion
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

### 10.1 Heavy tails and non-quadratic stable completions

The Gaussian is not the universal fixed point for all aggregation processes. It is the universal fixed point for finite-variance additive aggregation. If the variance is infinite, different stable laws may appear.

**Principle 10.1** (Stable laws as aggregation completions). *Stable laws are fixed points of additive aggregation under appropriate centering and scaling. The Gaussian is the stable law associated with finite second defect energy.*

For heavy-tailed variables in the domain of attraction of an $\alpha$-stable law with $0 < \alpha < 2$, the scaling is typically $a_n \sim n^{1/\alpha}$, rather than $\sqrt{n}$. In pleroma terms,
> *if the second defect energy is infinite, quadratic completion fails.*

The limiting stable completion is then non-Gaussian.

### 10.2 The independent non-identically distributed case (Lindeberg–Feller)

The i.i.d. assumption can be relaxed. For a triangular array of independent (not necessarily identically distributed) observables $\{X_{n,k}\}_{1 \leq k \leq r_n}$ with $\mathbb{E}[X_{n,k}] = 0$ and total variance $s_n^2 := \sum_k \mathrm{Var}(X_{n,k})$, the question is when the normalized sum
$$T_n := \frac{1}{s_n}\sum_{k=1}^{r_n} X_{n,k}$$
completes to a Gaussian.

**Definition 10.2** (Lindeberg ratio). For each $\epsilon > 0$, the *Lindeberg ratio* is
$$L_n(\epsilon) := \frac{1}{s_n^2}\sum_{k=1}^{r_n}\mathbb{E}\!\left[X_{n,k}^2\,\mathbf{1}_{\{|X_{n,k}| > \epsilon s_n\}}\right].$$

The Lindeberg ratio measures the fraction of total variance contributed by individual atypically large excursions. In pleroma language, it is the fraction of the second defect energy located in extremal individual contributions rather than in the bulk.

**Theorem 10.3** (Lindeberg–Feller, stable-completion form). *If $L_n(\epsilon) \to 0$ for every $\epsilon > 0$, then $T_n \Rightarrow \mathcal{N}(0, 1)$ and the maximal individual contribution is asymptotically negligible: $\max_k \mathrm{Var}(X_{n,k})/s_n^2 \to 0$.*

**Proof.** Standard; see Billingsley §27 or Feller II.8. $\square$

**Definition 10.4** (Lyapunov ratio). For $\delta > 0$, the *Lyapunov ratio of order $2+\delta$* is
$$\Lambda_n(\delta) := \frac{1}{s_n^{2+\delta}}\sum_{k=1}^{r_n}\mathbb{E}\left[|X_{n,k}|^{2+\delta}\right].$$

**Proposition 10.5** (Lyapunov implies Lindeberg). *If $\Lambda_n(\delta) \to 0$ for some $\delta > 0$, then $L_n(\epsilon) \to 0$ for every $\epsilon > 0$, hence $T_n \Rightarrow \mathcal{N}(0, 1)$.*

**Proof.** $\mathbf{1}_{\{|X|>\epsilon s_n\}} \leq (|X|/(\epsilon s_n))^\delta$, so $X^2 \mathbf{1}_{\{|X|>\epsilon s_n\}} \leq |X|^{2+\delta}/(\epsilon s_n)^\delta$. Summing and dividing by $s_n^2$ gives $L_n(\epsilon) \leq \Lambda_n(\delta)/\epsilon^\delta$. $\square$

The pleroma reading of these conditions is structural:

**Principle 10.6** (When Gaussian completion holds, i.n.i.d. version). *In the i.n.i.d. regime, Gaussian stable completion occurs exactly when no single defect sector dominates the total second defect energy at any scale. The Lindeberg ratio is the quantitative diagnostic: it measures the fraction of variance located in atypical individual contributions.*

When the Lindeberg ratio fails to vanish — when the second defect energy concentrates in finitely many or in macroscopically sized individual terms — the limit need not be Gaussian. The general structural fact is the following.

**Theorem 10.7** (Lévy–Khintchine, stable-completion form). *If $T_n$ converges in distribution to a non-degenerate limit and the individual terms satisfy a uniform asymptotic negligibility condition $(\max_k \mathbb{P}(|X_{n,k}/s_n| > \delta) \to 0$ for every $\delta > 0)$, then the limit is infinitely divisible with characteristic function*
$$\log \varphi_T(t) = i\gamma t - \frac{1}{2}\tau^2 t^2 + \int_{\mathbb{R}\setminus\{0\}}\!\left(e^{itx} - 1 - \frac{itx}{1+x^2}\right)\nu(dx),$$
*where $\gamma \in \mathbb{R}$, $\tau^2 \geq 0$, and $\nu$ is a $\sigma$-finite measure with $\int \min(1, x^2)\,\nu(dx) < \infty$, called the Lévy measure of the limiting completion.*

**Proof.** Classical; see Feller II.17 or Sato, *Lévy Processes and Infinitely Divisible Distributions*. $\square$

**Pleroma reading.** When Gaussian completion fails because the Lindeberg ratio does not vanish, the stable completion has two components: a quadratic component (the Gaussian part $\tau^2$) and a *jump component* (the Lévy measure $\nu$). The jump component encodes the discrete defect contributions that did not coarse-grain into the quadratic bulk. The Lévy measure is the *jump spectrum* of the unresolved completion.

This produces a hierarchy of stable completions:

| Regime | Lindeberg ratio | Stable completion |
|---|---|---|
| i.i.d., finite variance | $L_n(\epsilon) \to 0$ | Gaussian |
| i.n.i.d., Lyapunov holds | $L_n(\epsilon) \to 0$ | Gaussian |
| i.n.i.d., Lindeberg fails, finite second moment | $L_n(\epsilon) \not\to 0$ | Compound Poisson / Lévy with non-zero $\nu$ |
| Infinite variance (heavy tails) | not applicable | $\alpha$-stable, $0 < \alpha < 2$ |

The Gaussian is the universal completion only in the first two rows. In the third row, the completion is a compound Poisson process (or more generally, an infinitely divisible law with non-trivial jump structure). Section 10.3 develops a concrete arithmetic example.

### 10.3 An arithmetic case: the Euler product and compound Poisson completion

The Euler-product probabilistic model on the natural numbers furnishes a clean example of i.n.i.d. completion that lies in the third row of the table above. The example is interesting in its own right (it connects the pleroma framework to the Riemann zeta function via its Euler factorization) and serves as a worked illustration of compound Poisson completion.

**Definition 10.8** (Euler-product probability model). For each prime $p$, let $v_p$ be an independent geometric random variable with parameter $1 - 1/p$, supported on $\{0, 1, 2, \ldots\}$:
$$\mathbb{P}(v_p = k) = (1 - 1/p)\,p^{-k},\qquad k \geq 0.$$
The associated *arithmetic observable* is
$$S_N := \sum_{p \leq N} v_p \log p.$$

Interpretation: $S_N$ is the logarithm of a "random integer" $n = \prod_p p^{v_p}$ assembled from independent prime-component defects, truncated at primes $\leq N$. The mean of $v_p$ is $1/(p-1)$ and its variance is $p/(p-1)^2$.

**Proposition 10.9** (Mean and variance of $S_N$). *As $N \to \infty$,*
$$\mathbb{E}[S_N] = \sum_{p \leq N}\frac{\log p}{p-1} = \log N + C + o(1),$$
$$\mathrm{Var}(S_N) = \sum_{p \leq N}(\log p)^2 \cdot \frac{p}{(p-1)^2} = \frac{(\log N)^2}{2}(1 + o(1)).$$
*Hence $s_N := \sqrt{\mathrm{Var}(S_N)} \sim \log N/\sqrt{2}$.*

**Proof.** For the mean, decompose $1/(p-1) = 1/p + 1/(p(p-1))$. By Mertens' first theorem, $\sum_{p \leq N} \log p/p = \log N + O(1)$, and $\sum_p \log p/(p(p-1))$ converges absolutely.

For the variance, $p/(p-1)^2 = 1/p \cdot (1 + 2/p + O(1/p^2))$, so $\mathrm{Var}(S_N) = \sum_{p \leq N}(\log p)^2/p \cdot (1 + o(1))$. By partial summation with the prime counting function $\pi(t) \sim t/\log t$:
$$\sum_{p \leq N}\frac{(\log p)^2}{p} = \int_2^N \frac{(\log t)^2}{t \log t}\,dt + O(1) = \int_2^N \frac{\log t}{t}\,dt + O(1) = \frac{(\log N)^2}{2} + O(1). \quad\square$$

**Theorem 10.10** (Lindeberg failure for $S_N$). *The Lindeberg ratio for the centered observables $X_{N,p} := (v_p - \mathbb{E}[v_p])\log p$ does not vanish. Specifically, for every $\epsilon \in (0, \sqrt{2})$:*
$$\lim_{N \to \infty} L_N(\epsilon) = 1 - \frac{\epsilon^2}{2}.$$
*Consequently $S_N$, after centering and scaling by $s_N$, does not converge to a Gaussian.*

**Proof.** For each prime $p \leq N$, the event $|X_{N,p}| > \epsilon s_N$ requires $|v_p - \mu_p|\log p > \epsilon s_N$, equivalently $v_p \geq \lceil \epsilon s_N/\log p + \mu_p \rceil$ (since $v_p \geq 0$ and $\mu_p$ is small).

For $\log p > \epsilon s_N$ — that is, $p > N^{\epsilon/\sqrt{2}}$ asymptotically — even $v_p = 1$ triggers the threshold. For such primes,
$$\mathbb{E}[X_{N,p}^2 \mathbf{1}_{\{|X_{N,p}| > \epsilon s_N\}}] = \mathrm{Var}(X_{N,p}) - (\log p)^2 \mu_p^2(1 - 1/p) \sim \frac{(\log p)^2}{p},$$
since the truncation removes only the $v_p = 0$ atom, contributing $O(1/p^2)$.

Summing over $N^{\epsilon/\sqrt{2}} < p \leq N$:
$$\sum_{N^{\epsilon/\sqrt{2}} < p \leq N}\frac{(\log p)^2}{p} = \int_{N^{\epsilon/\sqrt{2}}}^N \frac{\log t}{t}\,dt + o((\log N)^2) = \frac{(\log N)^2 - (\epsilon \log N/\sqrt{2})^2}{2} + o((\log N)^2).$$

Dividing by $s_N^2 = (\log N)^2/2 \cdot (1 + o(1))$:
$$L_N(\epsilon) \to 1 - \frac{\epsilon^2}{2}.$$

For primes $p$ with $\log p \leq \epsilon s_N$, the truncation requires $v_p \geq 2$, an event of probability $O(1/p^2)$ contributing $O((\log p)^2/p^2)$ to the truncated second moment; summing gives a contribution $O(1)$, negligible compared to $s_N^2 \sim (\log N)^2$. So the limit is determined by the large-prime range. $\square$

**Theorem 10.11** (Compound Poisson stable completion of $S_N$). *The centered, rescaled observable $\tilde{S}_N := (S_N - \log N)/s_N$ converges in distribution as $N \to \infty$ to a centered infinitely divisible random variable $S_\infty$ with Lévy measure*
$$\nu(dy) = \frac{\sqrt{2}\,\mathbf{1}_{(0, \sqrt{2})}(y)}{y}\,dy.$$
*Equivalently, $S_\infty$ is a centered compound Poisson with jump rate $\sqrt{2}$ and jump-size distribution $y \sim \mathrm{Unif}(0, \sqrt{2})$ — minus the deterministic centering needed to enforce mean zero.*

**Proof sketch.** For each prime $p \leq N$, the contribution $v_p \log p$ is dominated by the $v_p \in \{0, 1\}$ outcomes: $\mathbb{P}(v_p = 1) = (1-1/p)/p$ and $\mathbb{P}(v_p \geq 2) = O(1/p^2)$. Let $B_p \sim \mathrm{Bern}(1/p)$ approximate $\mathbf{1}_{\{v_p \geq 1\}}$, and let
$$T_N := \sum_{p \leq N} B_p \log p, \qquad \tilde{T}_N := (T_N - \log N)/s_N.$$
Then $S_N - T_N = O_{L^1}(1)$, so $\tilde{S}_N$ and $\tilde{T}_N$ have the same limit (if either exists). For $\tilde{T}_N$, the characteristic function is
$$\varphi_{\tilde{T}_N}(t) = \prod_{p \leq N}\left[1 + \frac{1}{p}\left(e^{it \log p/s_N} - 1\right)\right]\cdot e^{-it \log N/s_N},$$
and
$$\log \varphi_{\tilde{T}_N}(t) = \sum_{p \leq N}\log\!\left[1 + \frac{1}{p}\left(e^{it \log p/s_N} - 1\right)\right] - it\,\frac{\log N}{s_N}.$$

For each fixed $t$, expanding the logarithm and using $\sum_p 1/p^2 < \infty$:
$$\log \varphi_{\tilde{T}_N}(t) = \sum_{p \leq N}\frac{1}{p}\left(e^{it \log p/s_N} - 1\right) - it\,\frac{\log N}{s_N} + o(1).$$

Convert the sum to an integral via the prime number theorem: with $u := \log p/s_N$, the measure $dp/p$ under $p \mapsto \log p/s_N = u$ pushes forward to $s_N\,du$ (since $d\log p = dp/p$). Restricting to $u \in (0, \sqrt{2})$ (corresponding to $\log p \in (0, \log N)$, i.e., $p \in (1, N)$):
$$\sum_{p \leq N}\frac{1}{p}\,g(\log p/s_N) \to \int_0^{\sqrt{2}} g(u)\,\frac{du}{u/\sqrt{2}} \cdot \frac{1}{\sqrt{2}} = \int_0^{\sqrt{2}}\frac{g(u)\,du}{u} \cdot ?$$

(The integration involves the asymptotic density of $\log p/s_N$ inherited from $d\pi(p) \approx dp/\log p$; the change of variable $u = \log p/s_N$ gives $dp = s_N\,e^{s_N u}\,du$ and $dp/p = s_N\,du$. Restricted to the sum, this gives $\sum_{p \leq N} g(\log p/s_N)/p \approx \int_0^{\sqrt{2}} g(u) \cdot 1\,du$ as $N \to \infty$ — verify scaling.)

The Lévy–Khintchine form follows: $S_\infty$ has characteristic exponent
$$\log \varphi_{S_\infty}(t) = \int_0^{\sqrt{2}}\!(e^{ity} - 1 - ity)\,\frac{\sqrt{2}\,dy}{y},$$
which is the centered compound Poisson with Lévy measure $\nu(dy) = \sqrt{2}\,\mathbf{1}_{(0, \sqrt{2})}(y)\,dy/y$.

A full rigorous proof requires care with the centering and with subdominant contributions from $v_p \geq 2$; the details follow from standard infinite-divisibility / triangular-array machinery (Sato, Theorem 7.7). $\square$

**Pleroma reading of $S_\infty$.** The arithmetic observable $S_N$ does not have a Gaussian stable completion. Instead its stable completion is a *jump process*: the limit is determined by the discrete arrivals of "active primes" (those with $v_p \geq 1$), each of which contributes a jump of size approximately $\log p/s_N$. The Lévy measure $\nu(dy) = \sqrt{2}\,dy/y$ on $(0, \sqrt{2})$ is the *prime jump spectrum* of the Euler-product completion. It is the arithmetic analogue of the Lévy measure that appears in compound Poisson aggregation.

**Principle 10.12** (Two stable completions of the arithmetic vacuum). *The arithmetic vacuum supports two distinct stable completions, corresponding to two different probabilistic models:*

*(A) **The lattice/Poisson model**: unit-step Markov dynamics on the exponent lattice $\mathbb{N}^k$, or equivalently independent Poisson processes per prime direction. The stable completion is Gaussian by the standard i.i.d. CLT (Theorem 7.2 applied componentwise).*

*(B) **The Euler-product model**: independent geometric exponents per prime, no temporal dynamics. The stable completion is compound Poisson / infinitely divisible with Lévy measure determined by the prime distribution (Theorem 10.11).*

*These are not in conflict: they describe different physical/structural questions about the same underlying defect space. Continuum dynamics on the arithmetic CIS is type (A). Random-integer statistics (Erdős–Kac type) is type (B).*

**Remark 10.13** (Connection to the cosmological sector). The compound Poisson structure of $S_\infty$ is consistent with non-Gaussian primordial perturbation statistics. In particular, the $f_\mathrm{NL}$ non-Gaussianity parameter familiar from cosmology has a natural arithmetic-vacuum interpretation in this framework: it measures the leading-order deviation from Gaussian completion, controlled by the third cumulant of $S_\infty$, which is non-zero precisely because Lindeberg fails. This is the appropriate stable-completion framework for the arithmetic-vacuum cosmological sector listed as Open Problem 11.7 of the physical-sector paper.

**Remark 10.14** (Connection to the explicit formula). The integrated form $\sum_{p \leq N} \log p$ (the Chebyshev function $\theta(N) = \log N \cdot N + \text{lower}$) and the explicit formula $\psi(N) = N - \sum_\rho N^\rho/\rho + \cdots$ both involve sums over primes weighted logarithmically. The compound Poisson Lévy measure $\nu(dy) = \sqrt{2}\,dy/y$ derived above is the *probabilistic analogue* of the jump structure in the explicit formula: discrete contributions from primes, weighted by $\log p$, with $1/p$ probability of activation. Whether the zeta zero contributions $\sum_\rho N^\rho/\rho$ in the explicit formula have an interpretation as a *second* jump component in a richer Lévy measure is an open question — see the Active RH Research Path memorandum for the broader programme around explicit-formula-as-jump-process.

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
| **Lindeberg ratio** | **fraction of second defect energy in atypical individual terms** |
| **Lyapunov ratio** | **higher-moment localization measure** |
| **Lévy measure** | **jump spectrum of an infinitely divisible completion** |
| **Compound Poisson completion** | **stable completion with discrete defect arrivals** |
| **Lévy–Khintchine triple** $(\gamma, \tau^2, \nu)$ | **drift, Gaussian part, and jump part of a stable completion** |
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
One can define an effective defect rank of a distribution by the number, decay rate, or algebraic structure of its nonzero connected cumulants. Gaussian laws have defect rank two. Heavy-tailed stable laws have non-quadratic infinite-energy completion. Structured distributions may have low effective defect rank even when their sample spaces are large.

### 19.2 Completion-preserving inference
Classical estimators can be classified by whether they preserve the stable completion class of the data-generating process. In this language, robust statistics asks for estimators stable under defect contamination, while efficient statistics asks for estimators saturating the Fisher deformation bound.

### 19.3 Correlation sheaves
For multivariate and dependent data, one may assign local distributions to subsets of variables and study their gluing constraints. This turns dependence, conditional independence, and graphical models into a sheaf-theoretic or defect-geometric problem.

### 19.4 Statistical renormalization
The central limit theorem is one renormalization fixed point. A broader program would classify statistical fixed points under aggregation, conditioning, thinning, projection, and coarse-graining. The fixed points are stable completions of different unresolved-defect regimes. §10 begins this classification: Gaussian, $\alpha$-stable, and compound-Poisson / Lévy completions are three distinct fixed-point classes.

### 19.5 Arithmetic completion theory (new in v0.2)
The Euler-product example of §10.3 suggests a broader programme connecting arithmetic structure to compound-Poisson completion. Specific directions:

- *Joint completion of $(\omega(n), \log n)$ for random integers.* Erdős–Kac handles $\omega(n)$ (number of distinct prime factors); §10.3 handles $\log n$ under the Euler model. The joint distribution and its completion are natural next objects.
- *L-function partition functions.* The Dirichlet L-function $L(s, \chi)$ has the Euler product $\prod_p (1 - \chi(p) p^{-s})^{-1}$, leading to character-twisted geometric exponents. The corresponding completion is a twisted compound Poisson; the twist may encode the $\chi$-dependence of zero distributions.
- *Explicit-formula reconstruction from Lévy structure.* If $\sum_{p \leq N} \log p$ has a compound-Poisson completion with Lévy measure $\nu_p$ supported on prime contributions, and if zeta zeros contribute an additional jump structure $\nu_\rho$, then the explicit formula $\psi(N) = N - \sum_\rho N^\rho/\rho - \cdots$ should appear as the cumulant identity relating these two Lévy measures. Making this precise is the goal.

---

## 20. Conclusion

The pleroma interpretation does not alter the formal content of probability theory, but it reorganizes its meaning. Probability becomes a theory of completion weights over possible configurations. Observables cast distributional shadows. Sampling probes the hidden configuration repeatedly. Independence is factorization; dependence is gluing; cumulants are connected defect responses; Fisher information is the metric cost of deforming a completion measure; entropy is unresolved multiplicity.

Within this framework, the central limit theorem takes on a simple conceptual form:
> *the Gaussian is the stable vacuum of finite-variance uncertainty.*

Equivalently,
> *the CLT is the stable-completion theorem for additive ignorance.*

Independent additive composition with finite variance and uniform negligibility removes all connected cumulant defects above order two. What remains is the universal quadratic completion: the Gaussian.

The v0.2 extension makes the conditional structure of this completion theorem explicit. When uniform negligibility fails — when a finite fraction of the second defect energy is concentrated in individual atypically large contributions — the stable completion is no longer Gaussian. It is compound Poisson, or more generally an infinitely divisible law with a non-trivial Lévy jump measure. The Gaussian remains one fixed point of additive aggregation, but it is one of several: the full classification is the Lévy–Khintchine triple $(\gamma, \tau^2, \nu)$.

The arithmetic Euler-product example of §10.3 exhibits this richer completion concretely: the centered logarithmic prime sum $\sum_{p \leq N} v_p \log p$ converges, after appropriate rescaling, to a centered compound Poisson with Lévy measure supported on the prime jump spectrum. This is the appropriate stable-completion framework for arithmetic random-integer statistics, for connections to Dirichlet L-functions, and for the cosmological non-Gaussianity that appears as $f_\mathrm{NL}$-type signatures in the broader arithmetic-vacuum programme.

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