# Expert Panel Discussion: MCMC Convergence in the WHO Tobacco Prevalence Model

## Participants

- **Andrew Gelman** (Columbia University) — Bayesian hierarchical models, Stan, convergence diagnostics
- **Donald Rubin** (Harvard University, Emeritus) — Multiple imputation, causal inference, Bayesian foundations
- **James Berger** (Duke University, Emeritus) — Bayesian decision theory, objective Bayesian analysis, model selection
- **David Dunson** (Duke University) — Bayesian nonparametrics, latent variable models, high-dimensional data

**Moderator**: The model developer

---

## Session 1: Problem Presentation

**Moderator**: Thank you all for joining this discussion. We have a WHO global tobacco prevalence projection model implemented in NIMBLE. The model uses a three-headed stick-breaking construction to jointly estimate cigarettes (CIG), other smoked tobacco (SMKEXTRA), and non-smoked tobacco (ANYEXTRA), enforcing the ordering constraint: P(cig) <= P(any smoked) <= P(any tobacco).

The CIG head has the most complex parameterization — country-level age and cohort splines on top of the hierarchical structure (global → regional → country). The SMKEXTRA and ANYEXTRA heads are simplified to regional spline means only, with country-level intercepts but no country-level splines.

Our convergence issue: the CIG head produces elevated R-hat values, particularly on country-level age and cohort spline coefficients. The SMKEXTRA and ANYEXTRA heads converge well. We run 4 chains, 10K burn-in, 30K iterations, thinning by 5. The model has approximately 200-400 monitored parameters depending on the number of countries with data for each sex.

The specific MCMC settings are:
- Natural cubic splines with 3 interior knots each (age: 25, 45, 65; cohort: 1945, 1965, 1985)
- This gives approximately 4 basis functions per spline
- Each CIG country gets: 1 intercept + 4 age spline coefficients + 4 cohort spline coefficients = 9 country-level parameters
- With ~100-150 countries per sex that have data, that's 900-1350 CIG country-level parameters alone
- Plus regional parameters, global means, variance components, and survey random effects
- dgamma(4,1) hyperpriors on between-region precisions, dgamma(3,1) on within-region precisions
- useConjugacy = FALSE in NIMBLE (defaulting to random walk Metropolis-Hastings samplers)
- Global intercept anchored with tight data-informed prior (sd=0.3)
- Truncated normal priors on age linear effects (forced negative)

**Gelman**: Before we get into solutions, I want to understand the diagnostic picture more precisely. When you say "elevated R-hat," what are the actual numbers? Are we talking 1.1, 1.5, 5.0? And is this concentrated in specific parameter groups, or is it spread across the board? The pattern of non-convergence tells us everything about the cause.

**Moderator**: The worst R-hat values are typically in the 1.2 to 2.0+ range, concentrated in the CIG country-level age spline and cohort spline coefficients — specifically `cig_age_spline[j, l]` and `cig_cohort_spline[j, m]`. The global intercepts, regional intercepts, and the SMKEXTRA/ANYEXTRA parameters generally have R-hat below 1.05. The variance/precision parameters (the dgamma hyperpriors) sometimes show mild elevation around 1.05-1.15. The survey random effects converge well.

**Gelman**: That's the classic signature. The country-level spline coefficients are the problem, and everything else is fine. That tells me immediately this is not a global non-identifiability issue — you've already dealt with that by anchoring the intercepts. This is a local mixing problem: the sampler is getting stuck in the high-dimensional space of correlated spline coefficients within each country.

---

## Session 2: Diagnosis — The Root Cause

**Gelman**: Let me lay out what I think is happening here, and then I want David to weigh in on the spline-specific aspects. You have a hierarchical model where each country's age spline coefficients are drawn from a regional mean:

```
cig_age_spline[j, l] ~ dnorm(cig_age_spline_region_mean[r, l],
                              sd = cig_age_spline_within_region_sd[r])
```

This is the classic "centered parameterization" — you're putting the prior directly on the parameter of interest. The problem is that when the within-region standard deviation is small (strong shrinkage), the posterior for `cig_age_spline[j, l]` is almost entirely determined by the data from country j. But when the standard deviation is large (weak shrinkage), the posterior is dominated by the regional mean. The sampler has to navigate between these two regimes, and random walk Metropolis-Hastings is terrible at this.

This is the exact problem that non-centered parameterization was designed to solve, and it's been the single most impactful reparameterization in the history of Bayesian computation. In Stan, we would write:

```
cig_age_spline_raw[j, l] ~ dnorm(0, 1)
cig_age_spline[j, l] = cig_age_spline_region_mean[r, l] +
                        cig_age_spline_within_region_sd[r] * cig_age_spline_raw[j, l]
```

Now the sampler explores the raw space, which has a much simpler geometry — the raw parameters are a priori independent of the variance component. This breaks the funnel geometry that causes the centered parameterization to mix poorly.

**Dunson**: I agree completely with Andrew's diagnosis, and I want to add a layer. The problem isn't just the centered parameterization — it's compounded by the fact that you have *correlated* spline coefficients within each country. Natural cubic spline basis functions from `ns()` are not orthogonal. The age spline coefficients `cig_age_spline[j, 1:4]` for country j are correlated in the posterior because the basis functions overlap. When basis function 2 increases, basis function 3 might need to decrease to maintain the same fitted curve. This creates ridges in the posterior that random walk samplers follow extremely slowly.

Moreover, you have *both* age and cohort splines for each country, and these interact through the data. Age and cohort are linearly related through period (Cohort = Period - Age), so even though you don't have an explicit period effect, the age and cohort spline coefficients can trade off against each other. Increasing the age effect at age 45 can be partially compensated by adjusting the cohort effect for the relevant birth cohorts. This is the APC identification problem manifesting as a computational problem.

**Berger**: Let me add a perspective on the prior specification. You're using `dgamma(3, 1)` for the within-region precisions. Let me think about what that implies. A Gamma(3, 1) on precision has mean 3 and variance 3, which translates to a prior on the standard deviation that's concentrated around 0.5 to 0.8 or so, with a mode around 0.7. That's actually a fairly informative prior. But the question is whether it's informative in the *right way*.

The issue with gamma priors on precisions — and this has been extensively discussed in the literature, including by Andrew — is that they can create a "precision trap." When the data for a particular country are sparse, the posterior for the precision can push toward very high values (very small variance), which traps the country-level parameters near the regional mean. The sampler then has difficulty moving the country parameters away from the regional mean because the precision is too high, but it can't reduce the precision because the few data points don't provide enough evidence against the prior.

I would recommend considering the half-Cauchy or half-t prior on the standard deviation instead:

```
cig_age_spline_within_region_sd[r] ~ T(dt(0, scale=1, df=3), 0, Inf)
```

This is less informative in the tails and doesn't create the same trapping behavior. It allows the standard deviation to be large when the data support it, which in turn allows the country-level parameters to move more freely.

**Rubin**: I want to step back and think about this from a design perspective. You have ~150 countries, each with 4 age spline coefficients and 4 cohort spline coefficients. That's 1,200 country-level CIG parameters. How many data points does a typical country have? Because the ratio of parameters to data is crucial for understanding why convergence is difficult.

**Moderator**: It varies enormously. Some countries have hundreds of observations from multiple surveys across many years and age groups. Others have as few as 10-20 observations from a single survey. The median is probably around 30-50 observations per country per sex.

**Rubin**: So for countries with 20 observations and 9 country-level CIG parameters (plus contributions to regional and global parameters), the data-to-parameter ratio at the country level is only about 2:1. That's desperately thin. The posterior for those countries is going to be dominated by the prior (the regional mean), and the likelihood surface is going to be nearly flat in many directions. The sampler is essentially doing a random walk in a 9-dimensional space with very little gradient to guide it.

This isn't just a computational problem — it's a modeling problem. You're asking the model to estimate country-specific nonlinear age and cohort patterns for countries where you simply don't have enough data to identify those patterns. The hierarchical structure is supposed to help by borrowing strength, but it can only do so much when the local likelihood is essentially uninformative.

**Gelman**: Don agrees with me on this, I think. This is what I call the "folk theorem of statistical computing": when you have convergence problems, there's almost always a problem with the model, not just the computation. In this case, the model is over-parameterized for data-poor countries. The right solution might not be purely computational — it might involve rethinking the model structure for those countries.

But let me be clear: I'm not saying to simplify the model uniformly. The data-rich countries benefit from having country-level splines. The question is how to handle the continuum from data-rich to data-poor gracefully. The hierarchical structure should handle this in principle — data-poor countries get shrunk toward the regional mean — but the *computation* struggles because the sampler doesn't know how to efficiently navigate between "shrunk to regional mean" and "informed by country data."

---

## Session 3: Solutions — Reparameterization

**Gelman**: Let me lay out the solutions in order of expected impact. The single most important change is the non-centered parameterization. In NIMBLE, you can't do this as cleanly as in Stan because NIMBLE doesn't support the same kind of deterministic transformations in the model definition. But you can approximate it.

The key idea: instead of sampling `cig_age_spline[j, l]` directly from `dnorm(regional_mean, sd)`, you sample a raw parameter and compute the actual parameter as a deterministic function:

```r
# In the NIMBLE model code:
for (j in 1:nCountry) {
  for (l in 1:nAgeSpline) {
    cig_age_spline_raw[j, l] ~ dnorm(0, 1)
    cig_age_spline[j, l] <- cig_age_spline_region_mean[Country_Region[j], l] +
                             cig_age_spline_within_region_sd[Country_Region[j]] *
                             cig_age_spline_raw[j, l]
  }
  for (m in 1:nCohortSpline) {
    cig_cohort_spline_raw[j, m] ~ dnorm(0, 1)
    cig_cohort_spline[j, m] <- cig_cohort_spline_region_mean[Country_Region[j], m] +
                                cig_cohort_spline_within_region_sd[Country_Region[j]] *
                                cig_cohort_spline_raw[j, m]
  }
}
```

Now the sampler explores `cig_age_spline_raw[j, l]`, which has a standard normal prior regardless of the variance component. This eliminates the funnel geometry. The actual `cig_age_spline[j, l]` values are computed deterministically and can still be monitored.

**Dunson**: I want to build on Andrew's point and add something specific to splines. The non-centered parameterization addresses the hierarchical correlation between levels, but it doesn't address the within-country correlation among spline coefficients. For that, we need to think about the basis.

One approach is to use an orthogonalized spline basis. Instead of natural cubic splines from `ns()`, you could use a B-spline basis and then apply a QR decomposition to orthogonalize it. This makes the basis functions uncorrelated in the prior:

```r
# In preprocessing:
B <- ns(ages, knots=c(25,45,65), Boundary.knots=c(15,80))
qr_B <- qr(B)
Q <- qr.Q(qr_B)  # Orthogonal basis
R_mat <- qr.R(qr_B)   # Transformation matrix
# Use Q as the basis in the model, then transform back afterward
```

With an orthogonal basis, the posterior for the coefficients tends to be more spherical, which random walk samplers handle much better.

However, I want to raise a more fundamental point about the spline specification. You're using free spline coefficients with a shared variance component across all basis functions within a country. This means the model treats the smoothness of each basis function contribution equally. A more principled approach would be to use penalized splines (P-splines) where you put a random walk or autoregressive prior on the spline coefficients:

```
cig_age_spline[j, 1] ~ dnorm(0, sd = tau_age[j])
cig_age_spline[j, l] ~ dnorm(cig_age_spline[j, l-1], sd = tau_age[j])  # for l > 1
```

This is equivalent to a first-difference penalty on the spline coefficients and automatically enforces smoothness. The single smoothing parameter `tau_age[j]` (or shared regionally) controls how wiggly the age curve is. This dramatically reduces the effective dimensionality of the spline space and makes MCMC much easier because adjacent coefficients are tied together.

**Berger**: David's point about penalized splines is well taken, but I want to point out that the current model already has substantial regularization through the hierarchical structure — the within-region precision effectively penalizes deviation from the regional mean. The question is whether we need *additional* smoothness constraints on top of that.

From a decision-theoretic perspective, I'd argue yes, because the regularization from the hierarchy is on the *level* of each coefficient (how far from the regional mean) but not on the *differences* between adjacent coefficients (how smooth the curve is). A country could have a very wiggly age curve that happens to have coefficients near the regional mean for each basis function individually, but whose *pattern* is quite different from the regional pattern. The P-spline prior addresses this by penalizing wiggliness directly.

**Gelman**: I want to caution against adding too much complexity to the prior simultaneously. If you do non-centered parameterization AND switch to P-splines AND change the hyperpriors, you won't know which change actually fixed the convergence. My recommendation is to implement changes sequentially:

1. First: non-centered parameterization (biggest expected impact)
2. If that's not enough: block sampling of correlated parameters
3. If still not enough: orthogonalized basis or P-spline prior
4. Last resort: reduce model complexity for data-poor countries

---

## Session 4: Solutions — Sampler Configuration

**Gelman**: Beyond reparameterization, the sampler configuration in NIMBLE matters enormously. You said `useConjugacy = FALSE`, which means NIMBLE defaults to scalar random walk Metropolis-Hastings for everything. That's the worst possible choice for correlated parameters. Each spline coefficient is being proposed independently, but they're highly correlated in the posterior.

NIMBLE has block samplers that propose multiple parameters jointly. For the CIG spline coefficients, you should at minimum block the age spline coefficients together and the cohort spline coefficients together for each country:

```r
mcmc_config <- configureMCMC(nimble_model, useConjugacy = FALSE)

# Remove default scalar samplers for CIG spline coefficients
# Add block samplers instead
for (j in 1:nCountry) {
  # Block all age spline coefficients for country j
  age_nodes <- paste0("cig_age_spline[", j, ", ", 1:nAgeSpline, "]")
  mcmc_config$removeSamplers(age_nodes)
  mcmc_config$addSampler(target = age_nodes,
                         type = "RW_block",
                         control = list(adaptInterval = 200))

  # Block all cohort spline coefficients for country j
  cohort_nodes <- paste0("cig_cohort_spline[", j, ", ", 1:nCohortSpline, "]")
  mcmc_config$removeSamplers(cohort_nodes)
  mcmc_config$addSampler(target = cohort_nodes,
                         type = "RW_block",
                         control = list(adaptInterval = 200))
}
```

The `RW_block` sampler proposes all coefficients in the block simultaneously using a multivariate normal proposal, and it adapts the proposal covariance matrix during the burn-in phase. This is vastly superior to scalar proposals for correlated parameters.

Even better, NIMBLE has the `AF_slice` sampler (automated factor slice sampler) which can handle correlations more robustly:

```r
mcmc_config$addSampler(target = age_nodes,
                       type = "AF_slice",
                       control = list(sliceAdaptFactorInterval = 200))
```

The AF_slice sampler identifies the principal components of the posterior and proposes along those directions. It's particularly good for funnel-shaped posteriors.

**Dunson**: I want to amplify Andrew's point about blocking. But I'd go further — don't just block the age spline coefficients together. Block the age spline coefficients *and* the cohort spline coefficients *and* the country intercept together. All 9 country-level CIG parameters for a country are correlated in the posterior, and they should be proposed jointly:

```r
for (j in 1:nCountry) {
  all_country_nodes <- c(
    paste0("cig_country_intercept[", j, "]"),
    paste0("cig_age_spline[", j, ", ", 1:nAgeSpline, "]"),
    paste0("cig_cohort_spline[", j, ", ", 1:nCohortSpline, "]")
  )
  mcmc_config$removeSamplers(all_country_nodes)
  mcmc_config$addSampler(target = all_country_nodes,
                         type = "AF_slice")
}
```

This creates a single block sampler for all 9 CIG parameters per country. The AF_slice sampler will learn the posterior covariance structure and propose efficient moves. The downside is that each proposal is more expensive (9-dimensional instead of 1-dimensional), but the improved mixing should more than compensate.

**Rubin**: I'm going to push back slightly on the aggressive blocking strategy. With 150 countries, each getting a 9-dimensional block sampler, you now have 150 expensive block proposals per MCMC iteration. Each one requires evaluating the likelihood for all data points in that country (and potentially the full likelihood if NIMBLE doesn't do efficient partial updates). This could make each iteration dramatically slower.

The right metric is effective samples per unit time, not just effective samples per iteration. If blocking makes each iteration 10x slower but only doubles the effective sample size per iteration, you've lost ground.

My recommendation: start with the non-centered parameterization and see if that alone fixes the problem. The non-centered parameterization is free — it doesn't change the computational cost per iteration at all — it just changes the geometry of the posterior. If it's not enough, then add blocking, but monitor the wall-clock time carefully.

**Gelman**: Don makes a good point about wall-clock time. In Stan, we don't have this issue because HMC naturally handles correlations through the gradient. But in NIMBLE with Metropolis-Hastings samplers, there's a real tradeoff between block size and computational cost.

I'd suggest a middle ground: block the age spline coefficients together (4-dimensional block) and the cohort spline coefficients together (4-dimensional block), but keep the country intercept as a separate scalar sampler. The intercept is less correlated with the spline coefficients because it shifts the whole curve up or down, while the spline coefficients change the shape. This gives you two 4-dimensional blocks per country instead of one 9-dimensional block, which is cheaper per proposal.

---

## Session 5: The APC Identification Problem

**Berger**: I want to devote some time to the age-period-cohort identification problem, because I think it's more central to the convergence issue than has been acknowledged so far. The fundamental constraint is:

```
Cohort = Period - Age
```

If you have both age effects and cohort effects in the model, and if the data contain implicit period information (which survey data does — each observation has a survey year), then the age and cohort effects are not jointly identifiable without constraints. Any shift in the age effect can be compensated by an equal and opposite shift in the cohort effect.

In your model, the period effect is absorbed by the survey random effects. But the survey random effects are centered at zero with a shared variance. If a genuine period trend exists, the survey random effects will try to capture it, but the zero-centering constraint means they can't capture a monotone trend perfectly. The leftover trend gets absorbed by some combination of the age and cohort effects, creating a ridge in the posterior.

Now, you've partially addressed this by using natural splines (which have linear constraints at the boundaries) and by anchoring the global intercept. But for data-poor countries, the identification problem is acute because the data don't provide enough information to separate age from cohort effects, and the hierarchical prior is the only thing keeping the parameters from wandering along the ridge.

**Gelman**: Jim is right that the APC identification issue is fundamental, but I want to be precise about how it manifests computationally. The APC identification problem creates a ridge in the posterior — a manifold of roughly equal density along which the age and cohort effects can trade off. The ridge isn't exactly flat (the priors and the curvature of the spline basis break the exact degeneracy), but it's very nearly flat.

Random walk samplers move slowly along ridges because they propose in random directions, most of which are nearly perpendicular to the ridge. Only the rare proposals that happen to align with the ridge direction get accepted. This is why the R-hat is elevated specifically for the spline coefficients — the sampler is slowly diffusing along the ridge.

The non-centered parameterization helps here because it changes the geometry of the ridge. In the non-centered parameterization, the ridge is "straightened out" in the raw parameter space, making it easier for the sampler to move along it.

But there's a more direct fix: explicitly break the identification by constraining the sum of cohort effects to be zero, or by using a contrast parameterization:

```
# Instead of:
cig_cohort_spline[j, m] ~ dnorm(regional_mean, sd)

# Use:
cig_cohort_spline_unconstrained[j, 1:(nCohortSpline-1)] ~ dnorm(regional_mean, sd)
cig_cohort_spline[j, nCohortSpline] <- -sum(cig_cohort_spline_unconstrained[j, 1:(nCohortSpline-1)])
```

This removes one degree of freedom from the cohort spline and breaks the exact linear dependence with the age spline. However, this constraint only makes sense if the spline basis is appropriately centered.

**Dunson**: I think the sum-to-zero constraint is too strong for this application. The birth cohort effect shouldn't necessarily sum to zero — there could be a genuine secular trend where later cohorts have uniformly lower smoking prevalence than earlier cohorts. Forcing a sum-to-zero constraint would bias the cohort effects.

A better approach for the APC identification problem in this specific model is to note that the model already has a partial solution: the linear age effect `cig_age_linear_smooth_effect * age_linear_smooth` captures the dominant linear relationship between age and prevalence. This effectively removes the linear component of the age effect from the spline, leaving only the nonlinear curvature. The APC identification problem is most severe for the linear component (any linear shift in age can be absorbed by a linear shift in cohort), so removing the linear component from the spline already helps enormously.

But I notice that you have the linear age effect as a global parameter, not a country-specific parameter. This means the model assumes the same linear age decline rate across all countries, and the country-specific splines have to capture both the nonlinear age pattern AND any deviation from the global linear trend. For countries where the true linear trend is substantially different from the global average, this could force the country-level splines to absorb a near-linear component, reintroducing the APC identification problem.

**Moderator**: That's an astute observation. The linear age effect is indeed global — `cig_age_linear_smooth_effect ~ T(dnorm(-0.02, sd = 0.05), -Inf, -0.001)` — shared across all countries. The truncation ensures it's negative (prevalence declines at older ages).

**Dunson**: Then consider making it country-specific, or at least regional. If the linear age decline is country-specific, then the country-level splines only need to capture the nonlinear residual, which is both lower-dimensional and less susceptible to the APC identification problem:

```r
for (j in 1:nCountry) {
  cig_age_linear_smooth_effect[j] ~ T(dnorm(cig_age_linear_global, sd = sigma_age_linear), -Inf, -0.001)
}
cig_age_linear_global ~ T(dnorm(-0.02, sd = 0.05), -Inf, -0.001)
sigma_age_linear ~ T(dt(0, scale=0.02, df=3), 0, Inf)
```

This adds ~150 parameters (one per country), but each has a very tight prior centered on the global estimate, and the truncation constraint ensures they stay negative. The country-level splines would then have less work to do, potentially improving convergence.

**Berger**: I'm sympathetic to David's suggestion, but I want to raise a caution. Every parameter you add to the model increases the computational burden and creates potential for new convergence issues. The country-specific linear age effect might converge well because it's a scalar with a strong prior, or it might create new correlations with the country intercept (both shift the overall level of the age curve).

My recommendation would be more conservative: keep the global linear age effect, but add a regional random effect on it:

```r
for (r in 1:nRegion) {
  cig_age_linear_region_effect[r] ~ T(dnorm(cig_age_linear_smooth_effect, sd = 0.01), -Inf, -0.001)
}
```

This allows modest regional variation in the linear age trend without adding 150 new parameters. The key regions where smoking epidemiology differs (e.g., South Asia with rising female smoking vs. Europe with declining) would get their own linear trend, reducing the burden on the country-level splines.

---

## Session 6: Practical NIMBLE Recommendations

**Gelman**: Let me shift to practical implementation. We've discussed several strategies, and I want to prioritize them for the NIMBLE implementation specifically:

**Priority 1: Non-Centered Parameterization**

This is the single highest-impact change. In NIMBLE, you implement it by:

```r
regional_hierarchical_global_ac_model_nimble <- nimbleCode({

  # ... [likelihood unchanged] ...

  # Country-level CIG parameters — NON-CENTERED
  for (j in 1:nCountry) {
    cig_country_intercept_raw[j] ~ dnorm(0, 1)
    cig_country_intercept[j] <- cig_intercept_within_region_sd[Country_Region[j]] *
                                 cig_country_intercept_raw[j]

    for (l in 1:nAgeSpline) {
      cig_age_spline_raw[j, l] ~ dnorm(0, 1)
      cig_age_spline[j, l] <- cig_age_spline_region_mean[Country_Region[j], l] +
                               cig_age_spline_within_region_sd[Country_Region[j]] *
                               cig_age_spline_raw[j, l]
    }

    for (m in 1:nCohortSpline) {
      cig_cohort_spline_raw[j, m] ~ dnorm(0, 1)
      cig_cohort_spline[j, m] <- cig_cohort_spline_region_mean[Country_Region[j], m] +
                                  cig_cohort_spline_within_region_sd[Country_Region[j]] *
                                  cig_cohort_spline_raw[j, m]
    }
  }
})
```

Monitor `cig_age_spline` and `cig_cohort_spline` (the derived parameters), not the raw ones. The sampler explores the raw space, but you can compute the actual spline coefficients for prediction.

In NIMBLE, deterministic nodes (`<-`) are automatically computed, so this should work without changes to the prediction code.

**Priority 2: Block Sampling**

After implementing non-centered parameterization, configure block samplers for the raw parameters:

```r
mcmc_config <- configureMCMC(nimble_model, useConjugacy = FALSE)

for (j in 1:nCountry) {
  # Block the raw age spline coefficients
  age_raw_nodes <- paste0("cig_age_spline_raw[", j, ", ", 1:nAgeSpline, "]")
  mcmc_config$removeSamplers(age_raw_nodes)
  mcmc_config$addSampler(target = age_raw_nodes, type = "RW_block")

  # Block the raw cohort spline coefficients
  cohort_raw_nodes <- paste0("cig_cohort_spline_raw[", j, ", ", 1:nCohortSpline, "]")
  mcmc_config$removeSamplers(cohort_raw_nodes)
  mcmc_config$addSampler(target = cohort_raw_nodes, type = "RW_block")
}
```

Note that you're blocking the *raw* parameters, not the original ones. In the non-centered parameterization, the raw parameters should be less correlated than the original parameters, so the blocks should mix better.

**Priority 3: Variance Component Priors**

If convergence is still unsatisfactory, consider changing the variance component priors from gamma-on-precision to half-Cauchy-on-standard-deviation. In NIMBLE:

```r
# Instead of:
cig_age_spline_within_region_precision[r] ~ dgamma(3, 1)
cig_age_spline_within_region_sd[r] <- 1 / sqrt(cig_age_spline_within_region_precision[r])

# Use:
cig_age_spline_within_region_sd[r] ~ T(dt(0, sigma=1, df=1), 0, Inf)  # Half-Cauchy
```

The half-Cauchy has heavier tails than the implied prior on SD from a Gamma(3,1) on precision, which means it's less likely to trap the variance component at small values. This gives data-poor countries more freedom to have country-specific patterns when the data support it.

Note: NIMBLE's `dt()` distribution uses a slightly different parameterization than base R. Check the NIMBLE manual for the exact syntax. If `dt()` with truncation is problematic in NIMBLE, an alternative is to use a Uniform(0, A) prior on the standard deviation with a generous upper bound A (like 5 or 10). This is the "Gelman-Rubin" recommendation from our 2006 paper with Jennifer Hill.

**Rubin**: I want to emphasize something practical that often gets lost in these discussions: the importance of initial values. You're running 4 chains with dispersed initial values generated by adding small Gaussian noise to fixed starting points. But for the problematic spline coefficients, the initial values might be placing all 4 chains in the same region of the posterior, giving a false sense of convergence.

For the CIG spline coefficients specifically, I would recommend:

1. Initialize each chain with spline coefficients computed from a simple country-specific regression (fit a simple GAM or GLM to each country's data separately, extract the coefficients, and use those as starting values for one chain). This gives at least one chain that starts near the likelihood mode.

2. For the other chains, use overdispersed starting values: set the spline coefficients to the regional means plus/minus large random perturbations. This ensures the chains start from genuinely different locations.

```r
generate_global_inits <- function(seed_offset = 0, emp_mean = -1.5) {
  set.seed(RANDOM_SEED + seed_offset)

  # For chains 2-4: overdispersed starts
  list(
    # ... other inits ...
    cig_age_spline = matrix(rnorm(nCountry * nAgeSpline, 0, 1.5),
                            nrow = nCountry, ncol = nAgeSpline),
    cig_cohort_spline = matrix(rnorm(nCountry * nCohortSpline, 0, 1.5),
                               nrow = nCountry, ncol = nCohortSpline)
  )
}
```

The key insight is that overdispersed starting values are essential for R-hat to be meaningful. If all chains start near the same point, R-hat will be close to 1.0 even if the chains haven't explored the full posterior. My original paper with Andrew (Gelman & Rubin, 1992) was explicit about this: overdispersion is *part of the diagnostic*, not just a convenience.

**Gelman**: I want to strongly second Don's point. And I'll add: our improved R-hat (Vehtari, Gelman, et al., 2021) uses rank normalization and fold-splitting, which makes R-hat more sensitive to within-chain non-stationarity. The improved R-hat also has a recommended threshold of 1.01, not 1.1 or even 1.05. If you're using the old R-hat from `coda::gelman.diag()`, you're getting the 1992 version, which is less sensitive. Consider using the `posterior` package in R which implements the improved R-hat:

```r
library(posterior)
draws <- as_draws_array(mcmc_samples)
summarise_draws(draws, Rhat = rhat, ESS_bulk = ess_bulk, ESS_tail = ess_tail)
```

The `ess_bulk` and `ess_tail` diagnostics from this package are also more informative than the simple ESS from `coda::effectiveSize()`. Bulk ESS measures sampling efficiency for the center of the distribution, while tail ESS measures it for the quantiles. You might find that some parameters have adequate bulk ESS but terrible tail ESS, which means the 95% credible intervals are unreliable even if the posterior mean is well-estimated.

---

## Session 7: Model Simplification and Adaptive Complexity

**Dunson**: I want to revisit a point Don raised earlier about the data-to-parameter ratio. For countries with 20 observations and 9 CIG parameters, the model is fundamentally over-parameterized at the country level. The hierarchical prior regularizes this, but it doesn't eliminate the computational difficulty.

I think a more elegant solution than uniform model complexity is to use an adaptive approach where the model complexity varies by country based on data availability. There are several ways to implement this:

**Option A: Spike-and-Slab on Country Deviations**

Instead of always estimating country-specific spline coefficients, use a spike-and-slab prior that allows some countries to use exactly the regional mean:

```r
for (j in 1:nCountry) {
  for (l in 1:nAgeSpline) {
    # Indicator: does country j need its own age spline coefficient l?
    gamma_age[j, l] ~ dbern(pi_age)
    cig_age_spline_deviation[j, l] ~ dnorm(0, sd = cig_age_spline_within_region_sd[Country_Region[j]])
    cig_age_spline[j, l] <- cig_age_spline_region_mean[Country_Region[j], l] +
                             gamma_age[j, l] * cig_age_spline_deviation[j, l]
  }
}
pi_age ~ dbeta(1, 1)  # Prior on inclusion probability
```

When `gamma_age[j, l] = 0`, the country uses exactly the regional mean. When `gamma_age[j, l] = 1`, the country has its own deviation. The data determines which countries need their own spline coefficients.

However, I should caution that spike-and-slab priors are notoriously difficult for MCMC — the sampler has to jump between the spike (gamma=0) and the slab (gamma=1), which creates its own convergence problems. In NIMBLE, you'd need special samplers for the binary indicators.

**Option B: Regularized Horseshoe Prior on Deviations**

A continuous alternative that achieves similar shrinkage:

```r
for (j in 1:nCountry) {
  for (l in 1:nAgeSpline) {
    lambda_age[j, l] ~ T(dt(0, 1, df=1), 0, Inf)  # Local shrinkage (half-Cauchy)
    cig_age_spline[j, l] ~ dnorm(cig_age_spline_region_mean[Country_Region[j], l],
                                  sd = lambda_age[j, l] * tau_age_global)
  }
}
tau_age_global ~ T(dt(0, 1, df=1), 0, Inf)  # Global shrinkage
```

This is the horseshoe prior adapted for hierarchical spline coefficients. The local shrinkage parameter `lambda_age[j, l]` can be very small (strong shrinkage to regional mean) or very large (no shrinkage), and the data determines which regime each coefficient falls into. The horseshoe prior avoids the discontinuous jump of spike-and-slab while achieving similar effective model selection.

**Berger**: The horseshoe prior is elegant, but I'm concerned about adding two new parameters per spline coefficient (the local shrinkage lambda and its auxiliary). That would roughly double the number of variance parameters in the model. For a model that's already struggling with convergence, adding more parameters seems counterproductive.

I would recommend a simpler approach: use the data availability to determine model complexity *before* fitting the model, not during. This is a pragmatic compromise between full Bayesian model selection and a one-size-fits-all model:

```r
# Before model fitting:
country_data_counts <- clean_data %>%
  group_by(wb_country_abv, sex) %>%
  summarise(n_obs = n(), n_years = n_distinct(year), n_ages = n_distinct(Age_Midpoint))

# Countries with enough data: full splines
# Countries with limited data: intercept only (splines fixed at regional mean)
data_rich_threshold <- 30  # observations
countries_with_splines <- country_data_counts %>%
  filter(n_obs >= data_rich_threshold) %>%
  pull(wb_country_abv)
```

Then in the model, only estimate country-level splines for data-rich countries. Data-poor countries use the regional spline means directly. This is technically a different model, but it's a practical and defensible approach. The threshold can be chosen by cross-validation or sensitivity analysis.

**Gelman**: Jim's pragmatic approach is the one I'd recommend if the goal is to get the model running reliably. The full Bayesian model selection (spike-and-slab, horseshoe) is theoretically beautiful but computationally fraught in MCMC. You'd be trading one convergence problem for another.

The pragmatic approach also has a clean interpretation: "Countries with at least 30 observations have enough data to estimate country-specific age and cohort patterns. Countries with fewer observations use the regional pattern." This is easy to communicate to stakeholders and reviewers.

You could even implement it in the existing model structure without changing the NIMBLE code. Simply fix the within-region standard deviation to zero for data-poor countries (effectively pinning their spline coefficients to the regional means). Or, more cleanly, remove those parameters from the model entirely for data-poor countries.

**Rubin**: I agree with Jim's pragmatic approach, but I want to add one refinement. Instead of a hard threshold (30 observations), use the effective number of parameters as a diagnostic *after* fitting the model. Fit the full model, compute the effective number of parameters for each country (e.g., using WAIC or LOO-CV), and check whether data-poor countries effectively use their 9 parameters. If the effective number of parameters for a country is close to 0 (meaning the posterior is close to the prior), that country doesn't need country-level splines and could be simplified.

This posterior-based approach is more principled than a pre-fitting threshold because it accounts for the information content of the data, not just the sample size. A country with 20 observations that span many ages and years might have more effective parameters than a country with 40 observations from a single survey.

---

## Session 8: NIMBLE-Specific Performance Tuning

**Gelman**: Let me address some NIMBLE-specific performance issues that compound the convergence problem.

**Conjugacy**: You have `useConjugacy = FALSE`. This is correct for the spline coefficients and the non-standard likelihood (the stick-breaking construction makes conjugate updates impossible). However, some parameters in the model might benefit from conjugate updates — specifically, the regional intercepts and regional spline means, which are normal with normal priors in a conditionally conjugate structure. I'd recommend:

```r
mcmc_config <- configureMCMC(nimble_model, useConjugacy = TRUE)
```

Let NIMBLE detect which parameters can use conjugate samplers (these will be exact Gibbs updates, which are always better than Metropolis-Hastings). Then manually override the samplers only for parameters that can't use conjugacy — the CIG spline coefficients, the global intercepts, etc.

**Adaptation**: NIMBLE's default adaptation interval is 200 iterations. For a model with this many parameters, the adaptation might not converge within 200 iterations. Increase it:

```r
mcmc_config$addSampler(target = nodes,
                       type = "RW_block",
                       control = list(adaptInterval = 500,
                                     adaptFactorExponent = 0.6))
```

The `adaptFactorExponent` controls how quickly the adaptation diminishes. The default of 0.8 might adapt too slowly; try 0.6 for more aggressive initial adaptation.

**Burn-in**: Your current burn-in of 10,000 iterations might be insufficient for a model this complex. The adaptation period alone might consume much of the burn-in. I'd recommend at least 20,000 burn-in iterations, and ideally 50,000 for the initial diagnostic run:

```r
NUMBER_OF_BURN <- 50000  # For diagnostic runs
NUMBER_OF_ITERATIONS <- 50000  # Also increase post-burn-in
THINNING_INTERVAL <- 10  # Increase thinning to manage memory
```

The increase in iterations is partly compensated by the improved mixing from non-centered parameterization and block sampling, but for the initial diagnostic run, you want to be generous with iterations.

**nimbleHMC**: NIMBLE now has an HMC (Hamiltonian Monte Carlo) sampler, which is the engine behind Stan's success. For the problematic CIG spline coefficients, HMC could be dramatically more efficient than random walk:

```r
# Replace RW_block with NUTS (the adaptive HMC variant)
mcmc_config$removeSamplers(cig_spline_nodes)
mcmc_config$addSampler(target = cig_spline_nodes,
                       type = "NUTS",
                       control = list(warmupMode = "default"))
```

This is available in recent versions of NIMBLE (via the `nimbleHMC` package). HMC uses gradient information to propose moves along the posterior surface, which is exactly what's needed for the correlated spline coefficients. The gradient naturally follows the ridge in the APC posterior, allowing efficient exploration.

However, HMC requires computing gradients of the log-posterior, which means the model must be differentiable everywhere. The truncated normal priors on the linear age effects (`T(dnorm(...), -Inf, -0.001)`) could be problematic because the truncation boundary creates a discontinuity in the gradient. You might need to replace the truncation with a soft constraint (e.g., a penalty term or a transformed parameterization that ensures negativity).

**Dunson**: Let me add one more NIMBLE-specific point. NIMBLE compiles the model to C++, and the compilation step caches the dependency structure. If you change the sampler configuration but don't recompile, you might get unexpected behavior. Always recompile after changing the MCMC configuration:

```r
compiled_mcmc <- compileNimble(mcmc_built, project = nimble_model, resetFunctions = TRUE)
```

Also, NIMBLE's `resetFunctions = TRUE` flag ensures that the compiled MCMC is rebuilt from scratch, which is important after configuration changes.

---

## Session 9: The Folk Theorem and Model Adequacy

**Gelman**: I want to return to a broader point. We've been discussing computational fixes — reparameterization, block sampling, HMC. These are all important. But the folk theorem of statistical computing says: if you're having trouble with MCMC, there's probably a problem with the model.

In your case, the "problem" might be that the model is trying to estimate more structure than the data can support. The CIG head has country-level age and cohort splines — 8 nonlinear parameters per country — while the SMKEXTRA and ANYEXTRA heads only have country-level intercepts. You made this asymmetry deliberately because the CIG data are richer. But "richer" is relative.

Consider this thought experiment: if you had unlimited data for every country, the CIG spline coefficients would converge perfectly because the data would pin down the coefficients with high precision. The posterior would be concentrated, and any sampler would work. The convergence problem arises because many countries have thin data that can't pin down 8 nonlinear parameters.

So the folk theorem's advice is: simplify the CIG head for data-poor countries. You've already done this for SMKEXTRA and ANYEXTRA — why not apply the same logic within the CIG head? Use country-level splines only where the data justify them, and use regional splines for the rest.

This isn't just a computational convenience — it's the right statistical decision. Estimating 8 country-level nonlinear parameters from 20 observations is asking the model to do something the data can't support, regardless of how good the hierarchical prior is.

**Dunson**: I want to offer a middle ground between full country splines and regional splines only. Instead of removing country-level splines entirely for data-poor countries, reduce the spline complexity. Use fewer knots for data-poor countries:

For data-rich countries (>50 observations): 3 interior knots (4 basis functions) — as currently.
For moderate-data countries (20-50 observations): 1 interior knot (2 basis functions).
For data-poor countries (<20 observations): no knots (linear only) — they use the regional spline pattern.

This creates a three-tier model complexity structure that respects the information content of the data. The challenge is implementing this in NIMBLE, because the spline basis dimension varies by country. One approach is to set up the model with the maximum number of basis functions but fix the extra ones to zero for countries that don't need them:

```r
for (j in 1:nCountry) {
  for (l in 1:nAgeSpline) {
    cig_age_spline[j, l] <- cig_age_spline_full[j, l] * spline_active[j, l]
  }
}
```

where `spline_active[j, l]` is a constant (0 or 1) set before model fitting based on data availability.

**Berger**: Or equivalently, set the within-region standard deviation to a very small value (effectively zero) for the unused spline coefficients of data-poor countries. This keeps the model structure uniform but effectively pins those parameters to their regional means:

```r
# In constants:
# For data-poor country j, set:
cig_age_spline_within_region_sd_effective[j] = 0.001  # Effectively zero
# For data-rich country j:
cig_age_spline_within_region_sd_effective[j] = cig_age_spline_within_region_sd[Country_Region[j]]
```

This is computationally cleaner because the model structure doesn't change — you're just modifying the prior to be extremely informative for data-poor countries. The sampler will quickly converge for those countries because the prior is essentially a point mass.

**Rubin**: I prefer Jim's approach because it's transparent and doesn't require new NIMBLE code. You're just setting different prior precision values for different countries, which is conceptually clean and easy to explain.

But I want to raise a philosophical point: if you simplify the model for data-poor countries, you should be explicit about this in any publication. The model effectively assumes that data-poor countries follow the regional pattern, which is a substantive assumption. In some cases, the countries with the least data might be the ones where smoking patterns are most unusual — precisely because their unusual patterns make survey data harder to collect or because they are low-income countries with different tobacco use patterns.

**Gelman**: Don raises an important concern. The solution is to present results for data-poor countries with appropriate uncertainty. The hierarchical model already does this in principle — data-poor countries have wider credible intervals. But explicitly acknowledging the model assumptions in the paper is essential for scientific integrity.

---

## Session 10: Penalized Complexity Priors

**Berger**: I want to return to the prior specification issue with a specific recommendation. The current priors on the variance components — `dgamma(3, 1)` on precision — have a specific problem that I've been thinking about.

The dgamma(3, 1) prior on precision implies a prior on the standard deviation of approximately:

```
E[sd] ≈ 0.6
mode[sd] ≈ 0.5
P(sd < 0.3) ≈ 0.15
P(sd > 1.0) ≈ 0.10
```

This prior is saying "I expect the between-country variation (on the logit scale) to be about 0.5." Is that reasonable? Let's think about what a logit-scale standard deviation of 0.5 means for prevalence. If the regional mean prevalence is 20% (logit = -1.39), then one standard deviation above the mean is logit(-1.39 + 0.5) = logit(-0.89) ≈ 29%, and one standard deviation below is logit(-1.39 - 0.5) = logit(-1.89) ≈ 13%. So the prior implies countries typically range from about 13% to 29% around a 20% regional mean. That seems plausible for intercepts.

But the same dgamma(3, 1) prior is used for the spline coefficient variances, and spline coefficients have a different interpretation. A spline coefficient of ±0.5 on the logit scale shifts the age-specific prevalence by about ±12 percentage points at the peak of the basis function. For a 4-dimensional spline with unconstrained coefficients, the total variation in the age-prevalence curve could be enormous. The dgamma(3, 1) prior might be too permissive for the spline coefficients.

I'd recommend using the penalized complexity (PC) prior framework from Simpson et al. (2017, Statistical Science). The PC prior is defined by specifying the base model (no country-level deviation) and the probability of exceeding a threshold:

```
# "I believe there's a 5% chance that the country-level deviation
#  exceeds 0.3 on the logit scale for any given spline coefficient"
# This implies:
cig_age_spline_within_region_sd[r] ~ dexp(lambda)
# where lambda = -log(0.05) / 0.3 ≈ 10
```

The exponential prior (which is the PC prior for a standard deviation parameter) has the property of penalizing complexity — it assigns the highest prior density to the simplest model (sd = 0) and decreasing density for more complex models (larger sd). This naturally implements Occam's razor without arbitrary threshold-based model selection.

In NIMBLE:

```r
for (r in 1:nRegion) {
  cig_age_spline_within_region_sd[r] ~ dexp(10)
  cig_age_spline_within_region_precision[r] <- pow(cig_age_spline_within_region_sd[r], -2)
}
```

**Gelman**: I'm not a huge fan of the PC prior for this application because the exponential prior has too much mass near zero — it strongly favors the base model. For a global tobacco model where we know country-level variation exists, I think the half-Cauchy is more appropriate because it allows substantial variation when the data support it while still regularizing toward zero:

```r
cig_age_spline_within_region_sd[r] ~ T(dt(0, scale=0.5, df=1), 0, Inf)
```

The half-Cauchy(0, 0.5) puts more mass on small values than the current dgamma(3,1)-on-precision but has heavier tails. It's a compromise between the strong regularization of the PC prior and the relative permissiveness of the current prior.

But I want to emphasize: the choice of prior on the variance components is less important than the reparameterization and blocking. If you implement non-centered parameterization and block sampling, the model will likely converge well with the current dgamma(3,1) prior. The prior choice affects the *results* but is secondary for *convergence*.

**Dunson**: I agree with Andrew that the prior choice is secondary for convergence. But I want to note that the PC prior has one practical advantage: it shrinks toward zero, which means data-poor countries will have their country-level deviations shrunk more aggressively toward the regional mean. This has the same practical effect as Jim's pragmatic approach of reducing complexity for data-poor countries, but it's data-adaptive and embedded in the model. The model automatically "turns off" country-level splines for countries that don't have enough data to justify them.

---

## Session 11: Diagnostics and Validation Strategy

**Gelman**: Before implementing any of these changes, I want to recommend a diagnostic strategy. Don't just look at R-hat — look at the full diagnostic picture:

1. **Trace plots for worst parameters**: Plot the 4 chains for the 10 worst R-hat parameters. Look for:
   - Chains stuck in different modes (label switching or multimodality)
   - Slow drift (chains moving slowly in one direction — needs more burn-in)
   - One chain behaving differently (initialization issue)

2. **Pairs plots for correlated parameters**: Plot the posterior samples of `cig_age_spline[j, 1]` vs `cig_age_spline[j, 2]`, and `cig_age_spline[j, 1]` vs `cig_cohort_spline[j, 1]` for the worst countries. Look for:
   - Banana-shaped distributions (strong nonlinear correlations — needs HMC or better proposal)
   - Funnel geometry (correlation between parameter and its variance — needs non-centered parameterization)
   - Bimodality (the model is actually multimodal — may need to constrain)

3. **ESS per wall-clock time**: After implementing each change, compute effective samples per second. A change that doubles ESS but triples runtime is a net loss.

4. **Posterior predictive checks**: After convergence, check that the model adequately fits the data. Generate posterior predictions for each country and overlay them on the observed data. If the model fits poorly even with good convergence, the model structure needs to be changed, not just the sampler.

5. **Prior-posterior overlap**: For each variance component, compare the prior and posterior distributions. If they're nearly identical, the data are not informative about that parameter, and you should consider whether it's needed. If the posterior is concentrated near the prior boundary (e.g., near zero for a variance), the model might be over-parameterized.

**Rubin**: I want to add one more diagnostic that I find invaluable: the **chain-split R-hat**. Split each chain in half and compute R-hat using the 8 half-chains instead of the 4 full chains. If the chain-split R-hat is much larger than the full-chain R-hat, it means the chains haven't reached stationarity within each half — there's a slow trend that doesn't show up in the between-chain comparison because all chains are trending in the same direction. This catches the case where all 4 chains start from similar initial values and slowly converge toward the same posterior mode.

Actually, the improved R-hat that Andrew mentioned (Vehtari et al., 2021) already does this split — it computes R-hat using the second half of each chain split into two. So if you switch to the `posterior` package, you get this for free.

**Berger**: From a model selection perspective, I recommend computing the WAIC (widely applicable information criterion) or LOO-CV before and after any model changes. This tells you whether the changes affect the model's predictive performance, not just its convergence. A model that converges beautifully but predicts poorly is worse than a model with some convergence issues that predicts well.

In NIMBLE, you can compute WAIC by monitoring the log-likelihood:

```r
mcmc_config$enableWAIC <- TRUE
# Or monitor node-level log-likelihoods manually
```

Compare WAIC for: (a) the current model, (b) the model with non-centered parameterization, (c) the simplified model for data-poor countries. The WAIC should be similar for (a) and (b) (they're the same model, just reparameterized) and slightly different for (c).

---

## Session 12: Summary of Recommendations

**Gelman**: Let me summarize what I think the implementation roadmap should be, in order of priority:

### Tier 1: Do These First (Biggest Impact)

1. **Non-centered parameterization** for all CIG country-level parameters (intercept + age splines + cohort splines). This is the single most impactful change. Implement by introducing `_raw` parameters with standard normal priors and computing the actual parameters as deterministic transformations.

2. **Block sampling** of CIG spline coefficients. Group the 4 age spline raw coefficients and the 4 cohort spline raw coefficients into blocks. Use `RW_block` or `AF_slice` samplers.

3. **Increase burn-in** to at least 20,000-50,000 iterations for diagnostic runs.

4. **Use improved R-hat** from the `posterior` package with a threshold of 1.01.

### Tier 2: If Tier 1 Isn't Enough

5. **Enable conjugacy** (`useConjugacy = TRUE`) and let NIMBLE use Gibbs updates where possible. Manually configure non-conjugate samplers for the remaining parameters.

6. **Variance component priors**: Switch from dgamma-on-precision to half-Cauchy-on-standard-deviation for the within-region variance components.

7. **Try nimbleHMC/NUTS** for the CIG spline blocks if RW_block mixing is still insufficient.

### Tier 3: Model Structure Changes (If Needed)

8. **Reduce spline complexity for data-poor countries**: Fix country-level spline coefficients to regional means for countries with fewer than ~30 observations. This can be done by setting the within-region SD to a very small value for those countries.

9. **Regional linear age effects**: Make the linear age decline rate regional instead of global, to reduce the burden on country-level splines.

10. **Penalized spline prior**: Add a first-difference penalty on adjacent spline coefficients to enforce smoothness.

**Dunson**: I agree with Andrew's prioritization. The non-centered parameterization should be the first thing you try, and it might be sufficient on its own. If the CIG spline R-hat values drop below 1.05 with just non-centered parameterization, you're done.

One additional recommendation: consider whether 4 chains are necessary. For a model this large, 3 chains might be sufficient for R-hat computation, and the saved compute time could be used for longer chains. The theoretical minimum for R-hat is 2 chains, but 3-4 gives more reliable diagnostics. With the improved R-hat (split chains), even 2 chains give 4 half-chains for the diagnostic.

**Berger**: I want to make one final point about the scientific context. This model is being used for WHO policy projections — the outputs directly inform tobacco control targets for 190+ countries. The convergence of the MCMC is not just a technical issue; it's a question of whether the uncertainty quantification is reliable. If the CIG spline coefficients haven't converged, then the credible intervals for country-level CIG prevalence projections are unreliable, and any policy conclusions based on those intervals are suspect.

Given the stakes, I would recommend running a "gold standard" diagnostic check for the final production run:
- 6 chains (not 4) with genuinely overdispersed starting values
- 100,000 burn-in iterations
- 200,000 post-burn-in iterations
- Thinning by 20
- Improved R-hat < 1.01 for ALL parameters

This will be expensive computationally, but for a WHO publication, the extra runtime is justified by the need for reliable results.

**Rubin**: I agree with Jim. For a model supporting international health policy, you cannot compromise on convergence. The CIG prevalence projections are the headline numbers — they're what journalists and policymakers will cite. If those numbers have unreliable uncertainty bounds, the entire exercise loses credibility.

My final recommendation: after implementing the computational fixes and confirming convergence, do a sensitivity analysis where you fit the model with two different prior specifications (e.g., dgamma(3,1) vs half-Cauchy vs PC prior on the variance components) and compare the country-level CIG prevalence projections. If the projections are robust to the prior choice, you can be confident in the results. If they're sensitive, it means the data are not informative enough to pin down the posterior for some countries, and the projections for those countries should be presented with appropriate caveats.

---

## Session 13: Advanced Approaches and Future Directions

**Dunson**: Before we close, I want to mention two advanced approaches that could be considered for future versions of the model, even if they're too complex for the current implementation:

**1. Gaussian Process Priors for Splines**

Instead of independent spline coefficients with a shared variance, model the age effect as a Gaussian process (GP) with a structured covariance function. The GP prior naturally enforces smoothness and handles the correlation between spline coefficients:

```
cig_age_effect[j, :] ~ GP(regional_mean_function, Matern_kernel(length_scale, amplitude))
```

The Matern kernel has a length-scale parameter that controls smoothness and an amplitude parameter that controls the magnitude of country-level deviations. The GP approach subsumes both the spline basis and the smoothness penalty into a single framework.

In NIMBLE, GP priors can be implemented using the `dmnorm()` distribution with a structured covariance matrix computed from the Matern kernel evaluated at the age knot points. This is a "spline GP" or "predictive process" approximation that combines the computational efficiency of splines with the principled regularization of GPs.

**2. Factor-Analytic Structure for Cross-Country Patterns**

If the country-level spline coefficients are highly correlated across countries (e.g., all countries in a region have similar age-prevalence shapes), you could use a factor-analytic structure to reduce dimensionality:

```
cig_age_spline[j, l] = regional_mean[r, l] + sum_k(loading[l, k] * factor[j, k])
```

where `factor[j, k]` is a K-dimensional latent factor for country j (K << nAgeSpline), and `loading[l, k]` maps the factor to the spline coefficients. This reduces the effective dimensionality from 4 per country to K per country, and with K = 1 or 2, the convergence would be dramatically better.

The factor-analytic structure also has the appealing interpretation that the latent factors represent "tobacco control archetypes" — different patterns of age-specific smoking behavior that countries can load onto to different degrees.

**Gelman**: David's suggestions are interesting for future work, but I want to be pragmatic. For the current model and the current publication timeline, the priority is to get the CIG head converging with the existing model structure. The non-centered parameterization + block sampling should get you there. The GP priors and factor-analytic structure are research directions for the next version of the model.

One more practical suggestion: after implementing the fixes, compare the posterior means from the fixed model to the current model's posterior means. If they're similar (within the credible intervals), then the current results are likely reliable despite the convergence issues — the posterior mean is more robust to mixing problems than the credible interval width. But if they differ substantially, it means the current results are biased by the poor mixing, and the fixed model should be used for the publication.

**Berger**: I'll close with a reminder of the fundamental principle: the goal of MCMC is to compute expectations under the posterior distribution. If you can verify that the posterior expectations of interest (country-level prevalence projections) are stable across independent runs of the MCMC, then the convergence is "good enough" for the purpose at hand, even if individual parameters have elevated R-hat. The question is always: does the convergence issue affect the quantities you actually care about?

In your case, the quantity of interest is the weighted prevalence for each country-year-sex combination, not the individual spline coefficients. The weighted prevalence integrates over the age distribution, which averages out the spline coefficient uncertainty. So the R-hat of the weighted prevalence might be much better than the R-hat of the individual spline coefficients. Check this explicitly:

```r
# In the monitors, add derived quantities:
monitors = c(
  # ... existing monitors ...
  "p_cig"  # Monitor the country-level CIG prevalence for each observation
)
```

Then compute R-hat for the `p_cig` values directly. If the fitted prevalence values have good R-hat even when the underlying spline coefficients don't, you can be more confident that the projection results are reliable.

---

## Appendix A: Quick Implementation Checklist

The following changes are listed in order of implementation priority:

### Step 1: Non-Centered Parameterization (05_models_nimble.R)

```r
# REPLACE the centered country-level CIG parameters:
for (j in 1:nCountry) {
  cig_country_intercept_raw[j] ~ dnorm(0, 1)
  cig_country_intercept[j] <- cig_intercept_within_region_sd[Country_Region[j]] *
                               cig_country_intercept_raw[j]

  for (l in 1:nAgeSpline) {
    cig_age_spline_raw[j, l] ~ dnorm(0, 1)
    cig_age_spline[j, l] <- cig_age_spline_region_mean[Country_Region[j], l] +
                             cig_age_spline_within_region_sd[Country_Region[j]] *
                             cig_age_spline_raw[j, l]
  }

  for (m in 1:nCohortSpline) {
    cig_cohort_spline_raw[j, m] ~ dnorm(0, 1)
    cig_cohort_spline[j, m] <- cig_cohort_spline_region_mean[Country_Region[j], m] +
                                cig_cohort_spline_within_region_sd[Country_Region[j]] *
                                cig_cohort_spline_raw[j, m]
  }
}
```

Monitor both `_raw` (for diagnostics) and the derived parameters (for prediction).

### Step 2: Block Sampling (06_run_global_model.R)

```r
# After configureMCMC:
for (j in 1:nimble_constants$nCountry) {
  # Block raw age spline coefficients
  age_nodes <- paste0("cig_age_spline_raw[", j, ", ", 1:nimble_constants$nAgeSpline, "]")
  mcmc_config$removeSamplers(age_nodes)
  mcmc_config$addSampler(target = age_nodes, type = "RW_block",
                         control = list(adaptInterval = 500))

  # Block raw cohort spline coefficients
  cohort_nodes <- paste0("cig_cohort_spline_raw[", j, ", ", 1:nimble_constants$nCohortSpline, "]")
  mcmc_config$removeSamplers(cohort_nodes)
  mcmc_config$addSampler(target = cohort_nodes, type = "RW_block",
                         control = list(adaptInterval = 500))
}
```

### Step 3: Increase MCMC Iterations (00_config.R)

```r
NUMBER_OF_BURN <- 30000       # Increased from 10000
NUMBER_OF_ITERATIONS <- 60000  # Increased from 30000
THINNING_INTERVAL <- 10       # Increased from 5
```

### Step 4: Use Improved R-hat (mcmc_diagnostics.R)

```r
# Add to diagnostics:
if (requireNamespace("posterior", quietly = TRUE)) {
  draws <- posterior::as_draws_array(mcmc_list)
  improved_diagnostics <- posterior::summarise_draws(draws,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
}
```

---

## Appendix B: Estimated Computational Impact

| Change | Wall-Clock Impact | ESS Improvement | Implementation Effort |
|--------|-------------------|-----------------|----------------------|
| Non-centered param. | +5-10% per iteration | 5-50x for spline params | Moderate (model code change) |
| Block sampling (RW_block) | +20-50% per iteration | 2-10x for spline params | Easy (config change only) |
| Block sampling (AF_slice) | +50-100% per iteration | 5-20x for spline params | Easy (config change only) |
| nimbleHMC/NUTS | +100-300% per iteration | 10-100x for spline params | Moderate (requires nimbleHMC) |
| Increased iterations | Linear increase | Linear increase | Trivial (config change) |
| Half-Cauchy priors | Negligible | Variable (prior-dependent) | Easy (model code change) |
| Data-poor simplification | -10-30% per iteration | Large for simplified countries | Moderate (pre-processing) |
| Orthogonal basis | Negligible | 2-5x for spline params | Easy (pre-processing) |

**Expected scenario**: Non-centered parameterization + block sampling should reduce max R-hat from 1.5-2.0 to below 1.05 for most parameters, with ~30-50% increase in wall-clock time per iteration. Net improvement in ESS/hour: 3-10x.

---

## Appendix C: Key References

1. **Gelman, A. & Rubin, D. B.** (1992). Inference from iterative simulation using multiple sequences. *Statistical Science*, 7(4), 457-472. — The original R-hat diagnostic.

2. **Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C.** (2021). Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667-718. — Improved R-hat with threshold 1.01.

3. **Papaspiliopoulos, O., Roberts, G. O., & Sköld, M.** (2007). A general framework for the parametrization of hierarchical models. *Statistical Science*, 22(1), 59-73. — Theory of centered vs non-centered parameterizations.

4. **Betancourt, M. & Girolami, M.** (2015). Hamiltonian Monte Carlo for hierarchical models. *Current Trends in Bayesian Methodology with Applications*, 79-101. — HMC for hierarchical models, including funnel geometry.

5. **Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sørbye, S. H.** (2017). Penalising model component complexity: A principled, practical approach to constructing priors. *Statistical Science*, 32(1), 1-28. — PC priors.

6. **de Valpine, P., Turek, D., Paciorek, C. J., Anderson-Bergman, C., Lang, D. T., & Bodik, R.** (2017). Programming with models: writing statistical algorithms for general model structures with NIMBLE. *Journal of Computational and Graphical Statistics*, 26(2), 403-413. — NIMBLE reference.

7. **Holford, T. R.** (1991). Understanding the effects of age, period, and cohort on incidence and mortality rates. *Annual Review of Public Health*, 12, 425-457. — APC models and identification.

8. **Lang, S. & Brezger, A.** (2004). Bayesian P-splines. *Journal of Computational and Graphical Statistics*, 13(1), 183-212. — Bayesian penalized splines.

9. **Bhatt, S., et al.** (2015). The effect of malaria control on Plasmodium falciparum in Africa between 2000 and 2015. *Nature*, 526, 207-211. — Example of hierarchical Bayesian geostatistical model with similar convergence challenges.

10. **Gelman, A.** (2006). Prior distributions for variance parameters in hierarchical models. *Bayesian Analysis*, 1(3), 515-534. — Half-Cauchy and other priors for variance components.

---

*Meeting concluded. Transcript prepared for the WHO Tobacco Control Prevalence Projection Model development team.*
