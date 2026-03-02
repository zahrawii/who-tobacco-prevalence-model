# Expert Panel: MCMC Convergence Failure in WHO Tobacco Prevalence Model

## Simulated Academic Discussion

**Panelists:**
- **Andrew Gelman** (Columbia University) — Hierarchical models, Bayesian workflow, Stan
- **Donald Rubin** (Harvard University) — Causal inference, missing data, multiple imputation
- **James Berger** (Duke University) — Bayesian decision theory, objective Bayes, model selection

**Date:** March 2, 2026
**Subject:** Convergence failure in a hierarchical Bayesian APC model for global tobacco prevalence projections (WHO FCTC)

---

## 1. MODEL PRESENTATION

### 1.1 Scientific Context

The WHO Framework Convention on Tobacco Control (FCTC) requires member states to achieve a 30% relative reduction in tobacco prevalence by 2025 (base year 2010). This model projects age-period-cohort (APC) trends for ~190 countries, both sexes, across 6 tobacco indicators.

**Data:** ~25,000 survey observations from national health surveys (GATS, STEPS, DHS, etc.) spanning 15 WHO regions, ~190 countries, ages 15-100, birth cohorts 1900-2010.

### 1.2 Model Specification (NIMBLE)

The model uses a **3-head stick-breaking construction** to ensure probabilistic coherence across nested tobacco product categories:

```
P(cigarettes) = logit^{-1}(mu_cig)
P(any smoked) = P(cig) + logit^{-1}(mu_smkextra) * (1 - P(cig))
P(any tobacco) = P(any smoked) + logit^{-1}(mu_anyextra) * (1 - P(any smoked))
```

This guarantees: P(cigarettes) <= P(any smoked) <= P(any tobacco).

### 1.3 Linear Predictor for Cigarettes (Head 1 — Full Complexity)

```
mu_cig[i] = cig_global_intercept                                           # (1 param)
           + cig_region_intercept[region(country(i))]                       # (15 params)
           + cig_country_intercept[country(i)]                              # (~190 params)
           + cig_def_code_shared * Def_Code_Binary[i]                       # (1 param)
           + w_spline[i] * sum(cig_age_spline[country(i), 1:6] * B_age[i, 1:6])   # (190*6 = 1140 params)
           + w_linear[i] * cig_age_linear_smooth_effect * age_linear[i]     # (1 param)
           + sum(cig_cohort_spline[country(i), 1:6] * B_cohort[i, 1:6])    # (190*6 = 1140 params)
           + sum(cig_age_cohort_interaction[1:36] * B_interaction[i, 1:36]) # (36 params)
           + survey_intercept[survey(i)]                                     # (~200+ params)
```

**Total CIG parameters: ~2,700+**

### 1.4 Linear Predictors for Heads 2 and 3 (Simplified)

```
mu_smkextra[i] = smkextra_global_intercept
               + smkextra_region_intercept[region]
               + smkextra_country_intercept[country]
               + 0.3 * cig_def_code_shared * Def_Code_Binary[i]
               + w_spline[i] * sum(smkextra_age_spline_REGIONAL[region, 1:6] * B_age[i, 1:6])
               + w_linear[i] * smkextra_age_linear_smooth_effect * age_linear[i]
               + sum(smkextra_cohort_spline_REGIONAL[region, 1:6] * B_cohort[i, 1:6])
```

Key difference: Heads 2 and 3 use **regional-level** splines (15*6=90 each) instead of country-level (190*6=1140). No age-cohort interactions.

### 1.5 Hierarchical Structure

```
GLOBAL:     cig_global_intercept ~ N(empirical_mean, sd=0.3)
  |
REGION:     cig_region_intercept[r] ~ N(0, sd=1/sqrt(precision_between))
  |                                     precision_between ~ Gamma(4,1)
COUNTRY:    cig_country_intercept[j] ~ N(0, sd=1/sqrt(precision_within[region(j)]))
                                        precision_within[r] ~ Gamma(4,1)
```

Same hierarchy for splines:
```
GLOBAL:     cig_age_spline_global_mean[l] ~ N(0, sd=1.5)
REGION:     cig_age_spline_region_mean[r,l] ~ N(global_mean[l], sd=1/sqrt(precision_between_age))
COUNTRY:    cig_age_spline[j,l] ~ N(region_mean[region(j),l], sd=1/sqrt(precision_within_age[region(j)]))
```

### 1.6 Key Priors

| Parameter | Prior | Rationale |
|-----------|-------|-----------|
| `cig_global_intercept` | N(empirical_mean, sd=0.3) | Data-informed anchor |
| `smkextra_global_intercept` | N(-3.0, sd=0.3) | ~5% conditional probability |
| `anyextra_global_intercept` | N(-4.0, sd=0.3) | ~2% conditional probability |
| `cig_def_code_shared` | N(0.3, sd=1.0) | Daily vs current definition |
| `residual_sd` | LogNormal(log(0.7), sdlog=0.5) | Observation noise |
| `cig_age_linear_smooth_effect` | TruncNormal(-0.02, sd=0.05, upper=-0.001) | Forced negative age decline |
| `cig_age_spline_global_mean[l]` | N(0, sd=1.5) | Weakly informative |
| `cig_age_cohort_interaction[k]` | N(0, sd=0.2) | Tight — interactions are weak |
| Between-region precisions | Gamma(4,1) | Informative — strong shrinkage |
| Within-region precisions | Gamma(3,1) or Gamma(4,1) | Informative |
| Survey intercepts | N(0, sd=1/sqrt(prec)), prec~Gamma(3,1) | Random effects |

### 1.7 MCMC Configuration

| Setting | Value |
|---------|-------|
| Sampler | NIMBLE default (scalar RW-MH, `useConjugacy=FALSE`) |
| Chains | 4 |
| Burn-in | 5,000 |
| Post-burn-in | 10,000 |
| Thinning | 5 |
| Retained samples/chain | 2,000 |
| Total retained | 8,000 |
| Block samplers | None configured |
| HMC/NUTS | Not used |

### 1.8 Likelihood

The response is logit-prevalence observed with heteroscedastic normal noise:
```
Prevalence[i] ~ N(logit(p_obs[i]), precision = residual_sd^{-2} * weight[i])
```

where `weight[i]` is an observation-level precision weight (e.g., from sample size).

---

## 2. DIAGNOSTICS SUMMARY

### 2.1 Males Model (Feb 27, 2026)

| Metric | Value |
|--------|-------|
| **Grade** | **F (Failed)** |
| Total parameters | 2,542 |
| Median R-hat | 1.061 |
| Max R-hat | **5.67** (`cig_global_intercept`) |
| Parameters R-hat > 1.05 | 1,377 (54.2%) |
| Parameters R-hat > 1.10 | 1,039 (40.9%) |
| Parameters R-hat > 2.00 | 50 (2.0%) |
| Mean ESS | 874 |
| Min ESS | **8** (`cig_global_intercept`) |
| Parameters ESS < 100 | 479 (18.8%) |

**By head type:**
| Head | N params | Median R-hat | % > 1.1 |
|------|----------|-------------|---------|
| CIG | 1,872 | 1.110 | 52.6% |
| SMKEXTRA | 335 | 1.013 | 9.3% |
| ANYEXTRA | 335 | 1.003 | 7.2% |

**Worst parameters (males):**

| Parameter | R-hat | ESS | Chain means |
|-----------|-------|-----|-------------|
| `cig_global_intercept` | 5.672 | 8 | [-0.66, -0.17, -0.75, -0.66] |
| `cig_region_intercept[13]` | 5.602 | 10 | [0.20, 0.60, 1.26, 1.35] |
| `cig_cohort_spline_region_mean[13,3]` | 4.530 | 30 | [-5.85, -7.54, -7.67, -7.98] |
| `cig_region_intercept[10]` | 4.074 | 19 | [0.97, 1.40, 1.86, 2.09] |
| `cig_region_intercept[3]` | 3.938 | 21 | [0.15, 0.23, 0.63, 0.70] |

**All 15 `cig_region_intercept` parameters:** 100% above R-hat 1.1, 73% above 2.0, 100% ESS below 100.

**Multimodality:** 3 modes detected for `cig_global_intercept`, 3-5 modes for regional intercepts.

**Key correlations (males):**
| Pair | Correlation |
|------|-------------|
| `cig_global_intercept` ↔ `cig_age_linear_smooth_effect` | -0.475 |
| `cig_global_intercept` ↔ `cig_region_intercept[3]` | -0.513 |
| `cig_region_intercept[13]` ↔ `cig_cohort_spline_region_mean[13,3]` | **-0.841** |
| `cig_region_intercept[10]` ↔ `cig_region_intercept[13]` | **0.843** |
| `cig_age_linear_smooth_effect` ↔ `cig_cohort_spline_region_mean[13,3]` | **0.697** |

### 2.2 Females Model (Feb 28, 2026)

| Metric | Value |
|--------|-------|
| **Grade** | **F (Failed)** |
| Total parameters | 2,553 |
| Median R-hat | 1.074 |
| Max R-hat | **6.77** (`cig_global_intercept`) |
| Parameters R-hat > 1.05 | 1,479 (57.9%) |
| Parameters R-hat > 1.10 | 1,134 (44.4%) |
| Parameters R-hat > 2.00 | 25 (1.0%) |
| Mean ESS | 829 |
| Min ESS | **9** (`cig_global_intercept`) |
| Parameters ESS < 100 | 458 (17.9%) |

**By head type:**
| Head | N params | Median R-hat | % > 1.1 |
|------|----------|-------------|---------|
| CIG | 1,881 | 1.116 | 53.7% |
| SMKEXTRA | 336 | 1.019 | 19.3% |
| ANYEXTRA | 336 | 1.006 | 17.3% |

**Worst parameters (females):**

| Parameter | R-hat | ESS | Chain means |
|-----------|-------|-----|-------------|
| `cig_global_intercept` | 6.774 | 9 | [-3.62, -2.75, -2.75, -2.58] |
| `cig_region_intercept[11]` | 4.340 | 19 | [0.51, -0.05, 0.40, 0.09] |
| `cig_region_intercept[7]` | 3.719 | 34 | [0.50, 0.08, 0.09, -0.54] |
| `cig_region_intercept[6]` | 3.247 | 19 | [0.71, -0.18, 0.83, 0.52] |
| `cig_region_intercept[9]` | 2.806 | 22 | [-0.46, -0.91, -0.76, -1.26] |

**All 15 `cig_region_intercept` parameters:** 100% above R-hat 1.1, 80% above 2.0, 100% ESS below 100.

**Key correlations (females):**
| Pair | Correlation |
|------|-------------|
| `cig_global_intercept` ↔ `cig_age_linear_smooth_effect` | **-0.783** |
| `cig_global_intercept` ↔ `cig_region_intercept[7]` | **-0.725** |
| `cig_global_intercept` ↔ `cig_region_intercept[9]` | -0.577 |
| `cig_region_intercept[7]` ↔ `cig_region_intercept[9]` | **0.729** |

---

## 3. EXPERT DISCUSSION

---

### ANDREW GELMAN (Columbia University)

**Opening Statement:**

Let me be direct: this model has a textbook non-identifiability problem, and R-hat values of 5-7 are not just "bad convergence" — they mean **the posterior is multimodal and the chains are not mixing**. An R-hat above 1.1 is concerning. An R-hat of 5.67 means the chains are sampling from effectively different distributions. You do not have a posterior; you have four separate point estimates from four separate modes.

**Diagnosis:**

The core issue is what I'd call a **hierarchical intercept funnel**. You have:

```
mu = global + region[r] + country[j] + splines + ...
```

The prior on `cig_global_intercept` is N(empirical_mean, sd=0.3). That sounds tight, but it's not tight enough when you have 15 regional intercepts that can absorb any shift. The likelihood only constrains the *sum* (global + region + country) — it cannot separately identify the three levels. Your sd=0.3 prior creates a soft constraint, but the chains find different "compromise" points between the global level and the regional deviations.

Look at the male model: Chain 2 has `cig_global_intercept` at -0.17 while the other chains are at -0.66 to -0.75. That's a 0.5 unit gap on the logit scale — which is *huge* for prevalence. Chain 2 has compensated by shifting all regional intercepts down. The likelihood is nearly identical; the chains are exploring different ridges of an elongated posterior.

The -0.841 correlation between `cig_region_intercept[13]` and `cig_cohort_spline_region_mean[13,3]` is even more damning. This means the model can't distinguish "region 13 has high baseline prevalence" from "region 13 has a particular cohort trend." These two explanations are observationally equivalent over the data range.

**The females model is worse** (R-hat=6.77 vs 5.67) because female smoking rates are lower in many regions, so there's even less data to identify the CIG parameters. Chain 1 is stuck at -3.62 while the others are near -2.6 to -2.75 — a full unit apart.

**My recommendations, in priority order:**

**1. Reparameterize to non-centered form.** This is the single most impactful change for hierarchical models in MCMC. Instead of:
```
region[r] ~ N(0, sigma_region)
```
Write:
```
region_raw[r] ~ N(0, 1)
region[r] = sigma_region * region_raw[r]
```
This decouples the variance parameter from the location parameters and eliminates the funnel geometry. I've written about this extensively — it's why Stan's non-centered parameterization works so much better than BUGS/JAGS for hierarchical models. NIMBLE supports this; you just have to write it manually.

**2. Use HMC, not random-walk Metropolis.** Your model has 2,500+ parameters with massive correlations. Scalar RW-MH is the worst possible sampler for this geometry. It proposes one parameter at a time and can never traverse the correlated ridges efficiently. HMC uses gradient information to propose moves *along* the ridges. NIMBLE has `nimbleHMC` — use it. Or better yet, rewrite in Stan. I know that's not what you want to hear, but Stan's NUTS sampler was designed precisely for this class of problem.

**3. If staying with RW-MH, you need block samplers.** At minimum, block-sample `(cig_global_intercept, cig_region_intercept[1:15])` together. Also block the regional cohort spline means with their intercepts — the -0.841 correlation demands it.

**4. The MCMC budget is laughably insufficient.** 2,000 post-burnin samples per chain for 2,500 parameters with multimodality? You need at minimum 10,000-20,000 *effective* samples for reliable inference. Given your current ESS of 8 for the worst parameters, you'd need roughly 250x more iterations to get ESS=2000 — which means millions of iterations. This is another argument for HMC, which would give you ESS close to the nominal sample size.

**5. Re-enable conjugacy.** You set `useConjugacy=FALSE`, which is actively harmful. Many of your precision parameters (the Gamma priors) have conjugate posteriors. Gibbs sampling those exactly is always better than approximating with RW-MH. There is no good reason to disable this.

**On the prior for `cig_global_intercept`:** sd=0.3 on the logit scale spans roughly a 15 percentage point range in prevalence (depending on the baseline). That's wide enough for the global mean to wander significantly while regions compensate. But I would *not* make this tighter — that would just mask the identifiability problem by making the prior dominate. Fix the parameterization instead.

**On the country-level splines:** You have 190 countries * 6 cohort basis functions = 1,140 cohort spline parameters for CIG alone. The hierarchical prior (through regional means) is supposed to regularize these, but the regional means themselves aren't converging. This is a cascading failure. I'd consider whether country-level cohort splines are justified epidemiologically — do you really believe that birth cohort effects differ by country within a region? If the evidence is weak (and the non-convergence suggests it is), move to regional cohort splines for CIG as well, matching the SMKEXTRA/ANYEXTRA specification.

---

### DONALD RUBIN (Harvard University)

**Opening Statement:**

I want to approach this differently from Andrew. Before we dive into MCMC mechanics, let me ask a more fundamental question: **Is this the right model for the scientific question?**

You're trying to project tobacco prevalence forward to 2025 and 2040. The model has ~2,500 parameters fit to ~25,000 observations — roughly 1 parameter per 10 observations. For a projection model, that's alarming. Overfitting is your enemy for prediction, and the convergence failure may actually be telling you something important about model complexity.

**On the missing data perspective:**

Each country-sex-age-cohort-year cell should, in principle, have an observation. Most don't. What you're really doing is a massive imputation exercise — filling in the unobserved cells using the hierarchical structure. The question is whether your model has the right structure for that imputation.

**The stick-breaking construction concerns me.** You're modeling:
```
P(cig) = logit^{-1}(mu_cig)
P(any smoked | not cig) = logit^{-1}(mu_smkextra)
P(any tobacco | not smoked) = logit^{-1}(mu_anyextra)
```

This creates an ordering dependency: the CIG head must be identified first, and errors propagate through. If `mu_cig` is poorly identified (which it clearly is), then even if `mu_smkextra` is well-identified *conditionally*, the marginal inference on P(any smoked) inherits the non-identifiability. Have you verified that the posterior predictive distributions for P(any smoked) and P(any tobacco) are well-behaved despite the CIG convergence failure?

**My primary concern is the causal structure.** Your linear predictor mixes age, period (through the intercept trend), and cohort effects. This is the classic APC identification problem — you cannot separately identify age, period, and cohort effects without constraints because `cohort = period - age`. Your constraints come implicitly from the spline basis functions and the hierarchical structure, but the non-convergence suggests these constraints are insufficient.

**The 0.3 scaling factor on `cig_def_code_shared` for SMKEXTRA and ANYEXTRA** — where does this come from? This is a strong assumption (the definition effect is 30% as large for other products) that is not estimated from data. If this is wrong, it could create systematic bias in the stick-breaking decomposition.

**My recommendations:**

**1. Conduct a posterior predictive check, even with the current non-converged samples.** If the *predictions* are stable across chains even though individual parameters aren't, you have a benign non-identifiability. The scientific question is about prevalence projections, not about the value of `cig_global_intercept`. Check whether the chains agree on P(cig) for key country-age-year combinations. If they do, the non-identifiability may not matter for your policy conclusions.

**2. If predictions differ across chains, reduce model complexity.** Start with the simplest model that could answer the scientific question — perhaps regional-level parameters for everything — and add complexity only where the data demand it. This is the principle of parsimony for prediction. Your current approach of starting with maximum complexity and trying to regularize it with priors is backwards for a projection exercise.

**3. Consider a two-stage approach.** First fit a simpler model that converges (perhaps with regional splines for all heads). Then use those results as fixed priors for a second stage with country-level refinements. This is essentially what your country-specific model does, but the first stage needs to converge properly.

**4. On the APC identifiability:** Consider whether you need both age splines AND cohort splines AND an age-cohort interaction for CIG. The 36 interaction parameters (6*6 tensor) are particularly suspicious — they give the model enormous flexibility to trade off between age and cohort effects. The N(0, sd=0.2) prior constrains the magnitude, but doesn't resolve the identifiability between the interaction terms and the main effects. I would start by removing the interaction and seeing if convergence improves.

**5. The `weight[i]` in the likelihood:** This heteroscedastic weighting is doing a lot of work. If some observations have much higher weight than others, they effectively anchor the parameters, but only for specific countries/ages. This could create "information gaps" where parameters are poorly identified. Examine whether the non-converging parameters correspond to regions/countries with low total weight.

---

### JAMES BERGER (Duke University)

**Opening Statement:**

Both Andrew and Donald raise valid points, but I want to focus on what I see as the deepest issue: **the prior specification reveals a fundamental tension between identifiability and scientific coherence.**

Let me explain. You've placed a N(empirical_mean, sd=0.3) prior on `cig_global_intercept`. This is what I would call a "pragmatic identifiability prior" — its purpose is computational (anchor the intercept) rather than scientific (express genuine prior knowledge). The sd=0.3 was presumably chosen because tighter priors felt too informative and looser ones didn't identify the model.

This is a red flag. When you need a specific prior width to achieve identification, you're in a region where the likelihood is flat or multimodal. The prior is doing the work that the model structure should be doing. I'd prefer to see identifiability achieved through *structural constraints* rather than informative priors on individual parameters.

**Formal identifiability analysis:**

Your linear predictor for CIG has:
```
mu = alpha_0 + alpha_r + alpha_j + f_age(x) + f_cohort(x) + g(age,cohort) + survey_s
```

where `alpha_0` is global, `alpha_r` is regional, `alpha_j` is country, `f_age` is a spline in age, `f_cohort` is a spline in birth cohort, `g` is the interaction, and `survey_s` is a random effect.

The non-identifiabilities are:

1. **Level non-identifiability:** `alpha_0 + c` and `alpha_r - c` give the same likelihood for any constant `c`. This is the classic random effects centering issue. Your prior breaks this softly, but the MCMC evidence shows it's not sufficient.

2. **Trend non-identifiability:** The constant basis function in the cohort spline overlaps with the intercept hierarchy. If the cohort spline basis includes an intercept-like component (e.g., B-splines sum to 1), then the cohort spline's "level" is confounded with the country intercept. This explains the -0.841 correlation between `cig_region_intercept[13]` and `cig_cohort_spline_region_mean[13,3]`.

3. **APC confounding:** As Donald noted, age + cohort = period. The interaction terms create additional confounding between main effects and interactions.

**My recommendations:**

**1. Impose sum-to-zero constraints on the regional intercepts.** This is the standard solution in the random effects literature and it works. Define:
```
cig_region_intercept[r] ~ N(0, sigma) for r = 1, ..., 14
cig_region_intercept[15] = -sum(cig_region_intercept[1:14])
```
This eliminates non-identifiability #1 *structurally* rather than relying on the prior. It's cleaner than tightening the prior because it doesn't inject artificial information. The global intercept then has a unique interpretation as the unweighted mean across regions.

Critically, apply the same constraint to country intercepts within each region:
```
cig_country_intercept[j] such that sum over countries in region r = 0
```

This makes every level of the hierarchy identifiable by construction.

**2. Orthogonalize the spline basis against the intercepts.** If your B-spline basis includes a constant component (or nearly constant component), the cohort spline "level" is confounded with the intercept. Center each spline basis function to have mean zero across the observed covariate values:
```
B_centered[i,l] = B[i,l] - mean(B[,l])
```
This ensures the spline captures *shape* variations only, while the intercept captures the *level*. This directly addresses non-identifiability #2 and should dramatically reduce the intercept-spline correlations.

**3. Use a reference prior analysis to calibrate your priors.** Rather than guessing sd=0.3, compute the reference prior (Jeffreys-style) for `cig_global_intercept` conditional on the other parameters at their MAP estimates. This tells you the maximum information the data can provide about this parameter. If the reference prior is very flat, it confirms the non-identifiability is fundamental and no proper prior of reasonable width will resolve it — you need structural constraints.

**4. On the Gamma(4,1) precision priors:** These imply a prior mean of 4 for the precision, which corresponds to sd = 0.5. The prior puts most mass on sd between 0.3 and 1.0. This is reasonable for the intercept hierarchy, but may be too constraining for the spline precisions. If the true within-region variation in cohort splines is larger than sd=1, your Gamma(3,1) prior is pulling the precision too high (variance too low), which forces the spline coefficients to be similar across countries within a region. When the data disagree, this creates tension that manifests as non-convergence.

Consider a half-Cauchy(0, 1) prior on the standard deviations instead. This has the heavy tail needed to accommodate unexpected variation while still regularizing. Andrew has written extensively about this — I believe his 2006 paper with the "folk theorem" about hierarchical variance parameters is directly relevant here.

**5. On the `useConjugacy=FALSE` setting:** I agree with Andrew that this is harmful, but for a deeper reason. When you disable conjugacy, you're replacing *exact* conditional distributions with *approximate* ones. For the precision parameters, the Gibbs update from the conjugate Gamma posterior is not just faster — it's *exact*. RW-MH proposals for these parameters will have lower acceptance rates and poorer mixing, which propagates to the parameters that depend on them (the intercepts and splines). Re-enabling conjugacy won't fix the fundamental identifiability issue, but it removes an unnecessary source of inefficiency.

---

## 4. DEBATE AND CROSS-EXAMINATION

---

### Gelman responds to Rubin:

**On posterior predictive checks:** Donald, I strongly disagree with the suggestion that if predictions are stable, the non-identifiability "may not matter." Yes, for point predictions it might be benign. But for **uncertainty quantification** — which is the whole point of a Bayesian model for policy — non-identifiability inflates posterior variance in unpredictable ways. The probability of achieving the 30% reduction target depends on the tails of the predictive distribution, not just the mean. If the chains disagree about parameter values, they'll disagree about uncertainty, which means your "probability of achieving target" numbers are unreliable.

Furthermore, with R-hat > 5, you don't even have a valid posterior to check. The chains haven't converged, so any posterior predictive check is comparing four different models, not four draws from one model.

**On model complexity:** I agree the model is complex, but simplifying to regional-level everything would lose the whole point — country-specific projections. The solution isn't to simplify the model; it's to parameterize it in a way that the sampler can explore. Non-centered parameterization + HMC can handle thousands of parameters routinely — I've seen Stan models with 50,000+ parameters converge beautifully with proper parameterization.

### Rubin responds to Gelman:

**On prediction vs. parameters:** Andrew, you're right about uncertainty quantification in principle. But let me push back: if the model is not identified for *parameters*, what does the uncertainty interval on prevalence projections mean? It's partly reflecting genuine uncertainty and partly reflecting the artificial variance from non-identifiability. These are mixed together and you can't separate them. A simpler model that converges would give you honest uncertainty intervals, even if wider, because they'd reflect only genuine uncertainty.

**On non-centered parameterization:** You and I have a deep disagreement about whether reparameterization is "fixing" a problem or "hiding" it. The non-centered parameterization works beautifully when the hierarchy is well-specified. But if the fundamental issue is that country-level cohort splines are not supported by the data, reparameterization will just give you faster convergence to a posterior that's dominated by the prior. You'll get good R-hat values and feel good about it, but the scientific content of the inference won't improve.

The question I'd want answered is: **for how many countries do you have enough data to estimate country-specific cohort trends?** If it's only 20 out of 190, then for the other 170, the "country-level" cohort spline is just the regional mean plus noise from the prior. That's not scientific; it's prior hallucination with extra steps.

### Berger responds to both:

**On reparameterization vs. structural constraints:** Andrew, I appreciate the pragmatic value of non-centered parameterization, but I side partially with Donald here. Non-centered form helps the *sampler* but doesn't resolve the *statistical* non-identifiability. If you achieve R-hat < 1.01 with non-centered HMC, you've converged to a posterior — but is it the right posterior?

Consider: with non-centered parameterization and HMC, the chains will explore the full ridge of the posterior efficiently. But the *marginal* posterior for `cig_global_intercept` will still be wide and prior-dominated, and the regional intercepts will show the compensating pattern. You'll have converged to a posterior that says "we're uncertain about whether the global mean is -0.2 or -0.7, but we know the sum is constrained." That's valid *given the model*, but the wide marginals might concern reviewers and policymakers.

With sum-to-zero constraints, the global intercept has a unique, interpretable value, and the marginal posteriors are all well-behaved. This is not just a computational convenience — it produces cleaner scientific communication.

**On the 36 age-cohort interaction terms:** I notice that nobody has defended these yet, and I think they should be the first thing to go. The interaction has 36 parameters (6 age * 6 cohort basis functions), each with an N(0, 0.2) prior. But these are a rank-36 tensor product. In practice, most of the interaction is low-rank — perhaps rank 2-3. You're fitting 36 parameters where 6-9 would suffice. This doesn't just waste parameters; it creates 36 dimensions of near-non-identifiability. Use a low-rank tensor decomposition or simply remove the interaction term. Your SMKEXTRA and ANYEXTRA heads don't have it and they converge fine.

### Gelman's rebuttal:

**On structural constraints:** Jim, I agree that sum-to-zero is cleaner for communication. But in my experience, it's harder to implement in NIMBLE than non-centered parameterization, and it can create its own problems (the "derived" 15th intercept has a different posterior geometry). In Stan, I'd use a soft sum-to-zero with a tight prior: `sum(region_intercept) ~ N(0, 0.001)`. In NIMBLE, you'd need a custom constraint.

But here's where I'll concede to both of you: **the right first step is to remove the age-cohort interaction AND move CIG cohort splines to regional level, matching the other heads.** This cuts ~1,200 parameters and removes the worst non-identifiabilities. Then reparameterize, then add HMC. If after that, you still see non-convergence, we can discuss structural constraints.

The key point is: **don't try to fix everything at once.** Fix the biggest problems first, refit, diagnose, and iterate. This is the Bayesian workflow.

### Rubin's closing:

I want to add one practical point. This model feeds into WHO policy decisions. The probability-of-achieving-target calculations will be used to classify countries as "on track" or "not on track." A country classified as 60% likely to achieve the target vs. 40% likely will receive different policy attention. If your posterior is multimodal, the probability calculation is essentially averaging over different models, which is not what the policymakers think they're getting.

**My final recommendation:** Before touching the MCMC, run the model on a simplified subset — perhaps 30 countries in 5 regions — and verify that you can achieve convergence. If you can, scale up carefully. If you can't, the model is fundamentally mis-specified for this problem.

### Berger's closing:

I want to address the elephant in the room: **the observation model.** You're fitting logit-prevalence with a normal likelihood:
```
logit(prevalence) ~ N(logit(p_obs), precision)
```

But `p_obs` itself is a nonlinear function of the parameters through the stick-breaking construction. This means `logit(p_obs)` is a nonlinear transform of multiple logistic functions, and the Jacobian of this transform varies across the parameter space. The normal likelihood on the logit scale is an approximation to what should be a binomial likelihood on the raw scale. This approximation may create artificial multimodality — the nonlinear transform maps multiple parameter combinations to similar `logit(p_obs)` values, creating ridges in parameter space that wouldn't exist with a proper binomial likelihood.

**If feasible, switch to a binomial likelihood** with the stick-breaking probabilities directly:
```
Y[i] ~ Binomial(n[i], p_obs[i])
```

This would remove one layer of nonlinearity from the model and may significantly improve identifiability.

---

## 5. CONSENSUS RECOMMENDATIONS

Despite disagreements on approach, the three panelists converge on the following ordered recommendations:

### Tier 1: Immediate (do before the next model run)

| # | Recommendation | Gelman | Rubin | Berger | Impact |
|---|---------------|--------|-------|--------|--------|
| 1 | Re-enable conjugacy (`useConjugacy=TRUE`) | Strong agree | Agree | Strong agree | Low effort, immediate ESS gains |
| 2 | Remove the 36 age-cohort interaction terms for CIG | Agree | Strong agree | Strong agree | Cuts 36 params, reduces confounding |
| 3 | Increase MCMC budget to 50k+ iterations post-burnin | Strong agree | Agree | Agree | Necessary regardless of other fixes |

### Tier 2: Structural (require model code changes)

| # | Recommendation | Gelman | Rubin | Berger | Impact |
|---|---------------|--------|-------|--------|--------|
| 4 | Move CIG cohort splines to regional level | Agree | Strong agree | Agree | Cuts ~1,100 params |
| 5 | Block-sample correlated intercept groups | Strong agree | Neutral | Agree | Breaks ridge geometry |
| 6 | Sum-to-zero constraint on regional intercepts | Cautious | Agree | Strong agree | Structural identifiability |
| 7 | Orthogonalize spline basis against intercepts | Neutral | Neutral | Strong agree | Removes level confounding |
| 8 | Non-centered parameterization | Strong agree | Cautious | Agree | Eliminates funnel geometry |

### Tier 3: Major Changes (consider if Tiers 1-2 insufficient)

| # | Recommendation | Gelman | Rubin | Berger | Impact |
|---|---------------|--------|-------|--------|--------|
| 9 | Switch to HMC/NUTS (nimbleHMC or Stan) | Strong agree | Neutral | Agree | Fundamental sampler improvement |
| 10 | Replace normal-logit with binomial likelihood | Neutral | Agree | Strong agree | More faithful observation model |
| 11 | Reduce to regional-only model first, add complexity | Disagree | Strong agree | Neutral | Simplification for projection |
| 12 | Run on subset (30 countries, 5 regions) to validate | Neutral | Strong agree | Agree | Diagnostic step |

### Tier 4: Diagnostic Steps (parallel to fixes)

| # | Recommendation | All agree |
|---|---------------|-----------|
| 13 | Posterior predictive check: do chains agree on *predictions* even if they disagree on *parameters*? | Yes |
| 14 | Count how many countries have enough data for country-level cohort estimation | Yes |
| 15 | Check whether non-converging parameters correspond to data-sparse regions | Yes |
| 16 | Examine whether the observation `weight[i]` values create information gaps | Yes |

---

## 6. IMPLEMENTATION PRIORITY

If the goal is to fix convergence with minimal model restructuring, the panelists suggest this sequence:

```
Step 1: Re-enable conjugacy + remove age-cohort interaction + increase iterations
         → Refit and check diagnostics
         → Expected: moderate improvement, interaction parameters converge

Step 2: Move CIG cohort splines to regional level + block samplers for intercepts
         → Refit and check diagnostics
         → Expected: substantial improvement, ~1100 fewer parameters

Step 3: Sum-to-zero constraint OR non-centered parameterization
         (choose one based on preference)
         → Refit and check diagnostics
         → Expected: resolves remaining non-identifiability

Step 4: If still not converging, consider HMC or Stan rewrite
```

**Estimated convergence after each step (expert prediction):**
- After Step 1: CIG median R-hat drops from ~1.1 to ~1.06, global intercept R-hat from ~6 to ~3
- After Step 2: CIG median R-hat drops to ~1.02, global intercept R-hat to ~1.5
- After Step 3: All parameters R-hat < 1.05, global intercept R-hat < 1.02
- After Step 4: All parameters R-hat < 1.01 (if needed)

---

## 7. KEY DISAGREEMENTS (Unresolved)

1. **Gelman vs. Rubin on model complexity:** Gelman believes proper parameterization can support 2,500+ parameters; Rubin believes the data cannot support this complexity for projection and advocates simplification. This is a philosophical difference about the role of hierarchical priors: are they scientific (representing real between-country variation) or regularization (preventing overfitting)?

2. **Gelman vs. Berger on structural constraints:** Gelman prefers soft computational solutions (reparameterization, better samplers) that let the posterior "be what it is." Berger prefers hard structural constraints (sum-to-zero) that eliminate non-identifiability by construction. Both are valid; the choice depends on whether you want the model to express uncertainty about the global mean (Gelman) or want the global mean to be uniquely defined (Berger).

3. **Rubin on prediction vs. inference:** Rubin's point about posterior predictive stability is important but unresolved — nobody has checked whether the predictions are actually stable across chains. This should be the very first diagnostic step.

---

*Note: This document presents a simulated academic discussion based on the published methodological positions and known expertise of the named researchers. It does not represent their actual opinions on this specific model.*
