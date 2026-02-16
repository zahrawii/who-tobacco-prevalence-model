# WHO Tobacco Control Prevalence Projection Model - NIMBLE v4 Documentation

## Document Metadata

**Document Version**: 2.3.2
**Code Version**: 2.3.2 (NIMBLE with vectorized prediction optimizations)
**Last Updated**: 2026-02-13
**Developed for**: World Health Organization

---

## Quick Reference Card

### Critical Equations

```
# Stick-Breaking Construction
P(cigarettes)  = logistic(μ_cig)
P(any smoked)  = P(cig) + logistic(μ_smkextra) × (1 - P(cig))
P(any tobacco) = P(smoked) + logistic(μ_anyextra) × (1 - P(smoked))

# Smooth Elderly Transition
w(age) = 1 / (1 + exp((age - 65) / 3))
age_effect = w × spline_effect + (1 - w) × linear_effect
```

### Key File Paths

| File | Purpose |
|------|---------|
| `R/pipeline_monolith.R` | Main pipeline (~4,800 lines) |
| `R/_main.R` | Orchestrator for modular execution |
| `R/mcmc_diagnostics.R` | Convergence analysis |
| `R/prediction_all_countries.R` | Prediction functions |
| `data/Prevalence RESHAPE *.xlsx` | Input survey data |
| `results/final_weighted_results_apc_post_selection.csv` | Final results |

### MCMC Configuration

| Parameter | Value | Effective Samples |
|-----------|-------|-------------------|
| Chains | 4 | - |
| Burn-in | 5,000 | (discarded) |
| Iterations | 10,000 | 40,000 total |
| Thinning | 5 | 8,000 retained |

### Common Failure Modes

| Problem | Solution |
|---------|----------|
| "maximal number of DLLs reached" | Set `R_MAX_NUM_DLLS = 600` |
| Initial logProb = -Inf | Check initial values, verify data bounds |
| High R-hat on intercepts | Use data-informed priors (empirical mean) |
| Memory leak in country loop | Use correct cleanup order (OPT-1) |

---

## Reading Guide

### If You Want to Run the Model
1. Section 1: Setup & Configuration
2. Section 2: Data Loading
3. Section 6: Usage (Running the Pipeline)
4. Section 19: Troubleshooting

### If You Want to Understand the Mathematics
1. Section 5: Mathematical Model Specification
2. Section 14: Mathematical Appendix
3. Section 7: NIMBLE Model Code (annotated)

### If You Want to Modify the Model
1. Section 17: Complete Variable Reference
2. Section 7: NIMBLE Model Definitions
3. Section 8: Global Model Fitting
4. Section 11: Country-Specific Models

### If Something Broke
1. Section 19: Troubleshooting
2. Section 18: Common Errors and Solutions
3. Section 11: MCMC Diagnostics

---

## Table of Contents

1. [Pipeline Architecture Overview](#1-pipeline-architecture-overview)
2. [Data Loading and Preparation](#2-data-loading-and-preparation)
3. [Regional Classification System](#3-regional-classification-system)
4. [Age-Cohort Transformation](#4-age-cohort-transformation)
5. [Mathematical Model Specification](#5-mathematical-model-specification)
6. [Global Model Fitting](#6-global-model-fitting)
7. [NIMBLE Model Code Reference](#7-nimble-model-code-reference)
8. [Prior Extraction for Country Models](#8-prior-extraction-for-country-models)
9. [Country-Specific Model Fitting](#9-country-specific-model-fitting)
10. [Vectorized Prediction System](#10-vectorized-prediction-system)
11. [MCMC Diagnostics](#11-mcmc-diagnostics)
12. [Target Evaluation](#12-target-evaluation)
13. [Model Selection (RMSE-Based)](#13-model-selection-rmse-based)
14. [Mathematical Appendix](#14-mathematical-appendix)
15. [Population Weights](#15-population-weights)
16. [Performance Optimizations](#16-performance-optimizations)
17. [Complete Variable Reference](#17-complete-variable-reference)
18. [Common Errors and Solutions](#18-common-errors-and-solutions)
19. [Troubleshooting Guide](#19-troubleshooting-guide)
20. [Output Files Reference](#20-output-files-reference)
21. [Glossary](#21-glossary)

---

## 1. Pipeline Architecture Overview

### 1.1 Pipeline Flow

```
Section 1: Setup & Configuration
       ↓
Section 2: Data Loading & Preparation  ──→  WHO Prevalence Data (Excel)
       ↓
Section 3: Regional Classification     ──→  Country-Region Mapping (14 regions)
       ↓
Section 4: Age-Cohort Transformation   ──→  Spline Bases, Logit Transform
       ↓
Section 5: Population Weights          ──→  UN WPP Age-Sex-Year Weights
       ↓
Section 6: Utility Functions
       ↓
Section 7: NIMBLE Model Definitions    ──→  Global + Country Models
       ↓
Section 8: Global Model Fitting (NIMBLE)
       ├───→ Section 8.13: Prior Extraction (for country models)
       ├───→ Section 8.15: Prediction Generation
       ↓
Section 9: Target Prevalence Calculation
       ↓
Section 10: Weighted Prevalence Trends
       ↓
Section 11: Country-Specific Models (360+ fits)
       ├───→ [OPT-2] Precompute age components
       ├───→ [OPT-3] Vectorized tensor product
       ├───→ [OPT-4] Matrix prediction
       ├───→ [OPT-5] rowMeans summaries
       ├───→ [OPT-6] Precompute cohort splines
       ↓
Section 12: Country-Specific Weighted Trends
       ↓
Section 13: Model Evaluation (RMSE)    ──→  Select best model per country
       ↓
Section 14: Post-Selection Combined Results
       ↓
Section 15-19: Aggregation, Visualization, Tables
```

### 1.2 Configuration Parameters — Complete Reference

#### 1.2.1 Environment Variables

| Variable | Value | Location | Purpose |
|----------|-------|----------|---------|
| `R_MAX_NUM_DLLS` | 600 | Line 101 | Prevents "maximal number of DLLs reached" error. Default R limit is 100; we need ~400 for 191 countries × 2 sexes |

**Why 600?**
- Each NIMBLE model compiles to a separate DLL
- 191 countries × 2 sexes = 382 potential models
- Add safety factor → 600

**Platform differences:**
- Windows: May require Rtools installed
- Linux/Mac: Uses system C++ compiler
- All platforms: Set before loading NIMBLE

#### 1.2.2 NIMBLE Options

| Option | Value | Purpose |
|--------|-------|---------|
| `verbose` | FALSE | Suppresses compilation output |
| `MCMCprogressBar` | TRUE | Shows sampling progress bar |
| `buildInterfacesForCompiledNestedNimbleFunctions` | FALSE | Memory optimization; prevents creation of unnecessary R interfaces for nested compiled functions |

#### 1.2.3 MCMC Parameters — Deep Dive

| Parameter | Value | Rationale | Memory Impact |
|-----------|-------|-----------|---------------|
| `NUMBER_OF_CHAINS` | 4 | Standard for Gelman-Rubin diagnostic (requires ≥2). 4 chains provide robust convergence assessment with reasonable parallelization overhead. | 4× model size |
| `NUMBER_OF_ADAPT` | 1000 | Informational only—NIMBLE handles adaptation internally. This is the period where NIMBLE tunes its proposal distributions. | Minimal |
| `NUMBER_OF_BURN` | 5000 | Empirically sufficient for this model. Can verify with trace plots—chains should stabilize well before 5000. Total discarded: 5000 × 4 = 20,000 samples. | Discarded after run |
| `NUMBER_OF_ITERATIONS` | 10000 | Post-burn-in samples per chain. With 4 chains and thinning of 5: 4 × 10000 / 5 = 8000 retained samples. ESS typically 50-80% of nominal. | Stored in memory |
| `THINNING_INTERVAL` | 5 | Reduces autocorrelation and storage. Rule of thumb: thin by ~autocorrelation time. For this model, lag-5 autocorrelation is typically < 0.1. | 5× storage reduction |

**When to Increase Iterations:**
- ESS < 400 for key parameters
- Trace plots show poor mixing
- R-hat > 1.10 for any parameter

#### 1.2.4 Analysis Parameters

| Parameter | Value | Rationale | WHO Reference |
|-----------|-------|-----------|---------------|
| `BASE_YEAR` | 2010 | WHO baseline year for 30% reduction target | WHA66.10 (2013) |
| `TARGET_YEAR` | 2025 | WHO commitment year | WHO FCTC |
| `PROJECTED_YEARS` | 20 | Extends projections to 2045 for long-term planning | - |
| `REDUCTION_PERCENTAGE` | 30 | WHO 30% relative reduction target | WHA66.10 |
| `RANDOM_SEED` | 42 | Reproducibility; any integer works | - |

#### 1.2.5 Age Transition Parameters

| Parameter | Value | Mathematical Role | Epidemiological Rationale |
|-----------|-------|-------------------|---------------------------|
| `TRANSITION_START` | 65 | Center of sigmoid | Retirement age; mortality selection accelerates; cessation rates increase |
| `TRANSITION_WIDTH` | 3 | Controls sigmoid steepness | Transition zone ≈ 59-71 years (±2σ). Narrower values create sharper transitions; wider values smooth the boundary. |
| `MAX_AGE_SPLINE` | 80 | Upper spline boundary | Beyond 80, pure linear extrapolation. Spline data too sparse above 80. |

**Sigmoid Weight Formula:**
```r
w(age) = 1 / (1 + exp((age - 65) / 3))
```

**Weight Values at Key Ages:**

| Age | Spline Weight | Linear Weight | Interpretation |
|-----|---------------|---------------|----------------|
| 15 | 1.000 | 0.000 | Pure spline |
| 25 | 1.000 | 0.000 | Pure spline |
| 45 | 0.999 | 0.001 | Essentially pure spline |
| 60 | 0.844 | 0.156 | Beginning transition |
| 65 | 0.500 | 0.500 | Midpoint |
| 70 | 0.156 | 0.844 | Approaching linear |
| 75 | 0.035 | 0.965 | Essentially pure linear |
| 85 | 0.001 | 0.999 | Pure linear |
| 100 | 0.000 | 1.000 | Pure linear |

#### 1.2.6 Endgame Target Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| `MANUAL_TARGET_PREVALENCE` | 4.0 | "Tobacco endgame" threshold used by New Zealand, UK |
| `MANUAL_TARGET_PROPORTION` | 0.04 | Derived: 4.0 / 100 |
| `MANUAL_TARGET_EVAL_YEAR` | 2040 | 15 years beyond 2025 target |

---

## 2. Data Loading and Preparation

### 2.1 Input File Specification

**File**: `data/Prevalence RESHAPE 28 Sep 2024.xlsx`
**Sheet**: `reshape.dta`
**Skip**: 1 row (header in Excel format)

#### 2.1.1 Required Columns

| Column | Type | Allowed Values | Example | Notes |
|--------|------|----------------|---------|-------|
| `un_name` | string | Country names | "United States" | Converted to lowercase |
| `wb_country_abv` | string | ISO3 codes | "usa", "gbr" | Primary key, 3 lowercase letters |
| `un_region` | string | UN regions | "Americas" | Replaced with custom classification |
| `survey_year` | numeric | 1980–2024 | 2018 | Must be valid year |
| `title_pri` | string | Survey name | "NHANES" | Creates survey random effect |
| `sex_name` | string | "male", "female" | "male" | Lowercase |
| `start_age` | numeric | 15–100 | 25 | Lower bound of age band |
| `end_age` | numeric | 15–100+ | 34 | Upper bound of age band |
| `def_code` | string | "cd.0101", "cd.0112" | "cd.0112" | Definition code |
| `type_code` | string | "tt.001", "tt.002", "tt.003" | "tt.003" | Product type code |
| `prevalence` | numeric | 0–100 | 23.5 | Percentage scale |

### 2.2 Data Cleaning Steps

#### 2.2.1 Column Renaming

| Original | Renamed To |
|----------|------------|
| `un_name` | `country` |
| `un_region` | `region` |
| `title_pri` | `survey` |
| `survey_year` | `year` |
| `sex_name` | `sex` |

#### 2.2.2 Type Conversions

```r
start_age  = as.numeric(start_age)
end_age    = as.numeric(end_age)
prevalence = as.numeric(prevalence)
year       = as.numeric(year)
```

#### 2.2.3 Indicator Creation

**Definition Code Binary** (`def_code_binary`):
| `def_code` | Binary Value | Meaning |
|------------|--------------|---------|
| `cd.0101` | 0 | Daily user |
| `cd.0112` | 1 | Current user |

**Why binary?** This indicator enters the linear predictor directly. Current users include daily users plus occasional users, so prevalence is higher when `def_code_binary = 1`.

**Product Type Indicators** (mutually exclusive):
| Indicator | Type Code | Product |
|-----------|-----------|---------|
| `type_any_tobacco` | `tt.001` | Any tobacco product |
| `type_smoked_tobacco` | `tt.002` | Any smoked tobacco |
| `type_cigarettes` | `tt.003` | Cigarettes only |

**Validation:** Exactly one indicator should equal 1 for each observation. Observations with all zeros (unrecognized type codes) are removed with a warning.

### 2.3 Age Midpoint Calculation

```r
Age_Midpoint = (start_age + end_age) / 2
```

**Filter:** `Age_Midpoint >= 15` (removes child observations)

### 2.4 Redundant Observation Removal

The data may contain both total-population observations (ages 15-99) AND age-stratified observations for the same country/year/survey. The model cannot use both without double-counting.

**Algorithm:**
```r
for each observation with (start_age = 15) AND (end_age >= 99):
  if other observations exist for same (sex, country, def_type, year):
    drop the total-population observation
```

**Rationale:**
- Age-stratified data provides more information
- Total-population observations have different variance structure
- Keeping both would double-count prevalence signal

### 2.5 Prevalence Transformation

```r
# Step 1: Convert from percentage to proportion
prevalence = prevalence / 100

# Step 2: Bound away from 0 and 1 (logit is undefined there)
prevalence = ifelse(prevalence == 0, 0.001, prevalence)
prevalence = ifelse(prevalence == 1, 0.999, prevalence)

# Step 3: Logit transform
prevalence = log(prevalence / (1 - prevalence))
```

**Why bound at 0.001 and 0.999?**
- logit(0) = -∞, logit(1) = +∞
- 0.001 → logit ≈ -6.9
- 0.999 → logit ≈ +6.9
- These bounds are wide enough to not distort estimates but prevent numerical issues

### 2.6 Country Name Mapping

```r
country_name_mapping <- clean_data %>%
  distinct(wb_country_abv, country) %>%
  { setNames(.$country, .$wb_country_abv) }
```

**Structure:** Named character vector
- Names: ISO3 country codes (e.g., "usa")
- Values: Full country names (e.g., "united states")

**Usage:** Convert between codes and names throughout pipeline

---

## 3. Regional Classification System

### 3.1 Complete Country-Region Mapping

The model uses 14 WHO-inspired regional groupings:

| Region | Countries (ISO3 codes) |
|--------|------------------------|
| **Central Asia** | tjk, tkm, kaz |
| **Eastern Asia** | chn, prk, jpn, mng, kor |
| **Eastern Europe** | blr, bgr, cze, hun, mda, pol, rou, rus, svk, ukr |
| **Latin America & Caribbean** | bhs, brb, blz, cub, dma, slv, grd, gtm, hti, hnd, jam, mex, nic, kna, lca, vct, tto, dom, cri, pan |
| **North Africa & Middle East** | dza, egy, lby, mar, sdn, tun, pse |
| **North America** | usa, can |
| **Oceania & Pacific** | kir, mhl, fsm, nru, niu, plw, png, wsm, slb, ton, vut, cok, tuv, nzl, aus, fji |
| **South America** | arg, bol, bra, chl, col, ecu, guy, per, sur, ury, ven, pry |
| **Southeastern Asia** | brn, idn, lao, mys, mmr, sgp, tha, vnm, phl, tls, khm |
| **Southern Asia** | afg, bgd, ind, irn, kgz, mdv, pak, lka, uzb, npl, btn |
| **Southern Europe** | alb, and, bih, grc, hrv, ita, mlt, mne, mkd, prt, smr, srb, svn, esp |
| **Sub-Saharan Africa** | ago, ben, bwa, bfa, bdi, cpv, cmr, caf, tcd, com, cog, cod, civ, dji, gnq, eri, eth, gab, gmb, gha, gin, gnb, ken, lso, lbr, mdg, mwi, mli, mrt, mus, moz, nam, ner, nga, rwa, stp, sen, syc, sle, zaf, swz, tza, tgo, uga, zmb, zwe |
| **Western Asia** | arm, aze, bhr, cyp, geo, irq, isr, jor, kwt, lbn, omn, qat, sau, syr, tur, are, yem |
| **Western Europe** | aut, bel, deu, fra, gbr, irl, nld, che, lux |
| **Northern Europe** | dnk, fin, isl, nor, swe, est, lva, ltu |

### 3.2 Region Numeric Encoding

For NIMBLE, regions must be integer-indexed:

```r
unique_regions <- sort(unique(country_region_mapping$region_consolidated))
region_to_num <- setNames(seq_along(unique_regions), unique_regions)
```

**Critical Fix (v2.3.1):** `Country_Region` must be length `nCountry`, not length `N`:

```r
# CORRECT: One region index per country
country_region_for_model <- country_lookup_df %>%
  arrange(Num_Country) %>%
  pull(Region_Num)

# This vector is indexed by Num_Country, not by observation i
Country_Region[Country[i]]  # ← How it's used in the model
```

---

## 4. Age-Cohort Transformation

### 4.1 Birth Cohort Calculation

```r
Birth_Cohort = year - Age_Midpoint
```

**Example:**
- Survey year 2015, age midpoint 30 → Birth cohort 1985
- Survey year 2000, age midpoint 45 → Birth cohort 1955

**Typical range in dataset:** ~1890 to ~2010

### 4.2 Centering Constants

| Constant | Formula | Stored Globally | Purpose |
|----------|---------|-----------------|---------|
| `COHORT_CENTER_CONSTANT` | `median(Birth_Cohort)` | Yes | Centers cohort spline |
| `AGE_LINEAR_CENTER_CONSTANT` | `mean(Age_Linear)` | Yes | Centers linear age term |

**Why median for cohort?** Robustness to outliers (very old surveys)

**Why mean for age linear?** Symmetric linear term; mean is natural center

**Critical:** These constants MUST be stored and reused for predictions. Using different centering in prediction vs. fitting will produce biased results.

### 4.3 Sigmoid Weight Calculation — Step by Step

```r
# Step 1: Compute sigmoid input
sigmoid_input = (Age_Midpoint - TRANSITION_START) / TRANSITION_WIDTH
# For age 65: sigmoid_input = (65 - 65) / 3 = 0
# For age 70: sigmoid_input = (70 - 65) / 3 = 1.67

# Step 2: Apply logistic function
spline_weight = 1 / (1 + exp(sigmoid_input))
# For age 65: spline_weight = 1 / (1 + exp(0)) = 0.5
# For age 70: spline_weight = 1 / (1 + exp(1.67)) = 0.16

# Step 3: Compute complement
linear_weight = 1 - spline_weight
```

### 4.4 Age Linear Term

```r
# Step 1: Linear term starts at transition point
Age_Linear = pmax(0, Age_Midpoint - TRANSITION_START)
# For age 50: Age_Linear = max(0, 50 - 65) = 0
# For age 75: Age_Linear = max(0, 75 - 65) = 10

# Step 2: Center for numerical stability
Age_Linear_Centered = Age_Linear - AGE_LINEAR_CENTER_CONSTANT
```

**Interpretation:** `Age_Linear` is zero for ages below 65, then increases linearly. The coefficient on this term (constrained negative) ensures prevalence declines at older ages.

### 4.5 Natural Spline Basis — Deep Dive

#### What is a Natural Cubic Spline?

A natural cubic spline with interior knots at $k_1, \ldots, k_K$ and boundary knots at $B_1, B_2$ is a function $f(x)$ such that:

1. **Cubic between knots:** $f(x)$ is a cubic polynomial on each interval $[k_i, k_{i+1}]$
2. **Continuous through second derivative:** $f$, $f'$, and $f''$ are continuous at all knots
3. **Linear beyond boundaries:** $f''(x) = 0$ for $x < B_1$ and $x > B_2$

The `ns()` function in R returns a **basis matrix** where each column is one basis function evaluated at the input values.

#### R Implementation

```r
age_spline_basis <- ns(
  clean_data$Age_For_Spline,    # Vector of ages
  knots = c(25, 45, 65),        # Interior knots
  Boundary.knots = c(15, 80)    # Boundary knots
)
```

**Output dimensions:** `N × (K + 1)` where K = number of interior knots
- With 3 interior knots: 4 basis functions
- Each column represents one basis function

#### Attribute Preservation

```r
# These attributes MUST be saved for prediction
age_spline_knots_attr    <- attr(age_spline_basis, "knots")
age_spline_boundary_attr <- attr(age_spline_basis, "Boundary.knots")

# In prediction, use the SAME attributes:
ns(new_ages,
   knots = age_spline_knots_attr,
   Boundary.knots = age_spline_boundary_attr)
```

### 4.6 Knot Placement — Literature Justification

#### Age Knots: 25, 45, 65

| Knot | Inflection Point | Evidence |
|------|------------------|----------|
| **25** | End of initiation | >80% of ever-smokers begin by age 25 (NHANES, GBD 2015). After 25, initiation is rare; prevalence stabilizes or begins declining. |
| **45** | Cessation onset | Cessation rates lowest between ages 45-64 ("hardening" hypothesis). After 45, net prevalence begins declining as cessation exceeds initiation. |
| **65** | Mortality selection | Differential mortality accelerates apparent prevalence decline. Smokers die younger, so surviving cohort has lower prevalence. |

**References:**
- GBD 2015 Lancet: "Male prevalence peaks between ages 25 and 35"
- NHANES 2025 Sci Rep: "Odds of smoking decrease after age ~27"

#### Cohort Knots: 1945, 1965, 1985

| Knot | Historical Context | Evidence |
|------|-------------------|----------|
| **1945** | Post-war peak | Baby boom generation. Female smoking peaked for 1935-1945 cohorts (Holford 2014). |
| **1965** | Post–Surgeon General | First cohort born after 1964 Surgeon General's Report. Reached adulthood amid growing health awareness. |
| **1985** | Modern tobacco control | Cohorts reaching adulthood from 2000 faced taxes, advertising bans, smoke-free policies. |

### 4.7 Age-Cohort Interaction Tensor

```r
n_age_splines    <- 4  # ncol(age_spline_basis)
n_cohort_splines <- 4  # ncol(cohort_spline_basis)
n_interactions   <- 16 # 4 × 4

# Build tensor product row by row
for (i in 1:N) {
  age_values    <- age_spline_basis[i, ]     # 1 × 4
  cohort_values <- cohort_spline_basis[i, ]  # 1 × 4

  # outer() creates 4×4 matrix, flatten to 1×16 vector
  age_cohort_interaction_matrix[i, ] <- as.vector(outer(age_values, cohort_values))
}
```

**Column ordering:** The tensor product is stored in column-major order:
- Columns 1-4: age_1 × (cohort_1, cohort_2, cohort_3, cohort_4)
- Columns 5-8: age_2 × (cohort_1, cohort_2, cohort_3, cohort_4)
- etc.

This ordering must match between fitting and prediction.

---

## 5. Mathematical Model Specification

### 5.1 Three-Head Stick-Breaking Construction

The model estimates three nested tobacco categories simultaneously:

```
Any Tobacco Product (tt.001)
       ↓
Any Smoked Tobacco (tt.002)
       ↓
Cigarettes (tt.003)
```

By definition: $P(\text{cig}) \leq P(\text{smoked}) \leq P(\text{any})$

The **stick-breaking** parameterization guarantees this ordering:

$$
\begin{aligned}
P(\text{cigarettes}) &= \text{logistic}(\mu_{\text{cig}}) \\
P(\text{any smoked}) &= P(\text{cig}) + \text{logistic}(\mu_{\text{smkextra}}) \cdot (1 - P(\text{cig})) \\
P(\text{any tobacco}) &= P(\text{smoked}) + \text{logistic}(\mu_{\text{anyextra}}) \cdot (1 - P(\text{smoked}))
\end{aligned}
$$

**Interpretation:**
- $\text{logistic}(\mu_{\text{cig}})$: Probability of smoking cigarettes
- $\text{logistic}(\mu_{\text{smkextra}})$: Probability of "other smoked tobacco" given "not cigarettes"
- $\text{logistic}(\mu_{\text{anyextra}})$: Probability of "smokeless tobacco" given "not smoked"

### 5.2 Linear Predictor — Full Equation

For observation $i$ from country $j$ in region $r$:

$$
\mu_{\text{cig}}[i] = \underbrace{\alpha^{(G)}}_{\text{global}} + \underbrace{\alpha^{(R)}_r}_{\text{region}} + \underbrace{\alpha^{(C)}_j}_{\text{country}} + \underbrace{\beta \cdot D_i}_{\text{definition}} + \underbrace{w_i \cdot \mathbf{s}_a(a_i)^\top \boldsymbol{\gamma}_j}_{\text{age spline}} + \underbrace{(1-w_i) \cdot \delta \cdot \tilde{a}_i}_{\text{linear tail}} + \underbrace{\mathbf{s}_c(c_i)^\top \boldsymbol{\phi}_j}_{\text{cohort spline}} + \underbrace{\mathbf{z}_i^\top \boldsymbol{\psi}}_{\text{age-cohort}} + \underbrace{\eta_k}_{\text{survey}}
$$

### 5.3 Model Complexity by Head

| Component | CIG | SMKEXTRA | ANYEXTRA |
|-----------|-----|----------|----------|
| Global intercept | ✓ | ✓ | ✓ |
| Regional intercept | ✓ | ✓ | ✓ |
| Country intercept | ✓ | ✓ | ✓ |
| Definition effect | Full | 0.3× | 0.3× |
| Age splines | Country | Regional | Regional |
| Cohort splines | Country | Regional | Regional |
| Age-cohort interaction | ✓ (16 params) | ✗ | ✗ |
| Age linear (elderly) | ✓ | ✓ | ✓ |

**Why simplified for SMKEXTRA/ANYEXTRA?**
- Signal is weaker (fewer observations of non-cigarette products)
- Country-specific splines don't converge well
- Regional pooling improves stability

### 5.4 Likelihood

```
logit_p[i] = log(p_obs[i] / (1 - p_obs[i]))
tau[i] = (1 / residual_sd²) × weight[i]
Prevalence[i] ~ Normal(logit_p[i], tau[i])
```

Where:
- `p_obs[i]` = appropriate probability based on product type
- `weight[i]` = inverse age-band width (downweights wide bands)

### 5.5 Prior Specifications — Complete Table

| Parameter | Prior | Hyperparameters | Rationale |
|-----------|-------|-----------------|-----------|
| `cig_global_intercept` | N(emp_mean, 0.3²) | emp_mean = mean(logit(prevalence)) | Data-informed anchor; breaks non-identifiability |
| `smkextra_global_intercept` | N(-3.0, 0.3²) | | P ≈ 5% on probability scale |
| `anyextra_global_intercept` | N(-4.0, 0.3²) | | P ≈ 2% on probability scale |
| `cig_def_code_shared` | N(0.3, 1.0²) | | Current > Daily; modest effect |
| `residual_sd` | LogNormal(log(0.7), 0.5²) | | Mode ≈ 0.7; right-skewed |
| `*_age_linear_smooth_effect` | TruncN(-0.02, 0.05², upper=-0.001) | | Constrained negative |
| `intercept_between_region_precision` | Gamma(4, 1) | | Mean = 4, mode = 3 |
| `*_age_spline_global_mean` | N(0, 1.5²) | | Weakly informative |
| `cig_age_cohort_interaction[k]` | N(0, 0.2²) | | Tight prior; regularization |
| `survey_intercept[s]` | N(0, survey_sd) | survey_sd ~ 1/√Gamma(3,1) | Absorbs survey-specific bias |

---

## 6. Global Model Fitting

### 6.1 NIMBLE Workflow Overview

```r
# 1. Build model (parse, validate)
nimble_model <- nimbleModel(code, constants, data, inits, name)

# 2. Configure MCMC (select samplers, set monitors)
mcmc_config <- configureMCMC(nimble_model, monitors, thin)

# 3. Build MCMC object
mcmc_built <- buildMCMC(mcmc_config)

# 4. Compile to C++ (slow, ~5-10 minutes)
compiled_model <- compileNimble(nimble_model)
compiled_mcmc <- compileNimble(mcmc_built, project = nimble_model)

# 5. Run sampling
samples <- runMCMC(compiled_mcmc, niter, nburnin, nchains, inits)
```

### 6.2 Constants vs Data Separation

**NIMBLE requires separating structure from observations:**

**Constants** (define model structure):
- Loop bounds: `N`, `nCountry`, `nRegion`, `nSurvey`
- Spline dimensions: `nAgeSpline`, `nCohortSpline`, `nAgeXCohortSplines`
- Index vectors: `Country[i]`, `Country_Region[j]`, `Survey[i]`
- Data-informed hyperparameters: `empirical_mean_cig`

**Data** (observed values):
- Response: `Prevalence[i]`
- Covariates: `Def_Code_Binary[i]`, `Type_Cig[i]`, `age_spline_matrix[i,l]`
- Weights: `weight[i]`

### 6.3 Initial Values Strategy

Initial values are "warm-started" near prior means with small jitter:

```r
generate_global_inits <- function(seed_offset, emp_mean) {
  set.seed(RANDOM_SEED + seed_offset)
  list(
    # Start at prior mean ± small noise
    cig_global_intercept = emp_mean + rnorm(1, 0, 0.05),
    smkextra_global_intercept = -3.0 + rnorm(1, 0, 0.05),
    anyextra_global_intercept = -4.0 + rnorm(1, 0, 0.05),

    # Constrained parameters: start in valid region
    cig_age_linear_smooth_effect = -0.02 + rnorm(1, 0, 0.005),

    # Precision parameters: start near prior mean
    intercept_between_region_precision = 4 + rnorm(1, 0, 0.3)
  )
}

# Different seed for each chain
inits_list <- lapply(1:NUMBER_OF_CHAINS, function(i) {
  generate_global_inits(i, emp_mean = empirical_mean_cig)
})
```

### 6.4 MCMC Configuration Options

| Option | Value | Purpose |
|--------|-------|---------|
| `monitors` | (list of parameters) | Which parameters to save |
| `thin` | 5 | Keep every 5th sample |
| `enableWAIC` | FALSE | WAIC computation disabled (slow) |
| `useConjugacy` | FALSE | Conjugate sampling disabled; use default Metropolis-Hastings |

**Why `useConjugacy = FALSE`?**
Conjugate samplers can be faster but don't exist for all distributions in this model. Disabling ensures consistent sampler selection across all parameters.

---

## 7. NIMBLE Model Code Reference

### 7.1 Global Model — Annotated Code

```r
regional_hierarchical_global_ac_model_nimble <- nimbleCode({

  # ==================================================================
  # LIKELIHOOD LOOP — Process each observation
  # ==================================================================
  for (i in 1:N) {

    # ----------------------------------------------------------------
    # HEAD 1: CIGARETTES (Full complexity)
    # ----------------------------------------------------------------
    mu_cig[i] <- cig_global_intercept +                           # Global mean
      cig_region_intercept[Country_Region[Country[i]]] +          # Region deviation
      cig_country_intercept[Country[i]] +                         # Country deviation
      cig_def_code_shared * Def_Code_Binary[i] +                  # Definition effect
      spline_weight_var[i] * inprod(                              # Age spline (weighted)
        cig_age_spline[Country[i], 1:nAgeSpline],
        age_spline_matrix[i, 1:nAgeSpline]
      ) +
      linear_weight_var[i] * cig_age_linear_smooth_effect *       # Linear tail
        age_linear_smooth[i] +
      inprod(                                                     # Cohort spline
        cig_cohort_spline[Country[i], 1:nCohortSpline],
        cohort_spline_matrix[i, 1:nCohortSpline]
      ) +
      inprod(                                                     # Age-cohort interaction
        cig_age_cohort_interaction[1:nAgeXCohortSplines],
        age_cohort_interaction_matrix[i, 1:nAgeXCohortSplines]
      ) +
      survey_intercept[Survey[i]]                                 # Survey random effect

    # ----------------------------------------------------------------
    # HEAD 2: OTHER SMOKED (Simplified — regional splines)
    # ----------------------------------------------------------------
    mu_smkextra[i] <- smkextra_global_intercept +
      smkextra_region_intercept[Country_Region[Country[i]]] +
      smkextra_country_intercept[Country[i]] +
      0.3 * cig_def_code_shared * Def_Code_Binary[i] +           # 30% of def effect
      spline_weight_var[i] * inprod(
        smkextra_age_spline_region_mean[Country_Region[Country[i]], 1:nAgeSpline],
        age_spline_matrix[i, 1:nAgeSpline]
      ) +
      linear_weight_var[i] * smkextra_age_linear_smooth_effect *
        age_linear_smooth[i] +
      inprod(
        smkextra_cohort_spline_region_mean[Country_Region[Country[i]], 1:nCohortSpline],
        cohort_spline_matrix[i, 1:nCohortSpline]
      )
      # NOTE: No age-cohort interaction

    # ----------------------------------------------------------------
    # HEAD 3: NON-SMOKED (Simplified — regional splines)
    # ----------------------------------------------------------------
    mu_anyextra[i] <- anyextra_global_intercept +
      anyextra_region_intercept[Country_Region[Country[i]]] +
      anyextra_country_intercept[Country[i]] +
      0.3 * cig_def_code_shared * Def_Code_Binary[i] +
      spline_weight_var[i] * inprod(
        anyextra_age_spline_region_mean[Country_Region[Country[i]], 1:nAgeSpline],
        age_spline_matrix[i, 1:nAgeSpline]
      ) +
      linear_weight_var[i] * anyextra_age_linear_smooth_effect *
        age_linear_smooth[i] +
      inprod(
        anyextra_cohort_spline_region_mean[Country_Region[Country[i]], 1:nCohortSpline],
        cohort_spline_matrix[i, 1:nCohortSpline]
      )
      # NOTE: No age-cohort interaction

    # ----------------------------------------------------------------
    # STICK-BREAKING CONSTRUCTION
    # ----------------------------------------------------------------
    p_cig[i]    <- ilogit(mu_cig[i])
    p_anysmk[i] <- p_cig[i] + ilogit(mu_smkextra[i]) * (1 - p_cig[i])
    p_anytob[i] <- p_anysmk[i] + ilogit(mu_anyextra[i]) * (1 - p_anysmk[i])

    # Select observed probability based on indicator
    p_obs[i] <- Type_Cig[i] * p_cig[i] +
                Type_Smoked[i] * p_anysmk[i] +
                Type_Any[i] * p_anytob[i]

    # ----------------------------------------------------------------
    # LIKELIHOOD
    # ----------------------------------------------------------------
    logit_p[i] <- log(p_obs[i] / (1 - p_obs[i]))
    tau[i] <- pow(residual_sd, -2) * weight[i]
    Prevalence[i] ~ dnorm(logit_p[i], tau = tau[i])
  }

  # [PRIORS omitted for brevity — see Section 5.5]
})
```

### 7.2 Key NIMBLE Functions

| NIMBLE | R Equivalent | Notes |
|--------|--------------|-------|
| `ilogit(x)` | `plogis(x)` | Inverse logit (logistic) |
| `inprod(a, b)` | `sum(a * b)` | Inner product |
| `dnorm(x, mean, tau)` | - | Parameterized by precision, not SD |
| `T(dist, lower, upper)` | - | Truncated distribution |
| `pow(x, y)` | `x^y` | Power function |

---

## 8. Prior Extraction for Country Models

### 8.1 Purpose

After fitting the global model, extract country-specific posterior distributions to serve as informative priors for country-specific models.

### 8.2 Total Intercept Calculation — Sample-Wise

```r
# CORRECT: Compute sum for each MCMC sample (preserves correlations)
compute_total_intercept_samplewise <- function(head_prefix) {
  global_col  <- paste0(head_prefix, "_global_intercept")
  region_col  <- paste0(head_prefix, "_region_intercept[", country_region, "]")
  country_col <- paste0(head_prefix, "_country_intercept[", country_num, "]")

  # Sum sample-by-sample
  total_samples <- combined_samples_matrix[, global_col] +
                   combined_samples_matrix[, region_col] +
                   combined_samples_matrix[, country_col]

  list(
    mean = mean(total_samples),
    sd   = sd(total_samples)
  )
}
```

**Why sample-wise?**
If we computed means separately and then summed, we would assume independence:
- $\text{Var}(A + B + C) = \text{Var}(A) + \text{Var}(B) + \text{Var}(C)$ (only if independent)

But these parameters are correlated. Sample-wise computation preserves the true posterior variance.

### 8.3 Shrinkage Factor

```r
shrinkage = 0.7
prior_sd = posterior_sd * shrinkage
```

**Why shrink?**
- Prevents country model from wandering too far from global fit
- Allows meaningful deviation when data warrant it
- Empirically tuned; 0.7 balances flexibility and stability

### 8.4 Prior File Format

Output: `country_priors/{gender}/{country_code}_regional_ac_priors_nested.csv`

| Column | Description |
|--------|-------------|
| `parameter` | Parameter name (e.g., "cig_intercept", "cig_age_spline_1") |
| `mean` | Posterior mean from global model |
| `sd` | Posterior SD × shrinkage factor |

---

## 9. Country-Specific Model Fitting

### 9.1 Loop Structure

```r
for (gender in c("males", "females")) {
  # Load regional structure (shared centering constants)
  regional_info <- readRDS(...)

  # [OPT-2] Precompute age-only components (shared across countries)

  for (country_code in countries_with_data) {
    # Load country-specific priors
    # Filter country data
    # Build NIMBLE model with informative priors
    # Run MCMC
    # Generate predictions
    # [OPT-1] Clean up compiled objects
  }
}
```

### 9.2 Model Structure Differences from Global

| Aspect | Global Model | Country Model |
|--------|--------------|---------------|
| Intercepts | Global + Region + Country | Single country intercept |
| Splines | Country-specific | Country-specific |
| Priors | Weakly informative | Informative (from global) |
| Regions | Modeled | Not needed (single country) |
| Parameters | ~3000+ | ~50 per country |

### 9.3 Precision Conversion

```r
sd_to_prec <- function(sd_val) {
  sd_safe <- pmax(sd_val, 0.01)    # Minimum SD = 0.01
  prec <- 1 / (sd_safe^2)
  return(pmin(prec, 10000))        # Maximum precision = 10000
}
```

**Why clamp?**
- Very small SD → very large precision → numerical instability
- Minimum SD of 0.01 corresponds to maximum precision of 10000
- Prevents NIMBLE from crashing on degenerate priors

---

## 10. Vectorized Prediction System

### 10.1 Overview

The prediction system generates prevalence estimates for all combinations of:
- Ages: 15, 16, ..., 100 (86 ages)
- Years: min_year to target_year + 20
- Definition codes: daily, current
- Product types: cigarettes, smoked, any

**Output:** Matrices of shape `(n_ages × n_samples)` saved as RDS files.

### 10.2 Precomputed Components — OPT-2

```r
# These are computed ONCE outside the country loop

# Age grid (shared across all countries)
age_midpoints_shared <- seq(15, 100, by = 1)  # 86 values

# Sigmoid weights (computed once)
sigmoid_input_shared <- (age_midpoints_shared - 65) / 3
spline_weight_shared <- 1 / (1 + exp(sigmoid_input_shared))
linear_weight_shared <- 1 - spline_weight_shared

# Age linear term (computed once)
age_linear_shared <- pmax(0, age_midpoints_shared - 65)
age_linear_centered_shared <- age_linear_shared - AGE_LINEAR_CENTER_CONSTANT

# Pre-multiply linear components
linear_age_product_shared <- linear_weight_shared * age_linear_centered_shared

# Age spline basis (computed once)
age_spline_mat_shared <- as.matrix(ns(
  age_midpoints_shared,
  knots = regional_info$age_spline_knots,
  Boundary.knots = regional_info$age_spline_boundary
))
```

### 10.3 Cohort Spline Precomputation — OPT-6

```r
# Compute ALL birth cohorts for ALL years at once
years_pred <- seq(min_year, target_year + 20, by = 1)
n_years <- length(years_pred)

# Flatten: (n_years × n_ages) vector
all_birth_cohorts <- rep(years_pred, each = n_ages) -
                     rep(age_midpoints_shared, times = n_years)

# Single call to ns() for all cohort values
all_cohort_spline_mat <- as.matrix(ns(
  all_birth_cohorts,
  knots = regional_info$cohort_spline_knots,
  Boundary.knots = regional_info$cohort_spline_boundary
))

# Index by year for later use
cohort_spline_by_year <- split(rows, year_indices)
```

### 10.4 Vectorized Tensor Product — OPT-3

```r
# Original (slow): row-by-row outer products
for (i in 1:n_ages) {
  int_mat[i, ] <- as.vector(outer(age_spline[i,], cohort_spline[i,]))
}

# Optimized: column multiplication
for (m in 1:nCohortSpline) {
  col_start <- (m - 1) * nAgeSpline + 1
  col_end   <- m * nAgeSpline
  int_mat[, col_start:col_end] <- age_spline_mat * cohort_mat[, m]
}
```

### 10.5 Vectorized Prediction — OPT-4

```r
# Full matrix algebra: n_ages × n_samples result

# CIG linear predictor
mu_cig <- ones_ages %o% cig_int_s +                              # Intercept broadcast
  def_code_binary * (ones_ages %o% def_shared_s) +               # Definition effect
  (age_spline_mat %*% t(cig_age_spline_s)) * spline_weight +     # Age spline
  outer(linear_age_product, cig_age_lin_s) +                     # Linear tail
  cohort_mat %*% t(cig_cohort_s) +                               # Cohort spline
  int_mat %*% t(cig_ac_s)                                        # Age-cohort

# Stick-breaking (element-wise, all matrices)
p_cig    <- plogis(mu_cig)
p_smoked <- p_cig + plogis(mu_smkextra) * (1 - p_cig)
p_any    <- p_smoked + plogis(mu_anyextra) * (1 - p_smoked)
```

### 10.6 Vectorized Summary Statistics — OPT-5

```r
# Use rowMeans instead of apply(..., mean) for ~5x speedup
summarise_predictions <- function(p_matrix) {
  ci <- t(apply(p_matrix, 1, quantile, probs = c(0.025, 0.975)))
  data.frame(
    Prevalence = rowMeans(p_matrix),  # ← FAST
    lower_ci   = ci[, 1],
    upper_ci   = ci[, 2]
  )
}
```

---

## 11. MCMC Diagnostics

### 11.1 R-hat Interpretation

| R-hat Range | Grade | Interpretation | Action |
|-------------|-------|----------------|--------|
| < 1.05 | Excellent | Full convergence | None needed |
| 1.05–1.10 | Good | Acceptable | Monitor but OK |
| 1.10–1.50 | Concerning | Possible non-convergence | Increase iterations |
| 1.50–2.00 | Poor | Results questionable | Investigate/reparameterize |
| > 2.00 | Failed | Results unreliable | Do not use |

### 11.2 ESS Interpretation

| ESS Range | Interpretation | Action |
|-----------|----------------|--------|
| > 400 | Good | Reliable posterior summaries |
| 200–400 | Adequate | Point estimates OK; CIs may be rough |
| 100–200 | Marginal | Consider more iterations |
| < 100 | Insufficient | Increase iterations or reparameterize |

### 11.3 Parameter Groups

The diagnostics module classifies parameters into 30 groups for structured analysis:

| Group | Prefix | Description |
|-------|--------|-------------|
| 01 | `cig_global_intercept` | Global mean cigarette prevalence |
| 05 | `cig_region_intercept` | Regional deviations (14 values) |
| 08 | `cig_country_intercept` | Country deviations (~180 values) |
| 11 | `cig_country_age_spline` | Country age splines (~180 × 4) |
| 19 | `cig_age_cohort_interaction` | Age-cohort interactions (16 values) |

---

## 12. Target Evaluation

### 12.1 WHO 30% Reduction Target

```r
# For each country and posterior sample
target_prevalence = base_year_prevalence * (1 - 0.30)
prob_achieving = mean(pred_2025 < target_prevalence)
```

**Output interpretation:**
- `prob_achieving_target = 0.85`: 85% posterior probability of achieving ≥30% reduction
- `prob_achieving_target = 0.12`: 12% probability; unlikely to achieve target

### 12.2 Population-Weighted Prevalence

```r
calculate_weighted_prevalence <- function(iterations_matrix_logit, weights_vector) {
  p_matrix <- plogis(iterations_matrix_logit)  # n_ages × n_samples

  weighted_probs <- numeric(ncol(p_matrix))
  for (s in 1:ncol(p_matrix)) {
    weighted_probs[s] <- sum(p_matrix[, s] * weights_vector) / sum(weights_vector)
  }

  return(weighted_probs)  # n_samples vector
}
```

---

## 13. Model Selection (RMSE-Based)

### 13.1 RMSE Calculation

```r
# Match predictions to observed data points
matched_data <- join(observed, predicted, by = c("Year", "Age_Midpoint", "Indicator"))

# Calculate RMSE
global_rmse <- sqrt(mean((global_pred - observed)^2))
country_rmse <- sqrt(mean((country_pred - observed)^2))

# Selection rule
selected_model <- if(country_rmse < global_rmse) "Country" else "Global"
```

### 13.2 Selection Summary

Typical results:
- ~60-70% of country-sex-indicator combinations select "Global"
- Country model wins when sufficient local data exists
- Global model wins for data-sparse countries (borrowing strength)

---

## 14. Mathematical Appendix

### 14.1 Natural Spline Properties

A natural cubic spline $f(x)$ with interior knots $\{k_1, \ldots, k_K\}$ satisfies:

1. $f(x)$ is a cubic polynomial on each interval $[k_i, k_{i+1}]$
2. $f$, $f'$, $f''$ are continuous everywhere
3. $f''(x) = 0$ for $x < k_1$ and $x > k_K$ (linear extrapolation)

**Basis dimension:** $K + 1$ basis functions for $K$ interior knots

### 14.2 Stick-Breaking Proof of Ordering

**Claim:** For any $\mu_1, \mu_2, \mu_3 \in \mathbb{R}$, the stick-breaking construction guarantees $p_1 \leq p_2 \leq p_3$.

**Proof:**
$$
p_2 = p_1 + \underbrace{\text{logistic}(\mu_2)}_{\in (0,1)} \cdot \underbrace{(1 - p_1)}_{\geq 0} \geq p_1
$$

Since logistic is always in $(0,1)$ and $1-p_1 \geq 0$, the second term is non-negative. Similarly for $p_3 \geq p_2$. ∎

### 14.3 Sigmoid Transition Properties

$$
w(a) = \frac{1}{1 + e^{(a-65)/3}}
$$

**Derivative at midpoint:**
$$
w'(65) = -\frac{1}{3} \cdot w(65) \cdot (1 - w(65)) = -\frac{1}{12}
$$

**Limiting behavior:**
- $\lim_{a \to -\infty} w(a) = 1$ (pure spline for young ages)
- $\lim_{a \to +\infty} w(a) = 0$ (pure linear for old ages)

---

## 15. Population Weights

### 15.1 Weight File Format

**File:** `data/weights_15_2022.csv`

| Column | Description |
|--------|-------------|
| `area` | ISO3 country code |
| `year` | Calendar year |
| `sex` | "male" or "female" |
| `15`, `16`, ..., `100` | Population count at each age |

**Source:** UN World Population Prospects

### 15.2 Observation Weights

Within the model, observations are weighted by inverse age-band width:

```r
age_range = end_age - start_age
weight = 1 / (age_range + 1)
weight = weight / mean(weight)  # Normalize
```

**Rationale:** Wide age bands (e.g., 15-99) have more averaging error than narrow bands (e.g., 25-29).

---

## 16. Performance Optimizations

### 16.1 Summary of Optimizations

| Code | Name | Before | After | Speedup |
|------|------|--------|-------|---------|
| OPT-1 | Cleanup order fix | DLL leak | DLL released | Memory stability |
| OPT-2 | Age precomputation | 86 ns() calls/country | 1 ns() call total | ~50× |
| OPT-3 | Tensor vectorization | Row-by-row outer() | Matrix column mult | ~5× |
| OPT-4 | Matrix prediction | Row-by-row loop | Matrix algebra | ~10× |
| OPT-5 | rowMeans vs apply | apply(..., mean) | rowMeans() | ~5× |
| OPT-6 | Cohort precomputation | ns() per year | Single ns() call | ~30× |

### 16.2 OPT-1: Fixed clearCompiled Order

**Bug:** Original code called `clearCompiled(nimble_model)` AFTER `rm(nimble_model)`, so the DLL was never released.

**Fix:**
```r
model_ref <- nimble_model                    # Save reference
rm(compiled_model, compiled_mcmc, ...)       # Remove other objects
nimble::clearCompiled(model_ref)             # Now works
rm(nimble_model, model_ref)
gc()
```

---

## 17. Complete Variable Reference

### 17.1 Configuration Constants

| Variable | Value | Set (Line) | Used For |
|----------|-------|------------|----------|
| `NUMBER_OF_CHAINS` | 4 | 108 | runMCMC() |
| `NUMBER_OF_ADAPT` | 1000 | 109 | Documentation only |
| `NUMBER_OF_BURN` | 5000 | 110 | runMCMC(nburnin) |
| `NUMBER_OF_ITERATIONS` | 10000 | 111 | runMCMC(niter) |
| `THINNING_INTERVAL` | 5 | 112 | configureMCMC(thin) |
| `BASE_YEAR` | 2010 | 115 | Target calculation |
| `TARGET_YEAR` | 2025 | 116 | Projection endpoint |
| `PROJECTED_YEARS` | 20 | 117 | Extends to 2045 |
| `REDUCTION_PERCENTAGE` | 30 | 118 | Target: 30% reduction |
| `RANDOM_SEED` | 42 | 119 | set.seed() |
| `TRANSITION_START` | 65 | 122 | Sigmoid center |
| `TRANSITION_WIDTH` | 3 | 123 | Sigmoid slope |
| `MAX_AGE_SPLINE` | 80 | 124 | Spline boundary |
| `MANUAL_TARGET_PREVALENCE` | 4.0 | 127 | Endgame threshold |
| `MANUAL_TARGET_EVAL_YEAR` | 2040 | 129 | Endgame year |

### 17.2 Derived Constants

| Variable | Formula | Computed At | Purpose |
|----------|---------|-------------|---------|
| `AGE_LINEAR_CENTER_CONSTANT` | `mean(Age_Linear)` | Section 4 | Centering linear age |
| `COHORT_CENTER_CONSTANT` | `median(Birth_Cohort)` | Section 4 | Centering cohort |
| `empirical_mean_cig` | `mean(logit(prevalence))` | Section 8.5 | Prior anchor |
| `age_knots` | `c(25, 45, 65)` | Section 4 | Spline knots |
| `cohort_knots` | `c(1945, 1965, 1985)` | Section 4 | Spline knots |

### 17.3 Data Vectors (per observation)

| Variable | Length | Range | Description |
|----------|--------|-------|-------------|
| `Prevalence` | N | (-∞, +∞) | Observed logit prevalence |
| `Def_Code_Binary` | N | {0, 1} | 0=daily, 1=current |
| `Type_Cig` | N | {0, 1} | Cigarette indicator |
| `Type_Smoked` | N | {0, 1} | Smoked tobacco indicator |
| `Type_Any` | N | {0, 1} | Any tobacco indicator |
| `weight` | N | (0, ∞) | Observation weight |
| `spline_weight_var` | N | [0, 1] | Sigmoid weight for spline |
| `linear_weight_var` | N | [0, 1] | Sigmoid weight for linear |
| `age_linear_smooth` | N | centered | Linear age term |

### 17.4 Index Vectors

| Variable | Length | Range | Description |
|----------|--------|-------|-------------|
| `Country` | N | 1:nCountry | Country index per observation |
| `Country_Region` | nCountry | 1:nRegion | Region per country (NOT length N!) |
| `Survey` | N | 1:nSurvey | Survey index per observation |

### 17.5 Spline Matrices

| Variable | Dimensions | Description |
|----------|------------|-------------|
| `age_spline_matrix` | N × 4 | Age basis evaluated at observations |
| `cohort_spline_matrix` | N × 4 | Cohort basis evaluated at observations |
| `age_cohort_interaction_matrix` | N × 16 | Tensor product |

### 17.6 Model Parameters (Stochastic Nodes)

| Parameter | Dimensions | Prior | Monitored |
|-----------|------------|-------|-----------|
| `cig_global_intercept` | 1 | N(emp_mean, 0.3²) | Yes |
| `cig_region_intercept` | 14 | N(0, between_sd²) | Yes |
| `cig_country_intercept` | ~180 | N(0, within_sd²) | Yes |
| `cig_age_spline` | ~180 × 4 | N(region_mean, region_sd²) | Yes |
| `cig_cohort_spline` | ~180 × 4 | N(region_mean, region_sd²) | Yes |
| `cig_age_cohort_interaction` | 16 | N(0, 0.2²) | Yes |
| `cig_age_linear_smooth_effect` | 1 | TruncN(-0.02, 0.05², <-0.001) | Yes |
| `smkextra_global_intercept` | 1 | N(-3.0, 0.3²) | Yes |
| `smkextra_age_spline_region_mean` | 14 × 4 | N(global, between_sd²) | Yes |
| `survey_intercept` | ~100s | N(0, survey_sd) | No |
| `residual_sd` | 1 | LogN(log(0.7), 0.5) | No |

---

## 18. Common Errors and Solutions

### 18.1 NIMBLE Errors

| Error Message | Cause | Solution |
|---------------|-------|----------|
| "maximal number of DLLs reached" | R's DLL limit (100) exceeded | `Sys.setenv(R_MAX_NUM_DLLS = 600)` before loading NIMBLE |
| "node X is not stochastic" | Missing prior for parameter | Add prior statement in model code |
| "initial log probability is -Inf" | Invalid initial values | Check inits are in valid range; verify data bounds |
| "model building failed" | Syntax error in nimbleCode | Check for missing commas, brackets |
| "could not find function X" | NIMBLE not loaded | `library(nimble)` |

### 18.2 Data Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| Zero countries processed | File path wrong | Check `data/` directory exists |
| NaN in predictions | Log of zero or negative | Check prevalence bounds (0.001, 0.999) |
| Negative prevalence | Transformation error | Verify logit/plogis direction |
| Country mismatch | Case sensitivity | Ensure all country codes lowercase |

### 18.3 Memory Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| R crashes during country loop | DLL memory leak | Use correct cleanup order (OPT-1) |
| Out of memory | Large MCMC samples | Increase RAM or reduce samples |
| Slow country loop | Accumulating DLLs | Ensure `clearCompiled()` succeeds |

---

## 19. Troubleshooting Guide

### 19.1 High R-hat on Intercepts

**Problem:** Global or regional intercepts have R-hat > 1.5

**Diagnosis:** Non-identifiability between hierarchical levels

**Solution:**
1. Use data-informed prior on global intercept: `cig_global_intercept ~ N(emp_mean, 0.3²)`
2. Verify `empirical_mean_cig` is computed correctly
3. Check that regional intercepts have zero-centered prior

### 19.2 Low ESS Everywhere

**Problem:** ESS < 100 for most parameters

**Diagnosis:** Poor MCMC mixing; high autocorrelation

**Solution:**
1. Increase `NUMBER_OF_ITERATIONS`
2. Reduce `THINNING_INTERVAL` (counterintuitive but preserves samples)
3. Check initial values are reasonable
4. Consider reparameterization (non-centered)

### 19.3 Ordering Violations

**Problem:** `check_ordering_per_draw()` returns warnings

**Diagnosis:** Stick-breaking producing p_cig > p_smoked

**Solution:**
1. This should not happen mathematically—check implementation
2. Verify `ilogit`/`plogis` are applied correctly
3. Check for NaN values in linear predictors

### 19.4 Countries Missing from Output

**Problem:** Some countries have no predictions

**Diagnosis:** Model fitting failed for those countries

**Solution:**
1. Check `warnings()` for error messages
2. Verify priors file exists: `country_priors/{gender}/{country}_*.csv`
3. Check country has sufficient data

---

## 20. Output Files Reference

### 20.1 Primary Results

| File | Description |
|------|-------------|
| `results/final_ac_predictions_nested.csv` | Global model predictions |
| `results/final_predictions_country_specific.csv` | Country model predictions |
| `results/final_weighted_results_apc_post_selection.csv` | Post-selection weighted prevalence |
| `results/target_prevalences_dual_evaluation.csv` | Target calculations |

### 20.2 Diagnostics

| File | Description |
|------|-------------|
| `diagnostics/logs/global_{gender}_diagnostics.log` | Human-readable convergence log |
| `diagnostics/tables/rhat_master.csv` | All R-hat values |
| `diagnostics/tables/ess_master.csv` | All ESS values |
| `diagnostics/tables/convergence_summary.csv` | Per-model summary |

### 20.3 Country Priors

| File | Description |
|------|-------------|
| `country_priors/{gender}/{country}_regional_ac_priors_nested.csv` | Extracted priors |

### 20.4 Predictions (Raw)

| Directory | Contents |
|-----------|----------|
| `processing/{gender}/{country}/{def_type}/{year}.rds` | Global model: n_ages × n_samples matrix |
| `results/country_specific_ac_nested/{gender}/{country}/{def_type}/{year}.rds` | Country model: n_ages × n_samples matrix |

---

## 21. Glossary

| Term | Definition |
|------|------------|
| **APC** | Age-Period-Cohort model decomposing prevalence trends |
| **Birth Cohort** | Year of birth; calculated as survey_year − age |
| **Credible Interval** | Bayesian analog of confidence interval; contains true value with stated probability |
| **ESS** | Effective Sample Size; accounts for autocorrelation in MCMC samples |
| **Hierarchical Model** | Model with nested random effects (global → region → country) |
| **ilogit** | Inverse logit function; ilogit(x) = 1/(1+exp(-x)) = plogis(x) |
| **Logit** | log(p/(1-p)); maps (0,1) to (-∞, +∞) |
| **MCMC** | Markov Chain Monte Carlo; sampling algorithm for Bayesian inference |
| **Natural Spline** | Cubic spline with linear extrapolation beyond boundary knots |
| **NIMBLE** | R package for Bayesian modeling with C++ compilation |
| **Posterior** | Distribution of parameters given data |
| **Prior** | Distribution of parameters before seeing data |
| **R-hat** | Gelman-Rubin convergence diagnostic; compares within- vs between-chain variance |
| **Shrinkage** | Pulling estimates toward a common mean; stronger for sparse data |
| **Stick-Breaking** | Construction ensuring nested probabilities sum correctly |
| **Tensor Product** | Outer product of basis functions; captures interactions |

---

*This documentation provides a complete technical reference for the WHO Tobacco Control Prevalence Projection Model implemented in R/pipeline_monolith.R v2.3.2.*
