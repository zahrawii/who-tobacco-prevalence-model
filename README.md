# WHO Tobacco Prevalence Projection Model

**Version 2.3.2** | Developed for the World Health Organization

A hierarchical Bayesian age-cohort model for projecting tobacco smoking prevalence across 191 WHO Member States, designed to evaluate national progress toward the WHO 2025 and 2030 reduction targets.

---

## Citation

If you use this software or methodology, please cite:

> Jamil, H., & Gilmour, S. (2025). Global Estimation of Tobacco Use Indicators from 1990 to 2040 Using Hierarchical Bayesian Modeling: Assessing Progress Toward WHO Targets. *[Journal]*, *[Volume]*, [Pages].

This work builds upon and updates previous global tobacco estimation methodology:

> Bilano, V., Gilmour, S., Moffiet, T., et al. (2015). Global trends and projections for tobacco use, 1990–2025: an analysis of smoking indicators from the WHO Comprehensive Information Systems for Tobacco Control. *The Lancet*, 385(9972), 966–76. https://doi.org/10.1016/S0140-6736(15)60264-1

---

## Quick Start

```r
# Set working directory to project root
setwd("path/to/who_modeling_r_code")

# Option 1: Run complete pipeline (8-12 hours)
source("R/pipeline_monolith.R")

# Option 2: Run test with 4 countries (~20 minutes)
source("R/test_pipeline.R")
run_test_pipeline()

# Option 3: Run modular pipeline phase-by-phase
source("R/_main.R")
run_modular_pipeline("prep")     # Phase 1: Data preparation
run_modular_pipeline("global")   # Phase 2: Global model
run_modular_pipeline("country")  # Phase 3: Country models
run_modular_pipeline("eval")     # Phase 4: Evaluation
run_modular_pipeline("agg")      # Phase 5: Aggregation
run_modular_pipeline("viz")      # Phase 6: Visualization
```

---

## Background and Motivation

Tobacco use remains the leading cause of preventable morbidity and mortality worldwide, responsible for over 8 million deaths annually (WHO, 2023). Despite decades of global tobacco control efforts—including the WHO Framework Convention on Tobacco Control (FCTC) and the MPOWER measures—tobacco consumption persists as a significant public health challenge, particularly in low- and middle-income countries where the majority of the world's 1.3 billion tobacco users reside.

In 2013, WHO Member States committed to achieving a **30% relative reduction** in tobacco prevalence by 2025, using 2010 as the baseline year. Tracking progress toward this target requires statistical models capable of handling several simultaneous challenges:

| Challenge | Description |
|-----------|-------------|
| **Sparse data** | Many countries have conducted only one or two national surveys; 6 countries have no survey data at all |
| **Inconsistent methodology** | Surveys vary in definitions (daily vs. current use), age groupings, and sampling frames |
| **Age structure** | Prevalence varies systematically across the lifespan—rising sharply in youth, peaking in middle age, and declining with mortality selection |
| **Generational shifts** | Birth cohort effects are substantial: someone born in 1945 has a fundamentally different smoking trajectory than someone born in 1985 |
| **Nested products** | Cigarette smoking is a subset of all smoked tobacco, which is itself a subset of any tobacco use—these logical constraints must be preserved |

The last comprehensive global tobacco estimation study (Bilano et al., 2015) is now a decade old and does not reflect recent trends, new survey data, or advances in Bayesian computation. This pipeline provides an updated analysis using data from **548 surveys** across **191 countries** compiled by the WHO Comprehensive Information System for Tobacco Control.

---

## Model Overview

### Hierarchical Structure

The model estimates tobacco prevalence within a four-level hierarchy that enables borrowing of strength across countries while allowing local deviation where data warrant:

```
Global (pooled across all data)
    │
    ├── Region (14 WHO regions)
    │       │
    │       └── Country (191 Member States)
    │               │
    │               └── Survey (548 surveys, absorbs measurement error)
```

Countries with extensive survey histories can deviate substantially from global patterns; countries with sparse data are partially pooled toward regional and global estimates, with uncertainty appropriately widened.

### Tobacco Use Indicators

The model jointly estimates three nested tobacco product categories:

| Indicator | Definition |
|-----------|------------|
| **Current cigarette smoking** | Smoked manufactured or hand-rolled cigarettes at least once in the past 30 days |
| **Current use of any smoked tobacco** | Used any smoked tobacco product (cigarettes, cigars, pipes, waterpipe, etc.) in the past 30 days |
| **Current use of any tobacco product** | Used any tobacco product (smoked or smokeless) in the past 30 days |

By definition, these indicators are nested: $P(\text{cigarettes}) \leq P(\text{any smoked}) \leq P(\text{any tobacco})$. The model enforces this constraint for every posterior draw using a stick-breaking construction.

---

## Mathematical Specification

### Likelihood

Survey observations are modeled on the logit scale. For observation $i$ from survey $k$ in country $j$ within region $r$:

$$
y_i \sim \text{Normal}(\mu_i, \sigma^2 / w_i)
$$

where $y_i = \text{logit}(\text{prevalence}_i)$, $w_i$ is a precision weight based on age-group width (narrower age bands receive higher weight), and $\mu_i$ is the linear predictor.

### Linear Predictor (Cigarettes)

The cigarette smoking component receives the most complex specification, as signal is strongest for this indicator:

$$
\mu_{\text{cig},i} = \underbrace{\alpha^{(G)}}_{\text{global}} + \underbrace{\alpha^{(R)}_r}_{\text{region}} + \underbrace{\alpha^{(C)}_j}_{\text{country}} + \underbrace{\beta \cdot D_i}_{\text{definition}} + \underbrace{f_{\text{age}}(a_i; \boldsymbol{\gamma}_j)}_{\text{age effect}} + \underbrace{f_{\text{cohort}}(c_i; \boldsymbol{\phi}_j)}_{\text{cohort effect}} + \underbrace{\mathbf{z}_i^\top \boldsymbol{\psi}}_{\text{age-cohort interaction}} + \underbrace{\eta_{k}}_{\text{survey}}
$$

| Term | Description |
|------|-------------|
| $\alpha^{(G)}$ | Global intercept (anchored at empirical logit-mean prevalence) |
| $\alpha^{(R)}_r$ | Regional deviation from global mean |
| $\alpha^{(C)}_j$ | Country deviation from regional mean |
| $\beta \cdot D_i$ | Adjustment for definition type (daily = 0, current = 1) |
| $f_{\text{age}}$ | Country-specific natural spline for age effects |
| $f_{\text{cohort}}$ | Country-specific natural spline for birth cohort effects |
| $\mathbf{z}_i^\top \boldsymbol{\psi}$ | Tensor product interaction (16 parameters: 4 age × 4 cohort basis functions) |
| $\eta_k$ | Survey-specific random effect (absorbs systematic survey bias) |

### Stick-Breaking Construction

To guarantee $P(\text{cig}) \leq P(\text{smoked}) \leq P(\text{any})$ for every posterior draw without post-hoc truncation:

$$
\begin{aligned}
P(\text{cigarettes}) &= \text{logistic}(\mu_{\text{cig}}) \\[0.5em]
P(\text{any smoked}) &= P(\text{cig}) + \text{logistic}(\mu_{\text{smkextra}}) \cdot (1 - P(\text{cig}}) \\[0.5em]
P(\text{any tobacco}) &= P(\text{smoked}) + \text{logistic}(\mu_{\text{anyextra}}) \cdot (1 - P(\text{smoked}})
\end{aligned}
$$

The "extra" components model the conditional probability of using additional tobacco products given non-use of the nested category. For example, $\text{logistic}(\mu_{\text{smkextra}})$ represents the probability of using other smoked tobacco (pipes, cigars, waterpipe) among those who do not smoke cigarettes.

**Numerical example:** If $P(\text{cigarettes}) = 0.20$ and $\text{logistic}(\mu_{\text{smkextra}}) = 0.05$, then:
$$P(\text{any smoked}) = 0.20 + 0.05 \times 0.80 = 0.24$$

### Smooth Elderly Transition

Survey data for ages beyond 70 are sparse, and survivorship bias complicates interpretation (smokers die younger, artificially deflating observed elderly prevalence). Rather than extrapolate flexible splines into data-sparse regions—which can produce implausible oscillations—the model transitions smoothly from spline-based estimation to constrained linear extrapolation:

$$
w(a) = \frac{1}{1 + \exp\left(\frac{a - 65}{3}\right)}
$$

$$
f_{\text{age}}(a) = w(a) \cdot \text{spline}(a) + (1 - w(a)) \cdot \delta \cdot (a - 65)
$$

The linear coefficient $\delta$ is constrained to be negative (truncated normal prior), enforcing the epidemiologically necessary property that prevalence declines with age beyond 65.

---

## Spline Configuration

Knot placement reflects established inflection points in smoking epidemiology, not arbitrary choices.

### Age Spline Knots: 25, 45, 65

| Knot | Epidemiological Rationale |
|------|---------------------------|
| **25** | End of initiation window. Over 80% of ever-smokers begin before age 25; initiation after this age is rare (USDHHS, 2012). |
| **45** | Onset of net decline. Cessation rates remain low through middle age—the "hardening" hypothesis—then begin to rise as health consequences manifest. |
| **65** | Mortality selection threshold. Differential survival accelerates apparent prevalence decline as the smoking-attributable mortality burden increases. |

Boundary knots: 15 (minimum age of policy interest) and 80 (transition to linear extrapolation).

### Cohort Spline Knots: 1945, 1965, 1985

| Knot | Epidemiological Rationale |
|------|---------------------------|
| **1945** | Post-war cohort peak. In many countries, cohorts born 1940–1950 exhibited the highest lifetime smoking prevalence, shaped by post-war economic conditions and pre-regulation tobacco marketing. |
| **1965** | Post–Surgeon General generation. The first cohort born after the 1964 U.S. Surgeon General's Report entered adulthood amid growing awareness of health risks and early advertising restrictions. |
| **1985** | Modern tobacco control era. Cohorts reaching adulthood after 2000 faced comprehensive MPOWER policies: taxes, advertising bans, smoke-free legislation, and graphic warning labels. |

Boundary knots are set at the observed cohort range in the data (approximately 1910–2005).

---

## Prior Specification

All priors are proper and weakly informative. The global cigarette intercept uses a data-informed prior (empirical Bayes anchoring) to resolve non-identifiability between hierarchical intercept levels.

| Parameter | Prior | Rationale |
|-----------|-------|-----------|
| $\alpha^{(G)}_{\text{cig}}$ | $\mathcal{N}(\bar{y}_{\text{logit}}, 0.3^2)$ | Anchored at empirical logit-mean; breaks intercept non-identifiability |
| $\alpha^{(G)}_{\text{smkextra}}$ | $\mathcal{N}(-3.0, 0.3^2)$ | Implies ~5% baseline for non-cigarette smoked tobacco |
| $\alpha^{(G)}_{\text{anyextra}}$ | $\mathcal{N}(-4.0, 0.3^2)$ | Implies ~2% baseline for smokeless tobacco |
| Regional/Country SD | $\text{Gamma}(4, 1) \to \sigma$ | Weakly informative; allows substantial heterogeneity |
| Age spline coefficients | $\mathcal{N}(0, 1)$ | Weakly informative; regularizes toward flat age profile |
| Cohort spline coefficients | $\mathcal{N}(0, 1)$ | Weakly informative; regularizes toward flat cohort profile |
| Linear elderly slope $\delta$ | $\mathcal{N}_{(-\infty, -0.001)}(-0.02, 0.05^2)$ | **Truncated** to enforce declining prevalence beyond age 65 |
| Survey effects | $\mathcal{N}(0, \sigma^2_{\text{survey}})$ | Absorbs survey-level systematic bias |
| Survey SD | $\text{Gamma}(2, 1) \to \sigma_{\text{survey}}$ | Weakly informative |

The truncation on the linear elderly slope deserves emphasis: without this constraint, posterior draws occasionally predict rising prevalence among the very old—biologically implausible given mortality selection.

---

## Two-Stage Estimation

The pipeline fits two complementary models:

### Stage 1: Global Hierarchical Model

A single model estimates all parameters using data from all 191 countries simultaneously. This maximizes borrowing of strength across countries but presents computational challenges (thousands of parameters, complex posterior geometry).

### Stage 2: Country-Specific Models

For each country with survey data, a simpler model is fitted using **informative priors** derived from the global model's posterior. This allows individual countries to deviate from global patterns when their data warrant it while remaining anchored to the global fit.

#### Prior Extraction

After fitting the global model, posterior distributions for each country's parameters are summarized (mean and SD). These become priors for the country-specific model, with a shrinkage factor applied:

$$
\sigma_{\text{prior}}^{(\text{country})} = 0.7 \times \sigma_{\text{posterior}}^{(\text{global})}
$$

The 0.7 shrinkage factor was selected based on preliminary analyses showing it produced stable country-specific fits without excessive deviation from the global pattern. Results are robust to moderate variation in this factor (0.5–0.9).

### Model Selection

For each country-sex stratum, both models generate predictions for observed survey years. The model with lower root mean squared error (RMSE) against observed data is selected for final predictions. This pragmatic approach avoids the computational burden of fully Bayesian model averaging while ensuring the best-fitting model is used for each country.

---

## Countries Without Survey Data

Six countries lack any tobacco survey data in the WHO database:

| Country | Region | Notes |
|---------|--------|-------|
| Taiwan | Eastern Asia | Uses regional hierarchical prior |
| Hong Kong SAR | Eastern Asia | Uses regional hierarchical prior |
| Macao SAR | Eastern Asia | Uses regional hierarchical prior |
| Kosovo | Southern Europe | Uses regional hierarchical prior |
| South Sudan | Sub-Saharan Africa | Uses regional hierarchical prior |
| Somalia | Sub-Saharan Africa | Uses regional hierarchical prior |

For these countries, predictions are generated by sampling from the regional hierarchical prior—statistically valid since, in the absence of likelihood contributions, the posterior equals the prior. These predictions have appropriately wider credible intervals reflecting the greater uncertainty.

**Important:** Results for these countries are flagged with `Has_Survey_Data = FALSE` in output files and should be interpreted with appropriate caution in any publication.

---

## Target Evaluation

### WHO 30% Reduction Target (2025)

The primary output is the posterior probability that a country achieves a 30% relative reduction in adult tobacco prevalence between 2010 (baseline) and 2025 (target year):

$$
P\left( \frac{\text{prev}_{2025}}{\text{prev}_{2010}} < 0.70 \right)
$$

This probability is computed across all posterior draws, preserving full uncertainty propagation through population-weighted age aggregation.

### Extended Targets

The pipeline also evaluates:

- **2030 Target**: 30% reduction by 2030 (for countries unlikely to meet 2025)
- **Tobacco Endgame**: Probability of reaching <5% prevalence by 2040

---

## Installation

### System Requirements

| Requirement | Specification |
|-------------|---------------|
| R version | 4.3.0 or later |
| NIMBLE | 1.0.0 or later (requires C++ compiler) |
| RAM | 16 GB minimum; 32 GB recommended |
| CPU | 8+ cores recommended for parallel chains |
| Disk space | ~5 GB for intermediate files |
| Operating system | Tested on Windows 10/11, Ubuntu 22.04, macOS 13+ |

### C++ Compiler Setup

NIMBLE requires a working C++ compiler:

- **Windows**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
- **macOS**: Install Xcode Command Line Tools (`xcode-select --install`)
- **Linux**: Install `build-essential` (`sudo apt install build-essential`)

### R Package Dependencies

```r
install.packages(c(
  # Core modeling
  "nimble", "coda", "splines", "Matrix",

  # Data manipulation
 "tidyverse", "readxl", "janitor",

  # Parallel processing
  "foreach", "doParallel", "parallel",

  # Visualization
  "ggplot2", "scales", "viridis", "patchwork",
  "sf", "rnaturalearth", "rnaturalearthdata",

  # Publication tables
  "gt", "knitr",

  # Optional: completion notification
  "beepr"
))
```

### Verify Installation

```r
# Test NIMBLE compilation
library(nimble)
nimbleOptions(verbose = FALSE)
testModel <- nimbleModel(
  nimbleCode({ y ~ dnorm(0, 1) }),
  data = list(y = 0)
)
# If no errors, NIMBLE is correctly configured
```

---

## Usage

### Running the Full Pipeline

```r
# Ensure working directory is project root
setwd("path/to/who_modeling_r_code")

# Run complete pipeline
source("R/pipeline_monolith.R")
```

**Expected console output:**
```
###########################################################################
#
#              WHO TOBACCO CONTROL PREVALENCE PROJECTION
#                         VERSION 2.3.2
#
###########################################################################

Loading required packages...
  Sourcing R/prediction_all_countries.R...
  Sourcing R/mcmc_diagnostics.R...

========== LOADING DATA ==========
  Reading: data/Prevalence RESHAPE 28 Sep 2024.xlsx
  Observations: 45,623
  Countries: 185
  Year range: 1990-2022
...
```

### Configuration

Key parameters at the top of `R/pipeline_monolith.R`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `NUMBER_OF_CHAINS` | 4 | MCMC chains (minimum 3 recommended) |
| `NUMBER_OF_BURN` | 5,000 | Burn-in iterations per chain |
| `NUMBER_OF_ITERATIONS` | 10,000 | Post-burn-in samples per chain |
| `THINNING_INTERVAL` | 5 | Keep every Nth sample |
| `BASE_YEAR` | 2010 | WHO baseline year |
| `TARGET_YEAR` | 2025 | Primary target year |
| `ENDGAME_YEAR` | 2040 | Long-term projection horizon |

### Modular Execution

For debugging, development, or running on limited computational resources:

```r
source("R/_main.R")

# Run phases sequentially (checkpoints saved between phases)
run_modular_pipeline("prep")     # ~2 minutes: load data, build splines
run_modular_pipeline("global")   # ~4 hours: fit global hierarchical model
run_modular_pipeline("country")  # ~3 hours: fit 185 country-specific models
run_modular_pipeline("eval")     # ~10 minutes: compute RMSE, select models
run_modular_pipeline("agg")      # ~20 minutes: regional aggregation
run_modular_pipeline("viz")      # ~30 minutes: generate figures

# Resume from checkpoint after interruption
load("checkpoints/phase2_global.RData")
run_modular_pipeline("country")
```

### Expected Runtime

| Component | Time (16-core machine) | Time (4-core machine) |
|-----------|------------------------|----------------------|
| Data preparation | 2 minutes | 2 minutes |
| Global model (males) | 1.5–2 hours | 4–6 hours |
| Global model (females) | 1.5–2 hours | 4–6 hours |
| Country models (all) | 2–3 hours | 6–8 hours |
| Evaluation & aggregation | 30 minutes | 1 hour |
| Visualization | 30 minutes | 1 hour |
| **Total** | **8–12 hours** | **20–30 hours** |

---

## Data Requirements

### Primary Input: Survey Data

An Excel file containing WHO tobacco prevalence estimates from national surveys.

**Required columns:**

| Column | Description | Example |
|--------|-------------|---------|
| `wb_country_abv` | ISO 3166-1 alpha-3 country code (lowercase) | `usa`, `gbr`, `ind` |
| `survey_year` | Year survey was conducted | `2018` |
| `sex_name` | Sex category | `male`, `female` |
| `start_age`, `end_age` | Age band boundaries | `25`, `34` |
| `prevalence` | Point estimate (0–100 scale) | `23.5` |
| `def_code` | Definition code | `cd.0101` (daily), `cd.0112` (current) |
| `type_code` | Product type code | `tt.003` (cigarettes), `tt.002` (smoked), `tt.001` (any) |

**Data source:** WHO Comprehensive Information System for Tobacco Control

### Secondary Input: Population Weights

A CSV file providing age-sex-year population estimates for computing weighted national prevalence.

```
data/weights_15_2022.csv
```

---

## Output Structure

```
results/
├── final_ac_predictions_nested.csv              # Global model predictions
├── final_predictions_country_specific.csv       # Country model predictions
├── final_weighted_results_apc_post_selection.csv # Post-selection combined results
├── target_prevalences_dual_evaluation.csv       # Target achievement calculations
└── regional_aggregation/
    └── regional_full_summary_lancet.csv         # Regional summary statistics

diagnostics/
├── logs/                           # Human-readable convergence logs per model
└── tables/
    ├── rhat_master.csv             # All R-hat values across all parameters
    ├── ess_master.csv              # All ESS values across all parameters
    └── convergence_summary.csv     # Summary: % converged, worst parameters

evaluation/
├── model_evaluation_by_country_sex.csv    # RMSE by country-sex stratum
├── model_evaluation_summary.csv           # Model selection summary
└── model_selection_summary.csv            # Which model selected per stratum

outputs/
├── figures/
│   ├── publication/               # Publication-ready figures
│   ├── maps/                      # Choropleth maps
│   ├── trends/                    # Country trend plots
│   └── validation/                # Observed vs. predicted plots
└── tables/
    └── publication/               # GT publication tables
```

### Key Output Columns

| Column | Description |
|--------|-------------|
| `Country` | ISO 3166-1 alpha-3 code |
| `Year` | Calendar year |
| `Sex` | `males` or `females` |
| `Age_Midpoint` | Midpoint of 5-year age band |
| `Prevalence` | Posterior mean (probability scale, 0–1) |
| `lower_ci`, `upper_ci` | 95% credible interval bounds |
| `prob_achieving_target` | $P(\text{prev}_{2025} < 0.70 \times \text{prev}_{2010})$ |
| `Selected_Model` | `Global` or `Country` (whichever had lower RMSE) |
| `Has_Survey_Data` | `TRUE` if country has survey data; `FALSE` for prior-driven predictions |

---

## Convergence Diagnostics

### Gelman-Rubin Statistic (R-hat)

| R-hat | Interpretation | Action |
|-------|----------------|--------|
| < 1.05 | Excellent convergence | Results reliable |
| 1.05–1.10 | Acceptable | Interpret with minor caution |
| 1.10–1.50 | Concerning | Increase iterations; check model specification |
| > 1.50 | Poor convergence | Do not trust results; reparameterize |

### Effective Sample Size (ESS)

| ESS | Interpretation | Action |
|-----|----------------|--------|
| > 400 | Reliable posterior summaries | Sufficient for publication |
| 100–400 | Adequate for point estimates | CIs may be rough; consider more iterations |
| < 100 | Insufficient | Increase iterations or thin less aggressively |

### Troubleshooting Poor Convergence

If R-hat > 1.10 for key parameters:

1. **Increase burn-in**: Set `NUMBER_OF_BURN = 10000`
2. **Increase iterations**: Set `NUMBER_OF_ITERATIONS = 20000`
3. **Check data**: Examine the flagged country for data anomalies
4. **Simplify model**: Consider removing age-cohort interactions for problematic countries
5. **Consult diagnostics**: See `diagnostics/logs/` for detailed trace plots

---

## Validation

### Internal Validation

- **RMSE-based model selection**: For each country-sex stratum, the model (global or country-specific) with lower RMSE against observed survey data is selected
- **Posterior predictive checks**: Model fit assessed by comparing observed data to simulated data from the posterior predictive distribution

### Known Limitations

1. **No external validation cohort.** All available survey data are used for fitting. Out-of-sample prediction accuracy is assessed only through cross-validation metrics, not truly held-out data.

2. **Survey representativeness assumed.** The model assumes surveys are nationally representative. Where they are not (e.g., urban-only samples in some countries), predictions may be biased.

3. **Prior-driven predictions for data-sparse countries.** For the 6 countries without survey data, predictions are entirely prior-driven. While statistically valid, these should be interpreted with substantial caution.

4. **Smooth cohort effects.** The model assumes generational change is smooth. Abrupt policy shocks (e.g., sudden large tax increases) may induce discontinuities the spline-based model cannot capture.

5. **No explicit period effects.** The model does not estimate calendar-year effects separate from age and cohort. This is a deliberate simplification; age-period-cohort models face well-known identification problems (the APC constraint: Age + Cohort = Period).

6. **Projection uncertainty grows with horizon.** Projections to 2040 carry substantially more uncertainty than estimates for 2025. Credible intervals appropriately widen, but structural model uncertainty is not fully captured.

---

## Directory Structure

```
who_modeling_r_code/
├── R/
│   ├── pipeline_monolith.R          # Complete pipeline (~4,800 lines)
│   ├── _main.R                      # Orchestrator for modular execution
│   ├── 00_config.R                  # Configuration and constants
│   ├── 01_data_prep.R               # Data loading and cleaning
│   ├── 02_regional_mapping.R        # WHO regional classification (14 regions)
│   ├── 03_splines.R                 # Age and cohort spline construction
│   ├── 04_utils.R                   # Utility functions
│   ├── 05_models_nimble.R           # NIMBLE model definitions
│   ├── 06_run_global_model.R        # Global hierarchical model fitting
│   ├── 07_run_country_model.R       # Country-specific model fitting
│   ├── 08_evaluation.R              # RMSE calculation and model selection
│   ├── 09_aggregation.R             # Population-weighted regional aggregation
│   ├── 10_visualization.R           # Plotting functions
│   ├── 11_publication_tables.R      # Table orchestration
│   ├── 12_publication_figures.R     # Publication figure generation
│   ├── 13_strata_figures.R          # Country-level diagnostic figures
│   ├── mcmc_diagnostics.R           # Comprehensive convergence analysis
│   ├── prediction_all_countries.R   # Prediction functions (including no-data countries)
│   ├── publication_tables.R         # GT table generation code
│   └── test_pipeline.R              # Validation with minimal settings
├── data/
│   ├── Prevalence RESHAPE *.xlsx    # Primary survey data (WHO)
│   └── weights_15_2022.csv          # Population weights by age-sex-year
├── results/                         # Model outputs and predictions
├── diagnostics/                     # MCMC convergence diagnostics
├── evaluation/                      # Model comparison metrics
├── processing/                      # Intermediate MCMC samples
├── country_priors/                  # Extracted priors for country models
├── checkpoints/                     # RData files for pipeline resumption
├── outputs/                         # Figures and tables
├── README.md                        # This file
└── NIMBLEV4_DOCUMENTATION.md        # Detailed technical documentation
```

---

## Reproducibility

### Random Seeds

Random number generator seeds are set at pipeline initialization (`set.seed(12345)`) to ensure reproducibility within a given computational environment. However, due to:

- MCMC stochasticity across runs
- Platform-specific numerical precision differences
- Parallel execution order variability

**Exact replication of point estimates is not guaranteed across different systems.** Qualitative conclusions (which countries achieve targets, regional patterns) should be robust.

### Session Information

For maximum reproducibility, the pipeline logs session information at completion. To manually capture:

```r
sessionInfo()
# Save to file for publication supplementary materials
writeLines(capture.output(sessionInfo()), "session_info.txt")
```

### Tested Configuration

This pipeline was developed and tested with:

- R 4.3.2
- NIMBLE 1.0.1
- Ubuntu 22.04 LTS / Windows 11
- 32 GB RAM, 16-core CPU

---

## Acknowledgments

This work was conducted at the Graduate School of Public Health, St. Luke's International University, Tokyo, Japan.

**Data source:** World Health Organization Comprehensive Information System for Tobacco Control (WHO CiT)

**Methodological foundation:** This pipeline builds upon and substantially extends the methodology developed in Bilano et al. (2015), incorporating advances in Bayesian computation, age-cohort modeling, and uncertainty quantification.

---

## References

1. World Health Organization. (2023). *WHO report on the global tobacco epidemic 2023: Protect people from tobacco smoke.* Geneva: WHO. https://www.who.int/publications/i/item/9789240077164

2. Bilano, V., Gilmour, S., Moffiet, T., et al. (2015). Global trends and projections for tobacco use, 1990–2025: an analysis of smoking indicators from the WHO Comprehensive Information Systems for Tobacco Control. *The Lancet*, 385(9972), 966–76. https://doi.org/10.1016/S0140-6736(15)60264-1

3. GBD 2019 Tobacco Collaborators. (2021). Spatial, temporal, and demographic patterns in prevalence of smoking tobacco use and attributable disease burden in 204 countries and territories, 1990–2019. *The Lancet*, 397(10292), 2337–60.

4. World Health Organization. (2013). *Global action plan for the prevention and control of noncommunicable diseases 2013-2020.* Geneva: WHO. https://www.who.int/publications/i/item/9789241506236

5. World Health Organization. (2003). *WHO Framework Convention on Tobacco Control.* Geneva: WHO. https://fctc.who.int/who-fctc/overview

---

## Contact

For questions about the methodology or code:

**Hasan Jamil**
Graduate School of Public Health
St. Luke's International University
Tokyo, Japan

**Supervisor: Stuart Gilmour, PhD**
