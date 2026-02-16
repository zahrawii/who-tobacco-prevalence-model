#########################################################################################
#
#                    WHO TOBACCO CONTROL PREVALENCE PROJECTION MODEL
#                      05_models_nimble.R - NIMBLE Model Definitions
#
#   Contains: Global hierarchical model, Country-specific model definitions
#   Requires: 00_config.R (for NIMBLE to be loaded)
#   Outputs: regional_hierarchical_global_ac_model_nimble,
#            regional_country_specific_ac_model_nimble
#
#########################################################################################

# ---- 5.1 Global Hierarchical Model (Stick-Breaking Construction) - NIMBLE ----
# UPDATED: Strong data-informed priors on global intercepts to fix non-identifiability
# while preserving the global mean interpretation

regional_hierarchical_global_ac_model_nimble <- nimbleCode({

  # ==================================================================
  # LIKELIHOOD LOOP
  # ==================================================================
  for (i in 1:N) {

    # ----------------------------------------------------------------
    # Head 1: Cigarettes (FULL COMPLEXITY - has strong signal)
    # ----------------------------------------------------------------
    mu_cig[i] <- cig_global_intercept +
      cig_region_intercept[Country_Region[Country[i]]] +
      cig_country_intercept[Country[i]] +
      cig_def_code_shared * Def_Code_Binary[i] +
      spline_weight_var[i] * inprod(cig_age_spline[Country[i], 1:nAgeSpline],
                                    age_spline_matrix[i, 1:nAgeSpline]) +
      linear_weight_var[i] * cig_age_linear_smooth_effect * age_linear_smooth[i] +
      inprod(cig_cohort_spline[Country[i], 1:nCohortSpline],
             cohort_spline_matrix[i, 1:nCohortSpline]) +
      inprod(cig_age_cohort_interaction[1:nAgeXCohortSplines],
             age_cohort_interaction_matrix[i, 1:nAgeXCohortSplines]) +
      survey_intercept[Survey[i]]

    # ----------------------------------------------------------------
    # Head 2: Other smoked tobacco (SIMPLIFIED - weak signal)
    # Uses: country intercept + REGIONAL splines (not country splines)
    # No age-cohort interactions (the worst convergence offenders)
    # ----------------------------------------------------------------
    mu_smkextra[i] <- smkextra_global_intercept +
      smkextra_region_intercept[Country_Region[Country[i]]] +
      smkextra_country_intercept[Country[i]] +
      0.3 * cig_def_code_shared * Def_Code_Binary[i] +
      spline_weight_var[i] * inprod(smkextra_age_spline_region_mean[Country_Region[Country[i]], 1:nAgeSpline],
                                    age_spline_matrix[i, 1:nAgeSpline]) +
      linear_weight_var[i] * smkextra_age_linear_smooth_effect * age_linear_smooth[i] +
      inprod(smkextra_cohort_spline_region_mean[Country_Region[Country[i]], 1:nCohortSpline],
             cohort_spline_matrix[i, 1:nCohortSpline])

    # ----------------------------------------------------------------
    # Head 3: Non-smoked tobacco (SIMPLIFIED - weak signal)
    # Uses: country intercept + REGIONAL splines (not country splines)
    # No age-cohort interactions
    # ----------------------------------------------------------------
    mu_anyextra[i] <- anyextra_global_intercept +
      anyextra_region_intercept[Country_Region[Country[i]]] +
      anyextra_country_intercept[Country[i]] +
      0.3 * cig_def_code_shared * Def_Code_Binary[i] +
      spline_weight_var[i] * inprod(anyextra_age_spline_region_mean[Country_Region[Country[i]], 1:nAgeSpline],
                                    age_spline_matrix[i, 1:nAgeSpline]) +
      linear_weight_var[i] * anyextra_age_linear_smooth_effect * age_linear_smooth[i] +
      inprod(anyextra_cohort_spline_region_mean[Country_Region[Country[i]], 1:nCohortSpline],
             cohort_spline_matrix[i, 1:nCohortSpline])

    # ----------------------------------------------------------------
    # STICK-BREAKING CONSTRUCTION (unchanged)
    # ----------------------------------------------------------------
    p_cig[i]    <- ilogit(mu_cig[i])
    p_anysmk[i] <- p_cig[i] + ilogit(mu_smkextra[i]) * (1 - p_cig[i])
    p_anytob[i] <- p_anysmk[i] + ilogit(mu_anyextra[i]) * (1 - p_anysmk[i])

    # Select observed probability based on product type
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

  # ==================================================================
  # GLOBAL INTERCEPTS - STRONG DATA-INFORMED PRIORS
  # KEY FIX: These tight priors anchor the global intercepts,
  # breaking the non-identifiability with regional intercepts
  # while preserving the "global mean" interpretation
  # ==================================================================

  # CIG: Tight prior centered at empirical mean (passed as constant)
  cig_global_intercept ~ dnorm(empirical_mean_cig, sd = 0.3)

  # SMKEXTRA: Small component - P(other smoked | not cig) ~ 3-5%
  smkextra_global_intercept ~ dnorm(-3.0, sd = 0.3)

  # ANYEXTRA: Very small component - P(smokeless | not smoked) ~ 1-3%
  anyextra_global_intercept ~ dnorm(-4.0, sd = 0.3)

  # ==================================================================
  # OTHER GLOBAL PRIORS (tightened for better convergence)
  # ==================================================================

  # Definition code effect
  cig_def_code_shared ~ dnorm(0.3, sd = 1.0)

  # Residual SD - slightly tighter
  residual_sd ~ dlnorm(log(0.7), sdlog = 0.5)

  # Age spline global means - TIGHTER
  for (l in 1:nAgeSpline) {
    cig_age_spline_global_mean[l] ~ dnorm(0, sd = 1.5)
    smkextra_age_spline_global_mean[l] ~ dnorm(0, sd = 0.75)
    anyextra_age_spline_global_mean[l] ~ dnorm(0, sd = 0.75)
  }

  # Cohort spline global means - TIGHTER
  for (m in 1:nCohortSpline) {
    cig_cohort_spline_global_mean[m] ~ dnorm(0, sd = 1.5)
    smkextra_cohort_spline_global_mean[m] ~ dnorm(0, sd = 0.75)
    anyextra_cohort_spline_global_mean[m] ~ dnorm(0, sd = 0.75)
  }

  # Constrained age linear effects (must decline with age) - tighter
  cig_age_linear_smooth_effect ~ T(dnorm(-0.02, sd = 0.05), -Inf, -0.001)
  smkextra_age_linear_smooth_effect ~ T(dnorm(-0.02, sd = 0.05), -Inf, -0.001)
  anyextra_age_linear_smooth_effect ~ T(dnorm(-0.02, sd = 0.05), -Inf, -0.001)

  # ==================================================================
  # REGIONAL HYPERPRIORS - MORE INFORMATIVE FOR BETTER SHRINKAGE
  # ==================================================================

  # Between-region precisions - more informative (stronger shrinkage)
  intercept_between_region_precision ~ dgamma(4, 1)
  age_spline_between_region_precision ~ dgamma(4, 1)
  cohort_spline_between_region_precision ~ dgamma(4, 1)

  # Convert precision to sd
  intercept_between_region_sd <- 1 / sqrt(intercept_between_region_precision)
  age_spline_between_region_sd <- 1 / sqrt(age_spline_between_region_precision)
  cohort_spline_between_region_sd <- 1 / sqrt(cohort_spline_between_region_precision)

  for (r in 1:nRegion) {
    # Regional intercepts - deviations from global mean (all 3 heads)
    # Now well-identified because global intercepts are anchored
    cig_region_intercept[r] ~ dnorm(0, sd = intercept_between_region_sd)
    smkextra_region_intercept[r] ~ dnorm(0, sd = intercept_between_region_sd)
    anyextra_region_intercept[r] ~ dnorm(0, sd = intercept_between_region_sd)

    # Within-region precisions (CIG ONLY) - more informative
    cig_intercept_within_region_precision[r] ~ dgamma(4, 1)
    cig_age_spline_within_region_precision[r] ~ dgamma(3, 1)
    cig_cohort_spline_within_region_precision[r] ~ dgamma(3, 1)

    # Convert to sd (CIG only)
    cig_intercept_within_region_sd[r] <- 1 / sqrt(cig_intercept_within_region_precision[r])
    cig_age_spline_within_region_sd[r] <- 1 / sqrt(cig_age_spline_within_region_precision[r])
    cig_cohort_spline_within_region_sd[r] <- 1 / sqrt(cig_cohort_spline_within_region_precision[r])

    # Regional age spline means (all 3 heads)
    for (l in 1:nAgeSpline) {
      cig_age_spline_region_mean[r, l] ~ dnorm(cig_age_spline_global_mean[l],
                                               sd = age_spline_between_region_sd)
      smkextra_age_spline_region_mean[r, l] ~ dnorm(smkextra_age_spline_global_mean[l],
                                                    sd = age_spline_between_region_sd)
      anyextra_age_spline_region_mean[r, l] ~ dnorm(anyextra_age_spline_global_mean[l],
                                                    sd = age_spline_between_region_sd)
    }

    # Regional cohort spline means (all 3 heads)
    for (m in 1:nCohortSpline) {
      cig_cohort_spline_region_mean[r, m] ~ dnorm(cig_cohort_spline_global_mean[m],
                                                  sd = cohort_spline_between_region_sd)
      smkextra_cohort_spline_region_mean[r, m] ~ dnorm(smkextra_cohort_spline_global_mean[m],
                                                       sd = cohort_spline_between_region_sd)
      anyextra_cohort_spline_region_mean[r, m] ~ dnorm(anyextra_cohort_spline_global_mean[m],
                                                       sd = cohort_spline_between_region_sd)
    }
  }

  # Shared within-region precision for smkextra/anyextra intercepts
  smkextra_intercept_within_region_precision ~ dgamma(4, 1)
  anyextra_intercept_within_region_precision ~ dgamma(4, 1)
  smkextra_intercept_within_region_sd <- 1 / sqrt(smkextra_intercept_within_region_precision)
  anyextra_intercept_within_region_sd <- 1 / sqrt(anyextra_intercept_within_region_precision)

  # ==================================================================
  # COUNTRY-LEVEL PARAMETERS
  # ==================================================================

  for (j in 1:nCountry) {
    # CIG: Full country-level parameters (intercept + splines)
    # Country intercept is deviation from regional intercept
    cig_country_intercept[j] ~ dnorm(0, sd = cig_intercept_within_region_sd[Country_Region[j]])

    for (l in 1:nAgeSpline) {
      cig_age_spline[j, l] ~ dnorm(cig_age_spline_region_mean[Country_Region[j], l],
                                   sd = cig_age_spline_within_region_sd[Country_Region[j]])
    }

    for (m in 1:nCohortSpline) {
      cig_cohort_spline[j, m] ~ dnorm(cig_cohort_spline_region_mean[Country_Region[j], m],
                                      sd = cig_cohort_spline_within_region_sd[Country_Region[j]])
    }

    # SMKEXTRA/ANYEXTRA: Only country intercepts (splines are regional)
    # Country intercept is deviation from regional intercept
    smkextra_country_intercept[j] ~ dnorm(0, sd = smkextra_intercept_within_region_sd)
    anyextra_country_intercept[j] ~ dnorm(0, sd = anyextra_intercept_within_region_sd)
  }

  # ==================================================================
  # SURVEY RANDOM EFFECTS
  # ==================================================================

  survey_intercept_precision ~ dgamma(3, 1)
  survey_intercept_sd <- 1 / sqrt(survey_intercept_precision)

  for (s in 1:nSurvey) {
    survey_intercept[s] ~ dnorm(0, sd = survey_intercept_sd)
  }

  # ==================================================================
  # AGE-COHORT INTERACTIONS (CIG ONLY - tighter prior)
  # ==================================================================

  for (k in 1:nAgeXCohortSplines) {
    cig_age_cohort_interaction[k] ~ dnorm(0, sd = 0.2)
  }
})


# ---- 5.2 Country-Specific Model (With Informative Priors) - NIMBLE ----
# CORRECTED v2.3.1:
# - Added truncation to age linear effects (must be negative)
# - Removed smkextra/anyextra age-cohort interactions (not in global model)

regional_country_specific_ac_model_nimble <- nimbleCode({

  for (i in 1:N) {
    # Head 1: Cigarettes (FULL COMPLEXITY)
    mu_cig[i] <- cig_intercept +
      cig_def_code_shared * Def_Code_Binary[i] +
      spline_weight_var[i] * inprod(cig_age_spline[1:nAgeSpline],
                                    age_spline_matrix[i, 1:nAgeSpline]) +
      linear_weight_var[i] * cig_age_linear_smooth_effect * age_linear_smooth[i] +
      inprod(cig_cohort_spline[1:nCohortSpline],
             cohort_spline_matrix[i, 1:nCohortSpline]) +
      inprod(cig_age_cohort_interaction[1:nAgeXCohortSplines],
             age_cohort_interaction_matrix[i, 1:nAgeXCohortSplines]) +
      survey_intercept[Survey[i]]

    # Head 2: Other smoked (SIMPLIFIED - matches global model structure)
    # NO age-cohort interactions (global model doesn't have them for smkextra)
    mu_smkextra[i] <- smkextra_intercept +
      0.3 * cig_def_code_shared * Def_Code_Binary[i] +
      spline_weight_var[i] * inprod(smkextra_age_spline[1:nAgeSpline],
                                    age_spline_matrix[i, 1:nAgeSpline]) +
      linear_weight_var[i] * smkextra_age_linear_smooth_effect * age_linear_smooth[i] +
      inprod(smkextra_cohort_spline[1:nCohortSpline],
             cohort_spline_matrix[i, 1:nCohortSpline])

    # Head 3: Non-smoked (SIMPLIFIED - matches global model structure)
    # NO age-cohort interactions (global model doesn't have them for anyextra)
    mu_anyextra[i] <- anyextra_intercept +
      0.3 * cig_def_code_shared * Def_Code_Binary[i] +
      spline_weight_var[i] * inprod(anyextra_age_spline[1:nAgeSpline],
                                    age_spline_matrix[i, 1:nAgeSpline]) +
      linear_weight_var[i] * anyextra_age_linear_smooth_effect * age_linear_smooth[i] +
      inprod(anyextra_cohort_spline[1:nCohortSpline],
             cohort_spline_matrix[i, 1:nCohortSpline])

    # Stick-breaking construction
    p_cig[i]    <- ilogit(mu_cig[i])
    p_anysmk[i] <- p_cig[i] + ilogit(mu_smkextra[i]) * (1 - p_cig[i])
    p_anytob[i] <- p_anysmk[i] + ilogit(mu_anyextra[i]) * (1 - p_anysmk[i])

    p_obs[i] <- Type_Cig[i] * p_cig[i] +
      Type_Smoked[i] * p_anysmk[i] +
      Type_Any[i] * p_anytob[i]

    logit_p[i] <- log(p_obs[i] / (1 - p_obs[i]))
    tau[i] <- pow(residual_sd, -2) * weight[i]
    Prevalence[i] ~ dnorm(logit_p[i], tau = tau[i])
  }

  # ==================================================================
  # INFORMATIVE PRIORS (from global model posteriors)
  # Prior means and PRECISIONS come from constants (Section 11)
  # ==================================================================

  # Definition code effect (shared across heads)
  cig_def_code_shared ~ dnorm(cig_def_code_shared_prior_mean,
                              tau = cig_def_code_shared_prior_prec)

  # ---- CIGARETTES HEAD ----
  cig_intercept ~ dnorm(cig_intercept_prior_mean, tau = cig_intercept_prior_prec)

  # CORRECTED: Add truncation to ensure negative (prevalence declines at old ages)
  cig_age_linear_smooth_effect ~ T(dnorm(cig_age_linear_smooth_effect_prior_mean,
                                         tau = cig_age_linear_smooth_effect_prior_prec),
                                   -Inf, -0.001)

  for (l in 1:nAgeSpline) {
    cig_age_spline[l] ~ dnorm(cig_age_spline_prior_means[l],
                              tau = cig_age_spline_prior_precs[l])
  }
  for (m in 1:nCohortSpline) {
    cig_cohort_spline[m] ~ dnorm(cig_cohort_spline_prior_means[m],
                                 tau = cig_cohort_spline_prior_precs[m])
  }
  for (k in 1:nAgeXCohortSplines) {
    cig_age_cohort_interaction[k] ~ dnorm(cig_age_cohort_prior_means[k],
                                          tau = cig_age_cohort_prior_precs[k])
  }

  # ---- SMKEXTRA HEAD ----
  smkextra_intercept ~ dnorm(smkextra_intercept_prior_mean, tau = smkextra_intercept_prior_prec)

  # CORRECTED: Add truncation
  smkextra_age_linear_smooth_effect ~ T(dnorm(smkextra_age_linear_smooth_effect_prior_mean,
                                              tau = smkextra_age_linear_smooth_effect_prior_prec),
                                        -Inf, -0.001)

  for (l in 1:nAgeSpline) {
    smkextra_age_spline[l] ~ dnorm(smkextra_age_spline_prior_means[l],
                                   tau = smkextra_age_spline_prior_precs[l])
  }
  for (m in 1:nCohortSpline) {
    smkextra_cohort_spline[m] ~ dnorm(smkextra_cohort_spline_prior_means[m],
                                      tau = smkextra_cohort_spline_prior_precs[m])
  }
  # REMOVED: smkextra_age_cohort_interaction (not in global model)

  # ---- ANYEXTRA HEAD ----
  anyextra_intercept ~ dnorm(anyextra_intercept_prior_mean, tau = anyextra_intercept_prior_prec)

  # CORRECTED: Add truncation
  anyextra_age_linear_smooth_effect ~ T(dnorm(anyextra_age_linear_smooth_effect_prior_mean,
                                              tau = anyextra_age_linear_smooth_effect_prior_prec),
                                        -Inf, -0.001)

  for (l in 1:nAgeSpline) {
    anyextra_age_spline[l] ~ dnorm(anyextra_age_spline_prior_means[l],
                                   tau = anyextra_age_spline_prior_precs[l])
  }
  for (m in 1:nCohortSpline) {
    anyextra_cohort_spline[m] ~ dnorm(anyextra_cohort_spline_prior_means[m],
                                      tau = anyextra_cohort_spline_prior_precs[m])
  }
  # REMOVED: anyextra_age_cohort_interaction (not in global model)

  # ---- SURVEY RANDOM EFFECTS ----
  survey_intercept_precision ~ dgamma(3, 1)  # Match global model
  survey_intercept_sd <- 1 / sqrt(survey_intercept_precision)

  for (s in 1:nSurvey) {
    survey_intercept[s] ~ dnorm(0, sd = survey_intercept_sd)
  }

  # ---- RESIDUAL SD ----
  residual_sd ~ dlnorm(log(0.7), sdlog = 1)
})

cat("NIMBLE model definitions loaded.\n")
cat("  - regional_hierarchical_global_ac_model_nimble (Global Model)\n")
cat("  - regional_country_specific_ac_model_nimble (Country-Specific Model)\n")
