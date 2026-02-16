#########################################################################################
#
#                    PREDICTION FOR ALL COUNTRIES (Including Those Without Data)
#
#   Statistical Approach:
#   - Countries WITH data: Use MCMC posterior samples of country-specific parameters
#   - Countries WITHOUT data: Sample from hierarchical prior using regional hyperparameters
#
#   This is statistically valid because:
#   - For countries without data, posterior = prior (no likelihood update)
#   - The prior hyperparameters (regional means, variances) ARE informed by data
#   - Uncertainty is properly propagated (larger for missing countries)
#
#   Hierarchical Structure (from NIMBLE model):
#   - Global → Regional → Country
#   - CIG: Full country-specific splines
#   - SMKEXTRA/ANYEXTRA: Country intercepts only, regional splines
#
#########################################################################################


#' Sample CIG parameters from hierarchical prior for a country without data
#'
#' For the cigarette head, countries have full country-specific parameters:
#' - country_intercept ~ N(0, cig_intercept_within_region_sd[region])
#' - age_spline[l] ~ N(cig_age_spline_region_mean[region,l], cig_age_spline_within_region_sd[region])
#' - cohort_spline[m] ~ N(cig_cohort_spline_region_mean[region,m], cig_cohort_spline_within_region_sd[region])
#'
#' @param region_num Region index for this country
#' @param mcmc_samples Matrix of MCMC samples (rows = iterations, cols = parameters)
#' @param n_samples Number of MCMC samples
#' @param nAgeSpline Number of age spline basis functions
#' @param nCohortSpline Number of cohort spline basis functions
#' @return List with sampled country intercept, age splines, cohort splines
sample_cig_params_from_prior <- function(region_num,
                                          mcmc_samples,
                                          n_samples,
                                          nAgeSpline,
                                          nCohortSpline) {

  # --- Country Intercept ---
  # cig_country_intercept[j] ~ dnorm(0, sd = cig_intercept_within_region_sd[Country_Region[j]])
  intercept_sd_col <- paste0("cig_intercept_within_region_sd[", region_num, "]")
  if (!intercept_sd_col %in% colnames(mcmc_samples)) {
    stop(sprintf("Required column not found: %s", intercept_sd_col))
  }
  within_region_sd <- mcmc_samples[, intercept_sd_col]
  country_intercept <- rnorm(n_samples, mean = 0, sd = within_region_sd)

  # --- Age Splines ---
  # cig_age_spline[j, l] ~ dnorm(cig_age_spline_region_mean[Country_Region[j], l],
  #                              sd = cig_age_spline_within_region_sd[Country_Region[j]])
  age_splines <- matrix(0, nrow = n_samples, ncol = nAgeSpline)

  age_sd_col <- paste0("cig_age_spline_within_region_sd[", region_num, "]")
  if (age_sd_col %in% colnames(mcmc_samples)) {
    age_within_sd <- mcmc_samples[, age_sd_col]
  } else {
    warning(sprintf("Age spline SD not found (%s), using default 0.5", age_sd_col))
    age_within_sd <- rep(0.5, n_samples)
  }

  for (l in 1:nAgeSpline) {
    region_mean_col <- paste0("cig_age_spline_region_mean[", region_num, ", ", l, "]")
    if (region_mean_col %in% colnames(mcmc_samples)) {
      region_mean <- mcmc_samples[, region_mean_col]
      age_splines[, l] <- rnorm(n_samples, mean = region_mean, sd = age_within_sd)
    } else {
      warning(sprintf("Regional mean not found: %s, using 0", region_mean_col))
      age_splines[, l] <- rnorm(n_samples, mean = 0, sd = age_within_sd * 1.5)
    }
  }

  # --- Cohort Splines ---
  # cig_cohort_spline[j, m] ~ dnorm(cig_cohort_spline_region_mean[Country_Region[j], m],
  #                                 sd = cig_cohort_spline_within_region_sd[Country_Region[j]])
  cohort_splines <- matrix(0, nrow = n_samples, ncol = nCohortSpline)

  cohort_sd_col <- paste0("cig_cohort_spline_within_region_sd[", region_num, "]")
  if (cohort_sd_col %in% colnames(mcmc_samples)) {
    cohort_within_sd <- mcmc_samples[, cohort_sd_col]
  } else {
    warning(sprintf("Cohort spline SD not found (%s), using default 0.5", cohort_sd_col))
    cohort_within_sd <- rep(0.5, n_samples)
  }

  for (m in 1:nCohortSpline) {
    region_mean_col <- paste0("cig_cohort_spline_region_mean[", region_num, ", ", m, "]")
    if (region_mean_col %in% colnames(mcmc_samples)) {
      region_mean <- mcmc_samples[, region_mean_col]
      cohort_splines[, m] <- rnorm(n_samples, mean = region_mean, sd = cohort_within_sd)
    } else {
      warning(sprintf("Regional mean not found: %s, using 0", region_mean_col))
      cohort_splines[, m] <- rnorm(n_samples, mean = 0, sd = cohort_within_sd * 1.5)
    }
  }

  return(list(
    country_intercept = country_intercept,
    age_splines = age_splines,
    cohort_splines = cohort_splines
  ))
}


#' Sample SMKEXTRA/ANYEXTRA parameters from hierarchical prior
#'
#' For smkextra/anyextra heads, countries only have intercepts (splines are regional):
#' - smkextra_country_intercept[j] ~ dnorm(0, sd = smkextra_intercept_within_region_sd)
#' - anyextra_country_intercept[j] ~ dnorm(0, sd = anyextra_intercept_within_region_sd)
#' Note: These use global (not per-region) within-region SD
#'
#' @param head_name Either "smkextra" or "anyextra"
#' @param mcmc_samples Matrix of MCMC samples
#' @param n_samples Number of MCMC samples
#' @return Sampled country intercept vector
sample_extra_intercept_from_prior <- function(head_name,
                                               mcmc_samples,
                                               n_samples) {

  # smkextra/anyextra use a single global within-region SD (not per-region)
  sd_col <- paste0(head_name, "_intercept_within_region_sd")

  if (sd_col %in% colnames(mcmc_samples)) {
    within_sd <- mcmc_samples[, sd_col]
  } else {
    warning(sprintf("%s SD not found, using default 0.5", sd_col))
    within_sd <- rep(0.5, n_samples)
  }

  country_intercept <- rnorm(n_samples, mean = 0, sd = within_sd)
  return(country_intercept)
}


#' Main prediction function for a single country
#'
#' Generates predictions for all three heads (cig, smkextra, anyextra) using
#' stick-breaking construction.
#'
#' @param country_code Country code (lowercase)
#' @param region_num Region index
#' @param has_data Logical: does this country have survey data in MCMC?
#' @param country_num_in_mcmc Country index in MCMC (if has_data = TRUE)
#' @param cig_samples List of CIG MCMC samples (extracted from combined_samples_matrix)
#' @param smkextra_samples List of SMKEXTRA MCMC samples
#' @param anyextra_samples List of ANYEXTRA MCMC samples
#' @param def_code_shared_samples Vector of def_code_shared samples
#' @param combined_samples_matrix Full MCMC samples matrix (for sampling from prior)
#' @param current_data Data frame with covariates for this country
#' @param n_samples Number of MCMC iterations
#' @param nAgeSpline Number of age spline basis
#' @param nCohortSpline Number of cohort spline basis
#' @return List with predictions_cig, predictions_smoked, predictions_any matrices
predict_single_country_all_heads <- function(country_code,
                                              region_num,
                                              has_data,
                                              country_num_in_mcmc = NULL,
                                              cig_samples,
                                              smkextra_samples,
                                              anyextra_samples,
                                              def_code_shared_samples,
                                              combined_samples_matrix,
                                              current_data,
                                              n_samples,
                                              nAgeSpline,
                                              nCohortSpline) {

  # ==================================================================
  # EXTRACT OR SAMPLE COUNTRY PARAMETERS
  # ==================================================================

  if (has_data && !is.null(country_num_in_mcmc)) {
    # Country HAS data - use fitted MCMC parameters

    # CIG parameters
    cig_country_intercept <- cig_samples$country_intercept[, paste0("cig_country_intercept[", country_num_in_mcmc, "]")]
    cig_age_spline <- cig_samples$age_spline[, paste0("cig_age_spline[", country_num_in_mcmc, ", ", 1:nAgeSpline, "]")]
    cig_cohort_spline <- cig_samples$cohort_spline[, paste0("cig_cohort_spline[", country_num_in_mcmc, ", ", 1:nCohortSpline, "]")]

    # SMKEXTRA/ANYEXTRA intercepts
    smkextra_country_intercept <- smkextra_samples$country_intercept[, paste0("smkextra_country_intercept[", country_num_in_mcmc, "]")]
    anyextra_country_intercept <- anyextra_samples$country_intercept[, paste0("anyextra_country_intercept[", country_num_in_mcmc, "]")]

  } else {
    # Country has NO data - sample from hierarchical prior

    # CIG: Sample full parameters (intercept + splines) from prior
    cig_prior_params <- sample_cig_params_from_prior(
      region_num = region_num,
      mcmc_samples = combined_samples_matrix,
      n_samples = n_samples,
      nAgeSpline = nAgeSpline,
      nCohortSpline = nCohortSpline
    )
    cig_country_intercept <- cig_prior_params$country_intercept
    cig_age_spline <- cig_prior_params$age_splines
    cig_cohort_spline <- cig_prior_params$cohort_splines

    # SMKEXTRA/ANYEXTRA: Sample only intercepts
    smkextra_country_intercept <- sample_extra_intercept_from_prior(
      head_name = "smkextra",
      mcmc_samples = combined_samples_matrix,
      n_samples = n_samples
    )
    anyextra_country_intercept <- sample_extra_intercept_from_prior(
      head_name = "anyextra",
      mcmc_samples = combined_samples_matrix,
      n_samples = n_samples
    )
  }

  # ==================================================================
  # GENERATE PREDICTIONS FOR EACH DATA POINT
  # ==================================================================

  n_pred <- nrow(current_data)
  predictions_cig <- matrix(0, nrow = n_pred, ncol = n_samples)
  predictions_smoked <- matrix(0, nrow = n_pred, ncol = n_samples)
  predictions_any <- matrix(0, nrow = n_pred, ncol = n_samples)

  for (j in 1:n_pred) {
    # Extract covariates
    def_code_binary <- current_data$Def_Code_Binary[j]
    age_spline_values <- as.numeric(current_data[j, grep("^age_spline_", names(current_data))])
    age_linear_smooth <- current_data$age_linear_smooth[j]
    spline_weight_value <- current_data$spline_weight_var[j]
    linear_weight_value <- current_data$linear_weight_var[j]
    cohort_spline_values <- as.numeric(current_data[j, grep("^cohort_spline_", names(current_data))])
    age_cohort_interaction_values <- as.numeric(current_data[j, grep("^age_cohort_", names(current_data))])

    # ------------------------------------------------------------------
    # CIG: Full model with country-specific splines
    # ------------------------------------------------------------------
    mu_cig <- cig_samples$global_intercept +
      cig_samples$region_intercept[, paste0("cig_region_intercept[", region_num, "]")] +
      cig_country_intercept +
      def_code_shared_samples * def_code_binary +
      spline_weight_value * as.matrix(cig_age_spline) %*% age_spline_values +
      linear_weight_value * cig_samples$age_linear_smooth * age_linear_smooth +
      as.matrix(cig_cohort_spline) %*% cohort_spline_values +
      as.matrix(cig_samples$age_cohort) %*% age_cohort_interaction_values

    # ------------------------------------------------------------------
    # SMKEXTRA: Country intercept + regional splines
    # ------------------------------------------------------------------
    mu_smkextra <- smkextra_samples$global_intercept +
      smkextra_samples$region_intercept[, paste0("smkextra_region_intercept[", region_num, "]")] +
      smkextra_country_intercept +
      0.3 * def_code_shared_samples * def_code_binary +
      spline_weight_value * as.matrix(smkextra_samples$age_spline_regional[, paste0("smkextra_age_spline_region_mean[", region_num, ", ", 1:length(age_spline_values), "]")]) %*% age_spline_values +
      linear_weight_value * smkextra_samples$age_linear_smooth * age_linear_smooth +
      as.matrix(smkextra_samples$cohort_spline_regional[, paste0("smkextra_cohort_spline_region_mean[", region_num, ", ", 1:length(cohort_spline_values), "]")]) %*% cohort_spline_values

    # ------------------------------------------------------------------
    # ANYEXTRA: Country intercept + regional splines
    # ------------------------------------------------------------------
    mu_anyextra <- anyextra_samples$global_intercept +
      anyextra_samples$region_intercept[, paste0("anyextra_region_intercept[", region_num, "]")] +
      anyextra_country_intercept +
      0.3 * def_code_shared_samples * def_code_binary +
      spline_weight_value * as.matrix(anyextra_samples$age_spline_regional[, paste0("anyextra_age_spline_region_mean[", region_num, ", ", 1:length(age_spline_values), "]")]) %*% age_spline_values +
      linear_weight_value * anyextra_samples$age_linear_smooth * age_linear_smooth +
      as.matrix(anyextra_samples$cohort_spline_regional[, paste0("anyextra_cohort_spline_region_mean[", region_num, ", ", 1:length(cohort_spline_values), "]")]) %*% cohort_spline_values

    # ------------------------------------------------------------------
    # Store CIG predictions
    # ------------------------------------------------------------------
    predictions_cig[j, ] <- mu_cig

    # ------------------------------------------------------------------
    # STICK-BREAKING: Compute smoked and any tobacco
    # ------------------------------------------------------------------
    p_cig <- plogis(mu_cig)
    p_smoked_full <- p_cig + plogis(mu_smkextra) * (1 - p_cig)
    predictions_smoked[j, ] <- log(p_smoked_full / (1 - p_smoked_full))

    p_any_full <- p_smoked_full + plogis(mu_anyextra) * (1 - p_smoked_full)
    predictions_any[j, ] <- log(p_any_full / (1 - p_any_full))
  }

  return(list(
    predictions_cig = predictions_cig,
    predictions_smoked = predictions_smoked,
    predictions_any = predictions_any
  ))
}


#' Create lookup for all countries (data + no-data)
#'
#' @param country_region_master Complete country-region mapping from country_region_manual
#' @param data_country_lookup Lookup table from countries with MCMC data
#' @param region_to_num Mapping from region name to region number
#' @return Data frame with all countries and their region mappings
create_all_country_lookup <- function(country_region_master,
                                       data_country_lookup,
                                       region_to_num) {

  # Start with master list
  all_countries <- unique(country_region_master$wb_country_abv)

  # Mark which have data
  data_countries <- unique(data_country_lookup$wb_country_abv)

  result <- data.frame(
    wb_country_abv = all_countries,
    stringsAsFactors = FALSE
  ) %>%
    left_join(country_region_master, by = "wb_country_abv") %>%
    mutate(
      has_data = wb_country_abv %in% data_countries,
      Region_Num = region_to_num[region_consolidated]
    ) %>%
    left_join(
      data_country_lookup %>% select(wb_country_abv, Num_Country),
      by = "wb_country_abv"
    )

  return(result)
}


#' Validate uncertainty is larger for countries without data
#'
#' This is a sanity check - countries without data should have wider credible intervals
#'
#' @param predictions_list Named list of prediction results with has_data flags
#' @param indicator Which indicator to summarize (cig, smoked, any)
#' @return Summary data frame
validate_uncertainty <- function(predictions_list, indicator = "cig") {

  summary_data <- data.frame(
    country = character(),
    has_data = logical(),
    region = character(),
    mean_pred = numeric(),
    ci_width = numeric(),
    stringsAsFactors = FALSE
  )

  pred_col <- switch(indicator,
                     "cig" = "predictions_cig",
                     "smoked" = "predictions_smoked",
                     "any" = "predictions_any")

  for (country_code in names(predictions_list)) {
    pred_info <- predictions_list[[country_code]]
    pred_matrix <- pred_info[[pred_col]]

    if (is.null(pred_matrix) || nrow(pred_matrix) == 0) next

    # Take mean across prediction points, then summarize across MCMC
    mean_across_points <- colMeans(plogis(pred_matrix))

    summary_data <- rbind(summary_data, data.frame(
      country = country_code,
      has_data = pred_info$has_data,
      region = pred_info$region,
      mean_pred = mean(mean_across_points),
      ci_width = quantile(mean_across_points, 0.975) - quantile(mean_across_points, 0.025),
      stringsAsFactors = FALSE
    ))
  }

  # Compare CI widths
  if (sum(summary_data$has_data) > 0 && sum(!summary_data$has_data) > 0) {
    data_ci <- mean(summary_data$ci_width[summary_data$has_data])
    nodata_ci <- mean(summary_data$ci_width[!summary_data$has_data])

    cat("\n  === Uncertainty Validation ===\n")
    cat(sprintf("  Mean CI width (countries WITH data):    %.4f\n", data_ci))
    cat(sprintf("  Mean CI width (countries WITHOUT data): %.4f\n", nodata_ci))
    cat(sprintf("  Ratio (should be > 1): %.2f\n", nodata_ci / data_ci))

    if (nodata_ci > data_ci) {
      cat("  [OK] Countries without data have larger uncertainty (as expected)\n")
    } else {
      cat("  [WARNING] Unexpected: countries without data have SMALLER uncertainty\n")
    }
  } else {
    cat("\n  Cannot validate uncertainty - need both data and no-data countries\n")
  }

  return(summary_data)
}


cat("\n  Module loaded: prediction_all_countries.R\n")
cat("  Functions available:\n")
cat("    - sample_cig_params_from_prior()\n")
cat("    - sample_extra_intercept_from_prior()\n")
cat("    - predict_single_country_all_heads()\n")
cat("    - create_all_country_lookup()\n")
cat("    - validate_uncertainty()\n\n")
