#########################################################################################
#
#                    WHO TOBACCO CONTROL PREVALENCE PROJECTION MODEL
#                    07_run_country_model.R - Country-Specific Model Fitting
#                                   VERSION 2.3.2
#
#   Contains: Section 11-12 from pipeline_monolith.R
#     - Country-specific model fitting loop (all 191 countries × 2 sexes)
#     - Vectorized prediction with OPT-1 through OPT-6 optimizations
#     - Country-specific weighted prevalence trends
#
#   PERFORMANCE OPTIMIZATIONS INCLUDED:
#     - OPT-1: Fixed clearCompiled order (prevents DLL memory leak)
#     - OPT-2: Precompute age components outside country loop
#     - OPT-3: Vectorized tensor product
#     - OPT-4: Matrix-based prediction
#     - OPT-5: rowMeans instead of apply
#     - OPT-6: Single ns() call for all cohort years
#
#   Requires: Global model fitted (06_run_global_model.R), country_priors/ directory
#   Outputs: final_predictions_country_specific.csv, country_specific_ac_nested/
#
#   EXTRACTED FROM: pipeline_monolith.R v2.3.2
#
#########################################################################################

cat("\n================================================================\n")
cat("  FITTING COUNTRY-SPECIFIC MODELS (NIMBLE)\n")
cat("================================================================\n")

if (!dir.exists("results/country_specific_ac_nested")) {
  dir.create("results/country_specific_ac_nested", recursive = TRUE)
}

# [MEM-OPT] Predictions are written to disk per-country (predictions.csv).
# After the loop, they are read back and combined. This avoids accumulating
# hundreds of data frames in RAM during the 370-iteration loop.

for (gender in genders) {
  cat(sprintf("\nCountry-specific models: %s\n", gender))

  gender_dir <- file.path("results/country_specific_ac_nested", gender)
  dir.create(gender_dir, showWarnings = FALSE)

  # Load regional structure with centering constants
  regional_info <- readRDS(file.path("processing", gender, "regional_structure.rds"))

  # ============================================================
  # [OPT-2] PRECOMPUTE AGE-ONLY COMPONENTS (shared across ALL countries)
  #
  # These matrices depend only on the fixed age grid (15-100).
  # Computing them once outside the country loop avoids redundant
  # ~86 ns() + outer() calls per country (~30,000 total for 360 countries).
  # ============================================================

  age_midpoints_shared <- seq(15, 100, by = 1)
  n_ages <- length(age_midpoints_shared)

  # Sigmoid weights for smooth elderly transition (per age)
  sigmoid_input_shared <- (age_midpoints_shared - TRANSITION_START) / TRANSITION_WIDTH
  spline_weight_shared <- 1 / (1 + exp(sigmoid_input_shared))  # length n_ages
  linear_weight_shared <- 1 - spline_weight_shared

  # Age linear centered (per age)
  age_linear_shared <- pmax(0, age_midpoints_shared - TRANSITION_START)
  age_linear_centered_shared <- age_linear_shared - AGE_LINEAR_CENTER_CONSTANT

  # Pre-multiply: linear_weight * age_linear_centered (used in every prediction row)
  linear_age_product_shared <- linear_weight_shared * age_linear_centered_shared  # length n_ages

  # Precompute age spline basis (n_ages x nAgeSpline)
  age_spline_mat_shared <- as.matrix(ns(
    age_midpoints_shared,
    knots = regional_info$age_spline_knots,
    Boundary.knots = regional_info$age_spline_boundary
  ))
  nAgeSpline_shared <- ncol(age_spline_mat_shared)

  # Helper vector for outer products
  ones_ages <- rep(1, n_ages)

  cat(sprintf("  Precomputed shared age components: %d ages, %d spline bases\n",
              n_ages, nAgeSpline_shared))

  countries <- list.files(file.path("country_priors", gender), pattern = "_regional_ac_priors_nested.csv")
  countries <- sub("_regional_ac_priors_nested.csv", "", countries)

  # [MEM-OPT] Counter for periodic memory monitoring
  country_counter <- 0L

  for (country_code in countries) {
    country_counter <- country_counter + 1L

    # [MEM-OPT] Periodic memory check (every 25 countries)
    if (country_counter %% 25 == 0) {
      cat(sprintf("  [Progress: country %d/%d]\n", country_counter, length(countries)))
      gc()  # Regular gc only — NOT gc(full=TRUE) which corrupts NIMBLE state
    }

    country_full_name <- country_name_mapping[country_code]
    cat(sprintf("  %s\n", country_full_name))

    country_dir <- file.path(gender_dir, country_full_name)
    dir.create(country_dir, showWarnings = FALSE)

    priors_file <- file.path("country_priors", gender, paste0(country_code, "_regional_ac_priors_nested.csv"))

    if (!file.exists(priors_file)) {
      warning(sprintf("No priors for %s, %s", country_code, gender))
      next
    }

    priors_df <- read.csv(priors_file)

    # Helper function to get prior values
    get_prior <- function(param_name) {
      row <- priors_df[priors_df$parameter == param_name, ]
      if (nrow(row) == 0) return(list(mean = 0, sd = 1))
      m <- row$mean[1]
      s <- row$sd[1]
      # Guard against NaN/NA from any source
      if (!is.finite(m)) m <- 0
      if (!is.finite(s) || s <= 0) s <- 1
      list(mean = m, sd = max(s, 1e-6))
    }

    # Helper function to convert SD to precision safely
    sd_to_prec <- function(sd_val) {
      sd_val[!is.finite(sd_val)] <- 1  # Replace NaN/NA/Inf with safe default
      sd_safe <- pmax(sd_val, 0.01)  # Minimum SD = 0.01
      prec <- 1 / (sd_safe^2)
      return(pmin(prec, 10000))  # Maximum precision = 10000
    }

    # Filter country data
    country_data <- clean_data %>%
      filter(wb_country_abv == country_code, sex == gender)

    if (nrow(country_data) == 0) {
      warning(sprintf("No data for %s, %s", country_code, gender))
      next
    }

    country_data <- country_data %>%
      mutate(
        Num_Survey = as.numeric(factor(survey)),
        Num_Year   = as.numeric(year),
        age_range  = end_age - start_age,
        weight     = 1 / (age_range + 1),
        weight     = weight / mean(weight)
      )

    # ---- Prepare NIMBLE Constants for Country-Specific Model ----

    nAgeSpline <- ncol(country_data[, grep("^age_spline_", names(country_data))])
    nCohortSpline <- ncol(country_data[, grep("^cohort_spline_", names(country_data))])
    nAgeXCohortSplines <- ncol(country_data[, grep("^age_cohort_", names(country_data))])
    nSurvey <- length(unique(country_data$Num_Survey))

    # Extract prior vectors with proper ordering
    get_prior_vector <- function(prefix, n) {
      means <- numeric(n)
      sds <- numeric(n)
      for (i in 1:n) {
        param_name <- paste0(prefix, i)
        row <- which(priors_df$parameter == param_name)
        if (length(row) == 1) {
          means[i] <- priors_df$mean[row]
          sds[i] <- priors_df$sd[row]
        } else {
          means[i] <- 0
          sds[i] <- 0.5
        }
      }
      list(means = means, sds = sds)
    }

    cig_age_priors <- get_prior_vector("cig_age_spline_", nAgeSpline)
    cig_cohort_priors <- get_prior_vector("cig_cohort_spline_", nCohortSpline)
    cig_ac_priors <- get_prior_vector("cig_age_cohort_interaction_", nAgeXCohortSplines)

    smkextra_age_priors <- get_prior_vector("smkextra_age_spline_", nAgeSpline)
    smkextra_cohort_priors <- get_prior_vector("smkextra_cohort_spline_", nCohortSpline)

    anyextra_age_priors <- get_prior_vector("anyextra_age_spline_", nAgeSpline)
    anyextra_cohort_priors <- get_prior_vector("anyextra_cohort_spline_", nCohortSpline)

    # ---- Prepare NIMBLE Data ----

    nimble_constants_country <- list(
      N = nrow(country_data),
      nAgeSpline = nAgeSpline,
      nCohortSpline = nCohortSpline,
      nAgeXCohortSplines = nAgeXCohortSplines,
      nSurvey = nSurvey,
      Survey = as.integer(country_data$Num_Survey)
    )

    nimble_data_country <- list(
      Prevalence = country_data$prevalence,
      Def_Code_Binary = country_data$def_code_binary,
      Type_Cig = country_data$Type_Cig,
      Type_Smoked = country_data$Type_Smoked,
      Type_Any = country_data$Type_Any,
      age_spline_matrix = as.matrix(country_data[, grep("^age_spline_", names(country_data))]),
      cohort_spline_matrix = as.matrix(country_data[, grep("^cohort_spline_", names(country_data))]),
      age_cohort_interaction_matrix = as.matrix(country_data[, grep("^age_cohort_", names(country_data))]),
      age_linear_smooth = country_data$age_linear_smooth,
      spline_weight_var = country_data$spline_weight_var,
      linear_weight_var = country_data$linear_weight_var,
      weight = country_data$weight,
      # Prior hyperparameters
      cig_intercept_prior_mean = get_prior("cig_intercept")$mean,
      cig_intercept_prior_prec = sd_to_prec(get_prior("cig_intercept")$sd),
      cig_def_code_prior_mean = get_prior("cig_def_code_shared")$mean,
      cig_def_code_prior_prec = sd_to_prec(get_prior("cig_def_code_shared")$sd),
      cig_age_linear_prior_mean = get_prior("cig_age_linear_smooth_effect")$mean,
      cig_age_linear_prior_prec = sd_to_prec(get_prior("cig_age_linear_smooth_effect")$sd),
      cig_age_spline_prior_means = cig_age_priors$means,
      cig_age_spline_prior_precs = sd_to_prec(cig_age_priors$sds),
      cig_cohort_spline_prior_means = cig_cohort_priors$means,
      cig_cohort_spline_prior_precs = sd_to_prec(cig_cohort_priors$sds),
      cig_age_cohort_prior_means = cig_ac_priors$means,
      cig_age_cohort_prior_precs = sd_to_prec(cig_ac_priors$sds),
      smkextra_intercept_prior_mean = get_prior("smkextra_intercept")$mean,
      smkextra_intercept_prior_prec = sd_to_prec(get_prior("smkextra_intercept")$sd),
      smkextra_age_linear_prior_mean = get_prior("smkextra_age_linear_smooth_effect")$mean,
      smkextra_age_linear_prior_prec = sd_to_prec(get_prior("smkextra_age_linear_smooth_effect")$sd),
      smkextra_age_spline_prior_means = smkextra_age_priors$means,
      smkextra_age_spline_prior_precs = sd_to_prec(smkextra_age_priors$sds),
      smkextra_cohort_spline_prior_means = smkextra_cohort_priors$means,
      smkextra_cohort_spline_prior_precs = sd_to_prec(smkextra_cohort_priors$sds),
      anyextra_intercept_prior_mean = get_prior("anyextra_intercept")$mean,
      anyextra_intercept_prior_prec = sd_to_prec(get_prior("anyextra_intercept")$sd),
      anyextra_age_linear_prior_mean = get_prior("anyextra_age_linear_smooth_effect")$mean,
      anyextra_age_linear_prior_prec = sd_to_prec(get_prior("anyextra_age_linear_smooth_effect")$sd),
      anyextra_age_spline_prior_means = anyextra_age_priors$means,
      anyextra_age_spline_prior_precs = sd_to_prec(anyextra_age_priors$sds),
      anyextra_cohort_spline_prior_means = anyextra_cohort_priors$means,
      anyextra_cohort_spline_prior_precs = sd_to_prec(anyextra_cohort_priors$sds)
    )

    # ---- Generate Initial Values ----

    generate_country_inits <- function(seed_offset = 0) {
      set.seed(RANDOM_SEED + seed_offset)
      list(
        cig_def_code_shared = rnorm(1, get_prior("cig_def_code_shared")$mean, 0.1),
        cig_intercept = rnorm(1, get_prior("cig_intercept")$mean, 0.2),
        smkextra_intercept = rnorm(1, get_prior("smkextra_intercept")$mean, 0.2),
        anyextra_intercept = rnorm(1, get_prior("anyextra_intercept")$mean, 0.2),
        cig_age_linear_smooth_effect = runif(1, -0.05, -0.01),
        smkextra_age_linear_smooth_effect = runif(1, -0.05, -0.01),
        anyextra_age_linear_smooth_effect = runif(1, -0.05, -0.01),
        residual_sd = runif(1, 0.5, 1.0),
        survey_intercept_precision = rgamma(1, 3, 1),
        survey_intercept = rnorm(nSurvey, 0, 0.1),
        cig_age_spline = rnorm(nAgeSpline, cig_age_priors$means, 0.05),
        cig_cohort_spline = rnorm(nCohortSpline, cig_cohort_priors$means, 0.05),
        cig_age_cohort_interaction = rnorm(nAgeXCohortSplines, cig_ac_priors$means, 0.02),
        smkextra_age_spline = rnorm(nAgeSpline, smkextra_age_priors$means, 0.05),
        smkextra_cohort_spline = rnorm(nCohortSpline, smkextra_cohort_priors$means, 0.05),
        anyextra_age_spline = rnorm(nAgeSpline, anyextra_age_priors$means, 0.05),
        anyextra_cohort_spline = rnorm(nCohortSpline, anyextra_cohort_priors$means, 0.05)
      )
    }

    inits_list <- lapply(1:NUMBER_OF_CHAINS, function(i) generate_country_inits(i))

    # ---- Fit Country Model ----

    tryCatch({
      nimble_model <- nimbleModel(
        code = regional_country_specific_ac_model_nimble,
        constants = nimble_constants_country,
        data = nimble_data_country,
        inits = inits_list[[1]],
        name = paste0("CountryModel_", country_code, "_", gender)
      )

      mcmc_config <- configureMCMC(
        nimble_model,
        monitors = c(
          "cig_intercept", "cig_def_code_shared", "cig_age_linear_smooth_effect",
          "cig_age_spline", "cig_cohort_spline", "cig_age_cohort_interaction",
          "smkextra_intercept", "smkextra_age_linear_smooth_effect",
          "smkextra_age_spline", "smkextra_cohort_spline",
          "anyextra_intercept", "anyextra_age_linear_smooth_effect",
          "anyextra_age_spline", "anyextra_cohort_spline"
        ),
        thin = THINNING_INTERVAL,
        enableWAIC = FALSE,
        useConjugacy = FALSE
      )

      mcmc_built <- buildMCMC(mcmc_config)
      compiled_model <- compileNimble(nimble_model)
      compiled_mcmc <- compileNimble(mcmc_built, project = nimble_model)

      samples <- runMCMC(
        compiled_mcmc,
        niter = NUMBER_OF_BURN + NUMBER_OF_ITERATIONS,
        nburnin = NUMBER_OF_BURN,
        nchains = NUMBER_OF_CHAINS,
        inits = inits_list,
        thin = THINNING_INTERVAL,
        samplesAsCodaMCMC = TRUE,
        progressBar = FALSE
      )

      combined_samples <- do.call(rbind, lapply(samples, as.matrix))
      n_samples <- nrow(combined_samples)

      # [MEM-OPT] Free raw MCMC samples immediately - only combined_samples needed
      rm(samples)

      # ============================================================
      # EXTRACT PARAMETER SAMPLES
      # ============================================================

      cig_int_s <- combined_samples[, "cig_intercept"]
      def_shared_s <- combined_samples[, "cig_def_code_shared"]
      cig_age_lin_s <- combined_samples[, "cig_age_linear_smooth_effect"]
      cig_age_spline_s <- combined_samples[, grep("^cig_age_spline\\[", colnames(combined_samples))]
      cig_cohort_s <- combined_samples[, grep("^cig_cohort_spline\\[", colnames(combined_samples))]
      cig_ac_s <- combined_samples[, grep("^cig_age_cohort_interaction\\[", colnames(combined_samples))]

      smk_int_s <- combined_samples[, "smkextra_intercept"]
      smk_age_lin_s <- combined_samples[, "smkextra_age_linear_smooth_effect"]
      smk_age_spline_s <- combined_samples[, grep("^smkextra_age_spline\\[", colnames(combined_samples))]
      smk_cohort_s <- combined_samples[, grep("^smkextra_cohort_spline\\[", colnames(combined_samples))]

      any_int_s <- combined_samples[, "anyextra_intercept"]
      any_age_lin_s <- combined_samples[, "anyextra_age_linear_smooth_effect"]
      any_age_spline_s <- combined_samples[, grep("^anyextra_age_spline\\[", colnames(combined_samples))]
      any_cohort_s <- combined_samples[, grep("^anyextra_cohort_spline\\[", colnames(combined_samples))]

      # ============================================================
      # [OPT-6] PRECOMPUTE COHORT SPLINES FOR ALL YEARS
      # ============================================================

      min_year <- min(country_data$year)
      max_year <- target_year + PROJECTED_YEARS
      years_pred <- seq(min_year, max_year, by = 1)
      n_years <- length(years_pred)

      # All birth cohorts flattened: (n_years * n_ages) length
      all_birth_cohorts <- rep(years_pred, each = n_ages) -
        rep(age_midpoints_shared, times = n_years)

      # Single ns() call for all cohorts
      all_cohort_spline_mat <- as.matrix(ns(
        all_birth_cohorts,
        knots = regional_info$cohort_spline_knots,
        Boundary.knots = regional_info$cohort_spline_boundary
      ))
      nCohortSpline_pred <- ncol(all_cohort_spline_mat)

      # [OPT-3] Vectorized tensor product for all years
      nA <- nAgeSpline_shared
      nC <- nCohortSpline_pred
      n_int <- nA * nC

      # Pre-allocate storage indexed by year
      cohort_spline_by_year <- vector("list", n_years)
      interaction_by_year   <- vector("list", n_years)
      names(cohort_spline_by_year) <- as.character(years_pred)
      names(interaction_by_year)   <- as.character(years_pred)

      for (y_idx in 1:n_years) {
        row_start <- (y_idx - 1) * n_ages + 1
        row_end   <- y_idx * n_ages
        cohort_mat <- all_cohort_spline_mat[row_start:row_end, , drop = FALSE]
        cohort_spline_by_year[[y_idx]] <- cohort_mat

        # Build tensor product: n_ages x (nA * nC)
        int_mat <- matrix(0, n_ages, n_int)
        for (m in 1:nC) {
          col_start <- (m - 1) * nA + 1
          col_end   <- m * nA
          int_mat[, col_start:col_end] <- age_spline_mat_shared * cohort_mat[, m]
        }
        interaction_by_year[[y_idx]] <- int_mat
      }

      rm(all_cohort_spline_mat, all_birth_cohorts)

      # ============================================================
      # [OPT-4] VECTORIZED PREDICTION LOOP
      # ============================================================
      def_codes    <- c("daily_user", "current_user")
      def_binaries <- c(0, 1)

      results_list <- list()

      for (y_idx in 1:n_years) {
        current_year <- years_pred[y_idx]
        cohort_mat   <- cohort_spline_by_year[[y_idx]]
        int_mat      <- interaction_by_year[[y_idx]]
        birth_cohorts_year <- current_year - age_midpoints_shared

        for (d in 1:2) {
          current_def     <- def_codes[d]
          def_code_binary <- def_binaries[d]

          cig_dir    <- file.path(country_dir, paste0(current_def, "_cigarettes"))
          smoked_dir <- file.path(country_dir, paste0(current_def, "_any_smoked_tobacco"))
          any_dir    <- file.path(country_dir, paste0(current_def, "_any_tobacco_product"))

          dir.create(cig_dir, showWarnings = FALSE)
          dir.create(smoked_dir, showWarnings = FALSE)
          dir.create(any_dir, showWarnings = FALSE)

          # ============================================================
          # VECTORIZED PREDICTION (full matrix algebra)
          # All operations produce n_ages x n_samples matrices.
          # ============================================================

          # CIG: Full model with age-cohort interactions
          mu_cig <- ones_ages %o% cig_int_s +
            def_code_binary * (ones_ages %o% def_shared_s) +
            (age_spline_mat_shared %*% t(cig_age_spline_s)) * spline_weight_shared +
            outer(linear_age_product_shared, cig_age_lin_s) +
            cohort_mat %*% t(cig_cohort_s) +
            int_mat %*% t(cig_ac_s)

          # SMKEXTRA: NO age-cohort interactions (0.3x def code scaling)
          mu_smkextra <- ones_ages %o% smk_int_s +
            (0.3 * def_code_binary) * (ones_ages %o% def_shared_s) +
            (age_spline_mat_shared %*% t(smk_age_spline_s)) * spline_weight_shared +
            outer(linear_age_product_shared, smk_age_lin_s) +
            cohort_mat %*% t(smk_cohort_s)

          # ANYEXTRA: NO age-cohort interactions (0.3x def code scaling)
          mu_anyextra <- ones_ages %o% any_int_s +
            (0.3 * def_code_binary) * (ones_ages %o% def_shared_s) +
            (age_spline_mat_shared %*% t(any_age_spline_s)) * spline_weight_shared +
            outer(linear_age_product_shared, any_age_lin_s) +
            cohort_mat %*% t(any_cohort_s)

          # STICK-BREAKING CONSTRUCTION (element-wise, all n_ages x n_samples)
          p_cig    <- plogis(mu_cig)
          p_smoked <- p_cig + plogis(mu_smkextra) * (1 - p_cig)
          p_any    <- p_smoked + plogis(mu_anyextra) * (1 - p_smoked)

          # Convert to logit scale for storage
          predictions_cig    <- mu_cig
          predictions_smoked <- log(p_smoked / (1 - p_smoked))
          predictions_any    <- log(p_any / (1 - p_any))

          # Check ordering constraints
          has_violations <- check_ordering_per_draw(predictions_cig, predictions_smoked, predictions_any)
          if (has_violations) {
            warning(sprintf("Ordering violations: %s, %s, year %d", country_code, current_def, current_year))
          }

          # Save prediction matrices (logit scale, n_ages x n_samples)
          saveRDS(predictions_cig, file = file.path(cig_dir, paste0(current_year, ".rds")), compress = TRUE)
          saveRDS(predictions_smoked, file = file.path(smoked_dir, paste0(current_year, ".rds")), compress = TRUE)
          saveRDS(predictions_any, file = file.path(any_dir, paste0(current_year, ".rds")), compress = TRUE)

          # ============================================================
          # [OPT-5] VECTORIZED SUMMARY STATISTICS
          # Use rowMeans instead of apply(..., mean) for ~5x speedup.
          # ============================================================

          summarise_predictions <- function(p_matrix) {
            ci <- t(apply(p_matrix, 1, quantile, probs = c(0.025, 0.975)))
            data.frame(
              Prevalence = rowMeans(p_matrix),
              lower_ci   = ci[, 1],
              upper_ci   = ci[, 2]
            )
          }

          base_df <- data.frame(
            Year         = current_year,
            Age_Midpoint = age_midpoints_shared,
            Birth_Cohort = birth_cohorts_year,
            Def_Code_Binary = def_code_binary,
            Def_Code     = current_def
          )

          # Cigarettes
          cig_summary <- summarise_predictions(p_cig)
          results_list[[length(results_list) + 1]] <- cbind(
            base_df, cig_summary,
            Def_Type_Code = paste0(current_def, "_cigarettes")
          )

          # Any smoked tobacco
          smoked_summary <- summarise_predictions(p_smoked)
          results_list[[length(results_list) + 1]] <- cbind(
            base_df, smoked_summary,
            Def_Type_Code = paste0(current_def, "_any_smoked_tobacco")
          )

          # Any tobacco
          any_summary <- summarise_predictions(p_any)
          results_list[[length(results_list) + 1]] <- cbind(
            base_df, any_summary,
            Def_Type_Code = paste0(current_def, "_any_tobacco_product")
          )

          # Free intermediate matrices
          rm(mu_cig, mu_smkextra, mu_anyextra,
             p_cig, p_smoked, p_any,
             predictions_cig, predictions_smoked, predictions_any)
        }
      }

      # Free precomputed year-specific matrices
      rm(cohort_spline_by_year, interaction_by_year)

      # Combine results
      result         <- do.call(rbind, results_list)
      result$Country <- country_code
      result$Sex     <- gender

      prediction_results <- result %>% select(
        Year, Age_Midpoint, Birth_Cohort, Def_Type_Code,
        Sex, Country, Prevalence, lower_ci, upper_ci
      )

      write.csv(prediction_results, file = file.path(country_dir, "predictions.csv"), row.names = FALSE)

      # ============================================================
      # [OPT-1] FIXED CLEANUP ORDER
      #
      # Original bug: clearCompiled(nimble_model) was called BEFORE
      # saving a reference, then rm() was called, so clearCompiled
      # always failed silently, leaking the compiled DLL. Over 360+
      # iterations this caused progressive slowdown.
      #
      # Fix: Save reference, remove other objects, THEN clearCompiled.
      # ============================================================

      model_ref <- nimble_model  # Save reference before removing anything
      rm(compiled_model, compiled_mcmc, mcmc_built, results_list,
         cig_int_s, cig_age_spline_s, cig_age_lin_s, cig_cohort_s, cig_ac_s,
         smk_int_s, smk_age_spline_s, smk_age_lin_s, smk_cohort_s,
         any_int_s, any_age_spline_s, any_age_lin_s, any_cohort_s,
         def_shared_s)
      nimble::clearCompiled(model_ref)  # Now actually works - releases the DLL
      rm(nimble_model, model_ref)
      gc()

    }, error = function(e) {
      warning(sprintf("Error fitting country model for %s, %s: %s", country_code, gender, e$message))
      cat(sprintf("    ERROR: %s\n", e$message))
    })

    # [MEM-OPT] Safety net: remove lingering NIMBLE objects after errors.
    # IMPORTANT: Do NOT call clearCompiled() here — on error paths the model
    # may never have been fully compiled, and clearCompiled on a half-built
    # model corrupts NIMBLE's internal DLL state, causing NaN in subsequent
    # country models. Just rm() and let gc() handle finalization.
    for (obj_name in c("nimble_model", "model_ref", "compiled_model",
                       "compiled_mcmc", "mcmc_built")) {
      if (exists(obj_name, inherits = FALSE)) {
        try(rm(list = obj_name), silent = TRUE)
      }
    }
    gc()
  }
}

# [MEM-OPT] Read back all country predictions from disk instead of growing list in RAM.
pred_files <- list.files("results/country_specific_ac_nested",
                         pattern = "^predictions\\.csv$",
                         recursive = TRUE, full.names = TRUE)

if (length(pred_files) > 0) {
  cat(sprintf("  Reading back %d country prediction files from disk...\n", length(pred_files)))
  final_predictions_country_specific_ac <- do.call(rbind,
    lapply(pred_files, function(f) read.csv(f, stringsAsFactors = FALSE)))
  rownames(final_predictions_country_specific_ac) <- NULL
} else {
  warning("No country-specific predictions were generated. All models may have failed.")
  final_predictions_country_specific_ac <- data.frame(
    Year = integer(), Age_Midpoint = numeric(), Birth_Cohort = numeric(),
    Def_Type_Code = character(), Sex = character(), Country = character(),
    Prevalence = numeric(), lower_ci = numeric(), upper_ci = numeric(),
    stringsAsFactors = FALSE
  )
}

write.csv(
  final_predictions_country_specific_ac,
  file = "results/final_predictions_country_specific.csv",
  row.names = FALSE
)

cat("\nCountry-specific model fitting complete (NIMBLE)\n")
