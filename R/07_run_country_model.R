#########################################################################################
#
#                    WHO TOBACCO CONTROL PREVALENCE PROJECTION MODEL
#                    07_run_country_model.R - Country-Specific Model Fitting
#                                   VERSION 2.4.0
#
#   Contains: Country-specific model fitting, prediction, weighted results
#     - Country-specific model fitting loop (all countries x 2 sexes)
#     - Vectorized prediction with OPT-2 through OPT-6 optimizations
#     - SUBPROCESS ISOLATION: Each country runs in a fresh R process via callr
#       to prevent DLL accumulation and progressive RAM growth
#
#   PERFORMANCE OPTIMIZATIONS INCLUDED:
#     - OPT-2: Precompute age components outside country loop
#     - OPT-3: Vectorized tensor product
#     - OPT-4: Matrix-based prediction
#     - OPT-5: rowMeans instead of apply
#     - OPT-6: Single ns() call for all cohort years
#     - SUBPROCESS: Each country in isolated R process (prevents DLL leak on Windows)
#
#   DLL ISOLATION:
#     Each NIMBLE model compiles 2 C++ DLLs (model + MCMC). On Windows,
#     dyn.unload() does not reliably release DLLs, causing "maximal number
#     of DLLs reached" after ~200 countries. Running each country in a fresh
#     R subprocess via callr::r() guarantees DLLs are released when the
#     subprocess exits (OS-level guarantee). Overhead is ~5-10 seconds
#     per country, negligible vs 20-30 min MCMC fitting time.
#
#   Requires: Global model fitted (06_run_global_model.R), country_priors/ directory
#   Outputs: final_predictions_country_specific.csv, country_specific_ac_nested/
#
#   Originally extracted from monolith, now the canonical source
#
#########################################################################################

# ============================================================
# WORKER FUNCTION: Runs in a fresh R subprocess via callr::r()
#
# Receives: country_code, gender, path to shared data RDS
# Produces: predictions.csv, posterior_samples.rds, prediction matrices
# Returns:  list(success = TRUE/FALSE)
# ============================================================

fit_single_country_model <- function(country_code, gender, shared_data_path) {

  # Load only packages needed for model fitting
  suppressPackageStartupMessages({
    library(nimble)
    library(coda)
    library(splines)
    library(dplyr)
  })
  nimbleOptions(verbose = FALSE)
  nimbleOptions(MCMCprogressBar = FALSE)
  nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = FALSE)

  # Load shared data prepared by the main process
  d <- readRDS(shared_data_path)

  # Extract config
  RANDOM_SEED          <- d$RANDOM_SEED
  NUMBER_OF_CHAINS     <- d$NUMBER_OF_CHAINS
  NUMBER_OF_BURN       <- d$NUMBER_OF_BURN
  NUMBER_OF_ITERATIONS <- d$NUMBER_OF_ITERATIONS
  THINNING_INTERVAL    <- d$THINNING_INTERVAL
  target_year          <- d$target_year
  PROJECTED_YEARS      <- d$PROJECTED_YEARS

  # Extract data
  regional_info        <- d$regional_info
  clean_data           <- d$clean_data
  country_name_mapping <- d$country_name_mapping
  model_code           <- d$model_code

  # Precomputed age components
  age_midpoints_shared       <- d$age_midpoints_shared
  n_ages                     <- d$n_ages
  spline_weight_shared       <- d$spline_weight_shared
  linear_weight_shared       <- d$linear_weight_shared
  linear_age_product_shared  <- d$linear_age_product_shared
  age_spline_mat_shared      <- d$age_spline_mat_shared
  nAgeSpline_shared          <- d$nAgeSpline_shared
  ones_ages                  <- d$ones_ages
  gender_dir                 <- d$gender_dir

  rm(d)

  # Helper: check ordering constraints (stick-breaking: cig <= smoked <= any)
  check_ordering_per_draw <- function(logit_cig, logit_smk, logit_any, tolerance = 1e-6) {
    pc <- plogis(logit_cig)
    ps <- plogis(logit_smk)
    pa <- plogis(logit_any)
    return(sum(pc > ps + tolerance) > 0 || sum(ps > pa + tolerance) > 0)
  }

  # --- Country Setup ---

  country_full_name <- country_name_mapping[country_code]
  country_dir <- file.path(gender_dir, country_full_name)
  dir.create(country_dir, showWarnings = FALSE, recursive = TRUE)

  priors_file <- file.path("country_priors", gender, paste0(country_code, "_regional_ac_priors_nested.csv"))

  if (!file.exists(priors_file)) {
    cat(sprintf("    WARNING: No priors for %s, %s - skipping\n", country_code, gender))
    return(list(success = FALSE, error = "No priors file"))
  }

  priors_df <- read.csv(priors_file)

  # Helper function to get prior values
  get_prior <- function(param_name) {
    row <- priors_df[priors_df$parameter == param_name, ]
    if (nrow(row) == 0) return(list(mean = 0, sd = 1))
    list(mean = row$mean, sd = max(row$sd, 1e-6))
  }

  # Helper function to convert SD to precision safely
  sd_to_prec <- function(sd_val) {
    sd_safe <- pmax(sd_val, 0.01)
    prec <- 1 / (sd_safe^2)
    return(pmin(prec, 10000))
  }

  # Filter country data
  country_data <- clean_data %>%
    filter(wb_country_abv == country_code, sex == gender)

  if (nrow(country_data) == 0) {
    cat(sprintf("    WARNING: No data for %s, %s - skipping\n", country_code, gender))
    return(list(success = FALSE, error = "No data"))
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
        sds[i] <- 1
      }
    }
    list(means = means, sds = sds)
  }

  cig_age_spline_priors <- get_prior_vector("cig_age_spline_", nAgeSpline)
  cig_cohort_spline_priors <- get_prior_vector("cig_cohort_spline_", nCohortSpline)
  smkextra_age_spline_priors <- get_prior_vector("smkextra_age_spline_", nAgeSpline)
  smkextra_cohort_spline_priors <- get_prior_vector("smkextra_cohort_spline_", nCohortSpline)

  anyextra_age_spline_priors <- get_prior_vector("anyextra_age_spline_", nAgeSpline)
  anyextra_cohort_spline_priors <- get_prior_vector("anyextra_cohort_spline_", nCohortSpline)

  # ---- Prepare NIMBLE Data ----

  # CORRECTED: No smkextra/anyextra age-cohort constants
  nimble_constants_country <- list(
    N = nrow(country_data),
    nSurvey = nSurvey,
    nAgeSpline = nAgeSpline,
    nCohortSpline = nCohortSpline,
    Survey = as.integer(country_data$Num_Survey),

    # CIG priors
    cig_def_code_shared_prior_mean = get_prior("cig_def_code_shared")$mean,
    cig_def_code_shared_prior_prec = sd_to_prec(get_prior("cig_def_code_shared")$sd),

    cig_intercept_prior_mean = get_prior("cig_intercept")$mean,
    cig_intercept_prior_prec = sd_to_prec(get_prior("cig_intercept")$sd),

    cig_age_linear_smooth_effect_prior_mean = get_prior("cig_age_linear_smooth_effect")$mean,
    cig_age_linear_smooth_effect_prior_prec = sd_to_prec(get_prior("cig_age_linear_smooth_effect")$sd),

    # CIG spline priors (using ordered extraction)
    cig_age_spline_prior_means = cig_age_spline_priors$means,
    cig_age_spline_prior_precs = sd_to_prec(cig_age_spline_priors$sds),

    cig_cohort_spline_prior_means = cig_cohort_spline_priors$means,
    cig_cohort_spline_prior_precs = sd_to_prec(cig_cohort_spline_priors$sds),

    # SMKEXTRA priors
    smkextra_intercept_prior_mean = get_prior("smkextra_intercept")$mean,
    smkextra_intercept_prior_prec = sd_to_prec(get_prior("smkextra_intercept")$sd),

    smkextra_age_linear_smooth_effect_prior_mean = get_prior("smkextra_age_linear_smooth_effect")$mean,
    smkextra_age_linear_smooth_effect_prior_prec = sd_to_prec(get_prior("smkextra_age_linear_smooth_effect")$sd),

    # SMKEXTRA spline priors (using ordered extraction)
    smkextra_age_spline_prior_means = smkextra_age_spline_priors$means,
    smkextra_age_spline_prior_precs = sd_to_prec(smkextra_age_spline_priors$sds),

    smkextra_cohort_spline_prior_means = smkextra_cohort_spline_priors$means,
    smkextra_cohort_spline_prior_precs = sd_to_prec(smkextra_cohort_spline_priors$sds),

    # ANYEXTRA priors (NO age-cohort)
    anyextra_intercept_prior_mean = get_prior("anyextra_intercept")$mean,
    anyextra_intercept_prior_prec = sd_to_prec(get_prior("anyextra_intercept")$sd),

    anyextra_age_linear_smooth_effect_prior_mean = get_prior("anyextra_age_linear_smooth_effect")$mean,
    anyextra_age_linear_smooth_effect_prior_prec = sd_to_prec(get_prior("anyextra_age_linear_smooth_effect")$sd),

    # ANYEXTRA spline priors (using ordered extraction)
    anyextra_age_spline_prior_means = anyextra_age_spline_priors$means,
    anyextra_age_spline_prior_precs = sd_to_prec(anyextra_age_spline_priors$sds),

    anyextra_cohort_spline_prior_means = anyextra_cohort_spline_priors$means,
    anyextra_cohort_spline_prior_precs = sd_to_prec(anyextra_cohort_spline_priors$sds)
  )

  # DATA
  nimble_data_country <- list(
    Prevalence = country_data$prevalence,
    Def_Code_Binary = country_data$def_code_binary,
    Type_Cig = country_data$Type_Cig,
    Type_Smoked = country_data$Type_Smoked,
    Type_Any = country_data$Type_Any,
    age_spline_matrix = as.matrix(country_data[, grep("^age_spline_", names(country_data))]),
    cohort_spline_matrix = as.matrix(country_data[, grep("^cohort_spline_", names(country_data))]),
    age_linear_smooth = country_data$age_linear_smooth,
    spline_weight_var = country_data$spline_weight_var,
    linear_weight_var = country_data$linear_weight_var,
    weight = country_data$weight
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
      # CIG spline arrays (using ordered prior means)
      cig_age_spline = rnorm(nAgeSpline, cig_age_spline_priors$means, 0.05),
      cig_cohort_spline = rnorm(nCohortSpline, cig_cohort_spline_priors$means, 0.05),
      # SMKEXTRA spline arrays
      smkextra_age_spline = rnorm(nAgeSpline, smkextra_age_spline_priors$means, 0.05),
      smkextra_cohort_spline = rnorm(nCohortSpline, smkextra_cohort_spline_priors$means, 0.05),

      # ANYEXTRA spline arrays (NO age-cohort)
      anyextra_age_spline = rnorm(nAgeSpline, anyextra_age_spline_priors$means, 0.05),
      anyextra_cohort_spline = rnorm(nCohortSpline, anyextra_cohort_spline_priors$means, 0.05)
    )
  }

  inits_list <- lapply(1:NUMBER_OF_CHAINS, function(i) generate_country_inits(i))

  # ---- Fit Country Model ----
  # No tryCatch needed: errors propagate to the parent process.
  # No cleanup needed: subprocess exit releases all DLLs and memory.

  nimble_model <- nimbleModel(
    code = model_code,
    constants = nimble_constants_country,
    data = nimble_data_country,
    inits = inits_list[[1]],
    name = paste0("CountryModel_", country_code, "_", gender)
  )

  # Check for uninitialized nodes
  uninit_nodes <- nimble_model$initializeInfo()$uninitializedNodes
  if (length(uninit_nodes) > 0) {
    cat(sprintf("    WARNING: Uninitialized nodes: %s\n",
                paste(head(uninit_nodes, 5), collapse = ", ")))
  }

  # Calculate initial log probability
  init_logprob <- try(nimble_model$calculate(), silent = TRUE)
  if (inherits(init_logprob, "try-error") || !is.finite(init_logprob)) {
    cat(sprintf("    WARNING: Initial logProb is invalid: %s\n", init_logprob))
    node_logprobs <- sapply(nimble_model$getNodeNames(stochOnly = TRUE), function(n) {
      tryCatch(nimble_model$calculate(n), error = function(e) NA)
    })
    bad_nodes <- names(node_logprobs)[!is.finite(node_logprobs)]
    if (length(bad_nodes) > 0) {
      cat(sprintf("    Problematic nodes: %s\n", paste(head(bad_nodes, 10), collapse = ", ")))
    }
  } else {
    cat(sprintf("    Initial logProb: %.2f\n", init_logprob))
  }

  # CORRECTED: Monitors - no smkextra/anyextra age-cohort
  mcmc_config <- configureMCMC(
    nimble_model,
    useConjugacy = FALSE,
    monitors = c(
      "cig_def_code_shared",
      "cig_intercept", "cig_age_spline", "cig_age_linear_smooth_effect",
      "cig_cohort_spline",
      "smkextra_intercept", "smkextra_age_spline", "smkextra_age_linear_smooth_effect",
      "smkextra_cohort_spline",
      "anyextra_intercept", "anyextra_age_spline", "anyextra_age_linear_smooth_effect",
      "anyextra_cohort_spline",
      "residual_sd"
    ),
    thin = THINNING_INTERVAL,
    enableWAIC = FALSE
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
    progressBar = FALSE,
    summary = TRUE
  )

  combined_samples <- do.call(rbind, lapply(samples$samples, as.matrix))

  # Check for NaN/Inf
  n_invalid <- sum(!is.finite(combined_samples))
  if (n_invalid > 0) {
    cat(sprintf("    WARNING: %d non-finite values in samples\n", n_invalid))
  }

  # R-hat diagnostic
  gelman_diag <- try({
    param_vars <- apply(combined_samples, 2, var, na.rm = TRUE)
    varying_params <- names(param_vars)[param_vars > 1e-10]
    if (length(varying_params) > 0) {
      samples_filtered <- lapply(samples$samples, function(s) s[, varying_params, drop = FALSE])
      class(samples_filtered) <- "mcmc.list"
      gelman.diag(samples_filtered, multivariate = FALSE)
    } else NULL
  }, silent = TRUE)

  if (!inherits(gelman_diag, "try-error") && !is.null(gelman_diag)) {
    max_rhat <- max(gelman_diag$psrf[, "Point est."], na.rm = TRUE)
    if (is.finite(max_rhat)) {
      cat(sprintf("    R-hat: %.3f %s\n", max_rhat,
                  ifelse(max_rhat < 1.1, "(Good)", "(Check convergence)")))
    }
  }

  rm(samples)

  saveRDS(combined_samples, file = file.path(country_dir, "posterior_samples.rds"), compress = TRUE)

  # Extract samples for prediction
  extract_head_samples_country <- function(prefix, combined_samples) {
    list(
      intercept     = combined_samples[, paste0(prefix, "_intercept")],
      age_linear    = combined_samples[, paste0(prefix, "_age_linear_smooth_effect")],
      age_spline    = combined_samples[, grep(paste0("^", prefix, "_age_spline\\["), colnames(combined_samples)), drop = FALSE],
      cohort_spline = combined_samples[, grep(paste0("^", prefix, "_cohort_spline\\["), colnames(combined_samples)), drop = FALSE]
    )
  }

  cig_samples_full      <- extract_head_samples_country("cig", combined_samples)
  smkextra_samples_full <- extract_head_samples_country("smkextra", combined_samples)
  anyextra_samples_full <- extract_head_samples_country("anyextra", combined_samples)
  def_code_shared_full  <- combined_samples[, "cig_def_code_shared"]

  # ============================================================
  # [OPT-4] Pre-subset to sampled_indices ONCE (avoids repeated indexing)
  # ============================================================
  n_samples <- min(1000, nrow(combined_samples))
  sampled_indices <- sort(sample(1:nrow(combined_samples), n_samples))

  cig_int_s          <- cig_samples_full$intercept[sampled_indices]
  cig_age_spline_s   <- cig_samples_full$age_spline[sampled_indices, , drop = FALSE]
  cig_age_lin_s      <- cig_samples_full$age_linear[sampled_indices]
  cig_cohort_s       <- cig_samples_full$cohort_spline[sampled_indices, , drop = FALSE]
  smk_int_s          <- smkextra_samples_full$intercept[sampled_indices]
  smk_age_spline_s   <- smkextra_samples_full$age_spline[sampled_indices, , drop = FALSE]
  smk_age_lin_s      <- smkextra_samples_full$age_linear[sampled_indices]
  smk_cohort_s       <- smkextra_samples_full$cohort_spline[sampled_indices, , drop = FALSE]

  any_int_s          <- anyextra_samples_full$intercept[sampled_indices]
  any_age_spline_s   <- anyextra_samples_full$age_spline[sampled_indices, , drop = FALSE]
  any_age_lin_s      <- anyextra_samples_full$age_linear[sampled_indices]
  any_cohort_s       <- anyextra_samples_full$cohort_spline[sampled_indices, , drop = FALSE]

  def_shared_s       <- def_code_shared_full[sampled_indices]

  # Free full sample matrices (no longer needed)
  rm(cig_samples_full, smkextra_samples_full, anyextra_samples_full,
     def_code_shared_full, combined_samples)

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
  # Pre-allocate storage indexed by year
  cohort_spline_by_year <- vector("list", n_years)
  names(cohort_spline_by_year) <- as.character(years_pred)

  for (y_idx in 1:n_years) {
    row_start <- (y_idx - 1) * n_ages + 1
    row_end   <- y_idx * n_ages
    cohort_mat <- all_cohort_spline_mat[row_start:row_end, , drop = FALSE]
    cohort_spline_by_year[[y_idx]] <- cohort_mat
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
    birth_cohorts_year <- current_year - age_midpoints_shared

    for (d_idx in 1:2) {
      current_def     <- def_codes[d_idx]
      def_code_binary <- def_binaries[d_idx]

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

      # CIG
      mu_cig <- ones_ages %o% cig_int_s +
        def_code_binary * (ones_ages %o% def_shared_s) +
        (age_spline_mat_shared %*% t(cig_age_spline_s)) * spline_weight_shared +
        outer(linear_age_product_shared, cig_age_lin_s) +
        cohort_mat %*% t(cig_cohort_s)

      # SMKEXTRA (0.3x def code scaling)
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
  rm(cohort_spline_by_year)

  # Combine results
  result         <- do.call(rbind, results_list)
  result$Country <- country_code
  result$Sex     <- gender

  prediction_results <- result %>% select(
    Year, Age_Midpoint, Birth_Cohort, Def_Type_Code,
    Sex, Country, Prevalence, lower_ci, upper_ci
  )

  write.csv(prediction_results, file = file.path(country_dir, "predictions.csv"), row.names = FALSE)

  # No cleanup needed: subprocess exit releases all DLLs and reclaims memory.
  return(list(success = TRUE))
}


# ============================================================
# MAIN SCRIPT: Launch subprocesses for each country
# ============================================================

cat("\n================================================================\n")
cat("  FITTING COUNTRY-SPECIFIC MODELS (NIMBLE)\n")
cat("  Using subprocess isolation (callr) for DLL management\n")
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

  # ============================================================
  # SAVE SHARED DATA FOR SUBPROCESS ISOLATION
  #
  # All data shared across countries (config, age components,
  # regional info, full dataset, model code) is saved to a
  # temporary RDS file. Each subprocess reads this file once
  # instead of having the data serialized per-call.
  # ============================================================

  shared_data <- list(
    # Config constants
    RANDOM_SEED          = RANDOM_SEED,
    NUMBER_OF_CHAINS     = NUMBER_OF_CHAINS,
    NUMBER_OF_BURN       = NUMBER_OF_BURN,
    NUMBER_OF_ITERATIONS = NUMBER_OF_ITERATIONS,
    THINNING_INTERVAL    = THINNING_INTERVAL,
    target_year          = target_year,
    PROJECTED_YEARS      = PROJECTED_YEARS,

    # Regional and dataset
    regional_info        = regional_info,
    clean_data           = clean_data,
    country_name_mapping = country_name_mapping,
    model_code           = regional_country_specific_ac_model_nimble,

    # Precomputed age components
    age_midpoints_shared       = age_midpoints_shared,
    n_ages                     = n_ages,
    spline_weight_shared       = spline_weight_shared,
    linear_weight_shared       = linear_weight_shared,
    linear_age_product_shared  = linear_age_product_shared,
    age_spline_mat_shared      = age_spline_mat_shared,
    nAgeSpline_shared          = nAgeSpline_shared,
    ones_ages                  = ones_ages,
    gender_dir                 = gender_dir
  )

  shared_file <- tempfile(pattern = "who_shared_", fileext = ".rds")
  saveRDS(shared_data, shared_file)
  rm(shared_data)
  cat(sprintf("  Shared data saved to: %s\n", basename(shared_file)))

  countries <- list.files(file.path("country_priors", gender), pattern = "_regional_ac_priors_nested.csv")
  countries <- sub("_regional_ac_priors_nested.csv", "", countries)

  n_countries <- length(countries)
  cat(sprintf("  Countries to fit: %d\n", n_countries))

  n_skipped <- 0
  for (i in seq_along(countries)) {
    country_code <- countries[i]
    country_full_name <- country_name_mapping[country_code]

    # Skip if predictions already exist on disk (enables resume after interruption)
    country_pred_file <- file.path(gender_dir, country_full_name, "predictions.csv")
    if (file.exists(country_pred_file)) {
      n_skipped <- n_skipped + 1
      next
    }

    cat(sprintf("  [%d/%d] %s\n", i, n_countries, country_full_name))

    # Pre-check: skip if no data (avoids subprocess startup overhead)
    has_data <- any(clean_data$wb_country_abv == country_code & clean_data$sex == gender)
    if (!has_data) {
      warning(sprintf("No data for %s, %s", country_code, gender))
      next
    }

    tryCatch({
      callr::r(
        func = fit_single_country_model,
        args = list(
          country_code = country_code,
          gender = gender,
          shared_data_path = shared_file
        ),
        show = TRUE,
        timeout = 7200  # 2 hours per country
      )
    }, error = function(e) {
      warning(sprintf("Error fitting country model for %s, %s: %s", country_code, gender, e$message))
      cat(sprintf("    ERROR: %s\n", e$message))
    })
  }

  if (n_skipped > 0) {
    cat(sprintf("  Skipped %d/%d countries with existing predictions on disk\n", n_skipped, n_countries))
    cat("  (Delete results/country_specific_ac_nested/ to force re-fitting)\n")
  }

  # Clean up shared data file
  unlink(shared_file)
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

#########################################################################################
#           TARGET PREVALENCE CALCULATION (COUNTRY-SPECIFIC MODEL)
#
#   Target prevalence calculation for country-specific models
#   Calculates base-year weighted prevalence and target thresholds
#   using the country-specific model predictions saved to disk.
#
#   Requires: weights_cleaned (from 04_utils.R),
#             calculate_weighted_prevalence (from 04_utils.R),
#             target_prevalence_df (from 06_run_global_model.R, Global-only)
#   Outputs: Updates target_prevalence_df with Country model entries
#########################################################################################

cat("\nCalculating target prevalences for country-specific models...\n")

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = RANDOM_SEED)
clusterExport(cl, "calculate_weighted_prevalence")

country_target_list_all <- list()

for (gender in genders) {
  cat(sprintf("  Target prevalences (Country model): %s\n", gender))

  countries_country <- list.dirs(
    file.path("results/country_specific_ac_nested", gender),
    full.names = FALSE, recursive = FALSE
  )

  country_results <- foreach(
    country = countries_country,
    .packages = c("dplyr", "tidyr"),
    .combine = rbind,
    .errorhandling = 'remove'
  ) %dopar% {

    country_target_list <- list()
    country_path <- file.path("results/country_specific_ac_nested", gender, country)
    def_types <- list.dirs(country_path, full.names = FALSE, recursive = FALSE)
    def_types <- def_types[grepl("(cigarettes|tobacco)", def_types)]

    country_code <- names(country_name_mapping)[country_name_mapping == country]
    if (length(country_code) == 0) country_code <- country

    base_year_weights <- weights_cleaned %>%
      filter(area == country_code, sex == gender, year == BASE_YEAR) %>%
      arrange(age)

    if (nrow(base_year_weights) == 0) return(NULL)

    for (def_type in def_types) {
      base_year_file <- file.path(country_path, def_type, paste0(BASE_YEAR, ".rds"))

      if (!file.exists(base_year_file)) next

      iterations_matrix_logit <- readRDS(base_year_file)

      if (nrow(iterations_matrix_logit) != nrow(base_year_weights)) next

      weighted_probs <- calculate_weighted_prevalence(
        iterations_matrix_logit,
        base_year_weights$weight
      )

      base_year_mean  <- mean(weighted_probs)
      base_year_lower <- quantile(weighted_probs, 0.025)
      base_year_upper <- quantile(weighted_probs, 0.975)

      calculated_target_prevalence <- base_year_mean * (1 - REDUCTION_PERCENTAGE / 100)

      base_record <- data.frame(
        Sex               = gender,
        Country           = country_code,
        Def_Type_Code     = def_type,
        BaseYear          = BASE_YEAR,
        BaseYearPrevalence = base_year_mean,
        BaseYearLower     = base_year_lower,
        BaseYearUpper     = base_year_upper,
        Model_Type        = "Country",
        stringsAsFactors  = FALSE
      )

      # 30% Reduction Target
      calculated_record <- base_record %>%
        mutate(
          TargetPrevalence    = calculated_target_prevalence,
          ReductionPercentage = REDUCTION_PERCENTAGE,
          TargetType          = "30% Reduction",
          TargetValue         = paste0(REDUCTION_PERCENTAGE, "% reduction")
        )

      country_target_list[[length(country_target_list) + 1]] <- calculated_record

      # Absolute Target
      manual_record <- base_record %>%
        mutate(
          TargetPrevalence    = MANUAL_TARGET_PROPORTION,
          ReductionPercentage = round((1 - MANUAL_TARGET_PROPORTION / base_year_mean) * 100, 1),
          TargetType          = paste0(MANUAL_TARGET_PREVALENCE, "% Target"),
          TargetValue         = paste0(MANUAL_TARGET_PREVALENCE, "% prevalence")
        )

      country_target_list[[length(country_target_list) + 1]] <- manual_record
    }

    if (length(country_target_list) > 0) {
      return(do.call(rbind, country_target_list))
    } else {
      return(NULL)
    }
  }

  country_target_list_all[[gender]] <- country_results
}

stopCluster(cl)

# Merge country model targets into the existing target_prevalence_df (which has Global entries)
country_targets <- do.call(rbind, country_target_list_all)
if (!is.null(country_targets) && nrow(country_targets) > 0) {
  rownames(country_targets) <- NULL
  target_prevalence_df <- bind_rows(target_prevalence_df, country_targets)
  cat(sprintf("  Added %d country model target entries to target_prevalence_df\n", nrow(country_targets)))
} else {
  cat("  WARNING: No country model target prevalences computed\n")
}

# Save combined targets
write.csv(target_prevalence_df, file = "results/target_prevalences_dual_evaluation.csv", row.names = FALSE)
cat("Target prevalence calculation complete (both models)\n")

#########################################################################################
#           WEIGHTED PREVALENCE TRENDS (COUNTRY-SPECIFIC MODEL)
#
#   Weighted prevalence trends for country-specific models
#   Calculates population-weighted prevalence and target achievement probabilities
#   using the country-specific model predictions saved to disk.
#
#   Requires: target_prevalence_df (from 06_run_global_model.R),
#             weights_cleaned (from 04_utils.R),
#             calculate_weighted_prevalence (from 04_utils.R)
#   Outputs: final_weighted_results_country
#########################################################################################

cat("\nCalculating weighted prevalence trends for country-specific models...\n")

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = RANDOM_SEED)
clusterExport(cl, "calculate_weighted_prevalence")

gender_weighted_results_country <- list()

for (gender in genders) {
  cat(sprintf("Processing weighted calculations (Country-specific model): %s\n", gender))

  country_dirs <- list.dirs(
    file.path("results/country_specific_ac_nested", gender),
    full.names = FALSE, recursive = FALSE
  )

  weighted_results <- foreach(
    country = country_dirs,
    .packages = c("dplyr", "tidyr"),
    .combine = rbind,
    .errorhandling = 'remove'
  ) %dopar% {

    country_path <- file.path("results/country_specific_ac_nested", gender, country)
    def_types <- list.dirs(country_path, full.names = FALSE, recursive = FALSE)
    def_types <- def_types[grepl("(cigarettes|tobacco)", def_types)]

    country_results <- list()

    country_code <- names(country_name_mapping)[country_name_mapping == country]
    if (length(country_code) == 0) country_code <- country

    for (def_type in def_types) {
      year_files <- list.files(file.path(country_path, def_type), pattern = "\\.rds$")
      years <- as.numeric(sub("\\.rds$", "", year_files))

      target_records <- target_prevalence_df %>%
        filter(Sex == gender, Country == country_code, Def_Type_Code == def_type,
               Model_Type == "Country")

      if (nrow(target_records) == 0) next

      for (year_filtered in years) {
        year_weights <- weights_cleaned %>%
          filter(area == country_code, sex == gender, year == year_filtered) %>%
          arrange(age)

        if (nrow(year_weights) == 0) next

        iterations_matrix_logit <- readRDS(
          file.path(country_path, def_type, paste0(year_filtered, ".rds"))
        )

        if (nrow(iterations_matrix_logit) != nrow(year_weights)) next

        weighted_probs <- calculate_weighted_prevalence(
          iterations_matrix_logit,
          year_weights$weight
        )

        weighted_mean     <- mean(weighted_probs)
        weighted_lower_ci <- quantile(weighted_probs, 0.025)
        weighted_upper_ci <- quantile(weighted_probs, 0.975)

        for (i in 1:nrow(target_records)) {
          target_record <- target_records[i, ]

          prob_achieving_target <- mean(weighted_probs < target_record$TargetPrevalence)

          result_record <- data.frame(
            Country              = country_code,
            Sex                  = gender,
            Year                 = year_filtered,
            Def_Type_Code        = def_type,
            weighted_mean        = weighted_mean,
            weighted_lower_ci    = weighted_lower_ci,
            weighted_upper_ci    = weighted_upper_ci,
            target_prevalence    = target_record$TargetPrevalence,
            prob_achieving_target = prob_achieving_target,
            target_type          = target_record$TargetType,
            target_value         = target_record$TargetValue,
            reduction_percentage = target_record$ReductionPercentage,
            Model_Type           = "Country",
            stringsAsFactors     = FALSE
          )

          country_results[[length(country_results) + 1]] <- result_record
        }
      }
    }

    if (length(country_results) > 0) {
      return(do.call(rbind, country_results))
    } else {
      return(NULL)
    }
  }

  gender_weighted_results_country[[gender]] <- weighted_results
}

stopCluster(cl)

final_weighted_results_country <- do.call(rbind, gender_weighted_results_country)
rownames(final_weighted_results_country) <- NULL

write.csv(
  final_weighted_results_country,
  file = "results/final_weighted_results_country_model.csv",
  row.names = FALSE
)

cat("Country-specific model weighted calculations complete\n")
