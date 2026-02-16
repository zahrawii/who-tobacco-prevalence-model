#########################################################################################
#
#                    WHO TOBACCO CONTROL PREVALENCE PROJECTION MODEL
#                      06_run_global_model.R - Global Model Fitting
#                                   VERSION 2.3.2
#
#   Contains: Sections 8-10 from pipeline_monolith.R
#     - Global hierarchical model fitting for males/females
#     - Country-specific prior extraction
#     - Prediction generation (with all OPT optimizations)
#     - Target prevalence calculation
#     - Weighted prevalence trends
#
#   Requires: All previous modules (00-05) to be sourced
#   Outputs: final_ac_predictions, final_weighted_results_global, gender_results
#
#   EXTRACTED FROM: pipeline_monolith.R v2.3.2
#
#########################################################################################

cat("\n================================================================\n")
cat("  GLOBAL MODEL FITTING MODULE\n")
cat("================================================================\n")

# ---- Initialize Gender Loop ----
genders <- c("males", "females")
gender_results <- list()

for (gender in genders) {

  cat("\n")
  cat("================================================================\n")
  cat("  FITTING GLOBAL HIERARCHICAL MODEL (NIMBLE):", toupper(gender), "\n")
  cat("================================================================\n")

  # ---- 6.1 Setup Directories ----

  gender_dir <- file.path("processing", gender)
  dir.create(gender_dir, showWarnings = FALSE)

  gender_priors_dir <- file.path("country_priors", gender)
  dir.create(gender_priors_dir, showWarnings = FALSE)

  # ---- 6.2 Filter Gender Data ----

  gender_data <- clean_data %>% filter(sex == gender)

  # Defensive check: Skip if no data for this gender
  if (nrow(gender_data) == 0) {
    cat(sprintf("  WARNING: No data for %s - skipping this gender\n", gender))
    next
  }

  # ---- 6.3 Create Region Mappings ----

  unique_regions <- sort(unique(country_region_mapping$region_consolidated))
  region_to_num <- setNames(seq_along(unique_regions), unique_regions)

  country_region_mapping_with_nums <- country_region_mapping %>%
    mutate(Region_Num = region_to_num[region_consolidated])

  gender_data <- gender_data %>%
    left_join(country_region_mapping_with_nums, by = "wb_country_abv") %>%
    mutate(
      Num_Country = as.numeric(factor(wb_country_abv)),
      Num_Survey  = as.numeric(factor(survey)),
      Num_Year    = as.numeric(year)
    )

  country_lookup_df <- gender_data %>%
    select(wb_country_abv, Num_Country, Region_Num, region_consolidated) %>%
    distinct() %>%
    arrange(Num_Country)

  country_to_region <- country_lookup_df$Region_Num
  names(country_to_region) <- country_lookup_df$Num_Country

  country_mapping <- setNames(
    unique(gender_data$wb_country_abv),
    unique(gender_data$Num_Country)
  )

  # ---- 6.4 Calculate Observation Weights ----

  gender_data <- gender_data %>%
    mutate(
      age_range = end_age - start_age,
      weight    = 1 / (age_range + 1),
      weight    = weight / mean(weight)
    )

  # ---- 6.5 Compute Empirical Mean for Informative Prior ----

  empirical_mean_cig <- mean(gender_data$prevalence)

  cat(sprintf("\n  Empirical mean prevalence (logit scale): %.3f\n", empirical_mean_cig))
  cat(sprintf("  Equivalent to approximately %.1f%% prevalence\n", 100 * plogis(empirical_mean_cig)))

  if (empirical_mean_cig < -4 || empirical_mean_cig > 0) {
    warning("  WARNING: Empirical mean outside expected range [-4, 0] - check data transformation")
  }

  # ---- 6.6 Prepare NIMBLE Constants ----
  # CRITICAL FIX (v2.3.1): Country_Region must be length nCountry, not length N

  country_region_for_model <- country_lookup_df %>%
    arrange(Num_Country) %>%
    pull(Region_Num)

  nimble_constants <- list(
    N = nrow(gender_data),
    nCountry = length(unique(gender_data$Num_Country)),
    nRegion = length(unique_regions),
    nSurvey = length(unique(gender_data$Num_Survey)),
    nAgeSpline = ncol(gender_data[, grep("^age_spline_", names(gender_data))]),
    nCohortSpline = ncol(gender_data[, grep("^cohort_spline_", names(gender_data))]),
    nAgeXCohortSplines = ncol(gender_data[, grep("^age_cohort_", names(gender_data))]),
    Country = as.integer(gender_data$Num_Country),
    Country_Region = as.integer(country_region_for_model),
    Survey = as.integer(gender_data$Num_Survey),
    empirical_mean_cig = empirical_mean_cig
  )

  cat(sprintf("  Country_Region length: %d (should equal nCountry: %d)\n",
              length(nimble_constants$Country_Region),
              nimble_constants$nCountry))

  # ---- 6.7 Prepare NIMBLE Data ----

  nimble_data <- list(
    Prevalence = gender_data$prevalence,
    Def_Code_Binary = gender_data$def_code_binary,
    Type_Cig = gender_data$Type_Cig,
    Type_Smoked = gender_data$Type_Smoked,
    Type_Any = gender_data$Type_Any,
    age_spline_matrix = as.matrix(gender_data[, grep("^age_spline_", names(gender_data))]),
    cohort_spline_matrix = as.matrix(gender_data[, grep("^cohort_spline_", names(gender_data))]),
    age_cohort_interaction_matrix = as.matrix(gender_data[, grep("^age_cohort_", names(gender_data))]),
    age_linear_smooth = gender_data$age_linear_smooth,
    spline_weight_var = gender_data$spline_weight_var,
    linear_weight_var = gender_data$linear_weight_var,
    weight = gender_data$weight
  )

  # ---- 6.8 Save Regional Structure ----

  regional_info <- list(
    country_to_region    = country_to_region,
    n_regions            = length(unique_regions),
    region_names         = unique_regions,
    country_lookup_df    = country_lookup_df,
    age_spline_knots     = age_spline_knots_attr,
    age_spline_boundary  = age_spline_boundary_attr,
    cohort_spline_knots  = cohort_spline_knots_attr,
    cohort_spline_boundary = cohort_spline_boundary_attr,
    age_linear_center_constant = AGE_LINEAR_CENTER_CONSTANT,
    cohort_center_constant     = COHORT_CENTER_CONSTANT
  )
  saveRDS(regional_info, file.path(gender_dir, "regional_structure.rds"), compress = "xz")

  # ---- 6.9 Generate Initial Values ----

  generate_global_inits <- function(seed_offset = 0, emp_mean = -1.5) {
    set.seed(RANDOM_SEED + seed_offset)
    list(
      cig_global_intercept = emp_mean + rnorm(1, 0, 0.05),
      smkextra_global_intercept = -3.0 + rnorm(1, 0, 0.05),
      anyextra_global_intercept = -4.0 + rnorm(1, 0, 0.05),
      cig_def_code_shared = rnorm(1, 0.3, 0.05),
      residual_sd = 0.7 + runif(1, -0.05, 0.05),
      cig_age_linear_smooth_effect = -0.02 + rnorm(1, 0, 0.005),
      smkextra_age_linear_smooth_effect = -0.02 + rnorm(1, 0, 0.005),
      anyextra_age_linear_smooth_effect = -0.02 + rnorm(1, 0, 0.005),
      intercept_between_region_precision = 4 + rnorm(1, 0, 0.3),
      age_spline_between_region_precision = 4 + rnorm(1, 0, 0.3),
      cohort_spline_between_region_precision = 4 + rnorm(1, 0, 0.3),
      survey_intercept_precision = 3 + rnorm(1, 0, 0.3),
      smkextra_intercept_within_region_precision = 4 + rnorm(1, 0, 0.3),
      anyextra_intercept_within_region_precision = 4 + rnorm(1, 0, 0.3)
    )
  }

  inits_list <- lapply(1:NUMBER_OF_CHAINS, function(i) {
    generate_global_inits(i, emp_mean = empirical_mean_cig)
  })

  # ---- 6.10 Build NIMBLE Model ----

  cat("  Building NIMBLE model...\n")
  start_build <- Sys.time()

  nimble_model <- nimbleModel(
    code = regional_hierarchical_global_ac_model_nimble,
    constants = nimble_constants,
    data = nimble_data,
    inits = inits_list[[1]],
    name = paste0("GlobalTobaccoModel_", gender)
  )

  build_time <- difftime(Sys.time(), start_build, units = "secs")
  cat(sprintf("  Model build time: %.1f seconds\n", build_time))

  # ---- 6.11 Configure MCMC ----

  cat("  Configuring MCMC...\n")

  mcmc_config <- configureMCMC(
    nimble_model,
    monitors = c(
      "cig_global_intercept", "cig_def_code_shared",
      "smkextra_global_intercept", "anyextra_global_intercept",
      "cig_region_intercept", "smkextra_region_intercept", "anyextra_region_intercept",
      "cig_age_spline_region_mean", "smkextra_age_spline_region_mean", "anyextra_age_spline_region_mean",
      "cig_cohort_spline_region_mean", "smkextra_cohort_spline_region_mean", "anyextra_cohort_spline_region_mean",
      "cig_country_intercept", "smkextra_country_intercept", "anyextra_country_intercept",
      "cig_age_spline", "cig_cohort_spline",
      "cig_age_cohort_interaction",
      "cig_age_linear_smooth_effect", "smkextra_age_linear_smooth_effect", "anyextra_age_linear_smooth_effect",
      "cig_age_spline_global_mean", "cig_cohort_spline_global_mean",
      "smkextra_age_spline_global_mean", "smkextra_cohort_spline_global_mean",
      "anyextra_age_spline_global_mean", "anyextra_cohort_spline_global_mean"
    ),
    thin = THINNING_INTERVAL,
    enableWAIC = FALSE,
    useConjugacy = FALSE
  )

  # ---- 6.12 Build and Compile MCMC ----

  cat("  Building MCMC...\n")
  mcmc_built <- buildMCMC(mcmc_config)

  cat("  Compiling model and MCMC (this may take several minutes)...\n")
  start_compile <- Sys.time()

  compiled_model <- compileNimble(nimble_model)
  compiled_mcmc <- compileNimble(mcmc_built, project = nimble_model)

  compile_time <- difftime(Sys.time(), start_compile, units = "mins")
  cat(sprintf("  Compilation time: %.1f minutes\n", compile_time))

  saveRDS(
    list(compile_time = compile_time, n_params = length(nimble_model$getNodeNames())),
    file.path("compiled_models", paste0("global_model_info_", gender, ".rds"))
  )

  # ---- 6.13 Run MCMC ----

  cat("  Running MCMC...\n")
  start_mcmc <- Sys.time()

  samples <- runMCMC(
    compiled_mcmc,
    niter = NUMBER_OF_BURN + NUMBER_OF_ITERATIONS,
    nburnin = NUMBER_OF_BURN,
    nchains = NUMBER_OF_CHAINS,
    inits = inits_list,
    thin = THINNING_INTERVAL,
    samplesAsCodaMCMC = TRUE,
    progressBar = TRUE,
    summary = TRUE
  )

  mcmc_time <- difftime(Sys.time(), start_mcmc, units = "mins")
  cat(sprintf("  MCMC time: %.1f minutes\n", mcmc_time))

  # ---- 6.14 Convergence Diagnostics ----

  combined_samples_matrix <- do.call(rbind, lapply(samples$samples, as.matrix))

  if (exists("generate_mcmc_diagnostics") && DIAGNOSTICS_AVAILABLE) {
    cat("  Running comprehensive convergence diagnostics...\n")

    global_diagnostics <- tryCatch({
      generate_mcmc_diagnostics(
        samples = samples,
        model_type = "global",
        model_name = "global",
        gender = gender,
        data_summary = gender_data,
        nimble_constants = nimble_constants,
        save_log = TRUE,
        save_tables = TRUE,
        verbose = FALSE
      )
    }, error = function(e) {
      cat(sprintf("  WARNING: Diagnostics failed: %s\n", e$message))
      NULL
    })

    if (!is.null(global_diagnostics)) {
      cat(sprintf("  Convergence Grade: %s\n",
                  global_diagnostics$convergence_summary$convergence_grade))
      cat(sprintf("  Max R-hat: %.3f (%s)\n",
                  global_diagnostics$convergence_summary$max_rhat,
                  global_diagnostics$convergence_summary$worst_param))
    }
  } else {
    gelman_diag <- try(gelman.diag(samples$samples, multivariate = FALSE), silent = TRUE)
    if (!inherits(gelman_diag, "try-error")) {
      max_rhat <- max(gelman_diag$psrf[, "Point est."], na.rm = TRUE)
      cat(sprintf("  Convergence: Max R-hat = %.3f %s\n", max_rhat,
                  ifelse(max_rhat < 1.1, "(Good)", "(Warning)")))
    }
  }

  # [MEM-OPT] Free raw MCMC samples immediately after diagnostics.
  # Only combined_samples_matrix is needed from here on.
  rm(samples)
  gc()

  # ---- 6.15 Extract Country-Specific Priors ----
  # (CORRECTED v2.3.1: Sample-wise computation, no double-counting)

  cat("  Extracting country priors (corrected method)...\n")

  for (country_num in names(country_mapping)) {
    country_code      <- country_mapping[country_num]
    country_full_name <- country_name_mapping[country_code]
    country_region    <- country_lookup_df$Region_Num[country_lookup_df$Num_Country == as.numeric(country_num)]
    region_name       <- unique_regions[country_region]

    # Helper: Compute total intercept SAMPLE-WISE (preserves correlations)
    compute_total_intercept_samplewise <- function(head_prefix) {
      global_col  <- paste0(head_prefix, "_global_intercept")
      region_col  <- paste0(head_prefix, "_region_intercept[", country_region, "]")
      country_col <- paste0(head_prefix, "_country_intercept[", country_num, "]")

      required_cols <- c(global_col, region_col, country_col)
      missing_cols <- required_cols[!required_cols %in% colnames(combined_samples_matrix)]

      if (length(missing_cols) > 0) {
        warning(sprintf("Missing columns for %s total intercept: %s",
                        head_prefix, paste(missing_cols, collapse = ", ")))
        return(list(mean = NA_real_, sd = NA_real_))
      }

      total_samples <- combined_samples_matrix[, global_col] +
        combined_samples_matrix[, region_col] +
        combined_samples_matrix[, country_col]

      list(mean = mean(total_samples), sd = sd(total_samples))
    }

    # CIG: Total intercept + country-specific splines
    cig_intercept_stats <- compute_total_intercept_samplewise("cig")

    cig_age_spline_means <- cig_age_spline_sds <- c()
    for (l in 1:nimble_constants$nAgeSpline) {
      country_col <- paste0("cig_age_spline[", country_num, ", ", l, "]")
      if (country_col %in% colnames(combined_samples_matrix)) {
        cig_age_spline_means[l] <- mean(combined_samples_matrix[, country_col])
        cig_age_spline_sds[l]   <- sd(combined_samples_matrix[, country_col])
      } else {
        regional_col <- paste0("cig_age_spline_region_mean[", country_region, ", ", l, "]")
        if (regional_col %in% colnames(combined_samples_matrix)) {
          cig_age_spline_means[l] <- mean(combined_samples_matrix[, regional_col])
          cig_age_spline_sds[l]   <- sd(combined_samples_matrix[, regional_col])
        } else {
          cig_age_spline_means[l] <- 0
          cig_age_spline_sds[l]   <- 1
        }
      }
    }

    cig_cohort_spline_means <- cig_cohort_spline_sds <- c()
    for (m in 1:nimble_constants$nCohortSpline) {
      country_col <- paste0("cig_cohort_spline[", country_num, ", ", m, "]")
      if (country_col %in% colnames(combined_samples_matrix)) {
        cig_cohort_spline_means[m] <- mean(combined_samples_matrix[, country_col])
        cig_cohort_spline_sds[m]   <- sd(combined_samples_matrix[, country_col])
      } else {
        regional_col <- paste0("cig_cohort_spline_region_mean[", country_region, ", ", m, "]")
        if (regional_col %in% colnames(combined_samples_matrix)) {
          cig_cohort_spline_means[m] <- mean(combined_samples_matrix[, regional_col])
          cig_cohort_spline_sds[m]   <- sd(combined_samples_matrix[, regional_col])
        } else {
          cig_cohort_spline_means[m] <- 0
          cig_cohort_spline_sds[m]   <- 1
        }
      }
    }

    # SMKEXTRA: Total intercept + regional splines
    smkextra_intercept_stats <- compute_total_intercept_samplewise("smkextra")

    smkextra_age_spline_means <- smkextra_age_spline_sds <- c()
    for (l in 1:nimble_constants$nAgeSpline) {
      col_name <- paste0("smkextra_age_spline_region_mean[", country_region, ", ", l, "]")
      if (col_name %in% colnames(combined_samples_matrix)) {
        smkextra_age_spline_means[l] <- mean(combined_samples_matrix[, col_name])
        smkextra_age_spline_sds[l]   <- sd(combined_samples_matrix[, col_name])
      } else {
        smkextra_age_spline_means[l] <- 0
        smkextra_age_spline_sds[l]   <- 0.5
      }
    }

    smkextra_cohort_spline_means <- smkextra_cohort_spline_sds <- c()
    for (m in 1:nimble_constants$nCohortSpline) {
      col_name <- paste0("smkextra_cohort_spline_region_mean[", country_region, ", ", m, "]")
      if (col_name %in% colnames(combined_samples_matrix)) {
        smkextra_cohort_spline_means[m] <- mean(combined_samples_matrix[, col_name])
        smkextra_cohort_spline_sds[m]   <- sd(combined_samples_matrix[, col_name])
      } else {
        smkextra_cohort_spline_means[m] <- 0
        smkextra_cohort_spline_sds[m]   <- 0.5
      }
    }

    # ANYEXTRA: Total intercept + regional splines
    anyextra_intercept_stats <- compute_total_intercept_samplewise("anyextra")

    anyextra_age_spline_means <- anyextra_age_spline_sds <- c()
    for (l in 1:nimble_constants$nAgeSpline) {
      col_name <- paste0("anyextra_age_spline_region_mean[", country_region, ", ", l, "]")
      if (col_name %in% colnames(combined_samples_matrix)) {
        anyextra_age_spline_means[l] <- mean(combined_samples_matrix[, col_name])
        anyextra_age_spline_sds[l]   <- sd(combined_samples_matrix[, col_name])
      } else {
        anyextra_age_spline_means[l] <- 0
        anyextra_age_spline_sds[l]   <- 0.5
      }
    }

    anyextra_cohort_spline_means <- anyextra_cohort_spline_sds <- c()
    for (m in 1:nimble_constants$nCohortSpline) {
      col_name <- paste0("anyextra_cohort_spline_region_mean[", country_region, ", ", m, "]")
      if (col_name %in% colnames(combined_samples_matrix)) {
        anyextra_cohort_spline_means[m] <- mean(combined_samples_matrix[, col_name])
        anyextra_cohort_spline_sds[m]   <- sd(combined_samples_matrix[, col_name])
      } else {
        anyextra_cohort_spline_means[m] <- 0
        anyextra_cohort_spline_sds[m]   <- 0.5
      }
    }

    # Age-cohort interactions (global, same for all countries)
    cig_age_cohort_cols <- grep("^cig_age_cohort_interaction\\[", colnames(combined_samples_matrix))
    if (length(cig_age_cohort_cols) > 0) {
      cig_age_cohort_means <- colMeans(combined_samples_matrix[, cig_age_cohort_cols, drop = FALSE])
      cig_age_cohort_sds   <- apply(combined_samples_matrix[, cig_age_cohort_cols, drop = FALSE], 2, sd)
    } else {
      n_interactions <- nimble_constants$nAgeXCohortSplines
      cig_age_cohort_means <- rep(0, n_interactions)
      cig_age_cohort_sds   <- rep(0.1, n_interactions)
    }

    # Global parameters
    cig_age_linear_mean <- mean(combined_samples_matrix[, "cig_age_linear_smooth_effect"])
    cig_age_linear_sd   <- sd(combined_samples_matrix[, "cig_age_linear_smooth_effect"])
    smkextra_age_linear_mean <- mean(combined_samples_matrix[, "smkextra_age_linear_smooth_effect"])
    smkextra_age_linear_sd   <- sd(combined_samples_matrix[, "smkextra_age_linear_smooth_effect"])
    anyextra_age_linear_mean <- mean(combined_samples_matrix[, "anyextra_age_linear_smooth_effect"])
    anyextra_age_linear_sd   <- sd(combined_samples_matrix[, "anyextra_age_linear_smooth_effect"])
    def_code_mean <- mean(combined_samples_matrix[, "cig_def_code_shared"])
    def_code_sd   <- sd(combined_samples_matrix[, "cig_def_code_shared"])

    # Build priors dataframe with shrinkage
    shrinkage <- 0.7
    min_sd_intercept <- 0.05; min_sd_spline <- 0.02
    min_sd_interaction <- 0.01; min_sd_linear <- 0.002

    all_priors <- data.frame(
      parameter = c(
        "cig_def_code_shared", "cig_intercept", "cig_age_linear_smooth_effect",
        paste0("cig_age_spline_", 1:length(cig_age_spline_means)),
        paste0("cig_cohort_spline_", 1:length(cig_cohort_spline_means)),
        paste0("cig_age_cohort_interaction_", 1:length(cig_age_cohort_means)),
        "smkextra_intercept", "smkextra_age_linear_smooth_effect",
        paste0("smkextra_age_spline_", 1:length(smkextra_age_spline_means)),
        paste0("smkextra_cohort_spline_", 1:length(smkextra_cohort_spline_means)),
        "anyextra_intercept", "anyextra_age_linear_smooth_effect",
        paste0("anyextra_age_spline_", 1:length(anyextra_age_spline_means)),
        paste0("anyextra_cohort_spline_", 1:length(anyextra_cohort_spline_means))
      ),
      mean = c(
        def_code_mean, cig_intercept_stats$mean, cig_age_linear_mean,
        cig_age_spline_means, cig_cohort_spline_means, cig_age_cohort_means,
        smkextra_intercept_stats$mean, smkextra_age_linear_mean,
        smkextra_age_spline_means, smkextra_cohort_spline_means,
        anyextra_intercept_stats$mean, anyextra_age_linear_mean,
        anyextra_age_spline_means, anyextra_cohort_spline_means
      ),
      sd = c(
        pmax(def_code_sd * shrinkage, min_sd_intercept),
        pmax(cig_intercept_stats$sd * shrinkage, min_sd_intercept),
        pmax(cig_age_linear_sd * shrinkage, min_sd_linear),
        pmax(cig_age_spline_sds * shrinkage, min_sd_spline),
        pmax(cig_cohort_spline_sds * shrinkage, min_sd_spline),
        pmax(cig_age_cohort_sds * shrinkage, min_sd_interaction),
        pmax(smkextra_intercept_stats$sd * shrinkage, min_sd_intercept),
        pmax(smkextra_age_linear_sd * shrinkage, min_sd_linear),
        pmax(smkextra_age_spline_sds * shrinkage, min_sd_spline),
        pmax(smkextra_cohort_spline_sds * shrinkage, min_sd_spline),
        pmax(anyextra_intercept_stats$sd * shrinkage, min_sd_intercept),
        pmax(anyextra_age_linear_sd * shrinkage, min_sd_linear),
        pmax(anyextra_age_spline_sds * shrinkage, min_sd_spline),
        pmax(anyextra_cohort_spline_sds * shrinkage, min_sd_spline)
      ),
      stringsAsFactors = FALSE
    )

    # Handle NAs
    na_means <- is.na(all_priors$mean)
    na_sds   <- is.na(all_priors$sd)
    if (any(na_means) || any(na_sds)) {
      all_priors$mean[na_means] <- 0
      all_priors$sd[na_sds]     <- 1
    }

    write.csv(
      all_priors,
      file = file.path(gender_priors_dir, paste0(country_code, "_regional_ac_priors_nested.csv")),
      row.names = FALSE
    )
  }

  cat("  Country priors saved.\n")

  # ---- 6.16 Generate Predictions ----

  n_samples <- min(1000, nrow(combined_samples_matrix))
  sampled_indices <- sort(sample(nrow(combined_samples_matrix), n_samples))

  cig_samples <- list(
    global_intercept  = combined_samples_matrix[sampled_indices, "cig_global_intercept"],
    region_intercept  = combined_samples_matrix[sampled_indices, grep("^cig_region_intercept\\[", colnames(combined_samples_matrix))],
    age_linear_smooth = combined_samples_matrix[sampled_indices, "cig_age_linear_smooth_effect"],
    country_intercept = combined_samples_matrix[sampled_indices, grep("^cig_country_intercept\\[", colnames(combined_samples_matrix))],
    age_spline        = combined_samples_matrix[sampled_indices, grep("^cig_age_spline\\[", colnames(combined_samples_matrix))],
    cohort_spline     = combined_samples_matrix[sampled_indices, grep("^cig_cohort_spline\\[", colnames(combined_samples_matrix))],
    age_cohort        = combined_samples_matrix[sampled_indices, grep("^cig_age_cohort_interaction\\[", colnames(combined_samples_matrix))]
  )

  smkextra_samples <- list(
    global_intercept     = combined_samples_matrix[sampled_indices, "smkextra_global_intercept"],
    region_intercept     = combined_samples_matrix[sampled_indices, grep("^smkextra_region_intercept\\[", colnames(combined_samples_matrix))],
    age_linear_smooth    = combined_samples_matrix[sampled_indices, "smkextra_age_linear_smooth_effect"],
    country_intercept    = combined_samples_matrix[sampled_indices, grep("^smkextra_country_intercept\\[", colnames(combined_samples_matrix))],
    age_spline_regional  = combined_samples_matrix[sampled_indices, grep("^smkextra_age_spline_region_mean\\[", colnames(combined_samples_matrix))],
    cohort_spline_regional = combined_samples_matrix[sampled_indices, grep("^smkextra_cohort_spline_region_mean\\[", colnames(combined_samples_matrix))]
  )

  anyextra_samples <- list(
    global_intercept     = combined_samples_matrix[sampled_indices, "anyextra_global_intercept"],
    region_intercept     = combined_samples_matrix[sampled_indices, grep("^anyextra_region_intercept\\[", colnames(combined_samples_matrix))],
    age_linear_smooth    = combined_samples_matrix[sampled_indices, "anyextra_age_linear_smooth_effect"],
    country_intercept    = combined_samples_matrix[sampled_indices, grep("^anyextra_country_intercept\\[", colnames(combined_samples_matrix))],
    age_spline_regional  = combined_samples_matrix[sampled_indices, grep("^anyextra_age_spline_region_mean\\[", colnames(combined_samples_matrix))],
    cohort_spline_regional = combined_samples_matrix[sampled_indices, grep("^anyextra_cohort_spline_region_mean\\[", colnames(combined_samples_matrix))]
  )

  def_code_shared_samples <- combined_samples_matrix[sampled_indices, "cig_def_code_shared"]

  age_midpoints        <- seq(15, 100, by = 1)
  years_for_prediction <- seq(min(gender_data$Num_Year), target_year + PROJECTED_YEARS, by = 1)
  unique_countries     <- unique(gender_data$Num_Country)

  # Create prediction grid
  new_data_age_period <- expand.grid(
    Year         = years_for_prediction,
    Def_Code_Binary = c(0, 1),
    Age_Midpoint = age_midpoints,
    Num_Country  = unique_countries
  )

  new_data_age_cohort <- new_data_age_period %>%
    mutate(
      Birth_Cohort        = Year - Age_Midpoint,
      Cohort_Centered     = Birth_Cohort - COHORT_CENTER_CONSTANT,
      sigmoid_input       = (Age_Midpoint - TRANSITION_START) / TRANSITION_WIDTH,
      spline_weight       = 1 / (1 + exp(sigmoid_input)),
      linear_weight       = 1 - spline_weight,
      Age_For_Spline      = Age_Midpoint,
      Age_Linear          = pmax(0, Age_Midpoint - TRANSITION_START),
      Age_Linear_Centered = Age_Linear - AGE_LINEAR_CENTER_CONSTANT
    )

  new_age_spline_basis <- ns(
    new_data_age_cohort$Age_For_Spline,
    knots = age_spline_knots_attr,
    Boundary.knots = age_spline_boundary_attr
  )
  new_age_spline_basis_df <- as.data.frame(new_age_spline_basis)
  colnames(new_age_spline_basis_df) <- paste0("age_spline_", 1:ncol(new_age_spline_basis_df))

  new_cohort_spline_basis <- ns(
    new_data_age_cohort$Birth_Cohort,
    knots = cohort_spline_knots_attr,
    Boundary.knots = cohort_spline_boundary_attr
  )
  new_cohort_spline_basis_df <- as.data.frame(new_cohort_spline_basis)
  colnames(new_cohort_spline_basis_df) <- paste0("cohort_spline_", 1:ncol(new_cohort_spline_basis_df))

  new_data_age_cohort$age_linear_smooth   <- new_data_age_cohort$Age_Linear_Centered
  new_data_age_cohort$spline_weight_var   <- new_data_age_cohort$spline_weight
  new_data_age_cohort$linear_weight_var   <- new_data_age_cohort$linear_weight

  new_n_interactions <- ncol(new_age_spline_basis_df) * ncol(new_cohort_spline_basis_df)
  new_age_cohort_interaction_matrix <- matrix(0, nrow = nrow(new_data_age_cohort), ncol = new_n_interactions)

  for (i in 1:nrow(new_data_age_cohort)) {
    age_values    <- as.numeric(new_age_spline_basis_df[i, ])
    cohort_values <- as.numeric(new_cohort_spline_basis_df[i, ])
    new_age_cohort_interaction_matrix[i, ] <- as.vector(outer(age_values, cohort_values))
  }

  new_age_cohort_interaction_df <- as.data.frame(new_age_cohort_interaction_matrix)
  colnames(new_age_cohort_interaction_df) <- paste0("age_cohort_", 1:ncol(new_age_cohort_interaction_df))

  new_data_apc <- cbind(new_data_age_cohort, new_age_spline_basis_df,
                        new_cohort_spline_basis_df, new_age_cohort_interaction_df)

  # ---- 6.17 Parallel Prediction Loop ----

  num_cores <- max(1, min(detectCores() - 2, 4))
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  clusterSetRNGStream(cl, iseed = RANDOM_SEED)
  clusterExport(cl, c("check_ordering_per_draw", "country_to_region"))

  country_results <- foreach(
    country = unique_countries,
    .packages = c("dplyr", "splines"),
    .combine = rbind,
    .errorhandling = 'remove'
  ) %dopar% {

    country_code      <- country_mapping[as.character(country)]
    country_full_name <- country_name_mapping[country_code]
    country_dir       <- file.path(gender_dir, country_full_name)
    dir.create(country_dir, showWarnings = FALSE)

    country_region <- country_to_region[as.character(country)]
    country_data <- new_data_apc[new_data_apc$Num_Country == country, ]
    country_data <- country_data %>%
      mutate(Def_Code = ifelse(Def_Code_Binary == 0, "daily_user", "current_user"))

    results_list    <- list()
    year_def_combos <- unique(country_data[, c("Year", "Def_Code")])

    for (i in 1:nrow(year_def_combos)) {
      current_year <- year_def_combos$Year[i]
      current_def  <- year_def_combos$Def_Code[i]

      cig_dir    <- file.path(country_dir, paste0(current_def, "_cigarettes"))
      smoked_dir <- file.path(country_dir, paste0(current_def, "_any_smoked_tobacco"))
      any_dir    <- file.path(country_dir, paste0(current_def, "_any_tobacco_product"))

      dir.create(cig_dir, showWarnings = FALSE)
      dir.create(smoked_dir, showWarnings = FALSE)
      dir.create(any_dir, showWarnings = FALSE)

      current_data <- country_data %>%
        filter(Year == current_year, Def_Code == current_def)

      predictions_cig    <- matrix(0, nrow = nrow(current_data), ncol = n_samples)
      predictions_smoked <- matrix(0, nrow = nrow(current_data), ncol = n_samples)
      predictions_any    <- matrix(0, nrow = nrow(current_data), ncol = n_samples)

      for (j in 1:nrow(current_data)) {
        def_code_binary               <- current_data$Def_Code_Binary[j]
        age_spline_values             <- as.numeric(current_data[j, grep("^age_spline_", names(current_data))])
        age_linear_smooth             <- current_data$age_linear_smooth[j]
        spline_weight_value           <- current_data$spline_weight_var[j]
        linear_weight_value           <- current_data$linear_weight_var[j]
        cohort_spline_values          <- as.numeric(current_data[j, grep("^cohort_spline_", names(current_data))])
        age_cohort_interaction_values <- as.numeric(current_data[j, grep("^age_cohort_", names(current_data))])

        mu_cig <- cig_samples$global_intercept +
          cig_samples$region_intercept[, paste0("cig_region_intercept[", country_region, "]")] +
          cig_samples$country_intercept[, paste0("cig_country_intercept[", country, "]")] +
          def_code_shared_samples * def_code_binary +
          spline_weight_value * as.matrix(cig_samples$age_spline[, paste0("cig_age_spline[", country, ", ", 1:length(age_spline_values), "]")]) %*% age_spline_values +
          linear_weight_value * cig_samples$age_linear_smooth * age_linear_smooth +
          as.matrix(cig_samples$cohort_spline[, paste0("cig_cohort_spline[", country, ", ", 1:length(cohort_spline_values), "]")]) %*% cohort_spline_values +
          as.matrix(cig_samples$age_cohort) %*% age_cohort_interaction_values

        mu_smkextra <- smkextra_samples$global_intercept +
          smkextra_samples$region_intercept[, paste0("smkextra_region_intercept[", country_region, "]")] +
          smkextra_samples$country_intercept[, paste0("smkextra_country_intercept[", country, "]")] +
          0.3 * def_code_shared_samples * def_code_binary +
          spline_weight_value * as.matrix(smkextra_samples$age_spline_regional[, paste0("smkextra_age_spline_region_mean[", country_region, ", ", 1:length(age_spline_values), "]")]) %*% age_spline_values +
          linear_weight_value * smkextra_samples$age_linear_smooth * age_linear_smooth +
          as.matrix(smkextra_samples$cohort_spline_regional[, paste0("smkextra_cohort_spline_region_mean[", country_region, ", ", 1:length(cohort_spline_values), "]")]) %*% cohort_spline_values

        mu_anyextra <- anyextra_samples$global_intercept +
          anyextra_samples$region_intercept[, paste0("anyextra_region_intercept[", country_region, "]")] +
          anyextra_samples$country_intercept[, paste0("anyextra_country_intercept[", country, "]")] +
          0.3 * def_code_shared_samples * def_code_binary +
          spline_weight_value * as.matrix(anyextra_samples$age_spline_regional[, paste0("anyextra_age_spline_region_mean[", country_region, ", ", 1:length(age_spline_values), "]")]) %*% age_spline_values +
          linear_weight_value * anyextra_samples$age_linear_smooth * age_linear_smooth +
          as.matrix(anyextra_samples$cohort_spline_regional[, paste0("anyextra_cohort_spline_region_mean[", country_region, ", ", 1:length(cohort_spline_values), "]")]) %*% cohort_spline_values

        predictions_cig[j, ] <- mu_cig

        p_cig         <- plogis(mu_cig)
        p_smoked_full <- p_cig + plogis(mu_smkextra) * (1 - p_cig)
        predictions_smoked[j, ] <- log(p_smoked_full / (1 - p_smoked_full))

        p_any_full <- p_smoked_full + plogis(mu_anyextra) * (1 - p_smoked_full)
        predictions_any[j, ] <- log(p_any_full / (1 - p_any_full))
      }

      has_violations <- check_ordering_per_draw(predictions_cig, predictions_smoked, predictions_any)
      if (has_violations) {
        stop("Ordering violations detected")
      }

      saveRDS(predictions_cig, file = file.path(cig_dir, paste0(current_year, ".rds")), compress = TRUE)
      saveRDS(predictions_smoked, file = file.path(smoked_dir, paste0(current_year, ".rds")), compress = TRUE)
      saveRDS(predictions_any, file = file.path(any_dir, paste0(current_year, ".rds")), compress = TRUE)

      # Summarize for results
      def_type_cig <- paste0(current_def, "_cigarettes")
      def_type_smoked <- paste0(current_def, "_any_smoked_tobacco")
      def_type_any <- paste0(current_def, "_any_tobacco_product")

      p_cig_matrix <- plogis(predictions_cig)
      p_smoked_matrix <- plogis(predictions_smoked)
      p_any_matrix <- plogis(predictions_any)

      for (k in 1:nrow(current_data)) {
        base_row <- data.frame(
          Year = current_year,
          Age_Midpoint = current_data$Age_Midpoint[k],
          Birth_Cohort = current_data$Birth_Cohort[k],
          Sex = gender,
          Country = country_code,
          Has_Survey_Data = TRUE,
          stringsAsFactors = FALSE
        )

        results_list[[length(results_list) + 1]] <- cbind(base_row, data.frame(
          Def_Type_Code = def_type_cig,
          Prevalence = mean(p_cig_matrix[k, ]),
          lower_ci = quantile(p_cig_matrix[k, ], 0.025),
          upper_ci = quantile(p_cig_matrix[k, ], 0.975),
          stringsAsFactors = FALSE
        ))

        results_list[[length(results_list) + 1]] <- cbind(base_row, data.frame(
          Def_Type_Code = def_type_smoked,
          Prevalence = mean(p_smoked_matrix[k, ]),
          lower_ci = quantile(p_smoked_matrix[k, ], 0.025),
          upper_ci = quantile(p_smoked_matrix[k, ], 0.975),
          stringsAsFactors = FALSE
        ))

        results_list[[length(results_list) + 1]] <- cbind(base_row, data.frame(
          Def_Type_Code = def_type_any,
          Prevalence = mean(p_any_matrix[k, ]),
          lower_ci = quantile(p_any_matrix[k, ], 0.025),
          upper_ci = quantile(p_any_matrix[k, ], 0.975),
          stringsAsFactors = FALSE
        ))
      }
    }

    if (length(results_list) > 0) {
      return(do.call(rbind, results_list))
    } else {
      return(NULL)
    }
  }

  # ==================================================================
  # HELPER FUNCTIONS FOR NO-DATA COUNTRIES
  # ==================================================================

  sample_cig_params_from_prior <- function(region_num, mcmc_samples, n_samples, nAgeSpline, nCohortSpline) {
    # Sample country intercept deviation from regional variance
    within_sd_col <- grep("intercept_within_region_sd", colnames(mcmc_samples))
    if (length(within_sd_col) > 0) {
      within_sd <- mcmc_samples[, within_sd_col[1]]
    } else {
      within_sd <- rep(0.3, n_samples)
    }
    country_intercept <- rnorm(n_samples, 0, within_sd)

    # Sample age splines from regional means with noise
    age_splines <- matrix(0, n_samples, nAgeSpline)
    for (l in 1:nAgeSpline) {
      regional_col <- paste0("cig_age_spline_region_mean[", region_num, ", ", l, "]")
      if (regional_col %in% colnames(mcmc_samples)) {
        age_splines[, l] <- mcmc_samples[, regional_col] + rnorm(n_samples, 0, 0.1)
      }
    }

    # Sample cohort splines from regional means
    cohort_splines <- matrix(0, n_samples, nCohortSpline)
    for (m in 1:nCohortSpline) {
      regional_col <- paste0("cig_cohort_spline_region_mean[", region_num, ", ", m, "]")
      if (regional_col %in% colnames(mcmc_samples)) {
        cohort_splines[, m] <- mcmc_samples[, regional_col] + rnorm(n_samples, 0, 0.1)
      }
    }

    list(country_intercept = country_intercept, age_splines = age_splines, cohort_splines = cohort_splines)
  }

  sample_extra_intercept_from_prior <- function(head_prefix, mcmc_samples, n_samples) {
    within_sd_col <- grep(paste0(head_prefix, "_intercept_within_region_sd"), colnames(mcmc_samples))
    if (length(within_sd_col) > 0) {
      within_sd <- mcmc_samples[, within_sd_col[1]]
    } else {
      within_sd <- rep(0.5, n_samples)
    }
    rnorm(n_samples, 0, within_sd)
  }

  # ==================================================================
  # ALL-COUNTRIES PREDICTION: Include countries without survey data
  # ==================================================================

  all_countries_master <- unique(country_region_manual$wb_country_abv)
  data_country_codes <- unique(country_lookup_df$wb_country_abv)

  all_country_lookup <- country_region_manual %>%
    distinct(wb_country_abv, region_consolidated) %>%
    mutate(
      has_data = wb_country_abv %in% data_country_codes,
      Region_Num = region_to_num[region_consolidated]
    ) %>%
    left_join(
      country_lookup_df %>% select(wb_country_abv, Num_Country),
      by = "wb_country_abv"
    )

  n_data_countries <- sum(all_country_lookup$has_data)
  n_nodata_countries <- sum(!all_country_lookup$has_data)

  cat(sprintf("\n  === ALL-COUNTRIES PREDICTION ===\n"))
  cat(sprintf("  Total countries in master list: %d\n", nrow(all_country_lookup)))
  cat(sprintf("  Countries WITH survey data:     %d\n", n_data_countries))
  cat(sprintf("  Countries WITHOUT survey data:  %d\n", n_nodata_countries))

  # ==================================================================
  # PREDICTION FOR COUNTRIES WITHOUT SURVEY DATA
  # ==================================================================

  nodata_countries <- all_country_lookup %>%
    filter(!has_data) %>%
    filter(!is.na(Region_Num))

  if (nrow(nodata_countries) > 0) {

    cat(sprintf("\n  Generating predictions for %d countries WITHOUT survey data...\n",
                nrow(nodata_countries)))

    nodata_prediction_grid <- expand.grid(
      Year = years_for_prediction,
      Def_Code_Binary = c(0, 1),
      Age_Midpoint = age_midpoints
    )

    nodata_prediction_grid <- nodata_prediction_grid %>%
      mutate(
        Birth_Cohort = Year - Age_Midpoint,
        Cohort_Centered = Birth_Cohort - COHORT_CENTER_CONSTANT,
        sigmoid_input = (Age_Midpoint - TRANSITION_START) / TRANSITION_WIDTH,
        spline_weight = 1 / (1 + exp(sigmoid_input)),
        linear_weight = 1 - spline_weight,
        Age_For_Spline = Age_Midpoint,
        Age_Linear = pmax(0, Age_Midpoint - TRANSITION_START),
        Age_Linear_Centered = Age_Linear - AGE_LINEAR_CENTER_CONSTANT
      )

    nodata_age_spline_basis <- ns(
      nodata_prediction_grid$Age_For_Spline,
      knots = age_spline_knots_attr,
      Boundary.knots = age_spline_boundary_attr
    )
    nodata_age_spline_df <- as.data.frame(nodata_age_spline_basis)
    colnames(nodata_age_spline_df) <- paste0("age_spline_", 1:ncol(nodata_age_spline_df))

    nodata_cohort_spline_basis <- ns(
      nodata_prediction_grid$Birth_Cohort,
      knots = cohort_spline_knots_attr,
      Boundary.knots = cohort_spline_boundary_attr
    )
    nodata_cohort_spline_df <- as.data.frame(nodata_cohort_spline_basis)
    colnames(nodata_cohort_spline_df) <- paste0("cohort_spline_", 1:ncol(nodata_cohort_spline_df))

    nodata_prediction_grid$age_linear_smooth <- nodata_prediction_grid$Age_Linear_Centered
    nodata_prediction_grid$spline_weight_var <- nodata_prediction_grid$spline_weight
    nodata_prediction_grid$linear_weight_var <- nodata_prediction_grid$linear_weight

    nodata_n_interactions <- ncol(nodata_age_spline_df) * ncol(nodata_cohort_spline_df)
    nodata_age_cohort_matrix <- matrix(0, nrow = nrow(nodata_prediction_grid), ncol = nodata_n_interactions)

    for (i in 1:nrow(nodata_prediction_grid)) {
      age_vals <- as.numeric(nodata_age_spline_df[i, ])
      cohort_vals <- as.numeric(nodata_cohort_spline_df[i, ])
      nodata_age_cohort_matrix[i, ] <- as.vector(outer(age_vals, cohort_vals))
    }

    nodata_age_cohort_df <- as.data.frame(nodata_age_cohort_matrix)
    colnames(nodata_age_cohort_df) <- paste0("age_cohort_", 1:ncol(nodata_age_cohort_df))

    nodata_grid_full <- cbind(nodata_prediction_grid, nodata_age_spline_df,
                              nodata_cohort_spline_df, nodata_age_cohort_df)

    mcmc_samples_subsetted <- combined_samples_matrix[sampled_indices, ]

    clusterExport(cl, c(
      "all_country_lookup", "nodata_grid_full", "country_name_mapping",
      "sample_cig_params_from_prior", "sample_extra_intercept_from_prior",
      "cig_samples", "smkextra_samples", "anyextra_samples",
      "def_code_shared_samples", "mcmc_samples_subsetted", "n_samples",
      "nimble_constants", "gender_dir", "check_ordering_per_draw"
    ), envir = environment())

    nodata_results <- foreach(
      country_code = nodata_countries$wb_country_abv,
      .packages = c("dplyr", "splines"),
      .combine = rbind,
      .errorhandling = 'remove'
    ) %dopar% {

      country_info <- all_country_lookup[all_country_lookup$wb_country_abv == country_code, ]
      country_region <- country_info$Region_Num

      if (is.na(country_region) || length(country_region) == 0) {
        return(NULL)
      }

      country_full_name <- country_name_mapping[country_code]
      if (is.na(country_full_name) || is.null(country_full_name)) {
        country_full_name <- toupper(country_code)
      }

      country_dir <- file.path(gender_dir, country_full_name)
      dir.create(country_dir, showWarnings = FALSE, recursive = TRUE)

      saveRDS(
        list(has_survey_data = FALSE, region = country_info$region_consolidated),
        file = file.path(country_dir, "_metadata.rds")
      )

      cig_prior_params <- sample_cig_params_from_prior(
        region_num = country_region,
        mcmc_samples = mcmc_samples_subsetted,
        n_samples = n_samples,
        nAgeSpline = nimble_constants$nAgeSpline,
        nCohortSpline = nimble_constants$nCohortSpline
      )
      cig_country_intercept <- cig_prior_params$country_intercept
      cig_age_spline_sampled <- cig_prior_params$age_splines
      cig_cohort_spline_sampled <- cig_prior_params$cohort_splines

      smkextra_country_intercept <- sample_extra_intercept_from_prior("smkextra", mcmc_samples_subsetted, n_samples)
      anyextra_country_intercept <- sample_extra_intercept_from_prior("anyextra", mcmc_samples_subsetted, n_samples)

      country_data <- nodata_grid_full %>%
        mutate(Def_Code = ifelse(Def_Code_Binary == 0, "daily_user", "current_user"))

      results_list <- list()
      year_def_combos <- unique(country_data[, c("Year", "Def_Code")])

      for (i in 1:nrow(year_def_combos)) {
        current_year <- year_def_combos$Year[i]
        current_def <- year_def_combos$Def_Code[i]

        cig_dir <- file.path(country_dir, paste0(current_def, "_cigarettes"))
        smoked_dir <- file.path(country_dir, paste0(current_def, "_any_smoked_tobacco"))
        any_dir <- file.path(country_dir, paste0(current_def, "_any_tobacco_product"))

        dir.create(cig_dir, showWarnings = FALSE)
        dir.create(smoked_dir, showWarnings = FALSE)
        dir.create(any_dir, showWarnings = FALSE)

        current_data <- country_data %>%
          filter(Year == current_year, Def_Code == current_def)

        predictions_cig <- matrix(0, nrow = nrow(current_data), ncol = n_samples)
        predictions_smoked <- matrix(0, nrow = nrow(current_data), ncol = n_samples)
        predictions_any <- matrix(0, nrow = nrow(current_data), ncol = n_samples)

        for (j in 1:nrow(current_data)) {
          def_code_binary <- current_data$Def_Code_Binary[j]
          age_spline_values <- as.numeric(current_data[j, grep("^age_spline_", names(current_data))])
          age_linear_smooth <- current_data$age_linear_smooth[j]
          spline_weight_value <- current_data$spline_weight_var[j]
          linear_weight_value <- current_data$linear_weight_var[j]
          cohort_spline_values <- as.numeric(current_data[j, grep("^cohort_spline_", names(current_data))])
          age_cohort_interaction_values <- as.numeric(current_data[j, grep("^age_cohort_", names(current_data))])

          mu_cig <- cig_samples$global_intercept +
            cig_samples$region_intercept[, paste0("cig_region_intercept[", country_region, "]")] +
            cig_country_intercept +
            def_code_shared_samples * def_code_binary +
            spline_weight_value * as.matrix(cig_age_spline_sampled) %*% age_spline_values +
            linear_weight_value * cig_samples$age_linear_smooth * age_linear_smooth +
            as.matrix(cig_cohort_spline_sampled) %*% cohort_spline_values +
            as.matrix(cig_samples$age_cohort) %*% age_cohort_interaction_values

          mu_smkextra <- smkextra_samples$global_intercept +
            smkextra_samples$region_intercept[, paste0("smkextra_region_intercept[", country_region, "]")] +
            smkextra_country_intercept +
            0.3 * def_code_shared_samples * def_code_binary +
            spline_weight_value * as.matrix(smkextra_samples$age_spline_regional[, paste0("smkextra_age_spline_region_mean[", country_region, ", ", 1:length(age_spline_values), "]")]) %*% age_spline_values +
            linear_weight_value * smkextra_samples$age_linear_smooth * age_linear_smooth +
            as.matrix(smkextra_samples$cohort_spline_regional[, paste0("smkextra_cohort_spline_region_mean[", country_region, ", ", 1:length(cohort_spline_values), "]")]) %*% cohort_spline_values

          mu_anyextra <- anyextra_samples$global_intercept +
            anyextra_samples$region_intercept[, paste0("anyextra_region_intercept[", country_region, "]")] +
            anyextra_country_intercept +
            0.3 * def_code_shared_samples * def_code_binary +
            spline_weight_value * as.matrix(anyextra_samples$age_spline_regional[, paste0("anyextra_age_spline_region_mean[", country_region, ", ", 1:length(age_spline_values), "]")]) %*% age_spline_values +
            linear_weight_value * anyextra_samples$age_linear_smooth * age_linear_smooth +
            as.matrix(anyextra_samples$cohort_spline_regional[, paste0("anyextra_cohort_spline_region_mean[", country_region, ", ", 1:length(cohort_spline_values), "]")]) %*% cohort_spline_values

          predictions_cig[j, ] <- mu_cig
          p_cig <- plogis(mu_cig)
          p_smoked_full <- p_cig + plogis(mu_smkextra) * (1 - p_cig)
          predictions_smoked[j, ] <- log(p_smoked_full / (1 - p_smoked_full))
          p_any_full <- p_smoked_full + plogis(mu_anyextra) * (1 - p_smoked_full)
          predictions_any[j, ] <- log(p_any_full / (1 - p_any_full))
        }

        saveRDS(predictions_cig, file = file.path(cig_dir, paste0(current_year, ".rds")), compress = TRUE)
        saveRDS(predictions_smoked, file = file.path(smoked_dir, paste0(current_year, ".rds")), compress = TRUE)
        saveRDS(predictions_any, file = file.path(any_dir, paste0(current_year, ".rds")), compress = TRUE)

        p_cig_matrix <- plogis(predictions_cig)
        p_smoked_matrix <- plogis(predictions_smoked)
        p_any_matrix <- plogis(predictions_any)

        for (k in 1:nrow(current_data)) {
          base_row <- data.frame(
            Year = current_year,
            Age_Midpoint = current_data$Age_Midpoint[k],
            Birth_Cohort = current_data$Birth_Cohort[k],
            Sex = gender,
            Country = country_code,
            Has_Survey_Data = FALSE,
            stringsAsFactors = FALSE
          )

          results_list[[length(results_list) + 1]] <- cbind(base_row, data.frame(
            Def_Type_Code = paste0(current_def, "_cigarettes"),
            Prevalence = mean(p_cig_matrix[k, ]),
            lower_ci = quantile(p_cig_matrix[k, ], 0.025),
            upper_ci = quantile(p_cig_matrix[k, ], 0.975)
          ))

          results_list[[length(results_list) + 1]] <- cbind(base_row, data.frame(
            Def_Type_Code = paste0(current_def, "_any_smoked_tobacco"),
            Prevalence = mean(p_smoked_matrix[k, ]),
            lower_ci = quantile(p_smoked_matrix[k, ], 0.025),
            upper_ci = quantile(p_smoked_matrix[k, ], 0.975)
          ))

          results_list[[length(results_list) + 1]] <- cbind(base_row, data.frame(
            Def_Type_Code = paste0(current_def, "_any_tobacco_product"),
            Prevalence = mean(p_any_matrix[k, ]),
            lower_ci = quantile(p_any_matrix[k, ], 0.025),
            upper_ci = quantile(p_any_matrix[k, ], 0.975)
          ))
        }
      }

      if (length(results_list) > 0) {
        return(do.call(rbind, results_list))
      } else {
        return(NULL)
      }
    }

    # Combine country results with nodata results
    # Handle case where country_results might be NULL (no countries with data for this gender)
    if (is.null(country_results) || nrow(country_results) == 0) {
      # No countries with survey data - use only nodata_results
      if (!is.null(nodata_results) && nrow(nodata_results) > 0) {
        country_results <- nodata_results
        cat(sprintf("  No survey data for %s - using regional priors for all %d countries\n",
                    gender, length(unique(nodata_results$Country))))
      } else {
        cat(sprintf("  WARNING: No predictions generated for %s\n", gender))
        country_results <- NULL
      }
    } else {
      # Has countries with survey data - add nodata results if any
      # Ensure Has_Survey_Data column exists (should already from base_row)
      if (!"Has_Survey_Data" %in% names(country_results)) {
        country_results <- country_results %>%
          mutate(Has_Survey_Data = TRUE)
      }

      if (!is.null(nodata_results) && nrow(nodata_results) > 0) {
        country_results <- bind_rows(country_results, nodata_results)
        cat(sprintf("  Added predictions for %d countries without survey data\n",
                    length(unique(nodata_results$Country))))
      }
    }
  }

  stopCluster(cl)

  gender_results[[gender]] <- country_results
  if (!is.null(country_results)) {
    cat(sprintf("  Generated predictions for %d total countries\n", length(unique(country_results$Country))))
  }

  # Clean up with OPT-1 pattern
  try(nimble::clearCompiled(nimble_model), silent = TRUE)
  rm(compiled_model, compiled_mcmc, mcmc_built, nimble_model)
  gc()
}

# ---- Final Output ----
final_ac_predictions <- do.call(rbind, gender_results)
rownames(final_ac_predictions) <- NULL

write.csv(final_ac_predictions, file = "results/final_ac_predictions_nested.csv", row.names = FALSE)
cat("\nGlobal model predictions saved to: results/final_ac_predictions_nested.csv\n")

# ==================================================================
# SECTION 9: TARGET PREVALENCE CALCULATION
# ==================================================================

cat("\n================================================================\n")
cat("  CALCULATING TARGET PREVALENCES\n")
cat("================================================================\n")

num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = RANDOM_SEED)
clusterExport(cl, "calculate_weighted_prevalence")

target_prevalence_list <- list()

for (gender in genders) {
  cat(sprintf("\nCalculating targets: %s\n", gender))

  countries_global <- list.dirs(file.path("processing", gender), full.names = FALSE, recursive = FALSE)

  global_results <- foreach(
    country = countries_global,
    .packages = c("dplyr", "tidyr"),
    .combine = rbind,
    .errorhandling = 'remove'
  ) %dopar% {

    country_target_list <- list()
    def_types <- list.dirs(file.path("processing", gender, country),
                           full.names = FALSE, recursive = FALSE)

    country_code <- names(country_name_mapping)[country_name_mapping == country]
    if (length(country_code) == 0) country_code <- country

    base_year_weights <- weights_cleaned %>%
      filter(area == country_code, sex == gender, year == BASE_YEAR) %>%
      arrange(age)

    if (nrow(base_year_weights) == 0) return(NULL)

    for (def_type in def_types) {
      base_year_file <- file.path("processing", gender, country, def_type,
                                  paste0(BASE_YEAR, ".rds"))

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
        Model_Type        = "Global",
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

      # Manual Target (e.g., 5%)
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

  target_prevalence_list[[paste0(gender, "_global")]] <- global_results
}

stopCluster(cl)

target_prevalence_df <- do.call(rbind, target_prevalence_list)
rownames(target_prevalence_df) <- NULL

write.csv(target_prevalence_df, file = "results/target_prevalences_global.csv", row.names = FALSE)
cat("Target prevalences saved to: results/target_prevalences_global.csv\n")

# ==================================================================
# SECTION 10: WEIGHTED PREVALENCE TRENDS (GLOBAL MODEL)
# ==================================================================

cat("\n================================================================\n")
cat("  CALCULATING WEIGHTED PREVALENCE TRENDS\n")
cat("================================================================\n")

num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = RANDOM_SEED)
clusterExport(cl, "calculate_weighted_prevalence")

gender_weighted_results <- list()

for (gender in genders) {
  cat(sprintf("Processing weighted calculations: %s\n", gender))

  country_dirs <- list.dirs(file.path("processing", gender),
                            full.names = FALSE, recursive = FALSE)

  weighted_results <- foreach(
    country = country_dirs,
    .packages = c("dplyr", "tidyr"),
    .combine = rbind,
    .errorhandling = 'remove'
  ) %dopar% {

    country_path <- file.path("processing", gender, country)
    def_types <- list.dirs(country_path, full.names = FALSE, recursive = FALSE)

    country_results <- list()

    country_code <- names(country_name_mapping)[country_name_mapping == country]
    if (length(country_code) == 0) country_code <- country

    target_records <- target_prevalence_df %>%
      filter(Sex == gender, Country == country_code)

    for (def_type in def_types) {
      year_files <- list.files(file.path(country_path, def_type), pattern = "\\.rds$")
      years <- as.numeric(sub("\\.rds$", "", year_files))

      def_targets <- target_records %>% filter(Def_Type_Code == def_type)
      if (nrow(def_targets) == 0) next

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

        for (i in 1:nrow(def_targets)) {
          target_record <- def_targets[i, ]

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
            Model_Type           = "Global",
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

  gender_weighted_results[[gender]] <- weighted_results
}

stopCluster(cl)

final_weighted_results_global <- do.call(rbind, gender_weighted_results)
rownames(final_weighted_results_global) <- NULL

write.csv(final_weighted_results_global, file = "results/final_weighted_results_global.csv", row.names = FALSE)
cat("Weighted prevalence trends saved to: results/final_weighted_results_global.csv\n")

cat("\n================================================================\n")
cat("  GLOBAL MODEL MODULE COMPLETE\n")
cat("================================================================\n")
