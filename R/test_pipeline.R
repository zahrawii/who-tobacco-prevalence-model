#########################################################################################
#
#                    WHO TOBACCO CONTROL PREVALENCE PROJECTION MODEL
#                           test_pipeline.R - Module Validation Script
#                                    VERSION 2.3.2
#
#   PURPOSE: Deep validation of the modular pipeline
#     - Runs with minimal settings (4 countries, low MCMC iterations)
#     - Validates each module output step-by-step
#     - Checks data structures, dimensions, and values
#     - Produces detailed diagnostic output
#
#   USAGE:
#     source("R/test_pipeline.R")
#     run_test_pipeline()
#
#########################################################################################

# ============================================================================
# TEST CONFIGURATION (Override production settings)
# ============================================================================

TEST_CONFIG <- list(
  # MCMC Settings - MINIMAL for testing
  NUMBER_OF_CHAINS     = 2,
  NUMBER_OF_BURN       = 50,
  NUMBER_OF_ITERATIONS = 100,
  THINNING_INTERVAL    = 1,

  # Test Countries - 4 from different regions with good data
  # Selected: USA (North America), GBR (Western Europe), CHN (Eastern Asia), BRA (South America)
  TEST_COUNTRIES = c("usa", "gbr", "chn", "bra"),

  # Test with only males to speed up
  TEST_GENDERS = c("males"),  # Set to c("males", "females") for full test

  # Verbosity
  VERBOSE = TRUE
)

# ============================================================================
# VALIDATION HELPER FUNCTIONS
# ============================================================================

validate_df <- function(df, name, required_cols = NULL, expected_rows = NULL) {
  cat(sprintf("\n--- Validating: %s ---\n", name))

  if (!is.data.frame(df) && !is.matrix(df)) {
    cat(sprintf("  ERROR: %s is not a data.frame/matrix (is %s)\n", name, class(df)[1]))
    return(FALSE)
  }

  cat(sprintf("  Dimensions: %d rows x %d cols\n", nrow(df), ncol(df)))

  if (!is.null(expected_rows) && nrow(df) != expected_rows) {
    cat(sprintf("  WARNING: Expected %d rows, got %d\n", expected_rows, nrow(df)))
  }

  if (!is.null(required_cols)) {
    missing <- setdiff(required_cols, names(df))
    if (length(missing) > 0) {
      cat(sprintf("  ERROR: Missing columns: %s\n", paste(missing, collapse = ", ")))
      return(FALSE)
    } else {
      cat(sprintf("  OK: All %d required columns present\n", length(required_cols)))
    }
  }

  # Check for NAs
  na_counts <- colSums(is.na(df))
  cols_with_na <- names(na_counts[na_counts > 0])
  if (length(cols_with_na) > 0) {
    cat(sprintf("  NOTE: %d columns have NA values: %s\n",
                length(cols_with_na),
                paste(head(cols_with_na, 5), collapse = ", ")))
  }

  # Print sample
  cat("  First few rows:\n")
  print(head(df, 3))

  return(TRUE)
}

validate_list <- function(lst, name, required_elements = NULL) {
  cat(sprintf("\n--- Validating: %s ---\n", name))

  if (!is.list(lst)) {
    cat(sprintf("  ERROR: %s is not a list (is %s)\n", name, class(lst)[1]))
    return(FALSE)
  }

  cat(sprintf("  Length: %d elements\n", length(lst)))
  cat(sprintf("  Elements: %s\n", paste(names(lst), collapse = ", ")))

  if (!is.null(required_elements)) {
    missing <- setdiff(required_elements, names(lst))
    if (length(missing) > 0) {
      cat(sprintf("  ERROR: Missing elements: %s\n", paste(missing, collapse = ", ")))
      return(FALSE)
    }
  }

  # Print structure summary
  for (elem in names(lst)) {
    val <- lst[[elem]]
    if (is.numeric(val) && length(val) == 1) {
      cat(sprintf("    %s: %s\n", elem, val))
    } else if (is.vector(val)) {
      cat(sprintf("    %s: vector[%d]\n", elem, length(val)))
    } else if (is.matrix(val)) {
      cat(sprintf("    %s: matrix[%d x %d]\n", elem, nrow(val), ncol(val)))
    } else {
      cat(sprintf("    %s: %s\n", elem, class(val)[1]))
    }
  }

  return(TRUE)
}

validate_nimble_samples <- function(samples, name) {
  cat(sprintf("\n--- Validating MCMC Samples: %s ---\n", name))

  if (is.list(samples)) {
    # Multiple chains
    n_chains <- length(samples)
    cat(sprintf("  Number of chains: %d\n", n_chains))

    for (i in 1:n_chains) {
      chain <- samples[[i]]
      if (is.matrix(chain)) {
        cat(sprintf("  Chain %d: %d samples x %d parameters\n", i, nrow(chain), ncol(chain)))
        if (i == 1) {
          cat(sprintf("  Parameters: %s...\n", paste(head(colnames(chain), 5), collapse = ", ")))
        }
      }
    }
  } else if (is.matrix(samples)) {
    cat(sprintf("  Single matrix: %d samples x %d parameters\n", nrow(samples), ncol(samples)))
    cat(sprintf("  Parameters: %s...\n", paste(head(colnames(samples), 10), collapse = ", ")))
  }

  return(TRUE)
}

# ============================================================================
# MODULE 00: CONFIG TEST
# ============================================================================

test_module_00_config <- function() {
  cat("\n")
  cat("################################################################\n")
  cat("# MODULE 00: CONFIGURATION\n")
  cat("################################################################\n")

  # Check DLL limit
  dll_limit <- Sys.getenv("R_MAX_NUM_DLLS")
  cat(sprintf("  R_MAX_NUM_DLLS: %s\n", dll_limit))

  # Override with test config
  NUMBER_OF_CHAINS     <<- TEST_CONFIG$NUMBER_OF_CHAINS
  NUMBER_OF_BURN       <<- TEST_CONFIG$NUMBER_OF_BURN
  NUMBER_OF_ITERATIONS <<- TEST_CONFIG$NUMBER_OF_ITERATIONS
  THINNING_INTERVAL    <<- TEST_CONFIG$THINNING_INTERVAL

  cat(sprintf("\n  TEST MCMC Settings Applied:\n"))
  cat(sprintf("    Chains:     %d\n", NUMBER_OF_CHAINS))
  cat(sprintf("    Burn-in:    %d\n", NUMBER_OF_BURN))
  cat(sprintf("    Iterations: %d\n", NUMBER_OF_ITERATIONS))
  cat(sprintf("    Thinning:   %d\n", THINNING_INTERVAL))

  # Check directories
  test_dirs <- c("processing", "country_priors", "results", "checkpoints")
  for (d in test_dirs) {
    if (dir.exists(d)) {
      cat(sprintf("  Directory OK: %s\n", d))
    } else {
      cat(sprintf("  Creating directory: %s\n", d))
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
  }

  return(TRUE)
}

# ============================================================================
# MODULE 01: DATA PREP TEST
# ============================================================================

test_module_01_data_prep <- function() {
  cat("\n")
  cat("################################################################\n")
  cat("# MODULE 01: DATA PREPARATION\n")
  cat("################################################################\n")

  if (!exists("clean_data")) {
    cat("  ERROR: clean_data does not exist!\n")
    return(FALSE)
  }

  # Validate structure
  required_cols <- c(
    "country", "wb_country_abv", "year", "sex", "prevalence",
    "start_age", "end_age", "Age_Midpoint", "def_code", "type_code",
    "def_code_binary", "def_type_code"
  )

  if (!validate_df(clean_data, "clean_data", required_cols)) {
    return(FALSE)
  }

  # Country summary
  country_counts <- clean_data %>%
    group_by(wb_country_abv) %>%
    summarise(n_obs = n(), .groups = "drop") %>%
    arrange(desc(n_obs))

  cat(sprintf("\n  Total countries in data: %d\n", nrow(country_counts)))
  cat(sprintf("  Total observations: %d\n", nrow(clean_data)))

  # Check test countries exist
  test_countries <- TEST_CONFIG$TEST_COUNTRIES
  available <- test_countries[test_countries %in% clean_data$wb_country_abv]
  missing <- setdiff(test_countries, available)

  cat(sprintf("\n  Test countries available: %s\n", paste(available, collapse = ", ")))
  if (length(missing) > 0) {
    cat(sprintf("  WARNING: Test countries NOT in data: %s\n", paste(missing, collapse = ", ")))
  }

  # Sex distribution
  cat("\n  Sex distribution:\n")
  print(table(clean_data$sex))

  # Prevalence range check
  cat(sprintf("\n  Prevalence range: [%.3f, %.3f]\n",
              min(clean_data$prevalence, na.rm = TRUE),
              max(clean_data$prevalence, na.rm = TRUE)))

  return(TRUE)
}

# ============================================================================
# MODULE 02: REGIONAL MAPPING TEST
# ============================================================================

test_module_02_regional_mapping <- function() {
  cat("\n")
  cat("################################################################\n")
  cat("# MODULE 02: REGIONAL MAPPING\n")
  cat("################################################################\n")

  if (!exists("country_region_manual")) {
    cat("  ERROR: country_region_manual does not exist!\n")
    return(FALSE)
  }

  if (!validate_df(country_region_manual, "country_region_manual",
                   c("wb_country_abv", "region_consolidated"))) {
    return(FALSE)
  }

  # Region summary
  region_counts <- country_region_manual %>%
    group_by(region_consolidated) %>%
    summarise(n_countries = n(), .groups = "drop") %>%
    arrange(desc(n_countries))

  cat("\n  Regions and country counts:\n")
  print(region_counts)

  # Check test countries mapping
  test_mapping <- country_region_manual %>%
    filter(wb_country_abv %in% TEST_CONFIG$TEST_COUNTRIES)

  cat("\n  Test country regions:\n")
  print(test_mapping)

  # Check clean_data has region
  if (!"region" %in% names(clean_data)) {
    cat("  ERROR: clean_data missing 'region' column!\n")
    return(FALSE)
  }

  return(TRUE)
}

# ============================================================================
# MODULE 03: SPLINES TEST
# ============================================================================

test_module_03_splines <- function() {
  cat("\n")
  cat("################################################################\n")
  cat("# MODULE 03: SPLINE BASIS\n")
  cat("################################################################\n")

  # Check spline columns exist
  age_spline_cols <- grep("^age_spline_", names(clean_data), value = TRUE)
  cohort_spline_cols <- grep("^cohort_spline_", names(clean_data), value = TRUE)
  interaction_cols <- grep("^age_cohort_", names(clean_data), value = TRUE)

  cat(sprintf("\n  Age spline columns: %d\n", length(age_spline_cols)))
  cat(sprintf("  Cohort spline columns: %d\n", length(cohort_spline_cols)))
  cat(sprintf("  Interaction columns: %d (expected: %d)\n",
              length(interaction_cols),
              length(age_spline_cols) * length(cohort_spline_cols)))

  # Check required variables
  required <- c("Birth_Cohort", "Cohort_Centered", "spline_weight_var",
                "linear_weight_var", "age_linear_smooth", "Type_Cig",
                "Type_Smoked", "Type_Any")

  missing <- setdiff(required, names(clean_data))
  if (length(missing) > 0) {
    cat(sprintf("  ERROR: Missing columns: %s\n", paste(missing, collapse = ", ")))
    return(FALSE)
  }

  # Check centering constants
  if (!exists("AGE_LINEAR_CENTER_CONSTANT")) {
    cat("  ERROR: AGE_LINEAR_CENTER_CONSTANT not defined!\n")
    return(FALSE)
  }
  if (!exists("COHORT_CENTER_CONSTANT")) {
    cat("  ERROR: COHORT_CENTER_CONSTANT not defined!\n")
    return(FALSE)
  }

  cat(sprintf("\n  AGE_LINEAR_CENTER_CONSTANT: %.3f\n", AGE_LINEAR_CENTER_CONSTANT))
  cat(sprintf("  COHORT_CENTER_CONSTANT: %.3f\n", COHORT_CENTER_CONSTANT))

  # Check prevalence is on logit scale
  prev_range <- range(clean_data$prevalence, na.rm = TRUE)
  cat(sprintf("\n  Prevalence range (should be logit): [%.3f, %.3f]\n",
              prev_range[1], prev_range[2]))

  if (prev_range[2] > 0 && prev_range[1] > -10) {
    cat("  OK: Prevalence appears to be on logit scale\n")
  } else {
    cat("  WARNING: Prevalence may not be on logit scale\n")
  }

  # Check spline attribute objects
  if (!exists("age_spline_knots_attr")) {
    cat("  ERROR: age_spline_knots_attr not defined!\n")
    return(FALSE)
  }

  cat(sprintf("\n  Age spline knots: %s\n", paste(age_spline_knots_attr, collapse = ", ")))
  cat(sprintf("  Cohort spline knots: %s\n", paste(cohort_spline_knots_attr, collapse = ", ")))

  return(TRUE)
}

# ============================================================================
# MODULE 04: UTILS TEST
# ============================================================================

test_module_04_utils <- function() {
  cat("\n")
  cat("################################################################\n")
  cat("# MODULE 04: UTILITY FUNCTIONS\n")
  cat("################################################################\n")

  # Check weights data
  if (!exists("weights_cleaned")) {
    cat("  ERROR: weights_cleaned not loaded!\n")
    return(FALSE)
  }

  cat(sprintf("\n  weights_cleaned:\n"))
  cat(sprintf("    Rows: %d\n", nrow(weights_cleaned)))
  cat(sprintf("    Years: %d - %d\n", min(weights_cleaned$year), max(weights_cleaned$year)))
  cat(sprintf("    Countries: %d\n", length(unique(weights_cleaned$area))))

  # Check test countries have weights
  test_countries <- tolower(TEST_CONFIG$TEST_COUNTRIES)
  weights_countries <- unique(weights_cleaned$area)

  available_weights <- test_countries[test_countries %in% weights_countries]
  missing_weights <- setdiff(test_countries, weights_countries)

  cat(sprintf("    Test countries with weights: %s\n", paste(available_weights, collapse = ", ")))
  if (length(missing_weights) > 0) {
    cat(sprintf("    WARNING: No weights for: %s\n", paste(missing_weights, collapse = ", ")))
  }

  # Check if calculate_weighted_prevalence exists
  if (!exists("calculate_weighted_prevalence")) {
    cat("  ERROR: calculate_weighted_prevalence function not defined!\n")
    return(FALSE)
  }

  # Test the function
  test_matrix <- matrix(rnorm(100, mean = -2, sd = 0.5), nrow = 10, ncol = 10)
  test_weights <- rep(1/10, 10)

  result <- tryCatch({
    calculate_weighted_prevalence(test_matrix, test_weights)
  }, error = function(e) {
    cat(sprintf("  ERROR in calculate_weighted_prevalence: %s\n", e$message))
    return(NULL)
  })

  if (!is.null(result)) {
    cat(sprintf("\n  calculate_weighted_prevalence test: OK (returned %d values)\n", length(result)))
  }

  # Check other utility functions
  funcs_to_check <- c("check_ordering_per_draw", "prec_to_sd",
                      "create_prediction_splines", "create_transition_weights")
  for (f in funcs_to_check) {
    if (exists(f)) {
      cat(sprintf("  %s: defined\n", f))
    } else {
      cat(sprintf("  WARNING: %s not defined\n", f))
    }
  }

  return(TRUE)
}

# ============================================================================
# MODULE 05: NIMBLE MODELS TEST
# ============================================================================

test_module_05_models <- function() {
  cat("\n")
  cat("################################################################\n")
  cat("# MODULE 05: NIMBLE MODEL DEFINITIONS\n")
  cat("################################################################\n")

  # Check model objects exist
  models_to_check <- c(
    "regional_hierarchical_global_ac_model_nimble",
    "regional_country_specific_ac_model_nimble"
  )

  for (model in models_to_check) {
    if (exists(model)) {
      obj <- get(model)
      cat(sprintf("  %s: %s\n", model, class(obj)[1]))
    } else {
      cat(sprintf("  ERROR: %s not defined!\n", model))
      return(FALSE)
    }
  }

  return(TRUE)
}

# ============================================================================
# SUBSET DATA FOR TEST
# ============================================================================

subset_data_for_test <- function() {
  cat("\n")
  cat("################################################################\n")
  cat("# SUBSETTING DATA FOR TEST\n")
  cat("################################################################\n")

  # Store original
  clean_data_full <<- clean_data

  # Subset to test countries
  test_countries <- TEST_CONFIG$TEST_COUNTRIES
  test_genders <- TEST_CONFIG$TEST_GENDERS

  clean_data <<- clean_data %>%
    filter(wb_country_abv %in% test_countries,
           sex %in% test_genders)

  cat(sprintf("\n  Original data: %d rows\n", nrow(clean_data_full)))
  cat(sprintf("  Test subset: %d rows\n", nrow(clean_data)))
  cat(sprintf("  Countries: %s\n", paste(unique(clean_data$wb_country_abv), collapse = ", ")))
  cat(sprintf("  Genders: %s\n", paste(unique(clean_data$sex), collapse = ", ")))

  # Update regional mapping
  country_region_mapping <<- clean_data %>%
    select(wb_country_abv) %>%
    distinct() %>%
    left_join(country_region_manual, by = "wb_country_abv")

  cat(sprintf("  Regions in subset: %s\n",
              paste(unique(country_region_mapping$region_consolidated), collapse = ", ")))

  return(TRUE)
}

# ============================================================================
# MODULE 06: GLOBAL MODEL TEST (SIMPLIFIED)
# ============================================================================

test_module_06_global_model <- function() {
  cat("\n")
  cat("################################################################\n")
  cat("# MODULE 06: GLOBAL MODEL (TEST RUN)\n")
  cat("################################################################\n")

  cat("\n  Running global model with test settings...\n")
  cat(sprintf("  This will fit %d gender(s) with %d chains x %d iterations\n",
              length(TEST_CONFIG$TEST_GENDERS),
              TEST_CONFIG$NUMBER_OF_CHAINS,
              TEST_CONFIG$NUMBER_OF_ITERATIONS))

  # CRITICAL: Override the genders variable to match TEST_CONFIG
  # The global model uses a 'genders' variable in its loop
  genders <<- TEST_CONFIG$TEST_GENDERS

  start_time <- Sys.time()

  tryCatch({
    source("R/06_run_global_model.R")

    end_time <- Sys.time()
    cat(sprintf("\n  Global model completed in: %.2f minutes\n",
                as.numeric(difftime(end_time, start_time, units = "mins"))))

    # Validate outputs
    if (exists("final_ac_predictions")) {
      validate_df(final_ac_predictions, "final_ac_predictions",
                  c("Year", "Age_Midpoint", "Sex", "Country", "Prevalence"))

      # Additional validation checks
      cat(sprintf("\n  Predictions summary:\n"))
      cat(sprintf("    Total rows: %d\n", nrow(final_ac_predictions)))
      cat(sprintf("    Countries: %s\n", paste(unique(final_ac_predictions$Country), collapse = ", ")))
      cat(sprintf("    Years: %d - %d\n",
                  min(final_ac_predictions$Year),
                  max(final_ac_predictions$Year)))
      cat(sprintf("    Prevalence range: [%.4f, %.4f]\n",
                  min(final_ac_predictions$Prevalence),
                  max(final_ac_predictions$Prevalence)))

      # Check for Has_Survey_Data column (no-data countries)
      if ("Has_Survey_Data" %in% names(final_ac_predictions)) {
        n_data <- sum(final_ac_predictions$Has_Survey_Data)
        n_nodata <- sum(!final_ac_predictions$Has_Survey_Data)
        cat(sprintf("    With survey data: %d rows\n", n_data))
        cat(sprintf("    Without survey data: %d rows\n", n_nodata))
      }
    } else {
      cat("  ERROR: final_ac_predictions not created!\n")
      return(FALSE)
    }

    if (exists("gender_results")) {
      cat(sprintf("\n  gender_results has %d elements\n", length(gender_results)))
    }

    return(TRUE)

  }, error = function(e) {
    cat(sprintf("\n  ERROR in global model: %s\n", e$message))
    cat(sprintf("  Error location: %s\n", deparse(conditionCall(e))))
    return(FALSE)
  })
}

# ============================================================================
# MODULE 07: COUNTRY MODEL TEST (SIMPLIFIED)
# ============================================================================

test_module_07_country_model <- function() {
  cat("\n")
  cat("################################################################\n")
  cat("# MODULE 07: COUNTRY MODEL (TEST RUN)\n")
  cat("################################################################\n")

  # Check prerequisites
  if (!dir.exists("country_priors")) {
    cat("  ERROR: country_priors directory missing!\n")
    return(FALSE)
  }

  # Check for country priors files
  priors_dir <- file.path("country_priors", TEST_CONFIG$TEST_GENDERS[1])
  if (dir.exists(priors_dir)) {
    prior_files <- list.files(priors_dir, pattern = "\\.csv$")
    cat(sprintf("\n  Found %d country prior files in %s\n", length(prior_files), priors_dir))
    if (length(prior_files) > 0) {
      cat(sprintf("  Sample files: %s\n", paste(head(prior_files, 3), collapse = ", ")))
    }
  }

  cat("\n  Running country-specific model with test settings...\n")

  # CRITICAL: Override the genders variable to match TEST_CONFIG
  genders <<- TEST_CONFIG$TEST_GENDERS

  start_time <- Sys.time()

  tryCatch({
    source("R/07_run_country_model.R")

    end_time <- Sys.time()
    cat(sprintf("\n  Country model completed in: %.2f minutes\n",
                as.numeric(difftime(end_time, start_time, units = "mins"))))

    # Validate outputs
    if (exists("final_predictions_country_specific_ac")) {
      validate_df(final_predictions_country_specific_ac,
                  "final_predictions_country_specific_ac",
                  c("Year", "Age_Midpoint", "Sex", "Country", "Prevalence"))

      # Additional validation
      cat(sprintf("\n  Country-specific predictions summary:\n"))
      cat(sprintf("    Total rows: %d\n", nrow(final_predictions_country_specific_ac)))
      cat(sprintf("    Countries fitted: %s\n",
                  paste(unique(final_predictions_country_specific_ac$Country), collapse = ", ")))
      cat(sprintf("    Def types: %s\n",
                  paste(unique(final_predictions_country_specific_ac$Def_Type_Code), collapse = ", ")))
    } else {
      cat("  WARNING: final_predictions_country_specific_ac not created\n")
      cat("  (This may be OK if no countries qualified for country-specific fitting)\n")
    }

    return(TRUE)

  }, error = function(e) {
    cat(sprintf("\n  ERROR in country model: %s\n", e$message))
    return(FALSE)
  })
}

# ============================================================================
# MAIN TEST RUNNER
# ============================================================================

run_test_pipeline <- function(modules = c("config", "data", "region", "splines",
                                           "utils", "models", "subset",
                                           "global", "country")) {
  cat("\n")
  cat("################################################################\n")
  cat("#\n")
  cat("#  WHO TOBACCO MODEL - MODULAR PIPELINE TEST\n")
  cat("#\n")
  cat("#  Test Configuration:\n")
  cat(sprintf("#    Countries: %s\n", paste(TEST_CONFIG$TEST_COUNTRIES, collapse = ", ")))
  cat(sprintf("#    Genders: %s\n", paste(TEST_CONFIG$TEST_GENDERS, collapse = ", ")))
  cat(sprintf("#    MCMC: %d chains x %d iterations\n",
              TEST_CONFIG$NUMBER_OF_CHAINS, TEST_CONFIG$NUMBER_OF_ITERATIONS))
  cat("#\n")
  cat("################################################################\n")

  results <- list()

  # Phase 1: Source preparation modules
  if ("config" %in% modules || "data" %in% modules) {
    cat("\n=== PHASE 1: Sourcing Preparation Modules ===\n")

    tryCatch({
      source("R/00_config.R")
      source("R/01_data_prep.R")
      source("R/02_regional_mapping.R")
      source("R/03_splines.R")
      source("R/04_utils.R")
      source("R/05_models_nimble.R")
    }, error = function(e) {
      cat(sprintf("ERROR sourcing modules: %s\n", e$message))
      return(results)
    })
  }

  # Test each module
  if ("config" %in% modules) {
    results$config <- test_module_00_config()
  }

  if ("data" %in% modules) {
    results$data <- test_module_01_data_prep()
  }

  if ("region" %in% modules) {
    results$region <- test_module_02_regional_mapping()
  }

  if ("splines" %in% modules) {
    results$splines <- test_module_03_splines()
  }

  if ("utils" %in% modules) {
    results$utils <- test_module_04_utils()
  }

  if ("models" %in% modules) {
    results$models <- test_module_05_models()
  }

  # Subset data for testing
  if ("subset" %in% modules) {
    results$subset <- subset_data_for_test()
  }

  # Test model fitting (optional, slow)
  if ("global" %in% modules) {
    results$global <- test_module_06_global_model()
  }

  if ("country" %in% modules) {
    results$country <- test_module_07_country_model()
  }

  # Summary
  cat("\n")
  cat("################################################################\n")
  cat("#  TEST SUMMARY\n")
  cat("################################################################\n")

  for (name in names(results)) {
    status <- if (isTRUE(results[[name]])) "PASS" else "FAIL"
    cat(sprintf("  %-12s: %s\n", name, status))
  }

  n_pass <- sum(sapply(results, isTRUE))
  n_total <- length(results)

  cat(sprintf("\n  Total: %d/%d passed\n", n_pass, n_total))

  if (n_pass == n_total) {
    cat("\n  ALL TESTS PASSED!\n")
  } else {
    cat("\n  SOME TESTS FAILED - Review output above\n")
  }

  invisible(results)
}

# ============================================================================
# QUICK VALIDATION (No model fitting)
# ============================================================================

run_quick_validation <- function() {
  run_test_pipeline(modules = c("config", "data", "region", "splines", "utils", "models"))
}

# ============================================================================
# FULL VALIDATION (With model fitting)
# ============================================================================

run_full_validation <- function() {
  run_test_pipeline(modules = c("config", "data", "region", "splines",
                                "utils", "models", "subset", "global", "country"))
}

# ============================================================================
# Print usage
# ============================================================================

cat("\n")
cat("################################################################\n")
cat("#  TEST PIPELINE LOADED\n")
cat("################################################################\n")
cat("\nAvailable commands:\n")
cat("  run_quick_validation()   - Test modules 00-05 (no model fitting)\n")
cat("  run_full_validation()    - Full test including MCMC (slow)\n")
cat("  run_test_pipeline()      - Run specific modules\n")
cat("\nTest Configuration:\n")
cat(sprintf("  Countries: %s\n", paste(TEST_CONFIG$TEST_COUNTRIES, collapse = ", ")))
cat(sprintf("  MCMC: %d chains x %d iterations\n",
            TEST_CONFIG$NUMBER_OF_CHAINS, TEST_CONFIG$NUMBER_OF_ITERATIONS))
cat("\n")
