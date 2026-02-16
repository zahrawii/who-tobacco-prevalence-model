#########################################################################################
#
#                    WHO TOBACCO CONTROL PREVALENCE PROJECTION MODEL
#                      04_utils.R - Utility Functions and Weights
#
#   Contains: Population weight loading, utility functions for model operations
#   Requires: 00_config.R
#   Outputs: weights_data, weights_cleaned, utility functions
#
#########################################################################################

# ---- 4.1 Load Population Weight Data ----

weights_data <- read_csv("data/weights_15_2022.csv", show_col_types = FALSE)

weights_cleaned <- weights_data %>%
  rename(`100` = `100+`) %>%
  pivot_longer(cols = `15`:`100`, names_to = "age", values_to = "weight") %>%
  mutate(
    area = tolower(area),
    sex  = tolower(sex),
    year = as.numeric(year),
    age  = as.numeric(age),
    sex  = case_when(
      sex == "male"   ~ "males",
      sex == "female" ~ "females",
      TRUE ~ sex
    )
  ) %>%
  arrange(area, year, sex, age)

cat("\nPopulation weights loaded.\n")
cat(sprintf("  Years: %d-%d\n", min(weights_cleaned$year), max(weights_cleaned$year)))
cat(sprintf("  Countries: %d\n", length(unique(weights_cleaned$area))))

#########################################################################################
#                              UTILITY FUNCTIONS                                         #
#########################################################################################

# ---- 4.2 Check Ordering Constraints ----
# Verifies stick-breaking constraint: P(cig) <= P(smoked) <= P(any)

check_ordering_per_draw <- function(logit_cig, logit_smk, logit_any, tolerance = 1e-6) {
  pc <- plogis(logit_cig)
  ps <- plogis(logit_smk)
  pa <- plogis(logit_any)

  violation_1 <- pc > ps + tolerance
  violation_2 <- ps > pa + tolerance

  n_viol_1 <- sum(violation_1)
  n_viol_2 <- sum(violation_2)

  return(n_viol_1 > 0 || n_viol_2 > 0)
}

# ---- 4.3 Calculate Population-Weighted Prevalence ----
# Aggregates age-specific prevalence using population weights

calculate_weighted_prevalence <- function(iterations_matrix_logit, weights_vector) {
  p_matrix <- plogis(iterations_matrix_logit)
  weighted_probs <- numeric(ncol(p_matrix))

  for (s in 1:ncol(p_matrix)) {
    weighted_probs[s] <- sum(p_matrix[, s] * weights_vector) / sum(weights_vector)
  }

  return(weighted_probs)
}

# ---- 4.4 Precision to SD Conversion Helper ----
# Converts NIMBLE precision parameters to standard deviation

prec_to_sd <- function(precision) {
  1 / sqrt(pmax(precision, 1e-10))
}

# ---- 4.5 Create Prediction Spline Bases ----
# Helper function to create spline bases for new prediction data

create_prediction_splines <- function(ages, cohorts,
                                       age_knots_attr, age_boundary_attr,
                                       cohort_knots_attr, cohort_boundary_attr) {

  # Age splines
  age_basis <- ns(ages,
                  knots = age_knots_attr,
                  Boundary.knots = age_boundary_attr)
  age_basis_df <- as.data.frame(age_basis)
  colnames(age_basis_df) <- paste0("age_spline_", 1:ncol(age_basis_df))

  # Cohort splines
  cohort_basis <- ns(cohorts,
                     knots = cohort_knots_attr,
                     Boundary.knots = cohort_boundary_attr)
  cohort_basis_df <- as.data.frame(cohort_basis)
  colnames(cohort_basis_df) <- paste0("cohort_spline_", 1:ncol(cohort_basis_df))

  # Interaction tensor
  n_age <- ncol(age_basis_df)
  n_cohort <- ncol(cohort_basis_df)
  n_obs <- length(ages)
  n_interactions <- n_age * n_cohort

  interaction_matrix <- matrix(0, nrow = n_obs, ncol = n_interactions)
  for (i in 1:n_obs) {
    age_values <- as.numeric(age_basis_df[i, ])
    cohort_values <- as.numeric(cohort_basis_df[i, ])
    interaction_matrix[i, ] <- as.vector(outer(age_values, cohort_values))
  }
  interaction_df <- as.data.frame(interaction_matrix)
  colnames(interaction_df) <- paste0("age_cohort_", 1:ncol(interaction_df))

  return(list(
    age_basis = age_basis_df,
    cohort_basis = cohort_basis_df,
    interaction = interaction_df
  ))
}

# ---- 4.6 Sigmoid Transition Weights ----
# Creates smooth transition weights for age extrapolation

create_transition_weights <- function(ages, transition_start = TRANSITION_START,
                                        transition_width = TRANSITION_WIDTH) {
  sigmoid_input <- (ages - transition_start) / transition_width
  spline_weight <- 1 / (1 + exp(sigmoid_input))
  linear_weight <- 1 - spline_weight

  return(list(
    spline_weight = spline_weight,
    linear_weight = linear_weight
  ))
}

cat("Utility functions loaded.\n")
