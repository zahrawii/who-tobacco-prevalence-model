#########################################################################################
#
#                    WHO TOBACCO CONTROL PREVALENCE PROJECTION MODEL
#                         03_splines.R - Spline Basis Construction
#
#   Contains: Age-cohort transformation, spline basis creation, interaction tensors
#   Requires: 00_config.R, 01_data_prep.R, 02_regional_mapping.R
#   Outputs: Updates clean_data with spline columns, spline attribute objects
#
#########################################################################################

# ---- 3.1 Calculate Birth Cohorts ----

clean_data <- clean_data %>%
  mutate(
    Birth_Cohort    = year - Age_Midpoint,
    Cohort_Centered = Birth_Cohort - median(Birth_Cohort, na.rm = TRUE)
  )

# ---- 3.2 Create Smooth Sigmoid Transition Weights ----

clean_data <- clean_data %>%
  mutate(
    sigmoid_input       = (Age_Midpoint - TRANSITION_START) / TRANSITION_WIDTH,
    spline_weight       = 1 / (1 + exp(sigmoid_input)),
    linear_weight       = 1 - spline_weight,
    Age_For_Spline      = Age_Midpoint,
    Age_Linear          = pmax(0, Age_Midpoint - TRANSITION_START),
    Age_Linear_Centered = Age_Linear - mean(Age_Linear, na.rm = TRUE)
  )

# ---- 3.2.1 Store Global Centering Constants ----

AGE_LINEAR_CENTER_CONSTANT <- mean(clean_data$Age_Linear, na.rm = TRUE)
COHORT_CENTER_CONSTANT     <- median(clean_data$Birth_Cohort, na.rm = TRUE)

cat("\nGlobal centering constants stored:\n")
cat("  AGE_LINEAR_CENTER_CONSTANT:", AGE_LINEAR_CENTER_CONSTANT, "\n")
cat("  COHORT_CENTER_CONSTANT:    ", COHORT_CENTER_CONSTANT, "\n\n")

# ---- 3.3 Define Spline Knots ----
# ============================================================================
# LITERATURE-INFORMED KNOT SELECTION
# Based on deep review of smoking epidemiology literature
# ============================================================================
#
# AGE KNOT RATIONALE (from GBD 2015, NHANES APC analysis, CISNET):
# ----------------------------------------------------------------
# The age-prevalence curve has a characteristic shape:
#   - Rapid rise from ~15 to ~25 (initiation phase: 82.6% initiate by age 25)
#   - Peak at 25-34 years (CDC data; GBD shows males peak 25-35 globally)
#   - Plateau from ~30-45 (stable prevalence, modest cessation)
#   - Decline from ~45+ (cessation accelerates; U-shaped cessation by age)
#   - Accelerated decline ~65+ (cumulative cessation + differential mortality)
#
# Key inflection points identified:
#   25 = End of initiation / peak onset (after this, odds of smoking decline)
#   45 = Decline onset (cessation begins to dominate; lowest cessation rates 45-64)
#   65 = Acceleration (mortality selection kicks in; cessation rates rise again)
#
# COHORT KNOT RATIONALE (from Holford/CISNET, 50-Year Surgeon General Report):
# ----------------------------------------------------------------
# The cohort-prevalence curve shows:
#   - Peak for males: 1920-1930 birth cohorts
#   - Peak for females: 1935-1945 birth cohorts
#   - 1950-1970 cohorts: HIGHEST overall smoking likelihood (both sexes)
#   - 1964 SGR impact: 1940 cohort was age 25 (past peak initiation)
#   - Post-1970 cohorts: Declining due to tobacco control
#
# Key inflection points:
#   1945 = Post-war peak (esp. females); baby boom generation
#   1965 = SGR generation (first cohort born after 1964 report)
#   1985 = Modern tobacco control era (taxes, bans, media campaigns)
#
# References:
#   - GBD 2015 Lancet: "Male prevalence peaks between ages 25 and 35"
#   - NHANES 2025 Sci Rep: "Odds of smoking decrease after age ~27"
#   - Holford 2014 AJPM: CISNET smoking history methodology
#   - 50-Year Surgeon General Report (2014): Birth cohort patterns
# ============================================================================

age_knots           <- c(25, 45, 65)
age_boundary_knots  <- c(15, MAX_AGE_SPLINE)

cohort_range        <- range(clean_data$Birth_Cohort, na.rm = TRUE)
cohort_knots        <- c(1945, 1965, 1985)
# Ensure knots fall within data range
cohort_knots        <- cohort_knots[cohort_knots > cohort_range[1] &
                                      cohort_knots < cohort_range[2]]
cohort_boundary_knots <- cohort_range

cat("\nSpline Configuration:\n")
cat("  Age knots:        ", paste(age_knots, collapse = ", "), "\n")
cat("  Age boundaries:   ", paste(age_boundary_knots, collapse = ", "), "\n")
cat("  Cohort knots:     ", paste(cohort_knots, collapse = ", "), "\n")
cat("  Cohort boundaries:", paste(cohort_boundary_knots, collapse = ", "), "\n")
cat("  Note: Ages 80-100 use linear extrapolation\n")
cat("  Rationale: Knots placed at epidemiologically meaningful transitions\n")
cat("             (peak prevalence, cessation onset, mortality selection)\n\n")

# ---- 3.4 Create Spline Bases ----

age_spline_basis <- ns(
  clean_data$Age_For_Spline,
  knots = age_knots,
  Boundary.knots = age_boundary_knots
)
age_spline_basis_df <- as.data.frame(age_spline_basis)
colnames(age_spline_basis_df) <- paste0("age_spline_", 1:ncol(age_spline_basis_df))

cohort_spline_basis <- ns(
  clean_data$Birth_Cohort,
  knots = cohort_knots,
  Boundary.knots = cohort_boundary_knots
)
cohort_spline_basis_df <- as.data.frame(cohort_spline_basis)
colnames(cohort_spline_basis_df) <- paste0("cohort_spline_", 1:ncol(cohort_spline_basis_df))

# Store knot attributes for prediction
age_spline_knots_attr    <- attr(age_spline_basis, "knots")
age_spline_boundary_attr <- attr(age_spline_basis, "Boundary.knots")
cohort_spline_knots_attr <- attr(cohort_spline_basis, "knots")
cohort_spline_boundary_attr <- attr(cohort_spline_basis, "Boundary.knots")

clean_data <- cbind(clean_data, age_spline_basis_df, cohort_spline_basis_df)

# ---- 3.5 Create Weight Variables ----

clean_data$age_linear_smooth   <- clean_data$Age_Linear_Centered
clean_data$spline_weight_var   <- clean_data$spline_weight
clean_data$linear_weight_var   <- clean_data$linear_weight

# ---- 3.6 Build Age-Cohort Interaction Tensor ----

n_age_splines    <- ncol(age_spline_basis_df)
n_cohort_splines <- ncol(cohort_spline_basis_df)
n_interactions   <- n_age_splines * n_cohort_splines

age_cohort_interaction_matrix <- matrix(0, nrow = nrow(clean_data), ncol = n_interactions)
for (i in 1:nrow(clean_data)) {
  age_values    <- as.numeric(clean_data[i, grep("^age_spline_", names(clean_data))])
  cohort_values <- as.numeric(clean_data[i, grep("^cohort_spline_", names(clean_data))])
  age_cohort_interaction_matrix[i, ] <- as.vector(outer(age_values, cohort_values))
}

age_cohort_interaction_df <- as.data.frame(age_cohort_interaction_matrix)
colnames(age_cohort_interaction_df) <- paste0("age_cohort_", 1:ncol(age_cohort_interaction_df))
clean_data <- cbind(clean_data, age_cohort_interaction_df)

# ---- 3.7 Transform Prevalence to Logit Scale ----

clean_data$prevalence <- clean_data$prevalence / 100
clean_data$prevalence <- ifelse(clean_data$prevalence == 0, 0.001, clean_data$prevalence)
clean_data$prevalence <- ifelse(clean_data$prevalence == 1, 0.999, clean_data$prevalence)
clean_data$prevalence <- log(clean_data$prevalence / (1 - clean_data$prevalence))

# ---- 3.8 Create Product Type Indicators ----

clean_data <- clean_data %>%
  mutate(
    Type_Cig    = as.integer(type_cigarettes == 1),
    Type_Smoked = as.integer(type_smoked_tobacco == 1),
    Type_Any    = as.integer(type_any_tobacco == 1)
  )

cat("Spline basis construction complete.\n")
cat(sprintf("  Age spline columns: %d\n", n_age_splines))
cat(sprintf("  Cohort spline columns: %d\n", n_cohort_splines))
cat(sprintf("  Interaction columns: %d\n", n_interactions))
