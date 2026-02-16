#########################################################################################
#                                                                                       
#                    WHO TOBACCO CONTROL PREVALENCE PROJECTION MODEL                    
#                                                                                       
#   A Bayesian Age-Cohort hierarchical model for projecting tobacco smoking            
#   prevalence and evaluating WHO reduction targets                                     
#                                                                                       
#   Model Features:                                                                     
#   - Three nested tobacco categories via stick-breaking construction                  
#   - Regional hierarchical structure (Global → Region → Country)                      
#   - Smooth sigmoid transition for elderly age extrapolation                          
#   - Dual target evaluation (30% reduction + absolute prevalence)                     
#   - Population-weighted prevalence estimates with full uncertainty                   
#                                                                                       
#   Developed for: World Health Organization
#   Version: 2.3.2 (NIMBLE with vectorized prediction optimizations)
#
#   Changes from v2.3.1:
#   - [OPT-1] Fixed clearCompiled cleanup order (prevents DLL memory leak)
#   - [OPT-2] Precompute age-only components outside country loop
#   - [OPT-3] Vectorize tensor product (eliminates row-by-row outer())
#   - [OPT-4] Vectorize prediction loop (full matrix algebra, 5-10x speedup)
#   - [OPT-5] Vectorize summary statistics (rowMeans + batch quantiles)
#   - [OPT-6] Precompute cohort splines & interactions per-year (single ns() call)
#   - Increased R_MAX_NUM_DLLS to 600 to prevent DLL limit errors
#
#   Changes from v2.2:                                                                 
#   - Migrated from JAGS to NIMBLE for 3-5x performance improvement                   
#   - Compiled C++ backend for faster MCMC sampling                                   
#   - Maintained identical model specification and output structure                    
#                                                                                       
#########################################################################################


#########################################################################################
#                              SECTION 1: SETUP & CONFIGURATION
#########################################################################################

# Record pipeline start time for total runtime calculation
pipeline_start_time <- Sys.time()
cat(sprintf("Pipeline started: %s\n\n", format(pipeline_start_time, "%Y-%m-%d %H:%M:%S")))

# ---- 1.1 Load All Required Packages ----

install_and_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

required_packages <- c(
  # Data manipulation
  "tidyverse", "readxl", "janitor", "purrr", "stringr", "Matrix",
  # Bayesian modeling
  "nimble", "coda",
  # Splines
  "splines",
  # Parallel processing
  "foreach", "doParallel", "parallel",
  # Visualization
  "ggplot2", "scales", "viridis", "ggpubr", "ggrepel", "ggsci", "patchwork", "grid",
  # Tables
  "knitr", "gt",
  # Mapping
  "sf", "rnaturalearth", "rnaturalearthdata",
  # Notification
  "beepr"
)

invisible(lapply(required_packages, install_and_load))

# ---- 1.1.1 Source Custom Modules ----

# Source prediction module for all countries (including those without data)
# Note: Working directory should be project root (who_modeling_r_code/)
if (file.exists("R/prediction_all_countries.R")) {
  cat("  Sourcing R/prediction_all_countries.R...\n")
  source("R/prediction_all_countries.R")
} else {
  cat("  WARNING: R/prediction_all_countries.R not found.\n")
}

# Source MCMC diagnostics module for comprehensive convergence analysis
if (file.exists("R/mcmc_diagnostics.R")) {
  cat("  Sourcing R/mcmc_diagnostics.R...\n")
  source("R/mcmc_diagnostics.R")
} else {
  cat("  WARNING: R/mcmc_diagnostics.R not found. Diagnostics will be limited.\n")
  DIAGNOSTICS_AVAILABLE <- FALSE
}
if (!exists("DIAGNOSTICS_AVAILABLE")) DIAGNOSTICS_AVAILABLE <- TRUE

# ---- 1.2 NIMBLE Configuration ----

# Increase max DLLs to prevent "maximal number of DLLs reached" error
# when fitting 360+ country models (each NIMBLE model compiles to a DLL)
Sys.setenv(R_MAX_NUM_DLLS = 600)

nimbleOptions(verbose = FALSE)
nimbleOptions(MCMCprogressBar = TRUE)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = FALSE)

# Production settings for reliable convergence
NUMBER_OF_CHAINS     <- 4
NUMBER_OF_ADAPT      <- 1000   # Informational only - NIMBLE handles adaptation internally
NUMBER_OF_BURN       <- 5000
NUMBER_OF_ITERATIONS <- 10000
THINNING_INTERVAL    <- 5

# Analysis parameters
BASE_YEAR              <- 2010
TARGET_YEAR            <- 2025
PROJECTED_YEARS        <- 20
REDUCTION_PERCENTAGE   <- 30
RANDOM_SEED            <- 42

# Age transition parameters (smooth sigmoid for elderly)
TRANSITION_START  <- 65
TRANSITION_WIDTH  <- 3
MAX_AGE_SPLINE    <- 80

# Absolute prevalence target
MANUAL_TARGET_PREVALENCE  <- 4.0
MANUAL_TARGET_PROPORTION  <- MANUAL_TARGET_PREVALENCE / 100
MANUAL_TARGET_EVAL_YEAR   <- 2040

# Legacy variable (used throughout code)
target_year <- TARGET_YEAR

# ---- 1.3 Output Directory Structure ----

output_dirs <- c(
  # Model outputs
  "outputs/model_results",
  "outputs/model_results/global",
  "outputs/model_results/country_specific",
  "outputs/predictions",
  "outputs/regional_aggregation",
  "outputs/evaluation",
  "outputs/data_exports",
  # Visualization outputs
  "outputs/figures/world_maps",
  "outputs/figures/world_maps/diagnostics",
  "outputs/figures/publication",
  "outputs/figures/trends/panels",
  "outputs/figures/trends/heatmaps",
  "outputs/figures/validation",
  "outputs/figures/aggregated",
  "outputs/figures/cohorts",
  "outputs/figures/maps",
  "outputs/tables/publication",
  # Processing directories
  "processing",
  "country_priors/males",
  "country_priors/females",
  "compiled_models",
  # Regional aggregation
  "results/regional_aggregation/selected_models",
  # Publication outputs
  "tables_publication",
  "plots/world_maps",
  "plots/birth_cohort_analysis",
  "plots/age_curves",
  # Legacy directories
  "results",
  "evaluation"
)

invisible(lapply(output_dirs, function(d) dir.create(d, recursive = TRUE, showWarnings = FALSE)))

# ---- 1.4 Print Configuration Summary ----

cat("
#########################################################################
#                                                                       
#              WHO TOBACCO CONTROL PREVALENCE PROJECTION                 
#                        VERSION 2.3 (NIMBLE)                            
#                                                                       
#########################################################################

CONFIGURATION:
--------------
  Base Year:           ", BASE_YEAR, "
  Target Year:         ", TARGET_YEAR, "
  Projection Period:   ", PROJECTED_YEARS, " years
  Reduction Target:    ", REDUCTION_PERCENTAGE, "%
  Absolute Target:     ", MANUAL_TARGET_PREVALENCE, "% by ", MANUAL_TARGET_EVAL_YEAR, "

MCMC SETTINGS (NIMBLE):
-----------------------
  Chains:              ", NUMBER_OF_CHAINS, "
  Adaptation:          ", NUMBER_OF_ADAPT, " (built into NIMBLE)
  Burn-in:             ", NUMBER_OF_BURN, "
  Iterations:          ", NUMBER_OF_ITERATIONS, "
  Thinning:            ", THINNING_INTERVAL, "

AGE TRANSITION:
---------------
  Spline-to-Linear:    Age ", TRANSITION_START, " (width: ", TRANSITION_WIDTH, ")
  Max Spline Age:      ", MAX_AGE_SPLINE, "

ENGINE: NIMBLE (C++ compiled MCMC)

#########################################################################
")


#########################################################################################
#                              SECTION 2: DATA LOADING & PREPARATION                    
#########################################################################################

# ---- 2.1 Load WHO Prevalence Data ----

data <- read_excel(
  "data/Prevalence RESHAPE 28 Sep 2024.xlsx",
  sheet = "reshape.dta",
  skip = 1
) %>%
  clean_names()

# ---- 2.2 Clean and Standardize Variable Names ----

clean_data <- data %>%
  rename(
    country        = un_name,
    region         = un_region,
    survey         = title_pri,
    year           = survey_year,
    sex            = sex_name,
    start_age      = start_age,
    end_age        = end_age,
    def_code       = def_code,
    type_code      = type_code,
    type           = published_calculated_estimated_unpublished_microdata,
    prevalence     = prevalence,
    wb_country_abv = wb_country_abv
  ) %>%
  mutate(across(where(is.character), tolower))

# ---- 2.3 Create Indicator Labels ----

clean_data <- clean_data %>%
  mutate(
    def_code_label = case_when(
      def_code == "cd.0112" ~ "current_user",
      def_code == "cd.0101" ~ "daily_user",
      TRUE ~ def_code
    ),
    type_code_label = case_when(
      type_code == "tt.001" ~ "any_tobacco_product",
      type_code == "tt.002" ~ "any_smoked_tobacco",
      type_code == "tt.003" ~ "cigarettes",
      TRUE ~ type_code
    )
  ) %>%
  unite("def_type_code", def_code_label, type_code_label, sep = "_")

# ---- 2.4 Create Binary Indicators ----

clean_data <- clean_data %>%
  mutate(
    def_code_binary    = case_when(
      def_code == "cd.0101" ~ 0,
      def_code == "cd.0112" ~ 1,
      TRUE ~ NA_real_
    ),
    type_any_tobacco   = as.integer(type_code == "tt.001"),
    type_smoked_tobacco = as.integer(type_code == "tt.002"),
    type_cigarettes    = as.integer(type_code == "tt.003"),
    start_age          = as.numeric(start_age),
    end_age            = as.numeric(end_age),
    prevalence         = as.numeric(prevalence),
    year               = as.numeric(year)
  ) %>%
  filter(!is.na(prevalence)) %>%
  mutate(Age_Midpoint = (start_age + end_age) / 2) %>%
  filter(Age_Midpoint >= 15)

# ---- 2.4.1 Validate Tobacco Type Indicators ----

tobacco_type_check <- clean_data %>%
  mutate(sum_type_indicators = type_any_tobacco + type_smoked_tobacco + type_cigarettes) %>%
  group_by(sum_type_indicators) %>%
  summarise(n_obs = n(), .groups = "drop")

cat("\nTobacco type indicator validation:\n")
print(tobacco_type_check)

# Remove observations with unrecognized type codes (all indicators = 0)
n_invalid <- sum(clean_data$type_any_tobacco + clean_data$type_smoked_tobacco + clean_data$type_cigarettes == 0)
if (n_invalid > 0) {
  cat(sprintf("  WARNING: Removing %d observations with unrecognized type_code\n", n_invalid))
  clean_data <- clean_data %>%
    filter(type_any_tobacco + type_smoked_tobacco + type_cigarettes > 0)
}

# ---- 2.5 Remove Redundant Total Prevalence Observations ----

clean_data <- clean_data %>%
  mutate(strata_id = paste(sex, wb_country_abv, def_type_code, year, sep = "_"))

total_prevalence_obs <- clean_data %>%
  filter(start_age == 15 & (end_age >= 100 | end_age == 99))

obs_to_drop <- c()
for (i in 1:nrow(total_prevalence_obs)) {
  current_strata <- total_prevalence_obs$strata_id[i]
  other_obs_in_strata <- clean_data %>%
    filter(strata_id == current_strata) %>%
    filter(!(start_age == 15 & (end_age >= 100 | end_age == 99))) %>%
    nrow()
  
  if (other_obs_in_strata > 0) {
    obs_to_drop <- c(obs_to_drop, which(
      clean_data$strata_id == current_strata &
        clean_data$start_age == 15 &
        (clean_data$end_age >= 100 | clean_data$end_age == 99)
    ))
  }
}

if (length(obs_to_drop) > 0) {
  cat("Removing", length(obs_to_drop), "redundant total prevalence observations\n")
  clean_data <- clean_data[-obs_to_drop, ]
}

clean_data <- clean_data %>% select(-strata_id)

# ---- 2.6 Create Country Name Lookup ----

country_name_mapping <- clean_data %>%
  distinct(wb_country_abv, country) %>%
  arrange(wb_country_abv) %>%
  { setNames(.$country, .$wb_country_abv) }


#########################################################################################
#                              SECTION 3: REGIONAL CLASSIFICATION                       
#########################################################################################

# ---- 3.1 Define WHO Regional Groups ----

country_region_manual <- tribble(
  ~wb_country_abv, ~region_consolidated,
  
  # Central Asia
  "tjk", "Central Asia", "tkm", "Central Asia", "kaz", "Central Asia",
  
  # Eastern Asia
  "chn", "Eastern Asia", "prk", "Eastern Asia", "jpn", "Eastern Asia",
  "mng", "Eastern Asia", "kor", "Eastern Asia",
  # Countries without survey data (will use regional priors):
  "twn", "Eastern Asia",   # Taiwan
  "hkg", "Eastern Asia",   # Hong Kong SAR
  "mac", "Eastern Asia",   # Macao SAR

  # Eastern Europe
  "blr", "Eastern Europe", "bgr", "Eastern Europe", "cze", "Eastern Europe",
  "hun", "Eastern Europe", "mda", "Eastern Europe", "pol", "Eastern Europe",
  "rou", "Eastern Europe", "rus", "Eastern Europe", "svk", "Eastern Europe",
  "ukr", "Eastern Europe",
  
  # Latin America & Caribbean
  "bhs", "Latin America & Caribbean", "brb", "Latin America & Caribbean",
  "blz", "Latin America & Caribbean", "cub", "Latin America & Caribbean",
  "dma", "Latin America & Caribbean", "slv", "Latin America & Caribbean",
  "grd", "Latin America & Caribbean", "gtm", "Latin America & Caribbean",
  "hti", "Latin America & Caribbean", "hnd", "Latin America & Caribbean",
  "jam", "Latin America & Caribbean", "mex", "Latin America & Caribbean",
  "nic", "Latin America & Caribbean", "kna", "Latin America & Caribbean",
  "lca", "Latin America & Caribbean", "vct", "Latin America & Caribbean",
  "tto", "Latin America & Caribbean", "dom", "Latin America & Caribbean",
  "cri", "Latin America & Caribbean", "pan", "Latin America & Caribbean",
  
  # North Africa & Middle East
  "dza", "North Africa & Middle East", "egy", "North Africa & Middle East",
  "lby", "North Africa & Middle East", "mar", "North Africa & Middle East",
  "sdn", "North Africa & Middle East", "tun", "North Africa & Middle East",
  "pse", "North Africa & Middle East",
  
  # North America
  "usa", "North America", "can", "North America",
  
  # Oceania & Pacific
  "kir", "Oceania & Pacific", "mhl", "Oceania & Pacific", "fsm", "Oceania & Pacific",
  "nru", "Oceania & Pacific", "niu", "Oceania & Pacific", "plw", "Oceania & Pacific",
  "png", "Oceania & Pacific", "wsm", "Oceania & Pacific", "slb", "Oceania & Pacific",
  "ton", "Oceania & Pacific", "vut", "Oceania & Pacific", "cok", "Oceania & Pacific",
  "tuv", "Oceania & Pacific", "nzl", "Oceania & Pacific", "aus", "Oceania & Pacific",
  "fji", "Oceania & Pacific",
  
  # South America
  "arg", "South America", "bol", "South America", "bra", "South America",
  "chl", "South America", "col", "South America", "ecu", "South America",
  "guy", "South America", "per", "South America", "sur", "South America",
  "ury", "South America", "ven", "South America", "pry", "South America",
  
  # Southeastern Asia
  "brn", "Southeastern Asia", "idn", "Southeastern Asia", "lao", "Southeastern Asia",
  "mys", "Southeastern Asia", "mmr", "Southeastern Asia", "sgp", "Southeastern Asia",
  "tha", "Southeastern Asia", "vnm", "Southeastern Asia", "phl", "Southeastern Asia",
  "tls", "Southeastern Asia", "khm", "Southeastern Asia",
  
  # Southern Asia
  "afg", "Southern Asia", "bgd", "Southern Asia", "ind", "Southern Asia",
  "irn", "Southern Asia", "kgz", "Southern Asia", "mdv", "Southern Asia",
  "pak", "Southern Asia", "lka", "Southern Asia", "uzb", "Southern Asia",
  "npl", "Southern Asia", "btn", "Southern Asia",
  
  # Southern Europe
  "alb", "Southern Europe", "and", "Southern Europe", "bih", "Southern Europe",
  "grc", "Southern Europe", "hrv", "Southern Europe", "ita", "Southern Europe",
  "mlt", "Southern Europe", "mne", "Southern Europe", "mkd", "Southern Europe",
  "prt", "Southern Europe", "smr", "Southern Europe", "srb", "Southern Europe",
  "svn", "Southern Europe", "esp", "Southern Europe",
  # Countries without survey data (will use regional priors):
  "xkx", "Southern Europe",  # Kosovo

  # Sub-Saharan Africa
  "ago", "Sub-Saharan Africa", "ben", "Sub-Saharan Africa", "bwa", "Sub-Saharan Africa",
  "bfa", "Sub-Saharan Africa", "bdi", "Sub-Saharan Africa", "cpv", "Sub-Saharan Africa",
  "cmr", "Sub-Saharan Africa", "caf", "Sub-Saharan Africa", "tcd", "Sub-Saharan Africa",
  "com", "Sub-Saharan Africa", "cog", "Sub-Saharan Africa", "cod", "Sub-Saharan Africa",
  "civ", "Sub-Saharan Africa", "dji", "Sub-Saharan Africa", "gnq", "Sub-Saharan Africa",
  "eri", "Sub-Saharan Africa", "eth", "Sub-Saharan Africa", "gab", "Sub-Saharan Africa",
  "gmb", "Sub-Saharan Africa", "gha", "Sub-Saharan Africa", "gin", "Sub-Saharan Africa",
  "gnb", "Sub-Saharan Africa", "ken", "Sub-Saharan Africa", "lso", "Sub-Saharan Africa",
  "lbr", "Sub-Saharan Africa", "mdg", "Sub-Saharan Africa", "mwi", "Sub-Saharan Africa",
  "mli", "Sub-Saharan Africa", "mrt", "Sub-Saharan Africa", "mus", "Sub-Saharan Africa",
  "moz", "Sub-Saharan Africa", "nam", "Sub-Saharan Africa", "ner", "Sub-Saharan Africa",
  "nga", "Sub-Saharan Africa", "rwa", "Sub-Saharan Africa", "stp", "Sub-Saharan Africa",
  "sen", "Sub-Saharan Africa", "syc", "Sub-Saharan Africa", "sle", "Sub-Saharan Africa",
  "zaf", "Sub-Saharan Africa", "swz", "Sub-Saharan Africa", "tza", "Sub-Saharan Africa",
  "tgo", "Sub-Saharan Africa", "uga", "Sub-Saharan Africa", "zmb", "Sub-Saharan Africa",
  "zwe", "Sub-Saharan Africa",
  # Countries without survey data (will use regional priors):
  "ssd", "Sub-Saharan Africa",  # South Sudan
  "som", "Sub-Saharan Africa",  # Somalia

  # Western Asia
  "arm", "Western Asia", "aze", "Western Asia", "bhr", "Western Asia",
  "cyp", "Western Asia", "geo", "Western Asia", "irq", "Western Asia",
  "isr", "Western Asia", "jor", "Western Asia", "kwt", "Western Asia",
  "lbn", "Western Asia", "omn", "Western Asia", "qat", "Western Asia",
  "sau", "Western Asia", "syr", "Western Asia", "tur", "Western Asia",
  "are", "Western Asia", "yem", "Western Asia",
  
  # Western Europe
  "aut", "Western Europe", "bel", "Western Europe", "deu", "Western Europe",
  "fra", "Western Europe", "gbr", "Western Europe", "irl", "Western Europe",
  "nld", "Western Europe", "che", "Western Europe", "lux", "Western Europe",
  
  # Northern Europe
  "dnk", "Northern Europe", "fin", "Northern Europe", "isl", "Northern Europe",
  "nor", "Northern Europe", "swe", "Northern Europe", "est", "Northern Europe",
  "lva", "Northern Europe", "ltu", "Northern Europe"
)

# ---- 3.2 Map Countries to Regions ----

country_region_mapping <- clean_data %>%
  select(wb_country_abv) %>%
  distinct() %>%
  left_join(country_region_manual, by = "wb_country_abv") %>%
  mutate(region_consolidated = if_else(is.na(region_consolidated), "Other", region_consolidated))

clean_data <- clean_data %>%
  select(-region) %>%
  left_join(country_region_mapping, by = "wb_country_abv") %>%
  rename(region = region_consolidated)

# ---- 3.3 Export Regional Mapping ----

write.csv(
  country_region_manual %>%
    mutate(Country_Name = country_name_mapping[wb_country_abv]) %>%
    select(Country_Code = wb_country_abv, Country_Name, Region = region_consolidated) %>%
    arrange(Region, Country_Name),
  file = "outputs/data_exports/country_region_mapping.csv",
  row.names = FALSE
)


#########################################################################################
#                              SECTION 4: AGE-COHORT TRANSFORMATION                     
#########################################################################################

# ---- 4.1 Calculate Birth Cohorts ----

clean_data <- clean_data %>%
  mutate(
    Birth_Cohort    = year - Age_Midpoint,
    Cohort_Centered = Birth_Cohort - median(Birth_Cohort, na.rm = TRUE)
  )

# ---- 4.2 Create Smooth Sigmoid Transition Weights ----

clean_data <- clean_data %>%
  mutate(
    sigmoid_input       = (Age_Midpoint - TRANSITION_START) / TRANSITION_WIDTH,
    spline_weight       = 1 / (1 + exp(sigmoid_input)),
    linear_weight       = 1 - spline_weight,
    Age_For_Spline      = Age_Midpoint,
    Age_Linear          = pmax(0, Age_Midpoint - TRANSITION_START),
    Age_Linear_Centered = Age_Linear - mean(Age_Linear, na.rm = TRUE)
  )

# ---- 4.2.1 Store Global Centering Constants (FIX from v2.2) ----

AGE_LINEAR_CENTER_CONSTANT <- mean(clean_data$Age_Linear, na.rm = TRUE)
COHORT_CENTER_CONSTANT     <- median(clean_data$Birth_Cohort, na.rm = TRUE)

cat("\nGlobal centering constants stored:\n")
cat("  AGE_LINEAR_CENTER_CONSTANT:", AGE_LINEAR_CENTER_CONSTANT, "\n")
cat("  COHORT_CENTER_CONSTANT:    ", COHORT_CENTER_CONSTANT, "\n\n")

# ---- 4.3 Define Spline Knots ----
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

# ---- 4.4 Create Spline Bases ----

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

# ---- 4.5 Create Weight Variables ----

clean_data$age_linear_smooth   <- clean_data$Age_Linear_Centered
clean_data$spline_weight_var   <- clean_data$spline_weight
clean_data$linear_weight_var   <- clean_data$linear_weight

# ---- 4.6 Build Age-Cohort Interaction Tensor ----

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

# ---- 4.7 Transform Prevalence to Logit Scale ----

clean_data$prevalence <- clean_data$prevalence / 100
clean_data$prevalence <- ifelse(clean_data$prevalence == 0, 0.001, clean_data$prevalence)
clean_data$prevalence <- ifelse(clean_data$prevalence == 1, 0.999, clean_data$prevalence)
clean_data$prevalence <- log(clean_data$prevalence / (1 - clean_data$prevalence))

# ---- 4.8 Create Product Type Indicators ----

clean_data <- clean_data %>%
  mutate(
    Type_Cig    = as.integer(type_cigarettes == 1),
    Type_Smoked = as.integer(type_smoked_tobacco == 1),
    Type_Any    = as.integer(type_any_tobacco == 1)
  )


#########################################################################################
#                              SECTION 5: POPULATION WEIGHTS                            
#########################################################################################

# ---- 5.1 Load Population Weight Data ----

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


#########################################################################################
#                              SECTION 6: UTILITY FUNCTIONS                             
#########################################################################################

# ---- 6.1 Check Ordering Constraints ----

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

# ---- 6.2 Calculate Population-Weighted Prevalence ----

calculate_weighted_prevalence <- function(iterations_matrix_logit, weights_vector) {
  p_matrix <- plogis(iterations_matrix_logit)
  weighted_probs <- numeric(ncol(p_matrix))
  
  for (s in 1:ncol(p_matrix)) {
    weighted_probs[s] <- sum(p_matrix[, s] * weights_vector) / sum(weights_vector)
  }
  
  return(weighted_probs)
}

# ---- 6.3 Precision to SD Conversion Helper ----

prec_to_sd <- function(precision) {
  1 / sqrt(pmax(precision, 1e-10))
}


#########################################################################################
#                              SECTION 7: NIMBLE MODEL DEFINITIONS                      
#########################################################################################

# ---- 7.1 Global Hierarchical Model (Stick-Breaking Construction) - NIMBLE ----
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
  
  # SMKEXTRA: Small component - P(other smoked | not cig) ≈ 3-5%
  smkextra_global_intercept ~ dnorm(-3.0, sd = 0.3)
  
  # ANYEXTRA: Very small component - P(smokeless | not smoked) ≈ 1-3%
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


# ---- 7.2 Country-Specific Model (With Informative Priors) - NIMBLE ----
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

#########################################################################################
#                              SECTION 8: GLOBAL MODEL FITTING (NIMBLE)                 
#########################################################################################

genders <- c("males", "females")


gender_results <- list()

for (gender in genders) {
  
  cat("\n")
  cat("================================================================\n")
  cat("  FITTING GLOBAL HIERARCHICAL MODEL (NIMBLE):", toupper(gender), "\n")
  cat("================================================================\n")
  
  # ---- 8.1 Setup Directories ----
  
  gender_dir <- file.path("processing", gender)
  dir.create(gender_dir, showWarnings = FALSE)
  
  gender_priors_dir <- file.path("country_priors", gender)
  dir.create(gender_priors_dir, showWarnings = FALSE)
  
  # ---- 8.2 Filter Gender Data ----
  
  gender_data <- clean_data %>% filter(sex == gender)
  
  # ---- 8.3 Create Region Mappings ----
  
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
  
  # ---- 8.4 Calculate Observation Weights ----
  
  gender_data <- gender_data %>%
    mutate(
      age_range = end_age - start_age,
      weight    = 1 / (age_range + 1),
      weight    = weight / mean(weight)
    )
  
  # ---- 8.5 Prepare NIMBLE Constants and Data (KEY CHANGE FROM JAGS) ----
  
  # NIMBLE requires separating constants (structure) from data (observations)
  
  # ==================================================================
  # KEY FIX: Compute empirical mean for informative prior
  # This anchors the global intercept, breaking non-identifiability
  # ==================================================================
  empirical_mean_cig <- mean(gender_data$prevalence)
  
  cat(sprintf("\n  Empirical mean prevalence (logit scale): %.3f\n", empirical_mean_cig))
  cat(sprintf("  Equivalent to approximately %.1f%% prevalence\n", 100 * plogis(empirical_mean_cig)))
  
  # Sanity check
  
  if (empirical_mean_cig < -4 || empirical_mean_cig > 0) {
    warning("  WARNING: Empirical mean outside expected range [-4, 0] - check data transformation")
  }
  
  # CRITICAL FIX: Create Country_Region as length nCountry, not length N
  # This maps each country index to its region index
  country_region_for_model <- country_lookup_df %>%
    arrange(Num_Country) %>%
    pull(Region_Num)
  
  # CONSTANTS: Define model structure (loop bounds, indices)
  nimble_constants <- list(
    N = nrow(gender_data),
    nCountry = length(unique(gender_data$Num_Country)),
    nRegion = length(unique_regions),
    nSurvey = length(unique(gender_data$Num_Survey)),
    nAgeSpline = ncol(gender_data[, grep("^age_spline_", names(gender_data))]),
    nCohortSpline = ncol(gender_data[, grep("^cohort_spline_", names(gender_data))]),
    nAgeXCohortSplines = ncol(gender_data[, grep("^age_cohort_", names(gender_data))]),
    
    # Index vectors
    Country = as.integer(gender_data$Num_Country),
    # FIXED: Country_Region is now length nCountry, indexed by country number
    Country_Region = as.integer(country_region_for_model),
    Survey = as.integer(gender_data$Num_Survey),
    
    # NEW: Empirical mean for informative prior on cig_global_intercept
    empirical_mean_cig = empirical_mean_cig
  )
  
  # Verify the fix
  cat(sprintf("  Country_Region length: %d (should equal nCountry: %d)\n", 
              length(nimble_constants$Country_Region), 
              nimble_constants$nCountry))
  
  # DATA: Observed values and covariates
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
  
  # ---- 8.6 Save Regional Structure ----
  
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
  
  # ---- 8.7 Generate Initial Values for NIMBLE ----
  # UPDATED: Tighter initial values, chains start closer together
  
  generate_global_inits <- function(seed_offset = 0, emp_mean = -1.5) {
    set.seed(RANDOM_SEED + seed_offset)
    list(
      # Global intercepts - start NEAR the prior means with small jitter
      # This ensures chains start in similar regions for faster mixing
      cig_global_intercept = emp_mean + rnorm(1, 0, 0.05),
      smkextra_global_intercept = -3.0 + rnorm(1, 0, 0.05),
      anyextra_global_intercept = -4.0 + rnorm(1, 0, 0.05),
      
      # Def code effect
      cig_def_code_shared = rnorm(1, 0.3, 0.05),
      
      # Residual SD
      residual_sd = 0.7 + runif(1, -0.05, 0.05),
      
      # Age linear effects (negative, constrained) - tighter
      cig_age_linear_smooth_effect = -0.02 + rnorm(1, 0, 0.005),
      smkextra_age_linear_smooth_effect = -0.02 + rnorm(1, 0, 0.005),
      anyextra_age_linear_smooth_effect = -0.02 + rnorm(1, 0, 0.005),
      
      # Precision parameters - start near prior means
      intercept_between_region_precision = 4 + rnorm(1, 0, 0.3),
      age_spline_between_region_precision = 4 + rnorm(1, 0, 0.3),
      cohort_spline_between_region_precision = 4 + rnorm(1, 0, 0.3),
      survey_intercept_precision = 3 + rnorm(1, 0, 0.3),
      
      # Shared precision for smkextra/anyextra intercepts
      smkextra_intercept_within_region_precision = 4 + rnorm(1, 0, 0.3),
      anyextra_intercept_within_region_precision = 4 + rnorm(1, 0, 0.3)
    )
  }
  
  # Generate initial values using the computed empirical mean
  inits_list <- lapply(1:NUMBER_OF_CHAINS, function(i) {
    generate_global_inits(i, emp_mean = empirical_mean_cig)
  })

  # ---- 8.8 Build NIMBLE Model ----
  
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
  
  # ---- 8.9 Configure MCMC ----
  
  cat("  Configuring MCMC...\n")
  
  mcmc_config <- configureMCMC(
    nimble_model,
    monitors = c(
      # Global parameters
      "cig_global_intercept", "cig_def_code_shared",
      "smkextra_global_intercept", "anyextra_global_intercept",
      
      # Regional intercepts (all heads)
      "cig_region_intercept", "smkextra_region_intercept", "anyextra_region_intercept",
      
      # Regional splines (all heads - smkextra/anyextra use these directly)
      "cig_age_spline_region_mean", "smkextra_age_spline_region_mean", "anyextra_age_spline_region_mean",
      "cig_cohort_spline_region_mean", "smkextra_cohort_spline_region_mean", "anyextra_cohort_spline_region_mean",
      
      # Country intercepts (all heads)
      "cig_country_intercept", "smkextra_country_intercept", "anyextra_country_intercept",
      
      # Country splines (CIG ONLY)
      "cig_age_spline", "cig_cohort_spline",
      
      # Age-cohort interactions (CIG ONLY)
      "cig_age_cohort_interaction",
      
      # Linear age effects (all heads)
      "cig_age_linear_smooth_effect", "smkextra_age_linear_smooth_effect", "anyextra_age_linear_smooth_effect",
      
      # Global spline means (all heads)
      "cig_age_spline_global_mean", "cig_cohort_spline_global_mean",
      "smkextra_age_spline_global_mean", "smkextra_cohort_spline_global_mean",
      "anyextra_age_spline_global_mean", "anyextra_cohort_spline_global_mean"
    ),
    thin = THINNING_INTERVAL,
    enableWAIC = FALSE,
    useConjugacy = FALSE  # CHANGED: Enable conjugacy for faster sampling
  )
  
  # ---- 8.10 Build and Compile MCMC ----
  
  cat("  Building MCMC...\n")
  mcmc_built <- buildMCMC(mcmc_config)
  
  cat("  Compiling model and MCMC (this may take several minutes)...\n")
  start_compile <- Sys.time()
  
  compiled_model <- compileNimble(nimble_model)
  compiled_mcmc <- compileNimble(mcmc_built, project = nimble_model)
  
  compile_time <- difftime(Sys.time(), start_compile, units = "mins")
  cat(sprintf("  Compilation time: %.1f minutes\n", compile_time))
  
  # Save compiled model info
  saveRDS(
    list(compile_time = compile_time, n_params = length(nimble_model$getNodeNames())),
    file.path("compiled_models", paste0("global_model_info_", gender, ".rds"))
  )
  
  # ---- 8.11 Run MCMC ----
  
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

  # ---- 8.12 Comprehensive Convergence Diagnostics ----

  # Combine samples from all chains first (needed for both diagnostics and downstream)
  combined_samples_matrix <- do.call(rbind, lapply(samples$samples, as.matrix))

  # Run comprehensive diagnostics if module available
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
      cat(sprintf("  Parameters > 1.1: %d (%.1f%%)\n",
                  global_diagnostics$convergence_summary$n_rhat_above_1.1,
                  global_diagnostics$convergence_summary$pct_rhat_above_1.1))
      cat(sprintf("  Mean ESS: %.0f | Min ESS: %.0f\n",
                  global_diagnostics$convergence_summary$mean_ess,
                  global_diagnostics$convergence_summary$min_ess))
      cat(sprintf("  Diagnostics saved to: %s\n",
                  global_diagnostics$log_filepath))
    }
  } else {
    # Fallback to basic convergence check
    gelman_diag <- try(gelman.diag(samples$samples, multivariate = FALSE), silent = TRUE)
    if (!inherits(gelman_diag, "try-error")) {
      max_rhat <- max(gelman_diag$psrf[, "Point est."], na.rm = TRUE)
      cat(sprintf("  Convergence: Max R-hat = %.3f %s\n", max_rhat,
                  ifelse(max_rhat < 1.1, "(Good)", "(Warning)")))
    }
  }
  
  # ---- 8.13 Extract Country-Specific Priors (FIX: No double-counting) ----
  
  cig_global_intercept    <- mean(combined_samples_matrix[, "cig_global_intercept"])
  cig_global_intercept_sd <- sd(combined_samples_matrix[, "cig_global_intercept"])
  
  smkextra_global_intercept    <- mean(combined_samples_matrix[, "smkextra_global_intercept"])
  smkextra_global_intercept_sd <- sd(combined_samples_matrix[, "smkextra_global_intercept"])
  
  anyextra_global_intercept    <- mean(combined_samples_matrix[, "anyextra_global_intercept"])
  anyextra_global_intercept_sd <- sd(combined_samples_matrix[, "anyextra_global_intercept"])
  
  def_code_shared_mean <- mean(combined_samples_matrix[, "cig_def_code_shared"])
  def_code_shared_sd   <- sd(combined_samples_matrix[, "cig_def_code_shared"])
  
  cig_age_linear_smooth    <- mean(combined_samples_matrix[, "cig_age_linear_smooth_effect"])
  cig_age_linear_smooth_sd <- sd(combined_samples_matrix[, "cig_age_linear_smooth_effect"])
  
  smkextra_age_linear_smooth    <- mean(combined_samples_matrix[, "smkextra_age_linear_smooth_effect"])
  smkextra_age_linear_smooth_sd <- sd(combined_samples_matrix[, "smkextra_age_linear_smooth_effect"])
  
  anyextra_age_linear_smooth    <- mean(combined_samples_matrix[, "anyextra_age_linear_smooth_effect"])
  anyextra_age_linear_smooth_sd <- sd(combined_samples_matrix[, "anyextra_age_linear_smooth_effect"])
  
  # CIG age-cohort interactions (CIG ONLY - smkextra/anyextra no longer have these)
  cig_age_cohort_cols  <- grep("cig_age_cohort_interaction\\[", colnames(combined_samples_matrix))
  cig_age_cohort_means <- colMeans(combined_samples_matrix[, cig_age_cohort_cols])
  cig_age_cohort_sds   <- apply(combined_samples_matrix[, cig_age_cohort_cols], 2, sd)
  
  # REMOVED: smkextra_age_cohort_* and anyextra_age_cohort_* - these parameters no longer exist
  
  # ---- 8.14 Save Country Priors (CORRECTED v2.3.1) ----
  # 
  # FIXES APPLIED:
  # 1. Total intercept = global + region + country (was missing region)
  # 2. CIG splines use country-specific posteriors (was using regional means)
  # 3. Uncertainty computed sample-wise (was assuming independence)
  # 4. Removed smkextra/anyextra age-cohort (not in global model)
  
  cat("  Extracting country priors (corrected method)...\n")
  
  for (country_num in names(country_mapping)) {
    country_code      <- country_mapping[country_num]
    country_full_name <- country_name_mapping[country_code]
    country_region    <- country_lookup_df$Region_Num[country_lookup_df$Num_Country == as.numeric(country_num)]
    region_name       <- unique_regions[country_region]
    
    # ==================================================================
    # HELPER FUNCTION: Compute total intercept SAMPLE-WISE
    # This accounts for correlations between global/regional/country
    # Total = global + region + country
    # ==================================================================
    
    compute_total_intercept_samplewise <- function(head_prefix) {
      global_col  <- paste0(head_prefix, "_global_intercept")
      region_col  <- paste0(head_prefix, "_region_intercept[", country_region, "]")
      country_col <- paste0(head_prefix, "_country_intercept[", country_num, "]")
      
      # Validate columns exist
      required_cols <- c(global_col, region_col, country_col)
      missing_cols <- required_cols[!required_cols %in% colnames(combined_samples_matrix)]
      
      if (length(missing_cols) > 0) {
        warning(sprintf("Missing columns for %s total intercept: %s", 
                        head_prefix, paste(missing_cols, collapse = ", ")))
        return(list(mean = NA_real_, sd = NA_real_))
      }
      
      # SAMPLE-WISE computation: sum for each MCMC sample, then compute statistics
      # This preserves correlations!
      total_samples <- combined_samples_matrix[, global_col] +
        combined_samples_matrix[, region_col] +
        combined_samples_matrix[, country_col]
      
      list(
        mean = mean(total_samples),
        sd   = sd(total_samples)
      )
    }
    
    # ==================================================================
    # CIG: Total intercept (global + region + country)
    # ==================================================================
    cig_intercept_stats <- compute_total_intercept_samplewise("cig")
    
    # ==================================================================
    # CIG: COUNTRY-SPECIFIC splines (NOT regional!)
    # The global model has cig_age_spline[j, l] for each country j
    # ==================================================================
    cig_age_spline_means <- c()
    cig_age_spline_sds   <- c()
    for (l in 1:nimble_constants$nAgeSpline) {
      # Use country-specific posterior (this is the key fix!)
      country_col <- paste0("cig_age_spline[", country_num, ", ", l, "]")
      if (country_col %in% colnames(combined_samples_matrix)) {
        cig_age_spline_means[l] <- mean(combined_samples_matrix[, country_col])
        cig_age_spline_sds[l]   <- sd(combined_samples_matrix[, country_col])
      } else {
        # Fallback to regional mean only if country-specific not found
        regional_col <- paste0("cig_age_spline_region_mean[", country_region, ", ", l, "]")
        if (regional_col %in% colnames(combined_samples_matrix)) {
          cig_age_spline_means[l] <- mean(combined_samples_matrix[, regional_col])
          cig_age_spline_sds[l]   <- sd(combined_samples_matrix[, regional_col])
          warning(sprintf("Using regional fallback for cig_age_spline[%s, %d]", country_num, l))
        } else {
          cig_age_spline_means[l] <- 0
          cig_age_spline_sds[l]   <- 1
          warning(sprintf("No data for cig_age_spline[%s, %d], using default", country_num, l))
        }
      }
    }
    
    cig_cohort_spline_means <- c()
    cig_cohort_spline_sds   <- c()
    for (m in 1:nimble_constants$nCohortSpline) {
      # Use country-specific posterior
      country_col <- paste0("cig_cohort_spline[", country_num, ", ", m, "]")
      if (country_col %in% colnames(combined_samples_matrix)) {
        cig_cohort_spline_means[m] <- mean(combined_samples_matrix[, country_col])
        cig_cohort_spline_sds[m]   <- sd(combined_samples_matrix[, country_col])
      } else {
        # Fallback to regional mean
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
    
    # ==================================================================
    # SMKEXTRA: Total intercept (global + region + country)
    # ==================================================================
    smkextra_intercept_stats <- compute_total_intercept_samplewise("smkextra")
    
    # ==================================================================
    # SMKEXTRA: REGIONAL splines (correct - global model uses regional directly)
    # ==================================================================
    smkextra_age_spline_means <- c()
    smkextra_age_spline_sds   <- c()
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
    
    smkextra_cohort_spline_means <- c()
    smkextra_cohort_spline_sds   <- c()
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
    
    # ==================================================================
    # ANYEXTRA: Total intercept (global + region + country)
    # ==================================================================
    anyextra_intercept_stats <- compute_total_intercept_samplewise("anyextra")
    
    # ==================================================================
    # ANYEXTRA: REGIONAL splines (correct - global model uses regional directly)
    # ==================================================================
    anyextra_age_spline_means <- c()
    anyextra_age_spline_sds   <- c()
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
    
    anyextra_cohort_spline_means <- c()
    anyextra_cohort_spline_sds   <- c()
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
    
    # ==================================================================
    # CIG Age-cohort interactions (global parameters, same for all countries)
    # ==================================================================
    cig_age_cohort_cols <- grep("^cig_age_cohort_interaction\\[", colnames(combined_samples_matrix))
    if (length(cig_age_cohort_cols) > 0) {
      cig_age_cohort_means <- colMeans(combined_samples_matrix[, cig_age_cohort_cols, drop = FALSE])
      cig_age_cohort_sds   <- apply(combined_samples_matrix[, cig_age_cohort_cols, drop = FALSE], 2, sd)
    } else {
      n_interactions <- nimble_constants$nAgeXCohortSplines
      cig_age_cohort_means <- rep(0, n_interactions)
      cig_age_cohort_sds   <- rep(0.1, n_interactions)
    }
    
    # ==================================================================
    # Age linear effects (global parameters)
    # ==================================================================
    cig_age_linear_mean <- mean(combined_samples_matrix[, "cig_age_linear_smooth_effect"])
    cig_age_linear_sd   <- sd(combined_samples_matrix[, "cig_age_linear_smooth_effect"])
    
    smkextra_age_linear_mean <- mean(combined_samples_matrix[, "smkextra_age_linear_smooth_effect"])
    smkextra_age_linear_sd   <- sd(combined_samples_matrix[, "smkextra_age_linear_smooth_effect"])
    
    anyextra_age_linear_mean <- mean(combined_samples_matrix[, "anyextra_age_linear_smooth_effect"])
    anyextra_age_linear_sd   <- sd(combined_samples_matrix[, "anyextra_age_linear_smooth_effect"])
    
    # ==================================================================
    # Definition code effect (global parameter)
    # ==================================================================
    def_code_mean <- mean(combined_samples_matrix[, "cig_def_code_shared"])
    def_code_sd   <- sd(combined_samples_matrix[, "cig_def_code_shared"])
    
    # ==================================================================
    # BUILD PRIORS DATAFRAME
    # ==================================================================
    
    # Shrinkage factor: tightens priors slightly to give country model 
    # some room to adjust while regularizing toward global estimates
    shrinkage <- 0.7
    
    # Minimum SDs to prevent numerical issues (precision explosion)
    min_sd_intercept   <- 0.05
    min_sd_spline      <- 0.02
    min_sd_interaction <- 0.01
    min_sd_linear      <- 0.002
    
    # CORRECTED: No smkextra/anyextra age-cohort (not in global model)
    all_priors <- data.frame(
      parameter = c(
        # Definition code effect
        "cig_def_code_shared",
        # CIG parameters
        "cig_intercept", 
        "cig_age_linear_smooth_effect",
        paste0("cig_age_spline_", 1:length(cig_age_spline_means)),
        paste0("cig_cohort_spline_", 1:length(cig_cohort_spline_means)),
        paste0("cig_age_cohort_interaction_", 1:length(cig_age_cohort_means)),
        # SMKEXTRA parameters (NO age-cohort)
        "smkextra_intercept", 
        "smkextra_age_linear_smooth_effect",
        paste0("smkextra_age_spline_", 1:length(smkextra_age_spline_means)),
        paste0("smkextra_cohort_spline_", 1:length(smkextra_cohort_spline_means)),
        # ANYEXTRA parameters (NO age-cohort)
        "anyextra_intercept", 
        "anyextra_age_linear_smooth_effect",
        paste0("anyextra_age_spline_", 1:length(anyextra_age_spline_means)),
        paste0("anyextra_cohort_spline_", 1:length(anyextra_cohort_spline_means))
      ),
      mean = c(
        # Definition code
        def_code_mean,
        # CIG
        cig_intercept_stats$mean,
        cig_age_linear_mean,
        cig_age_spline_means,
        cig_cohort_spline_means,
        cig_age_cohort_means,
        # SMKEXTRA
        smkextra_intercept_stats$mean,
        smkextra_age_linear_mean,
        smkextra_age_spline_means,
        smkextra_cohort_spline_means,
        # ANYEXTRA
        anyextra_intercept_stats$mean,
        anyextra_age_linear_mean,
        anyextra_age_spline_means,
        anyextra_cohort_spline_means
      ),
      sd = c(
        # Definition code
        pmax(def_code_sd * shrinkage, min_sd_intercept),
        # CIG
        pmax(cig_intercept_stats$sd * shrinkage, min_sd_intercept),
        pmax(cig_age_linear_sd * shrinkage, min_sd_linear),
        pmax(cig_age_spline_sds * shrinkage, min_sd_spline),
        pmax(cig_cohort_spline_sds * shrinkage, min_sd_spline),
        pmax(cig_age_cohort_sds * shrinkage, min_sd_interaction),
        # SMKEXTRA
        pmax(smkextra_intercept_stats$sd * shrinkage, min_sd_intercept),
        pmax(smkextra_age_linear_sd * shrinkage, min_sd_linear),
        pmax(smkextra_age_spline_sds * shrinkage, min_sd_spline),
        pmax(smkextra_cohort_spline_sds * shrinkage, min_sd_spline),
        # ANYEXTRA
        pmax(anyextra_intercept_stats$sd * shrinkage, min_sd_intercept),
        pmax(anyextra_age_linear_sd * shrinkage, min_sd_linear),
        pmax(anyextra_age_spline_sds * shrinkage, min_sd_spline),
        pmax(anyextra_cohort_spline_sds * shrinkage, min_sd_spline)
      ),
      stringsAsFactors = FALSE
    )
    
    # Validate: Check for NA values and replace with defaults
    na_means <- is.na(all_priors$mean)
    na_sds   <- is.na(all_priors$sd)
    
    if (any(na_means) || any(na_sds)) {
      warning(sprintf("NA values in priors for %s: %d means, %d sds", 
                      country_code, sum(na_means), sum(na_sds)))
      all_priors$mean[na_means] <- 0
      all_priors$sd[na_sds]     <- 1
    }
    
    # Save to CSV
    write.csv(
      all_priors,
      file = file.path(gender_priors_dir, paste0(country_code, "_regional_ac_priors_nested.csv")),
      row.names = FALSE
    )
  }
  
  cat("  Country priors saved (corrected method).\n")
  
  # ---- 8.15 Generate Predictions (using global centering constants) ----
  # CORRECTED: Added regional intercept extraction
  
  # Number of posterior samples for predictions

  # Higher = smoother uncertainty estimates, but slower
  # 50 = fast testing, 1000 = publication quality
  n_samples <- min(1000, nrow(combined_samples_matrix))
  sampled_indices <- sort(sample(nrow(combined_samples_matrix), n_samples))
  
  # CIG: Full extraction (country splines + interactions + REGIONAL INTERCEPT)
  cig_samples <- list(
    global_intercept  = combined_samples_matrix[sampled_indices, "cig_global_intercept"],
    region_intercept  = combined_samples_matrix[sampled_indices, grep("^cig_region_intercept\\[", colnames(combined_samples_matrix))],
    age_linear_smooth = combined_samples_matrix[sampled_indices, "cig_age_linear_smooth_effect"],
    country_intercept = combined_samples_matrix[sampled_indices, grep("^cig_country_intercept\\[", colnames(combined_samples_matrix))],
    age_spline        = combined_samples_matrix[sampled_indices, grep("^cig_age_spline\\[", colnames(combined_samples_matrix))],
    cohort_spline     = combined_samples_matrix[sampled_indices, grep("^cig_cohort_spline\\[", colnames(combined_samples_matrix))],
    age_cohort        = combined_samples_matrix[sampled_indices, grep("^cig_age_cohort_interaction\\[", colnames(combined_samples_matrix))]
  )
  
  # SMKEXTRA: Simplified (country intercepts + REGIONAL splines + REGIONAL INTERCEPT)
  smkextra_samples <- list(
    global_intercept     = combined_samples_matrix[sampled_indices, "smkextra_global_intercept"],
    region_intercept     = combined_samples_matrix[sampled_indices, grep("^smkextra_region_intercept\\[", colnames(combined_samples_matrix))],
    age_linear_smooth    = combined_samples_matrix[sampled_indices, "smkextra_age_linear_smooth_effect"],
    country_intercept    = combined_samples_matrix[sampled_indices, grep("^smkextra_country_intercept\\[", colnames(combined_samples_matrix))],
    age_spline_regional  = combined_samples_matrix[sampled_indices, grep("^smkextra_age_spline_region_mean\\[", colnames(combined_samples_matrix))],
    cohort_spline_regional = combined_samples_matrix[sampled_indices, grep("^smkextra_cohort_spline_region_mean\\[", colnames(combined_samples_matrix))]
  )
  
  # ANYEXTRA: Simplified (country intercepts + REGIONAL splines + REGIONAL INTERCEPT)
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

  # ==================================================================
  # ALL-COUNTRIES PREDICTION: Include countries without survey data
  # ==================================================================

  # Create lookup for ALL countries (from country_region_manual)
  all_countries_master <- unique(country_region_manual$wb_country_abv)
  data_country_codes <- unique(country_lookup_df$wb_country_abv)

  # Create comprehensive lookup including countries without data
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

  # Report statistics
  n_data_countries <- sum(all_country_lookup$has_data)
  n_nodata_countries <- sum(!all_country_lookup$has_data)

  cat(sprintf("\n  === ALL-COUNTRIES PREDICTION ===\n"))
  cat(sprintf("  Total countries in master list: %d\n", nrow(all_country_lookup)))
  cat(sprintf("  Countries WITH survey data:     %d (use fitted MCMC parameters)\n", n_data_countries))
  cat(sprintf("  Countries WITHOUT survey data:  %d (sample from hierarchical prior)\n", n_nodata_countries))

  # List countries without data by region
  nodata_by_region <- all_country_lookup %>%
    filter(!has_data) %>%
    group_by(region_consolidated) %>%
    summarise(
      n = n(),
      countries = paste(wb_country_abv, collapse = ", "),
      .groups = "drop"
    )

  if (nrow(nodata_by_region) > 0) {
    cat("\n  Countries without data by region:\n")
    for (i in 1:nrow(nodata_by_region)) {
      cat(sprintf("    %s (%d): %s\n",
                  nodata_by_region$region_consolidated[i],
                  nodata_by_region$n[i],
                  nodata_by_region$countries[i]))
    }
  }

  # Check for countries with missing region mapping
  missing_region <- all_country_lookup %>%
    filter(is.na(Region_Num))

  if (nrow(missing_region) > 0) {
    warning(sprintf("  %d countries have unmapped regions and will be skipped: %s",
                    nrow(missing_region),
                    paste(missing_region$wb_country_abv, collapse = ", ")))
  }

  # Create prediction grid for DATA countries (these use Num_Country indexing)
  new_data_age_period <- expand.grid(
    Year         = years_for_prediction,
    Def_Code_Binary = c(0, 1),
    Age_Midpoint = age_midpoints,
    Num_Country  = unique_countries
  )
  
  # FIX: Use global centering constants
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
  
  # ---- 8.16 Parallel Prediction ----
  # CORRECTED: Added regional intercept to prediction formulas
  
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
    
    # Get region for this country (needed for regional intercepts and smkextra/anyextra splines)
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
        
        # CIG: CORRECTED - includes regional intercept
        mu_cig <- cig_samples$global_intercept +
          cig_samples$region_intercept[, paste0("cig_region_intercept[", country_region, "]")] +
          cig_samples$country_intercept[, paste0("cig_country_intercept[", country, "]")] +
          def_code_shared_samples * def_code_binary +
          spline_weight_value * as.matrix(cig_samples$age_spline[, paste0("cig_age_spline[", country, ", ", 1:length(age_spline_values), "]")]) %*% age_spline_values +
          linear_weight_value * cig_samples$age_linear_smooth * age_linear_smooth +
          as.matrix(cig_samples$cohort_spline[, paste0("cig_cohort_spline[", country, ", ", 1:length(cohort_spline_values), "]")]) %*% cohort_spline_values +
          as.matrix(cig_samples$age_cohort) %*% age_cohort_interaction_values
        
        # SMKEXTRA: CORRECTED - includes regional intercept
        mu_smkextra <- smkextra_samples$global_intercept +
          smkextra_samples$region_intercept[, paste0("smkextra_region_intercept[", country_region, "]")] +
          smkextra_samples$country_intercept[, paste0("smkextra_country_intercept[", country, "]")] +
          0.3 * def_code_shared_samples * def_code_binary +
          spline_weight_value * as.matrix(smkextra_samples$age_spline_regional[, paste0("smkextra_age_spline_region_mean[", country_region, ", ", 1:length(age_spline_values), "]")]) %*% age_spline_values +
          linear_weight_value * smkextra_samples$age_linear_smooth * age_linear_smooth +
          as.matrix(smkextra_samples$cohort_spline_regional[, paste0("smkextra_cohort_spline_region_mean[", country_region, ", ", 1:length(cohort_spline_values), "]")]) %*% cohort_spline_values
        
        # ANYEXTRA: CORRECTED - includes regional intercept
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
      
      saveRDS(predictions_cig, file = file.path(cig_dir, paste0(current_year, ".rds")), compress = "xz")
      saveRDS(predictions_smoked, file = file.path(smoked_dir, paste0(current_year, ".rds")), compress = "xz")
      saveRDS(predictions_any, file = file.path(any_dir, paste0(current_year, ".rds")), compress = "xz")
      
      current_data_cig <- current_data %>%
        mutate(
          Prevalence    = apply(predictions_cig, 1, function(x) mean(plogis(x))),
          lower_ci      = apply(predictions_cig, 1, function(x) quantile(plogis(x), 0.025)),
          upper_ci      = apply(predictions_cig, 1, function(x) quantile(plogis(x), 0.975)),
          Def_Type_Code = paste0(current_def, "_cigarettes")
        )
      
      current_data_smoked <- current_data %>%
        mutate(
          Prevalence    = apply(predictions_smoked, 1, function(x) mean(plogis(x))),
          lower_ci      = apply(predictions_smoked, 1, function(x) quantile(plogis(x), 0.025)),
          upper_ci      = apply(predictions_smoked, 1, function(x) quantile(plogis(x), 0.975)),
          Def_Type_Code = paste0(current_def, "_any_smoked_tobacco")
        )
      
      current_data_any <- current_data %>%
        mutate(
          Prevalence    = apply(predictions_any, 1, function(x) mean(plogis(x))),
          lower_ci      = apply(predictions_any, 1, function(x) quantile(plogis(x), 0.025)),
          upper_ci      = apply(predictions_any, 1, function(x) quantile(plogis(x), 0.975)),
          Def_Type_Code = paste0(current_def, "_any_tobacco_product")
        )
      
      results_list[[length(results_list) + 1]] <- current_data_cig
      results_list[[length(results_list) + 1]] <- current_data_smoked
      results_list[[length(results_list) + 1]] <- current_data_any
      
      rm(predictions_cig, predictions_smoked, predictions_any)
      gc()
    }
    
    result         <- do.call(rbind, results_list)
    result$Country <- country_code
    result$Sex     <- gender
    result$Has_Survey_Data <- TRUE

    return(result %>% select(
      Year, Age_Midpoint, Birth_Cohort, Def_Type_Code,
      Sex, Country, Prevalence, lower_ci, upper_ci, Has_Survey_Data
    ))
  }

  # ==================================================================
  # PREDICTION FOR COUNTRIES WITHOUT SURVEY DATA
  # These countries sample from the hierarchical prior
  # ==================================================================

  nodata_countries <- all_country_lookup %>%
    filter(!has_data) %>%
    filter(!is.na(Region_Num))  # Only include countries with valid region mapping

  if (nrow(nodata_countries) > 0) {

    cat(sprintf("\n  Generating predictions for %d countries WITHOUT survey data...\n",
                nrow(nodata_countries)))

    # Create prediction grid for no-data countries (no Num_Country needed)
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

    # Add spline bases
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

    # Age-cohort interactions
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

    # CRITICAL: Use the SAME sampled indices as data countries for consistency
    # This ensures all countries (data and no-data) use the same MCMC iterations
    mcmc_samples_subsetted <- combined_samples_matrix[sampled_indices, ]

    # Export additional objects to cluster
    clusterExport(cl, c(
      "all_country_lookup", "nodata_grid_full", "country_name_mapping",
      "sample_cig_params_from_prior", "sample_extra_intercept_from_prior",
      "cig_samples", "smkextra_samples", "anyextra_samples",
      "def_code_shared_samples", "mcmc_samples_subsetted", "n_samples",
      "nimble_constants", "gender_dir"
    ), envir = environment())

    nodata_results <- foreach(
      country_code = nodata_countries$wb_country_abv,
      .packages = c("dplyr", "splines"),
      .combine = rbind,
      .errorhandling = 'remove'
    ) %dopar% {

      # Get region for this country
      country_info <- all_country_lookup[all_country_lookup$wb_country_abv == country_code, ]
      country_region <- country_info$Region_Num

      # Safety check: skip if region not found (shouldn't happen due to filter, but defensive)
      if (is.na(country_region) || length(country_region) == 0) {
        warning(sprintf("Skipping %s: no valid region mapping", country_code))
        return(NULL)
      }

      country_full_name <- country_name_mapping[country_code]

      if (is.na(country_full_name) || is.null(country_full_name)) {
        country_full_name <- toupper(country_code)
      }

      country_dir <- file.path(gender_dir, country_full_name)
      dir.create(country_dir, showWarnings = FALSE, recursive = TRUE)

      # Save metadata indicating this country has no survey data
      saveRDS(
        list(has_survey_data = FALSE, region = country_info$region_consolidated),
        file = file.path(country_dir, "_metadata.rds")
      )

      # Sample country parameters from hierarchical prior
      # IMPORTANT: Use mcmc_samples_subsetted (same indices as data countries)
      # CIG: Full country-specific parameters
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

      # SMKEXTRA/ANYEXTRA: Only intercepts
      smkextra_country_intercept <- sample_extra_intercept_from_prior("smkextra", mcmc_samples_subsetted, n_samples)
      anyextra_country_intercept <- sample_extra_intercept_from_prior("anyextra", mcmc_samples_subsetted, n_samples)

      # Prepare data with Def_Code column
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

          # CIG: Use sampled country parameters from hierarchical prior
          mu_cig <- cig_samples$global_intercept +
            cig_samples$region_intercept[, paste0("cig_region_intercept[", country_region, "]")] +
            cig_country_intercept +
            def_code_shared_samples * def_code_binary +
            spline_weight_value * as.matrix(cig_age_spline_sampled) %*% age_spline_values +
            linear_weight_value * cig_samples$age_linear_smooth * age_linear_smooth +
            as.matrix(cig_cohort_spline_sampled) %*% cohort_spline_values +
            as.matrix(cig_samples$age_cohort) %*% age_cohort_interaction_values

          # SMKEXTRA: Country intercept + regional splines
          mu_smkextra <- smkextra_samples$global_intercept +
            smkextra_samples$region_intercept[, paste0("smkextra_region_intercept[", country_region, "]")] +
            smkextra_country_intercept +
            0.3 * def_code_shared_samples * def_code_binary +
            spline_weight_value * as.matrix(smkextra_samples$age_spline_regional[, paste0("smkextra_age_spline_region_mean[", country_region, ", ", 1:length(age_spline_values), "]")]) %*% age_spline_values +
            linear_weight_value * smkextra_samples$age_linear_smooth * age_linear_smooth +
            as.matrix(smkextra_samples$cohort_spline_regional[, paste0("smkextra_cohort_spline_region_mean[", country_region, ", ", 1:length(cohort_spline_values), "]")]) %*% cohort_spline_values

          # ANYEXTRA: Country intercept + regional splines
          mu_anyextra <- anyextra_samples$global_intercept +
            anyextra_samples$region_intercept[, paste0("anyextra_region_intercept[", country_region, "]")] +
            anyextra_country_intercept +
            0.3 * def_code_shared_samples * def_code_binary +
            spline_weight_value * as.matrix(anyextra_samples$age_spline_regional[, paste0("anyextra_age_spline_region_mean[", country_region, ", ", 1:length(age_spline_values), "]")]) %*% age_spline_values +
            linear_weight_value * anyextra_samples$age_linear_smooth * age_linear_smooth +
            as.matrix(anyextra_samples$cohort_spline_regional[, paste0("anyextra_cohort_spline_region_mean[", country_region, ", ", 1:length(cohort_spline_values), "]")]) %*% cohort_spline_values

          predictions_cig[j, ] <- mu_cig

          # STICK-BREAKING
          p_cig <- plogis(mu_cig)
          p_smoked_full <- p_cig + plogis(mu_smkextra) * (1 - p_cig)
          predictions_smoked[j, ] <- log(p_smoked_full / (1 - p_smoked_full))

          p_any_full <- p_smoked_full + plogis(mu_anyextra) * (1 - p_smoked_full)
          predictions_any[j, ] <- log(p_any_full / (1 - p_any_full))
        }

        # Check ordering constraint
        has_violations <- check_ordering_per_draw(predictions_cig, predictions_smoked, predictions_any)
        if (has_violations) {
          warning(sprintf("Ordering violations for %s - %s - %d", country_code, current_def, current_year))
        }

        # Save predictions
        saveRDS(predictions_cig, file = file.path(cig_dir, paste0(current_year, ".rds")), compress = "xz")
        saveRDS(predictions_smoked, file = file.path(smoked_dir, paste0(current_year, ".rds")), compress = "xz")
        saveRDS(predictions_any, file = file.path(any_dir, paste0(current_year, ".rds")), compress = "xz")

        current_data_cig <- current_data %>%
          mutate(
            Prevalence = apply(predictions_cig, 1, function(x) mean(plogis(x))),
            lower_ci = apply(predictions_cig, 1, function(x) quantile(plogis(x), 0.025)),
            upper_ci = apply(predictions_cig, 1, function(x) quantile(plogis(x), 0.975)),
            Def_Type_Code = paste0(current_def, "_cigarettes")
          )

        current_data_smoked <- current_data %>%
          mutate(
            Prevalence = apply(predictions_smoked, 1, function(x) mean(plogis(x))),
            lower_ci = apply(predictions_smoked, 1, function(x) quantile(plogis(x), 0.025)),
            upper_ci = apply(predictions_smoked, 1, function(x) quantile(plogis(x), 0.975)),
            Def_Type_Code = paste0(current_def, "_any_smoked_tobacco")
          )

        current_data_any <- current_data %>%
          mutate(
            Prevalence = apply(predictions_any, 1, function(x) mean(plogis(x))),
            lower_ci = apply(predictions_any, 1, function(x) quantile(plogis(x), 0.025)),
            upper_ci = apply(predictions_any, 1, function(x) quantile(plogis(x), 0.975)),
            Def_Type_Code = paste0(current_def, "_any_tobacco_product")
          )

        results_list[[length(results_list) + 1]] <- current_data_cig
        results_list[[length(results_list) + 1]] <- current_data_smoked
        results_list[[length(results_list) + 1]] <- current_data_any

        rm(predictions_cig, predictions_smoked, predictions_any)
        gc()  # Memory cleanup after each year-def combo
      }

      result <- do.call(rbind, results_list)
      result$Country <- country_code
      result$Sex <- gender
      result$Has_Survey_Data <- FALSE

      return(result %>% select(
        Year, Age_Midpoint, Birth_Cohort, Def_Type_Code,
        Sex, Country, Prevalence, lower_ci, upper_ci, Has_Survey_Data
      ))
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
      if (!is.null(nodata_results) && nrow(nodata_results) > 0) {
        country_results <- bind_rows(country_results, nodata_results)
        cat(sprintf("  Added predictions for %d countries without survey data\n",
                    length(unique(nodata_results$Country))))
      }
    }
  }

  stopCluster(cl)
  gender_results[[gender]] <- country_results
  
  # Clean up NIMBLE objects
  try(nimble::clearCompiled(nimble_model), silent = TRUE)
  rm(nimble_model, compiled_model, compiled_mcmc, mcmc_built, samples)
  gc()
  
  rm(country_results, combined_samples_matrix)
  gc()
}

final_ac_predictions <- do.call(rbind, gender_results)
write.csv(final_ac_predictions, file = "results/final_ac_predictions_nested.csv", row.names = FALSE)

# Summary of all-countries predictions
cat("\n  === PREDICTION SUMMARY ===\n")
if ("Has_Survey_Data" %in% names(final_ac_predictions)) {
  n_data_countries <- length(unique(final_ac_predictions$Country[final_ac_predictions$Has_Survey_Data]))
  n_nodata_countries <- length(unique(final_ac_predictions$Country[!final_ac_predictions$Has_Survey_Data]))
  cat(sprintf("  Total countries predicted: %d\n", n_data_countries + n_nodata_countries))
  cat(sprintf("    - Countries WITH survey data: %d\n", n_data_countries))
  cat(sprintf("    - Countries WITHOUT survey data: %d (sampled from hierarchical prior)\n", n_nodata_countries))
} else {
  cat(sprintf("  Total countries predicted: %d\n", length(unique(final_ac_predictions$Country))))
}

cat("\n  Global model fitting complete (NIMBLE)\n")


#########################################################################################
#                              SECTION 9: TARGET PREVALENCE CALCULATION                 
#########################################################################################

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = RANDOM_SEED)
clusterExport(cl, "calculate_weighted_prevalence")

target_prevalence_list <- list()

for (gender in genders) {
  cat(sprintf("\nCalculating targets: %s\n", gender))
  
  # ---- GLOBAL MODEL TARGETS ----
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
      
      # 4% Absolute Target
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
  
  # ---- COUNTRY-SPECIFIC MODEL TARGETS ----
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
      
      # 4% Absolute Target
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
  
  target_prevalence_list[[paste0(gender, "_country")]] <- country_results
}

stopCluster(cl)

target_prevalence_df <- do.call(rbind, target_prevalence_list)
rownames(target_prevalence_df) <- NULL

saveRDS(target_prevalence_df, file = "results/target_prevalences_dual_evaluation.rds", compress = "xz")
write.csv(target_prevalence_df, file = "results/target_prevalences_dual_evaluation.csv", row.names = FALSE)

cat("Target prevalence calculation complete (both models)\n")


#########################################################################################
#                              SECTION 10: WEIGHTED PREVALENCE TRENDS (GLOBAL MODEL)    
#########################################################################################

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = RANDOM_SEED)
clusterExport(cl, "calculate_weighted_prevalence")

gender_weighted_results_global <- list()

for (gender in genders) {
  cat(sprintf("Processing weighted calculations (Global model): %s\n", gender))
  
  countries <- list.dirs(file.path("processing", gender), full.names = FALSE, recursive = FALSE)
  
  weighted_results <- foreach(
    country = countries,
    .packages = c("dplyr", "tidyr"),
    .combine = rbind,
    .errorhandling = 'remove'
  ) %dopar% {
    
    def_types <- list.dirs(file.path("processing", gender, country),
                           full.names = FALSE, recursive = FALSE)
    
    country_results <- list()
    
    country_code <- names(country_name_mapping)[country_name_mapping == country]
    if (length(country_code) == 0) country_code <- country
    
    for (def_type in def_types) {
      year_files <- list.files(file.path("processing", gender, country, def_type),
                               pattern = "\\.rds$")
      years <- as.numeric(sub("\\.rds$", "", year_files))
      
      target_records <- target_prevalence_df %>%
        filter(Sex == gender, Country == country_code, Def_Type_Code == def_type)
      
      if (nrow(target_records) == 0) next
      
      for (year_filtered in years) {
        year_weights <- weights_cleaned %>%
          filter(area == country_code, sex == gender, year == year_filtered) %>%
          arrange(age)
        
        if (nrow(year_weights) == 0) next
        
        iterations_matrix_logit <- readRDS(file.path("processing", gender, country, def_type,
                                                     paste0(year_filtered, ".rds")))
        
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
  
  gender_weighted_results_global[[gender]] <- weighted_results
}

stopCluster(cl)

final_weighted_results_global <- do.call(rbind, gender_weighted_results_global)
rownames(final_weighted_results_global) <- NULL

write.csv(
  final_weighted_results_global,
  file = "results/final_weighted_results_global_model.csv",
  row.names = FALSE
)

cat("Global model weighted calculations complete\n")


#########################################################################################
#                              SECTION 11: COUNTRY-SPECIFIC MODELS (NIMBLE)             
#########################################################################################

cat("\n")
cat("================================================================\n")
cat("  FITTING COUNTRY-SPECIFIC MODELS (NIMBLE)\n")
cat("================================================================\n")

if (!dir.exists("results/country_specific_ac_nested")) {
  dir.create("results/country_specific_ac_nested", recursive = TRUE)
}

country_specific_ac_predictions <- list()

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

  for (country_code in countries) {
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
      list(mean = row$mean, sd = max(row$sd, 1e-6))
    }
    
    # Helper function to convert SD to precision safely
    sd_to_prec <- function(sd_val) {
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
    
    # Extract prior vectors with proper error handling and CORRECT ORDERING
    # FIX: Use explicit index lookup to ensure parameters are in correct order (1, 2, 3, 4, ...)
    # instead of relying on grep which may return rows in CSV file order

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
    cig_age_cohort_priors <- get_prior_vector("cig_age_cohort_interaction_", nAgeXCohortSplines)

    smkextra_age_spline_priors <- get_prior_vector("smkextra_age_spline_", nAgeSpline)
    smkextra_cohort_spline_priors <- get_prior_vector("smkextra_cohort_spline_", nCohortSpline)
    # REMOVED: smkextra_age_cohort (not in model)

    anyextra_age_spline_priors <- get_prior_vector("anyextra_age_spline_", nAgeSpline)
    anyextra_cohort_spline_priors <- get_prior_vector("anyextra_cohort_spline_", nCohortSpline)
    # REMOVED: anyextra_age_cohort (not in model)
    
    # CORRECTED: No smkextra/anyextra age-cohort constants
    nimble_constants <- list(
      N = nrow(country_data),
      nSurvey = nSurvey,
      nAgeSpline = nAgeSpline,
      nCohortSpline = nCohortSpline,
      nAgeXCohortSplines = nAgeXCohortSplines,
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

      cig_age_cohort_prior_means = cig_age_cohort_priors$means,
      cig_age_cohort_prior_precs = sd_to_prec(cig_age_cohort_priors$sds),
      
      # SMKEXTRA priors (NO age-cohort)
      smkextra_intercept_prior_mean = get_prior("smkextra_intercept")$mean,
      smkextra_intercept_prior_prec = sd_to_prec(get_prior("smkextra_intercept")$sd),
      
      smkextra_age_linear_smooth_effect_prior_mean = get_prior("smkextra_age_linear_smooth_effect")$mean,
      smkextra_age_linear_smooth_effect_prior_prec = sd_to_prec(get_prior("smkextra_age_linear_smooth_effect")$sd),
      
      # SMKEXTRA spline priors (using ordered extraction)
      smkextra_age_spline_prior_means = smkextra_age_spline_priors$means,
      smkextra_age_spline_prior_precs = sd_to_prec(smkextra_age_spline_priors$sds),

      smkextra_cohort_spline_prior_means = smkextra_cohort_spline_priors$means,
      smkextra_cohort_spline_prior_precs = sd_to_prec(smkextra_cohort_spline_priors$sds),
      # REMOVED: smkextra_age_cohort_* (not in model)
      
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
      # REMOVED: anyextra_age_cohort_* (not in model)
    )
    
    # DATA
    nimble_data <- list(
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
      weight = country_data$weight
    )
    
    # ---- Generate Initial Values (CORRECTED - no smkextra/anyextra age-cohort) ----
    
    generate_country_inits <- function(seed_offset = 0) {
      set.seed(RANDOM_SEED + seed_offset)

      # FIX: Use ordered prior means from get_prior_vector() results
      list(
        # Scalar parameters
        cig_def_code_shared = rnorm(1, get_prior("cig_def_code_shared")$mean, 0.1),
        cig_intercept = rnorm(1, get_prior("cig_intercept")$mean, 0.2),
        smkextra_intercept = rnorm(1, get_prior("smkextra_intercept")$mean, 0.2),
        anyextra_intercept = rnorm(1, get_prior("anyextra_intercept")$mean, 0.2),

        # Age linear effects (negative, constrained)
        cig_age_linear_smooth_effect = runif(1, -0.05, -0.01),
        smkextra_age_linear_smooth_effect = runif(1, -0.05, -0.01),
        anyextra_age_linear_smooth_effect = runif(1, -0.05, -0.01),

        # Residual SD
        residual_sd = runif(1, 0.5, 1.0),

        # Survey effects
        survey_intercept_precision = rgamma(1, 3, 1),
        survey_intercept = rnorm(nSurvey, 0, 0.1),

        # CIG spline arrays (using ordered prior means)
        cig_age_spline = rnorm(nAgeSpline, cig_age_spline_priors$means, 0.05),
        cig_cohort_spline = rnorm(nCohortSpline, cig_cohort_spline_priors$means, 0.05),
        cig_age_cohort_interaction = rnorm(nAgeXCohortSplines, cig_age_cohort_priors$means, 0.02),

        # SMKEXTRA spline arrays (NO age-cohort)
        smkextra_age_spline = rnorm(nAgeSpline, smkextra_age_spline_priors$means, 0.05),
        smkextra_cohort_spline = rnorm(nCohortSpline, smkextra_cohort_spline_priors$means, 0.05),

        # ANYEXTRA spline arrays (NO age-cohort)
        anyextra_age_spline = rnorm(nAgeSpline, anyextra_age_spline_priors$means, 0.05),
        anyextra_cohort_spline = rnorm(nCohortSpline, anyextra_cohort_spline_priors$means, 0.05)
      )
    }
    
    inits_list <- lapply(1:NUMBER_OF_CHAINS, function(i) generate_country_inits(i))
    
    # ---- Build, Configure, Compile, Run NIMBLE ----
    
    tryCatch({
      nimble_model <- nimbleModel(
        code = regional_country_specific_ac_model_nimble,
        constants = nimble_constants,
        data = nimble_data,
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
        monitors = c(
          "cig_def_code_shared",
          "cig_intercept", "cig_age_spline", "cig_age_linear_smooth_effect",
          "cig_cohort_spline", "cig_age_cohort_interaction",
          "smkextra_intercept", "smkextra_age_spline", "smkextra_age_linear_smooth_effect",
          "smkextra_cohort_spline",
          "anyextra_intercept", "anyextra_age_spline", "anyextra_age_linear_smooth_effect",
          "anyextra_cohort_spline",
          "residual_sd"
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
      
      saveRDS(combined_samples, file = file.path(country_dir, "posterior_samples.rds"), compress = "xz")
      
      # Extract samples for prediction
      def_code_shared_samples <- combined_samples[, "cig_def_code_shared"]
      
      # CORRECTED: No age-cohort for smkextra/anyextra
      extract_head_samples_country <- function(prefix, combined_samples, has_age_cohort = FALSE) {
        result <- list(
          intercept     = combined_samples[, paste0(prefix, "_intercept")],
          age_linear    = combined_samples[, paste0(prefix, "_age_linear_smooth_effect")],
          age_spline    = combined_samples[, grep(paste0("^", prefix, "_age_spline\\["), colnames(combined_samples)), drop = FALSE],
          cohort_spline = combined_samples[, grep(paste0("^", prefix, "_cohort_spline\\["), colnames(combined_samples)), drop = FALSE]
        )
        if (has_age_cohort) {
          result$age_cohort <- combined_samples[, grep(paste0("^", prefix, "_age_cohort_interaction\\["), colnames(combined_samples)), drop = FALSE]
        }
        result
      }
      
      cig_samples_full      <- extract_head_samples_country("cig", combined_samples, has_age_cohort = TRUE)
      smkextra_samples_full <- extract_head_samples_country("smkextra", combined_samples, has_age_cohort = FALSE)
      anyextra_samples_full <- extract_head_samples_country("anyextra", combined_samples, has_age_cohort = FALSE)
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
      cig_ac_s           <- cig_samples_full$age_cohort[sampled_indices, , drop = FALSE]

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
      # [OPT-6] Precompute year-specific cohort splines and interactions
      # ============================================================
      years_pred <- seq(min(as.numeric(country_data$year)), target_year + PROJECTED_YEARS, by = 1)
      n_years <- length(years_pred)

      # Compute ALL birth cohorts at once: flatten (n_years * n_ages) vector
      all_birth_cohorts <- rep(years_pred, each = n_ages) - rep(age_midpoints_shared, times = n_years)

      # Single call to ns() for all cohort values
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
          saveRDS(predictions_cig, file = file.path(cig_dir, paste0(current_year, ".rds")), compress = "xz")
          saveRDS(predictions_smoked, file = file.path(smoked_dir, paste0(current_year, ".rds")), compress = "xz")
          saveRDS(predictions_any, file = file.path(any_dir, paste0(current_year, ".rds")), compress = "xz")

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
      
      country_specific_ac_predictions[[paste(country_code, gender, sep = "_")]] <- prediction_results
      
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
      rm(compiled_model, compiled_mcmc, mcmc_built, samples, results_list,
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
  }
}

final_predictions_country_specific_ac <- do.call(rbind, country_specific_ac_predictions)

# Handle case where all country models failed (empty list returns NULL)
if (is.null(final_predictions_country_specific_ac) || nrow(final_predictions_country_specific_ac) == 0) {
  warning("No country-specific predictions were generated. All models may have failed.")
  final_predictions_country_specific_ac <- data.frame(
    Year = integer(), Age_Midpoint = numeric(), Birth_Cohort = numeric(),
    Def_Type_Code = character(), Sex = character(), Country = character(),
    Prevalence = numeric(), lower_ci = numeric(), upper_ci = numeric(),
    stringsAsFactors = FALSE
  )
} else {
  rownames(final_predictions_country_specific_ac) <- NULL
}

write.csv(
  final_predictions_country_specific_ac,
  file = "results/final_predictions_country_specific.csv",
  row.names = FALSE
)

cat("\nCountry-specific model fitting complete (NIMBLE)\n")

#########################################################################################
#                 SECTION 12: WEIGHTED PREVALENCE TRENDS (COUNTRY-SPECIFIC MODEL)       #
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


#########################################################################################
#                              SECTION 13: MODEL EVALUATION (RMSE)                      #
#########################################################################################

cat("\n")
cat("================================================================\n")
cat("  MODEL EVALUATION: RMSE-BASED SELECTION (BY COUNTRY AND SEX)\n")
cat("================================================================\n")

dir.create("evaluation", showWarnings = FALSE)

# ---- 13.1 Prepare Observed Data ----

observed_data <- clean_data %>%
  mutate(
    Year       = as.numeric(year),
    Age_Midpoint = (start_age + end_age) / 2,
    Prevalence = plogis(prevalence),
    Data_Type  = "Observed"
  ) %>%
  select(Year, Age_Midpoint, def_type_code, sex, wb_country_abv,
         Prevalence, Data_Type) %>%
  rename(Indicator = def_type_code, Sex = sex, Country = wb_country_abv)

# ---- 13.2 Prepare Model Predictions ----

global_ac_predictions <- final_ac_predictions %>%
  mutate(Year = as.numeric(Year), Data_Type = "Global_Model") %>%
  select(Year, Age_Midpoint, Def_Type_Code, Sex, Country,
         Prevalence, Data_Type) %>%
  rename(Indicator = Def_Type_Code)

country_ac_predictions <- final_predictions_country_specific_ac %>%
  mutate(Year = as.numeric(Year), Data_Type = "Country_Model") %>%
  select(Year, Age_Midpoint, Def_Type_Code, Sex, Country,
         Prevalence, Data_Type) %>%
  rename(Indicator = Def_Type_Code)

all_ac_predictions <- bind_rows(global_ac_predictions, country_ac_predictions)

# ---- 13.3 Calculate RMSE ----

calculate_rmse <- function(country_code, gender, indicator) {
  obs_data <- observed_data %>%
    filter(Country == country_code, Sex == gender, Indicator == indicator)
  
  if (nrow(obs_data) == 0) {
    return(data.frame(
      Country        = country_code,
      Sex            = gender,
      Indicator      = indicator,
      Global_RMSE    = NA_real_,
      Country_RMSE   = NA_real_,
      N_Observations = 0,
      N_Matched      = 0,
      Selected_Model = NA_character_
    ))
  }
  
  model_data <- all_ac_predictions %>%
    filter(Country == country_code, Sex == gender, Indicator == indicator)
  
  global_matches <- 0
  global_sse     <- 0
  country_matches <- 0
  country_sse     <- 0
  
  age_tolerance <- 2.5
  
  for (i in 1:nrow(obs_data)) {
    curr_obs <- obs_data[i, ]
    
    global_pred <- model_data %>%
      filter(Data_Type == "Global_Model", Year == curr_obs$Year) %>%
      mutate(age_diff = abs(Age_Midpoint - curr_obs$Age_Midpoint)) %>%
      filter(age_diff <= age_tolerance) %>%
      slice_min(age_diff, n = 1)
    
    country_pred <- model_data %>%
      filter(Data_Type == "Country_Model", Year == curr_obs$Year) %>%
      mutate(age_diff = abs(Age_Midpoint - curr_obs$Age_Midpoint)) %>%
      filter(age_diff <= age_tolerance) %>%
      slice_min(age_diff, n = 1)
    
    if (nrow(global_pred) > 0) {
      global_matches <- global_matches + 1
      global_sse     <- global_sse + (global_pred$Prevalence[1] - curr_obs$Prevalence)^2
    }
    
    if (nrow(country_pred) > 0) {
      country_matches <- country_matches + 1
      country_sse     <- country_sse + (country_pred$Prevalence[1] - curr_obs$Prevalence)^2
    }
  }
  
  global_rmse  <- if (global_matches > 0) sqrt(global_sse / global_matches) else NA_real_
  country_rmse <- if (country_matches > 0) sqrt(country_sse / country_matches) else NA_real_
  
  selected_model <- case_when(
    is.na(global_rmse) & is.na(country_rmse) ~ NA_character_,
    is.na(country_rmse)                      ~ "Global",
    is.na(global_rmse)                       ~ "Country",
    global_rmse <= country_rmse              ~ "Global",
    TRUE                                     ~ "Country"
  )
  
  return(data.frame(
    Country        = country_code,
    Sex            = gender,
    Indicator      = indicator,
    Global_RMSE    = global_rmse,
    Country_RMSE   = country_rmse,
    N_Observations = nrow(obs_data),
    N_Matched      = global_matches,
    Selected_Model = selected_model
  ))
}

unique_countries  <- unique(c(observed_data$Country, all_ac_predictions$Country))
unique_genders    <- c("males", "females")
unique_indicators <- unique(c(observed_data$Indicator, all_ac_predictions$Indicator))

rmse_results <- list()
for (country in unique_countries) {
  for (gender in unique_genders) {
    for (indicator in unique_indicators) {
      result <- calculate_rmse(country, gender, indicator)
      rmse_results[[length(rmse_results) + 1]] <- result
    }
  }
}

rmse_df <- bind_rows(rmse_results)

# ---- 13.4 Create Country-Sex Summary ----

country_sex_summary <- rmse_df %>%
  filter(N_Observations > 0, !is.na(Selected_Model)) %>%
  group_by(Country, Sex) %>%
  summarize(
    Total_Observations    = sum(N_Observations, na.rm = TRUE),
    N_Indicators_Compared = n(),
    Global_Better_Count   = sum(Selected_Model == "Global", na.rm = TRUE),
    Country_Better_Count  = sum(Selected_Model == "Country", na.rm = TRUE),
    Mean_Global_RMSE      = mean(Global_RMSE, na.rm = TRUE),
    Mean_Country_RMSE     = mean(Country_RMSE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Final_Selected_Model = case_when(
      N_Indicators_Compared == 0             ~ "Global",
      Country_Better_Count > Global_Better_Count ~ "Country",
      TRUE                                   ~ "Global"
    )
  )

# Handle country-sex combinations with no data
all_country_sex <- expand.grid(
  Country = unique_countries,
  Sex     = unique_genders,
  stringsAsFactors = FALSE
)

country_sex_summary <- all_country_sex %>%
  left_join(country_sex_summary, by = c("Country", "Sex")) %>%
  mutate(
    Total_Observations    = replace_na(Total_Observations, 0),
    N_Indicators_Compared = replace_na(N_Indicators_Compared, 0),
    Global_Better_Count   = replace_na(Global_Better_Count, 0),
    Country_Better_Count  = replace_na(Country_Better_Count, 0),
    Final_Selected_Model  = replace_na(Final_Selected_Model, "Global")
  )

write.csv(rmse_df, "evaluation/model_evaluation_detailed.csv", row.names = FALSE)
write.csv(country_sex_summary, "evaluation/model_evaluation_by_country_sex.csv", row.names = FALSE)

# Summary statistics
selection_summary <- country_sex_summary %>%
  group_by(Sex, Final_Selected_Model) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Final_Selected_Model, values_from = n, values_fill = 0)

cat("\nModel Selection Results by Sex:\n")
print(selection_summary)

# Create legacy country_summary
country_summary <- country_sex_summary %>%
  filter(N_Indicators_Compared > 0) %>%
  group_by(Country) %>%
  summarize(
    Total_Observations   = sum(Total_Observations, na.rm = TRUE),
    Global_Better_Count  = sum(Final_Selected_Model == "Global", na.rm = TRUE),
    Country_Better_Count = sum(Final_Selected_Model == "Country", na.rm = TRUE),
    Mean_Global_RMSE     = mean(Mean_Global_RMSE, na.rm = TRUE),
    Mean_Country_RMSE    = mean(Mean_Country_RMSE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Final_Selected_Model = case_when(
      Country_Better_Count > Global_Better_Count ~ "Country",
      TRUE                                       ~ "Global"
    )
  )

# Add countries with no data
countries_with_no_data <- setdiff(unique_countries, country_summary$Country)
if (length(countries_with_no_data) > 0) {
  no_data_rows <- data.frame(
    Country              = countries_with_no_data,
    Total_Observations   = 0,
    Global_Better_Count  = 0,
    Country_Better_Count = 0,
    Mean_Global_RMSE     = NA_real_,
    Mean_Country_RMSE    = NA_real_,
    Final_Selected_Model = "Global"
  )
  country_summary <- bind_rows(country_summary, no_data_rows)
}

write.csv(country_summary, "evaluation/model_evaluation_summary.csv", row.names = FALSE)

cat(sprintf("\n  Country-Sex combinations using global model: %d\n",
            sum(country_sex_summary$Final_Selected_Model == "Global")))
cat(sprintf("  Country-Sex combinations using country model: %d\n",
            sum(country_sex_summary$Final_Selected_Model == "Country")))


#########################################################################################
#                  SECTION 14: CREATE POST-SELECTION COMBINED RESULTS                   #
#########################################################################################

cat("\n  Creating post-selection combined weighted results...\n")

# Combine global and country weighted results based on model selection
# For each Country-Sex combination, select results from the chosen model

# Add Model_Type to distinguish sources
final_weighted_results_global <- final_weighted_results_global %>%
  mutate(Model_Type = "Global")

final_weighted_results_country <- final_weighted_results_country %>%
  mutate(Model_Type = "Country")

# Combine both datasets
all_weighted_results <- bind_rows(
  final_weighted_results_global,
  final_weighted_results_country
)

# Join with model selection to filter to selected model only
final_weighted_results_selected <- all_weighted_results %>%
  left_join(
    country_sex_summary %>% select(Country, Sex, Final_Selected_Model),
    by = c("Country", "Sex")
  ) %>%
  # Keep only rows where Model_Type matches the selected model
  filter(Model_Type == Final_Selected_Model) %>%
  select(-Final_Selected_Model)  # Remove helper column

# Write to CSV for later use
write.csv(
  final_weighted_results_selected,
  file = "results/final_weighted_results_apc_post_selection.csv",
  row.names = FALSE
)

cat(sprintf("  Post-selection results: %d rows\n", nrow(final_weighted_results_selected)))
cat("  Saved to: results/final_weighted_results_apc_post_selection.csv\n")


#########################################################################################
#                                                                                       #
#              SECTION 15: POST-PROCESSING & PUBLICATION PIPELINE                       #
#                                                                                       #
#########################################################################################

cat("\n")
cat("================================================================\n")
cat("  STARTING POST-PROCESSING & VISUALIZATION PIPELINE\n")
cat("================================================================\n")

# ---- 15.1 Load Model Selection Data ----
# Ensure we have the selection summary from Section 12
if (!exists("country_sex_summary")) {
  country_sex_summary <- read.csv("evaluation/model_evaluation_by_country_sex.csv")
}

# ---- 15.2 Prepare Regional Population Weights ----
cat("  Preparing regional population weights...\n")

regional_population_weights <- weights_cleaned %>%
  left_join(country_region_mapping %>% select(wb_country_abv, region_consolidated), 
            by = c("area" = "wb_country_abv")) %>%
  mutate(region_consolidated = ifelse(is.na(region_consolidated), "Other", region_consolidated)) %>%
  group_by(region_consolidated, year, sex, age) %>%
  summarise(
    regional_population = sum(weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(region_consolidated, year, sex) %>%
  mutate(regional_age_weight = regional_population / sum(regional_population, na.rm = TRUE)) %>%
  ungroup()

country_regional_weights <- weights_cleaned %>%
  left_join(country_region_mapping %>% select(wb_country_abv, region_consolidated), 
            by = c("area" = "wb_country_abv")) %>%
  mutate(region_consolidated = ifelse(is.na(region_consolidated), "Other", region_consolidated)) %>%
  group_by(region_consolidated, year, sex, age) %>%
  mutate(
    regional_total_pop = sum(weight, na.rm = TRUE),
    country_regional_weight = weight / regional_total_pop
  ) %>%
  ungroup()

#########################################################################################
#                          REGIONAL AGGREGATION (SAMPLE-BASED)                          #
#########################################################################################

# ---- 15.3 Define Aggregation Function ----
aggregate_selected_country_samples_to_region <- function(region_name, gender, year_range) {
  
  # 1. Identify countries in region
  countries_in_region <- country_region_mapping %>%
    filter(region_consolidated == region_name) %>%
    pull(wb_country_abv)
  
  if(length(countries_in_region) == 0) return(NULL)
  
  # 2. Get Model Selection for this Region/Sex
  selection_subset <- country_sex_summary %>%
    filter(Country %in% countries_in_region, Sex == gender) %>%
    select(Country, Final_Selected_Model)
  
  # 3. Identify available definition types from one country (assuming consistency)
  # We check the processing folder as a fallback source of truth for folder structure
  example_country <- countries_in_region[1]
  example_path <- file.path("processing", gender, country_name_mapping[example_country])
  if(!dir.exists(example_path)) return(NULL)
  
  def_types <- list.dirs(example_path, full.names = FALSE, recursive = FALSE)
  def_types <- def_types[def_types != ""]

  for (def_type in def_types) {
    
    # Create output directory for this region/def_type
    regional_out_dir <- file.path("results/regional_aggregation", "selected_models", gender, region_name, def_type)
    dir.create(regional_out_dir, recursive = TRUE, showWarnings = FALSE)
    
    for (year in year_range) {
      
      # Prepare containers
      country_samples_matrix_list <- list()
      country_weights_vector_list <- list()
      
      for (country_code in countries_in_region) {
        country_full_name <- country_name_mapping[country_code]
        
        # Determine Model Path based on Selection
        model_choice <- selection_subset$Final_Selected_Model[selection_subset$Country == country_code]
        if(length(model_choice) == 0) model_choice <- "Global" # Default
        
        if (model_choice == "Country") {
          sample_path <- file.path("results/country_specific_ac_nested", gender, country_full_name, def_type, paste0(year, ".rds"))
        } else {
          sample_path <- file.path("processing", gender, country_full_name, def_type, paste0(year, ".rds"))
        }
        
        if (!file.exists(sample_path)) next
        
        # Load Samples
        s_mat <- readRDS(sample_path)
        
        # Get Weights
        w_vec <- country_regional_weights %>%
          filter(area == country_code, year == !!year, sex == gender) %>%
          arrange(age) %>%
          pull(country_regional_weight)
        
        if (nrow(s_mat) == length(w_vec)) {
          country_samples_matrix_list[[country_code]] <- s_mat
          country_weights_vector_list[[country_code]] <- w_vec
        }
      }
      
      if (length(country_samples_matrix_list) == 0) next
      
      # Perform Aggregation (Random Draw Matching)
      # We assume 1000 samples. If mismatch, we sample indices.
      n_samples_out <- 1000
      n_ages <- nrow(country_samples_matrix_list[[1]])
      regional_mat <- matrix(0, nrow = n_ages, ncol = n_samples_out)
      
      for (s in 1:n_samples_out) {
        # Create a combined weighted prevalence for draw 's'
        numerator <- rep(0, n_ages)
        denominator <- rep(0, n_ages)
        
        for (cc in names(country_samples_matrix_list)) {
          mat <- country_samples_matrix_list[[cc]]
          w <- country_weights_vector_list[[cc]]
          
          # Handle different sample counts by random sampling index
          idx <- sample(ncol(mat), 1)
          
          # Convert logit to prob, weight it
          prob <- plogis(mat[, idx])
          numerator <- numerator + (prob * w)
          denominator <- denominator + w
        }
        
        # Normalize and convert back to logit
        regional_prob <- numerator / denominator
        regional_prob <- pmin(pmax(regional_prob, 0.0001), 0.9999) # Numerical stability
        regional_mat[, s] <- qlogis(regional_prob)
      }
      
      saveRDS(regional_mat, file = file.path(regional_out_dir, paste0(year, ".rds")))
    }
  }
  return(TRUE)
}

# ---- 15.4 Run Aggregation (Parallel) ----
cat("  Running regional aggregation (Parallel)...\n")

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
clusterExport(cl, c("country_region_mapping", "country_sex_summary", "country_regional_weights", 
                    "country_name_mapping", "aggregate_selected_country_samples_to_region", "plogis", "qlogis"))

invisible(foreach(r = unique(country_region_mapping$region_consolidated), .packages = c("dplyr")) %:%
            foreach(g = c("males", "females"), .packages = c("dplyr")) %dopar% {
              aggregate_selected_country_samples_to_region(r, g, seq(BASE_YEAR, TARGET_YEAR + PROJECTED_YEARS))
            })
stopCluster(cl)
cat("  Regional aggregation complete.\n")

#########################################################################################
#                       CALCULATE REGIONAL METRICS & ACHIEVEMENT                        #
#########################################################################################

cat("  Calculating regional targets and high-certainty achievement rates...\n")

# Helper to calculate metrics from saved regional samples
calc_regional_metrics <- function() {
  results_list <- list()
  
  regions <- list.dirs("results/regional_aggregation/selected_models", recursive = FALSE) 
  # Note: structure is selected_models/gender/region/def_type/year.rds
  
  for (gender in c("males", "females")) {
    gender_path <- file.path("results/regional_aggregation/selected_models", gender)
    if (!dir.exists(gender_path)) next
    
    region_names <- list.dirs(gender_path, full.names = FALSE, recursive = FALSE)
    
    for (reg in region_names) {
      reg_path <- file.path(gender_path, reg)
      def_types <- list.dirs(reg_path, full.names = FALSE, recursive = FALSE)
      
      # Get Base Year Mean for Target Calc
      base_means <- list()
      
      for (dt in def_types) {
        base_file <- file.path(reg_path, dt, paste0(BASE_YEAR, ".rds"))
        if (file.exists(base_file)) {
          mat <- readRDS(base_file)
          # Get regional weights for aggregation to scalar
          w <- regional_population_weights %>% 
            filter(region_consolidated == reg, sex == gender, year == BASE_YEAR) %>% 
            arrange(age) %>% pull(regional_age_weight)
          
          if(length(w) == nrow(mat)) {
            p_mat <- plogis(mat)
            scalar_p <- colSums(p_mat * w)
            base_means[[dt]] <- mean(scalar_p)
          }
        }
      }
      
      # Process Target Year
      for (dt in def_types) {
        if (is.null(base_means[[dt]])) next
        
        target_file <- file.path(reg_path, dt, paste0(TARGET_YEAR, ".rds"))
        if (file.exists(target_file)) {
          mat <- readRDS(target_file)
          w <- regional_population_weights %>% 
            filter(region_consolidated == reg, sex == gender, year == TARGET_YEAR) %>% 
            arrange(age) %>% pull(regional_age_weight)
          
          if(length(w) == nrow(mat)) {
            p_mat <- plogis(mat)
            scalar_p <- colSums(p_mat * w)
            
            target_val <- base_means[[dt]] * (1 - REDUCTION_PERCENTAGE/100)
            prob_success <- mean(scalar_p < target_val)
            
            # Save Record
            results_list[[length(results_list)+1]] <- data.frame(
              Region = reg, Sex = gender, Def_Type_Code = dt,
              Base_Prevalence = base_means[[dt]],
              Target_Prevalence = target_val,
              Projected_Prevalence = mean(scalar_p),
              Prob_Achievement = prob_success
            )
          }
        }
      }
    }
  }
  return(do.call(rbind, results_list))
}

regional_results_df <- calc_regional_metrics()

# ---- 15.5 Calculate High-Certainty Country Percentages ----
# This addresses the Lancet requirement: "What % of countries in Region X are likely to succeed?"

country_achievement_summary <- final_weighted_results_selected %>% # From Section 15
  filter(Year == TARGET_YEAR) %>%
  mutate(High_Certainty = prob_achieving_target > 0.60) %>%
  left_join(country_region_mapping, by = c("Country" = "wb_country_abv")) %>%
  group_by(region_consolidated, Sex, Def_Type_Code) %>%
  summarise(
    Total_Countries = n(),
    Countries_High_Certainty = sum(High_Certainty, na.rm = TRUE),
    High_Certainty_Pct = (sum(High_Certainty, na.rm = TRUE) / n()) * 100,
    .groups = "drop"
  ) %>%
  rename(Region = region_consolidated)

# Merge Weighted Regional stats with Country Count stats
regional_full_summary <- regional_results_df %>%
  left_join(country_achievement_summary, by = c("Region", "Sex", "Def_Type_Code"))

write.csv(regional_full_summary, "results/regional_aggregation/regional_full_summary_lancet.csv", row.names = FALSE)

#########################################################################################
#                          PUBLICATION TABLES (GT)                                      #
#########################################################################################

cat("  Generating Publication-Ready GT Tables...\n")

# ============================================================================
# NEW: Source the comprehensive publication tables module
# This generates 8 tables following PhD thesis design guidelines:
#   - Main Tables 1-3: Primary indicator (current_user_any_tobacco_product)
#   - Supplementary Tables S1-S5: All 6 indicators, stratum-level details
# ============================================================================

tryCatch({
  # The publication_tables.R script expects clean_data and country_region_mapping
  # to be available in the environment

  if (file.exists("R/publication_tables.R")) {
    cat("  Sourcing R/publication_tables.R...\n")
    source("R/publication_tables.R")
  } else {
    cat("  WARNING: R/publication_tables.R not found. Skipping table generation.\n")
  }
}, error = function(e) {
  cat(sprintf("  WARNING: Publication tables generation failed: %s\n", e$message))
  cat("  Continuing with pipeline...\n")
})

#########################################################################################
#                          VISUALIZATION (PUBLICATION FIGURES)                           #
#########################################################################################

cat("  Generating Publication Figures...\n")

# ============================================================================
# SOURCE THE COMPREHENSIVE VISUALIZATION MODULE
# This generates all publication figures following PhD thesis design guidelines:
#   - Figure 1: Stoplight World Map (target achievement probability)
#   - Figure 2: Caterpillar Plot (country ranking by reduction)
#   - Figure 3: Endgame Trajectories (spaghetti plot with <5% zone)
#   - Figure 4: Age-Cohort Heatmap (Lexis diagram)
#   - eFigure 1: Model Validation (observed vs predicted)
#   - eFigure 2: Country Age Profiles (fitted vs observed with age bands)
#   - eFigure 3: Weighted Trend Curves (time series with uncertainty)
# ============================================================================

tryCatch({
  # Source Module 12: Publication Figures (Maps, Caterpillar, Trajectories, Heatmaps)
  if (file.exists("R/12_publication_figures.R")) {
    cat("  Sourcing R/12_publication_figures.R...\n")
    source("R/12_publication_figures.R")
  } else {
    cat("  WARNING: R/12_publication_figures.R not found.\n")
  }

  # Source Module 13: Strata-Level Figures (Age Curves, Trends per Country/Gender/Year)
  if (file.exists("R/13_strata_figures.R")) {
    cat("  Sourcing R/13_strata_figures.R...\n")
    source("R/13_strata_figures.R")
  } else {
    cat("  WARNING: R/13_strata_figures.R not found.\n")
    cat("  Falling back to basic visualization...\n")

    # Basic fallback visualization
    if (exists("clean_data") && exists("final_ac_predictions")) {
      plot_countries <- unique(clean_data$wb_country_abv)

      for (cntry in plot_countries[1:min(5, length(plot_countries))]) {
        cntry_dir <- file.path("plots/age_curves", cntry)
        dir.create(cntry_dir, recursive = TRUE, showWarnings = FALSE)

        obs_data <- clean_data %>%
          filter(wb_country_abv == cntry, def_type_code == "current_user_cigarettes") %>%
          mutate(Prevalence = plogis(prevalence), Year = as.numeric(year))

        if (nrow(obs_data) == 0) next

        for (yr in unique(obs_data$Year)) {
          p_obs <- obs_data %>% filter(Year == yr)

          p <- ggplot(p_obs, aes(x = Age_Midpoint, y = Prevalence)) +
            geom_point(size = 2, shape = 21, fill = "white") +
            geom_vline(xintercept = 65, linetype = "dotted", color = "grey50") +
            facet_wrap(~sex) +
            theme_minimal() +
            labs(title = paste0(cntry, " (", yr, ")"),
                 y = "Prevalence", x = "Age",
                 subtitle = "Current Cigarette Use")

          ggsave(file.path(cntry_dir, paste0("AgeCurve_", yr, ".png")),
                 p, width = 8, height = 5)
        }
      }
    }
  }
}, error = function(e) {
  cat(sprintf("  WARNING: Figure generation failed: %s\n", e$message))
  cat("  Continuing with pipeline...\n")
})

cat("\n================================================================\n")
cat("  PIPELINE COMPLETE\n")
cat("  - Publication tables generated (8 tables: 3 main + 5 supplementary)\n")
cat("  - Post-selection predictions consolidated.\n")
cat("  - Probability maps generated:\n")
cat("      * 6 indicators x 2 genders x 2 targets = 24 maps\n")
cat("      * 2025 WHO Target (30% reduction from 2010)\n")
cat("      * 2040 Endgame (<5% prevalence)\n")
cat("  - Age-Cohort visualizations saved.\n")
cat("================================================================\n")

#########################################################################################
#                                                                                       #
#             SECTION 16: ADVANCED VISUALIZATION PIPELINE                               #
#             Replaces Python modules: trends, validation, aggregated, cohorts, maps    #
#                                                                                       #
#########################################################################################

# ---- 16.1 Style Configuration ----
# Note: All required packages loaded in Section 1.1

# Color Palettes (Aligned with WHO/Lancet style)
who_colors <- list(
  primary   = "#2C3E50",
  secondary = "#3498DB",
  accent    = "#E74C3C",
  success   = "#27AE60",
  warning   = "#F39C12",
  gray      = "#95A5A6"
)

# Region Colors
region_colors <- c(
  "Africa" = "#1f77b4", "Americas" = "#ff7f0e", "Eastern Mediterranean" = "#2ca02c",
  "Europe" = "#d62728", "South-East Asia" = "#9467bd", "Western Pacific" = "#8c564b",
  "Global" = "#2C3E50" # Dark slate for Global
)

# Indicator Labels Map
format_indicator <- function(x) {
  case_when(
    x == "daily_user_cigarettes" ~ "Daily Cigarettes",
    x == "current_user_cigarettes" ~ "Current Cigarettes",
    x == "daily_user_any_tobacco_product" ~ "Daily Any Tobacco",
    x == "current_user_any_tobacco_product" ~ "Current Any Tobacco",
    TRUE ~ gsub("_", " ", str_to_title(x))
  )
}

# The Master Theme Function
theme_who <- function(base_size = 11, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Typography
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0, margin = margin(b = 5)),
      plot.subtitle = element_text(color = "#555555", size = rel(1.0), hjust = 0, margin = margin(b = 10)),
      plot.caption = element_text(color = "#777777", size = rel(0.7), hjust = 1, margin = margin(t = 10)),
      axis.title = element_text(face = "bold", size = rel(0.9)),
      axis.text = element_text(color = "#333333", size = rel(0.85)),
      
      # Clean Layout (No top/right spines)
      panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "#333333", linewidth = 0.3),
      
      # Legend
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_text(face = "bold", size = rel(0.8)),
      
      # Facets
      strip.background = element_rect(fill = "#F5F5F5", color = NA),
      strip.text = element_text(face = "bold", size = rel(0.9), hjust = 0, margin = margin(4,4,4,4))
    )
}
theme_set(theme_who())

#########################################################################################
#                              SECTION 17: VISUALIZATION DATA LOADING                   #
#########################################################################################

cat("  Loading data for visualization pipeline...\n")

# 1. Load Predictions (Age-Specific) - Tier 1
df_preds_global <- tryCatch({
  read.csv("results/final_ac_predictions_nested.csv") %>% mutate(Model_Type = "Global")
}, error = function(e) {
  cat("    WARNING: Could not load final_ac_predictions_nested.csv\n")
  data.frame()
})

df_preds_country <- tryCatch({
  read.csv("results/final_predictions_country_specific.csv") %>% mutate(Model_Type = "Country")
}, error = function(e) {
  cat("    WARNING: Could not load final_predictions_country_specific.csv\n")
  data.frame()
})

df_preds_all <- if (nrow(df_preds_global) > 0 || nrow(df_preds_country) > 0) {
  bind_rows(df_preds_global, df_preds_country)
} else {
  # Try alternate sources
  if (exists("final_ac_predictions")) {
    final_ac_predictions %>% mutate(Model_Type = "Global")
  } else {
    data.frame()
  }
}

# 2. Load Trends (Weighted) - Tier 2
df_trends <- tryCatch({
  read.csv("results/final_weighted_results_apc_post_selection.csv")
}, error = function(e) {
  cat("    WARNING: Could not load final_weighted_results_apc_post_selection.csv\n")
  if (exists("final_weighted_results_selected")) {
    final_weighted_results_selected
  } else {
    data.frame()
  }
})

# 3. Load Regional - Tier 3
df_regional <- tryCatch({
  read.csv("results/regional_aggregation/regional_full_summary_lancet.csv")
}, error = function(e) {
  cat("    WARNING: Could not load regional_full_summary_lancet.csv\n")
  if (exists("regional_full_summary")) {
    regional_full_summary
  } else {
    data.frame()
  }
})

# 4. Load Observations (from clean_data if available)
df_obs <- if (exists("clean_data") && nrow(clean_data) > 0) {
  clean_data %>%
    mutate(
      Year = as.numeric(year),
      Prevalence = plogis(prevalence),
      Indicator = def_type_code,
      Country = wb_country_abv
    )
} else {
  data.frame()
}

cat(sprintf("  Loaded: %d predictions, %d trend rows, %d regional rows, %d observations\n",
            nrow(df_preds_all), nrow(df_trends), nrow(df_regional), nrow(df_obs)))

#########################################################################################
#                              SECTION 18: PLOT MODULES                                 #
#########################################################################################

# ---- Module 1: Trend Plots (Trends.R equivalent) ----

plot_age_panel <- function(country_code, gender, indicator) {
  # Filter Data
  dat_p <- df_preds_all %>%
    filter(Country == country_code, Sex == gender, Def_Type_Code == indicator) %>%
    filter(Age_Midpoint %in% seq(20, 80, by = 10)) # Select specific ages for panel
  
  if(nrow(dat_p) == 0) return(NULL)
  
  country_name <- country_name_mapping[country_code]
  
  p <- ggplot(dat_p, aes(x = Year, y = Prevalence)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Model_Type), alpha = 0.2) +
    geom_line(aes(color = Model_Type), linewidth = 0.8) +
    facet_wrap(~Age_Midpoint, labeller = label_both, ncol = 4) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_color_manual(values = c("Country" = who_colors$accent, "Global" = who_colors$secondary)) +
    scale_fill_manual(values = c("Country" = who_colors$accent, "Global" = who_colors$secondary)) +
    labs(
      title = paste0(country_name, ": Age-Specific Trends"),
      subtitle = paste0(format_indicator(indicator), " | ", tools::toTitleCase(gender)),
      caption = "Shaded area: 95% Credible Interval"
    )
  return(p)
}

plot_age_year_heatmap <- function(country_code, gender, indicator) {
  dat_p <- df_preds_all %>%
    filter(Country == country_code, Sex == gender, Def_Type_Code == indicator, 
           Model_Type == "Country") # Use Country model for heatmap
  
  if(nrow(dat_p) == 0) return(NULL)
  country_name <- country_name_mapping[country_code]
  
  p <- ggplot(dat_p, aes(x = Year, y = Age_Midpoint, fill = Prevalence)) +
    geom_tile() +
    scale_fill_viridis(option = "magma", labels = percent, direction = -1) +
    scale_x_continuous(expand = c(0,0), breaks = seq(2010, 2030, 5)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(
      title = paste0(country_name, ": Prevalence Heatmap"),
      subtitle = paste0(format_indicator(indicator), " | ", tools::toTitleCase(gender)),
      y = "Age"
    ) +
    coord_fixed(ratio = 0.5) # Adjust aspect ratio
  return(p)
}

# ---- Module 2: Validation Plots (Validation.R equivalent) ----

plot_validation_facet <- function(country_code, gender, year) {
  # Get predictions for specific year
  dat_p <- df_preds_all %>%
    filter(Country == country_code, Sex == gender, Year == year)

  # Get observations close to that year (+- 2 years)
  dat_o <- df_obs %>%
    filter(wb_country_abv == country_code, sex == gender, abs(Year - year) <= 2)

  if(nrow(dat_p) == 0 && nrow(dat_o) == 0) return(NULL)

  # Safely get country name
  country_name <- if (exists("country_name_mapping") && country_code %in% names(country_name_mapping)) {
    country_name_mapping[country_code]
  } else {
    toupper(country_code)
  }

  # Ensure start_age and end_age exist
  if (!"start_age" %in% names(dat_o)) {
    dat_o$start_age <- dat_o$Age_Midpoint - 5
  }
  if (!"end_age" %in% names(dat_o)) {
    dat_o$end_age <- dat_o$Age_Midpoint + 5
  }

  p <- ggplot() +
    # Background grid
    geom_hline(yintercept = seq(0, 1, 0.2), color = "grey95", linewidth = 0.2)

  # Add predictions if available
  if (nrow(dat_p) > 0) {
    p <- p +
      # Prediction Ribbon
      geom_ribbon(data = dat_p,
                  aes(x = Age_Midpoint, ymin = lower_ci, ymax = upper_ci, fill = Model_Type),
                  alpha = 0.15) +
      geom_line(data = dat_p,
                aes(x = Age_Midpoint, y = Prevalence, color = Model_Type, linetype = Model_Type),
                linewidth = 0.8)
  }

  # Add observed data if available
  if (nrow(dat_o) > 0) {
    p <- p +
      # Age Range bars (showing survey age bands)
      geom_segment(data = dat_o,
                   aes(x = start_age, xend = end_age,
                       y = Prevalence, yend = Prevalence),
                   color = "black", linewidth = 0.5) +
      # Observed midpoints
      geom_point(data = dat_o,
                 aes(x = Age_Midpoint, y = Prevalence),
                 size = 1.5, shape = 21, fill = "white", stroke = 0.5)
  }

  # Add transition marker
  p <- p +
    geom_vline(xintercept = 65, linetype = "dashed", color = "orange", alpha = 0.6) +
    facet_wrap(~Def_Type_Code, labeller = as_labeller(format_indicator), scales = "free_y") +
    scale_y_continuous(labels = percent, limits = c(0, NA)) +
    scale_x_continuous(limits = c(15, 85)) +
    scale_color_manual(values = c("Country" = who_colors$accent, "Global" = who_colors$secondary),
                       name = "Model") +
    scale_fill_manual(values = c("Country" = who_colors$accent, "Global" = who_colors$secondary),
                      name = "Model") +
    scale_linetype_manual(values = c("Country" = "dashed", "Global" = "solid"),
                          name = "Model") +
    theme_who() +
    labs(
      title = paste0("Model Validation: ", country_name, " (", year, ")"),
      subtitle = paste0("Fitted curves vs survey data | ", tools::toTitleCase(gender)),
      x = "Age", y = "Prevalence",
      caption = "Horizontal bars show survey age bands. Orange line = spline/linear transition."
    )
  return(p)
}

# ---- Module 3: Aggregated Plots (Aggregated.R equivalent) ----

plot_regional_comparison <- function(indicator, gender) {
  dat_r <- df_regional %>%
    filter(Sex == gender, Def_Type_Code == indicator) %>%
    mutate(Region_Name = tools::toTitleCase(Region))
  
  # Add Global Trend (Tier 2 Data aggregated)
  # (Simulated here based on structure, assumes you calculate global avg elsewhere or add it to dataframe)
  
  p <- ggplot(dat_r) +
    # Regional Lines
    geom_segment(aes(x = Base_Prevalence, xend = Projected_Prevalence, y = reorder(Region_Name, Projected_Prevalence), yend = reorder(Region_Name, Projected_Prevalence)), color = "grey") +
    geom_point(aes(x = Base_Prevalence, y = Region_Name, color = "2010 Base"), size = 3) +
    geom_point(aes(x = Projected_Prevalence, y = Region_Name, color = "2025 Projected"), size = 3) +
    geom_vline(xintercept = 0.30, linetype = "dashed", color = "grey") + # Example reference
    
    scale_x_continuous(labels = percent) +
    scale_color_manual(values = c("2010 Base" = who_colors$primary, "2025 Projected" = who_colors$success)) +
    labs(
      title = paste0("Regional Shift: ", format_indicator(indicator)),
      subtitle = paste0(tools::toTitleCase(gender), " | 2010 vs 2025"),
      x = "Prevalence", y = NULL, color = "Year"
    )
  return(p)
}

# ---- Module 4: Cohort Plots (Cohorts.R equivalent) ----

plot_lexis_diagram <- function(country_code, gender, indicator) {
  dat_p <- df_preds_all %>%
    filter(Country == country_code, Sex == gender, Def_Type_Code == indicator, Model_Type == "Country")
  
  if(nrow(dat_p) == 0) return(NULL)
  
  # Create Cohort Diagonals (y = x - c  => Age = Year - Cohort)
  cohorts <- seq(1940, 2000, by = 20)
  
  p <- ggplot(dat_p, aes(x = Year, y = Age_Midpoint)) +
    geom_tile(aes(fill = Prevalence)) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1, labels = percent) +
    
    # Add diagonal cohort lines
    geom_abline(intercept = -cohorts, slope = 1, color = "white", linetype = "dashed", alpha = 0.5) +
    annotate("text", x = 2015, y = 2015 - cohorts, label = paste0("Born ", cohorts), 
             angle = 45, color = "white", size = 3, vjust = -0.5) +
    
    coord_fixed(ratio = 1, xlim = c(2000, 2030), ylim = c(15, 80)) +
    labs(
      title = paste0("Lexis Diagram: ", country_name_mapping[country_code]),
      subtitle = paste0("Diagonal lines track specific birth cohorts"),
      fill = "Prev."
    )
  return(p)
}

# ---- Module 5: Maps (Maps.R equivalent) ----
# NOTE: Comprehensive map functions are in R/12_publication_figures.R
# This section provides inline fallback and basic map functionality

# ============================================================================
# COUNTRY CODE HARMONIZATION (WHO -> Natural Earth ISO alpha-3)
# ============================================================================
# Comprehensive mapping based on thorough comparison of WHO tobacco data
# with Natural Earth map polygons. Fixes all known mismatches.

COUNTRY_CODE_MAP <- c(
  # Verified mismatches
  "rom" = "ROU",   # Romania (old code)
  "tmp" = "TLS",   # Timor-Leste (old code)
  "tls" = "TLS",   # Timor-Leste
  "pse" = "PSE",   # Palestine
  "wbg" = "PSE",   # West Bank & Gaza
  "zar" = "COD",   # DR Congo (old Zaire)
  "cod" = "COD",   # DR Congo
  "ksv" = "XKX",   # Kosovo (NE uses XKX)
  "xkx" = "XKX",   # Kosovo

  # East Asia
  "kor" = "KOR", "prk" = "PRK", "twn" = "TWN", "hkg" = "HKG", "mac" = "MAC",

  # Southeast Asia
  "lao" = "LAO", "mmr" = "MMR", "vnm" = "VNM", "khm" = "KHM", "brn" = "BRN",

  # Africa
  "civ" = "CIV", "swz" = "SWZ", "ssd" = "SSD",

  # Europe
  "mkd" = "MKD", "srb" = "SRB", "mne" = "MNE", "bih" = "BIH",
  "and" = "AND", "smr" = "SMR",

  # Pacific Islands (many too small for NE medium scale)
  "mhl" = "MHL", "fsm" = "FSM", "nru" = "NRU", "niu" = "NIU", "plw" = "PLW",
  "wsm" = "WSM", "tuv" = "TUV", "cok" = "COK", "kir" = "KIR", "ton" = "TON",
  "vut" = "VUT", "slb" = "SLB",

  # Caribbean
  "kna" = "KNA", "lca" = "LCA", "vct" = "VCT", "grd" = "GRD", "dma" = "DMA",
  "atg" = "ATG", "brb" = "BRB",

  # Middle East
  "are" = "ARE", "bhr" = "BHR", "kwt" = "KWT", "qat" = "QAT", "omn" = "OMN",

  # Indian Ocean
  "mdv" = "MDV", "syc" = "SYC", "mus" = "MUS", "com" = "COM", "stp" = "STP"
)

# ============================================================================
# INDICATOR CONFIGURATION FOR MAPS
# ============================================================================

# All 6 indicators
ALL_INDICATORS <- c(
  "current_user_any_tobacco_product",
  "current_user_any_smoked_tobacco",
  "current_user_cigarettes",

"daily_user_any_tobacco_product",
  "daily_user_any_smoked_tobacco",
  "daily_user_cigarettes"
)

# Indicator labels for display
INDICATOR_LABELS_MAP <- c(
  "current_user_any_tobacco_product" = "Current Any Tobacco",
  "current_user_any_smoked_tobacco" = "Current Smoked Tobacco",
  "current_user_cigarettes" = "Current Cigarettes",
  "daily_user_any_tobacco_product" = "Daily Any Tobacco",
  "daily_user_any_smoked_tobacco" = "Daily Smoked Tobacco",
  "daily_user_cigarettes" = "Daily Cigarettes"
)

# Short codes for filenames
INDICATOR_CODES_MAP <- c(
  "current_user_any_tobacco_product" = "CU_ATP",
  "current_user_any_smoked_tobacco" = "CU_AST",
  "current_user_cigarettes" = "CU_CIG",
  "daily_user_any_tobacco_product" = "DU_ATP",
  "daily_user_any_smoked_tobacco" = "DU_AST",
  "daily_user_cigarettes" = "DU_CIG"
)

# Stoplight colors for probability categories
STOPLIGHT_COLORS <- c(
  "Extremely Low (<10%)" = "#B71C1C",
  "Low (10-40%)" = "#E53935",
  "Uncertain (40-60%)" = "#FFC107",
  "Moderate (60-90%)" = "#8BC34A",
  "High (>90%)" = "#1B5E20",
  "No Data" = "#BDBDBD"
)

# Target year configuration
TARGET_YEAR <- 2025
ENDGAME_YEAR <- 2040
BASE_YEAR <- 2010

harmonize_country_code_for_map <- function(country_code) {
  cc <- tolower(country_code)
  if (cc %in% names(COUNTRY_CODE_MAP)) {
    return(COUNTRY_CODE_MAP[cc])
  }
  return(toupper(cc))
}

plot_prevalence_map_dual <- function(year, indicator) {
  world <- ne_countries(scale = "medium", returnclass = "sf")

  # Fix date line issues
  sf::sf_use_s2(FALSE)
  world <- tryCatch({
    st_crop(world, xmin = -180, xmax = 180, ymin = -90, ymax = 90)
  }, error = function(e) world)
  sf::sf_use_s2(TRUE)

  # Fix Natural Earth iso_a3 = "-99" issue for France, Norway, etc.
  # Use adm0_a3 as fallback when iso_a3 is invalid
  world <- world %>%
    mutate(
      iso_a3 = case_when(
        is.na(iso_a3) | iso_a3 == "-99" ~ adm0_a3,
        TRUE ~ iso_a3
      )
    )

  # Convert WHO country codes to ISO alpha-3 with harmonization
  dat_map <- df_trends %>%
    filter(Year == year, Def_Type_Code == indicator) %>%
    mutate(iso_a3 = sapply(Country, harmonize_country_code_for_map)) %>%
    select(iso_a3, Sex, weighted_mean)

  # Get valid Natural Earth codes (now with fix applied)
  ne_iso3 <- unique(world$iso_a3)
  ne_iso3 <- ne_iso3[!is.na(ne_iso3) & ne_iso3 != "-99"]

  # Diagnostic: which WHO countries matched?
  who_codes <- unique(dat_map$iso_a3)
  matched_codes <- who_codes[who_codes %in% ne_iso3]
  unmatched_codes <- who_codes[!who_codes %in% ne_iso3]

  cat(sprintf("    Map: %d/%d countries matched to Natural Earth\n",
              length(matched_codes), length(who_codes)))

  if (length(unmatched_codes) > 0 && length(unmatched_codes) <= 10) {
    cat(sprintf("    Unmatched: %s\n", paste(unmatched_codes, collapse = ", ")))
  } else if (length(unmatched_codes) > 10) {
    cat(sprintf("    Unmatched: %d countries (mostly small islands)\n", length(unmatched_codes)))
  }

  # Join on iso_a3
  map_data <- world %>%
    left_join(dat_map, by = "iso_a3") %>%
    filter(!is.na(Sex))

  p <- ggplot(map_data) +
    geom_sf(aes(fill = weighted_mean), color = "white", linewidth = 0.1) +
    facet_wrap(~Sex, ncol = 1) +
    scale_fill_viridis(option = "rocket", direction = -1, labels = percent, name = "Prevalence") +
    theme_void() +
    theme(legend.position = "right") +
    labs(
      title = paste0("Global Prevalence: ", format_indicator(indicator)),
      subtitle = paste0("Year: ", year)
    )
  return(p)
}

#########################################################################################
#                              SECTION 19: MASTER EXECUTION LOOP                        #
#########################################################################################

generate_all_visualizations <- function(output_base = "outputs/figures") {

  # Create directory structure (includes new map subdirectories)
  dirs <- c(
    "trends/panels", "trends/heatmaps", "validation", "aggregated", "cohorts",
    "maps",
    "maps/2025_target/men", "maps/2025_target/women",
    "maps/2040_endgame/men", "maps/2040_endgame/women",
    "maps/composites"
  )
  for(d in dirs) {
    dir.create(file.path(output_base, d), recursive = TRUE, showWarnings = FALSE)
  }

  cat("  Starting batch generation...\n")

  # 1. PROBABILITY MAPS (All 6 indicators × 2 targets × 2 genders)
  cat("  Generating Probability Maps...\n")

  # Try to use the comprehensive map function from 12_publication_figures.R
  if (exists("generate_all_probability_maps") && exists("df_trends") && nrow(df_trends) > 0) {
    tryCatch({
      cat("    Using generate_all_probability_maps() for full map suite...\n")
      generate_all_probability_maps(
        weighted_results = df_trends,
        mcmc_samples = NULL,
        output_dir = file.path(output_base, "maps"),
        primary_only = FALSE  # Generate all 6 indicators
      )
    }, error = function(e) {
      cat(sprintf("    Probability maps error: %s\n", e$message))
      cat("    Falling back to basic prevalence maps...\n")
      # Fallback to basic maps
      for(ind in unique(df_trends$Def_Type_Code)) {
        tryCatch({
          p <- plot_prevalence_map_dual(TARGET_YEAR, ind)
          if (!is.null(p)) {
            ggsave(file.path(output_base, "maps", paste0("map_", ind, "_", TARGET_YEAR, ".pdf")),
                   p, width = 8, height = 10)
          }
        }, error = function(e2) {})
      }
    })
  } else if (exists("df_trends") && nrow(df_trends) > 0) {
    # Basic prevalence maps fallback
    cat("    Using basic prevalence maps (generate_all_probability_maps not available)...\n")
    for(ind in unique(df_trends$Def_Type_Code)) {
      tryCatch({
        p <- plot_prevalence_map_dual(TARGET_YEAR, ind)
        if (!is.null(p)) {
          ggsave(file.path(output_base, "maps", paste0("map_", ind, "_", TARGET_YEAR, ".pdf")),
                 p, width = 8, height = 10)
        }
      }, error = function(e) {
        cat(sprintf("    Map error (%s): %s\n", ind, e$message))
      })
    }
  } else {
    cat("    Skipping maps (no df_trends data)\n")
  }

  # 2. AGGREGATED (Regional)
  cat("  Generating Regional Plots...\n")
  if (exists("df_regional") && nrow(df_regional) > 0) {
    for(ind in unique(df_regional$Def_Type_Code)) {
      for(g in c("males", "females")) {
        tryCatch({
          p <- plot_regional_comparison(ind, g)
          if (!is.null(p)) {
            ggsave(file.path(output_base, "aggregated", paste0("regional_", ind, "_", g, ".pdf")),
                   p, width = 10, height = 6)
          }
        }, error = function(e) {})
      }
    }
  } else {
    cat("    Skipping regional plots (no df_regional data)\n")
  }

  # 3. COUNTRY LOOPS
  if (exists("df_preds_all") && nrow(df_preds_all) > 0) {
    countries <- unique(df_preds_all$Country)
    countries <- countries[!is.na(countries) & countries != ""]
    cat(sprintf("  Generating Country Plots for %d countries...\n", length(countries)))

    pb <- txtProgressBar(min = 0, max = length(countries), style = 3)

    for(i in seq_along(countries)) {
      cntry <- countries[i]

      # Safely get country name
      c_name <- if (exists("country_name_mapping") && cntry %in% names(country_name_mapping)) {
        gsub(" ", "_", country_name_mapping[cntry])
      } else {
        toupper(cntry)
      }

      # Create Country Subfolder
      c_dir <- file.path(output_base, "trends", c_name)
      dir.create(c_dir, showWarnings = FALSE, recursive = TRUE)

      for(g in c("males", "females")) {
        for(ind in unique(df_preds_all$Def_Type_Code)) {
          tryCatch({
            # Trend Panel
            p1 <- plot_age_panel(cntry, g, ind)
            if(!is.null(p1)) {
              ggsave(file.path(c_dir, paste0("panel_", g, "_", ind, ".pdf")),
                     p1, width = 12, height = 8)
            }

            # Heatmap
            p2 <- plot_age_year_heatmap(cntry, g, ind)
            if(!is.null(p2)) {
              ggsave(file.path(c_dir, paste0("heatmap_", g, "_", ind, ".png")),
                     p2, width = 8, height = 6)
            }

            # Cohort Lexis
            if(ind == "daily_user_cigarettes" || ind == "current_user_cigarettes") {
              p3 <- plot_lexis_diagram(cntry, g, ind)
              if(!is.null(p3)) {
                ggsave(file.path(output_base, "cohorts", paste0("lexis_", c_name, "_", g, ".pdf")),
                       p3, width = 8, height = 8)
              }
            }

          }, error = function(e) {})
        }

        # Validation Plot
        tryCatch({
          # Get available years for this country
          obs_years <- if(exists("df_obs")) {
            unique(df_obs$Year[df_obs$wb_country_abv == cntry | df_obs$Country == cntry])
          } else {
            c(2015)
          }
          if (length(obs_years) > 0) {
            yr <- obs_years[ceiling(length(obs_years)/2)]  # Middle year
            p_val <- plot_validation_facet(cntry, g, yr)
            if(!is.null(p_val)) {
              ggsave(file.path(output_base, "validation", paste0("val_", c_name, "_", g, "_", yr, ".pdf")),
                     p_val, width = 10, height = 6)
            }
          }
        }, error = function(e) {})
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)
  } else {
    cat("    Skipping country plots (no df_preds_all data)\n")
  }

  cat("\n  Visualization generation complete.\n")
}

# Run the Master Function (with error handling)
tryCatch({
  generate_all_visualizations()
}, error = function(e) {
  cat(sprintf("  WARNING: Visualization batch failed: %s\n", e$message))
  cat("  Continuing with pipeline...\n")
})

#########################################################################################
#                          FINAL DIAGNOSTICS SUMMARY                                     #
#########################################################################################

cat("\n")
cat("###########################################################################\n")
cat("#                                                                         #\n")
cat("#                    FINAL DIAGNOSTICS SUMMARY                            #\n")
cat("#                                                                         #\n")
cat("###########################################################################\n\n")

# Generate comprehensive diagnostics summary if module available
if (exists("generate_all_models_summary") && DIAGNOSTICS_AVAILABLE) {
  tryCatch({
    cat("=== CONVERGENCE SUMMARY ACROSS ALL MODELS ===\n")
    all_models_summary <- generate_all_models_summary()

    # Generate worst parameters report
    cat("\n=== WORST PARAMETERS ACROSS ALL MODELS ===\n")
    worst_params <- generate_worst_parameters_report(top_n = 50)
    if (!is.null(worst_params)) {
      cat("\nTop 20 worst parameters:\n")
      print(head(as.data.frame(worst_params), 20), row.names = FALSE)
    }

    # Generate elegant GT table if gt package available
    if (requireNamespace("gt", quietly = TRUE)) {
      tryCatch({
        gt_table <- generate_convergence_gt_table(
          output_file = file.path("diagnostics", "tables", "convergence_summary_table.html")
        )
        cat("\n  Elegant GT table saved to diagnostics/tables/convergence_summary_table.html\n")
      }, error = function(e) {
        cat(sprintf("  Note: GT table generation skipped: %s\n", e$message))
      })
    }

  }, error = function(e) {
    cat(sprintf("  WARNING: Final diagnostics summary failed: %s\n", e$message))
  })
} else {
  cat("  Diagnostics module not available. See individual model outputs.\n")
}

#########################################################################################
#                              PIPELINE COMPLETE                                         #
#########################################################################################

pipeline_end_time <- Sys.time()
if (exists("pipeline_start_time")) {
  total_time <- difftime(pipeline_end_time, pipeline_start_time, units = "hours")
  cat(sprintf("\n  Total pipeline runtime: %.2f hours\n", as.numeric(total_time)))
}

cat("\n")
cat("###########################################################################\n")
cat("#                                                                         #\n")
cat("#                    PIPELINE EXECUTION COMPLETE                          #\n")
cat("#                                                                         #\n")
cat("###########################################################################\n\n")

cat("OUTPUT DIRECTORIES:\n")
cat("  - diagnostics/logs/      MCMC convergence log files\n")
cat("  - diagnostics/tables/    Master R-hat and ESS tables\n")
cat("  - results/               Final predictions and aggregations\n")
cat("  - figures_publication/   Publication-quality figures\n")
cat("  - tables_publication/    Publication-quality tables\n")
cat("  - processing/            Per-country prediction files\n")
cat("  - checkpoints/           MCMC checkpoints for resumption\n")

cat("\nKEY FILES:\n")
cat("  - diagnostics/tables/rhat_master.csv           All R-hat values\n")
cat("  - diagnostics/tables/convergence_summary.csv   Per-model summary\n")
cat("  - results/final_weighted_results_apc_post_selection.csv\n")
cat("  - results/final_ac_predictions_nested.csv\n")

# Play completion sound if available
if (requireNamespace("beepr", quietly = TRUE)) {
  tryCatch(beepr::beep(sound = "complete"), error = function(e) {})
}

cat("\n  WHO Tobacco Pipeline complete.\n\n")