#########################################################################################
#
#                    WHO TOBACCO CONTROL PREVALENCE PROJECTION MODEL
#                              00_config.R - Configuration
#                                   VERSION 2.3.2
#
#   Contains: All constants, MCMC settings, and output directory structure
#   Source this file first before any other module
#
#   CRITICAL: This file must be sourced BEFORE loading nimble to set DLL limit
#
#########################################################################################

# ---- CRITICAL: Increase DLL Limit BEFORE Loading NIMBLE ----
# Each NIMBLE model compiles to a separate DLL. With 191 countries × 2 sexes = 382+ models,
# we exceed R's default limit of 100. Must be set before loading nimble.
Sys.setenv(R_MAX_NUM_DLLS = 600)

# ---- Package Loading Function ----

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

# ---- NIMBLE Configuration ----

nimbleOptions(verbose = FALSE)
nimbleOptions(MCMCprogressBar = TRUE)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = FALSE)

# ---- MCMC Settings ----
# Production settings for reliable convergence

NUMBER_OF_CHAINS     <- 4
NUMBER_OF_ADAPT      <- 1000   # Informational only - NIMBLE handles adaptation internally
NUMBER_OF_BURN       <- 5000
NUMBER_OF_ITERATIONS <- 10000
THINNING_INTERVAL    <- 5

# ---- Analysis Parameters ----

BASE_YEAR              <- 2010
TARGET_YEAR            <- 2025
PROJECTED_YEARS        <- 20
REDUCTION_PERCENTAGE   <- 30
RANDOM_SEED            <- 42

# ---- Age Transition Parameters ----
# Smooth sigmoid for elderly age extrapolation

TRANSITION_START  <- 65
TRANSITION_WIDTH  <- 3
MAX_AGE_SPLINE    <- 80

# ---- Absolute Prevalence Target ----

MANUAL_TARGET_PREVALENCE  <- 4.0
MANUAL_TARGET_PROPORTION  <- MANUAL_TARGET_PREVALENCE / 100
MANUAL_TARGET_EVAL_YEAR   <- 2040

# ---- Legacy Variable (used throughout code) ----

target_year <- TARGET_YEAR

# ---- Output Directory Structure ----

output_dirs <- c(
  # Diagnostics outputs
  "diagnostics",
  "diagnostics/logs",
  "diagnostics/tables",
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
  "figures_publication",
  "plots/world_maps",
  "plots/birth_cohort_analysis",
  "plots/age_curves",
  # Legacy directories
  "results",
  "evaluation",
  "checkpoints"
)

invisible(lapply(output_dirs, function(d) dir.create(d, recursive = TRUE, showWarnings = FALSE)))

# ---- Source MCMC Diagnostics Module ----

DIAGNOSTICS_AVAILABLE <- FALSE
if (file.exists("R/mcmc_diagnostics.R")) {
  tryCatch({
    source("R/mcmc_diagnostics.R")
    DIAGNOSTICS_AVAILABLE <- TRUE
  }, error = function(e) {
    cat(sprintf("  Note: Could not load mcmc_diagnostics.R: %s\n", e$message))
  })
} else if (file.exists("mcmc_diagnostics.R")) {
  tryCatch({
    source("mcmc_diagnostics.R")
    DIAGNOSTICS_AVAILABLE <- TRUE
  }, error = function(e) {
    cat(sprintf("  Note: Could not load mcmc_diagnostics.R: %s\n", e$message))
  })
}

# ---- Print Configuration Summary ----

print_config <- function() {
  cat("
#########################################################################
#
#              WHO TOBACCO CONTROL PREVALENCE PROJECTION
#                        VERSION 2.3.2 (NIMBLE)
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

DIAGNOSTICS:
------------
  Module loaded:       ", ifelse(DIAGNOSTICS_AVAILABLE, "Yes", "No"), "
  Log output:          diagnostics/logs/
  Table output:        diagnostics/tables/

ENGINE: NIMBLE (C++ compiled MCMC)

#########################################################################
")
}

# Print on source
print_config()
