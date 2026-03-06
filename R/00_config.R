#########################################################################################
#
#                    WHO TOBACCO CONTROL PREVALENCE PROJECTION MODEL
#                              00_config.R - Configuration
#                                   VERSION 2.4.0
#
#   Contains: All constants, MCMC settings, and output directory structure
#   Source this file first before any other module
#
#   CRITICAL: This file must be sourced BEFORE loading nimble to set DLL limit
#
#########################################################################################

# ---- CRITICAL: DLL Limit for NIMBLE ----
# Each NIMBLE model compiles to a separate DLL. With 191 countries × 2 sexes = 382+ models,
# we exceed R's default limit. The actual limit is set in .Renviron (project root) because
# R reads R_MAX_NUM_DLLS at startup — Sys.setenv() is too late on Windows.
# This call is kept as a fallback for systems that support dynamic changes.
Sys.setenv(R_MAX_NUM_DLLS = 1000)
dll_limit <- as.integer(Sys.getenv("R_MAX_NUM_DLLS", "100"))
cat(sprintf("  R_MAX_NUM_DLLS = %d (need ~800 for full pipeline)\n", dll_limit))
if (dll_limit < 800) {
  cat("  WARNING: DLL limit may be too low. Ensure .Renviron has R_MAX_NUM_DLLS=1000\n")
  cat("  WARNING: and RESTART R (not just re-source) for it to take effect.\n")
}

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
  # Parallel processing & subprocess isolation
  "foreach", "doParallel", "parallel", "callr",
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
NUMBER_OF_BURN       <- 30000  # Increased: block samplers need more adaptation time
NUMBER_OF_ITERATIONS <- 60000  # Increased: compensate for thinning and ensure stable ESS
THINNING_INTERVAL    <- 10     # Increased: block samplers may have higher initial autocorrelation

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

# ---- Publication Target Years ----

TARGET_YEAR_2030  <- 2030    # WHO 2030 target (used in publication tables)
ENDGAME_YEAR      <- 2040    # Tobacco endgame year

# ---- Target Thresholds (dual units for different contexts) ----

REDUCTION_TARGET          <- 0.30   # 30% relative reduction (proportion form)
# Proportion form: for data comparisons where prevalence is 0-1
ENDGAME_THRESHOLD_PROP    <- 0.05   # <5% = endgame achieved
NEAR_ENDGAME_PROP         <- 0.10   # 5-10% = near endgame
VIRTUAL_ELIMINATION_PROP  <- 0.02   # <2% = virtual elimination
ON_TRACK_MARGIN           <- 0.10   # Within 10% of target = on track
# Percentage form: for ggplot annotations where y-axis is 0-100
ENDGAME_THRESHOLD_PCT     <- 5
NEAR_ENDGAME_PCT          <- 10

# ---- Valid Indicators (Single Source of Truth) ----

VALID_INDICATORS <- c(
  "current_user_cigarettes",
  "current_user_any_smoked_tobacco",
  "current_user_any_tobacco_product",
  "daily_user_cigarettes",
  "daily_user_any_smoked_tobacco",
  "daily_user_any_tobacco_product"
)

PRIMARY_INDICATOR <- "current_user_cigarettes"

# ---- Indicator Display Labels & Codes ----

INDICATOR_LABELS <- c(
  "current_user_cigarettes"          = "Current Cigarettes",
  "current_user_any_smoked_tobacco"  = "Current Smoked Tobacco",
  "current_user_any_tobacco_product" = "Current Any Tobacco",
  "daily_user_cigarettes"            = "Daily Cigarettes",
  "daily_user_any_smoked_tobacco"    = "Daily Smoked Tobacco",
  "daily_user_any_tobacco_product"   = "Daily Any Tobacco"
)

INDICATOR_CODES <- c(
  "current_user_cigarettes"          = "CU_CIG",
  "current_user_any_smoked_tobacco"  = "CU_AST",
  "current_user_any_tobacco_product" = "CU_ATP",
  "daily_user_cigarettes"            = "DU_CIG",
  "daily_user_any_smoked_tobacco"    = "DU_AST",
  "daily_user_any_tobacco_product"   = "DU_ATP"
)

INDICATOR_COLORS <- c(
  "current_user_any_tobacco_product" = "#E41A1C",
  "current_user_any_smoked_tobacco"  = "#377EB8",
  "current_user_cigarettes"          = "#4DAF4A",
  "daily_user_any_tobacco_product"   = "#984EA3",
  "daily_user_any_smoked_tobacco"    = "#FF7F00",
  "daily_user_cigarettes"            = "#A65628"
)

# ---- Sex/Gender Standardization ----

VALID_SEX <- c("males", "females")
GENDER_DISPLAY <- c(males = "Men", females = "Women")

# ---- Display Helper Functions ----

format_gender <- function(x) {
  dplyr::case_when(
    tolower(x) %in% c("males", "male") ~ "Men",
    tolower(x) %in% c("females", "female") ~ "Women",
    TRUE ~ as.character(x)
  )
}

format_indicator <- function(x) {
  ifelse(x %in% names(INDICATOR_LABELS), INDICATOR_LABELS[x],
         gsub("_", " ", tools::toTitleCase(x)))
}

# ---- WHO Publication Theme ----

who_colors <- list(
  primary = "#2C3E50", secondary = "#3498DB", accent = "#E74C3C",
  success = "#27AE60", warning = "#F39C12", gray = "#95A5A6"
)

theme_who <- function(base_size = 11, base_family = "sans") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      plot.title    = element_text(face = "bold", size = rel(1.2), hjust = 0, margin = margin(b = 5)),
      plot.subtitle = element_text(color = "#555555", size = rel(1.0), hjust = 0, margin = margin(b = 10)),
      plot.caption  = element_text(color = "#777777", size = rel(0.7), hjust = 1, margin = margin(t = 10)),
      axis.title    = element_text(face = "bold", size = rel(0.9)),
      axis.text     = element_text(color = "#333333", size = rel(0.85)),
      panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      axis.line     = element_line(color = "#333333", linewidth = 0.3),
      legend.position      = "top",
      legend.justification = "left",
      legend.title  = element_text(face = "bold", size = rel(0.8)),
      strip.background = element_rect(fill = "#F5F5F5", color = NA),
      strip.text    = element_text(face = "bold", size = rel(0.9), hjust = 0, margin = margin(4, 4, 4, 4))
    )
}

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
#                        VERSION 2.4.0 (NIMBLE)
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
