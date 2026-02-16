#########################################################################################
#
#                    WHO TOBACCO PIPELINE - MODULE 11: PUBLICATION TABLES
#
#   This module generates publication-quality tables following PhD thesis guidelines.
#
#   MAIN MANUSCRIPT TABLES (Current Smoking Only):
#     - Table 1: Baseline Characteristics and Prevalence Trends
#     - Table 2: WHO 30% Reduction Target Achievement
#     - Table 3: Tobacco Endgame Achievement by 2040
#
#   SUPPLEMENTARY TABLES (All Indicators):
#     - Table S1: Regional WHO Target Achievement by Indicator (WIDE)
#     - Table S2: Regional Endgame Achievement by Indicator (WIDE)
#     - Table S3: Country-Level Target Achievement by Stratum (LONG)
#     - Table S4: Model Selection and RMSE by Stratum (LONG)
#     - Table S5: Data Sources by Stratum - FULL survey names (LONG)
#
#   KEY DESIGN PRINCIPLES:
#     - COUNT COUNTRIES, never population-weighted aggregation
#     - Row spanners for gender: Men first, Women second
#     - Each supplementary row = gender x country x indicator stratum
#
#########################################################################################

cat("\n")
cat("###########################################################################\n")
cat("#                                                                         #\n")
cat("#        MODULE 11: PUBLICATION TABLES                                    #\n")
cat("#                                                                         #\n")
cat("###########################################################################\n\n")

# ============================================================================
# CONFIGURATION
# ============================================================================

# Year parameters
BASELINE_YEAR <- 2010
INTERIM_YEAR <- 2025
TARGET_YEAR <- 2030
ENDGAME_YEAR <- 2040

# Target thresholds
REDUCTION_TARGET <- 0.30
ENDGAME_THRESHOLD <- 0.05
NEAR_ENDGAME_UPPER <- 0.10
VIRTUAL_ELIMINATION <- 0.02
ON_TRACK_MARGIN <- 0.10

# Indicator labels
INDICATOR_LABELS <- c(
  "current_user_any_tobacco_product" = "Current Any Tobacco",
  "current_user_any_smoked_tobacco" = "Current Smoked",
  "current_user_cigarettes" = "Current Cigarettes",
  "daily_user_any_tobacco_product" = "Daily Any Tobacco",
  "daily_user_any_smoked_tobacco" = "Daily Smoked",
  "daily_user_cigarettes" = "Daily Cigarettes"
)

# Short labels for column headers
INDICATOR_SHORT <- c(
  "current_user_any_tobacco_product" = "CUR_ANY_TOB",
  "current_user_any_smoked_tobacco" = "CUR_SMOKED",
  "current_user_cigarettes" = "CUR_CIG",
  "daily_user_any_tobacco_product" = "DAILY_ANY_TOB",
  "daily_user_any_smoked_tobacco" = "DAILY_SMOKED",
  "daily_user_cigarettes" = "DAILY_CIG"
)

# Primary indicator for main tables
PRIMARY_INDICATOR <- "current_user_any_tobacco_product"

# Color palette
COLORS <- list(
  excellent = "#1B5E20",
  very_good = "#4CAF50",
  good = "#A5D6A7",
  moderate = "#FFF59D",
  low = "#FFCDD2",
  none = "#E57373",
  row_group = "#FFEBEE",
  subtotal = "#E0E0E0",
  header = "#F5F5F5",
  positive_change = "#2E7D32",
  negative_change = "#C62828"
)

# ============================================================================
# MAIN FUNCTION: Generate All Publication Tables
# ============================================================================

#' Generate All Publication Tables
#'
#' @param clean_data The cleaned survey data
#' @param country_region_mapping Country to region mapping dataframe
#' @param model_selection Optional model selection results
#' @param output_dir Output directory for tables (default: "tables_publication")
#' @return List with generation status for each table
generate_publication_tables <- function(clean_data,
                                         country_region_mapping,
                                         model_selection = NULL,
                                         output_dir = "tables_publication") {

  # Check dependencies
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required for publication tables. Install with: install.packages('gt')")
  }

  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(gt)
  })

  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  results <- list()

  cat("\n========== GENERATING PUBLICATION TABLES ==========\n")

  # Source the main publication_tables.R which contains all table functions
  # We do this to avoid code duplication

  # Get the directory of this script
  script_dir <- dirname(sys.frame(1)$ofile)
  if (is.null(script_dir)) script_dir <- "R"

  main_tables_file <- file.path(script_dir, "publication_tables.R")

  if (file.exists(main_tables_file)) {
    cat(sprintf("  Sourcing table functions from: %s\n", main_tables_file))

    # Source only the function definitions (not the execution block)
    # We'll execute the table generation here with our parameters

    # For now, we'll execute the full file which handles everything
    source(main_tables_file, local = TRUE)

    results$status <- "complete"
    results$output_dir <- output_dir

  } else {
    cat("  WARNING: publication_tables.R not found. Using inline definitions.\n")

    # Inline minimal table generation would go here
    # For production, always use publication_tables.R

    results$status <- "failed"
    results$error <- "publication_tables.R not found"
  }

  return(results)
}

# ============================================================================
# STANDALONE EXECUTION
# ============================================================================

# When sourced as part of the pipeline, check if data is available
if (exists("clean_data") && exists("country_region_mapping")) {

  cat("\n  Data available - delegating to publication_tables.R\n")

  # The full table generation is in publication_tables.R
  # Source it directly to execute all table generation
  if (file.exists("R/publication_tables.R")) {
    source("R/publication_tables.R")
  } else if (file.exists("publication_tables.R")) {
    source("publication_tables.R")
  } else {
    cat("  ERROR: Cannot find publication_tables.R\n")
  }

} else {
  cat("\n  No data available. Tables will be generated when clean_data and\n")
  cat("  country_region_mapping are loaded.\n")
  cat("\n  To generate tables, either:\n")
  cat("    1. Run the full pipeline first, OR\n")
  cat("    2. Call generate_publication_tables(clean_data, country_region_mapping)\n")
}

cat("\n  Module 11 loaded.\n")
