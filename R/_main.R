#########################################################################################
#
#                    WHO TOBACCO CONTROL PREVALENCE PROJECTION MODEL
#                            _main.R - Orchestrator Script
#
#   Orchestrates the modular pipeline with phased execution:
#     - Phased debugging and testing
#     - Selective re-running of specific pipeline stages
#     - Checkpoint save/restore between phases
#
#   Usage:
#     source("R/_main.R"); run_modular_pipeline()          # Run all phases
#     source("R/_main.R"); run_modular_pipeline("viz")     # Run only visualization
#     source("R/_main.R"); run_modular_pipeline("global")  # Run only global model
#
#########################################################################################

cat("
#########################################################################
#
#              WHO TOBACCO CONTROL PREVALENCE PROJECTION
#                     MODULAR VERSION 2.4.0
#
#########################################################################
")

# ---- Set Working Directory ----
# Ensure we're in the project root
if (!file.exists("data")) {
  stop("Please run this script from the project root directory (where 'data/' exists)")
}

# ============================================================================
# Run modular pipeline (phases can be run independently)
# ============================================================================

run_modular_pipeline <- function(phases = c("prep", "global", "country", "eval", "agg", "viz")) {

  # ---- Phase 1: Configuration and Data Preparation ----
  if ("prep" %in% phases) {
    cat("\n======= PHASE 1: Configuration & Data Preparation =======\n")

    source("R/00_config.R")
    source("R/01_data_prep.R")
    source("R/02_regional_mapping.R")
    source("R/03_splines.R")
    source("R/04_utils.R")
    source("R/05_models_nimble.R")

    # Source diagnostics module for comprehensive MCMC convergence analysis
    if (file.exists("R/mcmc_diagnostics.R")) {
      cat("  Sourcing MCMC diagnostics module...\n")
      source("R/mcmc_diagnostics.R")
      DIAGNOSTICS_AVAILABLE <<- TRUE
    } else {
      cat("  Note: mcmc_diagnostics.R not found. Diagnostics will be limited.\n")
      DIAGNOSTICS_AVAILABLE <<- FALSE
    }

    cat("\nPhase 1 complete. Data prepared.\n")

    # Save checkpoint
    save.image("checkpoints/phase1_prep.RData")
    cat("Checkpoint saved: checkpoints/phase1_prep.RData\n")
  }

  # ---- Phase 2: Global Model Fitting ----
  if ("global" %in% phases) {
    cat("\n======= PHASE 2: Global Model Fitting =======\n")

    # Load checkpoint if needed
    if (!exists("clean_data")) {
      if (file.exists("checkpoints/phase1_prep.RData")) {
        load("checkpoints/phase1_prep.RData")
      } else {
        stop("Phase 1 must be run first. Run: run_modular_pipeline(phases = 'prep')")
      }
    }

    source("R/06_run_global_model.R")

    cat("\nPhase 2 complete. Global models fitted.\n")
    save.image("checkpoints/phase2_global.RData")
  }

  # ---- Phase 3: Country-Specific Model Fitting ----
  if ("country" %in% phases) {
    cat("\n======= PHASE 3: Country-Specific Model Fitting =======\n")

    if (!exists("final_ac_predictions")) {
      if (file.exists("checkpoints/phase2_global.RData")) {
        load("checkpoints/phase2_global.RData")
      } else {
        stop("Phase 2 must be run first.")
      }
    }

    source("R/07_run_country_model.R")

    cat("\nPhase 3 complete. Country models fitted.\n")
    save.image("checkpoints/phase3_country.RData")
  }

  # ---- Phase 4: Model Evaluation and Selection ----
  if ("eval" %in% phases) {
    cat("\n======= PHASE 4: Model Evaluation & Selection =======\n")

    if (!exists("final_predictions_country_specific_ac")) {
      if (file.exists("checkpoints/phase3_country.RData")) {
        load("checkpoints/phase3_country.RData")
      } else {
        stop("Phase 3 must be run first.")
      }
    }

    source("R/08_evaluation.R")

    cat("\nPhase 4 complete. Models evaluated and selected.\n")
    save.image("checkpoints/phase4_eval.RData")
  }

  # ---- Phase 5: Aggregation and Post-Processing ----
  if ("agg" %in% phases) {
    cat("\n======= PHASE 5: Aggregation & Post-Processing =======\n")

    if (!exists("country_sex_summary")) {
      if (file.exists("checkpoints/phase4_eval.RData")) {
        load("checkpoints/phase4_eval.RData")
      } else {
        stop("Phase 4 must be run first.")
      }
    }

    source("R/09_aggregation.R")

    cat("\nPhase 5 complete. Results aggregated.\n")
    save.image("checkpoints/phase5_agg.RData")
  }

  # ---- Phase 6: Visualization ----
  if ("viz" %in% phases) {
    cat("\n======= PHASE 6: Visualization =======\n")

    # Load checkpoint if needed
    if (!exists("final_weighted_results_selected")) {
      if (file.exists("checkpoints/phase5_agg.RData")) {
        load("checkpoints/phase5_agg.RData")
      } else {
        stop("Phase 5 must be run first.")
      }
    }

    source("R/10_visualization.R")

    # Generate all basic/exploratory plots
    generate_all_visualizations()

    # ---- Publication Tables ----
    cat("\n  Generating publication tables...\n")
    tryCatch({
      source("R/11_publication_tables.R")
    }, error = function(e) {
      cat(sprintf("  WARNING: Publication tables failed: %s\n", e$message))
      cat("  Continuing with pipeline...\n")
    })

    # ---- Publication Figures ----
    cat("\n  Generating publication figures...\n")
    tryCatch({
      if (file.exists("R/12_publication_figures.R")) {
        source("R/12_publication_figures.R")
        # Explicit call — the file's auto-execution depends on final_ac_predictions
        # which may not be in the checkpoint. Call directly with available variables.
        generate_all_figures(
          clean_data = clean_data,
          global_predictions = if (exists("final_ac_predictions")) final_ac_predictions else NULL,
          country_predictions = if (exists("final_predictions_country_specific_ac")) final_predictions_country_specific_ac else NULL,
          weighted_results = if (exists("final_weighted_results_selected")) final_weighted_results_selected else NULL,
          country_region_mapping = if (exists("country_region_mapping")) country_region_mapping else NULL,
          country_name_mapping = if (exists("country_name_mapping")) country_name_mapping else NULL
        )
      } else {
        cat("  WARNING: R/12_publication_figures.R not found.\n")
      }
    }, error = function(e) {
      cat(sprintf("  WARNING: Publication figures failed: %s\n", e$message))
      cat("  Continuing with pipeline...\n")
    })

    # ---- Strata-Level Figures ----
    cat("\n  Generating strata-level figures...\n")
    tryCatch({
      if (file.exists("R/13_strata_figures.R")) {
        source("R/13_strata_figures.R")
        # Explicit call — the file never auto-executes, only defines functions
        generate_all_strata_figures(
          clean_data = clean_data,
          predictions = if (exists("final_ac_predictions")) final_ac_predictions else NULL,
          trend_data = if (exists("final_weighted_results_selected")) final_weighted_results_selected else NULL,
          country_name_mapping = if (exists("country_name_mapping")) country_name_mapping else NULL
        )
      } else {
        cat("  WARNING: R/13_strata_figures.R not found.\n")
      }
    }, error = function(e) {
      cat(sprintf("  WARNING: Strata figures failed: %s\n", e$message))
      cat("  Continuing with pipeline...\n")
    })

    cat("\nPhase 6 complete. Visualizations generated.\n")
  }

  # ---- Completion ----
  cat("
#########################################################################
#
#              PIPELINE COMPLETE
#
#  Outputs saved to:
#    - diagnostics/logs/     MCMC convergence log files
#    - diagnostics/tables/   R-hat and ESS master tables
#    - results/              Model predictions and weighted results
#    - evaluation/           Model comparison metrics
#    - outputs/figures/      All visualizations
#    - processing/           Intermediate sample files
#    - checkpoints/          RData files for each phase
#    - figures_publication/  Publication-quality figures
#    - tables_publication/   Publication-quality tables
#
#########################################################################
")

  # Notification sound (if beepr is loaded)
  try(beepr::beep(3), silent = TRUE)
}

# ============================================================================
# Create checkpoints directory
# ============================================================================

dir.create("checkpoints", showWarnings = FALSE)

# ============================================================================
# Print usage instructions
# ============================================================================

cat("
Pipeline loaded. Available commands:

  run_modular_pipeline()            - Run all phases
  run_modular_pipeline('prep')      - Run only Phase 1 (data preparation)
  run_modular_pipeline('global')    - Run only Phase 2 (global model)
  run_modular_pipeline('country')   - Run only Phase 3 (country models)
  run_modular_pipeline('eval')      - Run only Phase 4 (evaluation)
  run_modular_pipeline('agg')       - Run only Phase 5 (aggregation)
  run_modular_pipeline('viz')       - Run only Phase 6 (visualization)

Individual module files:
  R/00_config.R            - Configuration and constants
  R/01_data_prep.R         - Data loading and cleaning
  R/02_regional_mapping.R  - Regional classification
  R/03_splines.R           - Spline basis construction
  R/04_utils.R             - Utility functions
  R/05_models_nimble.R     - NIMBLE model definitions
  R/06_run_global_model.R  - Global model fitting (complete with no-data countries)
  R/07_run_country_model.R - Country model fitting (complete with OPT-1 to OPT-6)
  R/08_evaluation.R        - Model evaluation
  R/09_aggregation.R       - Regional aggregation
  R/10_visualization.R     - Plotting functions
  R/11_publication_tables.R - Publication table generation
  R/12_publication_figures.R - Publication figure generation
  R/mcmc_diagnostics.R     - MCMC convergence diagnostics

Diagnostics functions (after sourcing mcmc_diagnostics.R):
  generate_mcmc_diagnostics()        - Full diagnostics for one model
  generate_worst_parameters_report() - Find worst parameters across models
  generate_all_models_summary()      - Summary across all models
  generate_convergence_gt_table()    - Elegant GT table for publication
  clear_diagnostic_tables()          - Reset tables for fresh run

")
