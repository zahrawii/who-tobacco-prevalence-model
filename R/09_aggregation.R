#########################################################################################
#
#                    WHO TOBACCO CONTROL PREVALENCE PROJECTION MODEL
#                       09_aggregation.R - Regional Aggregation
#
#   Contains: Section 15 from pipeline_monolith.R
#     - Regional population weight preparation
#     - Sample-based regional aggregation
#     - Regional metrics and achievement calculation
#     - High-certainty country percentages
#     - Publication tables (GT)
#
#   Requires: All previous modules (00-08) to be sourced
#   Outputs: regional_results_df, regional aggregation files, GT tables
#
#########################################################################################

cat("\n")
cat("================================================================\n")
cat("  STARTING POST-PROCESSING & REGIONAL AGGREGATION PIPELINE\n")
cat("================================================================\n")

# ---- 15.1 Load Model Selection Data ----
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

  # 3. Identify available definition types from one country
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
        if(length(model_choice) == 0) model_choice <- "Global"

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
      n_samples_out <- 1000
      n_ages <- nrow(country_samples_matrix_list[[1]])
      regional_mat <- matrix(0, nrow = n_ages, ncol = n_samples_out)

      for (s in 1:n_samples_out) {
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
        regional_prob <- pmin(pmax(regional_prob, 0.0001), 0.9999)
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

calc_regional_metrics <- function() {
  results_list <- list()

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

country_achievement_summary <- final_weighted_results_selected %>%
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

# Indicator Labels for Display
labels_map <- c(
  "daily_user_cigarettes" = "Daily Cigarettes",
  "current_user_cigarettes" = "Current Cigarettes",
  "daily_user_any_smoked_tobacco" = "Daily Any Smoked",
  "current_user_any_smoked_tobacco" = "Current Any Smoked",
  "daily_user_any_tobacco_product" = "Daily Any Tobacco",
  "current_user_any_tobacco_product" = "Current Any Tobacco"
)

# ---- Table 1: Comprehensive Country Performance ----

country_table_data <- final_weighted_results_selected %>%
  filter(Year == TARGET_YEAR) %>%
  left_join(
    final_weighted_results_selected %>% filter(Year == BASE_YEAR) %>% select(Country, Sex, Def_Type_Code, Base_Prev = weighted_mean),
    by = c("Country", "Sex", "Def_Type_Code")
  ) %>%
  mutate(
    Country_Name = tools::toTitleCase(country_name_mapping[Country]),
    Gender = case_when(tolower(Sex) == "males" ~ "Men", tolower(Sex) == "females" ~ "Women", TRUE ~ Sex),
    Indicator = labels_map[Def_Type_Code],
    Base_Fmt = sprintf("%.1f%%", Base_Prev * 100),
    Proj_Fmt = sprintf("%.1f%%", weighted_mean * 100),
    Prob_Fmt = sprintf("%.0f%%", prob_achieving_target * 100),
    Model_Short = ifelse(Model_Type == "Global", "Global", "Country")
  ) %>%
  filter(!is.na(Indicator)) %>%
  select(Country_Name, Gender, Indicator, Model_Short, Base_Fmt, Proj_Fmt, Prob_Fmt, prob_achieving_target) %>%
  pivot_wider(
    names_from = Indicator,
    values_from = c(Base_Fmt, Proj_Fmt, Prob_Fmt, prob_achieving_target),
    names_sep = "_"
  )

# Generate GT Table
tryCatch({
  gt_country <- country_table_data %>%
    select(Country_Name, Gender, Model_Short,
           contains("Daily Cigarettes"), contains("Current Cigarettes")) %>%
    gt() %>%
    tab_header(
      title = md("**WHO Tobacco Control Target Achievement**"),
      subtitle = md(paste0("Country-level projections for ", TARGET_YEAR, " based on RMSE-selected AC models"))
    ) %>%
    cols_label(
      Country_Name = "Country", Gender = "Gender", Model_Short = "Model",
      `Base_Fmt_Daily Cigarettes` = "2010 Base", `Proj_Fmt_Daily Cigarettes` = "2025 Proj", `Prob_Fmt_Daily Cigarettes` = "Prob.",
      `Base_Fmt_Current Cigarettes` = "2010 Base", `Proj_Fmt_Current Cigarettes` = "2025 Proj", `Prob_Fmt_Current Cigarettes` = "Prob."
    ) %>%
    tab_spanner(label = "Daily Cigarettes", columns = contains("Daily Cigarettes")) %>%
    tab_spanner(label = "Current Cigarettes", columns = contains("Current Cigarettes")) %>%
    cols_hide(columns = contains("prob_achieving_target")) %>%
    tab_options(table.font.size = 10, data_row.padding = 3)

  gtsave(gt_country, "tables_publication/Table1_Country_Performance.html")
  cat("    Table 1 saved.\n")
}, error = function(e) cat("    Warning: Could not generate Table 1 -", e$message, "\n"))

# ---- Table 2: Regional Summary (Lancet Style) ----

tryCatch({
  regional_table_data <- regional_full_summary %>%
    mutate(
      Region_Name = tools::toTitleCase(Region),
      Gender = case_when(tolower(Sex) == "males" ~ "Men", tolower(Sex) == "females" ~ "Women", TRUE ~ Sex),
      Indicator = labels_map[Def_Type_Code],
      Prev_Fmt = sprintf("%.1f%%", Projected_Prevalence * 100),
      Achieve_Count_Fmt = sprintf("%d / %d", Countries_High_Certainty, Total_Countries),
      Achieve_Pct_Fmt = sprintf("%.0f%%", High_Certainty_Pct)
    ) %>%
    filter(!is.na(Indicator)) %>%
    select(Region_Name, Gender, Indicator, Prev_Fmt, Achieve_Count_Fmt, Achieve_Pct_Fmt) %>%
    pivot_wider(
      names_from = Indicator,
      values_from = c(Prev_Fmt, Achieve_Count_Fmt, Achieve_Pct_Fmt),
      names_sep = "_"
    )

  gt_regional <- regional_table_data %>%
    select(Region_Name, Gender, contains("Daily Cigarettes")) %>%
    gt() %>%
    tab_header(
      title = md("**Regional Tobacco Control Status**"),
      subtitle = md(paste0("Population-weighted regional prevalence and country-level success rates (", TARGET_YEAR, ")"))
    ) %>%
    cols_label(
      Region_Name = "Region", Gender = "Gender",
      `Prev_Fmt_Daily Cigarettes` = "Regional Prevalence",
      `Achieve_Count_Fmt_Daily Cigarettes` = "Successful Countries (n/N)",
      `Achieve_Pct_Fmt_Daily Cigarettes` = "% Countries with High Certainty"
    ) %>%
    tab_footnote(
      footnote = "High Certainty defined as >60% probability of achieving the 30% reduction target.",
      locations = cells_column_labels(columns = `Achieve_Pct_Fmt_Daily Cigarettes`)
    ) %>%
    tab_options(heading.background.color = "#F0F0F0")

  gtsave(gt_regional, "tables_publication/Table2_Regional_Summary.html")
  cat("    Table 2 saved.\n")
}, error = function(e) cat("    Warning: Could not generate Table 2 -", e$message, "\n"))

cat("\nRegional aggregation module complete.\n")
