#########################################################################################
#
#                    WHO TOBACCO PIPELINE - PUBLICATION TABLES
#
#   This file contains the actual table generation code for publication-quality tables.
#   Sourced by 11_publication_tables.R module.
#
#   MAIN MANUSCRIPT TABLES (Current Smoking Only):
#     - Table 1: Baseline Characteristics and Prevalence Trends
#     - Table 2: WHO 30% Reduction Target Achievement (2025 & 2030)
#     - Table 3: Tobacco Endgame Achievement by 2040 (<5% Prevalence)
#
#   SUPPLEMENTARY TABLES (All Indicators):
#     - Table S1: Regional WHO Target Achievement by Indicator (WIDE)
#     - Table S2: Regional Endgame Achievement by Indicator (WIDE)
#     - Table S3: Country-Level Target Achievement by Stratum (LONG)
#     - Table S4: Model Selection and RMSE by Stratum (LONG)
#     - Table S5: Data Sources by Stratum (LONG)
#
#   KEY DESIGN PRINCIPLES:
#     - COUNT COUNTRIES, never population-weighted aggregation
#     - Row spanners for gender: Men first, Women second
#     - Use "Men" and "Women" (NOT "Males" / "Females")
#     - Use "Gender" (NOT "Sex") in column headers
#     - Each supplementary row = gender x country x indicator stratum
#
#########################################################################################

cat("\n========== GENERATING PUBLICATION TABLES ==========\n")

# ============================================================================
# LOAD DEPENDENCIES
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(gt)
})

# ============================================================================
# CONFIGURATION
# ============================================================================

# Year parameters
BASELINE_YEAR <- 2010
COMPARISON_YEAR <- 2024
INTERIM_YEAR <- 2025
TARGET_YEAR <- 2030
ENDGAME_YEAR <- 2040

# Target thresholds
REDUCTION_TARGET <- 0.30          # 30% relative reduction
ENDGAME_THRESHOLD <- 0.05         # <5% = endgame
NEAR_ENDGAME_UPPER <- 0.10        # 5-10% = near-endgame
VIRTUAL_ELIMINATION <- 0.02       # <2% = virtual elimination
ON_TRACK_MARGIN <- 0.10           # Within 10% of target = on track

# Primary indicator (WHO FCTC priority)
PRIMARY_INDICATOR <- "current_user_any_tobacco_product"

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
# HELPER FUNCTIONS
# ============================================================================

#' Format gender for display (Men/Women, NOT Males/Females)
format_gender <- function(sex) {
  case_when(
    tolower(sex) %in% c("males", "male") ~ "Men",
    tolower(sex) %in% c("females", "female") ~ "Women",
    TRUE ~ as.character(sex)
  )
}

#' Extract percentage from formatted string "n (X%)"
extract_pct <- function(x) {
  as.numeric(str_extract(x, "(?<=\\()\\d+(?=%)"))
}

#' Format count with percentage
format_count_pct <- function(n, total) {
  pct <- round(n / total * 100)
  sprintf("%d (%d%%)", n, pct)
}

#' Format prevalence with range
format_prev_range <- function(median_val, min_val, max_val) {
  sprintf("%.1f (%.1f-%.1f)", median_val, min_val, max_val)
}

#' Common table options for publication quality
get_table_options <- function() {
  list(
    table.width = pct(100),
    table.font.size = px(10),
    table.font.names = c("Arial", "Helvetica", "sans-serif"),
    heading.title.font.size = px(14),
    heading.title.font.weight = "bold",
    heading.subtitle.font.size = px(11),
    heading.background.color = "#FFFFFF",
    heading.border.bottom.color = "#BDBDBD",
    heading.border.bottom.width = px(2),
    column_labels.background.color = "#F5F5F5",
    column_labels.font.size = px(9),
    column_labels.font.weight = "bold",
    column_labels.border.top.color = "#BDBDBD",
    column_labels.border.bottom.color = "#BDBDBD",
    row_group.background.color = "#FFEBEE",
    row_group.font.size = px(11),
    row_group.font.weight = "bold",
    row_group.border.top.color = "#BDBDBD",
    row_group.border.bottom.color = "#BDBDBD",
    stub.border.color = "#BDBDBD",
    data_row.padding = px(6),
    summary_row.background.color = "#EEEEEE",
    summary_row.border.color = "#BDBDBD",
    grand_summary_row.background.color = "#E0E0E0",
    footnotes.font.size = px(9),
    footnotes.padding = px(4),
    source_notes.font.size = px(9)
  )
}

# ============================================================================
# TABLE 1: BASELINE CHARACTERISTICS AND PREVALENCE TRENDS
# ============================================================================

generate_table1 <- function(clean_data,
                            weighted_results,
                            country_region_mapping,
                            output_dir = "tables_publication") {

  cat("  Generating Table 1: Baseline Characteristics and Prevalence Trends...\n")

  # Detect column names
  country_col <- if ("Country" %in% names(weighted_results)) "Country" else "wb_country_abv"
  sex_col <- if ("Sex" %in% names(weighted_results)) "Sex" else "sex"
  indicator_col <- if ("Def_Type_Code" %in% names(weighted_results)) "Def_Type_Code" else "def_type_code"
  year_col <- if ("Year" %in% names(weighted_results)) "Year" else "year"
  prev_col <- if ("weighted_mean" %in% names(weighted_results)) "weighted_mean" else "prevalence"

  # STEP 1: Get prevalence at baseline and comparison years
  prevalence_comparison <- weighted_results %>%
    filter(!!sym(indicator_col) == PRIMARY_INDICATOR) %>%
    filter(!!sym(year_col) %in% c(BASELINE_YEAR, COMPARISON_YEAR)) %>%
    select(!!sym(country_col), !!sym(sex_col), !!sym(year_col), !!sym(prev_col)) %>%
    rename(Country = !!sym(country_col), Sex = !!sym(sex_col),
           Year = !!sym(year_col), Prevalence = !!sym(prev_col)) %>%
    pivot_wider(
      names_from = Year,
      values_from = Prevalence,
      names_prefix = "prev_"
    )

  # Handle column names dynamically
  baseline_col <- paste0("prev_", BASELINE_YEAR)
  comparison_col <- paste0("prev_", COMPARISON_YEAR)

  prevalence_comparison <- prevalence_comparison %>%
    mutate(
      trend_direction = case_when(
        !!sym(comparison_col) < !!sym(baseline_col) ~ "Decrease",
        !!sym(comparison_col) > !!sym(baseline_col) ~ "Increase",
        TRUE ~ "No Change"
      )
    )

  # STEP 2: Join with region mapping
  prevalence_with_region <- prevalence_comparison %>%
    left_join(
      country_region_mapping %>% select(wb_country_abv, region_consolidated),
      by = c("Country" = "wb_country_abv")
    ) %>%
    filter(!is.na(region_consolidated))

  # STEP 3: Aggregate by region - COUNT COUNTRIES
  table1_trends <- prevalence_with_region %>%
    group_by(region_consolidated, Sex) %>%
    summarise(
      n_countries = n(),
      n_decrease = sum(trend_direction == "Decrease", na.rm = TRUE),
      n_increase = sum(trend_direction == "Increase", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      decrease_fmt = format_count_pct(n_decrease, n_countries),
      increase_fmt = format_count_pct(n_increase, n_countries)
    )

  # STEP 4: Get survey characteristics
  survey_characteristics <- clean_data %>%
    filter(def_type_code == PRIMARY_INDICATOR) %>%
    left_join(
      country_region_mapping %>% select(wb_country_abv, region_consolidated),
      by = c("wb_country_abv")
    ) %>%
    group_by(region_consolidated, sex) %>%
    summarise(
      n_surveys = n_distinct(survey),
      n_observations = n(),
      year_min = min(year, na.rm = TRUE),
      year_max = max(year, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      year_range = paste0(year_min, "-", year_max)
    )

  # STEP 5: Get baseline prevalence (median with range across countries)
  baseline_prevalence <- prevalence_with_region %>%
    group_by(region_consolidated, Sex) %>%
    summarise(
      median_baseline = median(!!sym(baseline_col) * 100, na.rm = TRUE),
      min_baseline = min(!!sym(baseline_col) * 100, na.rm = TRUE),
      max_baseline = max(!!sym(baseline_col) * 100, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      baseline_fmt = sprintf("%.1f (%.1f-%.1f)",
                             round(median_baseline, 1),
                             round(min_baseline, 1),
                             round(max_baseline, 1))
    )

  # STEP 6: Combine all data
  table1_data <- table1_trends %>%
    left_join(survey_characteristics,
              by = c("region_consolidated", "Sex" = "sex")) %>%
    left_join(baseline_prevalence,
              by = c("region_consolidated", "Sex")) %>%
    mutate(
      Region_Name = tools::toTitleCase(region_consolidated),
      Gender = format_gender(Sex)
    ) %>%
    arrange(desc(Gender == "Men"), Region_Name) %>%
    select(Region_Name, Gender, n_countries, n_surveys, n_observations,
           year_range, baseline_fmt, decrease_fmt, increase_fmt)

  # STEP 7: Calculate subtotals
  subtotals <- table1_data %>%
    group_by(Gender) %>%
    summarise(
      Region_Name = "Subtotal",
      n_countries = sum(n_countries, na.rm = TRUE),
      n_surveys = sum(n_surveys, na.rm = TRUE),
      n_observations = sum(n_observations, na.rm = TRUE),
      year_range = "-",
      baseline_fmt = "-",
      n_decrease = sum(as.numeric(str_extract(decrease_fmt, "^\\d+")), na.rm = TRUE),
      n_increase = sum(as.numeric(str_extract(increase_fmt, "^\\d+")), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      decrease_fmt = format_count_pct(n_decrease, n_countries),
      increase_fmt = format_count_pct(n_increase, n_countries)
    ) %>%
    select(-n_decrease, -n_increase)

  # Combine with subtotals
  table1_final <- bind_rows(
    table1_data %>% filter(Gender == "Men"),
    subtotals %>% filter(Gender == "Men"),
    table1_data %>% filter(Gender == "Women"),
    subtotals %>% filter(Gender == "Women")
  )

  # STEP 8: Create GT table
  gt_table1 <- table1_final %>%
    select(-Gender) %>%
    gt(rowname_col = "Region_Name", groupname_col = NULL) %>%
    tab_header(
      title = md("**Table 1: Baseline Characteristics and Prevalence Trends by WHO Region**"),
      subtitle = md("Current Smoking Prevalence Among Adults Aged 15+ Years")
    ) %>%
    tab_row_group(
      label = md("**Men**"),
      rows = 1:7
    ) %>%
    tab_row_group(
      label = md("**Women**"),
      rows = 8:14
    ) %>%
    row_group_order(groups = c("Men", "Women")) %>%
    cols_label(
      n_countries = "Countries (n)",
      n_surveys = "Surveys (n)",
      n_observations = "Observations (n)",
      year_range = "Survey Years",
      baseline_fmt = md("2010 Baseline<br>Median (Range)"),
      decrease_fmt = "Decrease",
      increase_fmt = "Increase"
    ) %>%
    tab_spanner(
      label = md("**Direction of Trend, 2010-2024**"),
      columns = c(decrease_fmt, increase_fmt)
    ) %>%
    fmt_number(
      columns = n_observations,
      use_seps = TRUE,
      decimals = 0
    ) %>%
    tab_style(
      style = list(
        cell_fill(color = COLORS$row_group),
        cell_text(weight = "bold", size = px(11))
      ),
      locations = cells_row_groups()
    ) %>%
    tab_style(
      style = list(
        cell_fill(color = COLORS$subtotal),
        cell_text(weight = "bold"),
        cell_borders(sides = "top", color = "#BDBDBD", weight = px(1))
      ),
      locations = cells_body(rows = Region_Name == "Subtotal")
    ) %>%
    tab_footnote(
      footnote = md("Current smoking defined as any tobacco product use in the past 30 days (WHO FCTC indicator)."),
      locations = cells_title(groups = "subtitle")
    ) %>%
    tab_footnote(
      footnote = md("Baseline prevalence shows median across countries in region with range (minimum-maximum). NOT population-weighted."),
      locations = cells_column_labels(columns = baseline_fmt)
    ) %>%
    tab_footnote(
      footnote = md("Direction of trend determined by comparing model-predicted 2024 prevalence to 2010 baseline for each country."),
      locations = cells_column_spanners(spanners = "Direction of Trend, 2010-2024")
    ) %>%
    tab_source_note(
      source_note = md("*Data source: WHO Global Tobacco Surveillance System (GTSS). Predictions from Bayesian Age-Period-Cohort models.*")
    ) %>%
    tab_options(
      table.width = pct(100),
      table.font.size = px(10),
      heading.title.font.size = px(14),
      heading.title.font.weight = "bold",
      column_labels.background.color = COLORS$header,
      column_labels.font.weight = "bold",
      row_group.background.color = COLORS$row_group,
      data_row.padding = px(6)
    )

  # Save table
  tryCatch({
    gtsave(gt_table1, file.path(output_dir, "Table1_Baseline_Trends.html"))
    cat("    Table 1 saved to", file.path(output_dir, "Table1_Baseline_Trends.html"), "\n")
  }, error = function(e) {
    cat("    Warning: Could not save Table 1 -", e$message, "\n")
  })

  return(gt_table1)
}

# ============================================================================
# TABLE 2: WHO 30% REDUCTION TARGET ACHIEVEMENT
# ============================================================================

generate_table2 <- function(weighted_results,
                            country_region_mapping,
                            output_dir = "tables_publication") {

  cat("  Generating Table 2: WHO 30% Reduction Target Achievement...\n")

  # Detect column names
  country_col <- if ("Country" %in% names(weighted_results)) "Country" else "wb_country_abv"
  sex_col <- if ("Sex" %in% names(weighted_results)) "Sex" else "sex"
  indicator_col <- if ("Def_Type_Code" %in% names(weighted_results)) "Def_Type_Code" else "def_type_code"
  year_col <- if ("Year" %in% names(weighted_results)) "Year" else "year"
  prev_col <- if ("weighted_mean" %in% names(weighted_results)) "weighted_mean" else "prevalence"

  # STEP 1: Get baseline values
  baseline_data <- weighted_results %>%
    filter(!!sym(indicator_col) == PRIMARY_INDICATOR,
           !!sym(year_col) == BASELINE_YEAR) %>%
    select(!!sym(country_col), !!sym(sex_col), !!sym(prev_col)) %>%
    rename(Country = !!sym(country_col), Sex = !!sym(sex_col),
           BaseYearPrevalence = !!sym(prev_col)) %>%
    mutate(
      Target_Threshold = BaseYearPrevalence * (1 - REDUCTION_TARGET)
    )

  # STEP 2: Get predictions at assessment years
  predictions <- weighted_results %>%
    filter(!!sym(indicator_col) == PRIMARY_INDICATOR,
           !!sym(year_col) %in% c(INTERIM_YEAR, TARGET_YEAR)) %>%
    select(!!sym(country_col), !!sym(sex_col), !!sym(year_col), !!sym(prev_col)) %>%
    rename(Country = !!sym(country_col), Sex = !!sym(sex_col),
           Year = !!sym(year_col), Prevalence = !!sym(prev_col)) %>%
    pivot_wider(
      names_from = Year,
      values_from = Prevalence,
      names_prefix = "pred_"
    )

  # Column names
  pred_interim_col <- paste0("pred_", INTERIM_YEAR)
  pred_target_col <- paste0("pred_", TARGET_YEAR)

  # STEP 3: Determine achievement status
  achievement_status <- baseline_data %>%
    left_join(predictions, by = c("Country", "Sex")) %>%
    left_join(
      country_region_mapping %>% select(wb_country_abv, region_consolidated),
      by = c("Country" = "wb_country_abv")
    ) %>%
    filter(!is.na(region_consolidated)) %>%
    mutate(
      Achieved_2025 = !!sym(pred_interim_col) <= Target_Threshold,
      Achieved_2030 = !!sym(pred_target_col) <= Target_Threshold,
      OnTrack_2025 = !!sym(pred_interim_col) <= Target_Threshold * (1 + ON_TRACK_MARGIN)
    )

  # STEP 4: Aggregate by region - COUNT COUNTRIES
  table2_data <- achievement_status %>%
    group_by(region_consolidated, Sex) %>%
    summarise(
      n_countries = n(),
      achieved_2025 = sum(Achieved_2025, na.rm = TRUE),
      ontrack_2025 = sum(OnTrack_2025, na.rm = TRUE),
      achieved_2030 = sum(Achieved_2030, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      additional_achievers = achieved_2030 - achieved_2025,
      achieved_2025_fmt = format_count_pct(achieved_2025, n_countries),
      ontrack_2025_fmt = format_count_pct(ontrack_2025, n_countries),
      achieved_2030_fmt = format_count_pct(achieved_2030, n_countries),
      change_fmt = sprintf("%+d", additional_achievers),
      Region_Name = tools::toTitleCase(region_consolidated),
      Gender = format_gender(Sex)
    ) %>%
    arrange(desc(Gender == "Men"), Region_Name)

  # STEP 5: Calculate subtotals
  subtotals <- table2_data %>%
    group_by(Gender) %>%
    summarise(
      Region_Name = "Subtotal",
      n_countries = sum(n_countries),
      achieved_2025 = sum(achieved_2025),
      ontrack_2025 = sum(ontrack_2025),
      achieved_2030 = sum(achieved_2030),
      additional_achievers = sum(additional_achievers),
      .groups = "drop"
    ) %>%
    mutate(
      achieved_2025_fmt = format_count_pct(achieved_2025, n_countries),
      ontrack_2025_fmt = format_count_pct(ontrack_2025, n_countries),
      achieved_2030_fmt = format_count_pct(achieved_2030, n_countries),
      change_fmt = sprintf("%+d", additional_achievers)
    )

  # Combine
  table2_final <- bind_rows(
    table2_data %>% filter(Gender == "Men") %>%
      select(Region_Name, Gender, n_countries, achieved_2025_fmt, ontrack_2025_fmt,
             achieved_2030_fmt, change_fmt),
    subtotals %>% filter(Gender == "Men") %>%
      select(Region_Name, Gender, n_countries, achieved_2025_fmt, ontrack_2025_fmt,
             achieved_2030_fmt, change_fmt),
    table2_data %>% filter(Gender == "Women") %>%
      select(Region_Name, Gender, n_countries, achieved_2025_fmt, ontrack_2025_fmt,
             achieved_2030_fmt, change_fmt),
    subtotals %>% filter(Gender == "Women") %>%
      select(Region_Name, Gender, n_countries, achieved_2025_fmt, ontrack_2025_fmt,
             achieved_2030_fmt, change_fmt)
  )

  # STEP 6: Create GT table
  n_men_rows <- sum(table2_final$Gender == "Men")
  n_women_rows <- sum(table2_final$Gender == "Women")

  gt_table2 <- table2_final %>%
    select(-Gender) %>%
    gt(rowname_col = "Region_Name") %>%
    tab_header(
      title = md("**Table 2: WHO 30% Reduction Target Achievement by Region**"),
      subtitle = md("Countries Achieving >=30% Reduction in Current Smoking Prevalence")
    ) %>%
    tab_row_group(
      label = md("**Men**"),
      rows = 1:n_men_rows
    ) %>%
    tab_row_group(
      label = md("**Women**"),
      rows = (n_men_rows + 1):(n_men_rows + n_women_rows)
    ) %>%
    row_group_order(groups = c("Men", "Women")) %>%
    cols_label(
      n_countries = "Countries (n)",
      achieved_2025_fmt = "Achieved",
      ontrack_2025_fmt = "On Track",
      achieved_2030_fmt = "Achieved",
      change_fmt = "Change from 2025"
    ) %>%
    tab_spanner(
      label = md("**2025 Interim Assessment**"),
      columns = c(achieved_2025_fmt, ontrack_2025_fmt)
    ) %>%
    tab_spanner(
      label = md("**2030 Target Assessment**"),
      columns = c(achieved_2030_fmt, change_fmt)
    ) %>%
    tab_style(
      style = list(
        cell_fill(color = COLORS$row_group),
        cell_text(weight = "bold")
      ),
      locations = cells_row_groups()
    ) %>%
    tab_style(
      style = list(
        cell_fill(color = COLORS$subtotal),
        cell_text(weight = "bold"),
        cell_borders(sides = "top", color = "#BDBDBD", weight = px(2))
      ),
      locations = cells_body(rows = Region_Name == "Subtotal")
    ) %>%
    tab_style(
      style = cell_text(color = COLORS$positive_change, weight = "bold"),
      locations = cells_body(
        columns = change_fmt,
        rows = str_detect(change_fmt, "^\\+[1-9]")
      )
    ) %>%
    tab_footnote(
      footnote = md("**Target achievement**: Predicted prevalence <=70% of 2010 baseline (>=30% relative reduction). Based on point estimates."),
      locations = cells_column_labels(columns = achieved_2030_fmt)
    ) %>%
    tab_footnote(
      footnote = md("**On Track**: 2025 prevalence within 10% of target threshold."),
      locations = cells_column_labels(columns = ontrack_2025_fmt)
    ) %>%
    tab_footnote(
      footnote = md("Each country counts as one unit (simple count, NOT population-weighted)."),
      locations = cells_column_labels(columns = n_countries)
    ) %>%
    tab_source_note(
      source_note = md("*WHO FCTC targets. Predictions from Bayesian Age-Period-Cohort models. Current smoking = any tobacco use in past 30 days.*")
    ) %>%
    tab_options(
      table.width = pct(100),
      table.font.size = px(10),
      heading.title.font.size = px(14),
      column_labels.background.color = COLORS$header,
      column_labels.font.weight = "bold",
      row_group.background.color = COLORS$row_group,
      data_row.padding = px(6)
    )

  # Save table
  tryCatch({
    gtsave(gt_table2, file.path(output_dir, "Table2_WHO_Target_Achievement.html"))
    cat("    Table 2 saved to", file.path(output_dir, "Table2_WHO_Target_Achievement.html"), "\n")
  }, error = function(e) {
    cat("    Warning: Could not save Table 2 -", e$message, "\n")
  })

  return(gt_table2)
}

# ============================================================================
# TABLE 3: TOBACCO ENDGAME ACHIEVEMENT BY 2040
# ============================================================================

generate_table3 <- function(weighted_results,
                            country_region_mapping,
                            output_dir = "tables_publication") {

  cat("  Generating Table 3: Tobacco Endgame Achievement by 2040...\n")

  # Detect column names
  country_col <- if ("Country" %in% names(weighted_results)) "Country" else "wb_country_abv"
  sex_col <- if ("Sex" %in% names(weighted_results)) "Sex" else "sex"
  indicator_col <- if ("Def_Type_Code" %in% names(weighted_results)) "Def_Type_Code" else "def_type_code"
  year_col <- if ("Year" %in% names(weighted_results)) "Year" else "year"
  prev_col <- if ("weighted_mean" %in% names(weighted_results)) "weighted_mean" else "prevalence"

  # STEP 1: Get 2040 predictions
  predictions_2040 <- weighted_results %>%
    filter(!!sym(indicator_col) == PRIMARY_INDICATOR,
           !!sym(year_col) == ENDGAME_YEAR) %>%
    select(!!sym(country_col), !!sym(sex_col), !!sym(prev_col)) %>%
    rename(Country = !!sym(country_col), Sex = !!sym(sex_col),
           pred_2040 = !!sym(prev_col))

  # STEP 2: Classify each country
  endgame_status <- predictions_2040 %>%
    left_join(
      country_region_mapping %>% select(wb_country_abv, region_consolidated),
      by = c("Country" = "wb_country_abv")
    ) %>%
    filter(!is.na(region_consolidated)) %>%
    mutate(
      achieved_endgame = pred_2040 < ENDGAME_THRESHOLD,
      near_endgame = pred_2040 >= ENDGAME_THRESHOLD & pred_2040 < NEAR_ENDGAME_UPPER,
      virtual_elimination = pred_2040 < VIRTUAL_ELIMINATION
    )

  # STEP 3: Aggregate by region - COUNT COUNTRIES
  table3_data <- endgame_status %>%
    group_by(region_consolidated, Sex) %>%
    summarise(
      n_countries = n(),
      n_endgame = sum(achieved_endgame, na.rm = TRUE),
      n_near_endgame = sum(near_endgame, na.rm = TRUE),
      n_virtual_elim = sum(virtual_elimination, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      n_below_10 = n_endgame + n_near_endgame,
      endgame_fmt = format_count_pct(n_endgame, n_countries),
      near_endgame_fmt = format_count_pct(n_near_endgame, n_countries),
      virtual_elim_fmt = format_count_pct(n_virtual_elim, n_countries),
      below_10_fmt = format_count_pct(n_below_10, n_countries),
      Region_Name = tools::toTitleCase(region_consolidated),
      Gender = format_gender(Sex)
    ) %>%
    arrange(desc(Gender == "Men"), Region_Name)

  # STEP 4: Calculate subtotals
  subtotals <- table3_data %>%
    group_by(Gender) %>%
    summarise(
      Region_Name = "Subtotal",
      n_countries = sum(n_countries),
      n_endgame = sum(n_endgame),
      n_near_endgame = sum(n_near_endgame),
      n_virtual_elim = sum(n_virtual_elim),
      n_below_10 = sum(n_below_10),
      .groups = "drop"
    ) %>%
    mutate(
      endgame_fmt = format_count_pct(n_endgame, n_countries),
      near_endgame_fmt = format_count_pct(n_near_endgame, n_countries),
      virtual_elim_fmt = format_count_pct(n_virtual_elim, n_countries),
      below_10_fmt = format_count_pct(n_below_10, n_countries)
    )

  # Combine
  table3_final <- bind_rows(
    table3_data %>% filter(Gender == "Men") %>%
      select(Region_Name, Gender, n_countries, endgame_fmt, near_endgame_fmt,
             virtual_elim_fmt, below_10_fmt),
    subtotals %>% filter(Gender == "Men") %>%
      select(Region_Name, Gender, n_countries, endgame_fmt, near_endgame_fmt,
             virtual_elim_fmt, below_10_fmt),
    table3_data %>% filter(Gender == "Women") %>%
      select(Region_Name, Gender, n_countries, endgame_fmt, near_endgame_fmt,
             virtual_elim_fmt, below_10_fmt),
    subtotals %>% filter(Gender == "Women") %>%
      select(Region_Name, Gender, n_countries, endgame_fmt, near_endgame_fmt,
             virtual_elim_fmt, below_10_fmt)
  )

  # STEP 5: Create GT table
  n_men_rows <- sum(table3_final$Gender == "Men")
  n_women_rows <- sum(table3_final$Gender == "Women")

  gt_table3 <- table3_final %>%
    select(-Gender) %>%
    gt(rowname_col = "Region_Name") %>%
    tab_header(
      title = md("**Table 3: Tobacco Endgame Achievement by 2040**"),
      subtitle = md("Countries Projected to Achieve <5% Current Smoking Prevalence")
    ) %>%
    tab_row_group(
      label = md("**Men**"),
      rows = 1:n_men_rows
    ) %>%
    tab_row_group(
      label = md("**Women**"),
      rows = (n_men_rows + 1):(n_men_rows + n_women_rows)
    ) %>%
    row_group_order(groups = c("Men", "Women")) %>%
    cols_label(
      n_countries = "Countries (n)",
      endgame_fmt = "Endgame (<5%)",
      near_endgame_fmt = "Near-Endgame (5-10%)",
      virtual_elim_fmt = "Virtual Elimination (<2%)",
      below_10_fmt = "Total <10%"
    ) %>%
    tab_spanner(
      label = md("**2040 Endgame Assessment**"),
      columns = c(endgame_fmt, near_endgame_fmt, virtual_elim_fmt, below_10_fmt)
    ) %>%
    tab_style(
      style = list(
        cell_fill(color = COLORS$row_group),
        cell_text(weight = "bold")
      ),
      locations = cells_row_groups()
    ) %>%
    tab_style(
      style = list(
        cell_fill(color = COLORS$subtotal),
        cell_text(weight = "bold"),
        cell_borders(sides = "top", color = "#BDBDBD", weight = px(2))
      ),
      locations = cells_body(rows = Region_Name == "Subtotal")
    ) %>%
    tab_footnote(
      footnote = md("**Endgame** defined as <5% adult current smoking prevalence (New Zealand, Ireland, Scotland frameworks)."),
      locations = cells_column_labels(columns = endgame_fmt)
    ) %>%
    tab_footnote(
      footnote = md("**Near-Endgame**: 5-10% prevalence. **Virtual Elimination**: <2% prevalence."),
      locations = cells_column_labels(columns = c(near_endgame_fmt, virtual_elim_fmt))
    ) %>%
    tab_footnote(
      footnote = md("Each country counts as one unit (simple count, NOT population-weighted)."),
      locations = cells_column_labels(columns = n_countries)
    ) %>%
    tab_source_note(
      source_note = md("*Projections from Bayesian Age-Period-Cohort models extended to 2040. Current smoking = any tobacco use in past 30 days.*")
    ) %>%
    tab_options(
      table.width = pct(100),
      table.font.size = px(10),
      heading.title.font.size = px(14),
      column_labels.background.color = COLORS$header,
      column_labels.font.weight = "bold",
      row_group.background.color = COLORS$row_group,
      data_row.padding = px(6)
    )

  # Save table
  tryCatch({
    gtsave(gt_table3, file.path(output_dir, "Table3_Endgame_2040.html"))
    cat("    Table 3 saved to", file.path(output_dir, "Table3_Endgame_2040.html"), "\n")
  }, error = function(e) {
    cat("    Warning: Could not save Table 3 -", e$message, "\n")
  })

  return(gt_table3)
}

# ============================================================================
# TABLE S1: REGIONAL WHO TARGET ACHIEVEMENT BY INDICATOR (WIDE)
# ============================================================================

generate_table_s1 <- function(weighted_results,
                              country_region_mapping,
                              output_dir = "tables_publication") {

  cat("  Generating Table S1: Regional WHO Target Achievement by Indicator...\n")

  # Detect column names
  country_col <- if ("Country" %in% names(weighted_results)) "Country" else "wb_country_abv"
  sex_col <- if ("Sex" %in% names(weighted_results)) "Sex" else "sex"
  indicator_col <- if ("Def_Type_Code" %in% names(weighted_results)) "Def_Type_Code" else "def_type_code"
  year_col <- if ("Year" %in% names(weighted_results)) "Year" else "year"
  prev_col <- if ("weighted_mean" %in% names(weighted_results)) "weighted_mean" else "prevalence"

  # Get all indicators
  all_indicators <- unique(weighted_results[[indicator_col]])
  all_indicators <- all_indicators[all_indicators %in% names(INDICATOR_LABELS)]

  # For each indicator, calculate achievement
  results_list <- list()

  for (ind in all_indicators) {
    # Baseline
    baseline <- weighted_results %>%
      filter(!!sym(indicator_col) == ind,
             !!sym(year_col) == BASELINE_YEAR) %>%
      select(!!sym(country_col), !!sym(sex_col), !!sym(prev_col)) %>%
      rename(Country = !!sym(country_col), Sex = !!sym(sex_col),
             Baseline = !!sym(prev_col)) %>%
      mutate(Target = Baseline * (1 - REDUCTION_TARGET))

    # 2030 prediction
    pred_2030 <- weighted_results %>%
      filter(!!sym(indicator_col) == ind,
             !!sym(year_col) == TARGET_YEAR) %>%
      select(!!sym(country_col), !!sym(sex_col), !!sym(prev_col)) %>%
      rename(Country = !!sym(country_col), Sex = !!sym(sex_col),
             Predicted = !!sym(prev_col))

    # Join and determine achievement
    achievement <- baseline %>%
      left_join(pred_2030, by = c("Country", "Sex")) %>%
      left_join(
        country_region_mapping %>% select(wb_country_abv, region_consolidated),
        by = c("Country" = "wb_country_abv")
      ) %>%
      filter(!is.na(region_consolidated)) %>%
      mutate(Achieved = Predicted <= Target)

    # Aggregate by region
    region_summary <- achievement %>%
      group_by(region_consolidated, Sex) %>%
      summarise(
        n_countries = n(),
        n_achieved = sum(Achieved, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        Indicator = ind,
        achieved_fmt = format_count_pct(n_achieved, n_countries)
      )

    results_list[[ind]] <- region_summary
  }

  # Combine all results
  all_results <- bind_rows(results_list)

  # Pivot to wide format (one column per indicator)
  table_s1_wide <- all_results %>%
    mutate(
      Region_Name = tools::toTitleCase(region_consolidated),
      Gender = format_gender(Sex),
      Indicator_Label = INDICATOR_LABELS[Indicator]
    ) %>%
    select(Region_Name, Gender, Indicator_Label, achieved_fmt) %>%
    pivot_wider(
      names_from = Indicator_Label,
      values_from = achieved_fmt
    ) %>%
    arrange(desc(Gender == "Men"), Region_Name)

  # Create GT table
  gt_table_s1 <- table_s1_wide %>%
    gt() %>%
    tab_header(
      title = md("**Table S1: Regional WHO Target Achievement by Indicator**"),
      subtitle = md("Countries Achieving >=30% Reduction by 2030, n (%)")
    ) %>%
    cols_label(
      Region_Name = "WHO Region",
      Gender = "Gender"
    ) %>%
    tab_style(
      style = list(
        cell_fill(color = COLORS$header),
        cell_text(weight = "bold")
      ),
      locations = cells_column_labels()
    ) %>%
    tab_footnote(
      footnote = md("Each cell shows n (%) of countries achieving >=30% reduction from 2010 baseline by 2030."),
      locations = cells_title(groups = "subtitle")
    ) %>%
    tab_source_note(
      source_note = md("*Predictions from Bayesian Age-Period-Cohort models. Simple country counts (NOT population-weighted).*")
    ) %>%
    tab_options(
      table.width = pct(100),
      table.font.size = px(10),
      heading.title.font.size = px(14),
      column_labels.background.color = COLORS$header,
      data_row.padding = px(6)
    )

  # Save table
  tryCatch({
    gtsave(gt_table_s1, file.path(output_dir, "TableS1_Regional_Target_by_Indicator.html"))
    cat("    Table S1 saved to", file.path(output_dir, "TableS1_Regional_Target_by_Indicator.html"), "\n")
  }, error = function(e) {
    cat("    Warning: Could not save Table S1 -", e$message, "\n")
  })

  return(gt_table_s1)
}

# ============================================================================
# TABLE S3: COUNTRY-LEVEL TARGET ACHIEVEMENT (LONG FORMAT)
# ============================================================================

generate_table_s3 <- function(weighted_results,
                              country_region_mapping,
                              output_dir = "tables_publication") {

  cat("  Generating Table S3: Country-Level Target Achievement Details...\n")

  # Detect column names
  country_col <- if ("Country" %in% names(weighted_results)) "Country" else "wb_country_abv"
  sex_col <- if ("Sex" %in% names(weighted_results)) "Sex" else "sex"
  indicator_col <- if ("Def_Type_Code" %in% names(weighted_results)) "Def_Type_Code" else "def_type_code"
  year_col <- if ("Year" %in% names(weighted_results)) "Year" else "year"
  prev_col <- if ("weighted_mean" %in% names(weighted_results)) "weighted_mean" else "prevalence"

  # Get baseline
  baseline <- weighted_results %>%
    filter(!!sym(indicator_col) == PRIMARY_INDICATOR,
           !!sym(year_col) == BASELINE_YEAR) %>%
    select(!!sym(country_col), !!sym(sex_col), !!sym(prev_col)) %>%
    rename(Country = !!sym(country_col), Sex = !!sym(sex_col),
           Baseline = !!sym(prev_col)) %>%
    mutate(Target = Baseline * (1 - REDUCTION_TARGET))

  # Get 2030 prediction with CI if available
  ci_lower_col <- if ("weighted_lower" %in% names(weighted_results)) "weighted_lower" else NULL
  ci_upper_col <- if ("weighted_upper" %in% names(weighted_results)) "weighted_upper" else NULL

  pred_cols <- c(country_col, sex_col, prev_col)
  if (!is.null(ci_lower_col)) pred_cols <- c(pred_cols, ci_lower_col, ci_upper_col)

  pred_2030 <- weighted_results %>%
    filter(!!sym(indicator_col) == PRIMARY_INDICATOR,
           !!sym(year_col) == TARGET_YEAR) %>%
    select(all_of(pred_cols)) %>%
    rename(Country = !!sym(country_col), Sex = !!sym(sex_col),
           Predicted = !!sym(prev_col))

  if (!is.null(ci_lower_col)) {
    pred_2030 <- pred_2030 %>%
      rename(CI_Lower = !!sym(ci_lower_col), CI_Upper = !!sym(ci_upper_col))
  }

  # Join and calculate
  country_details <- baseline %>%
    left_join(pred_2030, by = c("Country", "Sex")) %>%
    left_join(
      country_region_mapping %>%
        select(wb_country_abv, region_consolidated, country_name) %>%
        distinct(),
      by = c("Country" = "wb_country_abv")
    ) %>%
    filter(!is.na(region_consolidated)) %>%
    mutate(
      Achieved = Predicted <= Target,
      Baseline_pct = sprintf("%.1f%%", Baseline * 100),
      Target_pct = sprintf("%.1f%%", Target * 100),
      Predicted_pct = sprintf("%.1f%%", Predicted * 100),
      Status = ifelse(Achieved, "Achieved", "Not Achieved"),
      Region_Name = tools::toTitleCase(region_consolidated),
      Country_Name = tools::toTitleCase(country_name),
      Gender = format_gender(Sex)
    )

  # Add CI if available
  if ("CI_Lower" %in% names(country_details)) {
    country_details <- country_details %>%
      mutate(
        CI_fmt = sprintf("(%.1f%%-%.1f%%)", CI_Lower * 100, CI_Upper * 100)
      )
  } else {
    country_details <- country_details %>%
      mutate(CI_fmt = "-")
  }

  # Select and order columns
  table_s3_data <- country_details %>%
    select(Region_Name, Country_Name, Gender, Baseline_pct, Target_pct,
           Predicted_pct, CI_fmt, Status) %>%
    arrange(Region_Name, Country_Name, desc(Gender == "Men"))

  # Create GT table
  gt_table_s3 <- table_s3_data %>%
    gt() %>%
    tab_header(
      title = md("**Table S3: Country-Level Target Achievement Details**"),
      subtitle = md("Current Smoking Prevalence - 2030 Assessment")
    ) %>%
    cols_label(
      Region_Name = "WHO Region",
      Country_Name = "Country",
      Gender = "Gender",
      Baseline_pct = "2010 Baseline",
      Target_pct = "Target (70%)",
      Predicted_pct = "2030 Predicted",
      CI_fmt = "95% CI",
      Status = "Achievement"
    ) %>%
    tab_style(
      style = cell_fill(color = COLORS$good),
      locations = cells_body(
        columns = Status,
        rows = Status == "Achieved"
      )
    ) %>%
    tab_style(
      style = cell_fill(color = COLORS$low),
      locations = cells_body(
        columns = Status,
        rows = Status == "Not Achieved"
      )
    ) %>%
    tab_style(
      style = list(
        cell_fill(color = COLORS$header),
        cell_text(weight = "bold")
      ),
      locations = cells_column_labels()
    ) %>%
    tab_footnote(
      footnote = md("Target = 70% of 2010 baseline (30% relative reduction)."),
      locations = cells_column_labels(columns = Target_pct)
    ) %>%
    tab_source_note(
      source_note = md("*Current smoking = any tobacco product use in past 30 days. Predictions from Bayesian APC models.*")
    ) %>%
    tab_options(
      table.width = pct(100),
      table.font.size = px(9),
      heading.title.font.size = px(14),
      column_labels.background.color = COLORS$header,
      data_row.padding = px(4)
    )

  # Save table
  tryCatch({
    gtsave(gt_table_s3, file.path(output_dir, "TableS3_Country_Achievement_Details.html"))
    cat("    Table S3 saved to", file.path(output_dir, "TableS3_Country_Achievement_Details.html"), "\n")
  }, error = function(e) {
    cat("    Warning: Could not save Table S3 -", e$message, "\n")
  })

  return(gt_table_s3)
}

# ============================================================================
# TABLE S4: MODEL SELECTION AND RMSE
# ============================================================================

generate_table_s4 <- function(model_selection_results,
                              country_region_mapping,
                              output_dir = "tables_publication") {

  cat("  Generating Table S4: Model Selection and RMSE...\n")

  if (is.null(model_selection_results) || nrow(model_selection_results) == 0) {
    cat("    Warning: No model selection data available. Skipping Table S4.\n")
    return(NULL)
  }

  # Detect column names
  country_col <- if ("Country" %in% names(model_selection_results)) "Country" else "wb_country_abv"
  sex_col <- if ("Sex" %in% names(model_selection_results)) "Sex" else "sex"

  # Process model selection data
  table_s4_data <- model_selection_results %>%
    left_join(
      country_region_mapping %>%
        select(wb_country_abv, region_consolidated, country_name) %>%
        distinct(),
      by = c(setNames("wb_country_abv", country_col))
    ) %>%
    filter(!is.na(region_consolidated)) %>%
    mutate(
      Region_Name = tools::toTitleCase(region_consolidated),
      Country_Name = tools::toTitleCase(country_name),
      Gender = format_gender(!!sym(sex_col)),
      RMSE_Global = ifelse(!is.na(RMSE_Global), sprintf("%.4f", RMSE_Global), "-"),
      RMSE_Country = ifelse(!is.na(RMSE_Country), sprintf("%.4f", RMSE_Country), "-"),
      Selected_Model = ifelse(Final_Selected_Model == "Global", "Global APC", "Country APC")
    ) %>%
    select(Region_Name, Country_Name, Gender, n_observations,
           RMSE_Global, RMSE_Country, Selected_Model) %>%
    arrange(Region_Name, Country_Name, desc(Gender == "Men"))

  # Create GT table
  gt_table_s4 <- table_s4_data %>%
    gt() %>%
    tab_header(
      title = md("**Table S4: Model Selection and Performance Metrics**"),
      subtitle = md("RMSE-Based Model Selection by Country and Gender")
    ) %>%
    cols_label(
      Region_Name = "WHO Region",
      Country_Name = "Country",
      Gender = "Gender",
      n_observations = "Observations (n)",
      RMSE_Global = "Global APC RMSE",
      RMSE_Country = "Country APC RMSE",
      Selected_Model = "Selected Model"
    ) %>%
    tab_spanner(
      label = md("**RMSE Performance**"),
      columns = c(RMSE_Global, RMSE_Country)
    ) %>%
    tab_style(
      style = list(
        cell_fill(color = COLORS$header),
        cell_text(weight = "bold")
      ),
      locations = cells_column_labels()
    ) %>%
    tab_footnote(
      footnote = md("RMSE = Root Mean Square Error. Lower values indicate better fit."),
      locations = cells_column_spanners(spanners = "RMSE Performance")
    ) %>%
    tab_footnote(
      footnote = md("Country APC model selected when RMSE is lower than Global APC and sufficient data available."),
      locations = cells_column_labels(columns = Selected_Model)
    ) %>%
    tab_source_note(
      source_note = md("*'-' indicates insufficient data for country-specific model estimation.*")
    ) %>%
    tab_options(
      table.width = pct(100),
      table.font.size = px(9),
      heading.title.font.size = px(14),
      column_labels.background.color = COLORS$header,
      data_row.padding = px(4)
    )

  # Save table
  tryCatch({
    gtsave(gt_table_s4, file.path(output_dir, "TableS4_Model_Selection_RMSE.html"))
    cat("    Table S4 saved to", file.path(output_dir, "TableS4_Model_Selection_RMSE.html"), "\n")
  }, error = function(e) {
    cat("    Warning: Could not save Table S4 -", e$message, "\n")
  })

  return(gt_table_s4)
}

# ============================================================================
# TABLE S5: DATA SOURCES BY COUNTRY
# ============================================================================

generate_table_s5 <- function(clean_data,
                              country_region_mapping,
                              output_dir = "tables_publication") {

  cat("  Generating Table S5: Data Sources by Country...\n")

  # Summarize surveys by country
  survey_summary <- clean_data %>%
    left_join(
      country_region_mapping %>%
        select(wb_country_abv, region_consolidated, country_name) %>%
        distinct(),
      by = "wb_country_abv"
    ) %>%
    filter(!is.na(region_consolidated)) %>%
    group_by(region_consolidated, wb_country_abv, country_name) %>%
    summarise(
      n_observations = n(),
      n_surveys = n_distinct(survey),
      year_min = min(year, na.rm = TRUE),
      year_max = max(year, na.rm = TRUE),
      surveys = paste(sort(unique(survey)), collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(
      Region_Name = tools::toTitleCase(region_consolidated),
      Country_Name = tools::toTitleCase(country_name),
      Year_Range = paste0(year_min, "-", year_max)
    ) %>%
    select(Region_Name, Country_Name, n_observations, n_surveys, Year_Range, surveys) %>%
    arrange(Region_Name, Country_Name)

  # Create GT table
  gt_table_s5 <- survey_summary %>%
    gt() %>%
    tab_header(
      title = md("**Table S5: Data Sources by Country**"),
      subtitle = md("Survey Names, Years, and Sample Sizes")
    ) %>%
    cols_label(
      Region_Name = "WHO Region",
      Country_Name = "Country",
      n_observations = "Observations (n)",
      n_surveys = "Surveys (n)",
      Year_Range = "Survey Years",
      surveys = "Survey Names"
    ) %>%
    fmt_number(
      columns = n_observations,
      use_seps = TRUE,
      decimals = 0
    ) %>%
    tab_style(
      style = list(
        cell_fill(color = COLORS$header),
        cell_text(weight = "bold")
      ),
      locations = cells_column_labels()
    ) %>%
    tab_style(
      style = cell_text(size = px(8)),
      locations = cells_body(columns = surveys)
    ) %>%
    tab_source_note(
      source_note = md("*Data source: WHO Global Tobacco Surveillance System (GTSS), including GATS, STEPS, and national surveys.*")
    ) %>%
    tab_options(
      table.width = pct(100),
      table.font.size = px(9),
      heading.title.font.size = px(14),
      column_labels.background.color = COLORS$header,
      data_row.padding = px(4)
    )

  # Save table
  tryCatch({
    gtsave(gt_table_s5, file.path(output_dir, "TableS5_Data_Sources.html"))
    cat("    Table S5 saved to", file.path(output_dir, "TableS5_Data_Sources.html"), "\n")
  }, error = function(e) {
    cat("    Warning: Could not save Table S5 -", e$message, "\n")
  })

  return(gt_table_s5)
}

# ============================================================================
# MASTER EXECUTION
# ============================================================================

# Create output directory
dir.create("tables_publication", showWarnings = FALSE, recursive = TRUE)

# Check for required data
if (!exists("clean_data")) {
  cat("  ERROR: clean_data not found. Cannot generate tables.\n")
} else if (!exists("country_region_mapping") && !exists("country_region_manual")) {
  cat("  ERROR: country_region_mapping not found. Cannot generate tables.\n")
} else {

  # Use appropriate region mapping
  if (!exists("country_region_mapping") && exists("country_region_manual")) {
    country_region_mapping <- country_region_manual
  }

  # Detect weighted results
  weighted_results <- NULL
  if (exists("final_weighted_results_selected")) {
    weighted_results <- final_weighted_results_selected
    cat("  Using final_weighted_results_selected for tables.\n")
  } else if (exists("combined_apc_weighted_results")) {
    weighted_results <- combined_apc_weighted_results
    cat("  Using combined_apc_weighted_results for tables.\n")
  }

  if (is.null(weighted_results)) {
    cat("  WARNING: No weighted results found. Some tables may be skipped.\n")
  }

  # Generate tables
  tables_generated <- list()

  # Table 1: Baseline and Trends
  if (!is.null(weighted_results)) {
    tryCatch({
      tables_generated$table1 <- generate_table1(
        clean_data = clean_data,
        weighted_results = weighted_results,
        country_region_mapping = country_region_mapping,
        output_dir = "tables_publication"
      )
    }, error = function(e) {
      cat("  ERROR generating Table 1:", e$message, "\n")
    })
  }

  # Table 2: WHO Target Achievement
  if (!is.null(weighted_results)) {
    tryCatch({
      tables_generated$table2 <- generate_table2(
        weighted_results = weighted_results,
        country_region_mapping = country_region_mapping,
        output_dir = "tables_publication"
      )
    }, error = function(e) {
      cat("  ERROR generating Table 2:", e$message, "\n")
    })
  }

  # Table 3: Endgame 2040
  if (!is.null(weighted_results)) {
    tryCatch({
      tables_generated$table3 <- generate_table3(
        weighted_results = weighted_results,
        country_region_mapping = country_region_mapping,
        output_dir = "tables_publication"
      )
    }, error = function(e) {
      cat("  ERROR generating Table 3:", e$message, "\n")
    })
  }

  # Table S1: Regional by Indicator
  if (!is.null(weighted_results)) {
    tryCatch({
      tables_generated$table_s1 <- generate_table_s1(
        weighted_results = weighted_results,
        country_region_mapping = country_region_mapping,
        output_dir = "tables_publication"
      )
    }, error = function(e) {
      cat("  ERROR generating Table S1:", e$message, "\n")
    })
  }

  # Table S3: Country Details
  if (!is.null(weighted_results)) {
    tryCatch({
      tables_generated$table_s3 <- generate_table_s3(
        weighted_results = weighted_results,
        country_region_mapping = country_region_mapping,
        output_dir = "tables_publication"
      )
    }, error = function(e) {
      cat("  ERROR generating Table S3:", e$message, "\n")
    })
  }

  # Table S4: Model Selection
  if (exists("country_sex_summary")) {
    tryCatch({
      tables_generated$table_s4 <- generate_table_s4(
        model_selection_results = country_sex_summary,
        country_region_mapping = country_region_mapping,
        output_dir = "tables_publication"
      )
    }, error = function(e) {
      cat("  ERROR generating Table S4:", e$message, "\n")
    })
  } else {
    cat("  Skipping Table S4: country_sex_summary not available.\n")
  }

  # Table S5: Data Sources
  tryCatch({
    tables_generated$table_s5 <- generate_table_s5(
      clean_data = clean_data,
      country_region_mapping = country_region_mapping,
      output_dir = "tables_publication"
    )
  }, error = function(e) {
    cat("  ERROR generating Table S5:", e$message, "\n")
  })

  # Summary
  n_generated <- sum(!sapply(tables_generated, is.null))
  cat("\n========== TABLE GENERATION COMPLETE ==========\n")
  cat(sprintf("  Generated %d tables in tables_publication/\n", n_generated))
  cat("  Files:\n")
  for (f in list.files("tables_publication", pattern = "\\.html$")) {
    cat(sprintf("    - %s\n", f))
  }
}

cat("\n  Publication tables module complete.\n")
