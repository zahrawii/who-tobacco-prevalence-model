#########################################################################################
#
#                    WHO TOBACCO CONTROL PREVALENCE PROJECTION MODEL
#                       08_evaluation.R - Model Evaluation & Selection
#
#   Contains: Section 13-14 from pipeline_monolith.R
#     - RMSE calculation for global vs country models
#     - Model selection per country-sex combination
#     - Post-selection combined results creation
#
#   Requires: All previous modules (00-07) to be sourced
#   Outputs: country_sex_summary, final_weighted_results_selected
#
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
  filter(Model_Type == Final_Selected_Model) %>%
  select(-Final_Selected_Model)

# Write to CSV for later use
write.csv(
  final_weighted_results_selected,
  file = "results/final_weighted_results_apc_post_selection.csv",
  row.names = FALSE
)

cat(sprintf("  Post-selection results: %d rows\n", nrow(final_weighted_results_selected)))
cat("  Saved to: results/final_weighted_results_apc_post_selection.csv\n")

cat("\nModel evaluation module complete.\n")
