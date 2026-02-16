#########################################################################################
#
#                    WHO TOBACCO CONTROL PREVALENCE PROJECTION MODEL
#                          01_data_prep.R - Data Loading & Cleaning
#
#   Contains: Data loading, variable standardization, indicator creation
#   Requires: 00_config.R to be sourced first
#   Outputs: clean_data, country_name_mapping
#
#########################################################################################

# ---- 1.1 Load WHO Prevalence Data ----

data <- read_excel(
  "data/Prevalence RESHAPE 28 Sep 2024.xlsx",
  sheet = "reshape.dta",
  skip = 1
) %>%
  clean_names()

# ---- 1.2 Clean and Standardize Variable Names ----

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

# ---- 1.3 Create Indicator Labels ----

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

# ---- 1.4 Create Binary Indicators ----

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

# ---- 1.4.1 Validate Tobacco Type Indicators ----

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

# ---- 1.5 Remove Redundant Total Prevalence Observations ----

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

# ---- 1.6 Create Country Name Lookup ----

country_name_mapping <- clean_data %>%
  distinct(wb_country_abv, country) %>%
  arrange(wb_country_abv) %>%
  { setNames(.$country, .$wb_country_abv) }

cat("\nData preparation complete.\n")
cat(sprintf("  Total observations: %d\n", nrow(clean_data)))
cat(sprintf("  Countries: %d\n", length(unique(clean_data$wb_country_abv))))
