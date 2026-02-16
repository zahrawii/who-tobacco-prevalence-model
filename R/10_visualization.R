#########################################################################################
#
#                    WHO TOBACCO CONTROL PREVALENCE PROJECTION MODEL
#                       10_visualization.R - All Plotting Functions
#
#   Contains: Sections 16-19 from pipeline_monolith.R
#     - WHO color scheme and master theme
#     - Helper functions (format_indicator)
#     - Data loading for visualization
#     - Plot modules (trends, validation, aggregated, cohorts, maps)
#     - Master execution loop
#
#   Requires: All previous modules (00-09) to be sourced
#   Outputs: All visualization files in outputs/figures/
#
#########################################################################################

cat("\n================================================================\n")
cat("  VISUALIZATION MODULE\n")
cat("================================================================\n")

#########################################################################################
#                              SECTION 16: STYLE CONFIGURATION                          #
#########################################################################################

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
  "Global" = "#2C3E50"
)

# Indicator Labels Map
format_indicator <- function(x) {
  case_when(
    x == "daily_user_cigarettes" ~ "Daily Cigarettes",
    x == "current_user_cigarettes" ~ "Current Cigarettes",
    x == "daily_user_any_tobacco_product" ~ "Daily Any Tobacco",
    x == "current_user_any_tobacco_product" ~ "Current Any Tobacco",
    x == "daily_user_any_smoked_tobacco" ~ "Daily Any Smoked",
    x == "current_user_any_smoked_tobacco" ~ "Current Any Smoked",
    TRUE ~ gsub("_", " ", stringr::str_to_title(x))
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
df_preds_global <- NULL
df_preds_country <- NULL
df_preds_all <- NULL

if (file.exists("results/final_ac_predictions_nested.csv")) {
  df_preds_global <- read.csv("results/final_ac_predictions_nested.csv") %>% mutate(Model_Type = "Global")
  cat("    - Global predictions loaded\n")
}
if (file.exists("results/final_predictions_country_specific.csv")) {
  df_preds_country <- read.csv("results/final_predictions_country_specific.csv") %>% mutate(Model_Type = "Country")
  cat("    - Country predictions loaded\n")
}
if (!is.null(df_preds_global) && !is.null(df_preds_country)) {
  df_preds_all <- bind_rows(df_preds_global, df_preds_country)
} else if (!is.null(df_preds_global)) {
  df_preds_all <- df_preds_global
} else if (!is.null(df_preds_country)) {
  df_preds_all <- df_preds_country
}

# 2. Load Trends (Weighted) - Tier 2
df_trends <- NULL
if (file.exists("results/final_weighted_results_apc_post_selection.csv")) {
  df_trends <- read.csv("results/final_weighted_results_apc_post_selection.csv")
  cat("    - Weighted trends loaded\n")
}

# 3. Load Regional - Tier 3
df_regional <- NULL
if (file.exists("results/regional_aggregation/regional_full_summary_lancet.csv")) {
  df_regional <- read.csv("results/regional_aggregation/regional_full_summary_lancet.csv")
  cat("    - Regional data loaded\n")
}

# 4. Load Observations
df_obs <- NULL
if (exists("clean_data")) {
  df_obs <- clean_data %>%
    mutate(
      Year = as.numeric(year),
      Prevalence = plogis(prevalence),
      Indicator = def_type_code
    )
  cat("    - Observations prepared\n")
}

#########################################################################################
#                              SECTION 18: PLOT MODULES                                 #
#########################################################################################

# ---- Module 1: Trend Plots (Trends.R equivalent) ----

plot_age_panel <- function(country_code, gender, indicator) {
  if (is.null(df_preds_all)) return(NULL)

  dat_p <- df_preds_all %>%
    filter(Country == country_code, Sex == gender, Def_Type_Code == indicator) %>%
    filter(Age_Midpoint %in% seq(20, 80, by = 10))

  if(nrow(dat_p) == 0) return(NULL)

  country_name <- if (exists("country_name_mapping")) country_name_mapping[country_code] else country_code

  p <- ggplot(dat_p, aes(x = Year, y = Prevalence)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Model_Type), alpha = 0.2) +
    geom_line(aes(color = Model_Type), linewidth = 0.8) +
    facet_wrap(~Age_Midpoint, labeller = label_both, ncol = 4) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
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
  if (is.null(df_preds_all)) return(NULL)

  dat_p <- df_preds_all %>%
    filter(Country == country_code, Sex == gender, Def_Type_Code == indicator,
           Model_Type == "Country")

  if(nrow(dat_p) == 0) return(NULL)
  country_name <- if (exists("country_name_mapping")) country_name_mapping[country_code] else country_code

  p <- ggplot(dat_p, aes(x = Year, y = Age_Midpoint, fill = Prevalence)) +
    geom_tile() +
    scale_fill_viridis(option = "magma", labels = scales::percent, direction = -1) +
    scale_x_continuous(expand = c(0,0), breaks = seq(2010, 2030, 5)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(
      title = paste0(country_name, ": Prevalence Heatmap"),
      subtitle = paste0(format_indicator(indicator), " | ", tools::toTitleCase(gender)),
      y = "Age"
    ) +
    coord_fixed(ratio = 0.5)
  return(p)
}

# ---- Module 2: Validation Plots (Validation.R equivalent) ----

plot_validation_facet <- function(country_code, gender, year) {
  if (is.null(df_preds_all) || is.null(df_obs)) return(NULL)

  dat_p <- df_preds_all %>%
    filter(Country == country_code, Sex == gender, Year == year)

  dat_o <- df_obs %>%
    filter(wb_country_abv == country_code, sex == gender, abs(Year - year) <= 2)

  if(nrow(dat_p) == 0) return(NULL)
  country_name <- if (exists("country_name_mapping")) country_name_mapping[country_code] else country_code

  p <- ggplot() +
    geom_ribbon(data = dat_p, aes(x = Age_Midpoint, ymin = lower_ci, ymax = upper_ci, fill = Model_Type), alpha = 0.2) +
    geom_line(data = dat_p, aes(x = Age_Midpoint, y = Prevalence, color = Model_Type)) +
    geom_point(data = dat_o, aes(x = Age_Midpoint, y = Prevalence, shape = "Survey Data"), size = 2) +
    geom_segment(data = dat_o, aes(x = start_age, xend = end_age, y = Prevalence, yend = Prevalence),
                 linetype = "dotted", alpha = 0.5) +
    facet_wrap(~Def_Type_Code, labeller = as_labeller(format_indicator), scales = "free_y") +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values = c("Country" = who_colors$accent, "Global" = who_colors$secondary)) +
    scale_fill_manual(values = c("Country" = who_colors$accent, "Global" = who_colors$secondary)) +
    labs(
      title = paste0("Model Validation: ", country_name, " (", year, ")"),
      subtitle = paste0("Comparing Model vs Survey Data | ", tools::toTitleCase(gender)),
      x = "Age", y = "Prevalence"
    )
  return(p)
}

# ---- Module 3: Aggregated Plots (Aggregated.R equivalent) ----

plot_regional_comparison <- function(indicator, gender) {
  if (is.null(df_regional)) return(NULL)

  dat_r <- df_regional %>%
    filter(Sex == gender, Def_Type_Code == indicator) %>%
    mutate(Region_Name = tools::toTitleCase(Region))

  p <- ggplot(dat_r) +
    geom_segment(aes(x = Base_Prevalence, xend = Projected_Prevalence,
                     y = reorder(Region_Name, Projected_Prevalence),
                     yend = reorder(Region_Name, Projected_Prevalence)), color = "grey") +
    geom_point(aes(x = Base_Prevalence, y = Region_Name, color = "2010 Base"), size = 3) +
    geom_point(aes(x = Projected_Prevalence, y = Region_Name, color = "2025 Projected"), size = 3) +
    geom_vline(xintercept = 0.30, linetype = "dashed", color = "grey") +
    scale_x_continuous(labels = scales::percent) +
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
  if (is.null(df_preds_all)) return(NULL)

  dat_p <- df_preds_all %>%
    filter(Country == country_code, Sex == gender, Def_Type_Code == indicator, Model_Type == "Country")

  if(nrow(dat_p) == 0) return(NULL)

  cohorts <- seq(1940, 2000, by = 20)
  country_name <- if (exists("country_name_mapping")) country_name_mapping[country_code] else country_code

  p <- ggplot(dat_p, aes(x = Year, y = Age_Midpoint)) +
    geom_tile(aes(fill = Prevalence)) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1, labels = scales::percent) +
    geom_abline(intercept = -cohorts, slope = 1, color = "white", linetype = "dashed", alpha = 0.5) +
    annotate("text", x = 2015, y = 2015 - cohorts, label = paste0("Born ", cohorts),
             angle = 45, color = "white", size = 3, vjust = -0.5) +
    coord_fixed(ratio = 1, xlim = c(2000, 2030), ylim = c(15, 80)) +
    labs(
      title = paste0("Lexis Diagram: ", country_name),
      subtitle = paste0("Diagonal lines track specific birth cohorts"),
      fill = "Prev."
    )
  return(p)
}

# ---- Module 5: Maps (Maps.R equivalent) ----

plot_prevalence_map_dual <- function(year, indicator) {
  if (is.null(df_trends)) return(NULL)

  world <- ne_countries(scale = "medium", returnclass = "sf")

  # Fix Natural Earth iso_a3 = "-99" issue for France, Norway, etc.
  # Use adm0_a3 as fallback when iso_a3 is invalid
  world <- world %>%
    mutate(
      iso_a3 = case_when(
        is.na(iso_a3) | iso_a3 == "-99" ~ adm0_a3,
        TRUE ~ iso_a3
      )
    )

  dat_map <- df_trends %>%
    filter(Year == year, Def_Type_Code == indicator) %>%
    mutate(
      iso_a3 = toupper(Country),
      Gender = case_when(tolower(Sex) == "males" ~ "Men", tolower(Sex) == "females" ~ "Women", TRUE ~ Sex)
    ) %>%
    select(iso_a3, Gender, weighted_mean)

  map_data <- world %>%
    left_join(dat_map, by = "iso_a3") %>%
    filter(!is.na(Gender))

  p <- ggplot(map_data) +
    geom_sf(aes(fill = weighted_mean), color = "white", linewidth = 0.1) +
    facet_wrap(~Gender, ncol = 1) +
    scale_fill_viridis(option = "rocket", direction = -1, labels = scales::percent, name = "Prevalence") +
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

  # Create directory structure
  dirs <- c("trends/panels", "trends/heatmaps", "validation", "aggregated", "cohorts", "maps")
  purrr::walk(dirs, ~dir.create(file.path(output_base, .x), recursive = TRUE, showWarnings = FALSE))

  cat("  Starting batch generation...\n")

  # Check if data is available
  if (is.null(df_trends) && is.null(df_preds_all)) {
    cat("  No data available for visualization. Run model first.\n")
    return(invisible(NULL))
  }

  # 1. MAPS (Global)
  if (!is.null(df_trends)) {
    cat("  Generating Maps...\n")
    for(ind in unique(df_trends$Def_Type_Code)) {
      tryCatch({
        p <- plot_prevalence_map_dual(TARGET_YEAR, ind)
        if (!is.null(p)) {
          ggsave(file.path(output_base, "maps", paste0("map_", ind, "_", TARGET_YEAR, ".pdf")), p, width = 8, height = 10)
        }
      }, error = function(e) message("  Error generating map: ", e$message))
    }
  }

  # 2. AGGREGATED (Regional)
  if (!is.null(df_regional)) {
    cat("  Generating Regional Plots...\n")
    for(ind in unique(df_regional$Def_Type_Code)) {
      for(g in c("males", "females")) {
        tryCatch({
          p <- plot_regional_comparison(ind, g)
          if (!is.null(p)) {
            ggsave(file.path(output_base, "aggregated", paste0("regional_", ind, "_", g, ".pdf")), p, width = 10, height = 6)
          }
        }, error = function(e) {})
      }
    }
  }

  # 3. COUNTRY LOOPS
  if (!is.null(df_preds_all)) {
    countries <- unique(df_preds_all$Country)
    cat(sprintf("  Generating Country Plots for %d countries...\n", length(countries)))

    pb <- txtProgressBar(min = 0, max = length(countries), style = 3)

    for(i in seq_along(countries)) {
      cntry <- countries[i]
      c_name <- gsub(" ", "_", if (exists("country_name_mapping")) country_name_mapping[cntry] else cntry)

      # Create Country Subfolder
      c_dir <- file.path(output_base, "trends", c_name)
      dir.create(c_dir, showWarnings = FALSE)

      for(g in c("males", "females")) {
        for(ind in unique(df_preds_all$Def_Type_Code)) {
          tryCatch({
            # Trend Panel
            p1 <- plot_age_panel(cntry, g, ind)
            if(!is.null(p1)) ggsave(file.path(c_dir, paste0("panel_", g, "_", ind, ".pdf")), p1, width = 12, height = 8)

            # Heatmap
            p2 <- plot_age_year_heatmap(cntry, g, ind)
            if(!is.null(p2)) ggsave(file.path(c_dir, paste0("heatmap_", g, "_", ind, ".png")), p2, width = 8, height = 6)

            # Cohort Lexis (Only needed once per country/gender usually)
            if(ind == "daily_user_cigarettes") {
              p3 <- plot_lexis_diagram(cntry, g, ind)
              if(!is.null(p3)) ggsave(file.path(output_base, "cohorts", paste0("lexis_", c_name, "_", g, ".pdf")), p3, width = 8, height = 8)
            }

          }, error = function(e) {})
        }

        # Validation Plot (Per Year)
        tryCatch({
          p_val <- plot_validation_facet(cntry, g, 2015)
          if(!is.null(p_val)) ggsave(file.path(output_base, "validation", paste0("val_", c_name, "_", g, "_2015.pdf")), p_val, width = 10, height = 6)
        }, error = function(e) {})
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }

  cat("\n  Visualization generation complete.\n")
}

cat("  Visualization module loaded.\n")
cat("  Run generate_all_visualizations() to create plots.\n")
