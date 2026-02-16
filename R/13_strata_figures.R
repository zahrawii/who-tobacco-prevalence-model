#########################################################################################
#
#                    WHO TOBACCO PIPELINE - MODULE 13: STRATA-LEVEL FIGURES
#
#   Comprehensive visualization with proper loop structure:
#
#   AGE CURVES (eFigure 2):
#     - Folder: age_curves/{country}/{gender}/{year}/
#     - Plot A: Faceted by indicator (6 panels)
#     - Plot B: Overlay all indicators (colored lines)
#
#   WEIGHTED TRENDS (eFigure 3):
#     - Folder: weighted_trends/{country}/{gender}/
#     - Plot A: Faceted by indicator
#     - Plot B: Overlay all indicators
#     - Plot C: Individual indicator files
#
#   VALIDATION (eFigure 1):
#     - Global scatter (obs vs pred)
#     - Faceted by age group
#     - Faceted by indicator
#     - Residuals by age (spline quality)
#
#########################################################################################

cat("\n")
cat("###########################################################################\n")
cat("#                                                                         #\n")
cat("#        MODULE 13: STRATA-LEVEL PUBLICATION FIGURES                      #\n")
cat("#                                                                         #\n")
cat("###########################################################################\n\n")

# ============================================================================
# DEPENDENCIES
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(viridis)
  library(patchwork)
})

# ============================================================================
# CONFIGURATION
# ============================================================================

# Base output directory
FIG_BASE_DIR <- "outputs/figures"

# Year parameters
BASE_YEAR <- 2010
TARGET_YEAR <- 2025
ENDGAME_YEAR <- 2040
ENDGAME_THRESHOLD <- 5  # 5%

# Spline transition
TRANSITION_AGE <- 65

# Indicator colors (colorblind-friendly)
INDICATOR_COLORS <- c(
  "current_user_any_tobacco_product" = "#E41A1C",
  "current_user_any_smoked_tobacco" = "#377EB8",
  "current_user_cigarettes" = "#4DAF4A",
  "daily_user_any_tobacco_product" = "#984EA3",
  "daily_user_any_smoked_tobacco" = "#FF7F00",
  "daily_user_cigarettes" = "#A65628"
)

# Indicator labels
INDICATOR_LABELS <- c(
  "current_user_any_tobacco_product" = "Current Any Tobacco",
  "current_user_any_smoked_tobacco" = "Current Smoked Tobacco",
  "current_user_cigarettes" = "Current Cigarettes",
  "daily_user_any_tobacco_product" = "Daily Any Tobacco",
  "daily_user_any_smoked_tobacco" = "Daily Smoked Tobacco",
  "daily_user_cigarettes" = "Daily Cigarettes"
)

# Model colors
MODEL_COLORS <- c(
  "Global" = "#2166AC",
  "Country" = "#D55E00",
  "Observed" = "#000000"
)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Format gender for display
format_gender_label <- function(sex) {
  case_when(
    tolower(sex) %in% c("males", "male") ~ "Men",
    tolower(sex) %in% c("females", "female") ~ "Women",
    TRUE ~ as.character(sex)
  )
}

#' Format indicator for display
format_indicator_label <- function(indicator) {
  if (indicator %in% names(INDICATOR_LABELS)) {
    return(INDICATOR_LABELS[indicator])
  }
  return(gsub("_", " ", tools::toTitleCase(indicator)))
}

#' Create age group bins
create_age_groups <- function(ages, breaks = c(15, 25, 35, 45, 55, 65, 75, 85, 100)) {
  cut(ages, breaks = breaks, labels = paste0(head(breaks, -1), "-", tail(breaks, -1) - 1),
      include.lowest = TRUE, right = FALSE)
}

#' WHO publication theme
theme_who_pub <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.1), hjust = 0),
      plot.subtitle = element_text(color = "#555555", size = rel(0.9)),
      plot.caption = element_text(color = "#777777", size = rel(0.7), hjust = 1),
      axis.title = element_text(face = "bold", size = rel(0.9)),
      axis.text = element_text(color = "#333333"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.3),
      strip.background = element_rect(fill = "#F5F5F5", color = NA),
      strip.text = element_text(face = "bold", size = rel(0.85)),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = rel(0.8))
    )
}

#' Create directory structure
create_figure_dirs <- function(base_dir = FIG_BASE_DIR) {
  dirs <- c(
    "age_curves",
    "weighted_trends",
    "validation",
    "validation/by_age_group",
    "validation/by_indicator",
    "validation/by_region",
    "validation/by_country"
  )
  for (d in dirs) {
    dir.create(file.path(base_dir, d), recursive = TRUE, showWarnings = FALSE)
  }
}

# ============================================================================
# AGE CURVE FUNCTIONS
# ============================================================================

#' Create age curve plot - Faceted by Indicator
#'
#' @param obs_data Observed data for this country-gender-year
#' @param pred_data Predicted data (can include Global and/or Country models)
#' @param country_name Display name for country
#' @param gender Gender label
#' @param year Year
plot_age_curve_faceted <- function(obs_data, pred_data, country_name, gender, year) {

  # Get gender label
  gender_label <- format_gender_label(gender)

  # Create indicator labels for faceting
  if (nrow(pred_data) > 0) {
    pred_data <- pred_data %>%
      mutate(Indicator_Label = format_indicator_label(Def_Type_Code))
  }
  if (nrow(obs_data) > 0) {
    obs_data <- obs_data %>%
      mutate(Indicator_Label = format_indicator_label(def_type_code))
  }

  p <- ggplot()

  # Add predictions if available
  if (nrow(pred_data) > 0) {
    # Confidence ribbon
    if ("lower_ci" %in% names(pred_data) && "upper_ci" %in% names(pred_data)) {
      p <- p + geom_ribbon(
        data = pred_data,
        aes(x = Age_Midpoint, ymin = lower_ci * 100, ymax = upper_ci * 100,
            fill = Model_Type),
        alpha = 0.15
      )
    }

    # Fitted line
    p <- p + geom_line(
      data = pred_data,
      aes(x = Age_Midpoint, y = Prevalence * 100, color = Model_Type),
      linewidth = 0.8
    )
  }

  # Add observed data if available
  if (nrow(obs_data) > 0) {
    # Age-band segments
    if ("start_age" %in% names(obs_data) && "end_age" %in% names(obs_data)) {
      p <- p + geom_segment(
        data = obs_data,
        aes(x = start_age, xend = end_age,
            y = plogis(prevalence) * 100, yend = plogis(prevalence) * 100),
        color = "black", linewidth = 0.4, alpha = 0.7
      )
    }

    # Observed points
    p <- p + geom_point(
      data = obs_data,
      aes(x = Age_Midpoint, y = plogis(prevalence) * 100),
      size = 1.5, shape = 21, fill = "white", stroke = 0.5
    )
  }

  # Transition marker
  p <- p +
    geom_vline(xintercept = TRANSITION_AGE, linetype = "dashed",
               color = "orange", alpha = 0.6, linewidth = 0.5)

  # Facet and styling
  p <- p +
    facet_wrap(~ Indicator_Label, ncol = 3, scales = "free_y") +
    scale_color_manual(values = MODEL_COLORS, name = "Model") +
    scale_fill_manual(values = MODEL_COLORS, name = "Model") +
    scale_x_continuous(limits = c(15, 85), breaks = seq(20, 80, 20)) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, NA)) +
    theme_who_pub() +
    labs(
      title = paste0("Age-Specific Prevalence: ", country_name),
      subtitle = paste0(gender_label, " | Year: ", year),
      x = "Age",
      y = "Prevalence (%)",
      caption = "Horizontal bars show survey age bands. Orange dashed line = spline/linear transition (65)."
    )

  return(p)
}

#' Create age curve plot - Overlay All Indicators
plot_age_curve_overlay <- function(obs_data, pred_data, country_name, gender, year) {

  gender_label <- format_gender_label(gender)

  # Add indicator labels
  if (nrow(pred_data) > 0) {
    pred_data <- pred_data %>%
      mutate(Indicator_Label = format_indicator_label(Def_Type_Code))
  }
  if (nrow(obs_data) > 0) {
    obs_data <- obs_data %>%
      mutate(Indicator_Label = format_indicator_label(def_type_code))
  }

  p <- ggplot()

  # Predicted lines (no ribbons - too cluttered)
  if (nrow(pred_data) > 0) {
    p <- p + geom_line(
      data = pred_data,
      aes(x = Age_Midpoint, y = Prevalence * 100, color = Indicator_Label),
      linewidth = 0.8
    )
  }

  # Observed points
  if (nrow(obs_data) > 0) {
    p <- p + geom_point(
      data = obs_data,
      aes(x = Age_Midpoint, y = plogis(prevalence) * 100, color = Indicator_Label),
      size = 1.2, alpha = 0.6
    )
  }

  # Transition marker
  p <- p +
    geom_vline(xintercept = TRANSITION_AGE, linetype = "dashed",
               color = "grey50", alpha = 0.5)

  # Styling
  p <- p +
    scale_color_manual(values = INDICATOR_COLORS, name = "Indicator",
                       labels = function(x) gsub("Current |Daily ", "", x)) +
    scale_x_continuous(limits = c(15, 85)) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, NA)) +
    theme_who_pub() +
    theme(legend.position = "right") +
    labs(
      title = paste0("Age-Specific Prevalence: ", country_name),
      subtitle = paste0(gender_label, " | Year: ", year, " | All Indicators"),
      x = "Age",
      y = "Prevalence (%)"
    )

  return(p)
}

#' Create all age curve plots for a specific country-gender-year
create_age_curve_plots <- function(country_code, gender, year,
                                    observed_data, predictions,
                                    country_name_mapping = NULL,
                                    base_dir = FIG_BASE_DIR) {

  # Get country name
  country_name <- if (!is.null(country_name_mapping) && country_code %in% names(country_name_mapping)) {
    tools::toTitleCase(country_name_mapping[country_code])
  } else {
    toupper(country_code)
  }

  # Create output directory
  out_dir <- file.path(base_dir, "age_curves", country_code, gender, year)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Filter observed data
  obs <- observed_data %>%
    filter(wb_country_abv == country_code,
           sex == gender,
           year == !!year)

  # Filter predictions
  pred <- predictions %>%
    filter(tolower(Country) == country_code | Country == toupper(country_code),
           tolower(Sex) == gender,
           Year == !!year)

  # Add Model_Type if missing
  if (!"Model_Type" %in% names(pred)) {
    pred$Model_Type <- "Global"
  }

  # Check if we have any data
  if (nrow(obs) == 0 && nrow(pred) == 0) {
    return(list(success = FALSE, message = "No data"))
  }

  # Create faceted plot
  tryCatch({
    p_faceted <- plot_age_curve_faceted(obs, pred, country_name, gender, year)
    ggsave(file.path(out_dir, "faceted_by_indicator.pdf"), p_faceted,
           width = 10, height = 7)
    ggsave(file.path(out_dir, "faceted_by_indicator.png"), p_faceted,
           width = 10, height = 7, dpi = 150)
  }, error = function(e) {
    warning(sprintf("Faceted plot failed for %s/%s/%s: %s", country_code, gender, year, e$message))
  })

  # Create overlay plot
  tryCatch({
    p_overlay <- plot_age_curve_overlay(obs, pred, country_name, gender, year)
    ggsave(file.path(out_dir, "overlay_all_indicators.pdf"), p_overlay,
           width = 10, height = 6)
    ggsave(file.path(out_dir, "overlay_all_indicators.png"), p_overlay,
           width = 10, height = 6, dpi = 150)
  }, error = function(e) {
    warning(sprintf("Overlay plot failed for %s/%s/%s: %s", country_code, gender, year, e$message))
  })

  return(list(success = TRUE, out_dir = out_dir))
}

# ============================================================================
# WEIGHTED TREND FUNCTIONS
# ============================================================================

#' Create trend plot - Faceted by Indicator
plot_trend_faceted <- function(trend_data, country_name, gender) {

  gender_label <- format_gender_label(gender)

  trend_data <- trend_data %>%
    mutate(Indicator_Label = format_indicator_label(Def_Type_Code))

  p <- ggplot(trend_data, aes(x = Year))

  # Endgame zone
  p <- p +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = ENDGAME_THRESHOLD,
             fill = "#E8F5E9", alpha = 0.6)

  # Confidence ribbon if available
  if ("lower_ci" %in% names(trend_data) && "upper_ci" %in% names(trend_data)) {
    p <- p + geom_ribbon(
      aes(ymin = lower_ci * 100, ymax = upper_ci * 100),
      fill = "#2166AC", alpha = 0.2
    )
  }

  # Mean line
  mean_col <- if ("weighted_mean" %in% names(trend_data)) "weighted_mean" else "Prevalence"
  p <- p + geom_line(
    aes(y = !!sym(mean_col) * 100),
    color = "#2166AC", linewidth = 0.8
  )

  # Key year markers
  p <- p +
    geom_vline(xintercept = BASE_YEAR, linetype = "dotted", color = "grey50") +
    geom_vline(xintercept = TARGET_YEAR, linetype = "dotted", color = "grey50") +
    geom_hline(yintercept = ENDGAME_THRESHOLD, linetype = "dashed",
               color = "#1B5E20", linewidth = 0.5)

  # Facet and styling
  p <- p +
    facet_wrap(~ Indicator_Label, ncol = 3, scales = "free_y") +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, NA)) +
    theme_who_pub() +
    labs(
      title = paste0("Prevalence Trends: ", country_name),
      subtitle = paste0(gender_label, " | All Indicators"),
      x = "Year",
      y = "Age-Standardized Prevalence (%)",
      caption = paste0("Green zone: <", ENDGAME_THRESHOLD, "% endgame. ",
                       "Dotted lines: 2010 baseline, 2025 target.")
    )

  return(p)
}

#' Create trend plot - Overlay All Indicators
plot_trend_overlay <- function(trend_data, country_name, gender) {

  gender_label <- format_gender_label(gender)

  trend_data <- trend_data %>%
    mutate(Indicator_Label = format_indicator_label(Def_Type_Code))

  mean_col <- if ("weighted_mean" %in% names(trend_data)) "weighted_mean" else "Prevalence"

  p <- ggplot(trend_data, aes(x = Year, color = Indicator_Label, fill = Indicator_Label))

  # Endgame zone
  p <- p +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = ENDGAME_THRESHOLD,
             fill = "#E8F5E9", alpha = 0.6)

  # Light ribbons if available
  if ("lower_ci" %in% names(trend_data) && "upper_ci" %in% names(trend_data)) {
    p <- p + geom_ribbon(
      aes(ymin = lower_ci * 100, ymax = upper_ci * 100),
      alpha = 0.1, color = NA
    )
  }

  # Trend lines
  p <- p + geom_line(aes(y = !!sym(mean_col) * 100), linewidth = 0.8)

  # Key markers
  p <- p +
    geom_vline(xintercept = BASE_YEAR, linetype = "dotted", color = "grey50") +
    geom_vline(xintercept = TARGET_YEAR, linetype = "dotted", color = "grey50") +
    geom_hline(yintercept = ENDGAME_THRESHOLD, linetype = "dashed",
               color = "#1B5E20", linewidth = 0.5)

  # Styling
  p <- p +
    scale_color_manual(values = INDICATOR_COLORS, name = "Indicator") +
    scale_fill_manual(values = INDICATOR_COLORS, name = "Indicator") +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, NA)) +
    theme_who_pub() +
    theme(legend.position = "right") +
    labs(
      title = paste0("Prevalence Trends: ", country_name),
      subtitle = paste0(gender_label, " | All Indicators Overlaid"),
      x = "Year",
      y = "Age-Standardized Prevalence (%)"
    )

  return(p)
}

#' Create single indicator trend plot
plot_trend_single <- function(trend_data, indicator, country_name, gender) {

  gender_label <- format_gender_label(gender)
  indicator_label <- format_indicator_label(indicator)

  data <- trend_data %>% filter(Def_Type_Code == indicator)

  if (nrow(data) == 0) return(NULL)

  mean_col <- if ("weighted_mean" %in% names(data)) "weighted_mean" else "Prevalence"

  p <- ggplot(data, aes(x = Year))

  # Endgame zone
  p <- p +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = ENDGAME_THRESHOLD,
             fill = "#E8F5E9", alpha = 0.6) +
    annotate("text", x = min(data$Year) + 1, y = ENDGAME_THRESHOLD / 2,
             label = "Endgame Zone", hjust = 0, size = 3, color = "#1B5E20")

  # Confidence ribbon
  if ("lower_ci" %in% names(data) && "upper_ci" %in% names(data)) {
    p <- p + geom_ribbon(
      aes(ymin = lower_ci * 100, ymax = upper_ci * 100),
      fill = INDICATOR_COLORS[indicator], alpha = 0.2
    )
  }

  # Trend line
  p <- p + geom_line(
    aes(y = !!sym(mean_col) * 100),
    color = INDICATOR_COLORS[indicator], linewidth = 1
  )

  # Key markers
  p <- p +
    geom_vline(xintercept = BASE_YEAR, linetype = "dotted", color = "grey50") +
    geom_vline(xintercept = TARGET_YEAR, linetype = "dotted", color = "grey50") +
    geom_hline(yintercept = ENDGAME_THRESHOLD, linetype = "dashed",
               color = "#1B5E20", linewidth = 0.5)

  # Styling
  p <- p +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, NA)) +
    theme_who_pub() +
    labs(
      title = paste0(indicator_label, ": ", country_name),
      subtitle = gender_label,
      x = "Year",
      y = "Age-Standardized Prevalence (%)",
      caption = "Shaded area: 95% credible interval"
    )

  return(p)
}

#' Create all trend plots for a country-gender
create_trend_plots <- function(country_code, gender, trend_data,
                                country_name_mapping = NULL,
                                base_dir = FIG_BASE_DIR) {

  # Get country name
  country_name <- if (!is.null(country_name_mapping) && country_code %in% names(country_name_mapping)) {
    tools::toTitleCase(country_name_mapping[country_code])
  } else {
    toupper(country_code)
  }

  # Create output directory
  out_dir <- file.path(base_dir, "weighted_trends", country_code, gender)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Filter data
  # Handle different column name conventions
  country_col <- if ("Country" %in% names(trend_data)) "Country" else
    if ("wb_country_abv" %in% names(trend_data)) "wb_country_abv" else "area"
  sex_col <- if ("Sex" %in% names(trend_data)) "Sex" else "sex"
  indicator_col <- if ("Def_Type_Code" %in% names(trend_data)) "Def_Type_Code" else "def_type_code"

  data <- trend_data %>%
    filter(tolower(!!sym(country_col)) == country_code,
           tolower(!!sym(sex_col)) == gender)

  # Standardize column names
  if (!"Def_Type_Code" %in% names(data) && indicator_col != "Def_Type_Code") {
    data <- data %>% rename(Def_Type_Code = !!sym(indicator_col))
  }

  if (nrow(data) == 0) {
    return(list(success = FALSE, message = "No trend data"))
  }

  # Faceted plot
  tryCatch({
    p_faceted <- plot_trend_faceted(data, country_name, gender)
    ggsave(file.path(out_dir, "all_indicators_faceted.pdf"), p_faceted,
           width = 10, height = 7)
    ggsave(file.path(out_dir, "all_indicators_faceted.png"), p_faceted,
           width = 10, height = 7, dpi = 150)
  }, error = function(e) {
    warning(sprintf("Trend faceted plot failed: %s", e$message))
  })

  # Overlay plot
  tryCatch({
    p_overlay <- plot_trend_overlay(data, country_name, gender)
    ggsave(file.path(out_dir, "all_indicators_overlay.pdf"), p_overlay,
           width = 10, height = 6)
    ggsave(file.path(out_dir, "all_indicators_overlay.png"), p_overlay,
           width = 10, height = 6, dpi = 150)
  }, error = function(e) {
    warning(sprintf("Trend overlay plot failed: %s", e$message))
  })

  # Individual indicator plots
  for (ind in unique(data$Def_Type_Code)) {
    tryCatch({
      p_single <- plot_trend_single(data, ind, country_name, gender)
      if (!is.null(p_single)) {
        ind_short <- gsub("_user_|_tobacco_|_product", "", ind)
        ggsave(file.path(out_dir, paste0(ind_short, "_trend.pdf")), p_single,
               width = 8, height = 5)
        ggsave(file.path(out_dir, paste0(ind_short, "_trend.png")), p_single,
               width = 8, height = 5, dpi = 150)
      }
    }, error = function(e) {
      warning(sprintf("Single trend plot failed for %s: %s", ind, e$message))
    })
  }

  return(list(success = TRUE, out_dir = out_dir))
}

# ============================================================================
# VALIDATION FUNCTIONS
# ============================================================================

#' Create global validation scatter plot
plot_validation_global <- function(validation_data) {

  # Calculate metrics
  r2 <- cor(validation_data$Observed, validation_data$Predicted, use = "complete.obs")^2
  rmse <- sqrt(mean((validation_data$Observed - validation_data$Predicted)^2, na.rm = TRUE))
  n <- sum(!is.na(validation_data$Observed) & !is.na(validation_data$Predicted))

  p <- ggplot(validation_data, aes(x = Observed, y = Predicted)) +
    # Perfect fit line
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B71C1C", linewidth = 0.8) +
    # Points
    geom_point(alpha = 0.3, size = 0.8, color = "#2166AC") +
    # Linear fit
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
    # Metrics annotation
    annotate("text", x = min(validation_data$Observed, na.rm = TRUE) + 2,
             y = max(validation_data$Predicted, na.rm = TRUE) * 0.95,
             label = sprintf("R\u00B2 = %.3f\nRMSE = %.2f%%\nn = %s",
                             r2, rmse, format(n, big.mark = ",")),
             hjust = 0, size = 3.5, fontface = "bold") +
    coord_equal() +
    theme_who_pub() +
    labs(
      title = "Model Validation: Observed vs Predicted Prevalence",
      subtitle = "Global APC Model - All Countries and Indicators",
      x = "Observed Prevalence (%)",
      y = "Predicted Prevalence (%)",
      caption = "Red dashed = perfect fit. Black = linear regression."
    )

  return(p)
}

#' Create validation plot faceted by age group
plot_validation_by_age <- function(validation_data) {

  # Create age groups
  validation_data <- validation_data %>%
    mutate(Age_Group = create_age_groups(Age))

  # Calculate R² per group
  r2_by_age <- validation_data %>%
    group_by(Age_Group) %>%
    summarise(
      r2 = cor(Observed, Predicted, use = "complete.obs")^2,
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(label = sprintf("R\u00B2=%.2f\nn=%d", r2, n))

  p <- ggplot(validation_data, aes(x = Observed, y = Predicted)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B71C1C", linewidth = 0.5) +
    geom_point(alpha = 0.2, size = 0.5, color = "#2166AC") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.4) +
    geom_text(data = r2_by_age, aes(label = label),
              x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, size = 2.5) +
    facet_wrap(~ Age_Group, ncol = 4) +
    coord_equal() +
    theme_who_pub() +
    labs(
      title = "Model Validation by Age Group",
      subtitle = "Spline fit quality across age bands",
      x = "Observed (%)",
      y = "Predicted (%)"
    )

  return(p)
}

#' Create validation plot faceted by indicator
plot_validation_by_indicator <- function(validation_data) {

  validation_data <- validation_data %>%
    mutate(Indicator_Label = format_indicator_label(Indicator))

  # Calculate R² per indicator
  r2_by_ind <- validation_data %>%
    group_by(Indicator_Label) %>%
    summarise(
      r2 = cor(Observed, Predicted, use = "complete.obs")^2,
      rmse = sqrt(mean((Observed - Predicted)^2, na.rm = TRUE)),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(label = sprintf("R\u00B2=%.3f\nRMSE=%.1f%%", r2, rmse))

  p <- ggplot(validation_data, aes(x = Observed, y = Predicted)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B71C1C", linewidth = 0.5) +
    geom_point(alpha = 0.2, size = 0.5, color = "#2166AC") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.4) +
    geom_text(data = r2_by_ind, aes(label = label),
              x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, size = 2.5) +
    facet_wrap(~ Indicator_Label, ncol = 3) +
    coord_equal() +
    theme_who_pub() +
    labs(
      title = "Model Validation by Indicator",
      subtitle = "Comparing fit quality across tobacco product types",
      x = "Observed (%)",
      y = "Predicted (%)"
    )

  return(p)
}

#' Create residual plot by age (spline quality check)
plot_residuals_by_age <- function(validation_data) {

  validation_data <- validation_data %>%
    mutate(Residual = Observed - Predicted)

  # Calculate summary stats by age
  residual_summary <- validation_data %>%
    group_by(Age = round(Age)) %>%
    summarise(
      mean_resid = mean(Residual, na.rm = TRUE),
      sd_resid = sd(Residual, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )

  p <- ggplot(validation_data, aes(x = Age, y = Residual)) +
    # Zero line
    geom_hline(yintercept = 0, linetype = "dashed", color = "#1B5E20", linewidth = 0.5) +
    # Transition marker
    geom_vline(xintercept = TRANSITION_AGE, linetype = "dashed",
               color = "orange", alpha = 0.6) +
    annotate("text", x = TRANSITION_AGE + 1, y = max(validation_data$Residual, na.rm = TRUE) * 0.9,
             label = "Spline\u2192Linear", hjust = 0, size = 2.5, color = "orange") +
    # Points
    geom_point(alpha = 0.15, size = 0.5, color = "#2166AC") +
    # LOESS smooth
    geom_smooth(method = "loess", se = TRUE, color = "#B71C1C", fill = "#FFCDD2",
                linewidth = 0.8, alpha = 0.3) +
    # Summary ribbon
    geom_ribbon(data = residual_summary,
                aes(x = Age, ymin = mean_resid - sd_resid, ymax = mean_resid + sd_resid),
                alpha = 0.1, fill = "grey50") +
    scale_x_continuous(limits = c(15, 85)) +
    theme_who_pub() +
    labs(
      title = "Residuals by Age: Spline Fit Quality",
      subtitle = "Residual = Observed - Predicted. Ideal: flat at zero.",
      x = "Age",
      y = "Residual (Observed - Predicted, %)",
      caption = "Red line: LOESS smooth. Grey band: \u00B11 SD. Orange: spline/linear transition."
    )

  return(p)
}

#' Prepare validation data by merging observed and predicted
prepare_validation_data <- function(observed_data, predictions) {

  # Prepare observed
  obs <- observed_data %>%
    mutate(
      Year = as.numeric(year),
      Age = Age_Midpoint,
      Country = wb_country_abv,
      Sex = sex,
      Indicator = def_type_code,
      Observed = plogis(prevalence) * 100
    ) %>%
    select(Year, Age, Country, Sex, Indicator, Observed)

  # Prepare predictions
  pred <- predictions %>%
    mutate(
      Year = as.numeric(Year),
      Age = Age_Midpoint,
      Country = tolower(Country),
      Sex = tolower(Sex),
      Indicator = Def_Type_Code,
      Predicted = if(max(Prevalence, na.rm = TRUE) <= 1) Prevalence * 100 else Prevalence
    ) %>%
    select(Year, Age, Country, Sex, Indicator, Predicted)

  # Merge
  validation <- obs %>%
    inner_join(pred, by = c("Year", "Age", "Country", "Sex", "Indicator"))

  return(validation)
}

#' Create all validation plots
create_validation_plots <- function(observed_data, predictions, base_dir = FIG_BASE_DIR) {

  cat("  Creating validation plots...\n")

  # Prepare data
  validation <- prepare_validation_data(observed_data, predictions)

  if (nrow(validation) == 0) {
    warning("No matched data for validation plots")
    return(list(success = FALSE))
  }

  cat(sprintf("    Matched %d observation-prediction pairs\n", nrow(validation)))

  # Create directories
  val_dir <- file.path(base_dir, "validation")
  dir.create(val_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(val_dir, "by_age_group"), showWarnings = FALSE)
  dir.create(file.path(val_dir, "by_indicator"), showWarnings = FALSE)

  # Global scatter
  tryCatch({
    p_global <- plot_validation_global(validation)
    ggsave(file.path(val_dir, "global_scatter.pdf"), p_global, width = 8, height = 8)
    ggsave(file.path(val_dir, "global_scatter.png"), p_global, width = 8, height = 8, dpi = 150)
    cat("    Saved: global_scatter\n")
  }, error = function(e) warning(e$message))

  # By age group
  tryCatch({
    p_age <- plot_validation_by_age(validation)
    ggsave(file.path(val_dir, "by_age_group", "validation_age_facets.pdf"), p_age,
           width = 12, height = 8)
    ggsave(file.path(val_dir, "by_age_group", "validation_age_facets.png"), p_age,
           width = 12, height = 8, dpi = 150)
    cat("    Saved: by_age_group\n")
  }, error = function(e) warning(e$message))

  # By indicator
  tryCatch({
    p_ind <- plot_validation_by_indicator(validation)
    ggsave(file.path(val_dir, "by_indicator", "validation_indicator_facets.pdf"), p_ind,
           width = 10, height = 7)
    ggsave(file.path(val_dir, "by_indicator", "validation_indicator_facets.png"), p_ind,
           width = 10, height = 7, dpi = 150)
    cat("    Saved: by_indicator\n")
  }, error = function(e) warning(e$message))

  # Residuals by age
  tryCatch({
    p_resid <- plot_residuals_by_age(validation)
    ggsave(file.path(val_dir, "residuals_by_age.pdf"), p_resid, width = 10, height = 6)
    ggsave(file.path(val_dir, "residuals_by_age.png"), p_resid, width = 10, height = 6, dpi = 150)
    cat("    Saved: residuals_by_age\n")
  }, error = function(e) warning(e$message))

  return(list(success = TRUE, n_pairs = nrow(validation)))
}

# ============================================================================
# MASTER LOOP FUNCTIONS
# ============================================================================

#' Generate all age curve figures
generate_age_curve_figures <- function(clean_data, predictions,
                                        country_name_mapping = NULL,
                                        base_dir = FIG_BASE_DIR) {

  cat("\n========== GENERATING AGE CURVE FIGURES ==========\n")

  countries <- unique(clean_data$wb_country_abv)
  genders <- c("males", "females")

  total <- length(countries)
  success_count <- 0
  error_count <- 0

  for (i in seq_along(countries)) {
    country <- countries[i]

    for (gender in genders) {
      # Get years with observed data
      obs_years <- clean_data %>%
        filter(wb_country_abv == country, sex == gender) %>%
        pull(year) %>%
        unique() %>%
        sort()

      for (year in obs_years) {
        result <- tryCatch({
          create_age_curve_plots(
            country_code = country,
            gender = gender,
            year = year,
            observed_data = clean_data,
            predictions = predictions,
            country_name_mapping = country_name_mapping,
            base_dir = base_dir
          )
        }, error = function(e) {
          list(success = FALSE, message = e$message)
        })

        if (result$success) {
          success_count <- success_count + 1
        } else {
          error_count <- error_count + 1
        }
      }
    }

    # Progress
    cat(sprintf("\r  Progress: %d/%d countries (plots: %d success, %d failed)",
                i, total, success_count, error_count))
  }

  cat("\n")
  cat(sprintf("  Age curves complete: %d plots, %d errors\n", success_count, error_count))

  return(list(success = success_count, errors = error_count))
}

#' Generate all trend figures
generate_trend_figures <- function(trend_data, country_name_mapping = NULL,
                                    base_dir = FIG_BASE_DIR) {

  cat("\n========== GENERATING TREND FIGURES ==========\n")

  # Get countries from trend data
  country_col <- if ("Country" %in% names(trend_data)) "Country" else
    if ("wb_country_abv" %in% names(trend_data)) "wb_country_abv" else "area"

  countries <- unique(tolower(trend_data[[country_col]]))
  genders <- c("males", "females")

  total <- length(countries)
  success_count <- 0
  error_count <- 0

  for (i in seq_along(countries)) {
    country <- countries[i]

    for (gender in genders) {
      result <- tryCatch({
        create_trend_plots(
          country_code = country,
          gender = gender,
          trend_data = trend_data,
          country_name_mapping = country_name_mapping,
          base_dir = base_dir
        )
      }, error = function(e) {
        list(success = FALSE, message = e$message)
      })

      if (result$success) {
        success_count <- success_count + 1
      } else {
        error_count <- error_count + 1
      }
    }

    cat(sprintf("\r  Progress: %d/%d countries", i, total))
  }

  cat("\n")
  cat(sprintf("  Trends complete: %d country-gender pairs, %d errors\n",
              success_count, error_count))

  return(list(success = success_count, errors = error_count))
}

#' Generate all figures (master function)
generate_all_strata_figures <- function(clean_data,
                                         predictions,
                                         trend_data = NULL,
                                         country_name_mapping = NULL,
                                         base_dir = FIG_BASE_DIR) {

  cat("\n")
  cat("###########################################################################\n")
  cat("#                  GENERATING ALL STRATA-LEVEL FIGURES                    #\n")
  cat("###########################################################################\n")

  start_time <- Sys.time()

  # Create directory structure
  create_figure_dirs(base_dir)

  results <- list()

  # 1. Validation plots
  results$validation <- create_validation_plots(clean_data, predictions, base_dir)

  # 2. Age curve figures
  results$age_curves <- generate_age_curve_figures(
    clean_data, predictions, country_name_mapping, base_dir
  )

  # 3. Trend figures (if trend data available)
  if (!is.null(trend_data) && nrow(trend_data) > 0) {
    results$trends <- generate_trend_figures(trend_data, country_name_mapping, base_dir)
  } else {
    cat("\n  Skipping trend figures (no trend_data provided)\n")
  }

  # Summary
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")

  cat("\n")
  cat("###########################################################################\n")
  cat("#                           GENERATION SUMMARY                            #\n")
  cat("###########################################################################\n")
  cat(sprintf("  Time elapsed: %.1f minutes\n", as.numeric(elapsed)))
  cat(sprintf("  Validation pairs: %d\n", results$validation$n_pairs))
  cat(sprintf("  Age curve plots: %d success, %d errors\n",
              results$age_curves$success, results$age_curves$errors))
  if (!is.null(results$trends)) {
    cat(sprintf("  Trend plots: %d success, %d errors\n",
                results$trends$success, results$trends$errors))
  }
  cat("\n")

  return(results)
}

# ============================================================================
# STANDALONE EXECUTION
# ============================================================================

if (exists("clean_data") && exists("final_ac_predictions")) {
  cat("\n  Data available - ready to generate figures.\n")
  cat("  Call: generate_all_strata_figures(clean_data, final_ac_predictions, ...)\n")
} else {
  cat("\n  No data loaded. Load clean_data and predictions first.\n")
}

cat("\n  Module 13 loaded.\n")
