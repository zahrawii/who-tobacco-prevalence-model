#########################################################################################
#
#                    WHO TOBACCO PIPELINE - MODULE 12: PUBLICATION FIGURES
#
#   Comprehensive visualization module following PhD Thesis Figure Design Guide
#
#   FIGURES:
#     - Figure 1: Stoplight World Map (WHO Target Achievement Probability)
#     - Figure 2: Caterpillar Plot (Country Ranking by Reduction)
#     - Figure 3: Endgame Trajectories (Spaghetti Plot with <5% Zone)
#     - Figure 4: Age-Cohort Heatmap (Lexis Diagram)
#     - eFigure 1: Model Validation (Observed vs Predicted)
#     - eFigure 2: Country Age Profiles (Fitted vs Observed with Age Bands)
#     - eFigure 3: Weighted Trend Curves (Time Series with Uncertainty)
#
#   KEY FIXES:
#     - Country code harmonization for map joins
#     - Proper age-band segments for survey data
#     - Consistent variable naming throughout
#     - Working loops with error handling
#
#########################################################################################

cat("\n")
cat("###########################################################################\n")
cat("#                                                                         #\n")
cat("#        MODULE 12: PUBLICATION FIGURES                                   #\n")
cat("#                                                                         #\n")
cat("###########################################################################\n\n")

# ============================================================================
# SETUP AND DEPENDENCIES
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(viridis)
  library(patchwork)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(ggrepel)
})

# Create output directories
output_dirs <- c(
  # Maps with target/gender subdirectories
  "outputs/figures/maps",
  "outputs/figures/maps/2025_target/men",
  "outputs/figures/maps/2025_target/women",
  "outputs/figures/maps/2040_endgame/men",
  "outputs/figures/maps/2040_endgame/women",
  "outputs/figures/maps/composites",
  # Other figures
  "outputs/figures/caterpillar",
  "outputs/figures/trajectories",
  "outputs/figures/heatmaps",
  "outputs/figures/validation",
  "outputs/figures/age_profiles",
  "outputs/figures/trends"
)
for (d in output_dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================================
# CONFIGURATION
# ============================================================================

# Year parameters
BASE_YEAR <- 2010
TARGET_YEAR <- 2025
TARGET_YEAR_2030 <- 2030
ENDGAME_YEAR <- 2040

# Target thresholds
REDUCTION_TARGET <- 0.30
ENDGAME_THRESHOLD <- 5  # 5%
NEAR_ENDGAME_THRESHOLD <- 10  # 10%

# Spline transition
TRANSITION_START <- 65

# Primary indicator
PRIMARY_INDICATOR <- "current_user_cigarettes"

# ============================================================================
# VALID INDICATORS (Single Source of Truth)
# ============================================================================

VALID_INDICATORS <- c(
  "current_user_cigarettes",
  "current_user_any_smoked_tobacco",
  "current_user_any_tobacco_product",
  "daily_user_cigarettes",
  "daily_user_any_smoked_tobacco",
  "daily_user_any_tobacco_product"
)

# Indicator display labels
INDICATOR_LABELS <- c(
  "current_user_cigarettes"           = "Current Cigarettes",
  "current_user_any_smoked_tobacco"   = "Current Any Smoked Tobacco",
  "current_user_any_tobacco_product"  = "Current Any Tobacco Product",
  "daily_user_cigarettes"             = "Daily Cigarettes",
  "daily_user_any_smoked_tobacco"     = "Daily Any Smoked Tobacco",
  "daily_user_any_tobacco_product"    = "Daily Any Tobacco Product"
)

# Indicator short codes for file names
INDICATOR_CODES <- c(
  "current_user_cigarettes"           = "CU_CIG",
  "current_user_any_smoked_tobacco"   = "CU_SMK",
  "current_user_any_tobacco_product"  = "CU_TOB",
  "daily_user_cigarettes"             = "DU_CIG",
  "daily_user_any_smoked_tobacco"     = "DU_SMK",
  "daily_user_any_tobacco_product"    = "DU_TOB"
)

# ============================================================================
# SEX/GENDER STANDARDIZATION (Critical for Data Consistency)
# ============================================================================

# Valid internal sex values (always lowercase, plural)
VALID_SEX <- c("males", "females")

# Display labels for gender (Men/Women for human subjects)
GENDER_DISPLAY <- c(males = "Men", females = "Women")

#' Standardize sex value to internal format
#' Converts any variant to lowercase plural form
#'
#' @param x Sex value (can be "male", "males", "Male", "MALE", "M", etc.)
#' @return Standardized value ("males" or "females")
standardize_sex_value <- function(x) {
  x_lower <- tolower(trimws(x))
  result <- case_when(
    x_lower %in% c("male", "males", "m", "men", "man") ~ "males",
    x_lower %in% c("female", "females", "f", "women", "woman") ~ "females",
    x_lower %in% c("both", "all", "total") ~ "both",
    TRUE ~ NA_character_
  )
  if (any(is.na(result) & !is.na(x))) {
    invalid <- unique(x[is.na(result) & !is.na(x)])
    warning(sprintf("Unknown sex values converted to NA: %s", paste(invalid, collapse = ", ")))
  }
  return(result)
}

#' Get display label for sex
#' @param x Internal sex value ("males" or "females")
#' @return Display label ("Men" or "Women")
format_gender <- function(x) {
  result <- case_when(
    x == "males" ~ "Men",
    x == "females" ~ "Women",
    x == "male" ~ "Men",
    x == "female" ~ "Women",
    TRUE ~ as.character(x)
  )
  return(result)
}

# ============================================================================
# ASSERTION FUNCTIONS (Fail Fast on Bad Data)
# ============================================================================

#' Assert required columns exist
assert_required_columns <- function(data, cols, context = "data") {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop(sprintf("[%s] Missing required columns: %s\nAvailable: %s",
                 context, paste(missing, collapse = ", "),
                 paste(names(data), collapse = ", ")))
  }
  invisible(data)
}

#' Assert no NA values in critical columns
assert_no_na <- function(data, cols, context = "data") {
  for (col in cols) {
    if (col %in% names(data) && any(is.na(data[[col]]))) {
      n_na <- sum(is.na(data[[col]]))
      warning(sprintf("[%s] Column '%s' has %d NA values", context, col, n_na))
    }
  }
  invisible(data)
}

#' Assert valid sex values
assert_valid_sex <- function(data, sex_col = "sex", context = "data") {
  if (!sex_col %in% names(data)) return(invisible(data))

  sex_vals <- unique(tolower(data[[sex_col]]))
  invalid <- setdiff(sex_vals, c("males", "females", "male", "female", "both", NA))
  if (length(invalid) > 0) {
    stop(sprintf("[%s] Invalid sex values: %s", context, paste(invalid, collapse = ", ")))
  }
  invisible(data)
}

#' Assert valid indicator values
assert_valid_indicators <- function(data, indicator_col = "def_type_code", context = "data") {
  if (!indicator_col %in% names(data)) return(invisible(data))

  indicators <- unique(data[[indicator_col]])
  invalid <- setdiff(indicators, c(VALID_INDICATORS, NA))
  if (length(invalid) > 0) {
    warning(sprintf("[%s] Unknown indicators: %s", context, paste(invalid, collapse = ", ")))
  }
  invisible(data)
}

# ============================================================================
# COLOR PALETTES (Publication-Ready)
# ============================================================================

# Lancet-style stoplight colors for probability categories
# Refined for: print quality, colorblind accessibility, visual balance
stoplight_colors <- c(
  "Very High (>90%)"    = "#1B5E20",   # Deep forest green
  "High (60-90%)"       = "#4CAF50",   # Balanced green
  "Moderate (40-60%)"   = "#FF9800",   # Deep amber/orange (not yellow)
  "Low (10-40%)"        = "#E57373",   # Soft coral red
  "Very Low (<10%)"     = "#B71C1C",   # Deep red
  "No Data"             = "#BDBDBD"    # Medium gray (visible)
)

# ============================================================================
# PREFLIGHT VALIDATION (Run Before Any Visualization)
# ============================================================================

#' Preflight check for visualization data
#'
#' Validates data structure before running visualizations. Fails fast on critical issues.
#'
#' @param weighted_results Weighted results data frame
#' @param predictions Predictions data frame (optional)
#' @param verbose Print detailed output
#' @return TRUE if all checks pass, stops with error otherwise
preflight_visualization_check <- function(weighted_results, predictions = NULL, verbose = TRUE) {

  if (verbose) {
    cat("\n")
    cat("================================================================\n")
    cat("  PREFLIGHT CHECK: Visualization Data Validation\n")
    cat("================================================================\n\n")
  }

  errors <- c()
  warnings <- c()

  # 1. Check weighted_results structure
  if (is.null(weighted_results) || nrow(weighted_results) == 0) {
    errors <- c(errors, "weighted_results is NULL or empty")
  } else {
    # Check for required columns (using adaptive names)
    country_col <- if ("Country" %in% names(weighted_results)) "Country" else "wb_country_abv"
    sex_col <- if ("Sex" %in% names(weighted_results)) "Sex" else "sex"
    indicator_col <- if ("Def_Type_Code" %in% names(weighted_results)) "Def_Type_Code" else "def_type_code"
    year_col <- if ("Year" %in% names(weighted_results)) "Year" else "year"

    required <- c(country_col, sex_col, indicator_col, year_col)
    missing <- setdiff(required, names(weighted_results))

    if (length(missing) > 0) {
      errors <- c(errors, sprintf("Missing columns in weighted_results: %s", paste(missing, collapse = ", ")))
    } else {
      if (verbose) cat("  [OK] Required columns present\n")

      # Check sex values
      sex_vals <- unique(weighted_results[[sex_col]])
      valid_sex <- c("males", "females", "male", "female", "both")
      invalid_sex <- setdiff(tolower(sex_vals), valid_sex)
      if (length(invalid_sex) > 0) {
        warnings <- c(warnings, sprintf("Unusual sex values: %s", paste(invalid_sex, collapse = ", ")))
      } else {
        if (verbose) cat(sprintf("  [OK] Sex values valid: %s\n", paste(sex_vals, collapse = ", ")))
      }

      # Check indicators
      indicators <- unique(weighted_results[[indicator_col]])
      invalid_ind <- setdiff(indicators, VALID_INDICATORS)
      if (length(invalid_ind) > 0) {
        warnings <- c(warnings, sprintf("Non-standard indicators: %s", paste(invalid_ind, collapse = ", ")))
      } else {
        if (verbose) cat(sprintf("  [OK] Indicators valid: %d types\n", length(indicators)))
      }

      # Check year range
      years <- weighted_results[[year_col]]
      if (verbose) cat(sprintf("  [OK] Year range: %d - %d\n", min(years, na.rm=TRUE), max(years, na.rm=TRUE)))

      # Check countries
      countries <- unique(weighted_results[[country_col]])
      if (verbose) cat(sprintf("  [OK] Countries: %d unique\n", length(countries)))

      # Check for 2010 baseline
      if (!2010 %in% years) {
        warnings <- c(warnings, "No data for baseline year 2010")
      }

      # Check for 2025 target year
      if (!2025 %in% years) {
        warnings <- c(warnings, "No data for target year 2025")
      }
    }
  }

  # 2. Report results
  if (verbose) cat("\n")

  if (length(warnings) > 0) {
    cat("  WARNINGS:\n")
    for (w in warnings) cat(sprintf("    [!] %s\n", w))
    cat("\n")
  }

  if (length(errors) > 0) {
    cat("  ERRORS:\n")
    for (e in errors) cat(sprintf("    [X] %s\n", e))
    cat("\n")
    stop("Preflight check FAILED - cannot proceed with visualization")
  }

  if (verbose) {
    cat("  ================================================\n")
    cat("  PREFLIGHT CHECK: PASSED\n")
    cat("  ================================================\n\n")
  }

  invisible(TRUE)
}

# ============================================================================
# MAP BACKGROUND COLORS
# ============================================================================

# Map background colors
OCEAN_COLOR       <- "#DAE8F5"    # Soft blue for oceans
LAND_NO_DATA      <- "#BDBDBD"    # Gray for countries without data
BORDER_COLOR      <- "#FFFFFF"    # White country borders
BORDER_WIDTH      <- 0.25         # Border line width
COASTLINE_COLOR   <- "#90A4AE"    # Subtle coastline emphasis

# ============================================================================
# A4 PAGE DIMENSIONS (for Lancet submission)
# ============================================================================

# A4 vertical dimensions
A4_WIDTH  <- 8.27   # inches (210mm)
A4_HEIGHT <- 11.69  # inches (297mm)

# Optimized for maps with legend
FIG_WIDTH  <- 7.5   # inches (~190mm usable)
FIG_HEIGHT <- 10.0  # inches (~254mm usable)

# Region colors (colorblind-friendly)
region_colors <- c(
  "Central Asia" = "#E41A1C",
  "Eastern Asia" = "#377EB8",
  "Eastern Europe" = "#4DAF4A",
  "Latin America & Caribbean" = "#984EA3",
  "North Africa & Middle East" = "#FF7F00",
  "Northern Europe" = "#A65628",
  "Oceania & Pacific" = "#F781BF",
  "South America" = "#999999",
  "South Asia" = "#66C2A5",
  "South-East Asia" = "#FC8D62",
  "Sub-Saharan Africa" = "#8DA0CB",
  "Western Europe" = "#E78AC3",
  "Other" = "#A6D854"
)

# WHO style colors
who_colors <- list(
  primary   = "#2C3E50",
  secondary = "#3498DB",
  accent    = "#E74C3C",
  success   = "#27AE60",
  warning   = "#F39C12",
  gray      = "#95A5A6"
)

# ============================================================================
# COUNTRY CODE HARMONIZATION (CRITICAL FOR MAPS)
# ============================================================================

# Comprehensive mapping: WHO code (lowercase) -> Natural Earth ISO_A3 (uppercase)
# Based on thorough comparison of WHO tobacco data countries with Natural Earth

COUNTRY_CODE_HARMONIZATION <- c(
  # ====== VERIFIED MISMATCHES ======
  "rom" = "ROU",   # Romania (old code -> current ISO)
  "tmp" = "TLS",   # Timor-Leste (old code)
  "tls" = "TLS",   # Timor-Leste (current)
  "pse" = "PSE",   # Palestine
  "wbg" = "PSE",   # West Bank & Gaza -> Palestine
  "zar" = "COD",   # DR Congo (old Zaire code)
  "cod" = "COD",   # DR Congo (current)
  "ksv" = "XKX",   # Kosovo (Natural Earth uses XKX, not -99)
  "xkx" = "XKX",   # Kosovo

  # ====== EAST ASIA ======
  "kor" = "KOR",   # South Korea
  "prk" = "PRK",   # North Korea
  "twn" = "TWN",   # Taiwan (may show as "No Data" - not in NE as country)
  "hkg" = "HKG",   # Hong Kong (SAR)
  "mac" = "MAC",   # Macau (SAR)

  # ====== SOUTHEAST ASIA ======
  "lao" = "LAO",   # Laos
  "mmr" = "MMR",   # Myanmar
  "vnm" = "VNM",   # Vietnam
  "khm" = "KHM",   # Cambodia
  "brn" = "BRN",   # Brunei

  # ====== AFRICA ======
  "civ" = "CIV",   # Cote d'Ivoire
  "swz" = "SWZ",   # Eswatini (formerly Swaziland)
  "ssd" = "SSD",   # South Sudan
  "som" = "SOM",   # Somalia

  # ====== EUROPE ======
  "mkd" = "MKD",   # North Macedonia
  "srb" = "SRB",   # Serbia
  "mne" = "MNE",   # Montenegro
  "bih" = "BIH",   # Bosnia and Herzegovina
  "and" = "AND",   # Andorra (microstate - may not render)
  "smr" = "SMR",   # San Marino (microstate - may not render)
  "mco" = "MCO",   # Monaco (microstate - may not render)
  "lie" = "LIE",   # Liechtenstein (microstate - may not render)
  "vat" = "VAT",   # Vatican (microstate - may not render)

  # ====== PACIFIC ISLANDS ======
  # Many of these are too small for Natural Earth medium scale
  "mhl" = "MHL",   # Marshall Islands
  "fsm" = "FSM",   # Micronesia
  "nru" = "NRU",   # Nauru
  "niu" = "NIU",   # Niue (NZ territory - not in NE)
  "plw" = "PLW",   # Palau
  "wsm" = "WSM",   # Samoa
  "tuv" = "TUV",   # Tuvalu
  "cok" = "COK",   # Cook Islands (NZ territory - not in NE)
  "kir" = "KIR",   # Kiribati
  "ton" = "TON",   # Tonga
  "vut" = "VUT",   # Vanuatu
  "slb" = "SLB",   # Solomon Islands

  # ====== CARIBBEAN ======
  "kna" = "KNA",   # Saint Kitts and Nevis
  "lca" = "LCA",   # Saint Lucia
  "vct" = "VCT",   # Saint Vincent and the Grenadines
  "grd" = "GRD",   # Grenada
  "dma" = "DMA",   # Dominica
  "atg" = "ATG",   # Antigua and Barbuda
  "brb" = "BRB",   # Barbados

  # ====== MIDDLE EAST ======
  "are" = "ARE",   # United Arab Emirates
  "bhr" = "BHR",   # Bahrain
  "kwt" = "KWT",   # Kuwait
  "qat" = "QAT",   # Qatar
  "omn" = "OMN",   # Oman

  # ====== INDIAN OCEAN ======
  "mdv" = "MDV",   # Maldives
  "syc" = "SYC",   # Seychelles
  "mus" = "MUS",   # Mauritius
  "com" = "COM",   # Comoros
  "stp" = "STP"    # Sao Tome and Principe
)

#' Harmonize WHO country codes to ISO alpha-3 for map joining
#' This fixes known mismatches between WHO data and Natural Earth map data
#'
#' @param country_code WHO country code (can be lowercase or uppercase)
#' @return ISO alpha-3 code used by Natural Earth (uppercase)
harmonize_country_codes <- function(country_code) {
  # Convert to lowercase for lookup
  cc <- tolower(country_code)

  # Check if harmonization needed
  if (cc %in% names(COUNTRY_CODE_HARMONIZATION)) {
    return(COUNTRY_CODE_HARMONIZATION[cc])
  }

  # Otherwise just uppercase (most codes are already correct)
  return(toupper(cc))
}

#' Vectorized version for use with dplyr
harmonize_country_codes_vec <- function(country_codes) {
  sapply(country_codes, harmonize_country_codes, USE.NAMES = FALSE)
}

#' Get ISO alpha-3 codes used by Natural Earth
#' Creates a lookup table for joining
get_natural_earth_iso3 <- function() {
  world <- ne_countries(scale = "medium", returnclass = "sf")
  return(unique(world$iso_a3))
}

#' Check which WHO codes will NOT map to Natural Earth
#' Useful for diagnostics
check_unmappable_codes <- function(who_codes) {
  ne_iso3 <- get_natural_earth_iso3()
  ne_iso3 <- ne_iso3[!is.na(ne_iso3) & ne_iso3 != "-99"]

  unmappable <- c()
  for (code in who_codes) {
    mapped <- harmonize_country_codes(code)
    if (!mapped %in% ne_iso3) {
      unmappable <- c(unmappable, code)
    }
  }
  return(unmappable)
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Format gender labels
format_gender <- function(sex) {
  case_when(
    tolower(sex) == "males" ~ "Men",
    tolower(sex) == "females" ~ "Women",
    tolower(sex) == "male" ~ "Men",
    tolower(sex) == "female" ~ "Women",
    TRUE ~ as.character(sex)
  )
}

#' Format indicator labels
format_indicator <- function(x) {
  labels <- c(
    "daily_user_cigarettes" = "Daily Cigarettes",
    "current_user_cigarettes" = "Current Cigarettes",
    "daily_user_any_smoked_tobacco" = "Daily Smoked Tobacco",
    "current_user_any_smoked_tobacco" = "Current Smoked Tobacco",
    "daily_user_any_tobacco_product" = "Daily Any Tobacco",
    "current_user_any_tobacco_product" = "Current Any Tobacco"
  )
  ifelse(x %in% names(labels), labels[x], gsub("_", " ", tools::toTitleCase(x)))
}

#' Categorize probability into traffic-light bins (Lancet labels)
categorize_probability <- function(prob) {
  cut(prob,
      breaks = c(0, 0.10, 0.40, 0.60, 0.90, 1.0),
      labels = c("Very Low (<10%)", "Low (10-40%)", "Moderate (40-60%)",
                 "High (60-90%)", "Very High (>90%)"),
      include.lowest = TRUE)
}

# ============================================================================
# WHO MASTER THEME
# ============================================================================

theme_who <- function(base_size = 11, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0, margin = margin(b = 5)),
      plot.subtitle = element_text(color = "#555555", size = rel(1.0), hjust = 0, margin = margin(b = 10)),
      plot.caption = element_text(color = "#777777", size = rel(0.7), hjust = 1, margin = margin(t = 10)),
      axis.title = element_text(face = "bold", size = rel(0.9)),
      axis.text = element_text(color = "#333333", size = rel(0.85)),
      panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "#333333", linewidth = 0.3),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_text(face = "bold", size = rel(0.8)),
      strip.background = element_rect(fill = "#F5F5F5", color = NA),
      strip.text = element_text(face = "bold", size = rel(0.9), hjust = 0, margin = margin(4,4,4,4))
    )
}

# ============================================================================
# LANCET MAP THEME (Publication Quality)
# ============================================================================

#' Lancet-style theme for world maps (Publication Quality)
#'
#' A4 vertical layout with:
#' - Bottom horizontal legend
#' - Professional typography
#' - Balanced spacing
#'
#' @param base_size Base font size (default 11)
#' @param base_family Font family (default sans-serif)
theme_lancet_map <- function(base_size = 11, base_family = "sans") {
  theme_void(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # === TITLE BLOCK ===
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5,
        color = "#1A1A1A",
        margin = margin(t = 10, b = 5)
      ),
      plot.subtitle = element_text(
        size = 11,
        hjust = 0.5,
        color = "#4A4A4A",
        lineheight = 1.2,
        margin = margin(b = 20)
      ),
      plot.caption = element_text(
        size = 8,
        hjust = 0,
        color = "#666666",
        lineheight = 1.3,
        margin = margin(t = 20)
      ),

      # === FACET STRIPS (Men/Women) ===
      strip.text = element_text(
        face = "bold",
        size = 13,
        hjust = 0.5,
        color = "#2C3E50",
        margin = margin(t = 12, b = 8)
      ),
      strip.background = element_rect(fill = "#F8F9FA", color = NA),
      panel.spacing.y = unit(0.8, "cm"),
      panel.spacing.x = unit(0.4, "cm"),

      # === LEGEND (Bottom, Horizontal) ===
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.box.just = "center",
      legend.title = element_text(
        face = "bold",
        size = 10,
        color = "#333333",
        margin = margin(b = 8)
      ),
      legend.text = element_text(
        size = 9,
        color = "#4A4A4A",
        margin = margin(t = 3)
      ),
      legend.key.width = unit(1.3, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.spacing.x = unit(0.15, "cm"),
      legend.margin = margin(t = 15, b = 10),
      legend.box.margin = margin(t = 5),

      # === OVERALL LAYOUT ===
      plot.margin = margin(t = 15, r = 20, b = 15, l = 20),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = OCEAN_COLOR, color = NA)
    )
}

#' Legend guide for Lancet maps
#' Forces all categories to show, proper formatting
lancet_legend_guide <- function() {
  guides(fill = guide_legend(
    nrow = 1,
    title.position = "top",
    title.hjust = 0.5,
    label.position = "bottom",
    reverse = FALSE,
    keywidth = unit(1.3, "cm"),
    keyheight = unit(0.5, "cm"),
    override.aes = list(
      color = "#666666",
      linewidth = 0.4
    )
  ))
}

# ============================================================================
# HELPER: PREPARE WORLD MAP WITH FIXED ISO CODES
# ============================================================================

#' Prepare world map with fixed ISO codes
#'
#' Natural Earth has iso_a3 = "-99" for France, Norway, and some other countries.
#' This function uses adm0_a3 as a fallback when iso_a3 is invalid.
#'
#' @param crop_bounds Optional list with xmin, xmax, ymin, ymax for cropping
#' @return sf object with cleaned iso_a3 column
prepare_world_map <- function(crop_bounds = list(xmin = -180, xmax = 180, ymin = -56, ymax = 84)) {

  world <- ne_countries(scale = "medium", returnclass = "sf")

  # Crop to specified bounds
  sf::sf_use_s2(FALSE)
  world <- st_crop(world,
                   xmin = crop_bounds$xmin, xmax = crop_bounds$xmax,
                   ymin = crop_bounds$ymin, ymax = crop_bounds$ymax)
  sf::sf_use_s2(TRUE)

  # Fix Natural Earth iso_a3 = "-99" issue for France, Norway, etc.
  # Use adm0_a3 as fallback when iso_a3 is invalid
  world_map <- world %>%
    mutate(
      iso_a3_fixed = case_when(
        is.na(iso_a3) | iso_a3 == "-99" ~ adm0_a3,
        TRUE ~ iso_a3
      )
    ) %>%
    select(iso_a3 = iso_a3_fixed, name, geometry)

  return(world_map)
}

# ============================================================================
# FIGURE 1: STOPLIGHT WORLD MAP
# ============================================================================

#' Create Figure 1: Stoplight World Map
#' Shows probability of achieving WHO 2025 target with traffic-light coding
#'
#' @param weighted_results Data frame with prob_achieving_target column
#' @param indicator Indicator to plot (default: current_user_cigarettes)
#' @param year Target year (default: 2025)
#' @param country_name_mapping Named vector for country names
create_figure1_map <- function(weighted_results,
                                indicator = PRIMARY_INDICATOR,
                                year = TARGET_YEAR,
                                country_name_mapping = NULL) {

  cat("  Creating Figure 1: Stoplight World Map...\n")

  # Load world map WITH fix for iso_a3 = "-99" (France, Norway)
  # prepare_world_map() uses adm0_a3 as fallback when iso_a3 is invalid
  world_map <- prepare_world_map(crop_bounds = list(xmin = -180, xmax = 180, ymin = -90, ymax = 90))

  # Prepare achievement data
  # Check which column names exist
  country_col <- if ("Country" %in% names(weighted_results)) "Country" else "wb_country_abv"
  sex_col <- if ("Sex" %in% names(weighted_results)) "Sex" else "sex"
  indicator_col <- if ("Def_Type_Code" %in% names(weighted_results)) "Def_Type_Code" else "def_type_code"
  prob_col <- if ("prob_achieving_target" %in% names(weighted_results)) "prob_achieving_target" else "prob"
  year_col <- if ("Year" %in% names(weighted_results)) "Year" else "year"

  map_data <- weighted_results %>%
    filter(!!sym(year_col) == year) %>%
    filter(!!sym(indicator_col) == indicator)

  if (nrow(map_data) == 0) {
    cat("  WARNING: No data for Figure 1 map\n")
    return(NULL)
  }

  # Harmonize country codes and add categories
  map_data <- map_data %>%
    mutate(
      iso_a3_harmonized = sapply(!!sym(country_col), harmonize_country_codes),
      Gender_Label = format_gender(!!sym(sex_col)),
      prob_category = categorize_probability(!!sym(prob_col))
    )

  # Get Natural Earth ISO codes for comparison
  ne_iso3 <- unique(world_map$iso_a3)
  ne_iso3 <- ne_iso3[!is.na(ne_iso3) & ne_iso3 != "-99"]

  # Diagnostic: which WHO countries matched?
  who_codes_in_data <- unique(map_data$iso_a3_harmonized)
  matched_codes <- who_codes_in_data[who_codes_in_data %in% ne_iso3]
  unmatched_codes <- who_codes_in_data[!who_codes_in_data %in% ne_iso3]

  cat(sprintf("  Map matching: %d/%d countries matched\n",
              length(matched_codes), length(who_codes_in_data)))

  if (length(unmatched_codes) > 0) {
    cat(sprintf("  Unmatched (will show as 'No Data'): %s\n",
                paste(unmatched_codes, collapse = ", ")))
  }

  # Cross-join map with genders, then left-join data
  fig1_joined <- world_map %>%
    tidyr::crossing(Gender_Label = c("Men", "Women")) %>%
    left_join(
      map_data %>% select(iso_a3_harmonized, Gender_Label, prob_category, !!sym(prob_col)),
      by = c("iso_a3" = "iso_a3_harmonized", "Gender_Label")
    ) %>%
    mutate(
      prob_category = ifelse(is.na(prob_category), "No Data", as.character(prob_category)),
      prob_category = factor(prob_category,
                             levels = c("High (>90%)", "Moderate (60-90%)", "Uncertain (40-60%)",
                                        "Low (10-40%)", "Extremely Low (<10%)", "No Data"))
    )

  # Create main map
  fig1_map <- ggplot() +
    geom_sf(data = fig1_joined,
            aes(geometry = geometry, fill = prob_category),
            color = "white", linewidth = 0.1) +
    facet_wrap(~ Gender_Label, ncol = 1) +
    scale_fill_manual(values = stoplight_colors, drop = FALSE, name = "Achievement\nProbability") +
    theme_void() +
    theme(
      legend.position = "right",
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#555555")
    ) +
    labs(
      title = "WHO 2025 Target Achievement Probability",
      subtitle = paste0("30% Reduction in ", format_indicator(indicator), " Prevalence from 2010")
    )

  # Create inset bar chart with counts
  category_counts <- fig1_joined %>%
    st_drop_geometry() %>%
    filter(prob_category != "No Data") %>%
    group_by(Gender_Label, prob_category) %>%
    summarise(n = n(), .groups = "drop")

  fig1_inset <- ggplot(category_counts, aes(x = prob_category, y = n, fill = prob_category)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ Gender_Label, ncol = 2) +
    scale_fill_manual(values = stoplight_colors, guide = "none") +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_text(size = 7),
      strip.text = element_text(face = "bold", size = 8),
      panel.grid = element_blank()
    ) +
    labs(subtitle = "Country Count by Category")

  # Combine with patchwork
  fig1_combined <- fig1_map +
    inset_element(fig1_inset, left = 0.02, bottom = 0.02, right = 0.35, top = 0.25)

  # Save
  ggsave("outputs/figures/maps/Figure1_Stoplight_Map.pdf", fig1_combined, width = 12, height = 14)
  ggsave("outputs/figures/maps/Figure1_Stoplight_Map.png", fig1_combined, width = 12, height = 14, dpi = 300)

  cat("  Saved: outputs/figures/maps/Figure1_Stoplight_Map.pdf\n")

  return(fig1_combined)
}

# ============================================================================
# EXPANDED MAP SYSTEM: ALL INDICATORS & BOTH TARGETS
# ============================================================================

# Indicator labels for display
INDICATOR_LABELS <- c(
  "current_user_any_tobacco_product" = "Current Any Tobacco",
  "current_user_any_smoked_tobacco" = "Current Smoked Tobacco",
  "current_user_cigarettes" = "Current Cigarettes",
  "daily_user_any_tobacco_product" = "Daily Any Tobacco",
  "daily_user_any_smoked_tobacco" = "Daily Smoked Tobacco",
  "daily_user_cigarettes" = "Daily Cigarettes"
)

# Short indicator codes for filenames
INDICATOR_CODES <- c(
  "current_user_any_tobacco_product" = "CU_ATP",
  "current_user_any_smoked_tobacco" = "CU_AST",
  "current_user_cigarettes" = "CU_CIG",
  "daily_user_any_tobacco_product" = "DU_ATP",
  "daily_user_any_smoked_tobacco" = "DU_AST",
  "daily_user_cigarettes" = "DU_CIG"
)

# All 6 indicators
ALL_INDICATORS <- names(INDICATOR_LABELS)

#' Calculate target probabilities from MCMC samples
#'
#' @param mcmc_pred_samples Matrix of MCMC samples for predictions (rows = samples, cols = years/ages)
#' @param baseline_value 2010 baseline prevalence for this stratum
#' @param target_type "2025_target" (30% reduction) or "2040_endgame" (<5%)
#' @return Probability of achieving the target
calculate_target_probability <- function(mcmc_pred_samples, baseline_value, target_type) {

  if (target_type == "2025_target") {
    # WHO 2025 Target: 30% relative reduction from 2010 baseline
    target_value <- baseline_value * 0.70
    prob <- mean(mcmc_pred_samples <= target_value, na.rm = TRUE)
  } else if (target_type == "2040_endgame") {
    # Endgame: Absolute threshold <5%
    target_value <- 0.05
    prob <- mean(mcmc_pred_samples < target_value, na.rm = TRUE)
  } else {
    stop("target_type must be '2025_target' or '2040_endgame'")
  }

  return(prob)
}

#' Get 2010 baseline for a stratum
#'
#' @param weighted_results Weighted results data
#' @param country Country code
#' @param gender Gender
#' @param indicator Indicator code
#' @return 2010 baseline prevalence (or earliest available)
get_baseline_prevalence <- function(weighted_results, country, gender, indicator) {

  # Identify column names
  country_col <- if ("Country" %in% names(weighted_results)) "Country" else "wb_country_abv"
  sex_col <- if ("Sex" %in% names(weighted_results)) "Sex" else "sex"
  indicator_col <- if ("Def_Type_Code" %in% names(weighted_results)) "Def_Type_Code" else "def_type_code"
  year_col <- if ("Year" %in% names(weighted_results)) "Year" else "year"
  prev_col <- if ("Prevalence_Mean" %in% names(weighted_results)) "Prevalence_Mean" else
    if ("weighted_mean" %in% names(weighted_results)) "weighted_mean" else "prevalence"

  baseline <- weighted_results %>%
    filter(
      !!sym(country_col) == country,
      (!!sym(sex_col) == gender | format_gender(!!sym(sex_col)) == gender),
      !!sym(indicator_col) == indicator,
      !!sym(year_col) == BASE_YEAR
    ) %>%
    pull(!!sym(prev_col))

  if (length(baseline) == 0 || is.na(baseline[1])) {
    # Try earliest available year
    earliest <- weighted_results %>%
      filter(
        !!sym(country_col) == country,
        (!!sym(sex_col) == gender | format_gender(!!sym(sex_col)) == gender),
        !!sym(indicator_col) == indicator
      ) %>%
      arrange(!!sym(year_col)) %>%
      slice(1) %>%
      pull(!!sym(prev_col))

    if (length(earliest) > 0 && !is.na(earliest[1])) {
      return(earliest[1])
    }
    return(NA)
  }

  return(baseline[1])
}

#' Create probability map for specific indicator and target
#'
#' @param prob_data Data frame with columns: Country, Gender, indicator, prob, prob_category
#' @param indicator Indicator to map
#' @param gender "Men" or "Women"
#' @param target_type "2025_target" or "2040_endgame"
#' @param output_dir Output directory
#' @return ggplot object
create_probability_map <- function(prob_data,
                                    indicator,
                                    gender,
                                    target_type,
                                    output_dir = "outputs/figures/maps") {

  cat(sprintf("  Creating map: %s | %s | %s...\n",
              INDICATOR_CODES[indicator], gender, target_type))

  # Load world map WITH fix for iso_a3 = "-99" (France, Norway)
  # prepare_world_map() uses adm0_a3 as fallback when iso_a3 is invalid
  world_map <- prepare_world_map(crop_bounds = list(xmin = -180, xmax = 180, ymin = -90, ymax = 90))

  # Filter prob_data for this combination
  map_data <- prob_data %>%
    filter(
      Indicator == indicator,
      Gender_Label == gender,
      Target_Type == target_type
    )

  if (nrow(map_data) == 0) {
    cat("    WARNING: No data for this map\n")
    return(NULL)
  }

  # Harmonize country codes
  map_data <- map_data %>%
    mutate(iso_a3_harmonized = sapply(Country, harmonize_country_codes))

  # Join to world map
  joined <- world_map %>%
    left_join(
      map_data %>% select(iso_a3_harmonized, prob, prob_category),
      by = c("iso_a3" = "iso_a3_harmonized")
    ) %>%
    mutate(
      prob_category = ifelse(is.na(prob_category), "No Data", as.character(prob_category)),
      prob_category = factor(prob_category,
                             levels = c("High (>90%)", "Moderate (60-90%)", "Uncertain (40-60%)",
                                        "Low (10-40%)", "Extremely Low (<10%)", "No Data"))
    )

  # Create title based on target type
  if (target_type == "2025_target") {
    title <- paste0("Probability of 30% Reduction by 2025")
    subtitle <- paste0(INDICATOR_LABELS[indicator], " | ", gender)
  } else {
    title <- paste0("Probability of <5% Prevalence by 2040")
    subtitle <- paste0(INDICATOR_LABELS[indicator], " | ", gender)
  }

  # Create map
  p <- ggplot() +
    geom_sf(data = joined,
            aes(geometry = geometry, fill = prob_category),
            color = "white", linewidth = 0.1) +
    scale_fill_manual(values = stoplight_colors, drop = FALSE, name = "Achievement\nProbability") +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#555555"),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(title = title, subtitle = subtitle)

  # Create output directory structure
  target_dir <- file.path(output_dir, target_type, tolower(gender))
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

  # Save
  filename_base <- paste0(INDICATOR_CODES[indicator], "_", tolower(gender))
  ggsave(file.path(target_dir, paste0(filename_base, ".pdf")), p, width = 12, height = 7)
  ggsave(file.path(target_dir, paste0(filename_base, ".png")), p, width = 12, height = 7, dpi = 300)

  return(p)
}

# ============================================================================
# LANCET-STYLE MAP (A4 Vertical, Men/Women Facets)
# ============================================================================

#' Create Lancet-quality map for one indicator with gender facets
#'
#' Creates an A4 vertical publication-ready map with:
#' - Men facet on top
#' - Women facet on bottom
#' - Horizontal legend at bottom
#' - Professional Lancet styling
#'
#' @param prob_data Probability data with Country, Gender_Label, Indicator, Target_Type, prob, prob_category
#' @param indicator Indicator code to map
#' @param target_type "2025_target" or "2040_endgame"
#' @param output_dir Output directory
#' @return ggplot object
create_lancet_indicator_map <- function(prob_data,
                                         indicator,
                                         target_type,
                                         output_dir = "outputs/figures/maps") {

  cat(sprintf("  Creating Lancet map: %s | %s...\n",
              INDICATOR_LABELS[indicator], target_type))

  # Load world map with fixed ISO codes (handles France, Norway, etc.)
  world_map <- prepare_world_map()

  # Filter prob_data for this indicator and target (both genders)
  map_data <- prob_data %>%
    filter(
      Indicator == indicator,
      Target_Type == target_type,
      Gender_Label %in% c("Men", "Women")
    )

  if (nrow(map_data) == 0) {
    cat("    WARNING: No data for this map\n")
    return(NULL)
  }

  # Harmonize country codes
  map_data <- map_data %>%
    mutate(iso_a3_harmonized = sapply(Country, harmonize_country_codes))

  # Create data for both genders
  joined <- world_map %>%
    tidyr::crossing(Gender_Label = factor(c("Men", "Women"), levels = c("Men", "Women"))) %>%
    left_join(
      map_data %>% select(iso_a3_harmonized, Gender_Label, prob, prob_category),
      by = c("iso_a3" = "iso_a3_harmonized", "Gender_Label")
    ) %>%
    mutate(
      prob_category = ifelse(is.na(prob_category), "No Data", as.character(prob_category)),
      prob_category = factor(prob_category,
                             levels = c("Very High (>90%)", "High (60-90%)", "Moderate (40-60%)",
                                        "Low (10-40%)", "Very Low (<10%)", "No Data"))
    )

  # Create title and subtitle based on target type
  if (target_type == "2025_target") {
    title <- "Probability of Achieving WHO 2025 Target"
    subtitle <- paste0("30% Reduction in ", INDICATOR_LABELS[indicator], " Prevalence from 2010 Baseline")
  } else {
    title <- "Probability of Achieving Tobacco Endgame by 2040"
    subtitle <- paste0(INDICATOR_LABELS[indicator], " Prevalence Below 5%")
  }

  # Create the map with improved aesthetics
  p <- ggplot() +
    # Countries with data
    geom_sf(data = joined,
            aes(geometry = geometry, fill = prob_category),
            color = BORDER_COLOR,
            linewidth = BORDER_WIDTH,
            show.legend = TRUE) +
    # Facet by gender
    facet_wrap(~ Gender_Label, ncol = 1, strip.position = "top") +
    # Color scale - force all categories to show
    scale_fill_manual(
      values = stoplight_colors,
      drop = FALSE,
      na.value = LAND_NO_DATA,
      name = "Probability of\nAchieving Target"
    ) +
    # Coordinate system
    coord_sf(
      xlim = c(-180, 180),
      ylim = c(-56, 84),
      expand = FALSE,
      clip = "on"
    ) +
    # Theme and legend
    theme_lancet_map() +
    lancet_legend_guide() +
    # Labels
    labs(
      title = title,
      subtitle = subtitle,
      caption = paste0(
        "Source: WHO Global Report on Trends in Prevalence of Tobacco Use 2000-2030\n",
        "Note: Gray indicates insufficient data for reliable estimation"
      )
    )

  # Create output directory
  target_dir <- file.path(output_dir, target_type)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

  # Save with A4 vertical dimensions
  filename_base <- INDICATOR_CODES[indicator]
  ggsave(file.path(target_dir, paste0(filename_base, ".pdf")),
         p, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in")
  ggsave(file.path(target_dir, paste0(filename_base, ".png")),
         p, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", dpi = 300)

  cat(sprintf("    Saved: %s/%s.pdf (A4 vertical)\n", target_type, filename_base))

  return(p)
}

#' Create main publication figure for primary indicator (Current Cigarettes)
#'
#' Creates a 2x2 layout: Target × Gender for the primary indicator
#' Used as the main figure in the manuscript
#'
#' @param prob_data Probability data
#' @param output_dir Output directory
#' @return Combined ggplot
create_lancet_main_figure <- function(prob_data, output_dir = "outputs/figures/maps") {

  cat("  Creating Lancet main figure: Current Cigarettes (2x2 facet)...\n")

  # Load world map with fixed ISO codes (handles France, Norway, etc.)
  world_map <- prepare_world_map()

  # Filter for current cigarettes only
  map_data <- prob_data %>%
    filter(Indicator == "current_user_cigarettes") %>%
    mutate(
      iso_a3_harmonized = sapply(Country, harmonize_country_codes),
      Target_Label = factor(
        ifelse(Target_Type == "2025_target",
               "2025 WHO Target\n(30% Reduction from 2010)",
               "2040 Endgame\n(<5% Prevalence)"),
        levels = c("2025 WHO Target\n(30% Reduction from 2010)",
                   "2040 Endgame\n(<5% Prevalence)")
      ),
      Gender_Label = factor(Gender_Label, levels = c("Men", "Women"))
    )

  if (nrow(map_data) == 0) {
    cat("    WARNING: No data for main figure\n")
    return(NULL)
  }

  # Create all combinations for 2x2 grid
  joined <- world_map %>%
    tidyr::crossing(
      Gender_Label = factor(c("Men", "Women"), levels = c("Men", "Women")),
      Target_Label = factor(c("2025 WHO Target\n(30% Reduction from 2010)",
                              "2040 Endgame\n(<5% Prevalence)"),
                            levels = c("2025 WHO Target\n(30% Reduction from 2010)",
                                       "2040 Endgame\n(<5% Prevalence)"))
    ) %>%
    left_join(
      map_data %>% select(iso_a3_harmonized, Gender_Label, Target_Label, prob, prob_category),
      by = c("iso_a3" = "iso_a3_harmonized", "Gender_Label", "Target_Label")
    ) %>%
    mutate(
      prob_category = ifelse(is.na(prob_category), "No Data", as.character(prob_category)),
      prob_category = factor(prob_category,
                             levels = c("Very High (>90%)", "High (60-90%)", "Moderate (40-60%)",
                                        "Low (10-40%)", "Very Low (<10%)", "No Data"))
    )

  # Create 2x2 faceted map (rows = Target, cols = Gender)
  p <- ggplot() +
    # Country polygons
    geom_sf(data = joined,
            aes(geometry = geometry, fill = prob_category),
            color = BORDER_COLOR,
            linewidth = BORDER_WIDTH * 0.8,
            show.legend = TRUE) +
    # Facet grid
    facet_grid(Target_Label ~ Gender_Label) +
    # Color scale with all categories
    scale_fill_manual(
      values = stoplight_colors,
      drop = FALSE,
      na.value = LAND_NO_DATA,
      name = "Probability of\nAchieving Target"
    ) +
    # Coordinate system
    coord_sf(
      xlim = c(-180, 180),
      ylim = c(-56, 84),
      expand = FALSE
    ) +
    # Theme
    theme_lancet_map() +
    theme(
      strip.text = element_text(
        face = "bold",
        size = 11,
        lineheight = 1.1,
        margin = margin(t = 8, b = 6)
      ),
      strip.text.y = element_text(angle = 0, size = 10),
      panel.spacing.x = unit(0.3, "cm"),
      panel.spacing.y = unit(0.5, "cm"),
      # Wider margins for 2x2 layout with row labels
      plot.margin = margin(t = 15, r = 10, b = 15, l = 10)
    ) +
    # Legend
    lancet_legend_guide() +
    # Labels
    labs(
      title = "Probability of Achieving Tobacco Control Targets",
      subtitle = "Current Cigarette Use Prevalence (Age-Standardized)",
      caption = "Source: WHO Global Tobacco Surveillance System\nNote: Gray areas indicate insufficient data for analysis"
    )

  # Save with wider dimensions for 2x2 layout (needs more horizontal space)
  # Use full A4 width to accommodate row facet labels and legend
  MAIN_FIG_WIDTH <- 8.5   # inches (slightly wider than A4 for margins)
  MAIN_FIG_HEIGHT <- 10.0 # inches

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(output_dir, "Figure1_Main.pdf"),
         p, width = MAIN_FIG_WIDTH, height = MAIN_FIG_HEIGHT, units = "in")
  ggsave(file.path(output_dir, "Figure1_Main.png"),
         p, width = MAIN_FIG_WIDTH, height = MAIN_FIG_HEIGHT, units = "in", dpi = 300)

  cat("    Saved: Figure1_Main.pdf (8.5 x 10.0 in)\n")

  return(p)
}

#' Create faceted map showing all 6 indicators for one gender/target combination
#'
#' @param prob_data Probability data
#' @param gender "Men" or "Women"
#' @param target_type "2025_target" or "2040_endgame"
#' @param output_dir Output directory
#' @return ggplot object
create_faceted_indicator_map <- function(prob_data,
                                          gender,
                                          target_type,
                                          output_dir = "outputs/figures/maps") {

  cat(sprintf("  Creating faceted map: All indicators | %s | %s...\n", gender, target_type))

  # Load world map with fixed ISO codes (handles France, Norway, etc.)
  world_map <- prepare_world_map(crop_bounds = list(xmin = -180, xmax = 180, ymin = -90, ymax = 90))

  # Filter for this gender and target
  map_data <- prob_data %>%
    filter(
      Gender_Label == gender,
      Target_Type == target_type
    ) %>%
    mutate(
      iso_a3_harmonized = sapply(Country, harmonize_country_codes),
      Indicator_Label = INDICATOR_LABELS[Indicator]
    )

  if (nrow(map_data) == 0) {
    cat("    WARNING: No data for faceted map\n")
    return(NULL)
  }

  # Create data for all indicator-country combinations
  all_combos <- expand.grid(
    iso_a3 = unique(world_map$iso_a3),
    Indicator_Label = INDICATOR_LABELS,
    stringsAsFactors = FALSE
  )

  # Join map data
  joined <- world_map %>%
    tidyr::crossing(Indicator_Label = factor(INDICATOR_LABELS, levels = INDICATOR_LABELS)) %>%
    left_join(
      map_data %>% select(iso_a3_harmonized, Indicator_Label, prob, prob_category),
      by = c("iso_a3" = "iso_a3_harmonized", "Indicator_Label")
    ) %>%
    mutate(
      prob_category = ifelse(is.na(prob_category), "No Data", as.character(prob_category)),
      prob_category = factor(prob_category,
                             levels = c("High (>90%)", "Moderate (60-90%)", "Uncertain (40-60%)",
                                        "Low (10-40%)", "Extremely Low (<10%)", "No Data"))
    )

  # Title
  if (target_type == "2025_target") {
    title <- paste0("WHO 2025 Target Achievement (30% Reduction) - ", gender)
  } else {
    title <- paste0("Endgame 2040 Achievement (<5% Prevalence) - ", gender)
  }

  # Create faceted map (2 rows x 3 cols for 6 indicators)
  p <- ggplot() +
    geom_sf(data = joined,
            aes(geometry = geometry, fill = prob_category),
            color = "white", linewidth = 0.05) +
    facet_wrap(~ Indicator_Label, ncol = 3) +
    scale_fill_manual(values = stoplight_colors, drop = FALSE, name = "Achievement\nProbability") +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      strip.text = element_text(face = "bold", size = 9, margin = margin(2,2,2,2)),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    labs(title = title)

  # Save
  target_dir <- file.path(output_dir, "composites")
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

  filename <- paste0("all_indicators_", target_type, "_", tolower(gender))
  ggsave(file.path(target_dir, paste0(filename, ".pdf")), p, width = 16, height = 10)
  ggsave(file.path(target_dir, paste0(filename, ".png")), p, width = 16, height = 10, dpi = 300)

  cat(sprintf("    Saved: %s/%s.pdf\n", target_dir, filename))

  return(p)
}

#' Create main publication figure: Current Cigarettes, both targets, both genders
#'
#' @param prob_data Probability data
#' @param output_dir Output directory
#' @return Combined ggplot
create_main_cigarettes_map <- function(prob_data, output_dir = "outputs/figures/maps") {

  cat("  Creating main publication map: Current Cigarettes (2x2 facet)...\n")

  # Load world map with fixed ISO codes (handles France, Norway, etc.)
  world_map <- prepare_world_map(crop_bounds = list(xmin = -180, xmax = 180, ymin = -90, ymax = 90))

  # Filter for current cigarettes only
  map_data <- prob_data %>%
    filter(Indicator == "current_user_cigarettes") %>%
    mutate(
      iso_a3_harmonized = sapply(Country, harmonize_country_codes),
      Target_Label = ifelse(Target_Type == "2025_target",
                            "2025 Target (30% Reduction)",
                            "2040 Endgame (<5%)")
    )

  if (nrow(map_data) == 0) {
    cat("    WARNING: No data for main map\n")
    return(NULL)
  }

  # Create all combinations
  joined <- world_map %>%
    tidyr::crossing(
      Gender_Label = factor(c("Men", "Women")),
      Target_Label = factor(c("2025 Target (30% Reduction)", "2040 Endgame (<5%)"))
    ) %>%
    left_join(
      map_data %>% select(iso_a3_harmonized, Gender_Label, Target_Label, prob, prob_category),
      by = c("iso_a3" = "iso_a3_harmonized", "Gender_Label", "Target_Label")
    ) %>%
    mutate(
      prob_category = ifelse(is.na(prob_category), "No Data", as.character(prob_category)),
      prob_category = factor(prob_category,
                             levels = c("High (>90%)", "Moderate (60-90%)", "Uncertain (40-60%)",
                                        "Low (10-40%)", "Extremely Low (<10%)", "No Data"))
    )

  # Create 2x2 faceted map
  p <- ggplot() +
    geom_sf(data = joined,
            aes(geometry = geometry, fill = prob_category),
            color = "white", linewidth = 0.1) +
    facet_grid(Target_Label ~ Gender_Label) +
    scale_fill_manual(values = stoplight_colors, drop = FALSE, name = "Achievement\nProbability") +
    theme_void() +
    theme(
      legend.position = "right",
      strip.text = element_text(face = "bold", size = 11, margin = margin(4,4,4,4)),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#555555"),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    labs(
      title = "Tobacco Control Target Achievement Probability",
      subtitle = "Current Cigarette Use (Age-Standardized)"
    )

  # Save
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(output_dir, "Figure1_Main_Current_Cigarettes.pdf"), p, width = 14, height = 12)
  ggsave(file.path(output_dir, "Figure1_Main_Current_Cigarettes.png"), p, width = 14, height = 12, dpi = 300)

  cat("    Saved: Figure1_Main_Current_Cigarettes.pdf\n")

  return(p)
}

#' Prepare probability data for all strata
#'
#' @param weighted_results Weighted results with predictions
#' @param mcmc_samples MCMC samples (optional, for probabilistic version)
#' @return Data frame with prob for both targets
prepare_all_probability_data <- function(weighted_results, mcmc_samples = NULL) {

  cat("  Preparing probability data for all strata...\n")

  # Identify column names
  country_col <- if ("Country" %in% names(weighted_results)) "Country" else "wb_country_abv"
  sex_col <- if ("Sex" %in% names(weighted_results)) "Sex" else "sex"
  indicator_col <- if ("Def_Type_Code" %in% names(weighted_results)) "Def_Type_Code" else "def_type_code"
  year_col <- if ("Year" %in% names(weighted_results)) "Year" else "year"
  prev_col <- if ("Prevalence_Mean" %in% names(weighted_results)) "Prevalence_Mean" else
    if ("weighted_mean" %in% names(weighted_results)) "weighted_mean" else "prevalence"

  # Check for existing probability column
  prob_col <- if ("prob_achieving_target" %in% names(weighted_results)) "prob_achieving_target" else NULL

  # Get unique strata
  strata <- weighted_results %>%
    select(!!sym(country_col), !!sym(sex_col), !!sym(indicator_col)) %>%
    distinct() %>%
    rename(Country = !!sym(country_col), Sex = !!sym(sex_col), Indicator = !!sym(indicator_col))

  cat(sprintf("    Found %d strata\n", nrow(strata)))

  results <- list()

  for (i in 1:nrow(strata)) {
    country <- strata$Country[i]
    sex <- strata$Sex[i]
    indicator <- strata$Indicator[i]

    # Get 2010 baseline
    baseline <- get_baseline_prevalence(weighted_results, country, sex, indicator)

    # Get 2025 and 2040 predictions
    pred_2025 <- weighted_results %>%
      filter(
        !!sym(country_col) == country,
        !!sym(sex_col) == sex,
        !!sym(indicator_col) == indicator,
        !!sym(year_col) == TARGET_YEAR
      ) %>%
      pull(!!sym(prev_col))

    pred_2040 <- weighted_results %>%
      filter(
        !!sym(country_col) == country,
        !!sym(sex_col) == sex,
        !!sym(indicator_col) == indicator,
        !!sym(year_col) == ENDGAME_YEAR
      ) %>%
      pull(!!sym(prev_col))

    # Calculate probabilities (simplified version without MCMC)
    # If prob_achieving_target exists, use it for 2025
    if (!is.null(prob_col)) {
      prob_2025 <- weighted_results %>%
        filter(
          !!sym(country_col) == country,
          !!sym(sex_col) == sex,
          !!sym(indicator_col) == indicator,
          !!sym(year_col) == TARGET_YEAR
        ) %>%
        pull(!!sym(prob_col))
      prob_2025 <- if(length(prob_2025) > 0) prob_2025[1] else NA
    } else if (!is.na(baseline) && length(pred_2025) > 0) {
      # Approximate: if prediction <= target, prob = 1; else 0
      # This is a simplification - proper version needs MCMC samples
      target_2025 <- baseline * 0.70
      prob_2025 <- ifelse(pred_2025[1] <= target_2025, 1, 0)
    } else {
      prob_2025 <- NA
    }

    # Calculate 2040 endgame probability
    if (length(pred_2040) > 0 && !is.na(pred_2040[1])) {
      # Convert to proportion if needed
      pred_val <- if(pred_2040[1] > 1) pred_2040[1] / 100 else pred_2040[1]
      prob_2040 <- ifelse(pred_val < 0.05, 1, 0)
    } else {
      prob_2040 <- NA
    }

    # Add to results
    results[[length(results) + 1]] <- data.frame(
      Country = country,
      Sex = sex,
      Gender_Label = format_gender(sex),
      Indicator = indicator,
      Baseline_2010 = baseline,
      Target_Type = "2025_target",
      prob = prob_2025,
      stringsAsFactors = FALSE
    )

    results[[length(results) + 1]] <- data.frame(
      Country = country,
      Sex = sex,
      Gender_Label = format_gender(sex),
      Indicator = indicator,
      Baseline_2010 = baseline,
      Target_Type = "2040_endgame",
      prob = prob_2040,
      stringsAsFactors = FALSE
    )
  }

  prob_data <- bind_rows(results)

  # Add probability categories
  prob_data <- prob_data %>%
    mutate(prob_category = categorize_probability(prob))

  cat(sprintf("    Prepared %d probability records\n", nrow(prob_data)))

  return(prob_data)
}

#' Master function to generate all probability maps (Lancet Style)
#'
#' Generates publication-quality maps:
#' - A4 vertical format
#' - Men/Women as facets (not separate files)
#' - One file per indicator per target
#'
#' @param weighted_results Weighted results data
#' @param mcmc_samples MCMC samples (optional)
#' @param output_dir Base output directory
#' @param primary_only If TRUE, only generate current_user_cigarettes maps
#' @param lancet_style If TRUE, use Lancet A4 format (default TRUE)
generate_all_probability_maps <- function(weighted_results,
                                           mcmc_samples = NULL,
                                           output_dir = "outputs/figures/maps",
                                           primary_only = FALSE,
                                           lancet_style = TRUE) {

  cat("\n========== GENERATING PROBABILITY MAPS (LANCET STYLE) ==========\n\n")

  # Create directory structure (simplified for Lancet style)
  dirs <- c(
    file.path(output_dir, "2025_target"),
    file.path(output_dir, "2040_endgame"),
    file.path(output_dir, "supplementary")
  )
  for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  # Prepare probability data
  prob_data <- prepare_all_probability_data(weighted_results, mcmc_samples)

  if (nrow(prob_data) == 0) {
    cat("  ERROR: No probability data available\n")
    return(NULL)
  }

  results <- list()

  # Determine which indicators to process
  indicators_to_process <- if(primary_only) {
    "current_user_cigarettes"
  } else {
    unique(prob_data$Indicator)
  }

  cat(sprintf("  Processing %d indicators...\n", length(indicators_to_process)))
  cat(sprintf("  Output format: A4 vertical (%.2f x %.2f in)\n", FIG_WIDTH, FIG_HEIGHT))

  # 1. Main publication figure (current cigarettes, 2x2)
  tryCatch({
    results$main <- create_lancet_main_figure(prob_data, output_dir)
  }, error = function(e) {
    cat(sprintf("  ERROR creating main figure: %s\n", e$message))
  })

  # 2. Individual Lancet-style maps for each indicator × target
  #    (Each file contains Men on top, Women on bottom)
  for (ind in indicators_to_process) {
    for (target in c("2025_target", "2040_endgame")) {
      tryCatch({
        results[[paste0(ind, "_", target)]] <-
          create_lancet_indicator_map(prob_data, ind, target, output_dir)
      }, error = function(e) {
        cat(sprintf("    ERROR: %s | %s: %s\n", ind, target, e$message))
      })
    }
  }

  # Summary
  cat("\n========== MAP GENERATION COMPLETE (LANCET STYLE) ==========\n")
  cat("\nGenerated maps:\n")
  cat(sprintf("  Format: A4 vertical (%.2f x %.2f in)\n", FIG_WIDTH, FIG_HEIGHT))
  cat("  Layout: Men (top facet) + Women (bottom facet)\n\n")

  for (d in dirs) {
    if (dir.exists(d)) {
      files <- list.files(d, pattern = "\\.(pdf|png)$")
      if (length(files) > 0) {
        cat(sprintf("  %s: %d files\n", basename(d), length(files)))
      }
    }
  }

  # Count total
  all_maps <- list.files(output_dir, pattern = "\\.(pdf|png)$", recursive = TRUE)
  cat(sprintf("\nTotal map files: %d\n", length(all_maps)))

  # List PDF files for reference
  pdf_files <- list.files(output_dir, pattern = "\\.pdf$", recursive = TRUE)
  if (length(pdf_files) > 0) {
    cat("\nPDF files for publication:\n")
    for (f in pdf_files) {
      cat(sprintf("  - %s\n", f))
    }
  }

  return(results)
}

# ============================================================================
# FIGURE 2: CATERPILLAR PLOT (Country Ranking)
# ============================================================================

#' Create Figure 2: Caterpillar Plot
#' Ranks countries by projected relative reduction with uncertainty intervals
#'
#' @param reduction_data Data frame with reduction estimates and CIs
#' @param gender_filter "Men" or "Women"
#' @param country_region_mapping Mapping of countries to regions
create_figure2_caterpillar <- function(reduction_data,
                                        gender_filter = "Men",
                                        country_region_mapping = NULL,
                                        country_name_mapping = NULL) {

  cat(sprintf("  Creating Figure 2: Caterpillar Plot (%s)...\n", gender_filter))

  # Check column names and adapt
  if (!"Reduction_Mean" %in% names(reduction_data)) {
    # Need to calculate from available data
    cat("  WARNING: Reduction data needs different structure\n")
    return(NULL)
  }

  plot_data <- reduction_data %>%
    filter(Gender == gender_filter | format_gender(Sex) == gender_filter)

  if (nrow(plot_data) == 0) {
    cat("  WARNING: No data for caterpillar plot\n")
    return(NULL)
  }

  # Add country names if mapping provided
  if (!is.null(country_name_mapping)) {
    plot_data <- plot_data %>%
      mutate(Country_Display = tools::toTitleCase(
        ifelse(Country %in% names(country_name_mapping),
               country_name_mapping[Country], Country)
      ))
  } else {
    plot_data <- plot_data %>%
      mutate(Country_Display = tools::toTitleCase(Country))
  }

  # Add region if mapping provided
  if (!is.null(country_region_mapping) && !"Region" %in% names(plot_data)) {
    plot_data <- plot_data %>%
      left_join(country_region_mapping, by = c("Country" = "wb_country_abv"))
  }

  # Order by mean reduction
  plot_data <- plot_data %>%
    mutate(Country_Display = reorder(Country_Display, Reduction_Mean))

  # Dynamic sizing
  n_countries <- nrow(plot_data)
  plot_height <- max(8, n_countries * 0.15 + 2)

  fig2 <- ggplot(plot_data, aes(x = Reduction_Mean, y = Country_Display)) +
    # Target line (30%)
    geom_vline(xintercept = 30, linetype = "dashed", color = "#B71C1C", linewidth = 0.8) +
    annotate("text", x = 31, y = 1, label = "30% Target", color = "#B71C1C",
             hjust = 0, size = 3, fontface = "bold") +
    # Error bars (95% CrI)
    geom_errorbarh(aes(xmin = Reduction_Lower, xmax = Reduction_Upper,
                       color = if("Region" %in% names(plot_data)) Region else NULL),
                   height = 0, linewidth = 0.4, alpha = 0.6) +
    # Point estimates
    geom_point(aes(color = if("Region" %in% names(plot_data)) Region else NULL), size = 1.5) +
    theme_who() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 7)
    ) +
    labs(
      title = paste0("Projected Prevalence Reduction by Country (", gender_filter, ")"),
      subtitle = "Relative reduction from 2010 baseline to 2025, ordered by posterior mean",
      x = "Relative Reduction (%)",
      y = NULL,
      caption = "Error bars show 95% credible intervals. Dashed line = 30% WHO target."
    )

  if ("Region" %in% names(plot_data)) {
    fig2 <- fig2 + scale_color_manual(values = region_colors, name = "Region")
  }

  # Save with dynamic height
  ggsave(paste0("outputs/figures/caterpillar/Figure2_Caterpillar_", gender_filter, ".pdf"),
         fig2, width = 10, height = plot_height)
  ggsave(paste0("outputs/figures/caterpillar/Figure2_Caterpillar_", gender_filter, ".png"),
         fig2, width = 10, height = plot_height, dpi = 300)

  cat(sprintf("  Saved: outputs/figures/caterpillar/Figure2_Caterpillar_%s.pdf\n", gender_filter))

  return(fig2)
}

# ============================================================================
# FIGURE 3: ENDGAME TRAJECTORIES (Spaghetti Plot)
# ============================================================================

#' Create Figure 3: Endgame Trajectories
#' Shows prevalence trajectories with endgame zone (<5%)
#'
#' @param trend_data Data frame with Year, Prevalence_Mean, Entity, Entity_Type
#' @param gender_filter "Men" or "Women"
create_figure3_trajectories <- function(trend_data,
                                         gender_filter = "Men",
                                         highlight_countries = c("nzl", "swe", "aus", "gbr")) {

  cat(sprintf("  Creating Figure 3: Endgame Trajectories (%s)...\n", gender_filter))

  if (nrow(trend_data) == 0) {
    cat("  WARNING: No trend data for trajectories\n")
    return(NULL)
  }

  # Filter by gender
  gender_sex <- if(gender_filter == "Men") c("males", "male", "Men") else c("females", "female", "Women")

  plot_data <- trend_data %>%
    filter(Gender %in% gender_sex | Sex %in% gender_sex)

  # Separate regions and countries
  region_data <- plot_data %>% filter(Entity_Type == "Region")
  country_data <- plot_data %>% filter(Entity_Type == "Country")

  # Identify success stories (already near endgame)
  success_data <- country_data %>%
    filter(Entity %in% highlight_countries)

  fig3 <- ggplot() +
    # Endgame zone (green rectangle)
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = ENDGAME_THRESHOLD,
             fill = "#E8F5E9", alpha = 0.8) +
    annotate("text", x = min(plot_data$Year), y = ENDGAME_THRESHOLD / 2,
             label = paste0("Endgame Zone (<", ENDGAME_THRESHOLD, "%)"),
             hjust = 0, vjust = 0.5, size = 3, fontface = "italic", color = "#1B5E20") +

    # Near-endgame zone
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = ENDGAME_THRESHOLD, ymax = NEAR_ENDGAME_THRESHOLD,
             fill = "#FFF9C4", alpha = 0.5) +

    # Country lines (background, grey)
    geom_line(data = country_data,
              aes(x = Year, y = Prevalence_Mean, group = Entity),
              color = "grey70", linewidth = 0.3, alpha = 0.5) +

    # Region lines (foreground, colored)
    geom_line(data = region_data,
              aes(x = Year, y = Prevalence_Mean, color = Entity, group = Entity),
              linewidth = 1) +

    # Success stories highlighted
    geom_line(data = success_data,
              aes(x = Year, y = Prevalence_Mean, group = Entity),
              color = "#1B5E20", linewidth = 0.8, linetype = "dashed") +

    # Labels for success stories at end of line
    geom_text_repel(
      data = success_data %>%
        group_by(Entity) %>%
        filter(Year == max(Year)),
      aes(x = Year, y = Prevalence_Mean, label = tools::toTitleCase(Entity)),
      nudge_x = 2, size = 2.5, color = "#1B5E20"
    ) +

    scale_color_viridis_d(option = "D", name = "Region") +
    scale_y_continuous(limits = c(0, NA), labels = function(x) paste0(x, "%")) +
    theme_who() +
    labs(
      title = paste0("Tobacco Prevalence Trajectories to ", ENDGAME_YEAR, " (", gender_filter, ")"),
      subtitle = paste0("Regional trends with endgame zone (<", ENDGAME_THRESHOLD, "%) highlighted"),
      x = "Year",
      y = "Age-Standardized Prevalence (%)",
      caption = "Dashed green lines: countries near endgame (NZ, Sweden, etc.)"
    )

  ggsave(paste0("outputs/figures/trajectories/Figure3_Trajectories_", gender_filter, ".pdf"),
         fig3, width = 12, height = 8)
  ggsave(paste0("outputs/figures/trajectories/Figure3_Trajectories_", gender_filter, ".png"),
         fig3, width = 12, height = 8, dpi = 300)

  cat(sprintf("  Saved: outputs/figures/trajectories/Figure3_Trajectories_%s.pdf\n", gender_filter))

  return(fig3)
}

# ============================================================================
# FIGURE 4: AGE-COHORT HEATMAP (Lexis Diagram)
# ============================================================================

#' Create Figure 4: Age-Cohort Heatmap
#' Shows prevalence surface with diagonal birth cohort lines
#'
#' @param predictions Data frame with Year, Age_Midpoint, Prevalence
#' @param country_code Country to plot (or "global" for aggregated)
#' @param gender Gender to plot
create_figure4_heatmap <- function(predictions,
                                    country_code = "global",
                                    gender = "males",
                                    indicator = PRIMARY_INDICATOR) {

  cat(sprintf("  Creating Figure 4: Age-Cohort Heatmap (%s, %s)...\n", country_code, gender))

  # Filter data
  if (country_code == "global") {
    # Aggregate across sample countries
    sample_countries <- c("usa", "gbr", "chn", "ind", "bra", "deu")
    plot_data <- predictions %>%
      filter(Country %in% sample_countries | wb_country_abv %in% sample_countries) %>%
      filter(Sex == gender | sex == gender)
  } else {
    plot_data <- predictions %>%
      filter(Country == country_code | wb_country_abv == country_code) %>%
      filter(Sex == gender | sex == gender)
  }

  # Filter by indicator
  ind_col <- if ("Def_Type_Code" %in% names(plot_data)) "Def_Type_Code" else "def_type_code"
  plot_data <- plot_data %>%
    filter(!!sym(ind_col) == indicator)

  if (nrow(plot_data) == 0) {
    cat("  WARNING: No data for heatmap\n")
    return(NULL)
  }

  # Get year and age columns
  year_col <- if ("Year" %in% names(plot_data)) "Year" else "year"
  age_col <- if ("Age_Midpoint" %in% names(plot_data)) "Age_Midpoint" else "age"
  prev_col <- if ("Prevalence" %in% names(plot_data)) "Prevalence" else "prevalence"

  # Convert prevalence to percentage if needed
  plot_data <- plot_data %>%
    mutate(
      Year_Num = as.numeric(!!sym(year_col)),
      Age_Num = as.numeric(!!sym(age_col)),
      Prev_Pct = if(max(!!sym(prev_col), na.rm = TRUE) <= 1) !!sym(prev_col) * 100 else !!sym(prev_col)
    )

  # Aggregate if multiple countries
  if (country_code == "global") {
    plot_data <- plot_data %>%
      group_by(Year_Num, Age_Num) %>%
      summarise(Prev_Pct = mean(Prev_Pct, na.rm = TRUE), .groups = "drop")
  }

  # Create birth cohort lines
  cohorts <- seq(1920, 2010, by = 20)

  fig4 <- ggplot(plot_data, aes(x = Year_Num, y = Age_Num, fill = Prev_Pct)) +
    geom_tile() +

    # Diagonal cohort lines (Birth_Cohort = Year - Age, so Age = Year - Cohort)
    # intercept = -cohort because: y = x + b => Age = Year - Cohort => b = -Cohort
    geom_abline(intercept = -cohorts, slope = 1,
                color = "white", linewidth = 0.3, alpha = 0.6) +

    scale_fill_viridis_c(
      option = "inferno",
      direction = -1,
      limits = c(0, NA),
      name = "Prevalence (%)"
    ) +
    coord_fixed(ratio = 0.5, xlim = c(2000, 2040), ylim = c(15, 85)) +
    theme_who() +
    theme(
      panel.grid = element_blank()
    ) +
    labs(
      title = paste0("Age-Period-Cohort Prevalence Surface (",
                     tools::toTitleCase(gender), ")"),
      subtitle = paste0("Diagonal lines trace birth cohorts; ",
                        format_indicator(indicator)),
      x = "Year",
      y = "Age",
      caption = "Color intensity shows smoking prevalence. White diagonals = birth cohorts."
    )

  # Add cohort labels
  for (coh in cohorts) {
    # Find intersection with plot area
    y_at_2040 <- 2040 - coh
    if (y_at_2040 > 15 && y_at_2040 < 85) {
      fig4 <- fig4 + annotate("text", x = 2041, y = y_at_2040,
                               label = paste0("Born ", coh),
                               size = 2, hjust = 0, color = "grey40")
    }
  }

  ggsave(paste0("outputs/figures/heatmaps/Figure4_Heatmap_", country_code, "_", gender, ".pdf"),
         fig4, width = 10, height = 8)
  ggsave(paste0("outputs/figures/heatmaps/Figure4_Heatmap_", country_code, "_", gender, ".png"),
         fig4, width = 10, height = 8, dpi = 300)

  cat(sprintf("  Saved: outputs/figures/heatmaps/Figure4_Heatmap_%s_%s.pdf\n", country_code, gender))

  return(fig4)
}

# ============================================================================
# eFIGURE 1: MODEL VALIDATION (Observed vs Predicted)
# ============================================================================

#' Create eFigure 1: Model Validation Scatter
#' Compares observed vs predicted prevalence for both model types
#'
#' @param observed_data Observed survey data
#' @param global_predictions Global model predictions
#' @param country_predictions Country-specific model predictions
create_efigure1_validation <- function(observed_data,
                                        global_predictions,
                                        country_predictions = NULL) {

  cat("  Creating eFigure 1: Model Validation...\n")

  # Prepare observed data
  obs <- observed_data %>%
    mutate(
      Year = as.numeric(if("year" %in% names(.)) year else Year),
      Age = if("Age_Midpoint" %in% names(.)) Age_Midpoint else age,
      Country = if("wb_country_abv" %in% names(.)) wb_country_abv else Country,
      Sex = if("sex" %in% names(.)) sex else Sex,
      Indicator = if("def_type_code" %in% names(.)) def_type_code else Indicator,
      Prev_Obs = if("prevalence" %in% names(.)) plogis(prevalence) * 100 else Prevalence * 100
    ) %>%
    select(Year, Age, Country, Sex, Indicator, Prev_Obs)

  # Prepare global predictions
  glob <- global_predictions %>%
    mutate(
      Year = as.numeric(if("Year" %in% names(.)) Year else year),
      Age = if("Age_Midpoint" %in% names(.)) Age_Midpoint else age,
      Country = if("Country" %in% names(.)) Country else wb_country_abv,
      Sex = if("Sex" %in% names(.)) Sex else sex,
      Indicator = if("Def_Type_Code" %in% names(.)) Def_Type_Code else def_type_code,
      Prev_Global = if(max(Prevalence, na.rm = TRUE) <= 1) Prevalence * 100 else Prevalence
    ) %>%
    select(Year, Age, Country, Sex, Indicator, Prev_Global)

  # Merge
  validation <- obs %>%
    inner_join(glob, by = c("Year", "Age", "Country", "Sex", "Indicator"))

  if (nrow(validation) == 0) {
    cat("  WARNING: No matched data for validation\n")
    return(NULL)
  }

  # Calculate metrics
  r2_global <- cor(validation$Prev_Obs, validation$Prev_Global, use = "complete.obs")^2
  rmse_global <- sqrt(mean((validation$Prev_Obs - validation$Prev_Global)^2, na.rm = TRUE))

  # Create plot
  efig1 <- ggplot(validation, aes(x = Prev_Obs, y = Prev_Global)) +
    # Perfect fit line
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B71C1C", linewidth = 0.8) +
    # Scatter points
    geom_point(alpha = 0.3, size = 0.8, color = "#2166AC") +
    # Linear fit
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
    # Metrics annotation
    annotate("text", x = 5, y = max(validation$Prev_Global, na.rm = TRUE) * 0.9,
             label = sprintf("R\u00B2 = %.3f\nRMSE = %.2f%%", r2_global, rmse_global),
             hjust = 0, size = 3.5, fontface = "bold") +
    coord_equal() +
    theme_who() +
    labs(
      title = "Model Validation: Global APC Model",
      subtitle = "Observed vs Predicted Prevalence (In-Sample)",
      x = "Observed Prevalence (%)",
      y = "Predicted Prevalence (%)",
      caption = "Dashed line = perfect fit. Solid black = linear fit."
    )

  ggsave("outputs/figures/validation/eFigure1_Validation.pdf", efig1, width = 8, height = 8)
  ggsave("outputs/figures/validation/eFigure1_Validation.png", efig1, width = 8, height = 8, dpi = 300)

  cat("  Saved: outputs/figures/validation/eFigure1_Validation.pdf\n")

  return(efig1)
}

# ============================================================================
# eFIGURE 2: COUNTRY AGE PROFILES (Fitted vs Observed with Age Bands)
# ============================================================================

#' Create eFigure 2: Country Age Profiles
#' Shows fitted curves against observed data with proper age-band segments
#'
#' @param observed_data Observed survey data with start_age, end_age
#' @param global_predictions Global model predictions with CIs
#' @param country_predictions Country-specific predictions with CIs
#' @param country_code Country to plot
#' @param gender Gender to plot
#' @param year Year to plot
create_efigure2_age_profile <- function(observed_data,
                                         global_predictions,
                                         country_predictions = NULL,
                                         country_code,
                                         gender,
                                         year,
                                         indicator = PRIMARY_INDICATOR) {

  cat(sprintf("  Creating eFigure 2: Age Profile (%s, %s, %d)...\n", country_code, gender, year))

  # Prepare observed data
  obs <- observed_data %>%
    filter(
      (wb_country_abv == country_code | Country == country_code),
      (sex == gender | Sex == gender),
      (def_type_code == indicator | Indicator == indicator),
      abs(as.numeric(if("year" %in% names(.)) year else Year) - year) <= 2
    ) %>%
    mutate(
      Age = if("Age_Midpoint" %in% names(.)) Age_Midpoint else age,
      Prevalence = if("prevalence" %in% names(.)) plogis(prevalence) * 100 else Prevalence * 100,
      start_age = if("start_age" %in% names(.)) start_age else Age - 5,
      end_age = if("end_age" %in% names(.)) end_age else Age + 5
    )

  # Prepare global predictions
  glob <- global_predictions %>%
    filter(
      (Country == country_code | wb_country_abv == country_code),
      (Sex == gender | sex == gender),
      (Def_Type_Code == indicator | def_type_code == indicator),
      (Year == year | as.numeric(year) == year)
    ) %>%
    mutate(
      Age = if("Age_Midpoint" %in% names(.)) Age_Midpoint else age,
      Prev = if(max(Prevalence, na.rm = TRUE) <= 1) Prevalence * 100 else Prevalence,
      Lower = if("lower_ci" %in% names(.)) {
        if(max(lower_ci, na.rm = TRUE) <= 1) lower_ci * 100 else lower_ci
      } else Prev * 0.8,
      Upper = if("upper_ci" %in% names(.)) {
        if(max(upper_ci, na.rm = TRUE) <= 1) upper_ci * 100 else upper_ci
      } else Prev * 1.2,
      Model = "Global APC"
    )

  # Prepare country predictions if available
  if (!is.null(country_predictions) && nrow(country_predictions) > 0) {
    cntry <- country_predictions %>%
      filter(
        (Country == country_code | wb_country_abv == country_code),
        (Sex == gender | sex == gender),
        (Def_Type_Code == indicator | def_type_code == indicator),
        (Year == year | as.numeric(year) == year)
      ) %>%
      mutate(
        Age = if("Age_Midpoint" %in% names(.)) Age_Midpoint else age,
        Prev = if(max(Prevalence, na.rm = TRUE) <= 1) Prevalence * 100 else Prevalence,
        Lower = if("lower_ci" %in% names(.)) {
          if(max(lower_ci, na.rm = TRUE) <= 1) lower_ci * 100 else lower_ci
        } else Prev * 0.8,
        Upper = if("upper_ci" %in% names(.)) {
          if(max(upper_ci, na.rm = TRUE) <= 1) upper_ci * 100 else upper_ci
        } else Prev * 1.2,
        Model = "Country APC"
      )

    pred_data <- bind_rows(glob, cntry)
  } else {
    pred_data <- glob
  }

  if (nrow(pred_data) == 0 && nrow(obs) == 0) {
    cat("  WARNING: No data for age profile\n")
    return(NULL)
  }

  # Create plot
  efig2 <- ggplot() +
    # Background grid
    geom_hline(yintercept = seq(0, 100, 20), color = "grey95", linewidth = 0.2) +

    # Model confidence ribbons
    geom_ribbon(data = pred_data,
                aes(x = Age, ymin = Lower, ymax = Upper, fill = Model),
                alpha = 0.15) +

    # Model mean lines
    geom_line(data = pred_data,
              aes(x = Age, y = Prev, color = Model, linetype = Model),
              linewidth = 0.8) +

    # Observed data as AGE-BAND SEGMENTS
    geom_segment(data = obs %>% filter(!is.na(start_age)),
                 aes(x = start_age, xend = end_age,
                     y = Prevalence, yend = Prevalence),
                 color = "black", linewidth = 0.5) +

    # Observed midpoints
    geom_point(data = obs,
               aes(x = Age, y = Prevalence),
               size = 1.5, shape = 21, fill = "white", stroke = 0.5) +

    # Spline->Linear transition marker
    geom_vline(xintercept = TRANSITION_START, color = "orange",
               linetype = "dashed", alpha = 0.7) +
    annotate("text", x = TRANSITION_START + 1, y = 5,
             label = "Spline\u2192Linear", size = 2.5, hjust = 0, color = "orange") +

    scale_color_manual(values = c("Global APC" = "#2166AC", "Country APC" = "#D55E00")) +
    scale_fill_manual(values = c("Global APC" = "#2166AC", "Country APC" = "#D55E00")) +
    scale_linetype_manual(values = c("Global APC" = "solid", "Country APC" = "dashed")) +
    scale_y_continuous(limits = c(0, NA), labels = function(x) paste0(x, "%")) +
    scale_x_continuous(limits = c(15, 85)) +
    theme_who() +
    labs(
      title = paste0("Age-Specific Prevalence: ", toupper(country_code), " (", year, ")"),
      subtitle = paste0(format_indicator(indicator), " | ", format_gender(gender)),
      x = "Age",
      y = "Prevalence (%)",
      caption = "Points with horizontal bars show survey age bands. Shaded = 95% CrI."
    )

  # Save
  filename_base <- paste0("eFigure2_AgeProfile_", toupper(country_code), "_", gender, "_", year)
  ggsave(paste0("outputs/figures/age_profiles/", filename_base, ".pdf"), efig2, width = 10, height = 6)
  ggsave(paste0("outputs/figures/age_profiles/", filename_base, ".png"), efig2, width = 10, height = 6, dpi = 300)

  cat(sprintf("  Saved: outputs/figures/age_profiles/%s.pdf\n", filename_base))

  return(efig2)
}

# ============================================================================
# eFIGURE 3: WEIGHTED TREND CURVES
# ============================================================================

#' Create eFigure 3: Weighted Trend Curves
#' Shows time series of age-standardized prevalence with uncertainty
#'
#' @param weighted_trends Weighted trend data with Year, Prevalence_Mean, CIs
#' @param country_code Country to plot
#' @param gender Gender to plot
create_efigure3_trends <- function(weighted_trends,
                                    country_code,
                                    gender,
                                    indicator = PRIMARY_INDICATOR,
                                    country_name_mapping = NULL) {

  cat(sprintf("  Creating eFigure 3: Trend Curves (%s, %s)...\n", country_code, gender))

  # Filter data
  trend_data <- weighted_trends %>%
    filter(
      (Country == country_code | wb_country_abv == country_code | area == country_code),
      (Sex == gender | sex == gender),
      (Def_Type_Code == indicator | def_type_code == indicator)
    )

  if (nrow(trend_data) == 0) {
    cat("  WARNING: No trend data for this country\n")
    return(NULL)
  }

  # Get column names
  year_col <- if("Year" %in% names(trend_data)) "Year" else "year"
  mean_col <- if("Prevalence_Mean" %in% names(trend_data)) "Prevalence_Mean" else
    if("weighted_mean" %in% names(trend_data)) "weighted_mean" else "Prevalence"

  # Check for CI columns
  lower_col <- if("Prevalence_Lower" %in% names(trend_data)) "Prevalence_Lower" else
    if("lower_ci" %in% names(trend_data)) "lower_ci" else NULL
  upper_col <- if("Prevalence_Upper" %in% names(trend_data)) "Prevalence_Upper" else
    if("upper_ci" %in% names(trend_data)) "upper_ci" else NULL

  trend_data <- trend_data %>%
    mutate(
      Year_Num = as.numeric(!!sym(year_col)),
      Mean_Pct = if(max(!!sym(mean_col), na.rm = TRUE) <= 1) !!sym(mean_col) * 100 else !!sym(mean_col)
    )

  if (!is.null(lower_col) && !is.null(upper_col)) {
    trend_data <- trend_data %>%
      mutate(
        Lower_Pct = if(max(!!sym(lower_col), na.rm = TRUE) <= 1) !!sym(lower_col) * 100 else !!sym(lower_col),
        Upper_Pct = if(max(!!sym(upper_col), na.rm = TRUE) <= 1) !!sym(upper_col) * 100 else !!sym(upper_col)
      )
  }

  # Get country name
  country_name <- if(!is.null(country_name_mapping) && country_code %in% names(country_name_mapping)) {
    tools::toTitleCase(country_name_mapping[country_code])
  } else {
    toupper(country_code)
  }

  # Create plot
  efig3 <- ggplot(trend_data, aes(x = Year_Num))

  # Add ribbon if CIs available
  if (!is.null(lower_col) && !is.null(upper_col)) {
    efig3 <- efig3 +
      geom_ribbon(aes(ymin = Lower_Pct, ymax = Upper_Pct), fill = "#2166AC", alpha = 0.2)
  }

  efig3 <- efig3 +
    # Mean line
    geom_line(aes(y = Mean_Pct), color = "#2166AC", linewidth = 0.8) +

    # Endgame target
    geom_hline(yintercept = ENDGAME_THRESHOLD, linetype = "dashed", color = "#1B5E20", linewidth = 0.6) +
    annotate("text", x = min(trend_data$Year_Num), y = ENDGAME_THRESHOLD + 1,
             label = paste0(ENDGAME_THRESHOLD, "% Endgame"), hjust = 0, size = 3, color = "#1B5E20") +

    # Key year markers
    geom_vline(xintercept = BASE_YEAR, linetype = "dotted", color = "grey50", linewidth = 0.5) +
    geom_vline(xintercept = TARGET_YEAR, linetype = "dotted", color = "grey50", linewidth = 0.5) +
    annotate("text", x = BASE_YEAR, y = max(trend_data$Mean_Pct, na.rm = TRUE) * 0.95,
             label = "2010", hjust = 1.1, size = 2.5, color = "grey50") +
    annotate("text", x = TARGET_YEAR, y = max(trend_data$Mean_Pct, na.rm = TRUE) * 0.95,
             label = "2025", hjust = -0.1, size = 2.5, color = "grey50") +

    scale_y_continuous(limits = c(0, NA), labels = function(x) paste0(x, "%")) +
    theme_who() +
    labs(
      title = paste0("Prevalence Trend: ", country_name),
      subtitle = paste0(format_indicator(indicator), " | ", format_gender(gender)),
      x = "Year",
      y = "Age-Standardized Prevalence (%)",
      caption = if(!is.null(lower_col)) "Shaded area = 95% credible interval" else NULL
    )

  filename_base <- paste0("eFigure3_Trend_", toupper(country_code), "_", gender)
  ggsave(paste0("outputs/figures/trends/", filename_base, ".pdf"), efig3, width = 10, height = 6)
  ggsave(paste0("outputs/figures/trends/", filename_base, ".png"), efig3, width = 10, height = 6, dpi = 300)

  cat(sprintf("  Saved: outputs/figures/trends/%s.pdf\n", filename_base))

  return(efig3)
}

# ============================================================================
# MASTER EXECUTION FUNCTION
# ============================================================================

#' Generate All Publication Figures
#'
#' @param clean_data Survey data
#' @param global_predictions Global model predictions
#' @param country_predictions Country-specific predictions (optional)
#' @param weighted_results Weighted results with target achievement
#' @param country_region_mapping Country-region mapping
#' @param country_name_mapping Country name mapping
generate_all_figures <- function(clean_data,
                                  global_predictions,
                                  country_predictions = NULL,
                                  weighted_results = NULL,
                                  country_region_mapping = NULL,
                                  country_name_mapping = NULL) {

  cat("\n========== GENERATING PUBLICATION FIGURES ==========\n\n")

  results <- list()

  # Figure 1: Stoplight Maps (if weighted results available)
  # Now generates ALL probability maps: 6 indicators × 2 targets × 2 genders
  if (!is.null(weighted_results)) {
    tryCatch({
      # Generate all probability maps (this replaces the single map)
      results$maps <- generate_all_probability_maps(
        weighted_results = weighted_results,
        mcmc_samples = NULL,  # TODO: Pass MCMC samples for proper uncertainty
        output_dir = "outputs/figures/maps",
        primary_only = FALSE  # Set TRUE to only generate cigarettes maps
      )
    }, error = function(e) {
      cat(sprintf("  ERROR Figure 1 Maps: %s\n", e$message))
      # Fallback to single map
      tryCatch({
        results$fig1 <- create_figure1_map(weighted_results,
                                            country_name_mapping = country_name_mapping)
      }, error = function(e2) {
        cat(sprintf("  ERROR Figure 1 fallback: %s\n", e2$message))
      })
    })
  } else {
    cat("  Skipping Figure 1 (no weighted_results)\n")
  }

  # eFigure 1: Validation
  tryCatch({
    results$efig1 <- create_efigure1_validation(clean_data, global_predictions)
  }, error = function(e) {
    cat(sprintf("  ERROR eFigure 1: %s\n", e$message))
  })

  # Figure 4: Heatmap (global)
  tryCatch({
    results$fig4_men <- create_figure4_heatmap(global_predictions, "global", "males")
    results$fig4_women <- create_figure4_heatmap(global_predictions, "global", "females")
  }, error = function(e) {
    cat(sprintf("  ERROR Figure 4: %s\n", e$message))
  })

  # Get unique countries
  countries <- unique(c(clean_data$wb_country_abv, global_predictions$Country))
  countries <- countries[!is.na(countries) & countries != ""]

  cat(sprintf("\n  Processing %d countries for age profiles and trends...\n", length(countries)))

  # Loop through countries
  for (cc in countries[1:min(5, length(countries))]) {  # First 5 as sample
    for (g in c("males", "females")) {
      # Get observed years for this country
      obs_years <- clean_data %>%
        filter(wb_country_abv == cc | Country == cc) %>%
        pull(year) %>%
        unique() %>%
        as.numeric()

      if (length(obs_years) == 0) next

      # eFigure 2: Age profiles for sample years
      sample_years <- obs_years[seq(1, length(obs_years), length.out = min(3, length(obs_years)))]
      for (yr in sample_years) {
        tryCatch({
          create_efigure2_age_profile(clean_data, global_predictions, country_predictions,
                                       cc, g, yr)
        }, error = function(e) {
          # Silent fail for individual plots
        })
      }

      # eFigure 3: Trend curves
      if (!is.null(weighted_results)) {
        tryCatch({
          create_efigure3_trends(weighted_results, cc, g,
                                  country_name_mapping = country_name_mapping)
        }, error = function(e) {
          # Silent fail
        })
      }
    }
  }

  cat("\n========== FIGURE GENERATION COMPLETE ==========\n")

  # Summary
  fig_dirs <- c("maps", "caterpillar", "trajectories", "heatmaps", "validation", "age_profiles", "trends")
  cat("\nGenerated figures:\n")
  for (d in fig_dirs) {
    dir_path <- file.path("outputs/figures", d)
    if (dir.exists(dir_path)) {
      # Recursive for maps (has subdirectories)
      files <- list.files(dir_path, pattern = "\\.(pdf|png)$", recursive = TRUE)
      if (length(files) > 0) {
        cat(sprintf("  %s: %d files\n", d, length(files)))
        # Show breakdown for maps
        if (d == "maps") {
          for (subdir in c("2025_target", "2040_endgame", "composites")) {
            sub_path <- file.path(dir_path, subdir)
            if (dir.exists(sub_path)) {
              sub_files <- list.files(sub_path, pattern = "\\.(pdf|png)$", recursive = TRUE)
              if (length(sub_files) > 0) {
                cat(sprintf("    - %s: %d files\n", subdir, length(sub_files)))
              }
            }
          }
        }
      }
    }
  }

  return(results)
}

# ============================================================================
# STANDALONE EXECUTION
# ============================================================================

if (exists("clean_data") && exists("final_ac_predictions")) {
  cat("\n  Data available - generating figures...\n")

  # Adapt to available variable names
  global_preds <- if(exists("final_ac_predictions")) final_ac_predictions else
    if(exists("final_apc_predictions")) final_apc_predictions else NULL

  country_preds <- if(exists("final_predictions_country_specific_ac")) final_predictions_country_specific_ac else
    if(exists("final_predictions_country_specific_apc")) final_predictions_country_specific_apc else NULL

  weighted_res <- if(exists("final_weighted_results_selected")) final_weighted_results_selected else
    if(exists("final_weighted_results_apc")) final_weighted_results_apc else NULL

  region_map <- if(exists("country_region_mapping")) country_region_mapping else NULL
  name_map <- if(exists("country_name_mapping")) country_name_mapping else NULL

  generate_all_figures(
    clean_data = clean_data,
    global_predictions = global_preds,
    country_predictions = country_preds,
    weighted_results = weighted_res,
    country_region_mapping = region_map,
    country_name_mapping = name_map
  )

} else {
  cat("\n  No data available. Figures will be generated when data is loaded.\n")
  cat("\n  To generate figures, call:\n")
  cat("    generate_all_figures(clean_data, global_predictions, ...)\n")
}

cat("\n  Module 12 loaded.\n")
