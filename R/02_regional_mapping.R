#########################################################################################
#
#                    WHO TOBACCO CONTROL PREVALENCE PROJECTION MODEL
#                       02_regional_mapping.R - Regional Classification
#
#   Contains: WHO regional group definitions and country-to-region mapping
#   Requires: 00_config.R, 01_data_prep.R to be sourced first
#   Outputs: country_region_manual, country_region_mapping, updates clean_data
#
#########################################################################################

# ---- 2.1 Define WHO Regional Groups ----

country_region_manual <- tribble(
  ~wb_country_abv, ~region_consolidated,

  # Central Asia
  "tjk", "Central Asia", "tkm", "Central Asia", "kaz", "Central Asia",

  # Eastern Asia
  "chn", "Eastern Asia", "prk", "Eastern Asia", "jpn", "Eastern Asia",
  "mng", "Eastern Asia", "kor", "Eastern Asia",
  # Countries without survey data (will use regional priors):
  "twn", "Eastern Asia",   # Taiwan
  "hkg", "Eastern Asia",   # Hong Kong SAR
  "mac", "Eastern Asia",   # Macao SAR

  # Eastern Europe
  "blr", "Eastern Europe", "bgr", "Eastern Europe", "cze", "Eastern Europe",
  "hun", "Eastern Europe", "mda", "Eastern Europe", "pol", "Eastern Europe",
  "rou", "Eastern Europe", "rus", "Eastern Europe", "svk", "Eastern Europe",
  "ukr", "Eastern Europe",

  # Latin America & Caribbean
  "bhs", "Latin America & Caribbean", "brb", "Latin America & Caribbean",
  "blz", "Latin America & Caribbean", "cub", "Latin America & Caribbean",
  "dma", "Latin America & Caribbean", "slv", "Latin America & Caribbean",
  "grd", "Latin America & Caribbean", "gtm", "Latin America & Caribbean",
  "hti", "Latin America & Caribbean", "hnd", "Latin America & Caribbean",
  "jam", "Latin America & Caribbean", "mex", "Latin America & Caribbean",
  "nic", "Latin America & Caribbean", "kna", "Latin America & Caribbean",
  "lca", "Latin America & Caribbean", "vct", "Latin America & Caribbean",
  "tto", "Latin America & Caribbean", "dom", "Latin America & Caribbean",
  "cri", "Latin America & Caribbean", "pan", "Latin America & Caribbean",

  # North Africa & Middle East
  "dza", "North Africa & Middle East", "egy", "North Africa & Middle East",
  "lby", "North Africa & Middle East", "mar", "North Africa & Middle East",
  "sdn", "North Africa & Middle East", "tun", "North Africa & Middle East",
  "pse", "North Africa & Middle East",

  # North America
  "usa", "North America", "can", "North America",

  # Oceania & Pacific
  "kir", "Oceania & Pacific", "mhl", "Oceania & Pacific", "fsm", "Oceania & Pacific",
  "nru", "Oceania & Pacific", "niu", "Oceania & Pacific", "plw", "Oceania & Pacific",
  "png", "Oceania & Pacific", "wsm", "Oceania & Pacific", "slb", "Oceania & Pacific",
  "ton", "Oceania & Pacific", "vut", "Oceania & Pacific", "cok", "Oceania & Pacific",
  "tuv", "Oceania & Pacific", "nzl", "Oceania & Pacific", "aus", "Oceania & Pacific",
  "fji", "Oceania & Pacific",

  # South America
  "arg", "South America", "bol", "South America", "bra", "South America",
  "chl", "South America", "col", "South America", "ecu", "South America",
  "guy", "South America", "per", "South America", "sur", "South America",
  "ury", "South America", "ven", "South America", "pry", "South America",

  # Southeastern Asia
  "brn", "Southeastern Asia", "idn", "Southeastern Asia", "lao", "Southeastern Asia",
  "mys", "Southeastern Asia", "mmr", "Southeastern Asia", "sgp", "Southeastern Asia",
  "tha", "Southeastern Asia", "vnm", "Southeastern Asia", "phl", "Southeastern Asia",
  "tls", "Southeastern Asia", "khm", "Southeastern Asia",

  # Southern Asia
  "afg", "Southern Asia", "bgd", "Southern Asia", "ind", "Southern Asia",
  "irn", "Southern Asia", "kgz", "Southern Asia", "mdv", "Southern Asia",
  "pak", "Southern Asia", "lka", "Southern Asia", "uzb", "Southern Asia",
  "npl", "Southern Asia", "btn", "Southern Asia",

  # Southern Europe
  "alb", "Southern Europe", "and", "Southern Europe", "bih", "Southern Europe",
  "grc", "Southern Europe", "hrv", "Southern Europe", "ita", "Southern Europe",
  "mlt", "Southern Europe", "mne", "Southern Europe", "mkd", "Southern Europe",
  "prt", "Southern Europe", "smr", "Southern Europe", "srb", "Southern Europe",
  "svn", "Southern Europe", "esp", "Southern Europe",
  # Countries without survey data (will use regional priors):
  "xkx", "Southern Europe",  # Kosovo

  # Sub-Saharan Africa
  "ago", "Sub-Saharan Africa", "ben", "Sub-Saharan Africa", "bwa", "Sub-Saharan Africa",
  "bfa", "Sub-Saharan Africa", "bdi", "Sub-Saharan Africa", "cpv", "Sub-Saharan Africa",
  "cmr", "Sub-Saharan Africa", "caf", "Sub-Saharan Africa", "tcd", "Sub-Saharan Africa",
  "com", "Sub-Saharan Africa", "cog", "Sub-Saharan Africa", "cod", "Sub-Saharan Africa",
  "civ", "Sub-Saharan Africa", "dji", "Sub-Saharan Africa", "gnq", "Sub-Saharan Africa",
  "eri", "Sub-Saharan Africa", "eth", "Sub-Saharan Africa", "gab", "Sub-Saharan Africa",
  "gmb", "Sub-Saharan Africa", "gha", "Sub-Saharan Africa", "gin", "Sub-Saharan Africa",
  "gnb", "Sub-Saharan Africa", "ken", "Sub-Saharan Africa", "lso", "Sub-Saharan Africa",
  "lbr", "Sub-Saharan Africa", "mdg", "Sub-Saharan Africa", "mwi", "Sub-Saharan Africa",
  "mli", "Sub-Saharan Africa", "mrt", "Sub-Saharan Africa", "mus", "Sub-Saharan Africa",
  "moz", "Sub-Saharan Africa", "nam", "Sub-Saharan Africa", "ner", "Sub-Saharan Africa",
  "nga", "Sub-Saharan Africa", "rwa", "Sub-Saharan Africa", "stp", "Sub-Saharan Africa",
  "sen", "Sub-Saharan Africa", "syc", "Sub-Saharan Africa", "sle", "Sub-Saharan Africa",
  "zaf", "Sub-Saharan Africa", "swz", "Sub-Saharan Africa", "tza", "Sub-Saharan Africa",
  "tgo", "Sub-Saharan Africa", "uga", "Sub-Saharan Africa", "zmb", "Sub-Saharan Africa",
  "zwe", "Sub-Saharan Africa",
  # Countries without survey data (will use regional priors):
  "ssd", "Sub-Saharan Africa",  # South Sudan
  "som", "Sub-Saharan Africa",  # Somalia

  # Western Asia
  "arm", "Western Asia", "aze", "Western Asia", "bhr", "Western Asia",
  "cyp", "Western Asia", "geo", "Western Asia", "irq", "Western Asia",
  "isr", "Western Asia", "jor", "Western Asia", "kwt", "Western Asia",
  "lbn", "Western Asia", "omn", "Western Asia", "qat", "Western Asia",
  "sau", "Western Asia", "syr", "Western Asia", "tur", "Western Asia",
  "are", "Western Asia", "yem", "Western Asia",

  # Western Europe
  "aut", "Western Europe", "bel", "Western Europe", "deu", "Western Europe",
  "fra", "Western Europe", "gbr", "Western Europe", "irl", "Western Europe",
  "nld", "Western Europe", "che", "Western Europe", "lux", "Western Europe",

  # Northern Europe
  "dnk", "Northern Europe", "fin", "Northern Europe", "isl", "Northern Europe",
  "nor", "Northern Europe", "swe", "Northern Europe", "est", "Northern Europe",
  "lva", "Northern Europe", "ltu", "Northern Europe"
)

# ---- 2.2 Map Countries to Regions ----

country_region_mapping <- clean_data %>%
  select(wb_country_abv) %>%
  distinct() %>%
  left_join(country_region_manual, by = "wb_country_abv") %>%
  mutate(region_consolidated = if_else(is.na(region_consolidated), "Other", region_consolidated))

clean_data <- clean_data %>%
  select(-region) %>%
  left_join(country_region_mapping, by = "wb_country_abv") %>%
  rename(region = region_consolidated)

# ---- 2.3 Export Regional Mapping ----

write.csv(
  country_region_manual %>%
    mutate(Country_Name = country_name_mapping[wb_country_abv]) %>%
    select(Country_Code = wb_country_abv, Country_Name, Region = region_consolidated) %>%
    arrange(Region, Country_Name),
  file = "outputs/data_exports/country_region_mapping.csv",
  row.names = FALSE
)

cat("\nRegional mapping complete.\n")
cat(sprintf("  Regions: %d\n", length(unique(country_region_mapping$region_consolidated))))
