# Data Files

This folder should contain the following data files (not included in repository due to size/licensing):

## Required Files

### Primary Survey Data
- `Prevalence RESHAPE 28 Sep 2024.xlsx` - WHO tobacco prevalence survey data
  - Source: WHO Comprehensive Information System for Tobacco Control (CiT)
  - Contact WHO for access

### Population Weights
- `weights_15_2022.csv` - Age-sex-year population weights for 191 countries
  - Source: WHO CiT
  - Used for computing population-weighted national prevalence

## File Format Requirements

### Survey Data (Excel)
Required columns:
| Column | Description |
|--------|-------------|
| `wb_country_abv` | ISO 3166-1 alpha-3 country code (lowercase) |
| `survey_year` | Year survey was conducted |
| `sex_name` | `male` or `female` |
| `start_age`, `end_age` | Age band boundaries |
| `prevalence` | Point estimate (0-100 scale) |
| `def_code` | Definition code (e.g., `cd.0101` for daily, `cd.0112` for current) |
| `type_code` | Product type (e.g., `tt.003` for cigarettes) |

### Population Weights (CSV)
Required columns:
- Country identifier
- Year
- Age group
- Sex
- Population count

## Data Access

To obtain the required data files, contact:

**WHO Comprehensive Information System for Tobacco Control**
https://www.who.int/teams/health-promotion/tobacco-control/who-report-on-the-global-tobacco-epidemic
