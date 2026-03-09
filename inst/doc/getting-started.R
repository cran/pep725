## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----install, eval = FALSE----------------------------------------------------
# # Install from CRAN (when available)
# install.packages("pep725")
# 
# # Or install development version from GitHub
# # First install devtools if you don't have it:
# # install.packages("devtools")
# devtools::install_github("matthias-da/pep725")

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(pep725)

## ----download, message=FALSE--------------------------------------------------
# Download synthetic data (cached locally after first download)
pep <- pep_download()

## ----cachestatus, message=FALSE-----------------------------------------------
# Check cache status - shows file location and size
pep_cache_info()

# If you need to clear the cache (e.g., to get an updated version):
# pep_cache_clear()

## ----load-seed, eval = FALSE--------------------------------------------------
# # Load the included seed dataset (small, for quick tests only)
# data(pep_seed)
# 
# # This is a small subset (~1,300 rows) useful for:
# # - Testing code quickly
# # - Working offline
# # - Running examples in documentation

## ----simulate, eval = FALSE---------------------------------------------------
# # Generate synthetic data based on your real data
# # This preserves statistical properties while anonymizing observations
# pep_synthetic <- pep_simulate(your_real_pep_data)
# 
# # You can then share pep_synthetic freely

## ----import, eval = FALSE-----------------------------------------------------
# # Import from downloaded PEP725 files
# # Point to the folder containing your downloaded CSV files
# pep <- pep_import("path/to/pep725_data/")

## ----own-data, eval = FALSE---------------------------------------------------
# # Prepare your data with the required columns
# my_data <- data.frame(
#   s_id = c("MY001", "MY001", "MY002"),
#   lon = c(10.1, 10.1, 11.3),
#   lat = c(47.5, 47.5, 48.1),
#   genus = c("Malus", "Malus", "Malus"),
#   species = c("Malus domestica", "Malus domestica", "Malus domestica"),
#   phase_id = c(60, 65, 60),
#   year = c(2020, 2020, 2020),
#   day = c(110, 125, 108)
# )
# 
# # Convert to a pep object — validates required columns automatically
# my_pep <- as.pep(my_data)

## ----structure----------------------------------------------------------------
# View the data - the print method shows a summary
print(pep)

## ----first-look, fig.height = 4-----------------------------------------------
# Where are the stations?
plot(pep, type = "map")

# How are phenological events distributed across the year?
plot(pep, type = "histogram")

## ----class--------------------------------------------------------------------
# Check the class hierarchy
class(pep)

## ----summary------------------------------------------------------------------
# Default summary - shows top species by observation count
summary(pep)

## ----summary-phase------------------------------------------------------------
# Summary by phenological phase
summary(pep, by = "phase")

## ----summary-country----------------------------------------------------------
# Summary by country (showing top 5)
summary(pep, by = "country", top = 5)

## ----summary-year-------------------------------------------------------------
# Temporal coverage summary
summary(pep, by = "year")

## ----subset-------------------------------------------------------------------
# Filter to apple data using data.table syntax
# The syntax is: dataset[row_filter, column_selection]
apple <- pep[species == "Malus domestica"]
print(apple)

# Filter by year range
recent <- pep[year >= 2000 & year <= 2015]
print(recent)

# The summary method works on subsets too
summary(recent)

## ----bbch---------------------------------------------------------------------
# Get descriptions for BBCH codes present in your data
codes <- unique(pep$phase_id)
bbch_description(codes)

## ----filter-phase-------------------------------------------------------------
# Get all flowering observations (BBCH 60-69)
flowering <- pep[phase_id >= 60 & phase_id <= 69]
cat("Flowering observations:", nrow(flowering), "\n")

# Get only full flowering (BBCH 65)
full_flowering <- pep[phase_id == 65]
cat("Full flowering observations:", nrow(full_flowering), "\n")

# Get harvest-related observations
harvest <- pep[phase_id %in% c(87, 89)]
cat("Harvest observations:", nrow(harvest), "\n")

## ----coverage-all-------------------------------------------------------------
# Get a complete coverage report
pep_coverage(pep)

## ----coverage-temporal--------------------------------------------------------
# Temporal coverage with a visualization
pep_coverage(pep, kind = "temporal", plot = TRUE)

## ----coverage-geo-------------------------------------------------------------
# Geographical coverage - which countries have most data
pep_coverage(pep, kind = "geographical", top = 5)

## ----coverage-species---------------------------------------------------------
# Species coverage - which species are best represented
pep_coverage(pep, kind = "species", top = 5)

## ----coverage-by-country------------------------------------------------------
# Temporal coverage broken down by country
cov_by_country <- pep_coverage(pep, kind = "temporal", by = "country")
# Show observations by group
cov_by_country$temporal$obs_by_group

## ----quick-analysis-----------------------------------------------------------
# Step 1: Filter to apple flowering (BBCH 60 = first flowers)
apple_flowering <- pep[species == "Malus domestica" & phase_id == 60]

cat("Apple flowering observations:", nrow(apple_flowering), "\n")
cat("Year range:", min(apple_flowering$year), "-", max(apple_flowering$year), "\n")
cat("Countries:", length(unique(apple_flowering$country)), "\n")

# Step 2: Calculate annual mean DOY
# The data.table syntax here means:
#   - Group by year
#   - Calculate mean DOY and count observations
#   - Order by year
annual_mean <- apple_flowering[, .(
  mean_doy = mean(day, na.rm = TRUE),
  n_obs = .N
), by = year][order(year)]

# Step 3: Look at the results
print(annual_mean)

## ----quick-analysis-plot------------------------------------------------------
# Step 4: Plot the trend with ggplot2
library(ggplot2)

p <- ggplot(annual_mean, aes(x = year, y = mean_doy)) +
  geom_point(aes(size = n_obs), color = "steelblue", alpha = 0.6) +
  geom_line(color = "steelblue", alpha = 0.4) +
  geom_smooth(method = "lm", color = "red", linetype = 2, se = FALSE) +
  labs(
    title = "Apple Flowering Date Over Time",
    x = "Year", y = "Mean Day of Year",
    size = "Observations"
  ) +
  theme_minimal()
print(p)

# Step 5: Quantify the trend
trend <- lm(mean_doy ~ year, data = annual_mean)
slope <- coef(trend)[2]
cat("\nTrend:", round(slope * 10, 2), "days per decade\n")
if (slope < 0) {
  cat("Interpretation: Flowering is getting EARLIER\n")
} else {
  cat("Interpretation: Flowering is getting LATER\n")
}

## ----quick-analysis-vine------------------------------------------------------
# Step 1: Filter to grapevine flowering (BBCH 65 = full flowering)
vine_flowering <- pep[species == "Vitis vinifera" & phase_id == 65]

cat("Grapevine flowering observations:", nrow(vine_flowering), "\n")
cat("Year range:", min(vine_flowering$year), "-", max(vine_flowering$year), "\n")
cat("Countries:", length(unique(vine_flowering$country)), "\n")

# Step 2: Calculate annual mean DOY
vine_annual <- vine_flowering[, .(
  mean_doy = mean(day, na.rm = TRUE),
  n_obs = .N
), by = year][order(year)]

# Step 3: Plot the trend
p <- ggplot(vine_annual, aes(x = year, y = mean_doy)) +
  geom_point(aes(size = n_obs), color = "purple", alpha = 0.6) +
  geom_line(color = "purple", alpha = 0.4) +
  geom_smooth(method = "lm", color = "red", linetype = 2, se = FALSE) +
  labs(
    title = "Grapevine Flowering Date Over Time",
    x = "Year", y = "Mean Day of Year",
    size = "Observations"
  ) +
  theme_minimal()
print(p)

trend <- lm(mean_doy ~ year, data = vine_annual)
slope <- coef(trend)[2]
cat("\nTrend:", round(slope * 10, 2), "days per decade\n")

## ----eval = FALSE-------------------------------------------------------------
# vignette("phenological-analysis", package = "pep725")

## ----eval = FALSE-------------------------------------------------------------
# vignette("spatial-patterns", package = "pep725")

## ----eval = FALSE-------------------------------------------------------------
# vignette("data-quality", package = "pep725")

## ----session------------------------------------------------------------------
sessionInfo()

