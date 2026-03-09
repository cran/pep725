## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(pep725)
library(data.table)
library(ggplot2)

# Download the synthetic dataset
pep <- pep_download()

# Overview of the data
cat("Stations:", length(unique(pep$s_id)), "\n")
cat("Altitude range:", min(pep$alt, na.rm = TRUE), "-",
    max(pep$alt, na.rm = TRUE), "m\n")
cat("Latitude range:", round(min(pep$lat, na.rm = TRUE), 2), "-",
    round(max(pep$lat, na.rm = TRUE), 2), "N\n")

## ----gradient-altitude--------------------------------------------------------
# Calculate elevational gradient for apple flowering
grad_alt <- pheno_gradient(
  pep = pep,
  variable = "alt",
  species = "Malus domestica",
  phase_id = 60,
  method = "robust"
)

print(grad_alt)

## ----gradient-latitude--------------------------------------------------------
# Calculate latitude gradient
grad_lat <- pheno_gradient(
  pep = pep,
  variable = "lat",
  species = "Malus domestica",
  phase_id = 60,
  method = "robust"
)

print(grad_lat)

## ----gradient-methods---------------------------------------------------------
# Robust regression (default) - handles outliers
grad_robust <- pheno_gradient(pep, variable = "alt",
                              species = "Malus domestica",
                              phase_id = 60, method = "robust")

# Ordinary least squares
grad_ols <- pheno_gradient(pep, variable = "alt",
                           species = "Malus domestica",
                           phase_id = 60, method = "ols")

cat("Robust slope:", round(grad_robust$summary$slope, 2), "days/100m\n")
cat("OLS slope:", round(grad_ols$summary$slope, 2), "days/100m\n")

## ----gradient-vine------------------------------------------------------------
# Calculate gradient for grapevine
vine_data <- pep[species == "Vitis vinifera"]

# Verify data exists before analysis
if (nrow(vine_data) > 0) {
  grad_vine <- pheno_gradient(
    pep = pep,
    variable = "alt",
    species = "Vitis vinifera",
    phase_id = 65,  # Full flowering
    method = "robust"
  )

  cat("Comparison of elevation gradients:\n")
  cat("Apple:     ", round(grad_robust$summary$slope, 2), "days/100m",
      "(R² =", round(grad_robust$summary$r_squared, 3), ")\n")
  cat("Grapevine: ", round(grad_vine$summary$slope, 2), "days/100m",
      "(R² =", round(grad_vine$summary$r_squared, 3), ")\n")
} else {
  cat("No grapevine data available for gradient analysis.\n")
}

## ----gradient-by-region-------------------------------------------------------
# Calculate gradient by country
grad_by_country <- pheno_gradient(
  pep = pep,
  variable = "alt",
  species = "Malus domestica",
  phase_id = 60,
  by = "country",
  method = "robust"
)

print(grad_by_country$summary)

## ----gradient-plot, fig.height=5----------------------------------------------
# Plot the gradient relationship
if (!is.null(grad_alt$data) && nrow(grad_alt$data) > 2) {
  p <- plot(grad_alt)
  print(p)
}

## ----synchrony-basic----------------------------------------------------------
# Calculate synchrony by country and year
sync <- pheno_synchrony(
  pep = pep,
  species = "Malus domestica",
  phase_id = 60,
  by = c("country", "year"),
  min_stations = 3
)

print(sync)

## ----synchrony-summary--------------------------------------------------------
summary(sync)

## ----synchrony-trend----------------------------------------------------------
# Check if trend analysis was performed
if (!is.null(sync$trend) && nrow(sync$trend) > 0) {
  cat("Trend Analysis Results:\n")
  print(sync$trend)

  # Interpret significant trends
  sig_trends <- sync$trend[!is.na(significant) & significant == TRUE]
  if (nrow(sig_trends) > 0) {
    cat("\nSignificant trends detected:\n")
    print(sig_trends[, .(country, slope, direction, p_value)])
  }
}

## ----synchrony-plot, fig.height=5---------------------------------------------
# Plot synchrony time series
if (nrow(sync$data[!is.na(sd_doy)]) > 5) {
  p <- plot(sync)
  print(p)
}

## ----synchrony-simple---------------------------------------------------------
sync_simple <- pheno_synchrony(
  pep = pep,
  species = "Malus domestica",
  phase_id = 60,
  by = c("country", "year"),
  min_stations = 3,
  compute_trend = FALSE
)

# Just the synchrony data
head(sync_simple$data)

## ----synchrony-compare--------------------------------------------------------
# Calculate synchrony for grapevine
vine_check <- pep[species == "Vitis vinifera" & phase_id == 65]

if (nrow(vine_check) > 0) {
  sync_vine <- pheno_synchrony(
    pep = pep,
    species = "Vitis vinifera",
    phase_id = 65,
    by = c("country", "year"),
    min_stations = 3,
    compute_trend = FALSE
  )

  cat("Synchrony comparison (mean SD across years):\n")
  cat("Apple (flowering):    ", round(sync$overall$mean_sd_doy, 1), "days\n")
  cat("Grapevine (flowering):", round(sync_vine$overall$mean_sd_doy, 1), "days\n")
} else {
  cat("Insufficient grapevine data for synchrony comparison.\n")
}

## ----combined-analysis--------------------------------------------------------
# 1. Understand elevation gradient
gradient <- pheno_gradient(
  pep, variable = "alt",
  species = "Malus domestica",
  phase_id = 60
)

cat("Elevation Gradient Analysis:\n")
cat("  Slope:", round(gradient$summary$slope, 2), "days/100m\n")
cat("  R-squared:", round(gradient$summary$r_squared, 3), "\n\n")

# 2. Assess spatial synchrony
synchrony <- pheno_synchrony(
  pep,
  species = "Malus domestica",
  phase_id = 60,
  min_stations = 3
)

cat("Synchrony Analysis:\n")
cat("  Mean SD across stations:", round(synchrony$overall$mean_sd_doy, 1), "days\n")
cat("  Mean CV:", round(synchrony$overall$mean_cv_pct, 1), "%\n")

## ----leaflet-example, eval=FALSE----------------------------------------------
# # Launch interactive map (opens in viewer)
# # Draw rectangles or polygons to select stations
# selected_stations <- pheno_leaflet(pep, label_col = "species")
# 
# # The function returns a data.frame of selected stations
# print(selected_stations)

## ----leaflet-screenshot, echo=FALSE, out.width="100%", fig.cap="The pheno_leaflet() interface showing station locations with drawing tools for selection."----
# Display screenshot if it exists
screenshot_path <- "pheno_leaflet_screenshot.png"
if (file.exists(screenshot_path)) {
  knitr::include_graphics(screenshot_path)
} else {
  # Placeholder message when screenshot doesn't exist
  cat("*[Screenshot: Interactive map showing PEP725 stations across Europe with drawing toolbar for rectangle and polygon selection]*\n")
}

## ----leaflet-filtered, eval=FALSE---------------------------------------------
# # Filter to a specific species for faster loading (recommended)
# apple <- pep[species == "Malus domestica"]
# selected <- pheno_leaflet(apple, label_col = "s_id")
# 
# # Use selected stations for focused analysis
# apple_subset <- apple[s_id %in% selected$s_id]

## ----map-basic, fig.height=6, fig.width=8-------------------------------------
# Basic station map with country borders (no API key needed)
pheno_map(pep, background = "none", color_by = "none", point_size = 1.5)

## ----map-nobs, fig.height=6, fig.width=8--------------------------------------
# Color stations by number of observations
pheno_map(pep, background = "none", color_by = "n_obs", point_size = 1.5)

## ----map-nspecies, fig.height=6, fig.width=8----------------------------------
# Color by species diversity at each station
pheno_map(pep, background = "none", color_by = "n_species", point_size = 1.5)

## ----map-google, eval=FALSE---------------------------------------------------
# # Register your API key first
# ggmap::register_google(key = "your_api_key_here")
# 
# # Then use background = "google"
# pheno_map(pep, background = "google", color_by = "n_obs", zoom = 5)
# 
# # Regional focus with Google Maps
# pheno_map(
#   pep,
#   background = "google",
#   location = c(lon = 8.2, lat = 46.8),
#   zoom = 7,
#   color_by = "n_obs"
# )

## ----map-mean-doy, fig.height=6, fig.width=8----------------------------------
# Map mean phenological timing for flowering (phase 60)
# Earlier timing = darker colors, later = lighter (plasma scale)
pheno_map(pep, background = "none", color_by = "mean_doy",
        phase_id = 60, point_size = 1.5)

## ----map-trend, fig.height=6, fig.width=8-------------------------------------
# Map trends per station
# Blue = phenology getting earlier, Red = getting later
pheno_map(pep, background = "none", color_by = "trend",
        phase_id = 60, period = 1990:2020, min_years = 10, point_size = 1.5)

## ----map-species-cv, fig.height=6, fig.width=8--------------------------------
# Map species variation at each station
# Higher CV = more variation in timing among species
pheno_map(pep, background = "none", color_by = "species_cv",
        phase_id = 60, min_species = 3, point_size = 1.5)

## ----mapping-workflow, eval=FALSE---------------------------------------------
# # 1. Explore stations interactively
# selected <- pheno_leaflet(pep[genus == "Malus"])
# 
# # 2. Subset data to selected region
# pep_region <- pep[s_id %in% selected$s_id]
# 
# # 3. Analyze gradients in the selected region
# gradient <- pheno_gradient(pep_region, variable = "alt",
#                            species = "Malus domestica", phase_id = 60)
# 
# # 4. Assess synchrony
# synchrony <- pheno_synchrony(pep_region, species = "Malus domestica",
#                               phase_id = 60)
# 
# # 5. Create publication map of the region (no API needed)
# pheno_map(pep_region, background = "none", color_by = "mean_doy",
#         phase_id = 60, output_file = "regional_phenology.pdf")

## ----eval = FALSE-------------------------------------------------------------
# vignette("getting-started", package = "pep725")

## ----eval = FALSE-------------------------------------------------------------
# vignette("phenological-analysis", package = "pep725")

## ----eval = FALSE-------------------------------------------------------------
# vignette("data-quality", package = "pep725")

## ----session------------------------------------------------------------------
sessionInfo()

