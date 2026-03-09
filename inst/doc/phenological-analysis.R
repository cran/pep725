## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(pep725)

# Download the synthetic dataset
pep <- pep_download()

# For this vignette, we'll focus on two well-represented species:
# 1. Apple (Malus domestica) - excellent coverage across Europe
# 2. Grapevine (Vitis vinifera) - longest historical records

apple <- pep[species == "Malus domestica"]
vine <- pep[species == "Vitis vinifera"]

# Verify subsets contain data
stopifnot(nrow(apple) > 0, nrow(vine) > 0)

# Quick overview of our data
cat("=== Apple (Malus domestica) ===\n")
cat("Observations:", nrow(apple), "\n")
cat("Year range:", min(apple$year), "-", max(apple$year), "\n")
cat("Stations:", length(unique(apple$s_id)), "\n\n")

cat("=== Grapevine (Vitis vinifera) ===\n")
cat("Observations:", nrow(vine), "\n")
cat("Year range:", min(vine$year), "-", max(vine$year), "\n")
cat("Stations:", length(unique(vine$s_id)), "\n")

## ----normals-basic------------------------------------------------------------
# Calculate normals for apple across countries and phases
# Using 1990-2015 for illustration (ensures enough data in synthetic dataset)
normals <- pheno_normals(
  pep = apple,
  period = 1990:2015,
  by = c("country", "phase_id"),
  min_years = 10
)

print(normals)

## ----normals-summary----------------------------------------------------------
# Get a summary of the normals
summary(normals)

## ----normals-plot, fig.width=7, fig.height=5----------------------------------
plot(normals)

## ----normals-filtered---------------------------------------------------------
# Normals for apple first flowering (BBCH 60) only
flowering_normals <- pheno_normals(
  pep = pep,                          # Can use full dataset
  period = 1990:2015,
  species = "Malus domestica",         # Filter by species
  phase_id = 60,                      # Filter by phase (first flowers)
  by = "country",                     # Group only by country
  min_years = 5                       # Lower threshold for more countries
)

print(flowering_normals)

## ----normals-compare, warning=FALSE, message=FALSE----------------------------
# Calculate normals for two periods
# Period 1: Historical (1961-1990)
period1 <- pheno_normals(
  pep,
  period = 1961:1990,
  species = "Malus",
  phase_id = 60,
  by = "country",
  min_years = 10
)

# Period 2: Recent (1991-2020)
period2 <- pheno_normals(
  pep,
  period = 1991:2020,
  species = "Malus",
  phase_id = 60,
  by = "country",
  min_years = 10
)

# Compare the two periods by merging on country
# Negative shift = flowering got earlier
# Positive shift = flowering got later
if (nrow(period1) > 0 && nrow(period2) > 0) {
  comparison <- merge(
    period1[, .(country, mean_doy_early = mean_doy)],
    period2[, .(country, mean_doy_recent = mean_doy)],
    by = "country"
  )
  comparison[, shift := mean_doy_recent - mean_doy_early]

  cat("Countries compared:", nrow(comparison), "\n")
  cat("Mean shift:", round(mean(comparison$shift, na.rm = TRUE), 1), "days\n")
  cat("Range of shifts:", round(min(comparison$shift, na.rm = TRUE), 1), "to",
      round(max(comparison$shift, na.rm = TRUE), 1), "days\n")

  if (mean(comparison$shift, na.rm = TRUE) < 0) {
    cat("\nInterpretation: On average, flowering has shifted EARLIER\n")
  } else {
    cat("\nInterpretation: On average, flowering has shifted LATER\n")
  }
}

## ----normals-vine, warning=FALSE, message=FALSE-------------------------------
# Calculate normals for grapevine flowering
vine_normals <- pheno_normals(
  pep = vine,
  period = 1990:2015,
  phase_id = 65,  # Full flowering
  by = "country",
  min_years = 5
)

# Compare with apple for countries that have both
if (nrow(vine_normals) > 0) {
  cat("Grapevine flowering normals (BBCH 65):\n")
  print(vine_normals[!is.na(mean_doy), .(country, n_years, mean_doy, sd_doy)])
}

## ----anomalies-basic----------------------------------------------------------
# Calculate anomalies for apple
# Using 1990-2010 as baseline period
anomalies <- pheno_anomaly(
  pep = apple,
  baseline_period = 1990:2010,
  by = c("country", "phase_id"),
  min_years = 5
)

print(anomalies)

## ----anomalies-extreme--------------------------------------------------------
# Find extreme years (is_extreme == TRUE)
extreme <- anomalies[is_extreme == TRUE]

if (nrow(extreme) > 0) {
  cat("Total extreme years found:", nrow(extreme), "\n\n")

  # Show the most extreme early years
  cat("Most extreme EARLY years:\n")
  early <- extreme[direction == "early"][order(anomaly_days)]
  if (nrow(early) > 0) {
    print(early[1:min(5, nrow(early)),
          .(country, phase_id, year, anomaly_days, z_score)])
  }

  # Show the most extreme late years
  cat("\nMost extreme LATE years:\n")
  late <- extreme[direction == "late"][order(-anomaly_days)]
  if (nrow(late) > 0) {
    print(late[1:min(5, nrow(late)),
          .(country, phase_id, year, anomaly_days, z_score)])
  }
} else {
  cat("No extreme years detected in this dataset\n")
}

## ----anomalies-summary--------------------------------------------------------
summary(anomalies)

## ----anomalies-plot, fig.width=7, fig.height=5--------------------------------
plot(anomalies)

## ----quality-basic------------------------------------------------------------
# Assess quality at the station level
quality <- pep_quality(
  pep = apple,
  by = c("s_id", "phase_id")
)

print(quality)

## ----quality-summary----------------------------------------------------------
# Summary of quality across all stations
summary(quality)

## ----quality-plot, fig.height=5, fig.width=11---------------------------------
# Overview plot: grade distribution + station map (requires pep for coordinates)
plot(quality, pep = apple)

## ----quality-plot-grades, fig.height=4----------------------------------------
# Just the grade distribution (no pep data needed)
plot(quality, which = "grades")

## ----quality-plot-map, fig.height=5-------------------------------------------
# Just the map of station quality
plot(quality, which = "map", pep = apple)

## ----quality-filter-----------------------------------------------------------
# Get high-quality stations
high_quality <- quality[quality_grade %in% c("A", "B")]

cat("Quality distribution:\n")
print(table(quality$quality_grade))

cat("\nHigh quality (A or B):", nrow(high_quality), "of", nrow(quality),
    "(", round(100 * nrow(high_quality) / nrow(quality), 1), "%)\n")

# Filter original data to these stations
if (nrow(high_quality) > 0) {
  good_stations <- unique(high_quality$s_id)
  apple_hq <- apple[s_id %in% good_stations]

  cat("\nOriginal apple observations:", nrow(apple), "\n")
  cat("After quality filter:", nrow(apple_hq),
      "(", round(100 * nrow(apple_hq) / nrow(apple), 1), "%)\n")
}

## ----quality-country----------------------------------------------------------
# Country-level quality (coarser assessment)
country_quality <- pep_quality(
  pep = pep,
  by = c("country", "species", "phase_id")
)

print(country_quality)

## ----workflow-----------------------------------------------------------------
cat("=== STEP 1: Data Quality Assessment ===\n\n")

# Assess quality
quality <- pep_quality(apple, by = c("s_id", "phase_id"))
cat("Quality grades:\n")
print(table(quality$quality_grade))

cat("\n=== STEP 2: Filter to Good Quality Data ===\n\n")

# Keep grades A, B, and C (exclude only grade D)
good_stations <- quality[quality_grade %in% c("A", "B", "C"), s_id]
apple_clean <- apple[s_id %in% good_stations]

cat("Original observations:", nrow(apple), "\n")
cat("After quality filter:", nrow(apple_clean), "\n")
cat("Retained:", round(100 * nrow(apple_clean) / nrow(apple), 1), "%\n")

stopifnot(nrow(apple_clean) > 0)

cat("\n=== STEP 3: Calculate Baseline Normals ===\n\n")

normals <- pheno_normals(
  apple_clean,
  period = 1990:2010,
  by = c("country", "phase_id"),
  min_years = 5
)

cat("Normals calculated for", sum(!is.na(normals$mean_doy)), "groups\n")
cat("Mean DOY range:", round(min(normals$mean_doy, na.rm = TRUE), 0),
    "to", round(max(normals$mean_doy, na.rm = TRUE), 0), "\n")

cat("\n=== STEP 4: Calculate Anomalies ===\n\n")

anomalies <- pheno_anomaly(
  apple_clean,
  baseline_period = 1990:2010,
  by = c("country", "phase_id"),
  min_years = 5
)

# Convert to plain data.table so grouped aggregation works correctly
# (the pheno_anomaly S3 class can interfere with [.data.table dispatch)
anomalies_dt <- as.data.table(anomalies)

cat("Anomalies calculated for", sum(!is.na(anomalies_dt$anomaly_days)), "year-groups\n")
cat("Extreme years:", sum(anomalies_dt$is_extreme, na.rm = TRUE), "\n")

cat("\n=== STEP 5: Analyze Trends by Decade ===\n\n")

decade_summary <- anomalies_dt[
  !is.na(anomaly_days),
  .(
    mean_anomaly = round(mean(anomaly_days, na.rm = TRUE), 1),
    sd_anomaly = round(sd(anomaly_days, na.rm = TRUE), 1),
    n_obs = .N,
    pct_extreme = round(100 * sum(is_extreme, na.rm = TRUE) / .N, 1)
  ),
  by = .(decade = floor(year / 10) * 10)
][order(decade)]

print(decade_summary)

cat("\nInterpretation:\n")
first_decade <- decade_summary$mean_anomaly[1]
last_decade <- decade_summary$mean_anomaly[nrow(decade_summary)]
change <- last_decade - first_decade

if (change < -3) {
  cat("- Strong shift toward EARLIER phenology\n")
} else if (change < 0) {
  cat("- Moderate shift toward earlier phenology\n")
} else if (change > 3) {
  cat("- Strong shift toward LATER phenology\n")
} else {
  cat("- No clear directional change\n")
}

## ----timeseries-basic, fig.height=5-------------------------------------------
# Prepare aggregated data for visualization
apple_annual <- apple[, .(
  day = mean(day, na.rm = TRUE),
  n = .N
), by = .(year, species, phase_id)]

# Add human-readable phase labels
apple_annual[, phase := factor(phase_id,
  levels = c(60, 65, 100),
  labels = c("First flowers", "Full flowering", "Harvest"))]

# Remove rows where phase mapping failed
apple_annual <- apple_annual[!is.na(phase)]

# Plot time series with trend lines
if (nrow(apple_annual) > 0) {
  pheno_plot_timeseries(
    apple_annual,
    color_by = "phase",
    smooth = TRUE,
    title = "Apple Phenology Trends"
  )
}

## ----trends-turning, fig.height=5---------------------------------------------
# Get apple flowering data
apple_flowering <- pep[species == "Malus domestica" & phase_id == 60]

# Detect turning points
turning <- pheno_trend_turning(apple_flowering, min_years = 10)
print(turning)

# Visualize the analysis
plot(turning)

## ----kendall-simple-----------------------------------------------------------
# Aggregate to annual means
annual_doy <- apple_flowering[, .(
  mean_doy = mean(day, na.rm = TRUE)
), by = year][order(year)]

# Calculate Mann-Kendall test statistic
# Note: kendall_tau() returns the standardized z-statistic, not Kendall's
# tau correlation. Values beyond ±1.96 indicate significant trends.
tau <- kendall_tau(annual_doy$mean_doy)

cat("Mann-Kendall z-statistic:", round(tau, 3), "\n")
cat("\nInterpretation:\n")
if (tau < -2.58) {
  cat("- Highly significant earlier trend (p < 0.01)\n")
} else if (tau < -1.96) {
  cat("- Significant earlier trend (p < 0.05)\n")
} else if (tau > 2.58) {
  cat("- Highly significant later trend (p < 0.01)\n")
} else if (tau > 1.96) {
  cat("- Significant later trend (p < 0.05)\n")
} else {
  cat("- No statistically significant trend\n")
}

## ----regional-climate, eval=FALSE---------------------------------------------
# # Compile regional data and link to temperature anomalies
# regional_data <- pheno_regional(
#   pep = pep,
#   species_name = "Malus domestica",
#   year_min = 1961
# )
# 
# # Visualize phenology alongside temperature anomalies
# pheno_plot(regional_data)

## ----eval = FALSE-------------------------------------------------------------
# vignette("spatial-patterns", package = "pep725")

## ----eval = FALSE-------------------------------------------------------------
# vignette("data-quality", package = "pep725")

## ----session------------------------------------------------------------------
sessionInfo()

