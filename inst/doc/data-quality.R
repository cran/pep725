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

# Download the synthetic dataset
pep <- pep_download()

# For this vignette, we'll focus on flowering phases
flowering <- pep[phase_id %in% c(60, 65)]
cat("Flowering observations:", nrow(flowering), "\n")
cat("Species:", length(unique(flowering$species)), "\n")
cat("Year range:", min(flowering$year), "-", max(flowering$year), "\n")

## ----flag-basic---------------------------------------------------------------
# Flag outliers using the default 30-day rule
outliers <- pep_flag_outliers(
  pep = flowering,
  method = "30day",
  by = c("species", "phase_id")
)

print(outliers)

## ----flag-methods, eval=FALSE-------------------------------------------------
# # MAD method: flag if > 3 MAD from median (robust)
# outliers_mad <- pep_flag_outliers(flowering, method = "mad", threshold = 3)

## ----iqr-method, eval=FALSE---------------------------------------------------
# # IQR method: flag if outside 1.5 * IQR (standard boxplot rule)
# outliers_iqr <- pep_flag_outliers(flowering, method = "iqr", threshold = 1.5)

## ----zscore-method, eval=FALSE------------------------------------------------
# # Z-score method: flag if |z| > 3 (assumes normal distribution)
# outliers_zscore <- pep_flag_outliers(flowering, method = "zscore", threshold = 3)

## ----flag-summary-------------------------------------------------------------
summary(outliers)

## ----flag-compare-------------------------------------------------------------
# Check outliers for grapevine
vine_flowering <- pep[species == "Vitis vinifera" & phase_id %in% c(60, 65)]

if (nrow(vine_flowering) > 0) {
  outliers_vine <- pep_flag_outliers(
    pep = vine_flowering,
    method = "30day",
    by = c("species", "phase_id")
  )

  cat("Outlier comparison:\n")
  cat("All flowering species: ", round(100 * mean(outliers$is_outlier), 2), "%\n")
  cat("Grapevine only:        ", round(100 * mean(outliers_vine$is_outlier), 2), "%\n")
} else {
  cat("No grapevine flowering data available for comparison.\n")
}

## ----plot-overview, fig.height=6----------------------------------------------
# Overview of outlier patterns
pep_plot_outliers(outliers, type = "overview")

## ----plot-seasonal, fig.height=5----------------------------------------------
# When do outliers occur in the year?
pep_plot_outliers(outliers, type = "seasonal")

## ----plot-detail, fig.height=5------------------------------------------------
# See outliers in context of all observations
pep_plot_outliers(outliers, type = "detail", n_top = 15)

## ----plot-map, fig.height=5---------------------------------------------------
# Where are outliers located?
pep_plot_outliers(outliers, type = "map")

## ----completeness-basic-------------------------------------------------------
# Check completeness by station and phase
# Use year_range to focus on a specific period
completeness <- pep_completeness(
  pep = flowering,
  by = c("s_id", "phase_id"),
  year_range = c(1990, 2020)
)

print(completeness)

## ----completeness-summary-----------------------------------------------------
summary(completeness)

## ----completeness-filter------------------------------------------------------
# Get stations with good coverage (>= 70%)
good_coverage <- completeness[completeness_pct >= 70]
cat("Stations with >= 70% coverage:", nrow(good_coverage), "\n")

# Use these for trend analysis
good_stations <- unique(good_coverage$s_id)
flowering_complete <- flowering[s_id %in% good_stations]
cat("Observations from complete stations:", nrow(flowering_complete), "\n")

## ----completeness-plot, fig.height=5------------------------------------------
plot(completeness)

## ----check-basic--------------------------------------------------------------
# Check if expected phases are present for apple
apple <- pep[species == "Malus domestica"]

phase_check <- pep_check_phases(
  pep = apple,
  expected = c(60, 65, 87)  # flowering, full flowering, fruit maturity
)

print(phase_check)

## ----check-multi--------------------------------------------------------------
# Check phases for multiple species
multi_check <- pep_check_phases_multi(
  pep = pep,
  species_list = c("Malus domestica", "Vitis vinifera"),
  expected = c(60, 65, 87)
)

print(multi_check)

## ----workflow-----------------------------------------------------------------
# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: Assess temporal completeness
# ══════════════════════════════════════════════════════════════════════════════
# Why: Incomplete stations can bias trend estimates and normals

completeness <- pep_completeness(flowering, by = c("s_id", "phase_id"))
good_stations <- completeness[completeness_pct >= 50, s_id]
fl_filtered <- flowering[s_id %in% good_stations]

cat("Kept", length(good_stations), "stations with >= 50% completeness\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: Flag statistical outliers
# ══════════════════════════════════════════════════════════════════════════════
# Why: Identify observations that deviate from expected patterns

outliers_wf <- pep_flag_outliers(fl_filtered, method = "mad", threshold = 3)

cat("Flagged", sum(outliers_wf$is_outlier), "outliers",
    "(", round(100 * mean(outliers_wf$is_outlier), 1), "% )\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: Make informed decisions about exclusion
# ══════════════════════════════════════════════════════════════════════════════
# Key principle: Document your decisions!

# Option A: Strict cleaning (for normals calculation)
fl_strict <- outliers_wf[is_outlier == FALSE]

# Option B: Moderate cleaning (keep moderate outliers)
fl_moderate <- outliers_wf[is_outlier == FALSE | abs(deviation) < 60]

cat("Strict cleaning keeps:", nrow(fl_strict), "obs\n")
cat("Moderate cleaning keeps:", nrow(fl_moderate), "obs\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4: Proceed with analysis on cleaned data
# ══════════════════════════════════════════════════════════════════════════════

normals <- pheno_normals(fl_moderate, period = 1991:2020, min_years = 5)
cat("Normals calculated for", sum(!is.na(normals$mean_doy)), "groups\n")

## ----eval = FALSE-------------------------------------------------------------
# vignette("getting-started", package = "pep725")

## ----eval = FALSE-------------------------------------------------------------
# vignette("phenological-analysis", package = "pep725")

## ----eval = FALSE-------------------------------------------------------------
# vignette("spatial-patterns", package = "pep725")

## ----session------------------------------------------------------------------
sessionInfo()

