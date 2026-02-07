################################################################################
# Test: GBRS vs MUGA genoprobs — sync, interpolate, correlate
# Run from project root with renv active (e.g. in RStudio or R in qtl2utilz).
# Usage: source("scripts/test_gbrs_muga_correlation.R")
################################################################################

# Paths: use project-relative paths if running from project root
data_dir <- "data"
if (!dir.exists(data_dir)) data_dir <- "~/work/qtlviewer-project/qtl2utilz/data"

genoprobs.GBRS    <- readRDS(file.path(data_dir, "genoprobs.GBRS.rds"))
genoprobs.MUGA_ALL <- readRDS(file.path(data_dir, "genoprobs.MUGA_ALL.rds"))
K.GBRS            <- readRDS(file.path(data_dir, "K.GBRS.rds"))
K.MUGA_ALL        <- readRDS(file.path(data_dir, "K.MUGA_ALL.rds"))
map.GBRS          <- readRDS(file.path(data_dir, "map.GBRS.rds"))
map.MUGA_ALL      <- readRDS(file.path(data_dir, "map.MUGA_ALL.rds"))
markers.GBRS      <- readRDS(file.path(data_dir, "markers.GBRS.rds"))
markers.MUGA_ALL  <- readRDS(file.path(data_dir, "markers.MUGA_ALL.rds"))

# 1) Synchronize markers for each genoprobs (subset to marker tables, consistent order)
gbrs <- qtl2utilz::genoprobs_sync_markers(genoprobs.GBRS, markers.GBRS)
muga <- qtl2utilz::genoprobs_sync_markers(genoprobs.MUGA_ALL, markers.MUGA_ALL)

# 2) Interpolate GBRS onto MUGA marker grid (same positions/markers as muga$map)
gbrs_on_muga <- qtl2utilz::genoprobs_interpolate(gbrs$genoprobs, gbrs$map, muga$map)

# 3) Restrict to common samples (same order in both)
new <- qtl2utilz::genoprobs_sync_samples(gbrs_on_muga, muga$genoprobs)

# 4) Correlate GBRS-on-MUGA vs native MUGA per sample
cor_tbl <- qtl2utilz::genoprobs_correlate(new$genoprobs_1, new$genoprobs_2)

# --- Print results ---
cat("--- 1. Sync markers ---\n")
cat("gbrs: ", length(gbrs$genoprobs), " chrs, ", nrow(gbrs$markers), " markers\n", sep = "")
cat("muga: ", length(muga$genoprobs), " chrs, ", nrow(muga$markers), " markers\n", sep = "")

cat("\n--- 2. Interpolate GBRS onto MUGA grid ---\n")
cat("gbrs_on_muga: ", length(gbrs_on_muga), " chrs\n", sep = "")

cat("\n--- 3. Sync samples ---\n")
cat("common_samples: ", length(new$common_samples), "\n", sep = "")
cat("dropped_from_1: ", length(new$dropped_from_1), ", dropped_from_2: ", length(new$dropped_from_2), "\n", sep = "")

cat("\n--- 4. Correlate ---\n")
cat("cor_tbl: ", nrow(cor_tbl), " rows\n", sep = "")
cat("\nSummary of correlations:\n")
print(summary(cor_tbl$correlation))
cat("\nFlag_low (cor < 0.8): ", sum(cor_tbl$flag_low), " / ", nrow(cor_tbl), "\n", sep = "")
cat("\nFirst 10 rows of cor_tbl:\n")
print(head(cor_tbl, 10))
cat("\nLast 5 rows:\n")
print(tail(cor_tbl, 5))

# Optional: quick histogram
if (requireNamespace("graphics", quietly = TRUE)) {
  hist(cor_tbl$correlation, main = "GBRS vs MUGA genoprobs correlation", xlab = "Correlation")
}
