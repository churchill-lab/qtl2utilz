#' qtl2utilz: Utilities for R/qtl2 and GBRS Workflows
#'
#' Common functions and methods for working with R/qtl2 and GBRS (Genotype By
#' RNA-Seq) data. Enables importing GBRS genoprobs and read counts into R/qtl2
#' format, synchronizing and interpolating genoprobs across marker sets, and
#' comparing GBRS vs array-based genotypes (e.g. MUGA).
#'
#' @section GBRS workflow:
#' \code{\link{gbrs_find_files}}, \code{\link{gbrs_build_genoprobs}},
#' \code{\link{gbrs_build_counts}}
#'
#' @section Genoprobs utilities:
#' \code{\link{genoprobs_sync}}, \code{\link{genoprobs_sync_markers}},
#' \code{\link{genoprobs_sync_samples}}, \code{\link{genoprobs_intersect_samples}},
#' \code{\link{genoprobs_combine_samples}}, \code{\link{genoprobs_compare}},
#' \code{\link{genoprobs_interpolate}}, \code{\link{genoprobs_correlate}},
#' \code{\link{genoprobs_compute_similarity}}, \code{\link{genoprobs_detect_sample_swaps}}
#'
#' @section Other utilities:
#' \code{\link{positions_to_bp}}, \code{\link{markers_sort}}, \code{\link{rankZ}}
#'
"_PACKAGE"
