#' Detect duplicate samples within a single genoprobs object
#'
#' Computes pairwise similarity among all samples in a \code{calc_genoprob}
#' object and flags pairs above a threshold as potential duplicates.
#' This complements \code{qtl2::compare_geno()} by working directly on
#' allele-probability genoprobs (which may be all that is available when raw
#' marker genotypes have already been processed away).
#'
#' @section Why this matters:
#' Duplicate samples (two tubes containing DNA from the same mouse) cannot
#' be detected by the swap-detection pipeline, which compares a query dataset
#' against a reference. Duplicates appear as two samples that both match the
#' same reference identity — or as very-similar samples that the swap detector
#' marks \code{'match_ambiguous'} without raising an alarm. A dedicated
#' within-dataset pairwise check catches these cases.
#'
#' @param genoprobs A \code{calc_genoprob} object.
#' @param exclude_chr Chromosomes to exclude (e.g. \code{'X'} to avoid
#'   sex-chromosome effects inflating similarity for same-sex pairs).
#'   Default \code{NULL} uses all chromosomes.
#' @param threshold Minimum Pearson correlation to flag a pair as a potential
#'   duplicate. Default 0.9.
#' @param metric Similarity metric: \code{'pearson'} or \code{'cosine'}.
#'
#' @return A list with:
#'   \describe{
#'     \item{sim}{The full sample-by-sample similarity matrix (diagonal = 1).}
#'     \item{pairs}{data.frame of flagged pairs with columns:
#'       \code{sample_1}, \code{sample_2}, \code{similarity}. Sorted by
#'       descending similarity. Empty if no pairs exceed the threshold.}
#'     \item{threshold}{The threshold used.}
#'     \item{n_samples}{Number of samples.}
#'     \item{excluded_chr}{Chromosomes excluded (NULL if none).}
#'   }
#'
#' @export
genoprobs_detect_duplicates <- function(
    genoprobs,
    exclude_chr = NULL,
    threshold = 0.9,
    metric = c("pearson", "cosine")
) {
    metric <- match.arg(metric)

    # optionally drop chromosomes (e.g. sex chromosomes)
    gp <- genoprobs
    if (!is.null(exclude_chr)) {
        keep <- setdiff(names(gp), exclude_chr)
        if (length(keep) == 0) {
            stop("All chromosomes excluded; nothing left to compare.")
        }
        gp <- gp[keep]
        class(gp) <- class(genoprobs)
        attr(gp, "is_x_chr") <- attr(genoprobs, "is_x_chr")[keep]
        attr(gp, "crosstype") <- attr(genoprobs, "crosstype")
        attr(gp, "alleles") <- attr(genoprobs, "alleles")
        attr(gp, "alleleprobs") <- attr(genoprobs, "alleleprobs")
    }

    sim_result <- genoprobs_compute_similarity(gp, gp, metric = metric)
    sim <- sim_result$sim

    # use upper triangle only (avoid double-counting and self-matches)
    sim_upper <- sim
    diag(sim_upper) <- NA
    sim_upper[lower.tri(sim_upper)] <- NA

    above <- which(!is.na(sim_upper) & sim_upper >= threshold, arr.ind = TRUE)

    if (nrow(above) == 0) {
        pairs_df <- data.frame(
            sample_1   = character(0),
            sample_2   = character(0),
            similarity = numeric(0),
            stringsAsFactors = FALSE
        )
    } else {
        pairs_df <- data.frame(
            sample_1   = rownames(sim)[above[, 1]],
            sample_2   = colnames(sim)[above[, 2]],
            similarity = sim_upper[above],
            stringsAsFactors = FALSE
        )
        pairs_df <- pairs_df[order(pairs_df$similarity, decreasing = TRUE), ]
        rownames(pairs_df) <- NULL
    }

    list(
        sim          = sim_result$sim,
        pairs        = pairs_df,
        threshold    = threshold,
        n_samples    = nrow(sim),
        excluded_chr = exclude_chr
    )
}
