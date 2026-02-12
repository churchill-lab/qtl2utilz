#' Correlate genoprobs between two objects.
#'
#' Compute per-sample correlation between two qtl2-style genoprobs objects
#' (e.g. GBRS vs MUGA) over overlapping chromosomes and samples. For each
#' sample, all founder probabilities across all chromosomes are concatenated
#' and correlated. Use \code{genoprobs_sync_markers} and \code{genoprobs_sync_samples}
#' first so both objects share the same markers and samples if needed.
#'
#' @param genoprobs_1 First R/qtl2 calc_genoprobs object.
#' @param genoprobs_2 Second R/qtl2 calc_genoprobs object (same structure).
#' @param threshold Numeric; correlations below this value are flagged in
#'   \code{flag_low} (default 0.80).
#' @param method Correlation method passed to \code{\link[stats]{cor}}
#'   (e.g. \code{'pearson'}, \code{'spearman'}).
#'
#' @return A data frame with columns:
#'   \item{sample_id}{Sample ID.}
#'   \item{correlation}{Correlation between the two genoprobs for that sample.}
#'   \item{flag_low}{\code{TRUE} when correlation < \code{threshold}.}
#'
#' @export
genoprobs_correlate <- function(genoprobs_1,
                                genoprobs_2,
                                threshold = 0.80,
                                method = 'pearson') {

    # use only chromosomes present in both objects
    chrs <- intersect(names(genoprobs_1), names(genoprobs_2))
    if(!length(chrs)) {
        stop('No overlapping chromosomes.')
    }

    # sample IDs from first shared chromosome (same across chrs in qtl2 objects)
    s1 <- rownames(genoprobs_1[[chrs[1]]])
    s2 <- rownames(genoprobs_2[[chrs[1]]])
    common_samples <- intersect(s1, s2)
    if(!length(common_samples)) {
        stop('No overlapping samples.')
    }

    # one row per common sample; correlation filled in below
    results <- data.frame(
        sample_id = common_samples,
        correlation = NA_real_,
        flag_low = FALSE,
        stringsAsFactors = FALSE
    )

    # for each sample, concatenate all founder probs across chromosomes and correlate
    # Use [[chr]] only; do not use [chrs] (triggers qtl2 subset()).
    for(i in seq_along(common_samples)) {
        s <- common_samples[i]
        v1 <- c()
        v2 <- c()
        for(chr in chrs) {
            v1 <- c(v1, as.vector(genoprobs_1[[chr]][s, , ]))
            v2 <- c(v2, as.vector(genoprobs_2[[chr]][s, , ]))
        }
        results$correlation[i] <- cor(v1, v2, use = 'pairwise.complete.obs', method = method)
    }

    # flag samples whose correlation falls below threshold (e.g. QC)
    results$flag_low <- results$correlation < threshold
    results
}
