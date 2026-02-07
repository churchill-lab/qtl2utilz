#' Align samples across two genoprobs objects.
#'
#' Subset both genoprobs objects to their common samples (and optionally
#' sort sample order). Useful before \code{genoprobs_correlate} or other
#' comparisons so that both objects have the same samples in the same order.
#'
#' @param genoprobs_1 First qtl2-style genoprobs object.
#' @param genoprobs_2 Second qtl2-style genoprobs object.
#' @param sort If \code{TRUE}, order common samples alphabetically.
#'
#' @return A list with components:
#'   \item{genoprobs_1}{First genoprobs restricted to common chromosomes and common samples.}
#'   \item{genoprobs_2}{Second genoprobs restricted to common chromosomes and common samples.}
#'   \item{common_samples}{Character vector of sample IDs in both.}
#'   \item{dropped_from_1}{Sample IDs only in \code{genoprobs_1}.}
#'   \item{dropped_from_2}{Sample IDs only in \code{genoprobs_2}.}
#'
#' @details
#' Only overlapping chromosomes are considered. A sample is included only if it
#' appears on every common chromosome in both objects. Attributes and class of
#' each genoprobs are preserved in the returned subsets.
#'
#' @export
genoprobs_sync_samples <- function(genoprobs_1, genoprobs_2, sort = TRUE) {
    # use only chromosomes present in both objects
    common_chr <- intersect(names(genoprobs_1), names(genoprobs_2))
    if(!length(common_chr)) {
        stop('No overlapping chromosomes between the two genoprobs objects.')
    }

    # samples that appear on every common chromosome in each object
    # use [[chr]] only; avoid [common_chr] which triggers qtl2 subset() and can error
    s1_per_chr <- setNames(lapply(common_chr, function(chr) rownames(genoprobs_1[[chr]])), common_chr)
    s2_per_chr <- setNames(lapply(common_chr, function(chr) rownames(genoprobs_2[[chr]])), common_chr)
    if(any(vapply(s1_per_chr, is.null, NA)) || any(vapply(s2_per_chr, is.null, NA))) {
        stop('One or both genoprobs objects have NULL sample names on at least one chromosome.')
    }

    common_in_1 <- Reduce(intersect, s1_per_chr)
    common_in_2 <- Reduce(intersect, s2_per_chr)
    common <- intersect(common_in_1, common_in_2)
    if(!length(common)) {
        stop('No overlapping samples between the two genoprobs objects (after requiring samples to appear on all chromosomes).')
    }

    if(sort) common <- sort(common)

    # record which samples were dropped from each object for reporting
    all_in_1 <- Reduce(union, s1_per_chr)
    all_in_2 <- Reduce(union, s2_per_chr)
    dropped1 <- setdiff(all_in_1, common)
    dropped2 <- setdiff(all_in_2, common)

    # subset a genoprobs object to given samples on common chromosomes only
    # build result with new list and [[chr]] only; avoid gp[common_chr] which triggers qtl2 subset()
    subset_to <- function(gp, samples) {
        out <- setNames(vector('list', length(common_chr)), common_chr)
        for(chr in common_chr) {
            pr <- gp[[chr]]
            i <- match(samples, rownames(pr))
            if(anyNA(i)) {
                stop('Sample missing on chromosome ', chr, '; ensure samples appear on all chromosomes.')
            }
            out[[chr]] <- pr[i, , , drop = FALSE]
        }
        attributes(out) <- attributes(gp)
        class(out) <- class(gp)
        out
    }

    list(
        genoprobs_1      = subset_to(genoprobs_1, common),
        genoprobs_2      = subset_to(genoprobs_2, common),
        common_samples   = common,
        dropped_from_1   = dropped1,
        dropped_from_2   = dropped2
    )
}
