#' Intersect samples across two genoprobs objects.
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
genoprobs_intersect_samples <- function(genoprobs_1, genoprobs_2, sort = TRUE) {
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

#' Synchronize genoprobs to a samples data frame.
#'
#' Subset and reorder genoprobs so that it contains only samples that appear
#' in \code{samples_df}, optionally preserving the order from either the
#' genoprobs object or the samples data frame. Useful to align genoprobs with
#' a standard sample set (e.g. for filtering to specific samples or matching
#' with phenotype data).
#'
#' @param genoprobs qtl2-style genoprobs (named list of 3D arrays by chromosome).
#' @param samples_df Data frame with at least a sample ID column. Aliases are
#'   accepted and normalized via \code{resolve_col_samples()}.
#' @param sample_order How to order samples in the output:
#'   \code{'samples'} (default) preserves order from \code{samples_df},
#'   \code{'genoprobs'} preserves order from \code{genoprobs},
#'   \code{'alphabetical'} sorts samples alphabetically.
#'
#' @return A list with components:
#'   \item{genoprobs}{Genoprobs restricted to common chromosomes and
#'     samples present in both genoprobs and \code{samples_df}, in the
#'     order specified by \code{sample_order}.}
#'   \item{samples}{Sample data frame restricted to samples that appear in
#'     genoprobs; only samples present in both are included.}
#'   \item{common_samples}{Character vector of sample IDs present in both.}
#'   \item{dropped_from_genoprobs}{Sample IDs in genoprobs but not in \code{samples_df}.}
#'   \item{dropped_from_samples_df}{Sample IDs in \code{samples_df} but not in genoprobs.}
#'
#' @details
#' Only overlapping chromosomes are considered. A sample is included only if it
#' appears on every common chromosome in genoprobs. Attributes and class of
#' genoprobs are preserved in the returned subset.
#'
#' @export
genoprobs_sync_samples <- function(genoprobs,
                                   samples_df,
                                   sample_order = c('samples', 'genoprobs', 'alphabetical')) {
    sample_order <- match.arg(sample_order)
    samples_df <- resolve_col_samples(samples_df)

    if(!'sample_id' %in% names(samples_df)) {
        stop('samples_df must contain a sample_id column (or an accepted alias).')
    }

    # extract a plain list of chr arrays so subsetting [common_chr] is safe
    genoprobs_list <- extract_chr_list(genoprobs)

    # get sample IDs from samples_df
    samples_df_ids <- unique(samples_df$sample_id)

    # get sample IDs from genoprobs (must appear on all chromosomes)
    # use [[chr]] only; avoid [common_chr] which triggers qtl2 subset()
    common_chr <- names(genoprobs_list)
    if(!length(common_chr)) {
        stop('genoprobs has no chromosomes.')
    }

    gp_samples_per_chr <- setNames(
        lapply(common_chr, function(chr) rownames(genoprobs_list[[chr]])),
        common_chr
    )
    if(any(vapply(gp_samples_per_chr, is.null, NA))) {
        stop('genoprobs has NULL sample names on at least one chromosome.')
    }

    # samples that appear on every chromosome in genoprobs
    gp_samples <- Reduce(intersect, gp_samples_per_chr)

    # find intersection
    common_samples <- intersect(gp_samples, samples_df_ids)
    if(!length(common_samples)) {
        stop('No overlapping samples between genoprobs and samples_df.')
    }

    # determine final order based on sample_order parameter
    if(sample_order == 'samples') {
        # preserve order from samples_df (keep samples_df order, filter to common)
        common_samples <- samples_df_ids[samples_df_ids %in% common_samples]
    } else if(sample_order == 'genoprobs') {
        # preserve order from genoprobs
        common_samples <- gp_samples[gp_samples %in% common_samples]
    } else { # sample_order == 'alphabetical'
        common_samples <- sort(common_samples)
    }

    # record which samples were dropped
    dropped_from_genoprobs <- setdiff(gp_samples, samples_df_ids)
    dropped_from_samples_df <- setdiff(samples_df_ids, gp_samples)

    # subset genoprobs to common samples
    out <- setNames(vector('list', length(common_chr)), common_chr)
    for(chr in common_chr) {
        pr <- genoprobs_list[[chr]]
        if(length(dim(pr)) != 3) {
            stop(sprintf('Chr %s genoprobs is not a 3D array.', chr))
        }

        gp_chr_samples <- rownames(pr)
        if(is.null(gp_chr_samples)) {
            stop(sprintf('Chr %s genoprobs has NULL sample names.', chr))
        }

        # match to get indices in genoprobs order, then reorder to common_samples order
        idx <- match(common_samples, gp_chr_samples)
        if(anyNA(idx)) {
            stop(sprintf('Sample missing on chromosome %s; ensure samples appear on all chromosomes.', chr))
        }
        out[[chr]] <- pr[idx, , , drop = FALSE]
    }

    # restore original attributes/class so result behaves like input genoprobs
    attributes(out) <- attributes(genoprobs)
    class(out) <- class(genoprobs)

    # subset samples_df to only samples present in genoprobs
    samples_out <- samples_df[samples_df$sample_id %in% common_samples, , drop = FALSE]
    # remove duplicates if any, keeping first occurrence (preserves order)
    samples_out <- samples_out[!duplicated(samples_out$sample_id), , drop = FALSE]
    # reorder to match common_samples order
    samples_out <- samples_out[match(common_samples, samples_out$sample_id), , drop = FALSE]

    list(
        genoprobs                = out,
        samples                  = samples_out,
        common_samples           = common_samples,
        dropped_from_genoprobs   = dropped_from_genoprobs,
        dropped_from_samples_df  = dropped_from_samples_df
    )
}
