#' Synchronize a samples-by-features matrix to a samples data frame.
#'
#' Subset and reorder a matrix (rows = samples, columns = features) so that it
#' contains only samples that appear in \code{samples_df}, optionally preserving
#' the order from either the matrix or the samples data frame.  This is the
#' matrix analogue of \code{\link{genoprobs_sync_samples}} and is useful for
#' aligning expression counts, protein abundances, or any other samples-by-
#' features data with a canonical sample annotation table.
#'
#' @param mat A numeric matrix (or object coercible to one) with \code{rownames}
#'   set to sample IDs.
#' @param samples_df Data frame with at least a sample ID column.  Aliases are
#'   accepted and normalized via \code{resolve_col_samples()}.
#' @param sample_order How to order samples in the output:
#'   \code{'samples'} (default) preserves order from \code{samples_df},
#'   \code{'matrix'} preserves order from \code{mat},
#'   \code{'alphabetical'} sorts samples alphabetically.
#'
#' @return A list with components:
#'   \item{matrix}{Matrix restricted to samples present in both \code{mat} and
#'     \code{samples_df}, in the order specified by \code{sample_order}.}
#'   \item{samples}{Sample data frame restricted to samples present in
#'     \code{mat}; only samples present in both are included.}
#'   \item{common_samples}{Character vector of sample IDs present in both.}
#'   \item{dropped_from_matrix}{Sample IDs in \code{mat} but not in
#'     \code{samples_df}.}
#'   \item{dropped_from_samples_df}{Sample IDs in \code{samples_df} but not in
#'     \code{mat}.}
#'
#' @export
matrix_sync_samples <- function(mat,
                                samples_df,
                                sample_order = c('samples', 'matrix', 'alphabetical')) {
    sample_order <- match.arg(sample_order)

    if (!is.matrix(mat)) mat <- as.matrix(mat)
    if (is.null(rownames(mat))) {
        stop('mat must have rownames (sample IDs).')
    }

    samples_df <- resolve_col_samples(samples_df)
    if (!'sample_id' %in% names(samples_df)) {
        stop('samples_df must contain a sample_id column (or an accepted alias).')
    }

    mat_ids <- rownames(mat)
    df_ids  <- unique(samples_df$sample_id)

    common_samples <- intersect(mat_ids, df_ids)
    if (!length(common_samples)) {
        stop('No overlapping samples between mat and samples_df.')
    }

    if (sample_order == 'samples') {
        common_samples <- df_ids[df_ids %in% common_samples]
    } else if (sample_order == 'matrix') {
        common_samples <- mat_ids[mat_ids %in% common_samples]
    } else {
        common_samples <- sort(common_samples)
    }

    dropped_from_matrix     <- setdiff(mat_ids, df_ids)
    dropped_from_samples_df <- setdiff(df_ids, mat_ids)

    mat_out <- mat[common_samples, , drop = FALSE]

    samples_out <- samples_df[samples_df$sample_id %in% common_samples, , drop = FALSE]
    samples_out <- samples_out[!duplicated(samples_out$sample_id), , drop = FALSE]
    samples_out <- samples_out[match(common_samples, samples_out$sample_id), , drop = FALSE]

    list(
        matrix                  = mat_out,
        samples                 = samples_out,
        common_samples          = common_samples,
        dropped_from_matrix     = dropped_from_matrix,
        dropped_from_samples_df = dropped_from_samples_df
    )
}
