#' Synchronize genoprobs with samples and/or markers data frames.
#'
#' Convenience wrapper that synchronizes genoprobs to samples and/or markers
#' data frames in a single call. Calls \code{genoprobs_sync_samples()} and/or
#' \code{genoprobs_sync_markers()} internally, applying them in sequence.
#'
#' @param genoprobs qtl2-style genoprobs (named list of 3D arrays by chromosome).
#' @param samples_df Optional data frame with sample IDs. If provided, genoprobs
#'   will be synchronized to samples in this data frame. Aliases are accepted
#'   and normalized via \code{resolve_col_samples()}.
#' @param markers_df Optional data frame with marker information. If provided,
#'   genoprobs will be synchronized to markers in this data frame. Aliases are
#'   accepted and normalized via \code{resolve_col_markers()}.
#' @param sample_order How to order samples when \code{samples_df} is provided:
#'   \code{'samples'} (default) preserves order from \code{samples_df},
#'   \code{'genoprobs'} preserves order from \code{genoprobs},
#'   \code{'alphabetical'} sorts samples alphabetically.
#' @param unit Position unit for markers: \code{'auto'} (detect from values),
#'   \code{'Mb'}, or \code{'bp'}. Only used when \code{markers_df} is provided.
#'
#' @return A list with components:
#'   \item{genoprobs}{Genoprobs synchronized to samples and/or markers as specified.}
#'   \item{samples}{If \code{samples_df} was provided: sample data frame restricted
#'     to samples present in genoprobs.}
#'   \item{markers}{If \code{markers_df} was provided: marker data frame restricted
#'     to markers present in genoprobs, sorted by chr/pos.}
#'   \item{map}{If \code{markers_df} was provided: qtl2-style marker map.}
#'   \item{common_samples}{If \code{samples_df} was provided: character vector of
#'     sample IDs present in both.}
#'   \item{dropped_from_genoprobs_samples}{If \code{samples_df} was provided:
#'     sample IDs in genoprobs but not in \code{samples_df}.}
#'   \item{dropped_from_samples_df}{If \code{samples_df} was provided: sample IDs
#'     in \code{samples_df} but not in genoprobs.}
#'   \item{dropped_markers_not_in_gp}{If \code{markers_df} was provided:
#'     per-chromosome list of marker IDs in \code{markers_df} but not in genoprobs.}
#'   \item{dropped_markers_not_in_map}{If \code{markers_df} was provided:
#'     per-chromosome list of marker IDs in genoprobs but not in \code{markers_df}.}
#'
#' @details
#' When both \code{samples_df} and \code{markers_df} are provided, samples are
#' synchronized first, then markers. This ensures the final genoprobs object
#' contains only the intersection of samples and markers specified.
#'
#' @export
genoprobs_sync <- function(genoprobs,
                          samples_df = NULL,
                          markers_df = NULL,
                          sample_order = c('samples', 'genoprobs', 'alphabetical'),
                          unit = c('auto', 'Mb', 'bp')) {
    sample_order <- match.arg(sample_order)
    unit <- match.arg(unit)

    if(is.null(samples_df) && is.null(markers_df)) {
        stop('At least one of samples_df or markers_df must be provided.')
    }

    result <- list()

    # Sync samples first if provided
    if(!is.null(samples_df)) {
        sync_samples_result <- genoprobs_sync_samples(genoprobs, samples_df, sample_order = sample_order)
        genoprobs <- sync_samples_result$genoprobs
        result$genoprobs <- genoprobs
        result$samples <- sync_samples_result$samples
        result$common_samples <- sync_samples_result$common_samples
        result$dropped_from_genoprobs_samples <- sync_samples_result$dropped_from_genoprobs
        result$dropped_from_samples_df <- sync_samples_result$dropped_from_samples_df
    } else {
        result$genoprobs <- genoprobs
    }

    # Sync markers if provided
    if(!is.null(markers_df)) {
        sync_markers_result <- genoprobs_sync_markers(genoprobs, markers_df, unit = unit)
        result$genoprobs <- sync_markers_result$genoprobs
        result$markers <- sync_markers_result$markers
        result$map <- sync_markers_result$map
        result$dropped_markers_not_in_gp <- sync_markers_result$dropped_markers_not_in_gp
        result$dropped_markers_not_in_map <- sync_markers_result$dropped_markers_not_in_map
    }

    result
}
