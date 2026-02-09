#' Rank-based inverse normal transformation.
#'
#' Converts a numeric vector to approximate normality using ranks. Ranks are
#' scaled to (0, 1) and passed through \code{\link[stats]{qnorm}}. Useful for
#' transforming skewed expression or phenotype data before QTL mapping.
#'
#' @param x Numeric vector; \code{NA}s are preserved.
#'
#' @return Numeric vector of the same length as \code{x}; \code{NA} where
#'   \code{x} was \code{NA}, otherwise normal quantiles.
rankZ <- function(x) {
    x <- rank(x, na.last = 'keep', ties.method = 'average') / (sum(!is.na(x)) + 1)
    qnorm(x)
}


#' Read a single GBRS read-counts file.
#'
#' Reads a tab-delimited GBRS expected read-counts file. Expected columns:
#' \code{locus} (gene ID), founder columns A--H, \code{total}, and optionally
#' \code{notes}. Used internally by \code{gbrs_build_counts}; not exported.
#'
#' @param x Path to a GBRS read-counts TSV file.
#'
#' @return A data frame with at least \code{locus} and \code{total} columns.
read_counts_file <- function (x) {
    read.table(x, header=TRUE, sep = '\t')
}


#' Build expression matrix from GBRS read-counts files.
#'
#' Reads GBRS expected read-counts TSV files (one per sample), verifies that
#' all samples have identical locus IDs in the same order, then combines
#' \code{total} counts into a matrix with samples as rows and genes (loci)
#' as columns. Expects \code{gbrs_files_tbl} from \code{gbrs_find_files}.
#'
#' @param gbrs_files_tbl Tibble from \code{gbrs_find_files} with columns
#'   \code{sample_id} and \code{full_path_counts}.
#'
#' @return A matrix with rows = samples (rownames = sample IDs) and columns =
#'   genes (locus IDs). Values are total expected read counts.
#'
#' @export
gbrs_build_counts <- function(gbrs_files_tbl) {
    # read in the count files (one per sample)
    raw_counts <- sapply(
        gbrs_files_tbl$full_path_counts,
        read_counts_file,
        simplify = FALSE,
        USE.NAMES = TRUE
    )

    names(raw_counts) <- gbrs_files_tbl$sample_id

    # Sanity check: ensure all samples have the same locus IDs in the same order.
    locus_ref <- raw_counts[[1]]$locus
    locus_ok <- vapply(raw_counts[-1], function(x) identical(x$locus, locus_ref), logical(1))
    if (!all(locus_ok)) {
        bad <- names(which(!locus_ok))
        stop('All samples must have identical locus IDs in the same order. Mismatch in: ',
             paste(bad, collapse = ', '))
    }

    # extract locus and total columns, bind rows, pivot to samples x genes
    all_counts <- raw_counts
    for (i in 1:length(raw_counts)) {
        all_counts[[i]] <- dplyr::select(raw_counts[[i]], locus, total)
    }

    all_counts <- dplyr::bind_rows(all_counts, .id = 'sample_id')

    all_counts <- all_counts |>
        tidyr::pivot_wider(
            names_from  = locus,
            values_from = total,
            values_fill = NA
        ) |>
        tibble::column_to_rownames('sample_id')

    all_counts
}




