#' Get files in a directory matching a pattern.
#'
#' Recursively finds all files under \code{directory} whose names match
#' \code{file_pattern}. Extracts \code{sample_id} by stripping the pattern from
#' the basename (assumes sample ID is the part before the pattern). Used
#' internally by \code{gbrs_find_files}; not exported.
#'
#' @param directory Directory to search for files.
#' @param file_pattern Regular expression pattern to match file names (also
#'   used to derive sample_id from basename).
#'
#' @return A tibble with columns:
#'   \item{sample_id}{Character; derived by removing \code{file_pattern} from
#'     the file basename.}
#'   \item{full_path}{Character; full path to each file.}
#'   \item{file_name}{Character; basename of each file.}
get_files_tbl <- function(directory, file_pattern) {
    # recursively find all files
    file_list <- list.files(
        directory,
        pattern = file_pattern,
        recursive = TRUE,
        full.names = TRUE
    )

    file_list <- normalizePath(file_list, mustWork = TRUE)

    # create a tibble with the sample id, full path, file name
    # WARNING: Assuming sample id is embedded as first part of file name
    files_tbl <- tibble::tibble(
        sample_id  = basename(file_list) |>
            stringr::str_replace(file_pattern, ''),
        full_path = file_list,
        file_name = basename(file_list)
    )

    files_tbl
}

#' Find GBRS output files.
#'
#' Recursively finds GBRS genoprobs and read-counts files in a directory.
#' Ensures each sample has both file types; if any sample is missing one type,
#' reports the mismatches and stops. Returns a tibble with sample IDs and paths
#' to both file types for samples that have complete sets.
#'
#' @param directory Directory to search for GBRS output files.
#' @param genoprobs_file_pattern Regular expression to match genoprobs file
#'   names (default: \code{'.gbrs.interpolated.genoprobs.tsv'}).
#' @param counts_file_pattern Regular expression to match read-counts file
#'   names (default: \code{'.diploid.genes.expected_read_counts'}).
#'
#' @return A tibble with columns:
#'   \item{sample_id}{Character; sample identifier.}
#'   \item{full_path_genoprobs}{Character; full path to genoprobs file.}
#'   \item{file_name_genoprobs}{Character; basename of genoprobs file.}
#'   \item{full_path_counts}{Character; full path to counts file.}
#'   \item{file_name_counts}{Character; basename of counts file.}
#'
#' @export
gbrs_find_files <- function(
        directory,
        genoprobs_file_pattern = '\\.gbrs\\.interpolated\\.genoprobs\\.tsv$',
        counts_file_pattern = '\\.diploid\\.genes\\.expected_read_counts'
) {
    # get all genoprob files
    genoprobs_files_tbl <- get_files_tbl(directory, genoprobs_file_pattern)

    # get all count files
    counts_files_tbl <- get_files_tbl(directory, counts_file_pattern)

    #
    # WARNING: we are assuming that the genoprob files and counts files have
    # the same sample id
    #
    if (nrow(genoprobs_files_tbl) != nrow(counts_files_tbl)) {
        stop('ERROR: unequal number of genoprobs files and counts files')
    }

    # verify sample_ids are in both genoprobs and counts files
    tmp <- setdiff(genoprobs_files_tbl$sample_id, counts_files_tbl$sample_id)
    if (length(tmp) != 0) {
        message(
            'The following samples identifiers have a genoprobs file, but no ',
            'counts files: ', paste(tmp, collapse = ', ')
        )
    }

    # check if all sample identifiers in the sample file have a genoprobs file
    tmp2 <- setdiff(counts_files_tbl$sample_id, genoprobs_files_tbl$sample_id)
    if (length(tmp2) != 0) {
        message(
            'The following sample identifiers have a counts file, but no',
            'genoprobs file: ', paste(tmp2, collapse = ', ')
        )
    }

    if (length(tmp) + length(tmp2) != 0) {
        stop('Please validate the input directory')
    }

    # merge the sample file tibble with the genoprobs files tibble
    files_tbl <-
        genoprobs_files_tbl |>
        dplyr::inner_join(
            counts_files_tbl,
            suffix = c('_genoprobs', '_counts'),
            by = c('sample_id' = 'sample_id')
        )

    files_tbl
}
