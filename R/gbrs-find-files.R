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
gbrs_get_files_tbl <- function(directory, file_pattern) {
    # Recursively find all matching files under `directory`.
    file_list <- list.files(
        directory,
        pattern = file_pattern,
        recursive = TRUE,
        full.names = TRUE
    )

    file_list <- normalizePath(file_list, mustWork = TRUE)

    # Build a tibble with sample ID, full path, and base file name.
    # This assumes the sample ID is the leading part of the file name.
    files_tbl <- tibble::tibble(
        sample_id = basename(file_list) |>
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
        counts_file_pattern = '\\.diploid\\.genes\\.expected_read_counts$'
) {
    # Get all genoprobs files.
    genoprobs_files_tbl <- gbrs_get_files_tbl(directory, genoprobs_file_pattern)

    # Get all count files.
    counts_files_tbl <- gbrs_get_files_tbl(directory, counts_file_pattern)

    # Identify samples that are missing one of the two GBRS file types.
    missing_counts_samples <- setdiff(
        genoprobs_files_tbl$sample_id,
        counts_files_tbl$sample_id
    )
    if (length(missing_counts_samples) != 0) {
        message(
            'The following sample identifiers have a genoprobs file, but no ',
            'counts file: ', paste(missing_counts_samples, collapse = ', ')
        )
    }

    missing_genoprobs_samples <- setdiff(
        counts_files_tbl$sample_id,
        genoprobs_files_tbl$sample_id
    )
    if (length(missing_genoprobs_samples) != 0) {
        message(
            'The following sample identifiers have a counts file, but no ',
            'genoprobs file: ', paste(missing_genoprobs_samples, collapse = ', ')
        )
    }

    if (nrow(genoprobs_files_tbl) == 0 || nrow(counts_files_tbl) == 0) {
        stop(
            'No matching GBRS files found. ',
            'Found ', nrow(genoprobs_files_tbl), ' genoprobs files and ',
            nrow(counts_files_tbl), ' counts files. ',
            'Please verify directory and filename patterns.'
        )
    }

    if (length(missing_counts_samples) + length(missing_genoprobs_samples) != 0) {
        stop('Please validate the input directory')
    }

    # Join the two file tables once the sample sets match.
    files_tbl <-
        genoprobs_files_tbl |>
        dplyr::inner_join(
            counts_files_tbl,
            suffix = c('_genoprobs', '_counts'),
            by = c('sample_id' = 'sample_id')
        )

    files_tbl
}
