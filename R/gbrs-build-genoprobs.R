#' Read a single GBRS genoprobs file.
#'
#' Reads a tab-delimited GBRS genoprobs file and returns a matrix with
#' founder/haplotype columns labeled A--H. Used internally by
#' \code{gbrs_build_genoprobs}; not exported.
#'
#' @param x Path to a GBRS genoprobs TSV file.
#'
#' @return A matrix with rows = markers and columns = A--H (founder genotypes).
gbrs_read_genoprobs_file <- function (x) {
    x <- read.delim(x, header=TRUE, sep = '\t')
    colnames(x) <- LETTERS[1:8]
    x
}


#' Build qtl2 genoprobs from GBRS output files.
#'
#' Reads GBRS genoprobs TSV files (one per sample), combines them into a
#' 3D array (samples x founders x markers), and converts to R/qtl2 format
#' via \code{qtl2convert::probs_doqtl_to_qtl2}. Expects \code{gbrs_files_tbl}
#' from \code{gbrs_find_files} and a marker map with chr, pos, and marker IDs.
#'
#' @param gbrs_files_tbl Tibble from \code{gbrs_find_files} with columns
#'   \code{sample_id}, \code{full_path_genoprobs}, \code{file_name_genoprobs}.
#' @param markers Data frame with at least chromosome, position, and marker ID
#'   columns. Marker ID column may be named \code{marker}, \code{marker_id},
#'   \code{markerid}, \code{marker.id}, \code{marker_name}, or \code{markername}
#'   (first match wins). Chromosome and position columns must be \code{chr}
#'   and \code{pos}.
#'
#' @return A qtl2-style genoprobs object (as from \code{qtl2::calc_genoprob}).
#'
#' @export
gbrs_build_genoprobs <- function(gbrs_files_tbl, markers) {
    # read in the probabilities
    raw_probs = sapply(
        gbrs_files_tbl$full_path_genoprobs,
        gbrs_read_genoprobs_file,
        simplify = FALSE,
        USE.NAMES = TRUE
    )

    # get data dimensions
    samples   <- gbrs_files_tbl$sample_id
    n_samples <- length(samples)

    # detect marker column: accept marker, marker_id, markerid, marker.id, etc.
    marker_col_candidates <- c('marker', 'marker_id', 'markerid', 'marker.id', 'marker_name', 'markername')

    marker_col <- marker_col_candidates[marker_col_candidates %in% names(markers)][1]
    if (is.na(marker_col)) {
        stop('markers must have a column for marker IDs. Expected one of: ',
             paste(marker_col_candidates, collapse = ', '))
    }
    marker_ids <- markers[[marker_col]]
    n_markers <- length(marker_ids)

    # make an empty array (samples, haplotypes, markers) and fill it by sample
    probs <- array(
        NA,
        dim = c(n_samples, 8, n_markers),
        dimnames = list(samples, LETTERS[1:8], marker_ids)
    )

    # copy the data into the probs array
    for(i in seq_along(raw_probs)) {
        probs[i,,] = t(raw_probs[[i]])
    }

    # convert the probs array into a genoprobs object
    qtl2convert::probs_doqtl_to_qtl2(
        probs, 
        map = markers,
        marker_column = marker_col
    )
}
