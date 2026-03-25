#' Read a single GBRS genoprobs file.
#'
#' Reads a tab-delimited GBRS genoprobs file and returns a matrix with
#' founder/haplotype columns labeled A--H. Used internally by
#' \code{gbrs_build_genoprobs}; not exported.
#'
#' @param path Path to a GBRS genoprobs TSV file.
#'
#' @return A matrix with rows = markers and columns = A--H (founder genotypes).
#'
#' @keywords internal
gbrs_read_genoprobs_file <- function(path) {
    # read GBRS genoprobs TSV
    df <- read.delim(path, header = TRUE, sep = "\t", check.names = FALSE)

    # keep only numeric columns (founder probs); ignore any extra non-numeric cols
    num_cols <- vapply(df, is.numeric, logical(1))
    founders <- df[, num_cols, drop = FALSE]

    # We expect exactly 8 founder columns (A–H)
    if (ncol(founders) != 8L) {
        stop(
            "Expected 8 numeric founder columns in ", path,
            " but found ", ncol(founders), "."
        )
    }

    m <- as.matrix(founders)
    colnames(m) <- LETTERS[1:8]  # A–H
    m
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
    gbrs_files_tbl <- resolve_col_samples(gbrs_files_tbl)
    markers <- resolve_col_markers(markers)

    if (!all(c("sample_id", "full_path_genoprobs") %in% names(gbrs_files_tbl))) {
        stop("gbrs_files_tbl must contain 'sample_id' and 'full_path_genoprobs'")
    }
    if (nrow(gbrs_files_tbl) == 0) {
        stop("gbrs_files_tbl has 0 rows; no genoprobs files to read.")
    }
    if (anyDuplicated(gbrs_files_tbl$sample_id)) {
        stop("gbrs_files_tbl has duplicated sample_id values.")
    }
    if (!all(c("marker_id", "chr", "pos") %in% names(markers))) {
        stop("markers must contain 'marker_id', 'chr', and 'pos'")
    }
    if (nrow(markers) == 0) {
        stop("markers has 0 rows; no markers available.")
    }
    if (anyDuplicated(markers$marker_id)) {
        stop("markers has duplicated marker_id values.")
    }

    # read in the probabilities
    raw_probs <- sapply(
        gbrs_files_tbl$full_path_genoprobs,
        gbrs_read_genoprobs_file,
        simplify = FALSE,
        USE.NAMES = TRUE
    )

    # get data dimensions
    samples   <- gbrs_files_tbl$sample_id
    n_samples <- length(samples)
    marker_ids <- markers$marker_id
    n_markers <- length(marker_ids)

    # Validate GBRS file shape before array assignment.
    # Expected per file: rows = marker count, cols = 8 founders (A-H).
    bad_shape <- names(raw_probs)[!vapply(raw_probs, function(x) {
        is.matrix(x) &&
            nrow(x) == n_markers &&
            ncol(x) == 8
    }, logical(1))]
    if (length(bad_shape) > 0) {
        stop(
            "Genoprobs files have unexpected dimensions for samples: ",
            paste(bad_shape, collapse = ", "),
            ". Expected ", n_markers, " rows (markers) and 8 columns (founders A-H)."
        )
    }

    # make an empty array (samples, haplotypes, markers) and fill it by sample
    probs <- array(
        NA,
        dim = c(n_samples, 8, n_markers),
        dimnames = list(samples, LETTERS[1:8], marker_ids)
    )

    # copy the data into the probs array
    for(i in seq_along(raw_probs)) {
        probs[i,,] <- t(raw_probs[[i]])
    }

    # convert the probs array into a genoprobs object
    qtl2convert::probs_doqtl_to_qtl2(
        probs,
        map = markers,
        marker_column = "marker_id",
        chr_column = "chr",
        pos_column = "pos"
    )
}
