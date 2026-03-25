
.qtl2utilz_aliases <- list(
    sample_id = c("sample_id", "sample.id", "sampleID", "sample", "mouse_id", "mouse.id", "mouseID", "mouse"),
    marker_id = c("marker_id", "marker.id", "markerID", "markerid", "marker"),
    chr       = c("chr", "chrom", "chromosome", "chrom_id", "chromosome_id"),
    pos       = c("pos", "position")
)

# Normalize column tokens for flexible matching.
# Lowercase + collapse punctuation/whitespace to underscores.
.normalize_col_token <- function(x) {
    x <- tolower(as.character(x))
    x <- gsub("[^a-z0-9]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    x
}

# Resolve one canonical key to an existing column name.
# Internal utility used by the wrapper helpers below.
.resolve_col <- function(df, key, aliases = .qtl2utilz_aliases) {
    if (!is.data.frame(df)) {
        stop("Expected a data.frame for column resolution.")
    }
    if (!key %in% names(aliases)) {
        stop("Unknown alias key: ", key)
    }

    candidates <- unique(aliases[[key]])
    norm_candidates <- unique(.normalize_col_token(candidates))
    norm_names <- .normalize_col_token(names(df))
    hit_idx <- which(norm_names %in% norm_candidates)
    hits <- names(df)[hit_idx]

    if (length(hits) == 1) {
        hits
    }
    if (length(hits) == 0) {
        stop(
            "Missing column for '", key, "'. Accepted aliases: ",
            paste(candidates, collapse = ", "),
            " (case-insensitive; punctuation-insensitive)."
        )
    }

    stop(
        "Ambiguous columns for '", key, "': ",
        paste(hits, collapse = ", "),
        ". Keep only one of these aliases."
    )
}

#' Resolve and normalize sample key column to canonical sample_id
#'
#' Renames the detected sample column to \code{sample_id}. Common column name
#' variants (e.g. \code{sample}, \code{mouse_id}) are recognized case- and
#' punctuation-insensitively.
#'
#' @param df A data frame with a sample identifier column (or an accepted alias).
#'
#' @return \code{df} with the sample column renamed to \code{sample_id}.
#'
#' @export
resolve_col_samples <- function(df) {
    sample_col <- .resolve_col(df, "sample_id")
    names(df)[names(df) == sample_col] <- "sample_id"
    df
}

#' Resolve and normalize marker/map key columns to canonical names
#'
#' Renames columns to \code{marker_id}, \code{chr}, and \code{pos}. Common
#' name variants are recognized case- and punctuation-insensitively.
#'
#' @param df A data frame with marker ID, chromosome, and position columns (or accepted aliases).
#'
#' @return \code{df} with those columns renamed to \code{marker_id}, \code{chr}, and \code{pos}.
#'
#' @export
resolve_col_markers <- function(df) {
    marker_col <- .resolve_col(df, "marker_id")
    chr_col <- .resolve_col(df, "chr")
    pos_col <- .resolve_col(df, "pos")

    names(df)[names(df) == marker_col] <- "marker_id"
    names(df)[names(df) == chr_col] <- "chr"
    names(df)[names(df) == pos_col] <- "pos"
    df
}

#' Standardize column names with janitor
#'
#' Thin wrapper around \code{janitor::clean_names()}.
#'
#' @param df A data frame.
#'
#' @return \code{df} with syntactically consistent, lowercase snake_case names.
#'
#' @export
fix_cols <- function(df) {
    janitor::clean_names(df)
}

