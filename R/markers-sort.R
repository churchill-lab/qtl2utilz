#' Rank chromosomes for ordering.
#'
#' Returns a numeric rank for each chromosome name so that standard autosomes
#' (1–19) come first in numeric order, then X, Y, XY, MT, then any other
#' contigs in a stable order. Used by \code{markers_sort}. Not exported.
#'
#' @param chr Character (or coercible) vector of chromosome names (e.g.
#'   \code{'1'}, \code{'X'}, \code{'chr2'}).
#'
#' @return Numeric vector of ranks with the same length as \code{chr}.
chr_rank <- function(chr) {
    x <- toupper(as.character(chr))
    # drop leading 'chr' for consistent comparison
    x <- gsub('^CHR', '', x)

    # normalize mitochondria naming so M/MTDNA/MITO all sort together
    x[x %in% c('M', 'MTDNA', 'MITO')] <- 'MT'

    # numeric chromosomes (1-19) get their numeric value as rank
    is_num <- grepl('^[0-9]+$', x)
    r <- rep(NA_real_, length(x))
    r[is_num] <- as.numeric(x[is_num])

    # non-numeric chromosomes: fixed order after autosomes (X, Y, XY, MT)
    r[x == 'X']  <- 1e6 + 1
    r[x == 'Y']  <- 1e6 + 2
    r[x == 'XY'] <- 1e6 + 3
    r[x == 'MT'] <- 1e6 + 4

    # any other contigs/scaffolds: rank after MT, stable order by alphabet
    other <- is.na(r)
    if(any(other)) {
        ord_other <- order(x[other])
        r[other][ord_other] <- (1e6 + 100) + seq_len(sum(other))
    }

    r
}

#' Sort a data frame of genetic markers by chromosome and position
#'
#' Sorts rows by chromosome (using standard order 1–19, X, Y, XY, MT, then
#' others), then by position, then by marker ID for deterministic tie-breaking.
#' Marker/map key aliases are accepted and normalized via
#' \code{resolve_col_markers()}.
#'
#' @param markers Data frame containing marker information.
#'
#' @return The data frame \code{markers} with rows reordered. No columns are added
#'   or removed.
#'
#' @export
markers_sort <- function(markers) {
    markers <- resolve_col_markers(markers)

    if(!all(c('marker_id', 'chr', 'pos') %in% names(markers))) {
        stop('markers must contain "marker_id", "chr", and "pos"')
    }

    # order by chromosome rank, then position, then marker name (ties deterministic)
    o <- order(chr_rank(markers$chr),
               markers$pos,
               markers$marker_id)
    markers[o, , drop = FALSE]
}
