#' Convert genomic positions to base pairs
#'
#' Generic function to convert genomic positions to base pairs (bp). Supports
#' numeric vectors and data frames (updating a position column in place).
#' With \code{unit = "auto"}, units are inferred: if the maximum position is
#' less than 2000, values are treated as Mb and multiplied by 1e6; otherwise
#' treated as bp and returned unchanged.
#'
#' @param x Object containing genomic positions (numeric vector or data frame).
#' @param ... Arguments passed to methods (e.g. \code{pos_col}, \code{unit}).
#'
#' @return For \code{numeric}: a numeric vector in bp. For \code{data.frame}:
#'   \code{x} with the position column replaced by bp values.
#'
#' @seealso \code{\link{genoprobs_sync_markers}} (uses this for marker maps)
#'
#' @export
positions_to_bp <- function(x, ...) {
    UseMethod("positions_to_bp")
}

#' @rdname positions_to_bp
#' @param unit \code{"auto"} (infer from values), \code{"Mb"}, or \code{"bp"}.
#' @export
positions_to_bp.default <- function(x, ...) {
    # Only numeric and data.frame have methods; anything else is unsupported.
    stop(
        "Unsupported type for positions_to_bp(): ",
        paste(class(x), collapse = "/")
    )
}

#' @rdname positions_to_bp
#' @export
positions_to_bp.numeric <- function(x,
                                    unit = c("auto", "Mb", "bp"),
                                    ...) {
    unit <- match.arg(unit)
    # Delegate to single source of truth for vector conversion.
    positions_vec_to_bp(x, unit = unit)
}

#' @rdname positions_to_bp
#' @param pos_col Name of the column in \code{x} containing positions (for
#'   data frame method).
#' @export
positions_to_bp.data.frame <- function(x,
                                       pos_col = "pos",
                                       unit = c("auto", "Mb", "bp"),
                                       ...) {
    unit <- match.arg(unit)

    if(!pos_col %in% names(x)) {
        stop("Column not found: ", pos_col)
    }

    # Convert position column in place; leave other columns unchanged.
    x[[pos_col]] <- positions_vec_to_bp(
        x[[pos_col]],
        unit = unit
    )

    x
}
