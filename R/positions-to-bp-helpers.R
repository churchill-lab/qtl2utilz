#' Detect position units from a numeric vector
#'
#' If the maximum finite position is less than 2000, returns \code{"Mb"};
#' otherwise returns \code{"bp"}. Used internally for \code{unit = "auto"}.
#' Not exported.
#'
#' @param pos Numeric vector of genomic positions.
#'
#' @return Character string \code{"Mb"} or \code{"bp"}.
#'
#' @keywords internal
detect_position_units_vec <- function(pos) {
    pos <- as.numeric(pos)
    pos <- pos[is.finite(pos)]
    if(!length(pos)) stop('No finite positions found in vector')
    # If max position < 2000, assume Mb (mouse chr in Mb); else assume bp.
    if(max(pos, na.rm = TRUE) < 2000) 'Mb' else 'bp'
}

#' Detect position units from a data frame column
#'
#' Convenience wrapper: extracts \code{pos_col} and calls
#' \code{detect_position_units_vec}. Not exported.
#'
#' @param data Data frame with a position column.
#' @param pos_col Name of the position column.
#'
#' @return Character string \code{"Mb"} or \code{"bp"}.
#'
#' @keywords internal
detect_position_units_df <- function(data, pos_col = "pos") {
    if(!pos_col %in% names(data)) {
        stop("Column not found: ", pos_col)
    }
    # Use same heuristic as for a vector: infer from position values.
    detect_position_units_vec(data[[pos_col]])
}

#' Convert a position vector to base pairs
#'
#' Single source of truth for converting a numeric vector to bp. With
#' \code{unit = "auto"}, calls \code{detect_position_units_vec} to decide.
#' Not exported.
#'
#' @param pos Numeric vector of positions.
#' @param unit \code{"auto"}, \code{"Mb"}, or \code{"bp"}.
#'
#' @return Numeric vector of positions in bp.
#'
#' @keywords internal
positions_vec_to_bp <- function(pos,
                                unit = c("auto", "Mb", "bp")) {
    unit <- match.arg(unit)

    pos <- as.numeric(pos)

    # Auto: infer Mb vs bp from magnitude (max < 2000 => Mb).
    if(unit == "auto") {
        unit <- detect_position_units_vec(pos)
    }

    if(unit == "Mb") {
        pos <- pos * 1e6
    }
    # "bp" leaves pos unchanged

    pos
}
