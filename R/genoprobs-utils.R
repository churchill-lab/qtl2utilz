################################################################################
# Internal helpers
################################################################################

#' Extract genoprobs as a plain list by chromosome.
#'
#' Copies each chromosome array from a qtl2-style genoprobs object into a
#' plain list with the same names. Used to avoid subclass subsetting behavior
#' when modifying or combining genoprobs. Not exported.
#'
#' @param genoprobs R/qtl2 calc_genoprobs object.
#'
#' @return A plain list (no extra class) with names = chromosome names and
#'   elements = the corresponding 3D arrays.
extract_chr_list <- function(genoprobs) {
    chrs <- names(genoprobs)
    if(is.null(chrs) || !length(chrs)) {
        stop('genoprobs has no chromosome names; cannot extract.')
    }

    # build a plain list by extracting each chromosome with [[ ]] to avoid
    # subclass subsetting (e.g. [) that can alter structure or class
    out <- setNames(vector('list', length(chrs)), chrs)
    for(chr in chrs) {
        out[[chr]] <- genoprobs[[chr]]
    }

    # make sure result is a plain list with no extra class
    class(out) <- "list"
    out
}

#' Detect position units from a numeric vector.
#'
#' If the maximum finite position is less than 2000, returns 'Mb';
#' otherwise returns 'bp'. Used internally for 'unit = 'auto''.
#' Not exported.
#'
#' @param pos Numeric vector of genomic positions.
#'
#' @return Character string \code{'Mb'} or \code{'bp'}.
detect_position_units_vec <- function(pos) {
    pos <- as.numeric(pos)
    pos <- pos[is.finite(pos)]
    if(!length(pos)) stop('No finite positions found in vector')
    # If max position < 2000, assume Mb (mouse chr in Mb); else assume bp.
    if(max(pos, na.rm = TRUE) < 2000) 'Mb' else 'bp'
}

#' Detect position units from a data frame column.
#'
#' Convenience wrapper: extracts \code{pos_col} and calls
#' \code{detect_position_units_vec}. Not exported.
#'
#' @param data Data frame with a position column.
#' @param pos_col Name of the position column.
#'
#' @return Character string \code{'Mb'} or \code{'bp'}.
detect_position_units_df <- function(data, pos_col = 'pos') {
    if(!pos_col %in% names(data)) {
        stop('Column not found: ', pos_col)
    }
    
    detect_position_units_vec(data[[pos_col]])
}

#' Convert a position vector to base pairs.
#'
#' Single source of truth for converting a numeric vector to bp. With
#' \code{unit = 'auto'}, calls \code{detect_position_units_vec} to decide.
#' Not exported.
#'
#' @param pos Numeric vector of positions.
#' @param unit \code{'auto'}, \code{'Mb'}, or \code{'bp'}.
#'
#' @return Numeric vector of positions in bp.
positions_vec_to_bp <- function(pos,
                                unit = c('auto', 'Mb', 'bp')) {
    unit <- match.arg(unit)

    pos <- as.numeric(pos)

    # auto: infer Mb vs bp from magnitude (max < 2000 => Mb)
    if(unit == 'auto') {
        unit <- detect_position_units_vec(pos)
    }

    if(unit == 'Mb') {
        pos <- pos * 1e6
    }
    # 'bp' leaves pos unchanged

    pos
}
