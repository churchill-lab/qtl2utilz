#' Arrange rows in human/natural sort order across one or more columns
#'
#' Sorts a data frame using "human" ordering rather than strict
#' lexicographic ordering. For example, values like `sample2`, `sample10`,
#' and `sample100` sort in natural order rather than plain character order.
#'
#' This function accepts any number of columns and supports either:
#' - one direction applied to all columns, or
#' - one direction per column
#'
#' @param df A data frame or tibble.
#' @param ... Unquoted column names to sort by, in priority order.
#'   The first column is the primary sort key, the second breaks ties
#'   within the first, and so on.
#' @param .dir Sort direction. Either:
#'   - a single string, `"asc"` or `"desc"`, applied to all columns, or
#'   - a character vector with one value per column, each being
#'     `"asc"` or `"desc"`.
#'
#' @return A data frame with rows reordered according to the requested
#'   human-sort order.
#'
#' @details
#' The function works by repeatedly sorting the data frame from the
#' last requested column to the first requested column. This lets earlier
#' columns remain the highest-priority sort keys while still allowing
#' natural/human ordering within each column.
#'
#' Each selected column is converted to character before calling
#' [gtools::mixedorder()], so values like `M2`, `M10`, and `M100`
#' sort as a human would expect.
#'
#' @examples
#' df <- tibble::tibble(
#'   sex = c("F", "M", "F", "M"),
#'   sample_id = c("M10", "M2", "M1", "M11")
#' )
#'
#' # One-column natural sort
#' arrange_humanorder(df, sample_id)
#'
#' # Multiple columns, all ascending
#' arrange_humanorder(df, sex, sample_id)
#'
#' # Multiple columns, all descending
#' arrange_humanorder(df, sex, sample_id, .dir = "desc")
#'
#' # Mixed directions by column
#' arrange_humanorder(df, sex, sample_id, .dir = c("asc", "desc"))
#'
#' @export
arrange_humanorder <- function(df, ..., .dir = "asc") {
    # Capture unquoted column names from ...
    cols <- rlang::enquos(...)

    # If no columns were provided, return the input unchanged
    if (length(cols) == 0) {
        return(df)
    }

    # If one direction was supplied, recycle it across all columns
    if (length(.dir) == 1) {
        .dir <- rep(.dir, length(cols))
    }

    # Validate direction length
    if (length(.dir) != length(cols)) {
        stop("`.dir` must be length 1 or the same length as the number of columns.")
    }

    # Standardize case
    .dir <- tolower(.dir)

    # Validate direction values
    if (!all(.dir %in% c("asc", "desc"))) {
        stop("`.dir` must contain only 'asc' or 'desc'.")
    }

    # Sort from the last column to the first column so that the first
    # column remains the highest-priority sort key
    for (i in rev(seq_along(cols))) {
        values_i <- as.character(dplyr::pull(df, !!cols[[i]]))

        ord_i <- gtools::mixedorder(
            values_i,
            decreasing = (.dir[i] == "desc")
        )

        df <- dplyr::slice(df, ord_i)
    }

    df
}
