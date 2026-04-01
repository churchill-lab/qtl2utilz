#' Arrange rows in human/natural sort order across one or more columns
#'
#' Sorts a data frame using "human" ordering rather than strict
#' lexicographic ordering. For example, values like `sample2`, `sample10`,
#' and `sample100` will sort in natural order rather than plain character
#' order.
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
#' The function works by:
#' 1. Capturing the requested columns with tidy evaluation.
#' 2. Pulling each column from `df` and converting it to character.
#' 3. Converting each column to a mixed/natural rank using
#'    [gtools::mixedrank()].
#' 4. Negating the rank for descending columns.
#' 5. Combining all ranking vectors with [base::order()] and using the
#'    resulting row positions inside [dplyr::slice()].
#'
#' This is useful for columns containing embedded numbers such as:
#' - `M1`, `M2`, `M10`
#' - `sample_1`, `sample_2`, `sample_11`
#' - chromosome-like labels such as `chr1`, `chr2`, `chr10`
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
    # Capture the unquoted column names supplied through ...
    cols <- rlang::enquos(...)

    # If no columns were given, return the input unchanged
    if (length(cols) == 0) {
        return(df)
    }

    # If a single direction was supplied, recycle it across all columns
    if (length(.dir) == 1) {
        .dir <- rep(.dir, length(cols))
    }

    # Ensure .dir is either length 1 or matches the number of columns
    if (length(.dir) != length(cols)) {
        stop("`.dir` must be length 1 or the same length as the number of columns.")
    }

    # Standardize case before validation
    .dir <- tolower(.dir)

    # Validate accepted direction values
    if (!all(.dir %in% c("asc", "desc"))) {
        stop("`.dir` must contain only 'asc' or 'desc'.")
    }

    # Pull each selected column from df and coerce to character so
    # gtools::mixedrank() can apply natural/human ordering consistently
    order_inputs <- purrr::map(cols, ~ {
        as.character(dplyr::pull(df, !!.x))
    })

    # Convert each column to a ranking vector suitable for base::order().
    # For descending columns, negate the rank so larger values come first.
    for (i in seq_along(order_inputs)) {
        rank_i <- gtools::mixedrank(order_inputs[[i]])

        if (.dir[i] == "desc") {
            rank_i <- -rank_i
        }

        order_inputs[[i]] <- rank_i
    }

    # Combine all ranking vectors into one row order and
    # return df with rows rearranged accordingly
    df |>
        dplyr::slice(do.call(order, order_inputs))
}
