#' Sort a data frame by natural (human) order on one column
#'
#' Uses \code{gtools::mixedorder()} so identifiers like \code{F1}, \code{F11},
#' \code{F100} sort in intuitive order (not lexicographic).
#'
#' @param df A data frame.
#' @param col Unquoted column name (tidy-eval), as in \code{dplyr::arrange()}.
#'
#' @return \code{df} with rows reordered by mixed/natural order on \code{col}.
#'
#' @export
arrange_humanorder <- function(df, col) {
    col_quo <- rlang::enquo(col)

    df |>
        dplyr::slice(
            gtools::mixedorder(as.character(dplyr::pull(df, !!col_quo)))
        )
}
