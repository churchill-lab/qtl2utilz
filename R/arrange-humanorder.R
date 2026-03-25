#' Sort a data frame by an ID column like F1, F11, F100
#' so that F11 comes before F100
#' Sort a data frame by an ID column like F1, F11, F100
#' using gtools::mixedorder(), which is safer than mixedsort() for data frames
#' @export
arrange_humanorder <- function(df, col) {
    col_quo <- rlang::enquo(col)

    df |>
        dplyr::slice(
            gtools::mixedorder(as.character(dplyr::pull(df, !!col_quo)))
        )
}
