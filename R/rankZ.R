#' Rank-based inverse normal transformation.
#'
#' Converts a numeric vector to approximate normality using ranks. Ranks are
#' scaled to (0, 1) and passed through \code{\link[stats]{qnorm}}. Useful for
#' transforming skewed expression or phenotype data before QTL mapping.
#'
#' @param x Numeric vector; \code{NA}s are preserved.
#'
#' @return Numeric vector of the same length as \code{x}; \code{NA} where
#'   \code{x} was \code{NA}, otherwise normal quantiles.
#'
#' @export
rankZ <- function(x) {
    x <- rank(x, na.last = 'keep', ties.method = 'average') / (sum(!is.na(x)) + 1)
    qnorm(x)
}
