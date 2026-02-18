#' Thoroughly compare two genoprobs objects
#'
#' Performs a comprehensive comparison of two R/qtl2 calc_genoprob objects:
#' chromosomes, dimensions, dimension names (samples, founders, markers), and
#' numeric values. Does NOT assume dimension names are in the same order; aligns
#' by names before comparing values (e.g. samples F1, F2, M30 vs M30, F1, F2).
#'
#' @param genoprobs_1 First calc_genoprob object.
#' @param genoprobs_2 Second calc_genoprob object.
#' @param tol Numeric tolerance for value comparison; values are considered
#'   equal if \code{abs(a - b) <= tol} (default 1e-10).
#' @param compare_values If \code{TRUE}, compare numeric values on overlapping
#'   samples/founders/markers; set \code{FALSE} for structure-only comparison.
#'
#' @return A list with components:
#'   \item{chromosomes}{Summary of chromosome overlap and order.}
#'   \item{attributes}{Comparison of qtl2 attributes (crosstype, is_x_chr, alleles, alleleprobs).}
#'   \item{by_chr}{Per-chromosome details (samples, founders, markers, dims, values).}
#'   \item{summary}{Character vector of high-level summary messages.}
#'
#' @export
genoprobs_compare <- function(genoprobs_1,
                              genoprobs_2,
                              tol = 1e-10,
                              compare_values = TRUE) {

    gp1 <- extract_chr_list(genoprobs_1)
    gp2 <- extract_chr_list(genoprobs_2)

    chr1 <- names(gp1)
    chr2 <- names(gp2)
    chr_only_1 <- setdiff(chr1, chr2)
    chr_only_2 <- setdiff(chr2, chr1)
    chr_common <- intersect(chr1, chr2)
    chr_order_match <- length(chr_common) > 0 && identical(chr1[chr1 %in% chr_common], chr2[chr2 %in% chr_common])

    chromosomes <- list(
        in_1 = chr1,
        in_2 = chr2,
        only_in_1 = chr_only_1,
        only_in_2 = chr_only_2,
        common = chr_common,
        order_match = chr_order_match
    )

    # compare attributes
    attrs <- c('crosstype', 'is_x_chr', 'alleles', 'alleleprobs')
    attr_compare <- list()
    for(a in attrs) {
        v1 <- attr(genoprobs_1, a, exact = TRUE)
        v2 <- attr(genoprobs_2, a, exact = TRUE)
        attr_compare[[a]] <- list(
            in_1 = !is.null(v1),
            in_2 = !is.null(v2),
            match = if(is.null(v1) && is.null(v2)) TRUE else identical(v1, v2)
        )
        if(!attr_compare[[a]]$match && !is.null(v1) && !is.null(v2)) {
            attr_compare[[a]]$value_1 <- v1
            attr_compare[[a]]$value_2 <- v2
        }
    }
    attributes_out <- attr_compare

    by_chr <- setNames(vector('list', length(chr_common)), chr_common)
    summary_lines <- character()

    for(chr in chr_common) {
        pr1 <- gp1[[chr]]
        pr2 <- gp2[[chr]]

        if(length(dim(pr1)) != 3 || length(dim(pr2)) != 3) {
            by_chr[[chr]] <- list(error = paste0('Chr ', chr, ': one or both arrays are not 3D'))
            next
        }

        dn1 <- dimnames(pr1)
        dn2 <- dimnames(pr2)
        if(is.null(dn1) || is.null(dn2)) {
            by_chr[[chr]] <- list(error = paste0('Chr ', chr, ': one or both arrays have NULL dimnames'))
            next
        }
        s1 <- if(is.null(dn1[[1]])) character(0) else dn1[[1]]
        s2 <- if(is.null(dn2[[1]])) character(0) else dn2[[1]]
        f1 <- if(is.null(dn1[[2]])) character(0) else dn1[[2]]
        f2 <- if(is.null(dn2[[2]])) character(0) else dn2[[2]]
        m1 <- if(is.null(dn1[[3]])) character(0) else dn1[[3]]
        m2 <- if(is.null(dn2[[3]])) character(0) else dn2[[3]]

        s_only_1 <- setdiff(s1, s2)
        s_only_2 <- setdiff(s2, s1)
        s_common <- intersect(s1, s2)
        s_order_match <- length(s_common) > 0 && identical(s1[s1 %in% s_common], s2[s2 %in% s_common])

        f_only_1 <- setdiff(f1, f2)
        f_only_2 <- setdiff(f2, f1)
        f_common <- intersect(f1, f2)
        f_order_match <- length(f_common) > 0 && identical(f1[f1 %in% f_common], f2[f2 %in% f_common])

        m_only_1 <- setdiff(m1, m2)
        m_only_2 <- setdiff(m2, m1)
        m_common <- intersect(m1, m2)
        m_order_match <- length(m_common) > 0 && identical(m1[m1 %in% m_common], m2[m2 %in% m_common])

        dim_match <- identical(dim(pr1), dim(pr2))

        chr_result <- list(
            dimensions = list(
                dim_1 = dim(pr1),
                dim_2 = dim(pr2),
                match = dim_match
            ),
            samples = list(
                only_in_1 = s_only_1,
                only_in_2 = s_only_2,
                common = s_common,
                order_match = s_order_match,
                n_1 = length(s1),
                n_2 = length(s2)
            ),
            founders = list(
                only_in_1 = f_only_1,
                only_in_2 = f_only_2,
                common = f_common,
                order_match = f_order_match,
                n_1 = length(f1),
                n_2 = length(f2)
            ),
            markers = list(
                only_in_1 = m_only_1,
                only_in_2 = m_only_2,
                common = m_common,
                order_match = m_order_match,
                n_1 = length(m1),
                n_2 = length(m2)
            )
        )

        # value comparison: align by names, then compare
        values_summary <- NULL
        if(compare_values && length(s_common) > 0 && length(f_common) > 0 && length(m_common) > 0) {
            i_s1 <- match(s_common, s1)
            i_s2 <- match(s_common, s2)
            i_f1 <- match(f_common, f1)
            i_f2 <- match(f_common, f2)
            i_m1 <- match(m_common, m1)
            i_m2 <- match(m_common, m2)

            v1 <- pr1[i_s1, i_f1, i_m1, drop = FALSE]
            v2 <- pr2[i_s2, i_f2, i_m2, drop = FALSE]

            stopifnot(identical(dim(v1), dim(v2)))
            diff <- as.vector(v1) - as.vector(v2)
            finite_diff <- diff[is.finite(diff)]
            na1 <- is.na(v1)
            na2 <- is.na(v2)
            na_mismatch <- sum(na1 != na2)
            abs_diff <- abs(diff[is.finite(diff)])

            values_summary <- list(
                n_cells = length(v1),
                n_compared = length(finite_diff),
                na_mismatch = na_mismatch,
                max_abs_diff = if(length(abs_diff)) max(abs_diff) else NA_real_,
                mean_abs_diff = if(length(abs_diff)) mean(abs_diff) else NA_real_,
                identical = na_mismatch == 0 && all(abs_diff <= tol)
            )
            chr_result$values <- values_summary
        } else if(compare_values) {
            chr_result$values <- list(
                skipped = TRUE,
                reason = 'No overlapping samples, founders, and/or markers for value comparison'
            )
        }

        by_chr[[chr]] <- chr_result
    }

    summary_lines <- c(
        sprintf('Chromosomes: %d in both, %d only in 1, %d only in 2',
                length(chr_common), length(chr_only_1), length(chr_only_2)),
        sprintf('Chromosome order match: %s', chr_order_match)
    )
    if(length(chr_common) > 0) {
        order_diff <- vapply(by_chr, function(x) {
            if(!is.null(x$error)) return(FALSE)
            !x$samples$order_match || !x$founders$order_match || !x$markers$order_match
        }, FALSE)
        n_order_diff <- sum(order_diff)
        summary_lines <- c(summary_lines,
                          sprintf('Chr with dimension-order mismatch: %d', n_order_diff))
        if(compare_values) {
            n_ident <- sum(vapply(by_chr, function(x) {
                if(is.null(x$values) || isTRUE(x$values$skipped)) return(NA)
                isTRUE(x$values$identical)
            }, NA), na.rm = TRUE)
            n_total <- sum(vapply(by_chr, function(x) {
                if(is.null(x$values) || isTRUE(x$values$skipped)) return(0)
                1
            }, 0))
            summary_lines <- c(summary_lines,
                              sprintf('Chr with identical values (tol=%g): %d / %d', tol, n_ident, n_total))
        }
    }

    list(
        chromosomes = chromosomes,
        attributes = attributes_out,
        by_chr = by_chr,
        summary = summary_lines
    )
}
