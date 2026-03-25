#' Summarize a sample × sample correlation (similarity) matrix (row-centric view)
#'
#' @description
#' Given a correlation/similarity matrix where **rows represent samples from
#' dataset 1** and **columns represent samples from dataset 2**, this function
#' produces a tidy, *row-centric* summary. Each output row corresponds to one
#' row of the input matrix (one sample from dataset 1) and answers the question:
#'
#' > 'For this dataset‑1 sample (row), how well does it match its expected
#' > counterpart in dataset 2 (diagonal), and which dataset‑2 sample is its
#' > best overall match (row‑wise maximum)?'
#'
#' This is useful for quality control, especially for detecting potential sample
#' swaps or mixups from the **dataset‑1 perspective**.
#'
#' The input matrix is typically the output of something like:
#' \code{cor(t(X1), t(X2))}, where each row of X1 and X2 is a sample and
#' columns are features (e.g. flattened genotype probabilities).
#'
#' @details
#' Perspective and orientation:
#'
#' - Input: numeric matrix \code{cor_mat} with
#'   - \code{rownames(cor_mat)} = sample IDs from **dataset 1**
#'   - \code{colnames(cor_mat)} = sample IDs from **dataset 2**
#' - Output: tibble with **one row per input row** (i.e. one row per
#'   dataset‑1 sample).
#'
#' Concretely, for each row (sample from dataset 1), the function computes:
#'
#' \itemize{
#'   \item \strong{expected_score}:
#'     The diagonal value (correlation to the sample with the same ID in dataset 2).
#'     This represents how well a sample matches itself.
#'
#'   \item \strong{best_match} and \strong{best_score}:
#'     The highest similarity in the entire row. This may or may not be the
#'     expected (self) match.
#'
#'   \item \strong{second_match} and \strong{second_score}:
#'     The second highest similarity in the row (useful for context and tie detection).
#'
#'   \item \strong{best_nonself} and \strong{best_nonself_score}:
#'     The highest similarity excluding the diagonal (self match). This is
#'     particularly important for swap detection.
#'
#'   \item \strong{delta_best_vs_expected}:
#'     best_score - expected_score.
#'
#'     Interpretation:
#'       - 0   → self is best match
#'       - > 0 → another sample matches better than self (swap candidate)
#'
#'   \item \strong{delta_nonself_vs_expected}:
#'     best_nonself_score - expected_score.
#'
#'     Interpretation:
#'       - < 0 → self matches better than any other sample (good)
#'       - > 0 → some other sample matches better than self (suspicious)
#'
#'   \item \strong{margin_expected_vs_nonself}:
#'     expected_score - best_nonself_score.
#'
#'     Interpretation:
#'       - Large positive → strong separation, clean sample
#'       - Near zero → weak separation, possible contamination or relatedness
#'       - Negative → likely swap
#' }
#'
#' @param cor_mat Numeric similarity / correlation matrix with rownames and
#'   colnames. Rows are interpreted as 'dataset‑1 samples', columns as
#'   'dataset‑2 samples'. The matrix need not be symmetric.
#'
#' @return A tibble with one row per **row** in \code{cor_mat} (one per
#'   dataset‑1 sample).
#'
#' @export
summarize_correlation_matrix <- function(cor_mat) {
    if (!is.matrix(cor_mat)) {
        stop('cor_mat must be a matrix.')
    }

    row_ids <- rownames(cor_mat)
    col_ids <- colnames(cor_mat)

    if (is.null(row_ids) || is.null(col_ids)) {
        stop('cor_mat must have rownames and colnames.')
    }

    n_rows <- nrow(cor_mat)
    n_cols <- ncol(cor_mat)

    #
    # locate the expected (diagonal) column index for each row
    #
    # expected_col_index[i] gives the column index in cor_mat that
    # corresponds to row_ids[i], or NA if no such column exists
    #
    expected_col_index <- match(row_ids, col_ids)

    #
    # get best overall match (row-wise maximum)
    #
    # max.col() finds the column index of the maximum value for each row
    #
    best_col_index <- max.col(cor_mat, ties.method = 'first')

    best_score <- cor_mat[cbind(seq_len(n_rows), best_col_index)]
    best_match <- col_ids[best_col_index]

    #
    # get second-best match
    #
    # to obtain the second-highest value without sorting:
    #   - temporarily mask the best value in each row
    #   - recompute max.col() on the modified matrix
    #
    cor_masked_best <- cor_mat
    cor_masked_best[cbind(seq_len(n_rows), best_col_index)] <- -Inf

    second_col_index <- max.col(cor_masked_best, ties.method = 'first')
    second_score <- cor_masked_best[cbind(seq_len(n_rows), second_col_index)]
    second_match <- col_ids[second_col_index]

    # replace invalid (-Inf) second scores with NA
    invalid_second <- !is.finite(second_score)
    second_score[invalid_second] <- NA_real_
    second_match[invalid_second] <- NA_character_

    #
    # get expected (diagonal) score
    #
    # if the row name exists among column names, extract the diagonal value
    #
    expected_score <- rep(NA_real_, n_rows)

    has_expected <- !is.na(expected_col_index)
    expected_score[has_expected] <-
        cor_mat[cbind(which(has_expected), expected_col_index[has_expected])]

    #
    # get the best non-self match
    #
    # if a diagonal exists, mask it and compute the maximum among
    # remaining values. This identifies the strongest competitor
    #
    best_nonself_match <- rep(NA_character_, n_rows)
    best_nonself_score <- rep(NA_real_, n_rows)

    if (any(has_expected)) {
        cor_masked_diag <- cor_mat
        cor_masked_diag[cbind(which(has_expected),
                              expected_col_index[has_expected])] <- -Inf

        best_nonself_col_index <- max.col(cor_masked_diag,
                                          ties.method = 'first')

        best_nonself_score <-
            cor_masked_diag[cbind(seq_len(n_rows),
                                  best_nonself_col_index)]

        best_nonself_match <- col_ids[best_nonself_col_index]

        invalid_ns <- !is.finite(best_nonself_score)
        best_nonself_score[invalid_ns] <- NA_real_
        best_nonself_match[invalid_ns] <- NA_character_
    }

    #
    # set some derived QC metrics
    #
    delta_best_vs_expected <- best_score - expected_score
    delta_best_vs_second <- best_score - second_score
    delta_nonself_vs_expected <- best_nonself_score - expected_score
    margin_expected_vs_nonself <- expected_score - best_nonself_score

    # TRUE if the row's best match equals itself (diagonal wins)
    best_is_self <- !is.na(best_match) & (best_match == row_ids)

    # a simple swap suspicion rule is:
    #     best match is not self AND best beats expected
    swap_suspect <-
        !is.na(best_match) &
        (best_match != row_ids) &
        !is.na(delta_best_vs_expected) &
        (delta_best_vs_expected > 0)

    #
    # create the tibble
    #
    tibble::tibble(
        sample_id = row_ids,
        expected_match = row_ids,
        expected_score = expected_score,
        best_match = best_match,
        best_score = best_score,
        second_match = second_match,
        second_score = second_score,
        best_nonself = best_nonself_match,
        best_nonself_score = best_nonself_score,
        delta_best_vs_expected = delta_best_vs_expected,
        delta_best_vs_second = delta_best_vs_second,
        delta_nonself_vs_expected = delta_nonself_vs_expected,
        margin_expected_vs_nonself = margin_expected_vs_nonself,
        best_is_self = best_is_self,
        swap_suspect = swap_suspect
    )
}
