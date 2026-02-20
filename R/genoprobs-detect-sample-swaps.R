#' Detect likely sample swaps from a similarity matrix or two genoprob objects
#'
#' Generic function to detect sample swaps. Can accept either:
#' \itemize{
#'   \item A similarity matrix (from \code{genoprobs_compute_similarity}) - dispatches to \code{.matrix} method
#'   \item Two calc_genoprob objects - dispatches to \code{.calc_genoprob} method (computes similarity first)
#' }
#'
#' @param x Either a similarity matrix or a calc_genoprob object.
#' @param ... Additional arguments passed to methods. See method-specific documentation.
#'
#' @return data.frame with columns:
#' \describe{
#'   \item{col_sample}{Sample ID from the column set (genoprobs_2 or second object)}
#'   \item{expected_row_sample}{Expected matching sample from row set}
#'   \item{labeled_r}{Similarity to expected match}
#'   \item{best_row_sample}{Actual best matching sample from row set}
#'   \item{best_r}{Similarity to best match}
#'   \item{second_row_sample}{Second-best matching sample}
#'   \item{second_r}{Similarity to second-best match}
#'   \item{delta_best_second}{Difference: best - second}
#'   \item{delta_best_labeled}{Difference: best - labeled}
#'   \item{flag}{Match status: "match", "match_ambiguous", "mismatch_confident", "mismatch_ambiguous", or "unmapped"}
#'   \item{reciprocal_pair}{For swaps, identifier like "swap:A<->B"}
#' }
#'
#' @examples
#' \dontrun{
#' # Method 1: From similarity matrix
#' sim_result <- genoprobs_compute_similarity(gp_a, gp_b)
#' swaps <- genoprobs_detect_sample_swaps(sim_result$sim, sample_map)
#'
#' # Method 2: From genoprob objects (convenience wrapper)
#' swaps <- genoprobs_detect_sample_swaps(gp_a, gp_b, sample_map, metric = 'pearson')
#' }
#' @export
genoprobs_detect_sample_swaps <- function(x, ...) {
    UseMethod('genoprobs_detect_sample_swaps')
}

#' @rdname genoprobs_detect_sample_swaps
#' @param x For \code{.matrix} method: similarity matrix with rows = one set of samples
#'   (e.g. one platform/build), cols = the other set. Each cell x[i,j] = similarity
#'   between row sample i and col sample j. Must have rownames and colnames.
#' @param sample_map Optional data.frame with columns c("col_sample", "row_sample") giving
#'   expected pairing: for each col sample, which row sample it is expected to match.
#' @param min_delta_best_labeled Flag mismatch if best - labeled >= this (default 0.05)
#' @param min_delta_best_second Require best - second >= this for "confident" (default 0.02)
#' @export
genoprobs_detect_sample_swaps.matrix <- function(x,
                                                sample_map = NULL,
                                                min_delta_best_labeled = 0.05,
                                                min_delta_best_second  = 0.02,
                                                ...) {
  sim <- x

  # =========================================================================
  # STEP 1: Validate input and extract sample names
  # =========================================================================
  # The similarity matrix sim has:
  #   - ROWS = one set of samples (row set)
  #   - COLS = the other set (col set)
  # Each cell sim[i,j] = similarity between row sample i and col sample j.
  # We rely on rownames/colnames as the sample IDs.
  #
  # NOTE: This function assumes "larger = more similar". For Pearson/cosine,
  # larger is better (1 is perfect match). If you ever pass a distance matrix,
  # you must convert it (e.g., similarity = -distance) before using this.

  if (!is.matrix(sim)) {
    stop("x must be a matrix for the .matrix method.")
  }
  if (is.null(rownames(sim)) || is.null(colnames(sim))) {
    stop("sim must have rownames (row samples) and colnames (col samples).")
  }

  row_samples <- rownames(sim)
  col_samples <- colnames(sim)

  # Guard against duplicated sample IDs. Duplicates make mapping & swap logic ambiguous.
  if (anyDuplicated(row_samples)) stop("Row sample names (rownames(sim)) contain duplicates.")
  if (anyDuplicated(col_samples)) stop("Col sample names (colnames(sim)) contain duplicates.")

  # =========================================================================
  # STEP 2: Build or validate sample_map (col -> row expected pairing)
  # =========================================================================
  # sample_map tells us: for each col sample, which row sample it is EXPECTED to match.
  #
  # If sample_map is NULL, we assume "same IDs should match", but still return output
  # for all columns; columns without a same-named row become unmapped.

  if (is.null(sample_map)) {
    shared <- intersect(row_samples, col_samples)
    sample_map <- data.frame(col_sample = shared,
                             row_sample = shared,
                             stringsAsFactors = FALSE)
  } else {
    nm <- names(sample_map)

    # Only accept explicit col_sample/row_sample
    if (!all(c("col_sample", "row_sample") %in% nm)) {
      stop('sample_map must have columns "col_sample" and "row_sample".')
    }
    sample_map <- sample_map[, c("col_sample", "row_sample"), drop = FALSE]

    # Sanity checks: duplicates usually indicate metadata problems.
    if (anyDuplicated(sample_map$col_sample)) {
      stop("sample_map has duplicated col_sample entries; mapping must be 1-to-1 by col_sample.")
    }
    # Duplicated row_sample may or may not be allowed; warn by default because swaps become ambiguous.
    if (anyDuplicated(sample_map$row_sample)) {
      warning("sample_map has duplicated row_sample values. Reciprocal swap detection may be ambiguous.")
    }

    # Keep only mappings present in the matrix
    sample_map <- sample_map[sample_map$col_sample %in% col_samples &
                               sample_map$row_sample %in% row_samples, , drop = FALSE]
  }

  # expected_row[col] = expected row sample name (or NA if unmapped)
  expected_row <- setNames(sample_map$row_sample, sample_map$col_sample)
  expected_row <- expected_row[col_samples]  # reorder to match matrix columns; introduces NA for unmapped

  # =========================================================================
  # STEP 3: For each col sample, find best and second-best row matches
  # =========================================================================
  # We do this in a vectorized way:
  #   best_idx[j]   = which.max(sim[,j])  (ties broken by first)
  #   second_idx[j] = best of remaining after masking out best
  #
  # IMPORTANT EDGE CASE:
  # - If an entire column is NA (no similarities computable), best/second should be NA.

  # Identify columns that contain at least one finite value
  col_has_value <- vapply(seq_along(col_samples), function(j) any(is.finite(sim[, j])), logical(1))

  best_idx <- rep(NA_integer_, length(col_samples))
  best_r   <- rep(NA_real_,    length(col_samples))

  # max.col(m) returns column index of max in each ROW; we need row index of max in each COLUMN.
  # So use max.col(t(sim_sub)) to get, per column of sim_sub, which row has the max.
  if (any(col_has_value)) {
    sim_sub <- sim[, col_has_value, drop = FALSE]
    best_idx[col_has_value] <- max.col(t(sim_sub), ties.method = "first")
    best_r[col_has_value]   <- sim[cbind(best_idx[col_has_value], which(col_has_value))]
  }

  best_row_sample <- ifelse(is.na(best_idx), NA_character_, row_samples[best_idx])

  # Second-best: mask out the best entry per column, then take max again.
  second_idx <- rep(NA_integer_, length(col_samples))
  second_r   <- rep(NA_real_,    length(col_samples))

  if (any(col_has_value) && nrow(sim) >= 2) {
    sim2 <- sim[, col_has_value, drop = FALSE]
    # mask best entries (one per column)
    j_good <- which(col_has_value)
    sim2[cbind(best_idx[col_has_value], seq_along(j_good))] <- -Inf
    # row index of max in each column: max.col(t(sim2))
    second_idx[col_has_value] <- max.col(t(sim2), ties.method = "first")
    # if original column had only one finite entry, second could be -Inf; treat as NA
    tmp_second_r <- sim[cbind(second_idx[col_has_value], j_good)]
    tmp_second_r[!is.finite(tmp_second_r)] <- NA_real_
    second_r[col_has_value] <- tmp_second_r
  }

  second_row_sample <- ifelse(is.na(second_idx), NA_character_, row_samples[second_idx])

  # Optional tie diagnostics (cheap + useful): how many rows share the best score?
  # This helps interpret "ambiguous" cases.
  n_at_best <- vapply(seq_along(col_samples), function(j) {
    if (!col_has_value[j]) return(NA_integer_)
    mx <- max(sim[, j], na.rm = TRUE)
    sum(sim[, j] == mx, na.rm = TRUE)
  }, integer(1))

  # =========================================================================
  # STEP 4: Similarity of each col sample to its LABELED (expected) row
  # =========================================================================
  labeled_r <- vapply(seq_along(col_samples), function(j) {
    r <- expected_row[[j]]
    if (is.na(r)) return(NA_real_)
    sim[r, col_samples[j]]
  }, numeric(1))

  # =========================================================================
  # STEP 5: Compute deltas used for flagging
  # =========================================================================
  # delta_best_second:
  #   - Large means best is clearly better than runner-up (confident assignment)
  #   - Small means ambiguous (best ~ second)
  #
  # delta_best_labeled:
  #   - Large means best is much better than expected label -> possible swap/mislabel

  delta_best_second  <- best_r - second_r
  delta_best_labeled <- best_r - labeled_r

  # =========================================================================
  # STEP 6: Assign flags based on match quality and deltas
  # =========================================================================
  flag <- rep("match", length(col_samples))

  # If no expected mapping -> unmapped
  flag[is.na(expected_row)] <- "unmapped"

  # If the column has no usable similarity values, we also call it unmapped-like.
  # (You could choose a different label; "unmapped" is convenient and conservative.)
  flag[!col_has_value] <- "unmapped"

  # For mapped columns, apply mismatch/match ambiguity logic.
  mapped <- !is.na(expected_row) & col_has_value

  # Confident mismatch: best != expected AND best beats expected by enough AND best beats second by enough
  flag[mapped &
         (best_row_sample != expected_row) &
         (delta_best_labeled >= min_delta_best_labeled) &
         (delta_best_second  >= min_delta_best_second)] <- "mismatch_confident"

  # Ambiguous mismatch: best != expected AND best beats expected by enough BUT does NOT clearly beat second
  flag[mapped &
         (best_row_sample != expected_row) &
         (delta_best_labeled >= min_delta_best_labeled) &
         (!is.na(delta_best_second) & delta_best_second < min_delta_best_second)] <- "mismatch_ambiguous"

  # Ambiguous match: best == expected but second is close (or best is tied)
  flag[mapped &
         (best_row_sample == expected_row) &
         (
           (!is.na(delta_best_second) & delta_best_second < min_delta_best_second) |
             (!is.na(n_at_best) & n_at_best > 1L)
         )] <- "match_ambiguous"

  # =========================================================================
  # STEP 7: Build the output data.frame
  # =========================================================================
  out <- data.frame(
    col_sample           = col_samples,
    expected_row_sample  = unname(expected_row),
    labeled_r            = labeled_r,
    best_row_sample      = best_row_sample,
    best_r               = best_r,
    second_row_sample    = second_row_sample,
    second_r             = second_r,
    delta_best_second    = delta_best_second,
    delta_best_labeled   = delta_best_labeled,
    n_at_best            = n_at_best,   # NEW: helps interpret ambiguous/ties
    flag                 = flag,
    stringsAsFactors     = FALSE
  )

  # =========================================================================
  # STEP 8: Identify reciprocal swap pairs (A↔B swapped)
  # =========================================================================
  # A reciprocal swap means:
  #   - Col sample c1 best-matches row sample that is the expected row for c2
  #   - Col sample c2 best-matches row sample that is the expected row for c1
  #
  # NOTE: This logic assumes the mapping is roughly 1-to-1. If multiple columns map to
  # the same expected row, you can still get pairs but they may not be unique.

  reciprocal <- rep(NA_character_, nrow(out))

  idx_mis <- which(out$flag %in% c("mismatch_confident", "mismatch_ambiguous") &
                     !is.na(out$expected_row_sample) &
                     !is.na(out$best_row_sample))

  best_row <- setNames(out$best_row_sample, out$col_sample)
  exp_row  <- setNames(out$expected_row_sample, out$col_sample)

  for (i in idx_mis) {
    c1 <- out$col_sample[i]
    candidates <- names(exp_row)[exp_row == best_row[[c1]]]

    for (c2 in candidates) {
      if (!is.na(best_row[[c2]]) && best_row[[c2]] == exp_row[[c1]]) {
        rid <- paste0("swap:", paste(sort(c(c1, c2)), collapse = "<->"))
        reciprocal[match(c1, out$col_sample)] <- rid
        reciprocal[match(c2, out$col_sample)] <- rid
      }
    }
  }

  out$reciprocal_pair <- reciprocal
  out
}


#' @rdname genoprobs_detect_sample_swaps
#' @param x For \code{.calc_genoprob} method: first calc_genoprob object.
#' @param genoprobs_2 Second calc_genoprob object (required when x is a calc_genoprob object).
#' @param sample_map Optional data.frame with columns c("col_sample", "row_sample") giving
#'   expected pairing: for each col sample, which row sample it is expected to match.
#' @param metric Similarity metric to use when computing similarity matrix.
#'   One of \code{'pearson'} or \code{'cosine'} (default \code{'pearson'}).
#'   Only used when x is a calc_genoprob object.
#' @param min_delta_best_labeled Flag mismatch if best - labeled >= this (default 0.05)
#' @param min_delta_best_second Require best - second >= this for "confident" (default 0.02)
#' @export
genoprobs_detect_sample_swaps.calc_genoprob <- function(x,
                                                         genoprobs_2,
                                                         sample_map = NULL,
                                                         metric = c('pearson', 'cosine'),
                                                         min_delta_best_labeled = 0.05,
                                                         min_delta_best_second  = 0.02,
                                                         ...) {
    genoprobs_1 <- x
    metric <- match.arg(metric)

    # Compute similarity matrix first
    sim_result <- genoprobs_compute_similarity(genoprobs_1, genoprobs_2, metric = metric)

    # Then call the matrix method
    genoprobs_detect_sample_swaps.matrix(
        sim_result$sim,
        sample_map = sample_map,
        min_delta_best_labeled = min_delta_best_labeled,
        min_delta_best_second = min_delta_best_second
    )
}

#' @rdname genoprobs_detect_sample_swaps
#' @export
genoprobs_detect_sample_swaps.default <- function(x, ...) {
    stop(
        "genoprobs_detect_sample_swaps() does not support objects of class '",
        paste(class(x), collapse = "/"), "'.\n",
        "Expected either a similarity matrix or a calc_genoprob object.\n",
        "For calc_genoprob objects, also provide genoprobs_2 argument."
    )
}
