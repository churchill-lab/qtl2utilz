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
#'   Legacy columns "gbrs" and "array" are also accepted.
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
    #   - ROWS = one set of samples (e.g., genoprobs from platform/build A: MUGA GRCm38, etc.)
    #   - COLS = the other set (e.g., genoprobs from platform/build B: MUGA GRCm39, GBRS, etc.)
    # Each cell sim[i,j] = correlation/similarity between row sample i and col sample j.
    # Row and col can be any two comparable sources (e.g. two arrays, or array vs haplotype).

    # Require both row and column names to exist; we use them as sample IDs
    if (is.null(rownames(sim)) || is.null(colnames(sim))) {
        stop("sim must have rownames (row samples) and colnames (col samples).")
    }

    # Extract the ordered list of sample IDs from matrix dimensions
    row_samples <- rownames(sim)  # e.g., c("DO001", "DO002", ...)
    col_samples <- colnames(sim)  # e.g., c("DO001", "DO002", ...)

    # =========================================================================
    # STEP 2: Build or validate sample_map (col -> row expected pairing)
    # =========================================================================
    # sample_map tells us: "col sample X is EXPECTED to match row sample Y"
    # (e.g., from metadata: "DO001 from source B" should correspond to "DO001 from source A")

    if (is.null(sample_map)) {
        # No mapping provided: assume col and row use same IDs, and each col
        # sample is expected to match the row sample with the SAME name
        shared <- intersect(row_samples, col_samples)  # samples present in BOTH
        sample_map <- data.frame(col_sample = shared, row_sample = shared, stringsAsFactors = FALSE)
    } else {
        # User provided a mapping: validate it has required columns
        # Accept legacy names "gbrs"/"array" as well as generic "col_sample"/"row_sample"
        nm <- names(sample_map)
        if (all(c("col_sample", "row_sample") %in% nm)) {
            sample_map <- sample_map[, c("col_sample", "row_sample"), drop = FALSE]
        } else if (all(c("gbrs", "array") %in% nm)) {
            sample_map <- data.frame(col_sample = sample_map$gbrs, row_sample = sample_map$array, stringsAsFactors = FALSE)
        } else {
            stop("sample_map must have columns 'col_sample' and 'row_sample' (or legacy 'gbrs' and 'array').")
        }
        # Keep only rows where BOTH col and row sample exist in the sim matrix
        # (filter out any mappings to samples not in our data)
        sample_map <- sample_map[sample_map$col_sample %in% col_samples &
                                     sample_map$row_sample %in% row_samples, , drop = FALSE]
    }

    # Create a named vector: for each col sample, what row sample is it SUPPOSED to match?
    # expected_row["DO001"] = "DO001" means: "DO001 (col) should match DO001 (row)"
    expected_row <- setNames(sample_map$row_sample, sample_map$col_sample)
    # Reorder to match the column order of sim; samples not in sample_map become NA
    expected_row <- expected_row[col_samples]

    # =========================================================================
    # STEP 3: For each col sample, find the best and second-best row match
    # =========================================================================
    # We iterate over columns of sim. Column j = similarities of all row samples
    # to col sample j. We find which row sample has highest similarity (best)
    # and which has second-highest (second).

    best_row_sample   <- character(length(col_samples))  # best matching row sample name
    best_r            <- numeric(length(col_samples))     # similarity of best match
    second_row_sample <- character(length(col_samples))  # second-best row sample name
    second_r          <- numeric(length(col_samples))     # similarity of second-best match

    for (j in seq_along(col_samples)) {
        # Get column j: similarity of every row sample to this col sample
        col_vec <- sim[, j]

        # Order indices by similarity, descending (highest first), NAs last
        ord <- order(col_vec, decreasing = TRUE, na.last = TRUE)

        # First element = index of row sample with highest similarity
        best_row_sample[j] <- row_samples[ord[1]]
        best_r[j] <- col_vec[ord[1]]  # the actual similarity value

        if (length(ord) >= 2) {
            # Second element = second-best matching row sample
            second_row_sample[j] <- row_samples[ord[2]]
            second_r[j] <- col_vec[ord[2]]
        } else {
            # Only one row sample (or none); no "second" exists
            second_row_sample[j] <- NA_character_
            second_r[j] <- NA_real_
        }
    }

    # =========================================================================
    # STEP 4: Get the similarity of each col sample to its LABELED (expected) row
    # =========================================================================
    # For each col sample j: look up sim[expected_row[j], j]
    # i.e., the similarity between col j and the row sample we EXPECT it to match

    labeled_r <- vapply(seq_along(col_samples), function(j) {
        r <- expected_row[[j]]  # the row sample we expect to match
        if (is.na(r)) return(NA_real_)  # no expected match (unmapped)
        # sim[row_sample, col_sample] = similarity between that pair
        sim[r, col_samples[j]]
    }, numeric(1))

    # =========================================================================
    # STEP 5: Compute deltas (differences) used for flagging
    # =========================================================================
    # delta_best_second:  How much better is best than second?
    #   - Large = one clear winner (confident call)
    #   - Small = best and second are close (ambiguous)
    delta_best_second  <- best_r - second_r

    # delta_best_labeled: How much better is best than the labeled (expected) match?
    #   - Large = the best match is clearly NOT the expected one (suggests swap/mislabel)
    #   - Small or negative = best and labeled are similar (probably correct)
    delta_best_labeled <- best_r - labeled_r

    # =========================================================================
    # STEP 6: Assign flags based on match quality and deltas
    # =========================================================================
    # Start with "match" for everyone; then overwrite based on conditions
    flag <- rep("match", length(col_samples))

    # unmapped: no expected row sample for this col sample (can't check)
    flag[is.na(expected_row)] <- "unmapped"

    # mismatch_confident: best != expected AND we're confident it's a real mismatch
    #   - delta_best_labeled >= threshold: best is clearly better than expected
    #   - delta_best_second >= threshold: best is clearly better than second (not ambiguous)
    flag[!is.na(expected_row) & (best_row_sample != expected_row) &
             (delta_best_labeled >= min_delta_best_labeled) &
             (delta_best_second  >= min_delta_best_second)] <- "mismatch_confident"

    # mismatch_ambiguous: best != expected but second-best is close to best
    #   - We suspect a mismatch, but could be noise (best and second nearly tied)
    flag[!is.na(expected_row) & (best_row_sample != expected_row) &
             (delta_best_labeled >= min_delta_best_labeled) &
             (delta_best_second  <  min_delta_best_second)] <- "mismatch_ambiguous"

    # match_ambiguous: best == expected (correct!) but second is close to best
    #   - Probably correct, but if labels were wrong we couldn't tell (ambiguous)
    flag[!is.na(expected_row) & (best_row_sample == expected_row) &
             (delta_best_second  <  min_delta_best_second)] <- "match_ambiguous"

    # =========================================================================
    # STEP 7: Build the output data.frame
    # =========================================================================
    out <- data.frame(
        col_sample = col_samples,
        expected_row_sample = unname(expected_row),  # expected row match
        labeled_r = labeled_r,                       # similarity to expected
        best_row_sample = best_row_sample,            # actual best row match
        best_r = best_r,                              # similarity to best
        second_row_sample = second_row_sample,        # second-best row match
        second_r = second_r,                          # similarity to second-best
        delta_best_second = delta_best_second,        # best - second
        delta_best_labeled = delta_best_labeled,      # best - labeled
        flag = flag,
        stringsAsFactors = FALSE
    )

    # =========================================================================
    # STEP 8: Identify reciprocal swap pairs (A↔B swapped)
    # =========================================================================
    # A reciprocal swap: col A's best match is row B, and col B's best match is row A
    # i.e., best_row[A]==expected_row[B] and best_row[B]==expected_row[A]
    # We tag both A and B with the same reciprocal_pair ID (e.g., "swap:A<->B")

    reciprocal <- rep(NA_character_, nrow(out))  # will hold "swap:X<->Y" or NA

    # Only consider rows that are mismatches (confident or ambiguous) and have an expected sample
    idx_mis <- which(out$flag %in% c("mismatch_confident", "mismatch_ambiguous") &
                         !is.na(out$expected_row_sample))

    # Named vectors for quick lookup:
    # best_row[col] = row sample that best matches this col sample
    # exp_row[col]  = row sample we EXPECT to match this col sample (labeled)
    best_row <- setNames(out$best_row_sample, out$col_sample)
    exp_row  <- setNames(out$expected_row_sample, out$col_sample)

    for (i in idx_mis) {
        c1 <- out$col_sample[i]  # first col sample in potential swap pair

        # For a swap, we need c2 such that:
        #   best_row[c1] == exp_row[c2]  (c1's best match is c2's expected row)
        #   best_row[c2] == exp_row[c1]  (c2's best match is c1's expected row)
        # So: find all c2 whose EXPECTED row equals c1's BEST match
        candidates <- names(exp_row)[exp_row == best_row[[c1]]]

        for (c2 in candidates) {
            # Check the reverse: does c2's best match equal c1's expected row?
            if (!is.na(best_row[[c2]]) && best_row[[c2]] == exp_row[[c1]]) {
                # Yes! c1 and c2 form a reciprocal swap pair
                rid <- paste0("swap:", paste(sort(c(c1, c2)), collapse = "<->"))
                # Assign same ID to both rows (sort ensures consistent ordering)
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
#'   Legacy columns "gbrs" and "array" are also accepted.
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
