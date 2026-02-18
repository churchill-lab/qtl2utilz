#' Detect likely sample swaps from a similarity matrix
#'
#' @param sim Similarity matrix: rows=array/DNA samples, cols=GBRS/RNA samples
#' @param sample_map Optional data.frame with columns c("gbrs","array")
#' @param min_delta_best_labeled Flag mismatch if best - labeled >= this (default 0.05)
#' @param min_delta_best_second Require best - second >= this for "confident" (default 0.02)
#'
#' @return data.frame with best/second matches, deltas, flags, and reciprocal swap tags
#' @export
genoprobs_detect_sample_swaps <- function(sim,
                                          sample_map = NULL,
                                          min_delta_best_labeled = 0.05,
                                          min_delta_best_second  = 0.02) {

    # =========================================================================
    # STEP 1: Validate input and extract sample names
    # =========================================================================
    # The similarity matrix sim has:
    #   - ROWS = array samples (e.g., genotyping array / DNA samples)
    #   - COLS = GBRS samples (e.g., RNA-seq or haplotype-reconstructed samples)
    # Each cell sim[i,j] = correlation/similarity between array sample i and GBRS sample j

    # Require both row and column names to exist; we use them as sample IDs
    if (is.null(rownames(sim)) || is.null(colnames(sim))) {
        stop("sim must have rownames (array samples) and colnames (gbrs samples).")
    }

    # Extract the ordered list of sample IDs from matrix dimensions
    array_samples <- rownames(sim)  # e.g., c("DO001", "DO002", ...)
    gbrs_samples  <- colnames(sim)  # e.g., c("DO001", "DO002", ...)

    # =========================================================================
    # STEP 2: Build or validate sample_map (gbrs -> array expected pairing)
    # =========================================================================
    # sample_map tells us: "GBRS sample X is EXPECTED to match array sample Y"
    # (e.g., from metadata: "DO001 RNA" should correspond to "DO001 Array")

    if (is.null(sample_map)) {
        # No mapping provided: assume gbrs and array use same IDs, and each gbrs
        # sample is expected to match the array sample with the SAME name
        shared <- intersect(array_samples, gbrs_samples)  # samples present in BOTH
        sample_map <- data.frame(gbrs = shared, array = shared, stringsAsFactors = FALSE)
    } else {
        # User provided a mapping: validate it has required columns
        if (!all(c("gbrs", "array") %in% names(sample_map))) {
            stop("sample_map must have columns 'gbrs' and 'array'.")
        }
        # Keep only rows where BOTH gbrs and array sample exist in the sim matrix
        # (filter out any mappings to samples not in our data)
        sample_map <- sample_map[sample_map$gbrs %in% gbrs_samples &
                                     sample_map$array %in% array_samples, , drop = FALSE]
    }

    # Create a named vector: for each gbrs sample, what array sample is it SUPPOSED to match?
    # labeled_array["DO001"] = "DO001" means: "DO001 (gbrs) should match DO001 (array)"
    labeled_array <- setNames(sample_map$array, sample_map$gbrs)
    # Reorder to match the column order of sim; samples not in sample_map become NA
    labeled_array <- labeled_array[gbrs_samples]

    # =========================================================================
    # STEP 3: For each GBRS sample, find the best and second-best array match
    # =========================================================================
    # We iterate over columns of sim. Column j = similarities of all array samples
    # to GBRS sample j. We find which array sample has highest similarity (best)
    # and which has second-highest (second).

    best_array   <- character(length(gbrs_samples))  # best matching array sample name
    best_r       <- numeric(length(gbrs_samples))    # similarity of best match
    second_array <- character(length(gbrs_samples))  # second-best array sample name
    second_r     <- numeric(length(gbrs_samples))    # similarity of second-best match

    for (j in seq_along(gbrs_samples)) {
        # Get column j: similarity of every array sample to this GBRS sample
        col <- sim[, j]

        # Order indices by similarity, descending (highest first), NAs last
        ord <- order(col, decreasing = TRUE, na.last = TRUE)

        # First element = index of array sample with highest similarity
        best_array[j] <- array_samples[ord[1]]
        best_r[j] <- col[ord[1]]  # the actual similarity value

        if (length(ord) >= 2) {
            # Second element = second-best matching array sample
            second_array[j] <- array_samples[ord[2]]
            second_r[j] <- col[ord[2]]
        } else {
            # Only one array sample (or none); no "second" exists
            second_array[j] <- NA_character_
            second_r[j] <- NA_real_
        }
    }

    # =========================================================================
    # STEP 4: Get the similarity of each GBRS sample to its LABELED (expected) array
    # =========================================================================
    # For each gbrs sample j: look up sim[labeled_array[j], j]
    # i.e., the similarity between gbrs j and the array sample we EXPECT it to match

    labeled_r <- vapply(seq_along(gbrs_samples), function(j) {
        a <- labeled_array[[j]]  # the array sample we expect to match
        if (is.na(a)) return(NA_real_)  # no expected match (unmapped)
        # sim[array_sample, gbrs_sample] = similarity between that pair
        sim[a, gbrs_samples[j]]
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
    flag <- rep("match", length(gbrs_samples))

    # unmapped: no expected array sample for this gbrs sample (can't check)
    flag[is.na(labeled_array)] <- "unmapped"

    # mismatch_confident: best != labeled AND we're confident it's a real mismatch
    #   - delta_best_labeled >= threshold: best is clearly better than labeled
    #   - delta_best_second >= threshold: best is clearly better than second (not ambiguous)
    flag[!is.na(labeled_array) & (best_array != labeled_array) &
             (delta_best_labeled >= min_delta_best_labeled) &
             (delta_best_second  >= min_delta_best_second)] <- "mismatch_confident"

    # mismatch_ambiguous: best != labeled but second-best is close to best
    #   - We suspect a mismatch, but could be noise (best and second nearly tied)
    flag[!is.na(labeled_array) & (best_array != labeled_array) &
             (delta_best_labeled >= min_delta_best_labeled) &
             (delta_best_second  <  min_delta_best_second)] <- "mismatch_ambiguous"

    # match_ambiguous: best == labeled (correct!) but second is close to best
    #   - Probably correct, but if labels were wrong we couldn't tell (ambiguous)
    flag[!is.na(labeled_array) & (best_array == labeled_array) &
             (delta_best_second  <  min_delta_best_second)] <- "match_ambiguous"

    # =========================================================================
    # STEP 7: Build the output data.frame
    # =========================================================================
    out <- data.frame(
        gbrs_sample = gbrs_samples,
        labeled_array_sample = unname(labeled_array),  # expected array match
        labeled_r = labeled_r,                         # similarity to expected
        best_array_sample = best_array,                # actual best array match
        best_r = best_r,                               # similarity to best
        second_array_sample = second_array,            # second-best array match
        second_r = second_r,                           # similarity to second-best
        delta_best_second = delta_best_second,         # best - second
        delta_best_labeled = delta_best_labeled,       # best - labeled
        flag = flag,
        stringsAsFactors = FALSE
    )

    # =========================================================================
    # STEP 8: Identify reciprocal swap pairs (A↔B swapped)
    # =========================================================================
    # A reciprocal swap: gbrs A's best match is array B, and gbrs B's best match is array A
    # i.e., bestA[A]==labeledB and bestA[B]==labeledA
    # We tag both A and B with the same reciprocal_pair ID (e.g., "swap:A<->B")

    reciprocal <- rep(NA_character_, nrow(out))  # will hold "swap:X<->Y" or NA

    # Only consider rows that are mismatches (confident or ambiguous) and have a labeled sample
    idx_mis <- which(out$flag %in% c("mismatch_confident", "mismatch_ambiguous") &
                         !is.na(out$labeled_array_sample))

    # Named vectors for quick lookup:
    # bestA[gbrs] = array sample that best matches this gbrs
    # labA[gbrs]  = array sample we EXPECT to match this gbrs (labeled)
    bestA <- setNames(out$best_array_sample, out$gbrs_sample)
    labA  <- setNames(out$labeled_array_sample, out$gbrs_sample)

    for (i in idx_mis) {
        g1 <- out$gbrs_sample[i]  # first gbrs sample in potential swap pair

        # For a swap, we need g2 such that:
        #   bestA[g1] == labA[g2]  (g1's best match is g2's expected array)
        #   bestA[g2] == labA[g1]  (g2's best match is g1's expected array)
        # So: find all g2 whose LABELED array equals g1's BEST match
        candidates <- names(labA)[labA == bestA[[g1]]]

        for (g2 in candidates) {
            # Check the reverse: does g2's best match equal g1's labeled array?
            if (!is.na(bestA[[g2]]) && bestA[[g2]] == labA[[g1]]) {
                # Yes! g1 and g2 form a reciprocal swap pair
                rid <- paste0("swap:", paste(sort(c(g1, g2)), collapse = "<->"))
                # Assign same ID to both rows (sort ensures consistent ordering)
                reciprocal[match(g1, out$gbrs_sample)] <- rid
                reciprocal[match(g2, out$gbrs_sample)] <- rid
            }
        }
    }

    out$reciprocal_pair <- reciprocal
    out
}
