#' Detect sample swaps with robust, evidence-layered logic (column-centric view)
#'
#' @description
#' \code{genoprobs_detect_sample_swaps()} is a 'max robustness' reference
#' implementation for sample identity QC between two sample sets. It interprets
#' the input similarity matrix in a **column‑centric, query‑vs‑reference**
#' fashion:
#'
#' - **Rows** are the *reference* sample set (typically the more trusted set).
#' - **Columns** are the *query* sample set whose labels you are checking.
#'
#' For each **column sample** (query), it asks:
#'
#' > 'Which reference sample (row) is my best match, how different is that from
#' > my expected label, and is the evidence strong enough to confidently call
#' > a match, mismatch, or mixture?'
#'
#' The function is intentionally conservative and explicit:
#' \itemize{
#'   \item It keeps the familiar local evidence (best vs. second; best vs. expected).
#'   \item It adds empirical calibration (z-score and empirical tail probability).
#'   \item It adds one-to-one global assignment (greedy by default; Hungarian if \pkg{clue} is available).
#'   \item It marks low-information columns as \code{'no_call'} instead of forcing decisions.
#'   \item It annotates swap groups beyond simple reciprocal pairs.
#' }
#'
#' Why this shape? In practical bioinformatics QC pipelines, false certainty is
#' usually more harmful than an explicit inconclusive call. This function is
#' designed to make evidence and uncertainty visible, not hidden.
#'
#' @section Inputs and orientation:
#' You may provide either:
#' \itemize{
#'   \item a precomputed similarity matrix \code{x} with row/column sample
#'     names, interpreted as
#'     \itemize{
#'       \item \code{rownames(x)} = reference samples
#'       \item \code{colnames(x)} = query samples (one output row per column)
#'     }
#'   \item or two \code{calc_genoprob} objects (\code{x} and \code{genoprobs_2}),
#'     in which case similarity is computed internally via
#'     \code{\link{genoprobs_compute_similarity}} with:
#'     \itemize{
#'       \item \code{x}     → reference genoprobs (rows of the similarity matrix)
#'       \item \code{genoprobs_2} → query genoprobs (columns of the matrix)
#'     }
#' }
#'
#' In all cases, similarity must satisfy 'larger means more similar'
#' (e.g. Pearson/cosine), and the resulting matrix is treated as
#' \code{rows = reference}, \code{cols = query}.
#'
#' @section Decision strategy (high level):
#' For each **column sample** (i.e. each query sample) the function:
#' \enumerate{
#'   \item Finds the best and second-best **row** (reference) match.
#'   \item Compares the best match to the expected label (from \code{sample_map}
#'         or same‑ID default).
#'   \item Estimates whether the best is exceptional relative to that column's
#'         background (z‑score, empirical p‑value).
#'   \item Applies confidence-aware flags; may emit \code{'no_call'} instead of
#'         forcing a decision on weak evidence.
#'   \item Optionally compares the local best with a global one-to-one
#'         assignment (so that each reference is used at most once).
#' }
#'
#' @param x A similarity matrix (\code{rows = row samples}, \code{cols = col samples})
#'   or the first \code{calc_genoprob} object.
#' @param genoprobs_2 Optional second \code{calc_genoprob} object when \code{x} is a
#'   \code{calc_genoprob}.
#' @param sample_map Optional data.frame with columns
#'   \code{c('col_sample','row_sample')} describing expected pairings.
#'   If \code{NULL}, same-name mapping is used where available.
#' @param metric Similarity metric for \code{calc_genoprob} inputs:
#'   \code{'pearson'} or \code{'cosine'}.
#' @param assignment_method One-to-one assignment strategy:
#'   \code{'greedy'} (always available) or \code{'hungarian'}.
#'   If \code{'hungarian'} is requested but package \pkg{clue} is unavailable,
#'   the function falls back to \code{'greedy'} and emits a warning.
#' @param min_delta_best_labeled Minimum \code{best - labeled} to call a mismatch.
#' @param min_delta_best_second Minimum \code{best - second} for confident calls.
#' @param min_finite_per_col Minimum number of finite similarities required to
#'   consider a column evaluable; otherwise \code{'no_call'}.
#' @param min_z_best Minimum z-score of best-vs-background for confidence.
#' @param max_empirical_p_best Maximum empirical tail probability for confidence.
#' @param near_tie_gap Gap to define 'near ties' to the best value, used for
#'   mixture/ambiguity diagnostics.
#' @param allow_no_call If \code{TRUE}, low-information columns become \code{'no_call'}.
#'   If \code{FALSE}, they are labeled \code{'unmapped'} for compatibility with
#'   older workflows that do not expect \code{'no_call'}.
#' @param ... Reserved for forward compatibility.
#'
#' @return data.frame with one row per column sample and columns:
#' \describe{
#'   \item{col_sample}{Column sample ID.}
#'   \item{expected_row_sample}{Expected row sample ID (or \code{NA}).}
#'   \item{labeled_r}{Similarity to expected row sample.}
#'   \item{best_row_sample}{Best local row match (column-wise max).}
#'   \item{best_r}{Similarity to best local row match.}
#'   \item{second_row_sample}{Second-best local row match.}
#'   \item{second_r}{Similarity to second-best local row match.}
#'   \item{delta_best_second}{\code{best_r - second_r}.}
#'   \item{delta_best_labeled}{\code{best_r - labeled_r}.}
#'   \item{n_finite}{Count of finite similarities in the column.}
#'   \item{n_at_best}{How many rows tie for best similarity.}
#'   \item{n_near_best}{How many rows are within \code{near_tie_gap} of best.}
#'   \item{z_best}{Column-wise z-score of best relative to finite background.}
#'   \item{p_empirical_best}{Empirical upper-tail probability of best among
#'     finite row similarities in that column.}
#'   \item{global_row_sample}{Row assigned by global one-to-one assignment.}
#'   \item{assignment_consistent}{Whether global assignment equals local best.}
#'   \item{flag}{One of:
#'     \code{'match_confident'},
#'     \code{'match_ambiguous'},
#'     \code{'mismatch_confident'},
#'     \code{'mismatch_ambiguous'},
#'     \code{'possible_mixture'},
#'     \code{'no_call'},
#'     \code{'unmapped'}.}
#'   \item{swap_group_id}{ID of connected mismatch group (if any).}
#'   \item{swap_group_type}{\code{'pair'}, \code{'cycle'}, or \code{'complex'}.}
#' }
#'
#' @examples
#' \dontrun{
#' # Similarity matrix input
#' sim_res <- genoprobs_compute_similarity(gp_a, gp_b, metric = 'pearson')
#' out <- genoprobs_detect_sample_swaps(
#'   sim_res$sim,
#'   sample_map = sample_map,
#'   assignment_method = 'hungarian'
#' )
#'
#' # calc_genoprob input (computes similarity first)
#' out2 <- genoprobs_detect_sample_swaps(
#'   gp_a,
#'   genoprobs_2 = gp_b,
#'   sample_map = sample_map,
#'   metric = 'cosine'
#' )
#' }
#' @export
genoprobs_detect_sample_swaps <- function(
    x,
    genoprobs_2 = NULL,
    sample_map = NULL,
    metric = c('pearson', 'cosine'),
    assignment_method = c('greedy', 'hungarian'),
    min_delta_best_labeled = 0.05,
    min_delta_best_second = 0.02,
    min_finite_per_col = 25L,
    min_z_best = 3,
    max_empirical_p_best = 0.05,
    near_tie_gap = 0.01,
    allow_no_call = TRUE,
    ...
) {
    metric <- match.arg(metric)
    assignment_method <- match.arg(assignment_method)

    # -------------------------------------------------------------------------
    # STEP 1: Build or validate the similarity matrix.
    #
    # WHY: downstream logic assumes a row-by-column matrix where larger means
    # "more likely same sample". We normalize at the edge so internals can remain
    # simple and explicit.
    # -------------------------------------------------------------------------
    sim <- NULL

    if (is.matrix(x)) {
        sim <- x
    } else if (inherits(x, 'calc_genoprob')) {
        if (is.null(genoprobs_2)) {
            stop('When x is calc_genoprob, genoprobs_2 is required.')
        }
        sim <- genoprobs_compute_similarity(x, genoprobs_2, metric = metric)$sim
    } else {
        stop(
            'x must be either a similarity matrix or calc_genoprob.\n',
            'Received class: ', paste(class(x), collapse = '/')
        )
    }

    if (is.null(rownames(sim)) || is.null(colnames(sim))) {
        stop('sim must have rownames and colnames as sample IDs.')
    }
    if (anyDuplicated(rownames(sim))) {
        stop('rownames(sim) contain duplicates.')
    }
    if (anyDuplicated(colnames(sim))) {
        stop('colnames(sim) contain duplicates.')
    }

    row_samples <- rownames(sim)
    col_samples <- colnames(sim)
    n_col <- ncol(sim)

    # -------------------------------------------------------------------------
    # STEP 2: Build expected pairing map (col -> row).
    #
    # WHY: mismatch detection only has meaning relative to an expected label.
    # If no map is supplied, we default to same-name expectation where possible.
    # -------------------------------------------------------------------------
    if (is.null(sample_map)) {
        shared <- intersect(col_samples, row_samples)
        sample_map <- data.frame(
            col_sample = shared,
            row_sample = shared,
            stringsAsFactors = FALSE
        )
    } else {
        if (!all(c('col_sample', 'row_sample') %in% names(sample_map))) {
            stop('sample_map must contain columns \'col_sample\' and \'row_sample\'.')
        }
        sample_map <- sample_map[, c('col_sample', 'row_sample'), drop = FALSE]
        if (anyDuplicated(sample_map$col_sample)) {
            stop('sample_map has duplicated col_sample entries.')
        }
        if (anyDuplicated(sample_map$row_sample)) {
            warning(
                'sample_map has duplicated row_sample values. ',
                'One-to-one assignment and group labeling may be less interpretable.'
            )
        }
        sample_map <- sample_map[
            sample_map$col_sample %in% col_samples &
                sample_map$row_sample %in% row_samples,
            ,
            drop = FALSE
        ]
    }

    expected_row <- setNames(sample_map$row_sample, sample_map$col_sample)[col_samples]

    # -------------------------------------------------------------------------
    # STEP 3: Per-column local evidence.
    #
    # WHY: this preserves interpretability and backwards-compatibility with the
    # common "best/second/labeled" identity-QC pattern.
    # -------------------------------------------------------------------------
    n_finite <- colSums(is.finite(sim))
    col_has_value <- n_finite > 0L
    evaluable <- col_has_value & (n_finite >= as.integer(min_finite_per_col))

    best_idx <- rep(NA_integer_, n_col)
    best_r <- rep(NA_real_, n_col)
    second_idx <- rep(NA_integer_, n_col)
    second_r <- rep(NA_real_, n_col)
    n_at_best <- rep(NA_integer_, n_col)
    n_near_best <- rep(NA_integer_, n_col)
    z_best <- rep(NA_real_, n_col)
    p_empirical_best <- rep(NA_real_, n_col)

    for (j in seq_len(n_col)) {
        colv <- sim[, j]
        finite <- is.finite(colv)
        if (!any(finite)) next

        # Best local row for this column.
        finite_idx <- which(finite)
        local_vals <- colv[finite_idx]
        k_best <- which.max(local_vals)
        i_best <- finite_idx[k_best]
        best_idx[j] <- i_best
        best_r[j] <- colv[i_best]

        # Second best among finite values.
        if (length(finite_idx) >= 2L) {
            local_vals2 <- local_vals
            local_vals2[k_best] <- -Inf
            k_second <- which.max(local_vals2)
            i_second <- finite_idx[k_second]
            if (is.finite(colv[i_second])) {
                second_idx[j] <- i_second
                second_r[j] <- colv[i_second]
            }
        }

        # Tie/near-tie diagnostics help explain ambiguity.
        mx <- max(local_vals)
        n_at_best[j] <- sum(local_vals == mx)
        n_near_best[j] <- sum((mx - local_vals) <= near_tie_gap)

        # Empirical calibration: how exceptional is the best score within this column?
        #
        # WHY: fixed deltas alone can be brittle when score scales vary by sample,
        # platform, or depth. A local background gives a robust 'is this standout?'
        # signal without imposing a strict parametric model.
        mu <- mean(local_vals)
        sdv <- stats::sd(local_vals)
        if (is.finite(sdv) && sdv > 0) {
            z_best[j] <- (mx - mu) / sdv
        }
        p_empirical_best[j] <- (sum(local_vals >= mx) + 1) / (length(local_vals) + 1)
    }

    best_row_sample <- ifelse(is.na(best_idx), NA_character_, row_samples[best_idx])
    second_row_sample <- ifelse(is.na(second_idx), NA_character_, row_samples[second_idx])

    labeled_r <- vapply(seq_along(col_samples), function(j) {
        r <- expected_row[[j]]
        if (is.na(r)) {
            return(NA_real_)
        }
        sim[r, col_samples[j]]
    }, numeric(1))

    delta_best_second <- best_r - second_r
    delta_best_labeled <- best_r - labeled_r

    # -------------------------------------------------------------------------
    # STEP 4: Global one-to-one assignment (optional but recommended).
    #
    # WHY: local maxima are easy to interpret but can create many-to-one mappings.
    # A one-to-one assignment helps expose plate shifts/cycles and reduces local
    # greediness artifacts.
    # -------------------------------------------------------------------------
    global_row_for_col <- .swapmax_assign_rows(
        sim = sim,
        method = assignment_method
    )
    global_row_sample <- ifelse(
        is.na(global_row_for_col),
        NA_character_,
        row_samples[global_row_for_col]
    )
    assignment_consistent <- (global_row_sample == best_row_sample)
    assignment_consistent[is.na(global_row_sample) | is.na(best_row_sample)] <- NA

    # -------------------------------------------------------------------------
    # STEP 5: Confidence-aware flagging.
    #
    # WHY: explicit uncertainty categories are safer than binary forced calls.
    # -------------------------------------------------------------------------
    flag <- rep('match_ambiguous', n_col)

    # Missing mapping stays distinct from low-information no-call.
    flag[is.na(expected_row)] <- 'unmapped'

    if (allow_no_call) {
        flag[!evaluable] <- 'no_call'
    } else {
        flag[!evaluable] <- 'unmapped'
    }

    mapped_eval <- !is.na(expected_row) & evaluable

    confident <- mapped_eval &
        !is.na(delta_best_second) &
        !is.na(z_best) &
        !is.na(p_empirical_best) &
        (delta_best_second >= min_delta_best_second) &
        (z_best >= min_z_best) &
        (p_empirical_best <= max_empirical_p_best)

    is_match <- mapped_eval & (best_row_sample == expected_row)
    is_mismatch <- mapped_eval & (best_row_sample != expected_row)

    flag[is_match & confident] <- 'match_confident'
    flag[is_match & !confident] <- 'match_ambiguous'

    strong_disagreement <- is_mismatch &
        !is.na(delta_best_labeled) &
        (delta_best_labeled >= min_delta_best_labeled)

    flag[strong_disagreement & confident] <- 'mismatch_confident'
    flag[strong_disagreement & !confident] <- 'mismatch_ambiguous'

    # Mixture heuristic:
    # multiple near-best candidates + weak confidence + mismatch-like behavior.
    possible_mixture <- is_mismatch &
        (n_near_best >= 3L) &
        (
            is.na(delta_best_second) |
                delta_best_second < min_delta_best_second |
                is.na(z_best) |
                z_best < min_z_best
        )
    flag[possible_mixture] <- 'possible_mixture'

    # -------------------------------------------------------------------------
    # STEP 6: Group labeling for mismatches.
    #
    # WHY: identity failures are often not isolated pairs; they can be shifts or
    # cycle-like permutations. Group IDs help users inspect these patterns.
    # -------------------------------------------------------------------------
    grp <- .swapmax_label_groups(
        col_samples = col_samples,
        expected_row = unname(expected_row),
        assigned_row = global_row_sample,
        mismatch_mask = flag %in% c('mismatch_confident', 'mismatch_ambiguous', 'possible_mixture')
    )

    out <- data.frame(
        col_sample = col_samples,
        expected_row_sample = unname(expected_row),
        labeled_r = labeled_r,
        best_row_sample = best_row_sample,
        best_r = best_r,
        second_row_sample = second_row_sample,
        second_r = second_r,
        delta_best_second = delta_best_second,
        delta_best_labeled = delta_best_labeled,
        n_finite = n_finite,
        n_at_best = n_at_best,
        n_near_best = n_near_best,
        z_best = z_best,
        p_empirical_best = p_empirical_best,
        global_row_sample = global_row_sample,
        assignment_consistent = assignment_consistent,
        flag = flag,
        swap_group_id = grp$id,
        swap_group_type = grp$type,
        stringsAsFactors = FALSE
    )

    out
}


#' Greedy/Hungarian row assignment for columns
#'
#' Internal helper for \code{genoprobs_detect_sample_swaps()}.
#' Returns integer vector of row indices for each column (or \code{NA}).
#'
#' @keywords internal
.swapmax_assign_rows <- function(sim, method = c('greedy', 'hungarian')) {
    method <- match.arg(method)
    nr <- nrow(sim)
    nc <- ncol(sim)
    out <- rep(NA_integer_, nc)

    if (method == 'hungarian' && !requireNamespace('clue', quietly = TRUE)) {
        warning('assignment_method=\'hungarian\' requested but package \'clue\' is unavailable; using greedy.')
        method <- 'greedy'
    }

    if (method == 'greedy') {
        # WHY greedy fallback: always available, no external dependency.
        # We rank all finite pairs by similarity and keep the best non-conflicting
        # row/column matches.
        idx <- which(is.finite(sim), arr.ind = TRUE)
        if (!nrow(idx)) {
            return(out)
        }
        vals <- sim[idx]
        ord <- order(vals, decreasing = TRUE)
        row_used <- rep(FALSE, nr)
        col_used <- rep(FALSE, nc)

        for (k in ord) {
            i <- idx[k, 1]
            j <- idx[k, 2]
            if (!row_used[i] && !col_used[j]) {
                out[j] <- i
                row_used[i] <- TRUE
                col_used[j] <- TRUE
            }
        }
        return(out)
    }

    # Hungarian assignment (requires clue). We maximize similarity.
    #
    # `clue::solve_LSAP` solves a square problem. For rectangular matrices we pad
    # with very poor dummy scores; this permits unmatched assignments when needed.
    fin <- is.finite(sim)
    if (!any(fin)) {
        return(out)
    }

    finite_vals <- sim[fin]
    floor_val <- min(finite_vals) - (abs(min(finite_vals)) + 1)
    m <- sim
    m[!fin] <- floor_val

    n <- max(nr, nc)
    sq <- matrix(floor_val, nrow = n, ncol = n)
    sq[seq_len(nr), seq_len(nc)] <- m

    # solve_LSAP returns selected columns for each row.
    assign_col_for_row <- clue::solve_LSAP(sq, maximum = TRUE)
    assign_col_for_row <- as.integer(assign_col_for_row)

    # Invert row->col to col->row for real columns/rows only.
    for (i in seq_len(nr)) {
        j <- assign_col_for_row[i]
        if (j <= nc) out[j] <- i
    }

    # Drop assignments that correspond to non-finite original entries.
    assigned_cols <- which(!is.na(out))
    if (length(assigned_cols)) {
        assigned_rows <- out[assigned_cols]
        bad_assigned <- !is.finite(sim[cbind(assigned_rows, assigned_cols)])
        if (any(bad_assigned)) {
            out[assigned_cols[bad_assigned]] <- NA_integer_
        }
    }
    out
}


#' Label mismatch groups and broad topology type
#'
#' Internal helper for \code{genoprobs_detect_sample_swaps()}.
#'
#' @keywords internal
.swapmax_label_groups <- function(col_samples, expected_row, assigned_row, mismatch_mask) {
    n <- length(col_samples)
    id <- rep(NA_character_, n)
    type <- rep(NA_character_, n)

    # Map expected row -> column sample. If duplicated expectations exist, choose
    # first; this remains informative but less strictly interpretable.
    col_for_row <- setNames(col_samples, expected_row)

    # Directed edge: c1 -> c2 means c1 was assigned to row expected for c2.
    to <- rep(NA_integer_, n)
    for (i in which(mismatch_mask & !is.na(assigned_row))) {
        nxt <- col_for_row[[assigned_row[i]]]
        if (!is.null(nxt) && !is.na(nxt)) {
            j <- match(nxt, col_samples)
            if (!is.na(j) && j != i) to[i] <- j
        }
    }

    # Undirected adjacency for connected groups.
    adj <- vector('list', n)
    for (i in seq_len(n)) {
        if (!is.na(to[i])) {
            adj[[i]] <- unique(c(adj[[i]], to[i]))
            adj[[to[i]]] <- unique(c(adj[[to[i]]], i))
        }
    }

    seen <- rep(FALSE, n)
    comp_idx <- 0L
    candidates <- which(mismatch_mask & !is.na(to))
    for (start in candidates) {
        if (seen[start]) next
        comp_idx <- comp_idx + 1L
        q <- start
        members <- integer(0)
        seen[start] <- TRUE

        while (length(q)) {
            v <- q[1]
            q <- q[-1]
            members <- c(members, v)
            nbr <- adj[[v]]
            if (!length(nbr)) next
            for (u in nbr) {
                if (!seen[u]) {
                    seen[u] <- TRUE
                    q <- c(q, u)
                }
            }
        }

        members <- sort(unique(members))
        outdeg <- vapply(members, function(v) as.integer(!is.na(to[v])), integer(1))
        indeg <- vapply(members, function(v) sum(to[members] == v, na.rm = TRUE), integer(1))
        is_cycle <- all(outdeg == 1L) && all(indeg == 1L) && length(members) >= 3L
        is_pair <- length(members) == 2L && all(outdeg == 1L) && all(indeg == 1L)

        comp_type <- if (is_pair) {
            'pair'
        } else if (is_cycle) {
            'cycle'
        } else {
            'complex'
        }

        group_id <- paste0('group:', comp_idx)
        id[members] <- group_id
        type[members] <- comp_type
    }

    list(id = id, type = type)
}
