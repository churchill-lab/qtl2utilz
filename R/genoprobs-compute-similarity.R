#' Extract sample names from a calc_genoprob object
#' @keywords internal
get_gp_samples <- function(gp, chr = NULL) {
    if (is.null(chr)) chr <- names(gp)[1]
    dn <- dimnames(gp[[chr]])
    if (is.null(dn) || length(dn) < 1 || is.null(dn[[1]])) {
        stop('Could not find sample names in dimnames(genoprobs[[chr]])[[1]].')
    }
    dn[[1]]
}

#' Compute similarity matrix between two aligned calc_genoprob objects
#'
#' @section What this function does (big picture):
#'
#' Genoprobs are 3D arrays per chromosome: \code{dimnames[[1]]} = samples,
#' \code{dimnames[[2]]} = haplotypes (e.g. founders A–H), \code{dimnames[[3]]} =
#' markers. Each sample is thus a long vector of probabilities (one per
#' marker-haplotype combination). This function compares two such objects
#' (e.g. different platforms or builds) already aligned to the same grid (same
#' markers in the same order; e.g. via genoprobs_sync_markers). For every pair
#' (sample i from genoprobs_1, sample j from genoprobs_2) it computes one number: how
#' similar are those two probability vectors? Similarity is Pearson or cosine.
#' Result: matrix with rows = genoprobs_1 samples, cols = genoprobs_2 samples; cell [i,j] =
#' similarity of that pair. Downstream swap detection uses high similarity on
#' the diagonal (same mouse) and low off-diagonal.
#'
#' @section Sample order and sample names:
#'
#' Samples do **not** need to be in the same order in genoprobs_1 and genoprobs_2. The
#' function computes similarity for every (row of genoprobs_1, column of genoprobs_2) pair
#' and labels the output with each object's own sample names. Sample names are
#' **not** checked or required to match between the two objects; matching of
#' "which row sample is the same mouse as which column sample" is done
#' downstream (e.g. by genoprobs_detect_sample_swaps using a sample_map or
#' by inspecting the similarity matrix). Both objects must have the same
#' number of haplotypes and same number of markers per chromosome (enforced
#' internally); they may have different numbers of samples.
#'
#' @param genoprobs_1,genoprobs_2 calc_genoprob objects already aligned to the same map/grid
#'   (e.g. via genoprobs_sync_markers).
#' @param metric 'pearson' or 'cosine'
#' @param by_chr if TRUE, also return per-chromosome similarity matrices
#'
#' @return list(sim, sim_by_chr, n_features). sim has rownames = genoprobs_1 sample IDs,
#'   colnames = genoprobs_2 sample IDs; each cell is the similarity between that pair.
#' @export
genoprobs_compute_similarity <- function(genoprobs_1, genoprobs_2,
                                metric = c('pearson', 'cosine'),
                                by_chr = FALSE) {
    metric <- match.arg(metric)

    # =========================================================================
    # STEP 1: Validate shared chromosomes and get sample lists
    # =========================================================================
    # We need at least one chromosome in common so we have a comparable marker grid.
    # We do NOT require sample names to match between genoprobs_1 and genoprobs_2, or samples
    # to be in the same order; we just use each object's sample order and names
    # to label the output matrix.

    common_chr <- intersect(names(genoprobs_1), names(genoprobs_2))
    if (length(common_chr) == 0) {
        stop('No shared chromosomes between genoprobs_1 and genoprobs_2.')
    }

    # Get sample IDs and counts from each object (first shared chr is enough;
    # order/names can differ between genoprobs_1 and genoprobs_2).
    sa <- get_gp_samples(genoprobs_1, common_chr[1])
    sb <- get_gp_samples(genoprobs_2, common_chr[1])
    n_a <- length(sa)
    n_b <- length(sb)

    # =========================================================================
    # STEP 2: Initialize accumulators for whole-genome similarity
    # =========================================================================
    # Each sample is a long vector of probabilities (one per marker x haplotype).
    # Similarity between two samples = correlation or cosine of those vectors.
    # We accumulate sum(x), sum(x^2), sum(y), sum(y^2), and sum(x*y) over all chr
    # so we can compute Pearson or cosine at the end. Result: sim[i,j] =
    # similarity of genoprobs_1 sample i with genoprobs_2 sample j. Rows = genoprobs_1 (sa), cols = genoprobs_2 (sb).

    sumxy <- matrix(0, nrow = n_a, ncol = n_b, dimnames = list(sa, sb))
    sumx  <- numeric(n_a)
    sumx2 <- numeric(n_a)
    sumy  <- numeric(n_b)
    sumy2 <- numeric(n_b)
    p_tot <- 0L  # total features (markers * haplotypes) across all chr

    # Optional: per-chromosome similarity matrices (same layout as sim)
    sim_by_chr <- if (by_chr) setNames(vector('list', length(common_chr)), common_chr) else NULL

    # =========================================================================
    # STEP 3: Loop over chromosomes; accumulate cross-products and norms
    # =========================================================================
    # Layout: A and B are [samples x haplotypes x markers]. Aligned objects have
    # the same markers (dim 3) and haplotypes (dim 2) in the same order, so
    # column k in the flattened vector is the same marker-haplotype in both.
    for (ch in common_chr) {
        A <- genoprobs_1[[ch]]  # [samples x haplotypes x markers]
        B <- genoprobs_2[[ch]]
        dna <- dimnames(A)
        dnb <- dimnames(B)

        # Flattened length = haplotypes * markers. We need same length in both.
        if (dim(A)[2] != dim(B)[2] || dim(A)[3] != dim(B)[3]) {
            stop('Mismatch in haplotypes (dim 2) or markers (dim 3) for chr ', ch,
                 '. Ensure both objects are aligned (same founders, same marker grid).')
        }
        # Same identity and order of haplotypes and markers (fail fast if synced grid is violated).
        if (!is.null(dna[[2]]) && !is.null(dnb[[2]]) && !identical(dna[[2]], dnb[[2]])) {
            stop('Haplotype dimnames differ on chr ', ch, '.')
        }
        if (!is.null(dna[[3]]) && !is.null(dnb[[3]]) && !identical(dna[[3]], dnb[[3]])) {
            stop('Marker dimnames differ on chr ', ch, '.')
        }

        # Sample count and order must be consistent within each object across chromosomes.
        if (dim(A)[1] != n_a) {
            stop('Sample count differs on chr ', ch, ' in genoprobs_1.')
        }
        if (!is.null(dna[[1]]) && !identical(dna[[1]], sa)) {
            stop('Sample order/names differ on chr ', ch, ' in genoprobs_1.')
        }
        if (dim(B)[1] != n_b) {
            stop('Sample count differs on chr ', ch, ' in genoprobs_2.')
        }
        if (!is.null(dnb[[1]]) && !identical(dnb[[1]], sb)) {
            stop('Sample order/names differ on chr ', ch, ' in genoprobs_2.')
        }

        n_haplo  <- dim(A)[2]
        n_markers <- dim(A)[3]
        p_ch     <- n_haplo * n_markers  # features for this chromosome
        p_tot    <- p_tot + p_ch

        # Reshape to [samples x features]: each row = one sample's probability vector.
        # aperm(., c(1,3,2)) -> [samples, markers, haplotypes]; then flatten to [samples, markers*haplotypes]
        Xa <- matrix(aperm(A, c(1, 3, 2)), nrow = n_a)
        Xb <- matrix(aperm(B, c(1, 3, 2)), nrow = n_b)

        # Accumulate for whole-genome Pearson/cosine:
        # sumx/sumx2 = row sums and sum of squares for genoprobs_1; sumy/sumy2 for genoprobs_2
        # sumxy[i,j] = dot product of genoprobs_1 sample i with genoprobs_2 sample j
        sumx  <- sumx  + rowSums(Xa)
        sumx2 <- sumx2 + rowSums(Xa * Xa)
        sumy  <- sumy  + rowSums(Xb)
        sumy2 <- sumy2 + rowSums(Xb * Xb)
        sumxy <- sumxy + Xa %*% t(Xb)

        # Optionally compute and store per-chromosome similarity matrix
        if (by_chr) {
            if (metric == 'cosine') {
                # Cosine = dot(a,b) / (||a|| * ||b||)
                denom <- sqrt(outer(rowSums(Xa * Xa), rowSums(Xb * Xb), '*'))
                denom[denom == 0] <- NA_real_
                sim_by_chr[[ch]] <- (Xa %*% t(Xb)) / denom
            } else {
                # Pearson for this chr: (sxy - p_ch*mean(x)*mean(y)) / (sd(x)*sd(y)); p_ch = number of features
                sx  <- rowSums(Xa)
                sx2 <- rowSums(Xa * Xa)
                sy  <- rowSums(Xb)
                sy2 <- rowSums(Xb * Xb)
                sxy <- Xa %*% t(Xb)
                denom <- sqrt(outer(sx2 - (sx * sx) / p_ch, sy2 - (sy * sy) / p_ch, '*'))
                denom[denom == 0] <- NA_real_
                sim_by_chr[[ch]] <- (sxy - outer(sx, sy, '*') / p_ch) / denom
            }
            rownames(sim_by_chr[[ch]]) <- sa
            colnames(sim_by_chr[[ch]]) <- sb
        }

        rm(Xa, Xb, A, B)
    }

    # =========================================================================
    # STEP 4: Compute whole-genome similarity matrix from accumulators
    # =========================================================================
    # Same formula as per-chr but using totals (p_tot, sumx, sumy, sumx2, sumy2, sumxy).

    if (metric == 'cosine') {
        denom <- sqrt(outer(sumx2, sumy2, '*'))
        denom[denom == 0] <- NA_real_
        sim <- sumxy / denom
    } else {
        # Pearson: covariance / (sd_a * sd_b); variance = sumx2 - (sumx^2)/p_tot
        denom <- sqrt(outer(sumx2 - (sumx * sumx) / p_tot,
                          sumy2 - (sumy * sumy) / p_tot, '*'))
        denom[denom == 0] <- NA_real_
        sim <- (sumxy - outer(sumx, sumy, '*') / p_tot) / denom
    }

    list(sim = sim, sim_by_chr = sim_by_chr, n_features = p_tot)
}
