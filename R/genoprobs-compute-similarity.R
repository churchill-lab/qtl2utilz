#' Extract sample names from a calc_genoprob object
#' @keywords internal
get_gp_samples <- function(gp, chr = NULL) {
    if (is.null(chr)) chr <- names(gp)[1]
    dn <- dimnames(gp[[chr]])
    if (is.null(dn) || length(dn) < 1 || is.null(dn[[1]])) {
        stop("Could not find sample names in dimnames(genoprobs[[chr]])[[1]].")
    }
    dn[[1]]
}

#' Compute similarity matrix between two aligned calc_genoprob objects
#'
#' @param gp_a,gp_b calc_genoprob objects already aligned to the same map/grid.
#' @param metric "pearson" or "cosine"
#' @param by_chr if TRUE, also return per-chr similarity matrices
#'
#' @return list(sim, sim_by_chr, n_features)
#' @export
genoprob_compute_similarity <- function(gp_a, gp_b,
                                metric = c("pearson", "cosine"),
                                by_chr = FALSE) {
    metric <- match.arg(metric)

    common_chr <- intersect(names(gp_a), names(gp_b))
    if (length(common_chr) == 0) stop("No shared chromosomes between gp_a and gp_b.")

    sa <- get_gp_samples(gp_a, common_chr[1])
    sb <- get_gp_samples(gp_b, common_chr[1])
    n_a <- length(sa)
    n_b <- length(sb)

    # accumulators for whole-genome pearson/cosine
    sumxy <- matrix(0, nrow = n_a, ncol = n_b, dimnames = list(sa, sb))
    sumx  <- numeric(n_a); sumx2 <- numeric(n_a)
    sumy  <- numeric(n_b); sumy2 <- numeric(n_b)
    p_tot <- 0L

    sim_by_chr <- if (by_chr) setNames(vector("list", length(common_chr)), common_chr) else NULL

    for (ch in common_chr) {
        A <- gp_a[[ch]] # [ind x pos x state]
        B <- gp_b[[ch]]
        if (any(dim(A)[c(1,3)] != dim(B)[c(1,3)])) {
            stop("Mismatch in n_ind or n_state for chr ", ch,
                 ". Ensure both objects are comparable (same samples? same states?).")
        }
        n_pos   <- dim(A)[2]
        n_state <- dim(A)[3]
        p_ch <- n_pos * n_state
        p_tot <- p_tot + p_ch

        # reshape to [ind x (pos*state)]
        Xa <- matrix(aperm(A, c(1, 3, 2)), nrow = n_a)
        Xb <- matrix(aperm(B, c(1, 3, 2)), nrow = n_b)

        sumx  <- sumx  + rowSums(Xa)
        sumx2 <- sumx2 + rowSums(Xa * Xa)
        sumy  <- sumy  + rowSums(Xb)
        sumy2 <- sumy2 + rowSums(Xb * Xb)
        sumxy <- sumxy + Xa %*% t(Xb)

        if (by_chr) {
            if (metric == "cosine") {
                denom <- sqrt(outer(rowSums(Xa * Xa), rowSums(Xb * Xb), "*"))
                sim_by_chr[[ch]] <- (Xa %*% t(Xb)) / denom
            } else {
                sx  <- rowSums(Xa); sx2 <- rowSums(Xa * Xa)
                sy  <- rowSums(Xb); sy2 <- rowSums(Xb * Xb)
                sxy <- Xa %*% t(Xb)
                denom <- sqrt(outer(sx2 - (sx * sx) / p_ch, sy2 - (sy * sy) / p_ch, "*"))
                sim_by_chr[[ch]] <- (sxy - outer(sx, sy, "*") / p_ch) / denom
            }
            rownames(sim_by_chr[[ch]]) <- sa
            colnames(sim_by_chr[[ch]]) <- sb
        }

        rm(Xa, Xb, A, B)
    }

    if (metric == "cosine") {
        sim <- sumxy / sqrt(outer(sumx2, sumy2, "*"))
    } else {
        denom <- sqrt(outer(sumx2 - (sumx * sumx) / p_tot,
                            sumy2 - (sumy * sumy) / p_tot, "*"))
        sim <- (sumxy - outer(sumx, sumy, "*") / p_tot) / denom
    }

    list(sim = sim, sim_by_chr = sim_by_chr, n_features = p_tot)
}


