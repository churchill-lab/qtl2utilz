################################################################################
# Interpolate genotype probabilities from one marker set to another.
# Daniel Gatti, dan.gatti@jax.org, 2022-09-08
################################################################################

#' Interpolate genoprobs for one chromosome to new marker positions
#'
#' Linearly interpolate founder probabilities from the marker positions in
#' \code{pr1} to the positions in \code{mkr2}. Positions must be in base pairs.
#' Used internally by \code{genoprobs_interpolate}; not exported.
#'
#' @param pr1 3D numeric array: samples x founders x markers (dimnames for
#'   third dimension must match \code{names(mkr1)}).
#' @param mkr1 Named numeric vector of marker positions (bp) for \code{pr1}.
#' @param mkr2 Named numeric vector of target marker positions (bp) to
#'   interpolate to.
#'
#' @return 3D array with same samples and founders as \code{pr1}, third
#'   dimension = length(mkr2) with \code{dimnames[[3]] = names(mkr2)}.
#'
#' @keywords internal
interpolate_chromosome = function(pr1, mkr1, mkr2) {
    stopifnot(identical(dimnames(pr1)[[3]], names(mkr1)))
    stopifnot(min(mkr1) > 200)
    stopifnot(min(mkr2) > 200)

    # Integer positions for findInterval (bp; truncation is fine).
    x1 <- as.integer(as.numeric(mkr1))
    x2 <- as.integer(as.numeric(mkr2))
    n1 <- length(x1)

    # For each target position x2, get index of nearest source marker to the
    # left (proximal) and to the right (distal) for linear interpolation.
    proximal <- findInterval(x2, x1, left.open = TRUE)
    distal   <- findInterval(x2, x1, left.open = FALSE) + 1L

    # Clamp targets outside source range to first/last marker so we never
    # index out of bounds; mimics legacy NA handling.
    below <- proximal == 0L
    proximal[below] <- 1L
    distal[below]   <- 1L
    distal[distal > n1] <- n1
    proximal[proximal > n1] <- n1

    # Extract prob arrays at proximal and distal marker indices.
    prox_pr <- pr1[, , proximal, drop = FALSE]
    dist_pr <- pr1[, , distal,   drop = FALSE]

    # Linear interpolation weight: 0 at proximal position, 1 at distal.
    w <- (x2 - x1[proximal]) / (x1[distal] - x1[proximal])

    # Avoid NaN/Inf when proximal == distal (zero denominator) or invalid positions.
    w[is.nan(w)]      <- 0
    w[is.infinite(w)] <- 1
    w[!is.finite(w)]  <- 0

    # Broadcast weight to (samples x founders x target_markers) for vectorized calc.
    w_arr <- array(rep(w, each = nrow(pr1) * ncol(pr1)),
                   dim = c(nrow(pr1), ncol(pr1), length(x2)),
                   dimnames = list(rownames(pr1), colnames(pr1), names(mkr2)))

    out <- prox_pr + w_arr * (dist_pr - prox_pr)
    dimnames(out)[[3]] <- names(mkr2)
    out
}


#' Interpolate genoprobs to a new marker map
#'
#' Given qtl2-style genotype probabilities and their marker map, interpolate
#' the probabilities onto a second marker map (e.g. from MUGA or a different
#' array). Linear interpolation in base-pair position is used per chromosome.
#' If marker positions are in Mb (max < 200), they are converted to bp
#' internally.
#'
#' @param probs1 qtl2-style genoprobs: named list of 3D arrays (samples x
#'   founders x markers), one element per chromosome.
#' @param markers1 qtl2-style marker map for \code{probs1}: named list of
#'   named numeric vectors (marker positions), names = chromosome and marker.
#' @param markers2 qtl2-style marker map to interpolate onto: same list
#'   structure; \code{probs1} will be interpolated to these positions.
#'
#' @return A genoprobs object with the same structure and attributes as
#'   \code{probs1}, but with marker sets (and third dimension) from
#'   \code{markers2} for each chromosome.
#'
#' @details
#' \code{probs1} and \code{markers1} must have the same length (chromosomes).
#' For each chromosome, \code{interpolate_chromosome} is called. Attributes
#' (\code{crosstype}, \code{is_x_chr}, \code{alleles}, \code{alleleprobs}) and
#' class are copied from \code{probs1}.
#'
#' @export
genoprobs_interpolate = function(probs1, markers1, markers2) {

    if(length(probs1) != length(markers1)) {
        stop('Probs1 and markers1 must be the same length.')
    }

    all_chr <- names(probs1)
    missing_chr <- setdiff(all_chr, names(markers2))
    if(length(missing_chr) > 0L) {
        stop("markers2 is missing chromosomes: ", paste(missing_chr, collapse = ", "))
    }

    # If positions look like Mb (max < 200), convert to bp; interpolate_chromosome
    # expects bp. Chr 1 is < 200 Mb so this heuristic is safe for mouse.
    if(max(unlist(markers1)) < 200) {
        markers1 = lapply(markers1, '*', 1e6)
    }
    if(max(unlist(markers2)) < 200) {
        markers2 = lapply(markers2, '*', 1e6)
    }

    new_probs1        <- vector("list", length(probs1))
    names(new_probs1) <- names(probs1)

    # Interpolate each chromosome from markers1 positions to markers2 positions.
    for(chr in all_chr) {
        new_probs1[[chr]] = interpolate_chromosome(pr1  = probs1[[chr]],
                                                   mkr1 = markers1[[chr]],
                                                   mkr2 = markers2[[chr]])
    }

    # Preserve qtl2 attributes and class so result is a valid genoprobs object.
    attr(new_probs1, 'crosstype')    = attr(probs1, 'crosstype')
    attr(new_probs1, 'is_x_chr')     = attr(probs1, 'is_x_chr')
    attr(new_probs1, 'alleles')      = attr(probs1, 'alleles')
    attr(new_probs1, 'alleleprobs')  = attr(probs1, 'alleleprobs')
    class(new_probs1)                = class(probs1)

    new_probs1
}

