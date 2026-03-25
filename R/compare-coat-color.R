################################################################################
# Test inferred genotype at coat color loci versus reported coat color for
# Diversity Outbred mice.
# Users must report the coat color as 'agouti', 'albino', or 'black'.
# Indeterminate coat colors may be coded as NA.
# We only use the allele probs at the agouti locus on Chr 2 and the tyrosinase
# locus on Chr 7. Coordinates in GRCm39.
# Albino is epistatic (i.e. masks) to black. So test the albino locus
# first and the agouti locus second.
#
# 2025-05-23
################################################################################



#' Call the most probable marginal diplotype at a single marker.
#'
#' Originally written by Daniel Gatti (dan.gatti@jax.org).
#'
#' @param pr A samples-by-founders matrix of allele probabilities at one marker.
#' @return data.frame with sample ID, two-letter diplotype, and founder probs.
#' @keywords internal
call_diplotypes <- function(pr) {
    # Round probabilities to recover approximate homozygous/heterozygous states.
    pr2 <- round(2 * pr)

    # Collect the founder letters contributing to each diplotype call.
    diplotype_index <- apply(pr2 > 0, 1, which)
    diplotype_index <- lapply(diplotype_index, function(founders) {
        sort(names(founders))
    })
    diplotype_size <- sapply(diplotype_index, length)

    diplotype <- setNames(rep("", length(diplotype_index)), names(diplotype_index))

    # Set homozygotes.
    homozygous_rows <- which(diplotype_size == 1)
    diplotype[homozygous_rows] <- paste0(
        diplotype_index[homozygous_rows],
        diplotype_index[homozygous_rows]
    )

    # Set heterozygotes.
    heterozygous_rows <- which(diplotype_size == 2)
    diplotype[heterozygous_rows] <- sapply(
        diplotype_index[heterozygous_rows],
        paste0,
        collapse = ""
    )

    # Rows that do not sum to 2 are ambiguous after rounding. Fall back to the
    # top two founder probabilities so every sample still gets a best call.
    row_sum <- rowSums(pr2)
    ambiguous_rows <- which(row_sum != 2)

    for (row_idx in ambiguous_rows) {
        top_founders <- sort(pr[row_idx, ])[7:8]
        diplotype[row_idx] <- paste0(sort(names(top_founders)), collapse = "")
    }

    data.frame(id = names(diplotype), dt = diplotype, round(pr, digits = 3))
}


#' Compare genotype at coat color loci with reported coat color
#'
#' Infers coat color from allele probabilities at the albino (Tyrosinase,
#' Chr 7) and black (Agouti, Chr 2) loci, then compares with reported
#' phenotype.
#'
#' Originally written by Daniel Gatti (dan.gatti@jax.org).
#'
#' @param pheno data.frame with columns \code{id} (mouse ID) and \code{coat}
#'   (reported coat color: \code{'agouti'}, \code{'albino'}, \code{'black'},
#'   or \code{NA}).
#' @param probs qtl2-style allele probs. A named list with one element per
#'   chromosome. Each element is a 3D array (samples x founders x markers).
#' @param map qtl2-style marker map. A named list with one element per
#'   chromosome. Each element is a named numeric vector of marker positions
#'   in Mb.
#'
#' @return data.frame with columns: \code{id}, \code{coat} (reported),
#'   \code{albino_dt}, \code{black_dt}, \code{geno_coat} (inferred),
#'   \code{match} (logical), and per-founder allele probabilities at each
#'   locus.
#'
#' @export
compare_coat_color <- function(pheno, probs, map) {
    albino_locus <- list(chr = "7", pos = 87.108308)
    black_locus <- list(chr = "2", pos = 154.842726)

    # Verify that samples are aligned between `pheno` and `probs`.
    if (!identical(pheno$id, rownames(probs[[1]]))) {
        stop("pheno$id must match rownames(probs[[1]]) in the same order.")
    }

    # Verify that `probs` and `map` use the same chromosome order.
    if (!identical(names(map), names(probs))) {
        stop("names(map) must match names(probs) in the same order.")
    }

    # Verify that marker names are aligned between `map` and `probs` on Chr 2
    # and Chr 7.
    if (!identical(
        names(map[[black_locus$chr]]),
        dimnames(probs[[black_locus$chr]])[[3]]
    )) {
        stop("Black-locus marker names must align between map and probs.")
    }
    if (!identical(
        names(map[[albino_locus$chr]]),
        dimnames(probs[[albino_locus$chr]])[[3]]
    )) {
        stop("Albino-locus marker names must align between map and probs.")
    }

    # Get albino diplotypes.
    # Find the nearest marker to the albino locus.
    outer_diff <- abs(albino_locus$pos - map[[albino_locus$chr]])
    albino_marker <- names(outer_diff)[which.min(outer_diff)]

    # Get the allele probs at the albino marker.
    albino_probs <- probs[[albino_locus$chr]][, , albino_marker]

    # Call diplotypes.
    albino_dt <- call_diplotypes(albino_probs)
    colnames(albino_dt)[-1] <- paste0("albino_", colnames(albino_dt)[-1])

    # Get black diplotypes.
    # Find the nearest marker to the black locus.
    outer_diff <- abs(black_locus$pos - map[[black_locus$chr]])
    black_marker <- names(outer_diff)[which.min(outer_diff)]

    # Get the allele probs at the black marker.
    black_probs <- probs[[black_locus$chr]][, , black_marker]

    # Call diplotypes.
    black_dt <- call_diplotypes(black_probs)
    colnames(black_dt)[-1] <- paste0("black_", colnames(black_dt)[-1])

    # Create return value data.frame.
    result <- data.frame(
        id = pheno$id,
        coat = pheno$coat
    )
    result <- merge(result, albino_dt, by = "id")
    result <- merge(result, black_dt, by = "id")

    # Infer the coat color from the diplotype.
    # A/J (A) and NOD (D) contribute the albino allele.
    result$geno_coat <- ifelse(result$albino_dt %in% c("AA", "AD", "DD"),
                               "albino", "agouti")
    # A/J (A) and BL6 (B) contribute the black allele. Albino is epistatic
    # to black.
    result$geno_coat <- ifelse(result$black_dt %in% c("AA", "AB", "BB") &
                                   !result$albino_dt %in% c("AA", "AD", "DD"),
                               "black", result$geno_coat)

    # If the phenotype coat color and the inferred coat color match, set
    # `match` to TRUE.
    result$match <- result$coat == result$geno_coat

    # Reorder the columns before returning the final result.
    result[, c(
        "id", "coat", "albino_dt", "black_dt", "geno_coat", "match",
        paste0("albino_", LETTERS[1:8]),
        paste0("black_", LETTERS[1:8])
    )]
}
