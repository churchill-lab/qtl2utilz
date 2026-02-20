#' Correlate genoprobs between two objects (per sample)
#'
#' Compute per-sample similarity between two R/qtl2-style genotype probability
#' objects (a \code{calc_genoprob} list of 3D arrays, one per chromosome).
#'
#' @section Big picture:
#' For each chromosome, genoprobs are stored as a 3D array with dimensions:
#' \itemize{
#'   \item \code{dim 1}: samples (individuals)
#'   \item \code{dim 2}: states/haplotypes (e.g., founders A--H)
#'   \item \code{dim 3}: markers (positions along the chromosome)
#' }
#' For a given sample, the \code{haplotype x marker} probabilities can be
#' "flattened" into one long numeric vector. This function compares the two
#' objects sample-by-sample by correlating those flattened vectors across all
#' overlapping chromosomes.
#'
#' @section Why use accumulation instead of concatenation:
#' A straightforward implementation would concatenate probability vectors across
#' chromosomes and call \code{cor(v1, v2)}. That can be slow and memory-heavy
#' because repeated \code{c(v, new_values)} reallocates vectors many times.
#' For Pearson correlation, we can avoid building the full concatenated vectors:
#' Pearson correlation can be computed from running totals
#' \code{sum(x)}, \code{sum(y)}, \code{sum(x^2)}, \code{sum(y^2)}, \code{sum(xy)}
#' and the number of features. This version uses that approach for speed and
#' lower memory use. (For Spearman/Kendall, rank correlation requires the full
#' vectors, so this function falls back to concatenation for those methods.)
#'
#' @section Alignment assumptions (important):
#' This function assumes the two objects are comparable:
#' \itemize{
#'   \item same state/haplotype set (dim 2) per chromosome
#'   \item same marker grid (dim 3) per chromosome, in the same order
#' }
#' You can enforce this upstream with helpers such as
#' \code{genoprobs_sync_markers()} and \code{genoprobs_sync_samples()}.
#'
#' By default, this function also performs lightweight safety checks using
#' dimnames (when present) to ensure state and marker ordering matches before
#' comparing. This prevents silent "wrong-but-plausible" correlations when, for
#' example, marker order differs between objects.
#'
#' @param genoprobs_1 First R/qtl2 \code{calc_genoprob} object.
#' @param genoprobs_2 Second R/qtl2 \code{calc_genoprob} object.
#' @param threshold Numeric; correlations below this value are flagged in
#'   \code{flag_low}. Default is \code{0.80}.
#' @param method Correlation method. One of \code{"pearson"}, \code{"spearman"},
#'   or \code{"kendall"} (passed to \code{\link[stats]{cor}}). Default \code{"pearson"}.
#' @param check_dimnames Logical; if \code{TRUE} (default), verify that haplotype
#'   and marker dimnames match exactly (when available) for each chromosome. Set
#'   to \code{FALSE} only if you have already enforced alignment and want to skip
#'   these checks.
#'
#' @return A data frame with one row per overlapping sample:
#' \describe{
#'   \item{sample_id}{Sample ID (intersection of samples across objects).}
#'   \item{correlation}{Correlation between the two genoprobs vectors for that sample.}
#'   \item{flag_low}{\code{TRUE} when correlation is \code{NA} or < \code{threshold}.}
#' }
#'
#' @examples
#' \dontrun{
#' # Suppose gp_a and gp_b were produced from different sources (e.g. GBRS vs MUGA)
#' # and aligned to the same marker grid and sample set:
#' gp_b2 <- genoprobs_sync_markers(gp_b, gp_a)    # align marker grids
#' gp_b2 <- genoprobs_sync_samples(gp_b2, gp_a)   # align sample IDs/order (if desired)
#'
#' out <- genoprobs_correlate(gp_a, gp_b2, threshold = 0.9)
#' head(out)
#' subset(out, flag_low)
#' }
#'
#' @export
genoprobs_correlate <- function(genoprobs_1,
                                genoprobs_2,
                                threshold = 0.80,
                                method = c("pearson", "spearman", "kendall"),
                                check_dimnames = TRUE) {
  method <- match.arg(method)

  # ---------------------------------------------------------------------------
  # STEP 1: Determine overlapping chromosomes and samples
  # ---------------------------------------------------------------------------
  # We only compare chromosomes present in both objects.
  chrs <- intersect(names(genoprobs_1), names(genoprobs_2))
  if (!length(chrs)) stop("No overlapping chromosomes.")

  # We only compare samples present in both objects. Use the first shared
  # chromosome to obtain sample names (qtl2 objects should be consistent across
  # chromosomes, but we still validate per-chromosome sample consistency below).
  s1 <- rownames(genoprobs_1[[chrs[1]]])
  s2 <- rownames(genoprobs_2[[chrs[1]]])
  common_samples <- intersect(s1, s2)
  if (!length(common_samples)) stop("No overlapping samples.")

  # Output container
  res <- data.frame(
    sample_id = common_samples,
    correlation = NA_real_,
    flag_low = FALSE,
    stringsAsFactors = FALSE
  )

  # Helper: validate per-chromosome alignment in shape and (optionally) dimnames
  validate_chr <- function(A, B, chr) {
    # Expect [samples x states/haplotypes x markers]
    if (length(dim(A)) != 3L || length(dim(B)) != 3L) {
      stop("Expected 3D arrays for chr ", chr, ".")
    }
    # Must match states and markers (sample count may differ).
    if (dim(A)[2] != dim(B)[2] || dim(A)[3] != dim(B)[3]) {
      stop("Mismatch in states (dim 2) or markers (dim 3) for chr ", chr, ".")
    }

    if (check_dimnames) {
      dna <- dimnames(A)
      dnb <- dimnames(B)

      # If dimnames exist, enforce exact match of state and marker ordering.
      if (!is.null(dna) && !is.null(dnb)) {
        # states/haplotypes
        if (!is.null(dna[[2]]) && !is.null(dnb[[2]]) && !identical(dna[[2]], dnb[[2]])) {
          stop("State/haplotype dimnames differ on chr ", chr, ".")
        }
        # markers
        if (!is.null(dna[[3]]) && !is.null(dnb[[3]]) && !identical(dna[[3]], dnb[[3]])) {
          stop("Marker dimnames/order differ on chr ", chr, ".")
        }
      }
    }
    invisible(TRUE)
  }

  # ---------------------------------------------------------------------------
  # STEP 2A: Rank correlations (Spearman/Kendall) require full vectors
  # ---------------------------------------------------------------------------
  if (method != "pearson") {
    for (i in seq_along(common_samples)) {
      s <- common_samples[i]

      # Collect per-chromosome chunks then concatenate once
      chunks1 <- vector("list", length(chrs))
      chunks2 <- vector("list", length(chrs))

      for (k in seq_along(chrs)) {
        chr <- chrs[k]
        A <- genoprobs_1[[chr]]
        B <- genoprobs_2[[chr]]

        validate_chr(A, B, chr)

        # Flatten [states x markers] for this sample into a vector
        chunks1[[k]] <- as.vector(A[s, , ])
        chunks2[[k]] <- as.vector(B[s, , ])
      }

      v1 <- unlist(chunks1, use.names = FALSE)
      v2 <- unlist(chunks2, use.names = FALSE)

      res$correlation[i] <- stats::cor(v1, v2,
                                       use = "pairwise.complete.obs",
                                       method = method)
    }

  } else {
    # -------------------------------------------------------------------------
    # STEP 2B: Fast Pearson correlation using running totals (no concatenation)
    # -------------------------------------------------------------------------
    # Pearson correlation over many features can be computed from totals:
    #   r = (sumxy - sumx*sumy/n) / sqrt((sumx2 - sumx^2/n)*(sumy2 - sumy^2/n))
    # where n is the number of paired features included in the sums.
    for (i in seq_along(common_samples)) {
      s <- common_samples[i]

      sumx <- 0
      sumy <- 0
      sumx2 <- 0
      sumy2 <- 0
      sumxy <- 0
      n_tot <- 0L

      for (chr in chrs) {
        A <- genoprobs_1[[chr]]
        B <- genoprobs_2[[chr]]

        validate_chr(A, B, chr)

        x <- as.vector(A[s, , ])
        y <- as.vector(B[s, , ])

        # Use only paired finite values (similar in spirit to pairwise.complete.obs).
        ok <- is.finite(x) & is.finite(y)
        if (!any(ok)) next

        x <- x[ok]
        y <- y[ok]
        n <- length(x)
        n_tot <- n_tot + n

        sumx  <- sumx  + sum(x)
        sumy  <- sumy  + sum(y)
        sumx2 <- sumx2 + sum(x * x)
        sumy2 <- sumy2 + sum(y * y)
        sumxy <- sumxy + sum(x * y)
      }

      if (n_tot < 2L) {
        res$correlation[i] <- NA_real_
      } else {
        denom <- sqrt((sumx2 - (sumx * sumx) / n_tot) *
                      (sumy2 - (sumy * sumy) / n_tot))
        if (!is.finite(denom) || denom == 0) {
          res$correlation[i] <- NA_real_
        } else {
          res$correlation[i] <- (sumxy - (sumx * sumy) / n_tot) / denom
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # STEP 3: Flag low/failed correlations for QC
  # ---------------------------------------------------------------------------
  # Treat NA as flagged because it indicates insufficient variance or data.
  res$flag_low <- is.na(res$correlation) | (res$correlation < threshold)
  res
}
