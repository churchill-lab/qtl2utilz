#' Check sample sex from RNA-seq expression of sex-linked genes
#'
#' Uses Xist expression (high in females, near-zero in males) as the primary
#' signal. Optionally cross-checks against Y-linked genes (Kdm5d, Eif2s3y,
#' Uty, Ddx3y) if present in the counts matrix.
#'
#' @param counts A samples-by-genes matrix of expression counts (raw or
#'   normalized). Column names must be Ensembl gene IDs.
#' @param xist_id Ensembl ID for Xist (default \code{'ENSMUSG00000086503'}).
#' @param y_gene_ids Character vector of Y-linked gene Ensembl IDs to use
#'   for male-specific expression. Only genes present in
#'   \code{colnames(counts)} are used. Set to \code{NULL} to skip Y-gene
#'   check entirely.
#' @param ambig_margin Samples within this log2 distance of the cluster
#'   boundary are labeled \code{'Ambiguous'}.
#' @param seed Random seed for \code{kmeans}.
#'
#' @return data.frame with columns:
#'   \describe{
#'     \item{sample}{Sample ID from counts rownames.}
#'     \item{xist_log2}{log2(Xist count + 1).}
#'     \item{sex_xist}{Sex call from Xist alone.}
#'     \item{n_y_genes}{Number of Y-linked genes found in the counts matrix.}
#'     \item{y_mean_log2}{Mean log2(count + 1) across available Y genes.}
#'     \item{sex_y}{Sex call from Y genes alone (NA if none found).}
#'     \item{sex_call}{Consensus call combining Xist and Y-gene evidence.}
#'   }
#'
#' @export
expression_sex_check <- function(
        counts,
        xist_id = "ENSMUSG00000086503",
        y_gene_ids = c(
            "ENSMUSG00000056673",  # Kdm5d
            "ENSMUSG00000069049",  # Eif2s3y
            "ENSMUSG00000068457",  # Uty
            "ENSMUSG00000069045"   # Ddx3y
        ),
        ambig_margin = 0.5,
        seed = 1
) {
    counts <- as.matrix(counts)
    if (is.null(rownames(counts))) stop("counts must have rownames (sample IDs).")
    if (is.null(colnames(counts))) stop("counts must have colnames (gene IDs).")

    # --- Xist-based call ---
    if (!xist_id %in% colnames(counts)) {
        stop("Xist gene (", xist_id, ") not found in colnames(counts). ",
             "Column names must be Ensembl gene IDs.")
    }

    x_vals <- log2(counts[, xist_id] + 1)

    if (length(unique(x_vals[is.finite(x_vals)])) < 2) {
        stop("Xist expression has fewer than 2 distinct values; cannot separate sexes by kmeans.")
    }

    set.seed(seed)
    km <- kmeans(x_vals, centers = 2)

    centers <- tapply(x_vals, km$cluster, mean)
    male_cluster <- as.integer(names(which.min(centers)))

    sex_xist <- ifelse(km$cluster == male_cluster, "M", "F")
    xist_cutoff <- mean(centers)
    xist_margin <- abs(x_vals - xist_cutoff)
    sex_xist[xist_margin < ambig_margin] <- "Ambiguous"

    # --- Y-gene-based call ---
    y_present <- character(0)
    if (!is.null(y_gene_ids)) {
        y_present <- y_gene_ids[y_gene_ids %in% colnames(counts)]
    }
    n_y <- length(y_present)

    sex_y <- rep(NA_character_, nrow(counts))
    y_mean_log2 <- rep(NA_real_, nrow(counts))

    if (n_y > 0) {
        y_mat <- log2(counts[, y_present, drop = FALSE] + 1)
        y_mean_log2 <- rowMeans(y_mat, na.rm = TRUE)

        if (length(unique(y_mean_log2[is.finite(y_mean_log2)])) >= 2) {
            set.seed(seed)
            km_y <- kmeans(y_mean_log2, centers = 2)
            centers_y <- tapply(y_mean_log2, km_y$cluster, mean)
            male_cluster_y <- as.integer(names(which.max(centers_y)))

            sex_y <- ifelse(km_y$cluster == male_cluster_y, "M", "F")
            y_cutoff <- mean(centers_y)
            y_margin <- abs(y_mean_log2 - y_cutoff)
            sex_y[y_margin < ambig_margin] <- "Ambiguous"
        } else {
            sex_y[] <- NA_character_
        }
    }

    # --- Consensus ---
    sex_call <- sex_xist
    if (n_y > 0 && any(!is.na(sex_y))) {
        both_defined <- sex_xist %in% c("M", "F") & sex_y %in% c("M", "F")
        disagree <- both_defined & (sex_xist != sex_y)
        sex_call[disagree] <- "Ambiguous"
    }

    data.frame(
        sample      = rownames(counts),
        xist_log2   = x_vals,
        sex_xist    = sex_xist,
        n_y_genes   = n_y,
        y_mean_log2 = y_mean_log2,
        sex_y       = sex_y,
        sex_call    = sex_call,
        stringsAsFactors = FALSE
    )
}
