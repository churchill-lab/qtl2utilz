#' Check sample sex from X-chromosome genoprobs
#'
#' Females are heterozygous on chrX (high heterozygosity); males are
#' hemizygous (low heterozygosity). A 2-means clustering on per-sample
#' mean X-chromosome heterozygosity separates the two groups.
#'
#' @param genoprobs A \code{calc_genoprob} object (list of 3D arrays).
#' @param x_chr Name of the X chromosome element (default \code{'X'}).
#' @param ambig_margin Samples within this distance of the cluster boundary
#'   are labeled \code{'Ambiguous'} instead of M/F. Set to 0 to disable.
#' @param seed Random seed for \code{kmeans}.
#'
#' @return data.frame with columns:
#'   \describe{
#'     \item{sample}{Sample ID from genoprobs rownames.}
#'     \item{H_X}{Mean X-chromosome heterozygosity.}
#'     \item{sex_call}{Called sex: \code{'M'}, \code{'F'}, or \code{'Ambiguous'}.}
#'     \item{cutoff}{Midpoint between the two cluster centers.}
#'     \item{margin}{Distance from the cutoff (larger = more confident).}
#'   }
#'
#' @export
genoprobs_sex_check <- function(genoprobs, x_chr = "X",
                                ambig_margin = 0.03, seed = 1) {
    if (!x_chr %in% names(genoprobs)) {
        stop(
            "Chromosome '", x_chr, "' not found in genoprobs. ",
            "Available: ", paste(names(genoprobs), collapse = ", ")
        )
    }

    hx <- .hx_from_genoprobs(genoprobs[[x_chr]])

    set.seed(seed)
    km <- kmeans(as.numeric(hx), centers = 2)

    centers <- tapply(as.numeric(hx), km$cluster, mean, na.rm = TRUE)
    male_cluster <- as.integer(names(which.min(centers)))

    sex <- ifelse(km$cluster == male_cluster, "M", "F")
    cutoff <- mean(as.numeric(centers))
    margin <- abs(as.numeric(hx) - cutoff)

    if (ambig_margin > 0) {
        sex[margin < ambig_margin] <- "Ambiguous"
    }

    data.frame(
        sample   = names(hx),
        H_X      = as.numeric(hx),
        sex_call = sex,
        cutoff   = cutoff,
        margin   = margin,
        stringsAsFactors = FALSE
    )
}


#' Cross-check sex calls against sex labels from metadata
#'
#' Compares DNA-based (or expression-based) sex calls against the sex label
#' from sample metadata (e.g. an \code{annot.samples} data.frame). Note that
#' these labels may themselves be incorrect if samples were swapped or
#' mislabeled — that is exactly what this check is designed to detect.
#'
#' @param sex_df A data.frame with a \code{sample} column and a sex-call
#'   column (e.g. from \code{\link{genoprobs_sex_check}}).
#' @param labeled_sex Sex labels for each sample. Can be one of:
#'   \itemize{
#'     \item A named character vector where names are sample IDs and values
#'       are \code{'F'} or \code{'M'} (or \code{NA}).
#'     \item A data.frame with columns \code{sample} (or matching the value
#'       of \code{sample_col}) and \code{sex} (or matching \code{sex_label_col}).
#'   }
#' @param sex_col Name of the column in \code{sex_df} containing the
#'   computed sex call (default \code{'sex_call'}).
#' @param sample_col If \code{labeled_sex} is a data.frame, the column name
#'   containing sample IDs (default \code{'sample'}).
#' @param sex_label_col If \code{labeled_sex} is a data.frame, the column name
#'   containing the sex label (default \code{'sex'}).
#'
#' @return The input data.frame with added columns:
#'   \describe{
#'     \item{sex_label}{Sex label from metadata (\code{'F'}, \code{'M'}, or NA).}
#'     \item{sex_concordant}{Logical: TRUE if \code{sex_col} matches
#'       \code{sex_label}. NA if either is missing or Ambiguous.}
#'   }
#'
#' @export
genoprobs_sex_check_labels <- function(sex_df, labeled_sex,
                                       sex_col = "sex_call",
                                       sample_col = "sample",
                                       sex_label_col = "sex") {
    if (!sex_col %in% names(sex_df)) {
        stop("Column '", sex_col, "' not found in sex_df.")
    }
    if (!"sample" %in% names(sex_df)) {
        stop("Column 'sample' not found in sex_df.")
    }

    # Build a named lookup vector from whatever form labeled_sex takes
    if (is.data.frame(labeled_sex)) {
        if (!sample_col %in% names(labeled_sex)) {
            stop("Column '", sample_col, "' not found in labeled_sex data.frame.")
        }
        if (!sex_label_col %in% names(labeled_sex)) {
            stop("Column '", sex_label_col, "' not found in labeled_sex data.frame.")
        }
        sex_vec <- setNames(
            toupper(substr(as.character(labeled_sex[[sex_label_col]]), 1, 1)),
            as.character(labeled_sex[[sample_col]])
        )
    } else if (is.character(labeled_sex) || is.factor(labeled_sex)) {
        sex_vec <- toupper(substr(as.character(labeled_sex), 1, 1))
        if (!is.null(names(labeled_sex))) {
            names(sex_vec) <- names(labeled_sex)
        }
    } else {
        stop("labeled_sex must be a named character vector or a data.frame.")
    }

    # Normalize to single-letter codes
    sex_vec[!sex_vec %in% c("F", "M")] <- NA_character_

    # Match to sex_df samples
    label <- sex_vec[sex_df$sample]

    sex_df$sex_label <- as.character(label)
    called <- sex_df[[sex_col]]
    sex_df$sex_concordant <- ifelse(
        is.na(label) | is.na(called) | called == "Ambiguous",
        NA,
        label == called
    )

    sex_df
}


#' Compute per-sample heterozygosity on a single chromosome
#'
#' For an n x f x m genoprobs array (samples x founders x markers),
#' heterozygosity at each position is \code{1 - sum(p^2)} across founders.
#' Returns the per-sample mean across all positions.
#'
#' @param gp_chr A 3D array (samples x founders x markers).
#' @return Named numeric vector of mean heterozygosity per sample.
#' @keywords internal
.hx_from_genoprobs <- function(gp_chr) {
    if (length(dim(gp_chr)) != 3) {
        stop("gp_chr must be a 3D array: samples x founders x markers.")
    }
    # H[sample, pos] = 1 - sum_f p^2
    H_sp <- 1 - apply(gp_chr^2, c(1, 3), sum, na.rm = TRUE)
    rowMeans(H_sp, na.rm = TRUE)
}
