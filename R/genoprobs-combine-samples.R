#' Union genoprobs by sample
#'
#' Combine two qtl2-style genoprobs objects by taking the union of samples.
#' Chromosomes, founder/strain names, and marker names must match exactly
#' between the two objects. Useful for merging genoprobs from the same
#' cross computed in separate batches or from different marker sets that
#' have been synchronized first.
#'
#' @param genoprobs_1 First genoprobs object (list of 3D arrays by chromosome).
#' @param genoprobs_2 Second genoprobs object (same structure as \code{genoprobs_1}).
#' @param samples_1 Character vector of sample IDs to keep from \code{genoprobs_1};
#'   \code{NULL} means use all samples.
#' @param samples_2 Character vector of sample IDs to add from \code{genoprobs_2};
#'   \code{NULL} means use all samples.
#' @param sort If \code{TRUE}, order the union of samples alphabetically.
#' @param duplicates If \code{"first"}, samples present in both are taken from
#'   \code{genoprobs_1} and not duplicated from \code{genoprobs_2}. If
#'   \code{"error"}, an error is raised when any sample appears in both.
#'
#' @return A list with components:
#'   \item{genoprobs}{Combined genoprobs object (union of samples, same chromosomes and markers).}
#'   \item{union_samples}{Character vector of sample IDs in the result.}
#'   \item{duplicates}{Character vector of sample IDs that appeared in both inputs.}
#'   \item{dropped_from_2_due_to_duplicates}{When \code{duplicates = "first"}, samples from \code{genoprobs_2} not added because they were already in \code{genoprobs_1}.}
#'   \item{added_from_2}{Sample IDs that came from \code{genoprobs_2} only.}
#'
#' @details
#' Chromosome sets must be identical (same names). For each chromosome, founder
#' names (columns) and marker names (third dimension) must match. Attributes and
#' class of the result are taken from \code{genoprobs_1}.
#'
#' @export
genoprobs_combine_samples <- function(genoprobs_1,
                                    genoprobs_2,
                                    samples_1 = NULL,   # NULL = all
                                    samples_2 = NULL,   # NULL = all
                                    sort = TRUE,
                                    duplicates = c("first","error")) {
    duplicates <- match.arg(duplicates)

    # Work with plain lists to avoid subclass subsetting issues when we
    # reorder or subset chromosomes later.
    gp1 <- extract_chr_list(genoprobs_1)
    gp2 <- extract_chr_list(genoprobs_2)

    # Chromosomes must match exactly (same set and names) so we can
    # merge sample rows without ambiguity.
    chr1 <- names(gp1)
    chr2 <- names(gp2)
    if(!setequal(chr1, chr2)) {
        stop(sprintf("Chromosome sets differ.\nOnly in genoprobs_1: %s\nOnly in genoprobs_2: %s",
                     paste(setdiff(chr1, chr2), collapse = ", "),
                     paste(setdiff(chr2, chr1), collapse = ", ")))
    }
    # Use a consistent chromosome order (from genoprobs_1); align gp2 to it.
    chrs <- chr1
    gp2 <- gp2[chrs]

    # Markers and founders must match exactly for each chromosome so that
    # combined arrays have identical structure.
    for(chr in chrs) {
        if(!identical(colnames(gp1[[chr]]), colnames(gp2[[chr]]))) {
            stop(sprintf("Founder/strain names differ on chr %s.", chr))
        }
        if(!identical(dimnames(gp1[[chr]])[[3]], dimnames(gp2[[chr]])[[3]])) {
            stop(sprintf("Marker names/order differ on chr %s.", chr))
        }
    }

    # Select which samples to use from each object (NULL = all).
    chr0 <- chrs[1]
    s1_all <- rownames(gp1[[chr0]])
    s2_all <- rownames(gp2[[chr0]])
    if(is.null(s1_all) || is.null(s2_all)) stop("One or both genoprobs objects have NULL sample names.")

    s1 <- if(is.null(samples_1)) s1_all else intersect(s1_all, samples_1)
    s2 <- if(is.null(samples_2)) s2_all else intersect(s2_all, samples_2)
    if(!length(s1)) stop("No samples selected from genoprobs_1.")
    if(!length(s2)) stop("No samples selected from genoprobs_2.")

    # Detect duplicates; error if user requested strict no-overlap.
    dup <- intersect(s1, s2)
    if(length(dup) && duplicates == "error") {
        stop(sprintf("Duplicate samples found (%d). Example: %s",
                     length(dup), paste(head(dup, 5), collapse = ", ")))
    }

    # With duplicates = "first": keep all from 1, add from 2 only samples not in 1.
    keep1 <- s1
    add2  <- setdiff(s2, s1)

    union_samples <- c(keep1, add2)
    if(sort) union_samples <- sort(union_samples)

    # Build output list: one array per chromosome, filled from gp1 then gp2.
    out <- setNames(vector("list", length(chrs)), chrs)

    for(chr in chrs) {
        pr1 <- gp1[[chr]]
        pr2 <- gp2[[chr]]

        founders <- colnames(pr1)
        markers  <- dimnames(pr1)[[3]]

        # Pre-allocate result array for union of samples x founders x markers.
        res <- array(NA_real_,
                     dim = c(length(union_samples), length(founders), length(markers)),
                     dimnames = list(union_samples, founders, markers))

        # Fill from genoprobs_1 for all keep1 samples (rows in result).
        i_out_1 <- match(keep1, union_samples)
        i_in_1  <- match(keep1, rownames(pr1))
        res[i_out_1, , ] <- pr1[i_in_1, , , drop = FALSE]

        # Fill from genoprobs_2 only for samples that are not in genoprobs_1.
        if(length(add2)) {
            i_out_2 <- match(add2, union_samples)
            i_in_2  <- match(add2, rownames(pr2))
            res[i_out_2, , ] <- pr2[i_in_2, , , drop = FALSE]
        }

        out[[chr]] <- res
    }

    # Restore qtl2 attributes and class from first input so downstream
    # code (e.g. qtl2) still sees a valid genoprobs object.
    attributes(out) <- attributes(genoprobs_1)
    class(out) <- class(genoprobs_1)

    list(
        genoprobs = out,
        union_samples = union_samples,
        duplicates = dup,
        dropped_from_2_due_to_duplicates = if(duplicates == "first") dup else character(),
        added_from_2 = add2
    )
}
