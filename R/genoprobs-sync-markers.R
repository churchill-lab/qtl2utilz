#' Synchronize genoprobs to a marker data frame.
#'
#' Subset and reorder genoprobs so that each chromosome contains only markers
#' that appear in \code{markers_df}, in the order defined by the data frame
#' (after sorting by chromosome and position). Useful to align genoprobs with
#' a standard marker set (e.g. for correlation with another platform) or to
#' drop markers not in a given map. Positions are converted to bp if needed.
#'
#' @param genoprobs qtl2-style genoprobs (named list of 3D arrays by chromosome).
#' @param markers_df Data frame with at least marker ID, chromosome, and
#'   position columns. Aliases are accepted and normalized via
#'   \code{resolve_col_markers()}.
#' @param unit Position unit: \code{'auto'} (detect from values), \code{'Mb'}, or \code{'bp'}.
#'
#' @return A list with components:
#'   \item{genoprobs}{Genoprobs restricted to common chromosomes and
#'     markers present in both genoprobs and \code{markers_df}, in map order.}
#'   \item{markers}{Marker data frame with positions in bp, sorted by chr/pos;
#'     only markers that appear in \code{genoprobs} and \code{map} are included.}
#'   \item{map}{qtl2-style marker map (named list of named numeric position vectors).}
#'   \item{dropped_markers_not_in_gp}{Per-chromosome list of marker IDs in
#'     \code{markers_df} but not in genoprobs.}
#'   \item{dropped_markers_not_in_map}{Per-chromosome list of marker IDs in
#'     genoprobs but not in \code{markers_df}.}
#'
#' @export
genoprobs_sync_markers <- function(genoprobs,
                                   markers_df,
                                   unit       = c('auto','Mb','bp')) {
    unit <- match.arg(unit)
    markers_df <- resolve_col_markers(markers_df)

    if(!all(c('marker_id', 'chr', 'pos') %in% names(markers_df))) {
        stop('markers_df must contain "marker_id", "chr", and "pos"')
    }

    # extract a plain list of chr arrays so subsetting [common_chr] is safe
    genoprobs_list <- extract_chr_list(genoprobs)

    # normalize marker positions to bp and sort by chr, pos for consistent order
    markers_bp <- positions_to_bp(markers_df, unit = unit)
    markers_bp <- markers_sort(markers_bp)

    # build qtl2-style marker map (named list of named position vectors) for downstream use
    map <- qtl2convert::map_df_to_list(
        markers_bp,
        marker = 'marker_id',
        chr    = 'chr',
        pos    = 'pos'
    )

    # restrict to chromosomes present in both; plain list subsetting is safe here
    common_chr <- intersect(names(genoprobs_list), names(map))
    if(!length(common_chr)) {
        stop('No overlapping chromosome names between genoprobs and markers/map.')
    }

    out <- genoprobs_list[common_chr]

    dropped_markers_not_in_gp  <- list()
    dropped_markers_not_in_map <- list()

    # per chromosome: keep only markers that appear in both, in map order
    for(chr in common_chr) {

        pr <- out[[chr]]
        if(length(dim(pr)) != 3) {
            stop(sprintf('Chr %s genoprobs is not a 3D array.', chr))
        }

        gp_markers <- dimnames(pr)[[3]]
        if(is.null(gp_markers)) {
            stop(sprintf('Chr %s genoprobs has NULL marker dimnames.', chr))
        }

        map_markers <- names(map[[chr]])
        if(is.null(map_markers)) {
            stop(sprintf('Chr %s map has NULL names().', chr))
        }

        # keep markers in map order (sorted by chr/pos already); only those in genoprobs too
        keep <- map_markers[map_markers %in% gp_markers]

        dropped_markers_not_in_gp[[chr]]  <- setdiff(map_markers, gp_markers)
        dropped_markers_not_in_map[[chr]] <- setdiff(gp_markers, map_markers)

        if(!length(keep)) {
            stop(sprintf('Chr %s has zero overlapping markers between genoprobs and markers_df.', chr))
        }

        idx <- match(keep, gp_markers)
        out[[chr]] <- pr[, , idx, drop = FALSE]
        map[[chr]] <- map[[chr]][keep]

        stopifnot(identical(dimnames(out[[chr]])[[3]], names(map[[chr]])))
    }

    # restore original attributes/class so result behaves like input genoprobs
    attributes(out) <- attributes(genoprobs)
    class(out) <- class(genoprobs)

    # subset markers to only those present in genoprobs/map (same set per chr)
    kept_rows <- logical(nrow(markers_bp))
    for(chr in common_chr) {
        kept_rows <- kept_rows | (
            markers_bp$chr == chr &
            markers_bp$marker_id %in% names(map[[chr]])
        )
    }
    markers_out <- markers_bp[kept_rows, , drop = FALSE]

    list(
        genoprobs                  = out,
        markers                    = markers_out,
        map                        = map[common_chr],
        dropped_markers_not_in_gp  = dropped_markers_not_in_gp,
        dropped_markers_not_in_map = dropped_markers_not_in_map
    )
}

