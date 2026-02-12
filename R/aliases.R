.qtl2utilz_aliases <- list(
    sample_id = c('sample_id', 'sample.id', 'sampleID', 'sample', 'mouse_id', 'mouse.id', 'mouseID', 'mouse'),
    marker_id = c('marker_id', 'marker.id', 'markerID', 'markerid', 'marker'),
    chr       = c('chr', 'chrom', 'chromosome', 'chrom_id', 'chromosome_id'),
    pos       = c('pos', 'position')
)

# Normalize column tokens for flexible matching.
# Lowercase + collapse punctuation/whitespace to underscores.
.normalize_col_token <- function(x) {
    x <- tolower(as.character(x))
    x <- gsub('[^a-z0-9]+', '_', x)
    x <- gsub('^_+|_+$', '', x)
    x
}

# Resolve one canonical key to an existing column name.
# Internal utility used by the wrapper helpers below.
.resolve_col <- function(df, key, aliases = .qtl2utilz_aliases) {
    if(!is.data.frame(df)) {
        stop('Expected a data.frame for column resolution.')
    }
    if(!key %in% names(aliases)) {
        stop('Unknown alias key: ', key)
    }

    candidates <- unique(aliases[[key]])
    norm_candidates <- unique(.normalize_col_token(candidates))
    norm_names <- .normalize_col_token(names(df))
    hit_idx <- which(norm_names %in% norm_candidates)
    hits <- names(df)[hit_idx]

    if(length(hits) == 1) {
        return(hits)
    }
    if(length(hits) == 0) {
        stop(
            'Missing column for "', key, '". Accepted aliases: ',
            paste(candidates, collapse = ', '),
            ' (case-insensitive; punctuation-insensitive).'
        )
    }

    stop(
        'Ambiguous columns for "', key, '": ',
        paste(hits, collapse = ', '),
        '. Keep only one of these aliases.'
    )
}

# Resolve and normalize sample key column to canonical sample_id.
# Returns a data.frame with the renamed column.
# @keywords internal
resolve_col_samples <- function(df) {
    sample_col <- .resolve_col(df, 'sample_id')
    names(df)[names(df) == sample_col] <- 'sample_id'
    df
}

# Resolve and normalize marker/map key columns to canonical names.
# Returns a data.frame with marker_id, chr, and pos column names.
# @keywords internal
resolve_col_markers <- function(df) {
    marker_col <- .resolve_col(df, 'marker_id')
    chr_col <- .resolve_col(df, 'chr')
    pos_col <- .resolve_col(df, 'pos')

    names(df)[names(df) == marker_col] <- 'marker_id'
    names(df)[names(df) == chr_col] <- 'chr'
    names(df)[names(df) == pos_col] <- 'pos'
    df
}