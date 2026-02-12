# qtl2utilz

Utilities for working with [R/qtl2](https://kbroman.org/qtl2/) and [GBRS](https://github.com/churchill-lab/GBRS) (Genotype By RNA-Seq) data. Enables importing GBRS genoprobs and read counts into R/qtl2 format, synchronizing and interpolating genoprobs across marker sets, and comparing GBRS vs array-based genotypes (e.g. MUGA).

Developed in the [Churchill Lab](https://churchill-lab.jax.org/) at The Jackson Laboratory.

## Installation

```r
# Install from GitHub (requires remotes or devtools)
remotes::install_github("churchill-lab/qtl2utilz")
```

## Features

- **GBRS workflow**: Find GBRS output files, build genotype probabilities, and build expression matrices from read counts
- **Genoprobs utilities**: Sync markers/samples, interpolate to new marker grids, combine samples across batches
- **Comparison**: Correlate GBRS vs array-based genoprobs (e.g. MUGA) for QC
- **Flexible input**: Column aliases for `sample_id`, `marker_id`, `chr`, `pos` (e.g. mouse_id, marker, chrom, position)

## Quick start

### Import GBRS data into R/qtl2

```r
library(qtl2utilz)
library(qtl2)

# Find GBRS output files in a directory
files_tbl <- gbrs_find_files("/path/to/gbrs/output")

# Build genotype probabilities (requires a marker map)
markers <- read.csv("marker_map.csv")  # columns: marker_id, chr, pos (or aliases)
genoprobs <- gbrs_build_genoprobs(files_tbl, markers)

# Build expression matrix from read counts
expr <- gbrs_build_counts(files_tbl)
```

### Compare GBRS vs MUGA genoprobs

```r
# 1) Sync each genoprobs to its marker set
gbrs <- genoprobs_sync_markers(genoprobs_GBRS, markers_GBRS)
muga <- genoprobs_sync_markers(genoprobs_MUGA, markers_MUGA)

# 2) Interpolate GBRS onto MUGA marker grid
gbrs_on_muga <- genoprobs_interpolate(gbrs$genoprobs, gbrs$map, muga$map)

# 3) Restrict to common samples
aligned <- genoprobs_sync_samples(gbrs_on_muga, muga$genoprobs)

# 4) Per-sample correlation (QC)
cor_tbl <- genoprobs_correlate(aligned$genoprobs_1, aligned$genoprobs_2)
summary(cor_tbl$correlation)
```

## Documentation

- [Vignette: GBRS to R/qtl2](vignettes/gbrs-workflow.html) — full workflow with examples
- Function help: `?gbrs_build_genoprobs`, `?genoprobs_interpolate`, etc.

## Dependencies

- [R/qtl2](https://kbroman.org/qtl2/)
- [qtl2convert](https://github.com/rqtl/qtl2convert)
- dplyr, tidyr, tibble, stringr

## License

MIT
