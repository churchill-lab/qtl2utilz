test_that("resolve_col_samples normalizes sample column aliases", {
  df <- data.frame(mouse_id = 1:3, x = letters[1:3])
  out <- qtl2utilz:::resolve_col_samples(df)
  expect_named(out, c("sample_id", "x"))
  expect_identical(out$sample_id, 1:3)
})

test_that("resolve_col_markers normalizes marker columns", {
  df <- data.frame(marker = c("m1", "m2"), chrom = c("1", "1"), position = c(10, 20))
  out <- qtl2utilz:::resolve_col_markers(df)
  expect_named(out, c("marker_id", "chr", "pos"))
  expect_identical(out$marker_id, c("m1", "m2"))
  expect_identical(out$chr, c("1", "1"))
  expect_identical(out$pos, c(10, 20))
})

test_that("resolve_col_samples errors on missing column", {
  df <- data.frame(x = 1:3)
  expect_error(qtl2utilz:::resolve_col_samples(df), "Missing column")
})

test_that("resolve_col_markers errors on missing required columns", {
  df <- data.frame(marker = c("m1"), chr = c("1"))
  expect_error(qtl2utilz:::resolve_col_markers(df), "Missing column")
})
