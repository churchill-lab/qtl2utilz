test_that("markers_sort orders by chr then pos", {
  df <- data.frame(
    marker_id = c("m3", "m1", "m2"),
    chr = c("2", "1", "1"),
    pos = c(1, 3, 2)
  )
  out <- markers_sort(df)
  # chr 1 before chr 2; within chr 1, pos ascending (2 before 3)
  expect_identical(out$marker_id, c("m2", "m1", "m3"))
  expect_identical(out$chr, c("1", "1", "2"))
  expect_identical(out$pos, c(2, 3, 1))
})

test_that("markers_sort accepts chr alias", {
  df <- data.frame(marker = c("m1", "m2"), chrom = c("1", "1"), pos = c(1, 2))
  out <- markers_sort(df)
  expect_named(out, c("marker_id", "chr", "pos"))
  expect_identical(out$marker_id, c("m1", "m2"))
})

test_that("markers_sort puts X after autosomes", {
  df <- data.frame(marker_id = c("mx", "m1"), chr = c("X", "1"), pos = c(1, 1))
  out <- markers_sort(df)
  expect_identical(out$chr, c("1", "X"))
})
