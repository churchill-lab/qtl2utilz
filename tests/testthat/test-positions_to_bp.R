test_that("positions_to_bp.numeric converts Mb to bp", {
  expect_equal(positions_to_bp(c(1, 2, 3), unit = "Mb"), c(1e6, 2e6, 3e6))
})

test_that("positions_to_bp.numeric leaves bp unchanged", {
  x <- c(1e6, 2e6, 3e6)
  expect_equal(positions_to_bp(x, unit = "bp"), x)
})

test_that("positions_to_bp.numeric auto detects Mb when max < 2000", {
  expect_equal(positions_to_bp(c(1, 2), unit = "auto"), c(1e6, 2e6))
})

test_that("positions_to_bp.numeric auto detects bp when max >= 2000", {
  x <- c(1000, 2000, 3000)
  expect_equal(positions_to_bp(x, unit = "auto"), x)
})

test_that("positions_to_bp.data.frame converts position column", {
  df <- data.frame(chr = "1", pos = c(1, 2, 3), id = 1:3)
  out <- positions_to_bp(df, unit = "Mb")
  expect_equal(out$pos, c(1e6, 2e6, 3e6))
  expect_equal(out$chr, df$chr)
  expect_equal(out$id, df$id)
})

test_that("positions_to_bp errors on unsupported type", {
  expect_error(positions_to_bp("a"), "Unsupported type")
})
