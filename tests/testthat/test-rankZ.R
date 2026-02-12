test_that("rankZ returns same length as input", {
  x <- c(1, 2, 3, 4, 5)
  expect_length(rankZ(x), 5)
})

test_that("rankZ preserves NA", {
  x <- c(1, NA, 3, 4, 5)
  out <- rankZ(x)
  expect_true(is.na(out[2]))
  expect_true(all(!is.na(out[-2])))
})

test_that("rankZ produces roughly normal output", {
  x <- runif(100)
  out <- rankZ(x)
  expect_true(abs(mean(out, na.rm = TRUE)) < 0.5)
  expect_true(sd(out, na.rm = TRUE) > 0.5)
})
