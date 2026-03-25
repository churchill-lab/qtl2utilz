test_that("rank_z returns same length as input", {
  x <- c(1, 2, 3, 4, 5)
  expect_length(rank_z(x), 5)
})

test_that("rank_z preserves NA", {
  x <- c(1, NA, 3, 4, 5)
  out <- rank_z(x)
  expect_true(is.na(out[2]))
  expect_true(all(!is.na(out[-2])))
})

test_that("rank_z produces roughly normal output", {
  x <- runif(100)
  out <- rank_z(x)
  expect_true(abs(mean(out, na.rm = TRUE)) < 0.5)
  expect_true(sd(out, na.rm = TRUE) > 0.5)
})
