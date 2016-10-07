context("util")

test_that("assertions work", {
  expect_error(assert_scalar(NULL), "must be a scalar")
  expect_error(assert_scalar(numeric(0)), "must be a scalar")
  expect_error(assert_scalar(1:2), "must be a scalar")

  expect_error(assert_scalar_logical(1), "must be logical")
  expect_error(assert_scalar_logical(NA), "must not be NA")
  expect_error(assert_scalar_logical(c(TRUE, TRUE)), "must be a scalar")

  expect_error(assert_size(1.5), "must be integer")
  expect_error(assert_size(-2L), "must be nonnegative")
  expect_error(assert_size(NA_integer_), "must not be NA")
  expect_error(assert_size(c(1, 2)), "must be a scalar")
})
