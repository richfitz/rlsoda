context("interface")

## These all come directly from dde:tests/testthat/test-interface.R
test_that("invalid args", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 1, length.out = 11)
  expect_error(rlsoda(y, times, "lorenz", p, dllname = "lorenz",
                      unknown = 1),
               "Invalid dot arguments")
})

test_that("R function with dllname", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 1, length.out = 11)
  growth <- function(t, y, p) {
    y * p
  }
  total <- function(t, y, p) {
    sum(y)
  }

  expect_error(rlsoda(y, times, growth, p, dllname = "lorenz",
                      n_out = 1, output = total),
               "dllname must not be given")
})

test_that("R function with dllname", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 1, length.out = 11)
  growth <- function(t, y, p) {
    y * p
  }
  total <- function(t, y, p) {
    sum(y)
  }

  expect_error(rlsoda(y, times, "lorenz", p, dllname = "lorenz",
                      n_out = 1, output = total),
               "output must be a compiled function")
  expect_error(rlsoda(y, times, growth, p,
                      n_out = 2L, output = "lorenz_output"),
               "output must be an R function")
})

test_that("output with no n_out", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 1, length.out = 11)
  growth <- function(t, y, p) {
    y * p
  }
  total <- function(t, y, p) {
    sum(y)
  }

  expect_error(rlsoda(y, times, "lorenz", p, dllname = "lorenz",
                      output = "lorenz_output"),
               "n_out must be specified")
  expect_error(rlsoda(y, times, growth, p,
                      output = total),
               "n_out must be specified")
})

test_that("invalid function input", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 1, length.out = 11)
  expect_error(rlsoda(y, times, 1, p),
               "Invalid input for 'func'")
})

test_that("Missing output function", {
  growth <- function(t, y, p) y
  output <- function(t, y, p) sum(y)
  expect_error(rlsoda(1, 0:10, growth, NULL, n_out = 1),
               "Invalid input for 'output'")
})

test_that("Native Symbol interface", {
  p <- c(10, 28, 8 / 3)
  y0 <- c(10, 1, 1)
  times <- seq(0, 10, length.out = 101)
  func <- getNativeSymbolInfo("lorenz", PACKAGE = "lorenz")
  y <- rlsoda(y0, times, func, p)
  expect_identical(y, rlsoda(y0, times, "lorenz", p, dllname = "lorenz"))
})

test_that("NULL pointer safety", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 10, length.out = 101)

  func <- getNativeSymbolInfo("lorenz", PACKAGE = "lorenz")
  func <- unserialize(serialize(func, NULL))

  expect_error(rlsoda(y, times, func, p), "null pointer")
  expect_error(rlsoda(y, times, "lorenz", p, output = func, n_out = 2L),
               "null pointer")
})

test_that("time validation", {
  p <- runif(5)
  y <- runif(5)
  times <- seq(0, 1, length.out = 11)
  growth <- function(t, y, p) {
    y * p
  }

  expect_error(rlsoda(y, numeric(0), growth, p),
               "At least two times must be given")
  expect_error(rlsoda(y, c(0, 0), growth, p),
               "Beginning and end times are the same")
  expect_error(rlsoda(y, c(0, 2, 1), growth, p),
               "Times have inconsistent sign")
})
