context("rlsoda")

test_that("R", {
  p <- c(10, 28, 8 / 3)
  y0 <- c(10, 1, 1)

  tt <- seq(0, 1, length.out = 200)
  res <- rlsoda(y0, tt, lorenz, p, deSolve_compatible = TRUE)
  cmp <- deSolve::lsoda(y0, tt, function(...) list(lorenz(...)), p)
  expect_equal(res, unclass(cmp), check.attributes = FALSE, tolerance = 1e-5)

  res2 <- rlsoda(y0, tt, lorenz, p,
                n_out = 2L, output = lorenz_output,
                deSolve_compatible = TRUE)

  expect_equal(res2[, 1:4], res)
  expect_equal(res2[, 5:6],
               t(apply(res2[, 2:4], 1, range)),
               check.attributes = FALSE)
})

test_that("C", {
  skip_if_no_compiler()

  p <- c(10, 28, 8 / 3)
  y0 <- c(10, 1, 1)
  tt <- seq(0, 1, length.out = 200)

  res <- rlsoda(y0, tt, "lorenz", p, deSolve_compatible = TRUE, dllname = "lorenz")
  cmp <- rlsoda(y0, tt, lorenz, p, deSolve_compatible = TRUE)

  expect_equal(res, cmp, tolerance = 1e-10)
})

test_that("Critical times", {
  target <- function(t, y, p) {
    if (t <= 1) {
      y
    } else {
      -5 * y
    }
  }
  tt <- seq(0, 2, length.out = 200)
  res1 <- rlsoda(1, tt, target, numeric(0), return_statistics = TRUE)
  res2 <- rlsoda(1, tt, target, numeric(0), tcrit = 1, return_statistics = TRUE)
  res3 <- rlsoda(1, tt, target, numeric(0), tcrit = c(-1, 1),
                 return_statistics = TRUE)

  ## The tolerance is really low here because the non-tcrit version
  ## totally misses the big peak
  expect_equal(as.vector(res1), as.vector(res2), tolerance = 1e-3)

  s1 <- attr(res1, "statistics")
  s2 <- attr(res2, "statistics")
  s3 <- attr(res3, "statistics")
  expect_lt(s2[["n_step"]], s1[["n_step"]])
  expect_equal(s2, s3)
  expect_identical(res2, res3)

  expect_is(attr(res1, "step_size"), "numeric")

  ## A few more related cases:
  res4 <- rlsoda(1, tt, target, numeric(0), tcrit = c(0, 1),
                 return_statistics = TRUE)
  expect_identical(res4, res2)

  res5 <- rlsoda(1, tt, target, numeric(0), tcrit = c(1, 1, 1),
                 return_statistics = TRUE)
  expect_identical(res5, res2)

  ## No valid tcrit:
  res6 <- rlsoda(1, tt, target, numeric(0), tcrit = -1,
                 return_statistics = TRUE)
  expect_identical(res6, res1)
})

test_that("integration failure is caught", {
  p <- c(10, 28, 8 / 3)
  y0 <- c(10, 1, 1)
  tt <- seq(0, 5, length.out = 200)
  ## TODO: Actually too many steps may be taken; using max time here
  ## of 1 will complete!
  res <- rlsoda(y0, tt, "lorenz", p, dllname = "lorenz",
                return_statistics = TRUE)
  expect_error(rlsoda(y0, tt, "lorenz", p, dllname = "lorenz",
                      step_max_n = 10),
               "10 steps taken before reaching tout")
  gc() # because there's an early exit; confirm no crash
})

test_that("SEXP parameters", {
  tt <- seq(0, 10, length.out = 101)
  p <- c(10, 28, 8 / 3)
  y0 <- c(10, 1, 1)
  cmp <- rlsoda(y0, tt, "lorenz", p, dllname = "lorenz")
  res <- rlsoda(y0, tt, "lorenz", p, dllname = "lorenz2",
                parms_are_real = FALSE)
  expect_identical(res, cmp)
})

test_that("pointer parameters", {
  p <- c(10, 28, 8 / 3)
  ptr <- .Call("lorenz_init",  p, PACKAGE = "lorenz3")
  expect_is(ptr, "externalptr")

  tt <- seq(0, 10, length.out = 101)
  y0 <- c(10, 1, 1)

  cmp <- rlsoda(y0, tt, "lorenz", p, dllname = "lorenz")
  res <- rlsoda(y0, tt, "lorenz", ptr, dllname = "lorenz3",
                parms_are_real = FALSE)
  expect_identical(res, cmp)

  ## Argument is optional:
  res <- rlsoda(y0, tt, "lorenz", ptr, dllname = "lorenz3")
  expect_identical(res, cmp)
})
