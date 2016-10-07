context("rlsoda")

test_that("R", {
  p <- c(10, 28, 8 / 3)
  y0 <- c(10, 1, 1)

  tt <- seq(0, 1, length.out = 200)
  res <- lsoda(y0, tt, lorenz, p, deSolve_compatible = TRUE)
  cmp <- deSolve::lsoda(y0, tt, function(...) list(lorenz(...)), p)
  expect_equal(res, unclass(cmp), check.attributes = FALSE, tolerance = 1e-5)

  res2 <- lsoda(y0, tt, lorenz, p,
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

  res <- lsoda(y0, tt, "lorenz", p, deSolve_compatible = TRUE, dllname = "lorenz")
  cmp <- lsoda(y0, tt, lorenz, p, deSolve_compatible = TRUE)

  expect_equal(res, cmp, tolerance = 1e-10)
})
