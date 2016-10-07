context("rlsoda")

test_that("R", {
  p <- c(10, 28, 8 / 3)
  y0 <- c(10, 1, 1)

  lorenz <- function(t, y, p) {
    sigma <- p[[1L]]
    R <- p[[2L]]
    b <- p[[3L]]
    c(sigma * (y[[2L]] - y[[1L]]),
      R * y[[1L]] - y[[2L]] - y[[1L]] * y[[3L]],
      -b * y[[3L]] + y[[1L]] * y[[2L]])
  }

  tt <- seq(0, 1, length.out = 200)
  res <- lsoda(y0, tt, lorenz, p, deSolve_compatible = TRUE)
  cmp <- deSolve::lsoda(y0, tt, function(...) list(lorenz(...)), p)
  expect_equal(res, unclass(cmp), check.attributes = FALSE, tolerance = 1e-5)
})
