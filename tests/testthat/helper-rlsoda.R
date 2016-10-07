lorenz <- function(t, y, p) {
  sigma <- p[[1L]]
  R <- p[[2L]]
  b <- p[[3L]]
  c(sigma * (y[[2L]] - y[[1L]]),
    R * y[[1L]] - y[[2L]] - y[[1L]] * y[[3L]],
    -b * y[[3L]] + y[[1L]] * y[[2L]])
}
lorenz_output <- function(t, y, p) {
  c(min(y), max(y))
}

.first_time <- TRUE
.loaded_dlls <- FALSE
prepare_all <- function(reload = FALSE) {
  if (.first_time || reload) {
    cleanup_objects()
    .first_time <<- FALSE
  }
  if (requireNamespace("rcmdshlib", quietly = TRUE)) {
    files <- dir(pattern = "\\.c$")
    if (!reload) {
      files <- files[!(sub("\\.c$", "", files) %in% names(getLoadedDLLs()))]
    }
    for (f in files) {
      message("Compiling ", f)
      dyn.load(rcmdshlib::shlib(f, verbose = FALSE)$dll)
    }
    .loaded_dlls <<- TRUE
  } else {
    .loaded_dlls <<- FALSE
  }
}

cleanup_objects <- function() {
  files <- dir(pattern="\\.(o|so|dll)$")
  if (length(files) > 0L) {
    file.remove(files)
  }
  invisible()
}

try(prepare_all(), silent = !interactive())

skip_if_no_compiler <- function() {
  skip_if_not_installed("rcmdshlib")
  if (.loaded_dlls) {
    return()
  }
  testthat::skip("could not compile libraries")
}
