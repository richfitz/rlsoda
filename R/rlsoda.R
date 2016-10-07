##' Integrate an ODE with lsoda.
##'
##' @title Integrate ODE with lsoda
##'
##' @param y Initial conditions for the integration
##'
##' @param times Times where output is needed.  Unlike \code{deSolve}
##'   we won't actually stop at these times, but instead interpolate
##'   back to get the result.
##'
##' @param func Function to integrate.  Can be an R function of
##'   arguments \code{t, y, parms}, returning a numeric vector, or it
##'   can be the name or address of a C function with arguments
##'   \code{size_t n, double t, const double *y, double *dydt, void *data}.
##'
##' @param parms Parameters to pass through to the derivatives.
##'
##' @param ... Dummy arguments - nothing is allowed here, but this
##'   means that all further arguments \emph{must} be specified by
##'   name (not order) so I can easily reorder them later on.
##'
##' @param n_out Number of "output" variables (not differential
##'   equation variables) to compute via the routine \code{output}.
##'
##' @param output The output routine; either an R function taking
##'   arguments \code{t, y, parms} or the name/address of a C function
##'   taking arguments \code{size_t n, double t, const double *y,
##'   size_t n_out, double *out, void *data}.
##'
##' @param rtol The per-step relative tolerance.  The total accuracy
##'   will be less than this.
##'
##' @param atol The per-step absolute tolerance.
##'
##' @param step_size_min The minimum step size.  The actual minimum
##'   used will be the largest of the absolute value of this
##'   \code{step_size_min} or \code{.Machine$double.eps}.  If the
##'   integration attempts to make a step smaller than this, it will
##'   throw an error, stopping the integration (note that this differs
##'   from the treatment of \code{hmin} in \code{deSolve::lsoda}).
##'
##' @param step_size_max The largest step size.  By default there is
##'   no maximum step size (Inf) so the solver can take as large a
##'   step as it wants to.  If you have short events you want the
##'   solver to notice, then specify a smaller maximim step size here
##'   (or use \code{tcrit} below).
##'
##' @param step_size_initial The initial step size.  By default the
##'   integrator will guess the step size automatically, but one can
##'   be given here instead.
##'
##' @param step_max_n The maximum number of steps allowed.  If the
##'   solver takes more steps than this it will throw an error.  Note
##'   the number of evaluations of \code{func} will be a multiple of
##'   the number of steps.
##'
##' @param tcrit An optional vector of critical times that the solver
##'   must stop at (rather than interpolating over).  This can include
##'   an end time that we can't go past, or points within the
##'   integration that must be stopped at exactly (for example cases
##'   where the derivatives change abruptly).  Note that this differs
##'   from the interpretation of this parameter in deSolve; there
##'   \code{tcrit} is a single time that integration may not go past
##'   -- with dde we never go past the final time, and this is just
##'   for times that fall \emph{within} the range of times in
##'   \code{times}.
##'
##' @param dllname Name of the shared library (without extension) to
##'   find the function \code{func} (and \code{output} if given) in
##'   the case where \code{func} refers to compiled function.
##'
##' @param parms_are_real Logical, indicating if \code{parms} should
##'   be treated as vector of doubles by \code{func} (when it is a
##'   compiled function).  If \code{TRUE} (the default), then
##'   \code{REAL(parms)}, which is \code{double*} is passed through.
##'   If \code{FALSE} then if \code{params} is an externalptr type
##'   (\code{EXTPTRSXP}) we pass through the result of
##'   \code{R_ExternalPtrAddr}, otherwise we pass \code{params}
##'   through unmodified as a \code{SEXP}.  In the last case, in your
##'   target function you will need to include \code{<Rinternals.h>},
##'   cast to \code{SEXP} and then pull it apart using the R API (or
##'   Rcpp).
##'
##' @param ynames Logical, indicating if the output should be named
##'   following the names of the input vector \code{y}.
##'   Alternatively, if \code{ynames} is a character vector of the
##'   same length as \code{y}, these will be used as the output names.
##'
##' @param outnames An optional character vector, used when
##'   \code{n_out} is greater than 0, to name the model output matrix.
##'
##' @param by_column Logical, indicating if the output should be
##'   returned organised by column (rather than row).  This incurs a
##'   slight cost for transposing the matrices.  If you can work with
##'   matrices that are transposed relative to \code{deSolve}, then
##'   set this to \code{FALSE}.
##'
##' @param return_initial Logical, indicating if the output should
##'   include the initial conditions (like deSolve).
##'
##' @param return_statistics Logical, indicating if statistics about
##'   the run should be included.  If \code{TRUE}, then an integer
##'   vector containing the number of target evaluations, steps,
##'   accepted steps and rejected steps is returned (the vector is
##'   named).
##'
##' @param return_time Logical, indicating if a row (or column if
##'   \code{by_column} is \code{TRUE}) representing time is included
##'   (this matches deSolve).
##'
##' @param return_output_with_y Logical, indicating if the output
##'   should be bound together with the returned matrix \code{y} (as
##'   it is with \code{deSolve}).  Otherwise output will be returned
##'   as the attribute \code{output}.
##'
##' @param deSolve_compatible Logical, indicating if we should run in
##'   "deSolve compatible" output mode.  This enables the options
##'   \code{by_column}, \code{return_initial}, \code{return_time} and
##'   \code{return_output_with_y}.  This affects only some aspects of
##'   the returned value, and not the calculations themselves.
##'
##' @return At present the return value is transposed relative to
##'   deSolve.  This might change in future.
##'
##' @export
rlsoda <- function(y, times, func, parms, ...,
                  n_out = 0L, output = NULL,
                  rtol = 1e-6, atol = 1e-6,
                  step_size_min = 0, step_size_max = 0,
                  step_size_initial = 0, step_max_n = 100000L,
                  tcrit = NULL,
                  dllname = "",
                  parms_are_real = TRUE,
                  ynames = TRUE,
                  outnames = NULL,
                  by_column = FALSE, return_initial = FALSE,
                  return_statistics = FALSE, return_time = FALSE,
                  return_output_with_y = FALSE,
                  deSolve_compatible = FALSE) {
  DOTS <- list(...)
  if (length(DOTS) > 0L) {
    stop("Invalid dot arguments!")
  }
  if (deSolve_compatible) {
    by_column <- TRUE
    return_initial <- TRUE
    return_time <- TRUE
    return_output_with_y <- TRUE
  }

  func <- find_function_address(func, dllname)

  is_r_target <- is.function(func)

  if (is_r_target) {
    parms_are_real <- FALSE
    parms <- list(func = func, parms = parms, rho = parent.frame(),
                  n = as.integer(length(y)))
    func <- find_function_address("rlsoda_r_harness", "rlsoda")
    if (nzchar(dllname)) {
      stop("dllname must not be given when using an R function for 'func'")
    }
  }

  assert_scalar(rtol)
  assert_scalar(atol)
  assert_scalar(step_size_min)
  assert_scalar(step_size_max)
  assert_scalar(step_size_initial)
  assert_size(step_max_n)
  assert_scalar_logical(parms_are_real)
  assert_scalar_logical(by_column)
  assert_scalar_logical(return_initial)
  assert_scalar_logical(return_statistics)
  assert_scalar_logical(return_time)
  assert_scalar_logical(return_output_with_y)

  ## Needed here...
  atol <- rep_len(atol, length(y))
  rtol <- rep_len(rtol, length(y))

  ynames <- check_ynames(y, ynames, deSolve_compatible)

  assert_size(n_out)
  has_output <- n_out > 0L
  outnames <- check_outnames(n_out, outnames)

  if (has_output) {
    output <- find_function_address(output, dllname)
    ## NOTE: The same-typedness of output/func is not really
    ## necessary, but I think it's simplest to think about if we
    ## enforce it.  We should be able to put anything into a harness
    ## either way.
    if (is_r_target) {
      if (!is.function(output)) {
        stop("output must be an R function")
      }
      parms <- c(parms, list(output, n_out))
      output <- find_function_address("rlsoda_r_output_harness", "rlsoda")
    } else {
      if (is.function(output)) {
        stop("output must be a compiled function (name or address)")
      }
    }
  } else if (!is.null(output)) {
    stop("If 'output' is given, then n_out must be specified")
  }

  ret <- .Call(Clsoda, y, as.numeric(times), func, parms,
               as.integer(n_out), output,
               parms_are_real,
               ## Tolerance:
               rtol, atol,
               ## Step control:
               step_size_min, step_size_max,
               step_size_initial, as.integer(step_max_n),
               ## Other:
               tcrit,
               ## Return information:
               return_initial, return_statistics)

  prepare_output(ret, times, ynames, outnames, has_output,
                 by_column, return_initial, return_time,
                 return_output_with_y,
                 "time")
}
