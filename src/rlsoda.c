#include "rlsoda.h"
#include <common.h>
#include "util.h"
#include "stdbool.h"

// This is _heavily_ modelled off of the same interface in dde.  If
// the interpolation matrix can be exploited in lsoda then we'd be
// able to justify moving the two together, too.
SEXP r_lsoda(SEXP r_y_initial, SEXP r_times, SEXP r_func, SEXP r_data,
             SEXP r_n_out, SEXP r_output,
             SEXP r_data_is_real,
             // Tolerance:
             SEXP r_rtol, SEXP r_atol,
             // Step size control:
             SEXP r_step_size_min, SEXP r_step_size_max,
             SEXP r_step_size_initial, SEXP r_step_max_n,
             // Other:
             SEXP r_tcrit,
             // Return information:
             SEXP r_return_initial, SEXP r_return_statistics) {

  double *y_initial = REAL(r_y_initial);
  size_t n = length(r_y_initial);

  size_t n_times = LENGTH(r_times);
  double *times = REAL(r_times);

  if (n_times < 2) {
    Rf_error("At least two times must be given");
  }
  if (times[n_times - 1] == times[0]) {
    Rf_error("Initialisation failure: Beginning and end times are the same");
  }
  const double sign = copysign(1.0, times[n_times - 1] - times[0]);
  for (size_t i = 0; i < n_times - 1; ++i) {
    double t0 = times[i], t1 = times[i + 1];
    bool err = (sign > 0 && t1 < t0) || (sign < 0 && t1 > t0);
    if (err) {
      Rf_error("Initialisation failure: Times have inconsistent sign");
    }
  }

  size_t n_tcrit = 0;
  double *tcrit = NULL;
  if (r_tcrit != R_NilValue) {
    n_tcrit = LENGTH(r_tcrit);
    tcrit = REAL(r_tcrit);
  }

  deriv_func_ptr func = (deriv_func_ptr)R_ExternalPtrAddr(r_func);
  if (func == NULL) {
    Rf_error("Was passed null pointer for 'func'");
  }
  void *data = NULL;
  if (TYPEOF(r_data) == REALSXP && INTEGER(r_data_is_real)[0]) {
    data = (void*) REAL(r_data);
  } else if (TYPEOF(r_data) == EXTPTRSXP) {
    data = R_ExternalPtrAddr(r_data);
  } else {
    data = (void*) r_data;
  }

  bool return_initial = INTEGER(r_return_initial)[0];
  bool return_statistics = INTEGER(r_return_statistics)[0];
  size_t nt = return_initial ? n_times : n_times - 1;

  size_t n_out = INTEGER(r_n_out)[0];
  const bool has_output = n_out > 0;
  output_func_ptr output = NULL;
  double *out = NULL;
  SEXP r_out = R_NilValue;
  if (has_output) {
    output = (output_func_ptr)R_ExternalPtrAddr(r_output);
    if (output == NULL) {
      Rf_error("Was passed null pointer for 'output'");
    }
    r_out = PROTECT(allocMatrix(REALSXP, n_out, nt));
    out = REAL(r_out);
  }

  struct lsoda_opt_t opt = {0};
  opt.ixpr = 0;
  opt.itask = 1;
  opt.rtol = REAL(r_rtol);
  opt.atol = REAL(r_atol);
  opt.mxstep = INTEGER(r_step_max_n)[0];
  // These need careful treatment on default.
  opt.h0 = REAL(r_step_size_initial)[0];
  opt.hmin = fmax(fabs(REAL(r_step_size_min)[0]), DBL_EPSILON);
  opt.hmin = fmin(fabs(REAL(r_step_size_max)[0]), DBL_MAX);
  // Not setting here:
  // * mxhnil
  // * mxordn
  // * mxords
  // * hmxi

  struct lsoda_context_t ctx = {
    .function = func,
    .neq = n,
    .data = data,
    .state = 1,
  };

  SEXP ptr = PROTECT(rlsoda_ptr_create(&ctx));

  SEXP r_y = PROTECT(allocMatrix(REALSXP, n, nt));
  double *y = REAL(r_y);

  double *yi = (double*) R_alloc(n, sizeof(double));
  memcpy(yi, y_initial, n * sizeof(double));
  double t = times[0];

  size_t tcrit_idx = 0;
  bool has_tcrit = n_tcrit > 0;
  if (has_tcrit) {
    while (sign * tcrit[tcrit_idx] <= t && tcrit_idx < n_tcrit) {
      tcrit_idx++;
    }
    if (tcrit_idx < n_tcrit) {
      opt.tcrit = tcrit[tcrit_idx];
      opt.itask = 4;
    } else {
      has_tcrit = false;
    }
  }

  lsoda_prepare(&ctx, &opt);
  if (return_initial) {
    memcpy(y, yi, n * sizeof(double));
    y += n;
    if (has_output) {
      output(t, yi, out, data);
      out += n_out;
    }
  }

  double t_next = times[1];
  bool store_next = true;

  // This is *really* nasty book-keeping, but it should work OK.  This
  // interleaves the critical times with the storage times, organises
  // stopping at each of them and recording the output appropriately.
  for (size_t i = 1; i < n_times;) {
    lsoda(&ctx, yi, &t, t_next);
    if (ctx.error) {
      Rf_error("Integration error: %s", ctx.error);
    }
    if (store_next) {
      memcpy(y, yi, n * sizeof(double));
      y += n;
      if (has_output) {
        output(t, yi, out, data);
        out += n_out;
      }
      ++i;
      if (i < n_times) {
        t_next = times[i];
      }
    } else {
      store_next = true;
    }
    if (has_tcrit) {
      while (sign * tcrit[tcrit_idx] <= sign * t) {
        tcrit_idx++;
      }
      if (tcrit_idx < n_tcrit) {
        ctx.opt->tcrit = tcrit[tcrit_idx];
      } else {
        has_tcrit = false;
        ctx.opt->itask = 1;
      }
      if (sign * tcrit[tcrit_idx] <= sign * t_next) {
        t_next = tcrit[tcrit_idx];
        store_next = !(sign * tcrit[tcrit_idx] < sign * t_next);
      } else {
        store_next = true;
      }
    }
  }

  if (has_output) {
    setAttrib(r_y, install("output"), r_out);
    UNPROTECT(1);
  }

  if (return_statistics) {
    SEXP stats = PROTECT(allocVector(INTSXP, 5));
    SEXP stats_nms = PROTECT(allocVector(STRSXP, 5));
    INTEGER(stats)[0] = ctx.common->nfe;
    SET_STRING_ELT(stats_nms, 0, mkChar("n_eval"));
    INTEGER(stats)[1] = ctx.common->nst;
    SET_STRING_ELT(stats_nms, 1, mkChar("n_step"));
    INTEGER(stats)[2] = ctx.common->nje;
    SET_STRING_ELT(stats_nms, 2, mkChar("n_jacobian"));
    INTEGER(stats)[3] = ctx.common->nqu;
    SET_STRING_ELT(stats_nms, 3, mkChar("last_order"));
    INTEGER(stats)[4] = ctx.common->nq;
    SET_STRING_ELT(stats_nms, 4, mkChar("next_order"));
    setAttrib(stats, R_NamesSymbol, stats_nms);
    setAttrib(r_y, install("statistics"), stats);
    setAttrib(r_y, install("step_size"), ScalarReal(ctx.common->h));
    UNPROTECT(2);
  }

  // Deterministically clean up if we can, otherwise we clean up by R
  // running the finaliser for us when it garbage collects ptr above.
  rlsoda_ptr_finalizer(ptr);
  UNPROTECT(2);
  return r_y;
}

int rlsoda_r_harness(double t, double *y, double *dydt, void *data) {
  SEXP d = (SEXP)data;
  SEXP
    target = VECTOR_ELT(d, 0),
    parms = VECTOR_ELT(d, 1),
    rho = VECTOR_ELT(d, 2);
  const int n = INTEGER(VECTOR_ELT(d, 3))[0]; // NOTE: different to dde
  SEXP r_t = PROTECT(ScalarReal(t));
  SEXP r_y = PROTECT(allocVector(REALSXP, n));
  memcpy(REAL(r_y), y, n * sizeof(double));
  SEXP call = PROTECT(lang4(target, r_t, r_y, parms));
  SEXP ans = PROTECT(eval(call, rho));
  memcpy(dydt, REAL(ans), n * sizeof(double));
  UNPROTECT(4);
  return 0;
}

void rlsoda_r_output_harness(double t, const double *y,
                             double *out, void *data) {
  SEXP d = (SEXP)data;
  SEXP
    parms = VECTOR_ELT(d, 1),
    rho = VECTOR_ELT(d, 2),
    output = VECTOR_ELT(d, 4);
  const int
    n = INTEGER(VECTOR_ELT(d, 3))[0], // NOTE: different to dde
    n_out = INTEGER(VECTOR_ELT(d, 5))[0]; // NOTE: different to dde
  SEXP r_t = PROTECT(ScalarReal(t));
  SEXP r_y = PROTECT(allocVector(REALSXP, n));
  memcpy(REAL(r_y), y, n * sizeof(double));
  SEXP call = PROTECT(lang4(output, r_t, r_y, parms));
  SEXP ans = PROTECT(eval(call, rho));
  memcpy(out, REAL(ans), n_out * sizeof(double));
  UNPROTECT(4);
}

void rlsoda_ptr_finalizer(SEXP r_ptr) {
  void *obj = R_ExternalPtrAddr(r_ptr);
  if (obj) {
    lsoda_free((struct lsoda_context_t*) obj);
    R_ClearExternalPtr(r_ptr);
  }
}

SEXP rlsoda_ptr_create(struct lsoda_context_t *obj) {
  SEXP r_ptr = PROTECT(R_MakeExternalPtr(obj, R_NilValue, R_NilValue));
  R_RegisterCFinalizer(r_ptr, rlsoda_ptr_finalizer);
  UNPROTECT(1);
  return r_ptr;
}
