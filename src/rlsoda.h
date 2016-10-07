#include <R.h>
#include <Rinternals.h>
#include <lsoda.h>

// NOTE: These are different type definitions to those used by
// deSolve, dde, etc.
typedef _lsoda_f deriv_func_ptr;
typedef void (*output_func_ptr)(double t, const double *y,
                                double *out, const void *data);

SEXP rlsoda_ptr_create(struct lsoda_context_t *obj);
void rlsoda_ptr_finalizer(SEXP r_ptr);
int rlsoda_r_harness(double t, double *y, double *dydt, void *data);
void rlsoda_r_output_harness(double t, const double *y,
                             double *out, void *data);
SEXP r_lsoda(SEXP r_y_initial, SEXP r_times, SEXP r_func, SEXP r_data,
             SEXP r_n_out, SEXP r_output,
             SEXP r_data_is_real,
             // Tolerance:
             SEXP r_rtol, SEXP r_atol,
             // Step size control:
             SEXP r_step_size_min, SEXP r_step_size_max,
             SEXP r_step_size_initial, SEXP r_step_max_n,
             // Other:
             SEXP tcrit,
             // Return information:
             SEXP r_return_initial, SEXP r_return_statistics);
