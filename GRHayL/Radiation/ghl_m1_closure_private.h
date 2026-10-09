#ifndef GHL_M1_CLOSURE_PRIVATE_H
#define GHL_M1_CLOSURE_PRIVATE_H

#include "ghl_m1.h"

/* These declarations are shared only with the closure unit test. They are not
 * installed and are not part of the public M1 API. */
typedef struct {
  double P[3][3];
  double chi;
  double residual;
  double normalized_residual;
  double physical_xi;
} ghl_m1_closure_evaluation;

typedef struct {
  double xi;
  int iterations;
  ghl_m1_closure_solve_status_t status;
  ghl_m1_closure_evaluation evaluation;
} ghl_m1_closure_root_result;

typedef ghl_error_codes_t (*ghl_m1_closure_evaluator)(
      void *context,
      double xi,
      ghl_m1_closure_evaluation *evaluation);

ghl_error_codes_t ghl_m1_closure_private_evaluate_invariant(
      double J,
      double H2,
      double scale,
      double xi,
      double *H2_clipped,
      double *residual,
      double *normalized_residual,
      double *physical_xi);

ghl_error_codes_t ghl_m1_closure_private_finish_evaluation(
      const ghl_metric_quantities *metric,
      const double PDD[4][4],
      double energy_scale,
      double W,
      double J,
      double H2,
      double scale,
      double xi,
      double chi,
      ghl_m1_closure_evaluation *evaluation);

ghl_error_codes_t ghl_m1_closure_private_finalize_fallback(
      ghl_error_codes_t conversion_error,
      const ghl_m1_comoving *comoving,
      const ghl_m1_closure *candidate,
      ghl_m1_closure *closure);

ghl_error_codes_t ghl_m1_closure_private_construct_eulerian_pressure(
      double energy,
      const double gammaUU[3][3],
      double flux_factor,
      const double direction[3],
      double P[3][3],
      double *chi);

ghl_error_codes_t ghl_m1_closure_private_check_residual_gate(
      double normalized_residual,
      double residual_tolerance);

ghl_error_codes_t ghl_m1_closure_private_solve_root(
      const ghl_m1_parameters *params,
      ghl_m1_closure_evaluator evaluate,
      void *context,
      ghl_m1_closure_root_result *root);

#endif // GHL_M1_CLOSURE_PRIVATE_H
