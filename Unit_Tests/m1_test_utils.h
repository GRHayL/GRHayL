#ifndef UNIT_TESTS_M1_TEST_UTILS_H_
#define UNIT_TESTS_M1_TEST_UTILS_H_

#include "ghl_m1.h"
#include "ghl_unit_tests.h"
#include <float.h>
#include <math.h>

static inline void m1_setup_flat_metric(ghl_metric_quantities *restrict metric) {

  *metric = (ghl_metric_quantities){ 0 };
  metric->lapse = 1.0;
  metric->lapseinv = 1.0;
  metric->lapseinv2 = 1.0;
  metric->detgamma = 1.0;
  metric->sqrt_detgamma = 1.0;

  metric->gammaDD[0][0] = 1.0;
  metric->gammaDD[1][1] = 1.0;
  metric->gammaDD[2][2] = 1.0;
  metric->gammaUU[0][0] = 1.0;
  metric->gammaUU[1][1] = 1.0;
  metric->gammaUU[2][2] = 1.0;
}

static inline void m1_setup_metric_anchor_B(ghl_metric_quantities *restrict metric) {

  *metric = (ghl_metric_quantities){ 0 };
  metric->lapse = 0.91;
  metric->lapseinv = 1.0 / metric->lapse;
  metric->lapseinv2 = 1.0 / (metric->lapse * metric->lapse);
  metric->betaU[0] = 0.09;
  metric->betaU[1] = -0.06;
  metric->betaU[2] = 0.02;

  metric->gammaDD[0][0] = 1.00;
  metric->gammaDD[0][1] = 0.05;
  metric->gammaDD[0][2] = 0.00;
  metric->gammaDD[1][0] = 0.05;
  metric->gammaDD[1][1] = 1.15;
  metric->gammaDD[1][2] = -0.04;
  metric->gammaDD[2][0] = 0.00;
  metric->gammaDD[2][1] = -0.04;
  metric->gammaDD[2][2] = 0.93;

  metric->gammaUU[0][0] = 1.0021819205593225;
  metric->gammaUU[0][1] = -0.0436384111864486;
  metric->gammaUU[0][2] = -0.0018769209112451;
  metric->gammaUU[1][0] = -0.0436384111864486;
  metric->gammaUU[1][1] = 0.8727682237289727;
  metric->gammaUU[1][2] = 0.0375384182249020;
  metric->gammaUU[2][0] = -0.0018769209112451;
  metric->gammaUU[2][1] = 0.0375384182249020;
  metric->gammaUU[2][2] = 1.0768833728268774;

  metric->detgamma = 1.065575;
  metric->sqrt_detgamma = sqrt(metric->detgamma);
}

static inline void m1_permute_metric_xy(
      const ghl_metric_quantities *restrict input,
      ghl_metric_quantities *restrict output) {

  static const int p[3] = { 1, 0, 2 };

  *output = (ghl_metric_quantities){ 0 };
  output->lapse = input->lapse;
  output->lapseinv = input->lapseinv;
  output->lapseinv2 = input->lapseinv2;
  output->detgamma = input->detgamma;
  output->sqrt_detgamma = input->sqrt_detgamma;
  for(int i = 0; i < 3; i++) {
    output->betaU[i] = input->betaU[p[i]];
    for(int j = 0; j < 3; j++) {
      output->gammaDD[i][j] = input->gammaDD[p[i]][p[j]];
      output->gammaUU[i][j] = input->gammaUU[p[i]][p[j]];
    }
  }
}

static inline void m1_permute_state_xy(
      const ghl_m1_rad_state *restrict input,
      ghl_m1_rad_state *restrict output) {

  *output = *input;
  output->F[0] = input->F[1];
  output->F[1] = input->F[0];
}

static inline int m1_nearly_equal(
      const double a,
      const double b,
      const double rel_tol,
      const double abs_tol) {

  const double scale = fmax(1.0, fmax(fabs(a), fabs(b)));
  return fabs(a - b) <= fmax(abs_tol, rel_tol * scale);
}

static inline double m1_invariant_slack(const double scale) {

  return 128.0 * DBL_EPSILON * fmax(1.0, fabs(scale));
}

/* Compare two published closures by their named members.
 *
 * ghl_m1_closure mixes doubles with int/enum/bool, so it carries tail padding
 * whose bytes C does not define across a struct assignment (C11 6.2.6.1p6).
 * memcmp over the whole object therefore tests representation equality, not
 * the transactional contract these tests mean to assert: that no member was
 * published.  Compare the members instead.  Doubles are compared exactly on
 * purpose -- an unchanged closure must be bit-for-bit unchanged. */
static inline int m1_closure_identical(
      const ghl_m1_closure *restrict a,
      const ghl_m1_closure *restrict b) {

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      if(a->P[i][j] != b->P[i][j]) {
        return 0;
      }
    }
  }
  return a->chi == b->chi && a->xi == b->xi && a->root_residual == b->root_residual
         && a->root_iterations == b->root_iterations
         && a->solve_status == b->solve_status
         && a->four_point_compatibility == b->four_point_compatibility;
}

#endif // UNIT_TESTS_M1_TEST_UTILS_H_
