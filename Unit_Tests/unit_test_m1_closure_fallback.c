#include "ghl_m1.h"
#include "m1_test_utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Coverage for the Eulerian Minerbo admissibility fallback in
 * GRHayL/Radiation/ghl_m1_closure.c.
 *
 * Two distinct states publish four_point_compatibility == false: an exact-zero
 * Eulerian flux whose covariant thin dyad vanishes identically, and a
 * finite-flux candidate rejected by the physical PSD check. Both are separately
 * counted by ghl_m1_closure_counters. This test owns three properties:
 *
 *   1. Every published pressure tensor is exactly symmetric, on the fallback
 *      path as well as the primary path. The fallback builds its tensor from a
 *      direction dyad and the inverse metric; without explicit symmetrization
 *      those terms differ in their last bit under an index swap, and the shared
 *      tensor validator then rejects the tensor whenever an off-diagonal
 *      component is small compared with E.
 *   2. The fallback publishes rather than failing. A directional sweep through
 *      the PSD regime returns ghl_success for every state.
 *   3. The PSD regime stays where it is documented. Flux parallel or
 *      antiparallel to the fluid velocity does not reach it; transverse flux
 *      above roughly half light speed does. Pinning both directions makes a
 *      future closure change that widens the regime a test failure rather than
 *      a silent behavior change.
 */

static int failures = 0;

static void fail(const char *const what) {
  fprintf(stderr, "unit_test_m1_closure_fallback: %s\n", what);
  failures++;
}

static void setup_parameters(ghl_m1_parameters *restrict params) {
  if(ghl_m1_initialize(1.0e-8, 1.0e-12, 1.0e-12, 1.0e-6, 1.0e-8, 20, 1.0e-5, params)
     != ghl_success) {
    fail("initializer rejected the campaign control values");
    exit(1);
  }
}

/* Flat metric with zero shift and unit lapse, so prims->vU is the Eulerian
 * three-velocity and the reduced flux factor is the Euclidean |F|/E. */
static void setup_velocity(ghl_primitive_quantities *restrict prims, const double V[3]) {
  memset(prims, 0, sizeof(*prims));
  prims->vU[0] = V[0];
  prims->vU[1] = V[1];
  prims->vU[2] = V[2];
}

static int published_tensor_is_exactly_symmetric(const ghl_m1_closure *restrict closure) {
  for(int i = 0; i < 3; ++i) {
    for(int j = i + 1; j < 3; ++j) {
      if(closure->P[i][j] != closure->P[j][i]) {
        return 0;
      }
    }
  }
  return 1;
}

/* Sweep transverse flux directions at a speed inside the PSD regime. Every
 * state is admissible, so every call must publish, and every published tensor
 * must be exactly symmetric. */
static void check_psd_fallback_publishes_symmetric_tensors(void) {
  ghl_m1_parameters params;
  setup_parameters(&params);
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);

  const double V[3] = { 0.8, 0.0, 0.0 };
  ghl_primitive_quantities prims;
  setup_velocity(&prims, V);

  ghl_m1_reset_closure_counters();

  int sampled = 0;
  int fallbacks = 0;
  for(int step = 0; step < 720; ++step) {
    const double phi = 2.0 * M_PI * (double)step / 720.0;
    for(int k = 1; k <= 9; ++k) {
      const double flux_factor = 0.1 * (double)k;
      ghl_m1_rad_state rad_state;
      rad_state.E = 1.0;
      /* Transverse to V: no x component. */
      rad_state.F[0] = 0.0;
      rad_state.F[1] = flux_factor * cos(phi);
      rad_state.F[2] = flux_factor * sin(phi);

      ghl_m1_closure closure;
      memset(&closure, 0, sizeof(closure));
      const ghl_error_codes_t error
            = ghl_m1_compute_closure_minerbo(&params, &metric, &prims, &rad_state, &closure);
      if(error != ghl_success) {
        fail("admissible transverse-flux state failed to publish a closure");
        return;
      }
      if(!published_tensor_is_exactly_symmetric(&closure)) {
        fail("published pressure tensor is not exactly symmetric");
        return;
      }
      sampled++;
      if(!closure.four_point_compatibility) {
        fallbacks++;
      }
    }
  }

  if(sampled == 0) {
    fail("transverse sweep produced no samples");
    return;
  }
  if(fallbacks == 0) {
    fail("transverse sweep never reached the PSD admissibility fallback; this "
         "test no longer covers the path it owns");
    return;
  }

  ghl_m1_closure_counters counters;
  ghl_m1_get_closure_counters(&counters);
  if(counters.admissibility_fallback_psd != (unsigned long long)fallbacks) {
    fail("admissibility_fallback_psd does not match the observed substitutions");
  }
  if(counters.admissibility_fallback_zero_flux != 0ULL) {
    fail("finite-flux sweep incremented the zero-flux fallback counter");
  }
}

/* An exact-zero Eulerian flux with a moving fluid takes the same fallback for a
 * different, documented reason, and must be accounted separately. */
static void check_zero_flux_fallback_is_counted_separately(void) {
  ghl_m1_parameters params;
  setup_parameters(&params);
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);

  const double V[3] = { 0.3, 0.1, -0.2 };
  ghl_primitive_quantities prims;
  setup_velocity(&prims, V);

  ghl_m1_rad_state rad_state;
  rad_state.E = 1.0;
  rad_state.F[0] = 0.0;
  rad_state.F[1] = 0.0;
  rad_state.F[2] = 0.0;

  ghl_m1_reset_closure_counters();

  ghl_m1_closure closure;
  memset(&closure, 0, sizeof(closure));
  if(ghl_m1_compute_closure_minerbo(&params, &metric, &prims, &rad_state, &closure)
     != ghl_success) {
    fail("exact-zero-flux state failed to publish a closure");
    return;
  }
  if(!published_tensor_is_exactly_symmetric(&closure)) {
    fail("zero-flux published pressure tensor is not exactly symmetric");
  }
  if(closure.four_point_compatibility) {
    fail("exact-zero-flux state with a moving fluid did not take the fallback; "
         "this test no longer covers the path it owns");
    return;
  }

  ghl_m1_closure_counters counters;
  ghl_m1_get_closure_counters(&counters);
  if(counters.admissibility_fallback_zero_flux != 1ULL) {
    fail("zero-flux fallback was not counted");
  }
  if(counters.admissibility_fallback_psd != 0ULL) {
    fail("zero-flux fallback incremented the PSD fallback counter");
  }
}

/* Pin the documented regime boundary in both directions. */
static void check_psd_regime_boundary(void) {
  ghl_m1_parameters params;
  setup_parameters(&params);
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);

  const double flux_factors[] = { 0.05, 0.2, 0.4, 0.6, 0.8, 0.95 };
  const int flux_count = (int)(sizeof(flux_factors) / sizeof(flux_factors[0]));

  for(int speed_step = 1; speed_step <= 9; ++speed_step) {
    const double speed = 0.1 * (double)speed_step;
    const double V[3] = { speed, 0.0, 0.0 };
    ghl_primitive_quantities prims;
    setup_velocity(&prims, V);

    for(int k = 0; k < flux_count; ++k) {
      const double flux_factor = flux_factors[k];

      /* Flux aligned with and opposed to the flow must stay on the primary
       * construction at every speed tested. */
      for(int sign = -1; sign <= 1; sign += 2) {
        ghl_m1_rad_state rad_state;
        rad_state.E = 1.0;
        rad_state.F[0] = (double)sign * flux_factor;
        rad_state.F[1] = 0.0;
        rad_state.F[2] = 0.0;
        ghl_m1_closure closure;
        memset(&closure, 0, sizeof(closure));
        if(ghl_m1_compute_closure_minerbo(&params, &metric, &prims, &rad_state, &closure)
           != ghl_success) {
          fail("aligned-flux state failed to publish a closure");
          return;
        }
        if(!published_tensor_is_exactly_symmetric(&closure)) {
          fail("aligned-flux published pressure tensor is not exactly symmetric");
          return;
        }
        if(!closure.four_point_compatibility) {
          fail("flux aligned with the fluid velocity reached the PSD fallback; "
               "the documented regime has widened");
          return;
        }
      }

      /* Transverse flux: below the documented onset the primary construction
       * must still be published. */
      if(speed <= 0.4) {
        ghl_m1_rad_state rad_state;
        rad_state.E = 1.0;
        rad_state.F[0] = 0.0;
        rad_state.F[1] = flux_factor;
        rad_state.F[2] = 0.0;
        ghl_m1_closure closure;
        memset(&closure, 0, sizeof(closure));
        if(ghl_m1_compute_closure_minerbo(&params, &metric, &prims, &rad_state, &closure)
           != ghl_success) {
          fail("transverse-flux state below the onset failed to publish a closure");
          return;
        }
        if(!closure.four_point_compatibility) {
          fail("transverse flux below half light speed reached the PSD fallback; "
               "the documented regime has widened");
          return;
        }
      }
    }
  }
}

/* Crossing the E^2 overflow boundary must not change the dimensionless
 * closure. A nonidentity metric distinguishes F_i V^i from F^i V^i. */
static void check_large_energy_scaling(void) {
  ghl_m1_parameters params;
  setup_parameters(&params);
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  metric.gammaDD[0][0] = 4.0;
  metric.gammaUU[0][0] = 0.25;
  metric.detgamma = 4.0;
  metric.sqrt_detgamma = 2.0;
  const double V[3] = { 0.1, 0.0, 0.0 };
  ghl_primitive_quantities prims;
  setup_velocity(&prims, V);

  /* 2^520 has a finite pressure tensor but its square exceeds DBL_MAX. */
  const double energies[] = { 1.0, 0x1p520 };
  ghl_m1_closure closures[2];
  for(int k = 0; k < 2; ++k) {
    const ghl_m1_rad_state state
          = { .E = energies[k], .F = { 0.8 * energies[k], 0.0, 0.0 } };
    if(ghl_m1_compute_closure_with_primitives(
             &params, &metric, &prims, &state, &closures[k]) != ghl_success) {
      fail("energy-scaled closure failed to publish");
      return;
    }
  }
  const double tolerance = params.closure_root_residual_tolerance;
  if(fabs(closures[0].xi - closures[1].xi) > tolerance
     || fabs(closures[0].chi - closures[1].chi) > tolerance) {
    fail("dimensionless closure changed across the E^2 overflow boundary");
  }
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      if(fabs(closures[0].P[i][j] - closures[1].P[i][j] / energies[1])
         > tolerance) {
        fail("normalized pressure changed across the E^2 overflow boundary");
      }
    }
  }
}

int main(void) {
  check_large_energy_scaling();
  check_psd_fallback_publishes_symmetric_tensors();
  check_zero_flux_fallback_is_counted_separately();
  check_psd_regime_boundary();

  if(failures != 0) {
    fprintf(stderr, "unit_test_m1_closure_fallback: %d check(s) failed\n", failures);
    return 1;
  }

  ghl_info("unit_test_m1_closure_fallback: admissibility fallback symmetry, "
           "accounting, and regime checks passed\n");
  return 0;
}
