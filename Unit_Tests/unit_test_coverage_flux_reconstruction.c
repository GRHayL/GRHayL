#include "ghl_unit_tests.h"
#include <float.h>

typedef ghl_error_codes_t (*flux_function)(
      ghl_primitive_quantities *,
      ghl_primitive_quantities *,
      const ghl_eos_parameters *,
      const ghl_metric_quantities *,
      double,
      double,
      ghl_conservative_quantities *);

/* Every rejected speed must fail before EOS dispatch and preserve the output. */
static void check_wave_speed_rejection(void) {
  const flux_function functions[] = {
    ghl_calculate_HLLE_fluxes_dirn0_hybrid_checked,
    ghl_calculate_HLLE_fluxes_dirn1_hybrid_checked,
    ghl_calculate_HLLE_fluxes_dirn2_hybrid_checked,
    ghl_calculate_HLLE_fluxes_dirn0_hybrid_entropy_checked,
    ghl_calculate_HLLE_fluxes_dirn1_hybrid_entropy_checked,
    ghl_calculate_HLLE_fluxes_dirn2_hybrid_entropy_checked,
#ifndef GHL_DISABLE_HDF5
    ghl_calculate_HLLE_fluxes_dirn0_tabulated_checked,
    ghl_calculate_HLLE_fluxes_dirn1_tabulated_checked,
    ghl_calculate_HLLE_fluxes_dirn2_tabulated_checked,
    ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy_checked,
    ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy_checked,
    ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy_checked,
#endif
  };
  const double speeds[][2]
        = { { NAN, 1.0 },           { 1.0, NAN },        { INFINITY, 1.0 },
            { 1.0, INFINITY },      { -1.0, 1.0 },       { 1.0, -1.0 },
            { DBL_MAX, DBL_MAX },   { 0.0, 0.0 },        { 0.5 / DBL_MAX, 0.0 },
            { 0.0, 0.5 / DBL_MAX }, { 1.0e200, 1.0e200 } };
  for(size_t f = 0; f < sizeof(functions) / sizeof(functions[0]); ++f) {
    for(size_t c = 0; c < sizeof(speeds) / sizeof(speeds[0]); ++c) {
      ghl_conservative_quantities result;
      memset(&result, 0xa5, sizeof(result));
      unsigned char original[sizeof(result)];
      memcpy(original, &result, sizeof(result));
      const ghl_error_codes_t error
            = functions[f](NULL, NULL, NULL, NULL, speeds[c][0], speeds[c][1], &result);
      if(error != ghl_error_invalid_hlle_wavespeeds
         || memcmp(original, &result, sizeof(result)) != 0) {
        ghl_error("HLLE rejection changed output: function %zu, case %zu\n", f, c);
      }
    }
  }
}

static void check_maxmod_endpoints(void) {
  const double cases[][3]
        = { { 0.0, 1.0, 0.0 },    { 1.0, 0.0, 0.0 },  { -1.0, 0.0, 0.0 },
            { 0.0, -1.0, 0.0 },   { -1.0, 1.0, 0.0 }, { 1.0, -1.0, 0.0 },
            { 1.0, 2.0, 2.0 },    { 2.0, 1.0, 2.0 },  { -1.0, -2.0, -2.0 },
            { -2.0, -1.0, -2.0 }, { 1.0, 1.0, 1.0 },  { -1.0, -1.0, -1.0 } };
  for(size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); ++c) {
    if(ghl_maxmod(cases[c][0], cases[c][1]) != cases[c][2]) {
      ghl_error("maxmod endpoint %zu failed\n", c);
    }
  }
}

static void check_flattening(void) {
  const ghl_parameters params = { .ppm_flattening_epsilon = 0.33,
                                  .ppm_flattening_omega1 = 0.75,
                                  .ppm_flattening_omega2 = 10.0 };
  const double uniform[5] = { 2.0, 2.0, 2.0, 2.0, 2.0 };
  const double jump[5] = { 1.0, 1.0, 2.0, 3.0, 3.0 };
  const double equal_outer[5] = { 2.0, 1.0, 2.0, 3.0, 2.0 };
  const double compressing[5] = { 0.0, 1.0, 0.0, -1.0, 0.0 };
  const double expanding[5] = { 0.0, -1.0, 0.0, 1.0, 0.0 };
  const double stationary[5] = { 0.0 };
  const double signed_pressure[5] = { -1.0, -1.0, 0.0, 3.0, 3.0 };
  const double zero_outer[5] = { 0.0, 1.0, 2.0, 3.0, 0.0 };
  if(ghl_shock_detection_ftilde(&params, uniform, compressing) != 0.0
     || ghl_shock_detection_ftilde(&params, jump, compressing) != 1.0
     || ghl_shock_detection_ftilde(&params, equal_outer, compressing) != 1.0
     || ghl_shock_detection_ftilde(&params, jump, expanding) != 0.0
     || ghl_shock_detection_ftilde(&params, jump, stationary) != 0.0
     || ghl_shock_detection_ftilde(&params, zero_outer, compressing) != 1.0
     || ghl_shock_detection_ftilde(&params, signed_pressure, expanding) != 0.0
     || ghl_shock_detection_ftilde(&params, signed_pressure, compressing) != 0.0) {
    ghl_error("PPM flattening endpoint failed\n");
  }
}

static void check_steepening(void) {
  ghl_parameters params = { .ppm_shock_k0 = 0.1,
                            .ppm_shock_epsilon = 0.01,
                            .ppm_shock_eta1 = 20.0,
                            .ppm_shock_eta2 = 0.05 };
  const double pressure[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
  const double jump[5] = { 1.0, 1.0, 2.0, 3.0, 3.0 };
  double right = 2.5, left = 1.5;
  ghl_steepen_var(&params, jump, pressure, 2.0, &right, &left);
  if(right != 3.0 || left != 1.0) {
    ghl_error("PPM contact did not steepen to constant endpoint states\n");
  }
  /* A constant contact has no density jump: even with a zero trigger, its
   * zero slopes leave both faces at the constant value without division. */
  const double constant[5] = { 2.0, 2.0, 2.0, 2.0, 2.0 };
  params.ppm_shock_epsilon = 0.0;
  right = left = 2.0;
  ghl_steepen_var(&params, constant, pressure, 2.0, &right, &left);
  if(right != 2.0 || left != 2.0) {
    ghl_error("PPM zero-jump contact changed constant data\n");
  }
  /* Missing/nonpositive normalization inputs must retain the legacy product
   * comparison, not divide by zero. This constant stencil is unchanged. */
  const double zero_pressure[5] = { 0.0 };
  ghl_steepen_var(&params, constant, zero_pressure, 2.0, &right, &left);
  if(right != 2.0 || left != 2.0) {
    ghl_error("PPM zero-pressure contact changed constant data\n");
  }
  ghl_steepen_var(&params, constant, pressure, NAN, &right, &left);
  if(right != 2.0 || left != 2.0) {
    ghl_error("PPM nonfinite contact coefficient changed data\n");
  }
  /* Same-sign curvature forbids steepening, independently of the contact
   * normalization. Exercise the exceptional scaling paths without allowing
   * them to publish nonfinite face values. */
  const struct {
    double rho[5], press[5], gamma;
  } rejected[]
        = { { { 1.0, INFINITY, 2.0, INFINITY, 3.0 }, { 1, 1, 1, 1, 1 }, 2.0 },
            { { 1, 1, 2, 1, 1 }, { 1, INFINITY, 1, INFINITY, 1 }, 2.0 },
            { { 1.0, 1.0, 2.0, INFINITY, INFINITY }, { 1, 1, 1, 1, 1 }, 2.0 },
            { { 1, 1, 2, 1, 1 }, { 1.0, 1.0, 2.0, INFINITY, INFINITY }, 2.0 },
            { { DBL_MIN, DBL_MIN, DBL_MAX, DBL_MAX, INFINITY }, { 1, 1, 1, 1, 1 }, 2.0 },
            { { 1, 1, 2, 1, 1 }, { 1.0, DBL_MIN, DBL_MAX, DBL_MAX, DBL_MAX }, 2.0 },
            { { 1, 1, 4, 3, 3 }, { 1, 1, 1, 1, 1 }, DBL_MAX },
            { { 1, 1, 4, 3, 3 }, { 1, 1, 1, 1, 1 }, INFINITY },
            { { 1, 1, 4, 3, 3 }, { 1, 1, 1, 1, 1 }, 2.0 },
            { { 1, 1, 2, 3, 3 }, { 1, 1, 2, 3, 3 }, 2.0 },
            { { 1, 1, 2, 3, 3 }, { -1, -1, 0, 3, 3 }, 2.0 },
            { { 1, 1, 2, 1, 1 }, { 1, 1, 2, 3, 3 }, 2.0 },
            { { 1.0, 1.0, 2.0, 4.0, 4.0 }, { 1, 1, 1, 1, 1 }, -1.0 } };
  for(size_t c = 0; c < sizeof(rejected) / sizeof(rejected[0]); ++c) {
    right = 2.5;
    left = 1.5;
    ghl_steepen_var(
          &params, rejected[c].rho, rejected[c].press, rejected[c].gamma, &right, &left);
    if(right != 2.5 || left != 1.5) {
      ghl_error("PPM rejected contact %zu changed face data\n", c);
    }
  }
  params.ppm_shock_epsilon = 100.0;
  right = 2.5;
  left = 1.5;
  ghl_steepen_var(&params, jump, pressure, 2.0, &right, &left);
  if(right != 2.5 || left != 1.5) {
    ghl_error("PPM contact below the relative-jump trigger changed faces\n");
  }
}

int main(void) {
  check_wave_speed_rejection();
  check_maxmod_endpoints();
  check_flattening();
  check_steepening();
  ghl_info("Flux/reconstruction coverage boundary tests passed.\n");
  return 0;
}
