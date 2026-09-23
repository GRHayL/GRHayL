#ifndef GHL_M1_NEUTRINO_IMPLICIT_H
#define GHL_M1_NEUTRINO_IMPLICIT_H

/*
 * Generic private Radiation header for grey three-species neutrino M1
 * number-current, interaction-source, and implicit-solve helpers.
 *
 * This header declares helpers used by
 * ghl_m1_solve_neutrino_implicit_homogeneous_update across translation units.
 * Public residual, Jacobian, trial-state, and admissibility declarations live
 * in ghl_m1.h; the explicit-base variants and remaining helpers here remain
 * internal:
 *
 *   - The neutrino residual and Jacobian use frozen rates and frozen
 *     primitives.
 *   - The neutrino frozen-rate solve uses the shared Newton driver but keeps
 *     its rate residual, Jacobian, and solver helpers in this subdirectory.
 *
 * Its scope includes the shared number-current helper, interaction sources,
 * and the local implicit solve. It is listed in the neutrino build manifest's
 * private include line (#! INCS =).
 */

#include "ghl_m1.h"

/* The scaled-positive helpers below use frexp, scalbn, sqrt, and isfinite
 * directly, so this header states that dependency rather than relying on a
 * transitive include. */
#include <math.h>

/* Public single-species source boundaries additionally require that electron
 * rates contain no reaction requiring the missing partner state. */
ghl_error_codes_t ghl_m1_neutrino_validate_single_species_rates(
      const ghl_m1_neutrino_rates *restrict rates,
      ghl_m1_neutrino_diagnostics *restrict diagnostics);

typedef struct {
  double J;
  double h_n;
  double HU[3];
  double Gamma_N;
  double n_com;
  double number_flux[3];
  double number_transport_velocity[3];
} ghl_m1_neutrino_current;

/* Immutable inputs shared by all Newton residual and Jacobian evaluations in
 * one neutrino solve. The solve entry point validates this context before it
 * is handed to the private callbacks. */
typedef struct {
  const ghl_m1_parameters *restrict m1_params;
  const ghl_metric_quantities *restrict metric;
  const ghl_primitive_quantities *restrict prims_frozen;
  const ghl_m1_neutrino_rates *restrict rates;
} ghl_m1_neutrino_implicit_context;

/* Build the number current from moments and velocity quantities that have
 * already been evaluated for the same radiation state. */
ghl_error_codes_t ghl_m1_neutrino_build_current_from_moments(
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_comoving *restrict comoving,
      const double V_con[3],
      const double W,
      ghl_m1_neutrino_current *restrict current);

/* Radiation-private transactional exchange assembly. */
ghl_error_codes_t ghl_m1_neutrino_assemble_exchange(
      const ghl_m1_neutrino_state *restrict state_in,
      const ghl_m1_neutrino_state *restrict state_out,
      const ghl_m1_neutrino_rates *restrict rates,
      double dL_rad_cc,
      double sqrt_detgamma,
      double baryon_density_conserved,
      ghl_m1_neutrino_exchange *restrict exchange);

ghl_error_codes_t ghl_m1_neutrino_derive_current(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      ghl_m1_neutrino_current *restrict current);

/* Apply the shared endpoint-number policy used by source shortcuts and the
 * implicit fallback. A negative threshold preserves endpoint-Gamma
 * backward-Euler behavior; a nonnegative threshold selects the equilibrium
 * mean-energy projection when dt_alpha*kappa_a_N reaches it. The rate bundle
 * and endpoint current are expected to have crossed their owning validation
 * boundaries. */
ghl_error_codes_t ghl_m1_neutrino_update_endpoint_number_with_policy(
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double dt_alpha,
      const double thermalized_number_threshold,
      const ghl_m1_neutrino_state *restrict state_base,
      const ghl_m1_neutrino_current *restrict endpoint_current,
      double *restrict N_out,
      bool *restrict number_projected);

/* A nonnegative finite binary64 value represented without a range-limited
 * exponent.  The mantissa is normalized to [0.5, 1), except for zero. */
typedef struct {
  double mantissa;
  int exponent;
} ghl_m1_scaled_positive;

static inline bool ghl_m1_scaled_positive_from_double(
      const double value,
      ghl_m1_scaled_positive *restrict scaled) {
  if(scaled == NULL || !isfinite(value) || value < 0.0) {
    return false;
  }
  if(value == 0.0) {
    *scaled = (ghl_m1_scaled_positive){ .mantissa = 0.0, .exponent = 0 };
    return true;
  }
  scaled->mantissa = frexp(value, &scaled->exponent);
  return true;
}

static inline bool ghl_m1_scaled_positive_multiply(
      const ghl_m1_scaled_positive *restrict left,
      const ghl_m1_scaled_positive *restrict right,
      ghl_m1_scaled_positive *restrict product) {
  if(left == NULL || right == NULL || product == NULL) {
    return false;
  }
  if(left->mantissa == 0.0 || right->mantissa == 0.0) {
    *product = (ghl_m1_scaled_positive){ .mantissa = 0.0, .exponent = 0 };
    return true;
  }
  int normalization = 0;
  const double mantissa = frexp(left->mantissa * right->mantissa, &normalization);
  product->mantissa = mantissa;
  product->exponent = left->exponent + right->exponent + normalization;
  return true;
}

static inline bool ghl_m1_scaled_positive_add(
      const ghl_m1_scaled_positive *restrict left,
      const ghl_m1_scaled_positive *restrict right,
      ghl_m1_scaled_positive *restrict sum) {
  if(left == NULL || right == NULL || sum == NULL) {
    return false;
  }
  if(left->mantissa == 0.0) {
    *sum = *right; /* GCOVR_EXCL_LINE -- unreachable */
    return true;   /* GCOVR_EXCL_LINE -- unreachable */
  }
  if(right->mantissa == 0.0) {
    *sum = *left;
    return true;
  }

  const ghl_m1_scaled_positive *larger = left;
  const ghl_m1_scaled_positive *smaller = right;
  if(right->exponent > left->exponent) {
    larger = right;
    smaller = left;
  }
  const double mantissa
        = larger->mantissa
          + scalbn(smaller->mantissa, smaller->exponent - larger->exponent);
  int normalization = 0;
  sum->mantissa = frexp(mantissa, &normalization);
  sum->exponent = larger->exponent + normalization;
  return true;
}

static inline bool ghl_m1_scaled_positive_divide(
      const ghl_m1_scaled_positive *restrict numerator,
      const ghl_m1_scaled_positive *restrict denominator,
      double *restrict quotient) {
  if(numerator == NULL || denominator == NULL || quotient == NULL
     || denominator->mantissa == 0.0) {
    return false;
  }
  if(numerator->mantissa == 0.0) {
    *quotient = 0.0;
    return true;
  }
  *quotient = scalbn(
        numerator->mantissa / denominator->mantissa,
        numerator->exponent - denominator->exponent);
  return isfinite(*quotient) && *quotient >= 0.0;
}

static inline int ghl_m1_scaled_positive_compare(
      const ghl_m1_scaled_positive *restrict left,
      const ghl_m1_scaled_positive *restrict right) {
  if(left->mantissa == 0.0) {
    return right->mantissa == 0.0 ? 0 : -1;
  }
  if(right->mantissa == 0.0) {
    return 1;
  }
  if(left->exponent != right->exponent) {
    return left->exponent < right->exponent ? -1 : 1;
  }
  return left->mantissa < right->mantissa ? -1
                                          : (left->mantissa > right->mantissa ? 1 : 0);
}

static inline bool ghl_m1_scaled_positive_product(
      const double *restrict values,
      const int value_count,
      ghl_m1_scaled_positive *restrict product) {
  if(values == NULL || value_count <= 0 || product == NULL) {
    return false;
  }
  if(!ghl_m1_scaled_positive_from_double(1.0, product)) {
    return false;
  }
  for(int i = 0; i < value_count; ++i) {
    ghl_m1_scaled_positive factor;
    ghl_m1_scaled_positive next_product;
    if(!ghl_m1_scaled_positive_from_double(values[i], &factor)
       || !ghl_m1_scaled_positive_multiply(product, &factor, &next_product)) {
      return false;
    }
    *product = next_product;
  }
  return true;
}

static inline bool ghl_m1_scaled_positive_sqrt(
      const ghl_m1_scaled_positive *restrict value,
      ghl_m1_scaled_positive *restrict root) {
  if(value == NULL || root == NULL) {
    return false;
  }
  if(value->mantissa == 0.0) {
    *root = (ghl_m1_scaled_positive){ .mantissa = 0.0,
                                      .exponent
                                      = 0 }; /* GCOVR_EXCL_LINE -- zero product */
    return true;                             /* GCOVR_EXCL_LINE -- zero product */
  }

  double mantissa = value->mantissa;
  int exponent = value->exponent;
  if(exponent % 2 != 0) {
    mantissa *= 2.0;
    --exponent;
  }
  int normalization = 0;
  root->mantissa = frexp(sqrt(mantissa), &normalization);
  root->exponent = exponent / 2 + normalization;
  return true;
}

/* Evaluate (product of numerator_values) / (product of denominator_values)
 * without materializing an intermediate that overflows or underflows.  Returns
 * false when the true quotient is itself not representable as a finite
 * nonnegative binary64, so a genuinely nonrepresentable result stays an error.
 * Callers retain their direct expression for the normal range and fall back
 * here only when it is not representable, so ordinary results keep their
 * established evaluation order exactly. */
static inline bool ghl_m1_neutrino_scaled_ratio_of_products(
      const double *restrict numerator_values,
      const int numerator_count,
      const double *restrict denominator_values,
      const int denominator_count,
      double *restrict quotient) {
  ghl_m1_scaled_positive numerator;
  ghl_m1_scaled_positive denominator;
  if(quotient == NULL
     || !ghl_m1_scaled_positive_product(numerator_values, numerator_count, &numerator)
     || !ghl_m1_scaled_positive_product(
           denominator_values, denominator_count, &denominator)) {
    return false;
  }
  double candidate = 0.0;
  if(!ghl_m1_scaled_positive_divide(&numerator, &denominator, &candidate)) {
    return false;
  }
  /* A nonzero positive quotient that rounds to zero is below binary64's
   * representable range, not a valid zero result. */
  if(candidate == 0.0 && numerator.mantissa != 0.0) {
    return false;
  }
  *quotient = candidate;
  return true;
}

/* Compare two nonnegative finite products without forming a range-limited
 * intermediate. This is shared by source-policy selection and endpoint-number
 * policy selection; inclusive controls whether equality selects the policy. */
bool ghl_m1_neutrino_scaled_product_meets_threshold(
      const double *restrict left_values,
      int left_count,
      const double *restrict right_values,
      int right_count,
      bool inclusive);

/* Signed charged-current radiation lepton-number increment for validated
 * single-species rates. Ordinary backward Euler uses the un-repaired number
 * change; the optional thermalized projection uses dt_alpha times the
 * charged-current source at its physical endpoint. number_projected is the
 * decision returned by update_endpoint_number_with_policy. Endpoint Gamma_N
 * must be the normalization that selected physical_number_endpoint. */
ghl_error_codes_t ghl_m1_neutrino_charged_current_lepton_delta(
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt_alpha,
      const double number_initial,
      const bool number_projected,
      const double physical_number_endpoint,
      const double physical_number_gamma,
      double *restrict dL_rad_cc);

/* Private final-endpoint validation shared by all local source routes. */
ghl_error_codes_t ghl_m1_neutrino_check_EN_bounds(
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_m1_neutrino_current *restrict current);

/* Closure-supplied localization layer. Production callers use
 * ghl_m1_neutrino_derive_current, which first evaluates the configured
 * primitive-aware closure. Direct verifiers may use this boundary to compare
 * the downstream current algebra from one independently generated tensor. */
ghl_error_codes_t ghl_m1_neutrino_derive_current_from_closure(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_neutrino_current *restrict current);

ghl_error_codes_t ghl_m1_neutrino_physical_number_flux_from_current(
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_neutrino_current *restrict current,
      const ghl_m1_direction_t direction,
      double *restrict physical_number_flux);

/* Solver-observability variant. This private entry point only reports whether
 * the selected closure used a finite nonordinary fallback candidate. */
ghl_error_codes_t ghl_m1_neutrino_compute_EF_interaction_sources_diagnostics(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_neutrino_rates *restrict rates,
      bool *restrict closure_fallback_observed,
      ghl_m1_sources *restrict EF_sources);

/* Internal solve path. The metric and immutable M1 configuration are checked
 * once by the solve entry point; this variant retains trial-state and closure
 * output checks without repeating that boundary validation. */
ghl_error_codes_t ghl_m1_neutrino_compute_EF_interaction_sources_validated(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_neutrino_rates *restrict rates,
      bool *restrict closure_fallback_observed,
      ghl_m1_sources *restrict EF_sources);

/* Private production entry point for the branched source dispatcher. The
 * installed public solver remains the compatibility wrapper below in the
 * implementation and passes a negative threshold, preserving ordinary
 * backward-Euler number integration. */
ghl_error_codes_t ghl_m1_solve_neutrino_implicit_homogeneous_update_with_number_policy(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double n_b_cons,
      const double thermalized_number_threshold,
      const ghl_m1_neutrino_state *restrict state_in,
      ghl_m1_neutrino_state *restrict state_out,
      ghl_m1_neutrino_exchange *restrict exchange,
      ghl_m1_implicit_solve_diagnostics *restrict solve_diagnostics,
      ghl_m1_neutrino_diagnostics *restrict neutrino_diagnostics);

/**
 * Compute the same residual using an explicit densitized substep base U_base.
 *
 * This is used by the fallback substepping driver so substep k solves
 * U_{k+1} - U_k - dt_sub*S(U_{k+1}) = 0 instead of repeatedly subtracting
 * the original full-step input state. The public-like wrapper above passes
 * U_base derived from state_in for single-step tests.
 */
ghl_error_codes_t ghl_m1_neutrino_compute_implicit_residual_with_base(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double U_base[4],
      const double U[4],
      double residual[4]);

ghl_error_codes_t ghl_m1_neutrino_compute_implicit_residual_with_base_diagnostics(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double U_base[4],
      const double U[4],
      bool *restrict closure_fallback_observed,
      double residual[4]);

/* Private callback core. The context's metric/configuration/rates have
 * already crossed their public validation boundary. */
ghl_error_codes_t ghl_m1_neutrino_compute_implicit_residual_validated(
      const ghl_m1_neutrino_implicit_context *restrict context,
      const double dt,
      const double U_base[4],
      const double U[4],
      bool *restrict closure_fallback_observed,
      double residual[4]);

/**
 * Finite-difference Jacobian for the explicit-base residual variant.
 */
ghl_error_codes_t ghl_m1_neutrino_compute_implicit_jacobian_with_base(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double U_base[4],
      const double U[4],
      const double residual_0[4],
      double jacobian[4][4]);

ghl_error_codes_t ghl_m1_neutrino_compute_implicit_jacobian_validated(
      const ghl_m1_neutrino_implicit_context *restrict context,
      const double dt,
      const double U_base[4],
      const double U[4],
      const double residual_0[4],
      double jacobian[4][4]);

/* Run one validated E/F Newton substep. Pair-source code supplies an
 * internally constructed effective rate bundle after validating the original
 * provider rates; this bridge deliberately does not repeat public rate
 * validation. */
ghl_error_codes_t ghl_m1_neutrino_attempt_EF_newton_step(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt_sub,
      const double U_in[4],
      double U_out[4],
      ghl_m1_newton_diagnostics *restrict diagnostics,
      bool *restrict closure_fallback_observed);

/* Populate endpoint mean-energy observability fields after a paired solve. */
void ghl_m1_neutrino_populate_mean_energy_diagnostics(
      const ghl_m1_neutrino_state *restrict state_out,
      const ghl_m1_neutrino_current *restrict current,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_m1_neutrino_rates *restrict rates,
      ghl_m1_neutrino_diagnostics *restrict neutrino_diagnostics);

#endif // GHL_M1_NEUTRINO_IMPLICIT_H
