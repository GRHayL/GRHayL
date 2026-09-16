#include "m1_neutrino_seeded_test_utils.h"
#include "m1_thcm1_fixture_utils.h"
#include "m1_thcm1_stress_energy_fixture.h"
#include "m1_thcm1_source_fixture.h"
#include "../GRHayL/Radiation/ghl_m1_utils.h"

#include <string.h>

/*
 * Seeded local-contract coverage for grey three-species neutrino M1.
 *
 * The local-contract section generates reproducible
 * admissible inputs, runs the public local operators on a baseline and a
 * matched radiation perturbation, and checks finite-value, identity,
 * realizability, species, and transactional contracts. A separate fixture
 * section replays stored independent THC_M1 observations on current GRHayL.
 */

static void require_condition(
      const bool condition,
      const char *restrict label,
      const int case_index,
      const int species) {

  if(!condition)
    ghl_error("seeded neutrino case %d species %d: %s\n",
              case_index, species, label);
}

static void require_finite_value(
      const double value,
      const char *restrict label,
      const int case_index,
      const int species) {

  require_condition(isfinite(value), label, case_index, species);
}

static void require_finite_state(
      const ghl_m1_neutrino_state *restrict state,
      const char *restrict label,
      const int case_index,
      const int species) {

  require_finite_value(state->N, label, case_index, species);
  require_finite_value(state->E, label, case_index, species);
  for(int i = 0; i < 3; ++i)
    require_finite_value(state->F[i], label, case_index, species);
}

static void require_finite_sources(
      const ghl_m1_sources *restrict sources,
      const char *restrict label,
      const int case_index,
      const int species) {

  require_finite_value(sources->S_E, label, case_index, species);
  for(int i = 0; i < 3; ++i)
    require_finite_value(sources->S[i], label, case_index, species);
}

static void check_generated_state(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_neutrino_state *restrict state,
      const int case_index,
      const int species) {

  require_finite_state(state, "generated state is nonfinite", case_index, species);
  require_condition(state->N >= nu_params->N_floor,
                    "generated N is below the configured floor",
                    case_index, species);
  const double flux_factor = m1_neutrino_seeded_flux_factor(metric, state);
  require_finite_value(flux_factor, "generated flux factor is nonfinite",
                       case_index, species);
  require_condition(flux_factor <= 1.0 - m1_params->epsilon_c + 1.0e-12,
                    "generated state is outside the M1 cone",
                    case_index, species);
}

static void check_closure(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure,
      const bool nonzero_flux,
      const int case_index) {

  require_finite_value(closure->chi, "closure chi is nonfinite", case_index, -1);
  require_finite_value(closure->xi, "closure xi is nonfinite", case_index, -1);
  require_finite_value(closure->root_residual,
                       "closure residual is nonfinite", case_index, -1);
  require_condition(closure->chi >= 1.0 / 3.0 && closure->chi <= 1.0,
                    "closure chi is outside its physical interval",
                    case_index, -1);
  require_condition(closure->xi >= 0.0 && closure->xi <= 1.0,
                    "closure xi is outside its physical interval",
                    case_index, -1);
  require_condition(closure->root_residual >= 0.0 &&
                    closure->root_iterations >= 0 &&
                    closure->solve_status != ghl_m1_closure_solve_invalid,
                    "closure metadata is invalid", case_index, -1);

  long double trace = 0.0L;
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      require_finite_value(closure->P[i][j], "closure tensor is nonfinite",
                           case_index, -1);
      trace += (long double)metric->gammaDD[i][j] * closure->P[i][j];
    }
  }
  /* Directed zero-flux cases below check the trace with the production
   * invariant allowance; the generated nonzero-flux cases check it here. */
  if(nonzero_flux)
    require_condition(m1_nearly_equal((double)trace, state->E,
                                      5.0e-10, 5.0e-12),
                      "nonzero-flux closure trace identity failed",
                      case_index, -1);

  (void)m1_params;
}

static void check_closure_decomposition_diagnostic(
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure,
      const ghl_m1_closure_decomposition_diagnostic *restrict diagnostic,
      const int case_index) {

  require_finite_value(diagnostic->chi_minerbo_xi,
                       "decomposition diagnostic chi is nonfinite",
                       case_index, -1);
  require_finite_value(diagnostic->xi_HaHa_over_J2,
                       "decomposition diagnostic xi is nonfinite",
                       case_index, -1);
  require_finite_value(diagnostic->dthin_scalar,
                       "decomposition diagnostic thin weight is nonfinite",
                       case_index, -1);
  require_finite_value(diagnostic->dthick_scalar,
                       "decomposition diagnostic thick weight is nonfinite",
                       case_index, -1);
  require_condition(m1_nearly_equal(diagnostic->chi_minerbo_xi,
                                    closure->chi, 2.0e-12, 2.0e-14),
                    "decomposition diagnostic chi disagrees with closure",
                    case_index, -1);
  require_condition(m1_nearly_equal(diagnostic->xi_HaHa_over_J2,
                                    closure->xi * closure->xi,
                                    2.0e-11, 2.0e-13),
                    "decomposition diagnostic xi disagrees with closure",
                    case_index, -1);
  require_condition(m1_nearly_equal(
                        diagnostic->dthin_scalar,
                        0.5 * (3.0 * closure->chi - 1.0),
                        2.0e-12, 2.0e-14) &&
                    m1_nearly_equal(
                        diagnostic->dthick_scalar,
                        1.5 * (1.0 - closure->chi),
                        2.0e-12, 2.0e-14),
                    "decomposition diagnostic weights are inconsistent",
                    case_index, -1);

  double thick_trace = 0.0;
  for(int i = 0; i < 3; ++i) {
    require_finite_value(diagnostic->Pthick_minus_Pthin_dd[i],
                         "decomposition diagonal difference is nonfinite",
                         case_index, -1);
    require_condition(m1_nearly_equal(
                          diagnostic->Pthick_minus_Pthin_dd[i],
                          diagnostic->Pthick_dd[i][i] -
                              diagnostic->Pthin_dd[i][i],
                          2.0e-12, 2.0e-14),
                      "decomposition diagonal difference identity failed",
                      case_index, -1);
    for(int j = 0; j < 3; ++j) {
      require_finite_value(diagnostic->Pthin_dd[i][j],
                           "thin decomposition tensor is nonfinite",
                           case_index, -1);
      require_finite_value(diagnostic->Pthick_dd[i][j],
                           "thick decomposition tensor is nonfinite",
                           case_index, -1);
      require_condition(m1_nearly_equal(diagnostic->Pthin_dd[i][j],
                                        diagnostic->Pthin_dd[j][i],
                                        2.0e-12, 2.0e-14) &&
                        m1_nearly_equal(diagnostic->Pthick_dd[i][j],
                                        diagnostic->Pthick_dd[j][i],
                                        2.0e-12, 2.0e-14),
                        "decomposition tensors are not symmetric",
                        case_index, -1);
      thick_trace += metric->gammaUU[i][j] * diagnostic->Pthick_dd[i][j];
    }
  }
  require_finite_value(diagnostic->Pth_dd_3_3_UU,
                       "thick decomposition trace is nonfinite",
                       case_index, -1);
  require_finite_value(diagnostic->Pth_dd_0_0_DD,
                       "thick decomposition time component is nonfinite",
                       case_index, -1);
  require_condition(m1_nearly_equal(diagnostic->Pth_dd_3_3_UU,
                                    thick_trace / 3.0,
                                    2.0e-11, 2.0e-13) &&
                    m1_nearly_equal(diagnostic->Pth_dd_0_0_DD,
                                    diagnostic->Pthick_dd[0][0],
                                    2.0e-12, 2.0e-14),
                    "decomposition trace identities failed",
                    case_index, -1);

  /* The normal production closure is the weighted lower-index decomposition.
   * Exceptional admissibility fallbacks intentionally do not satisfy this
   * identity, but their diagnostic fields above remain well-defined. */
  if(closure->four_point_compatibility) {
    for(int i = 0; i < 3; ++i)
      for(int j = 0; j < 3; ++j) {
        double reconstructed = 0.0;
        for(int k = 0; k < 3; ++k)
          for(int l = 0; l < 3; ++l)
            reconstructed += metric->gammaUU[i][k] * metric->gammaUU[j][l]
                * (diagnostic->dthin_scalar * diagnostic->Pthin_dd[k][l]
                   + diagnostic->dthick_scalar * diagnostic->Pthick_dd[k][l]);
        require_condition(m1_nearly_equal(reconstructed, closure->P[i][j],
                                          3.0e-10, 3.0e-12),
                          "closure/decomposition reconstruction failed",
                          case_index, -1);
      }
  }

  (void)state;
}

static void check_comoving(
      const ghl_m1_comoving *restrict comoving,
      const int case_index) {

  require_finite_value(comoving->J, "comoving J is nonfinite", case_index, -1);
  require_condition(comoving->J >= 0.0, "comoving J is negative", case_index, -1);
  require_finite_value(comoving->Hn, "comoving Hn is nonfinite", case_index, -1);
  for(int i = 0; i < 3; ++i) {
    require_finite_value(comoving->HU[i], "comoving HU is nonfinite",
                         case_index, -1);
    require_finite_value(comoving->HD[i], "comoving HD is nonfinite",
                         case_index, -1);
  }
}

static void check_diagnostics(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_diagnostics *restrict diagnostics,
      const int case_index) {

  require_finite_value(diagnostics->closure_xi,
                       "diagnostic xi is nonfinite", case_index, -1);
  require_finite_value(diagnostics->closure_root_residual,
                       "diagnostic residual is nonfinite", case_index, -1);
  require_finite_value(diagnostics->r, "diagnostic r is nonfinite", case_index, -1);
  require_finite_value(diagnostics->chi_eddington,
                       "diagnostic chi is nonfinite", case_index, -1);
  require_condition(diagnostics->r >= 0.0 &&
                    diagnostics->r <= m1_params->one_minus_epsilon_c_sq + 1.0e-12,
                    "diagnostic r is outside its physical interval",
                    case_index, -1);
}

static void check_stress_energy(
      const ghl_stress_energy *restrict stress_energy,
      const int case_index) {

  for(int i = 0; i < 4; ++i) {
    for(int j = 0; j < 4; ++j) {
      require_finite_value(stress_energy->T4[i][j],
                           "stress-energy output is nonfinite", case_index, -1);
      require_condition(stress_energy->T4[i][j] == stress_energy->T4[j][i],
                        "stress-energy output is not symmetric",
                        case_index, -1);
    }
  }
}

static double compute_gamma_n(
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_comoving *restrict comoving) {

  double velocity[3];
  for(int i = 0; i < 3; ++i)
    velocity[i] = (prims->vU[i] + metric->betaU[i]) / metric->lapse;
  const double velocity_norm =
      m1_neutrino_seeded_metric_norm(metric->gammaDD, velocity);
  const double W = 1.0 / sqrt((1.0 - velocity_norm) * (1.0 + velocity_norm));
  return W - comoving->Hn / comoving->J;
}

static void check_transactional_contracts(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const m1_neutrino_seeded_case *restrict test_case,
      const ghl_m1_closure *restrict closure,
      const int case_index) {

  const int species = ghl_m1_neutrino_nue;

  /* Invalid number-flux input must not publish output. */
  ghl_m1_neutrino_state invalid_state = test_case->state;
  invalid_state.E = NAN;
  double number_flux[3] = {31.0, 32.0, 33.0};
  double transport_velocity[3] = {41.0, 42.0, 43.0};
  const ghl_error_codes_t flux_error = ghl_m1_compute_neutrino_number_flux(
      m1_params, nu_params, &test_case->metric, &test_case->prims,
      &invalid_state, number_flux, transport_velocity);
  require_condition(flux_error == ghl_error_m1_invalid_state,
                    "invalid number-flux state was accepted", case_index, species);
  for(int i = 0; i < 3; ++i) {
    require_condition(number_flux[i] == 31.0 + i,
                      "number flux changed after rejected input",
                      case_index, species);
    require_condition(transport_velocity[i] == 41.0 + i,
                      "transport velocity changed after rejected input",
                      case_index, species);
  }

  /* Invalid frozen rates must leave both source outputs untouched. */
  ghl_m1_neutrino_rates invalid_rates = test_case->rates[species];
  invalid_rates.kappa_a_E = -1.0;
  invalid_rates.eta_E = 1.0;
  ghl_m1_sources sources = {.S_E = 51.0, .S = {52.0, 53.0, 54.0}};
  const ghl_m1_sources sources_before = sources;
  double N_source = 55.0;
  const ghl_error_codes_t source_error =
      ghl_m1_compute_neutrino_interaction_sources_from_closure(
          m1_params, nu_params, &test_case->metric, &test_case->prims,
          &test_case->state, closure, &invalid_rates, &sources, &N_source);
  require_condition(source_error == ghl_error_m1_microphysics_failure,
                    "invalid rates were accepted by source operator",
                    case_index, species);
  require_condition(sources.S_E == sources_before.S_E &&
                    sources.S[0] == sources_before.S[0] &&
                    sources.S[1] == sources_before.S[1] &&
                    sources.S[2] == sources_before.S[2] &&
                    N_source == 55.0,
                    "source output changed after rejected rates",
                    case_index, species);

  /* Invalid supplied closure must not publish a stress-energy tensor. */
  ghl_m1_closure invalid_closure = *closure;
  invalid_closure.P[0][0] = NAN;
  ghl_stress_energy stress_energy;
  for(int i = 0; i < 4; ++i)
    for(int j = 0; j < 4; ++j)
      stress_energy.T4[i][j] = 61.0 + 4.0 * i + j;
  const ghl_stress_energy stress_before = stress_energy;
  const ghl_error_codes_t stress_error = ghl_m1_compute_neutrino_stress_energy(
      m1_params, &test_case->metric, &test_case->state, &invalid_closure,
      &stress_energy);
  require_condition(stress_error != ghl_success,
                    "invalid supplied closure was accepted", case_index, -1);
  for(int i = 0; i < 4; ++i)
    for(int j = 0; j < 4; ++j)
      require_condition(stress_energy.T4[i][j] == stress_before.T4[i][j],
                        "stress-energy changed after rejected closure",
                        case_index, -1);

  /* Nonfinite repair input must be rejected before the state is published. */
  ghl_m1_neutrino_state invalid_repair_state = test_case->state;
  invalid_repair_state.E = NAN;
  const double N_before = invalid_repair_state.N;
  const double F_before[3] = {invalid_repair_state.F[0],
                              invalid_repair_state.F[1],
                              invalid_repair_state.F[2]};
  const ghl_error_codes_t repair_error = ghl_m1_repair_neutrino_state(
      m1_params, nu_params, &test_case->metric, &invalid_repair_state, NULL);
  require_condition(repair_error == ghl_error_m1_invalid_state,
                    "nonfinite repair input was accepted", case_index, -1);
  require_condition(isnan(invalid_repair_state.E) &&
                    invalid_repair_state.N == N_before &&
                    invalid_repair_state.F[0] == F_before[0] &&
                    invalid_repair_state.F[1] == F_before[1] &&
                    invalid_repair_state.F[2] == F_before[2],
                    "repair state changed after rejected input", case_index, -1);
}

enum { pointwise_input_count = 68, pointwise_output_count = 38 };

static void evaluate_pointwise_fixture(
      const double input[pointwise_input_count],
      double output[pointwise_output_count],
      bool *admissibility_fallback) {
  const double alpha = input[4];
  const double *beta = input + 5;
  const double *gamma = input + 8;
  ghl_metric_quantities metric;
  ghl_initialize_metric(alpha, beta[0], beta[1], beta[2],
      gamma[0], gamma[1], gamma[2], gamma[4], gamma[5], gamma[8], &metric);
  ghl_primitive_quantities prims = {0};
  double v2 = 0.0;
  for(int i = 0; i < 3; ++i) {
    prims.vU[i] = input[17+i];
    for(int j = 0; j < 3; ++j)
      v2 += gamma[3*i+j] * (input[17+i]+beta[i]) *
            (input[17+j]+beta[j]) / (alpha*alpha);
  }
  if(!(alpha > 0.0) || !isfinite(v2) || v2 < 0.0 || v2 >= 1.0)
    ghl_error("pointwise fixture has invalid lapse or velocity\n");
  prims.u0 = 1.0 / (alpha * sqrt(1.0-v2));
  ghl_extrinsic_curvature curv = {0};
  ghl_metric_quantities deriv[3] = {0};
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      curv.K[i][j] = input[20+3*i+j];
  for(int d = 0; d < 3; ++d) {
    deriv[d].lapse = input[29+d];
    for(int i = 0; i < 3; ++i) {
      deriv[d].betaU[i] = input[32+3*d+i];
      for(int j = 0; j < 3; ++j)
        deriv[d].gammaDD[i][j] = input[41+9*d+3*i+j];
    }
  }
  /* Frozen controls of the retained pointwise producer. */
  ghl_m1_parameters params;
  ghl_error_codes_t error = ghl_m1_initialize(
      1e-8, 1e-12, 1.0, 1e-6, 1e-10, 100, 1e-10, &params);
  if(error != ghl_success) ghl_error("pointwise fixture initializer failed\n");
  const ghl_m1_rad_state rad = {.E=input[0], .F={input[1],input[2],input[3]}};
  ghl_m1_closure closure;
  ghl_m1_comoving moments;
  ghl_stress_energy stress;
  ghl_m1_sources geometry;
  error = ghl_m1_compute_closure_with_primitives(&params, &metric, &prims, &rad, &closure);
  if(error != ghl_success) ghl_error("pointwise fixture closure failed: %d\n", error);
  *admissibility_fallback =
      closure.solve_status == ghl_m1_closure_solve_endpoint_fallback &&
      !closure.four_point_compatibility;
  error = ghl_m1_compute_comoving_moments(&params, &metric, &prims, &rad, &closure, &moments);
  if(error != ghl_success) ghl_error("pointwise fixture moments failed: %d\n", error);
  error = ghl_m1_compute_stress_energy(&params, &metric, &rad, &closure, &stress);
  if(error != ghl_success) ghl_error("pointwise fixture stress failed: %d\n", error);
  error = ghl_m1_compute_geometry_sources(
      &params, &metric, deriv, deriv+1, deriv+2, &curv, &rad, &closure, &geometry);
  if(error != ghl_success) ghl_error("pointwise fixture geometry failed: %d\n", error);
  size_t next = 0;
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j) output[next++] = closure.P[i][j];
  output[next++] = closure.chi;
  output[next++] = closure.xi;
  output[next++] = moments.J;
  for(int i = 0; i < 3; ++i) output[next++] = moments.HD[i];
  for(int i = 0; i < 4; ++i)
    for(int j = 0; j < 4; ++j) output[next++] = stress.T4[i][j];
  output[next++] = geometry.S_E;
  for(int i = 0; i < 3; ++i) output[next++] = geometry.S[i];
  for(int d = 0; d < 3; ++d) {
    double sm, sp;
    error = ghl_m1_compute_wavespeeds(&metric, (ghl_m1_direction_t)d, &sm, &sp);
    if(error != ghl_success) ghl_error("pointwise fixture wavespeed failed: %d\n", error);
    output[next++] = fmax(-sm, sp);
  }
}

static bool zero_flux_reference_trace_difference(
      const double *input, const double *reference) {
  if(input[1] != 0.0 || input[2] != 0.0 || input[3] != 0.0) return false;
  long double trace = 0.0L;
  for(int i = 0; i < 9; ++i) trace += (long double)input[8+i] * reference[i];
  /* The published-tensor contract uses this bound in ghl_m1_utils.h. */
  return fabsl(trace-input[0]) >
      128.0L * DBL_EPSILON * fmaxl(fabsl(trace), fabsl(input[0]));
}

static void check_zero_flux_admissibility(
      const double *input, const double *actual, const bool fallback,
      const char *case_id) {
  long double trace = 0.0L;
  for(int i = 0; i < 9; ++i) trace += (long double)input[8+i] * actual[i];
  if(!fallback || !m1_thcm1_fixture_finite_vector(actual, pointwise_output_count) ||
     fabsl(trace-input[0]) >
       128.0L * DBL_EPSILON * fmaxl(fabsl(trace), fabsl(input[0])) ||
     fabs(actual[9]-1.0/3.0) > 128.0 * DBL_EPSILON)
    ghl_error("fixture %s did not preserve the zero-flux admissibility contract\n", case_id);
}

static void check_pointwise_fixtures(const char *fixture_dir) {
  static const char *const anchors[] = {
    "radiation-anchor-energy-floor-a01", "radiation-anchor-number-floor-a01",
    "radiation_N-anchor-number-floor-a01", "radiation-anchor-realizability-cone-a01",
    "radiation-anchor-closure-transition-a01", "matter-anchor-limiter-threshold-a01",
    "matter-anchor-eos-lower-a01", "matter-anchor-eos-upper-a01",
    "rates-anchor-wavespeed-cap-a01", "rates-anchor-stiff-a01",
    "rates-anchor-frozen-rate-endpoint-a01", "velocity-anchor-high-lorentz-a01",
    "metric-anchor-offdiagonal-spd-a01", "geometry-anchor-source-balance-a01"
  };
  static const char *const seeded_families[] = {
    "radiation", "radiation_N", "matter", "rates", "velocity", "metric", "geometry"
  };
  const char filename[] = "/pointwise_closure_moments.m1";
  const size_t path_size = strlen(fixture_dir) + sizeof(filename);
  char *path = malloc(path_size);
  if(path == NULL) {
    ghl_error("fixture path allocation failed\n");
    return;
  }
  snprintf(path, path_size, "%s%s", fixture_dir, filename);
  m1_thcm1_fixture_collection fixtures = {0};
  char error[M1_THCM1_FIXTURE_TEXT_MAX];
  if(!m1_thcm1_fixture_load(path, "pointwise_closure_moments",
        pointwise_input_count, pointwise_output_count, &fixtures, error, sizeof(error)))
    ghl_error("%s: %s\n", path, error);
  free(path);
  if(strcmp(fixtures.policy, "pointwise_a1_a2_v1") != 0)
    ghl_error("unrecognized pointwise fixture comparison policy\n");
  /* The retained effective-packet corpus contains 56 pairs. This independent
   * count prevents a shortened file/header from silently reducing coverage. */
  if(fixtures.record_count != 56)
    ghl_error("pointwise fixture corpus is incomplete\n");
  ghl_m1_reset_closure_counters();
  size_t sensitive = 0, controls = 0, admissibility_cases = 0;
  for(size_t i = 0; i < fixtures.record_count; ++i) {
    const m1_thcm1_fixture_record *record = &fixtures.records[i];
    /* Packet ordering and seed/amplitude axes are those of the retained
     * corpus; do not let unique-ID duplicates replace its physical cases. */
    if(i < sizeof(anchors)/sizeof(anchors[0])) {
      const char prefix[] = "rngpkt-v2-rd-";
      if(strncmp(record->case_id, prefix, sizeof(prefix)-1) != 0 ||
         strcmp(record->case_id + sizeof(prefix)-1, anchors[i]) != 0)
        ghl_error("pointwise directed packet inventory differs at %zu\n", i);
    } else {
      const size_t sample = i - sizeof(anchors)/sizeof(anchors[0]);
      char expected[sizeof("rngpkt-v2-rd-radiation_N-seed-v1-0000-a00")];
      snprintf(expected, sizeof(expected), "rngpkt-v2-rd-%s-seed-v1-%04zu-a%02zu",
               seeded_families[sample/6], (sample%6)/3, sample%3);
      if(strcmp(record->case_id, expected) != 0)
        ghl_error("pointwise seeded packet inventory differs at %zu\n", i);
    }
    const char *case_family = record->case_id + strlen("rngpkt-v2-rd-");
    const size_t family_length = strlen(record->family);
    if(strncmp(case_family, record->family, family_length) != 0 ||
       case_family[family_length] != '-' ||
       strcmp(record->origin, i < sizeof(anchors)/sizeof(anchors[0])
              ? "directed-anchor" : "seeded") != 0)
      ghl_error("pointwise packet family/origin differs at %zu\n", i);
    char expected_pair_id[M1_THCM1_FIXTURE_TEXT_MAX];
    snprintf(expected_pair_id, sizeof(expected_pair_id), "%s-pair",
             record->case_id);
    if(strcmp(record->pair_id, expected_pair_id) != 0)
      ghl_error("pointwise packet pair identity differs at %zu\n", i);
    double baseline[pointwise_output_count];
    double normalization[pointwise_output_count];
    const double scale = fmax(fabs(record->baseline_input[0]),
                               fabs(record->perturbed_input[0]));
    for(size_t j = 0; j < pointwise_output_count; ++j) {
      normalization[j] = (j == 9 || j == 10 || j >= 35) ? 1.0 : scale;
      if(normalization[j] != record->normalization[j])
        ghl_error("fixture %s normalization differs from its input\n", record->case_id);
    }
    bool baseline_fallback;
    evaluate_pointwise_fixture(record->baseline_input, baseline, &baseline_fallback);
    const bool baseline_difference = zero_flux_reference_trace_difference(
        record->baseline_input, record->baseline_output);
    const bool perturbed_difference = zero_flux_reference_trace_difference(
        record->perturbed_input, record->perturbed_output);
    /* Explicitly reviewed policy exceptions, not a runtime skip-on-mismatch.
     * A new discrepancy or a changed reference classification must fail. */
    const bool known_admissibility_case =
        strcmp(record->case_id, "rngpkt-v2-rd-radiation-anchor-energy-floor-a01") == 0 ||
        strcmp(record->case_id, "rngpkt-v2-rd-metric-anchor-offdiagonal-spd-a01") == 0;
    if(baseline_difference || perturbed_difference || known_admissibility_case) {
      if(!known_admissibility_case || !baseline_difference || !perturbed_difference)
        ghl_error("fixture %s has an unreviewed zero-flux policy classification\n", record->case_id);
      /* These two records are local fallback-policy assertions, not standard
       * cross-code replays. Preserve their existing two-state check while
       * keeping the ordinary replay baseline-only. */
      double perturbed[pointwise_output_count];
      bool perturbed_fallback;
      evaluate_pointwise_fixture(record->perturbed_input, perturbed,
                                 &perturbed_fallback);
      check_zero_flux_admissibility(record->baseline_input, baseline,
                                    baseline_fallback, record->case_id);
      check_zero_flux_admissibility(record->perturbed_input, perturbed,
                                    perturbed_fallback, record->case_id);
      ++admissibility_cases;
      continue;
    }
    m1_thcm1_fixture_comparison_report report;
    if(!m1_thcm1_fixture_compare_envelope(record, fixtures.policy, normalization,
                                          baseline, &report, error, sizeof(error)))
      ghl_error("%s\n", error);
    if(record->sensitivity_count) ++sensitive;
    else ++controls;
  }
  /* Exact inventory of the retained directed corpus, not sampling quotas. */
  if(sensitive != 29 || controls != 25 || admissibility_cases != 2)
    ghl_error("pointwise fixture corpus changed its sensitivity/control/policy inventory\n");
  ghl_m1_closure_counters counters;
  ghl_m1_get_closure_counters(NULL); /* Documented no-op. */
  ghl_m1_get_closure_counters(&counters);
  if(counters.ordinary_convergence + counters.endpoint_fallback +
       counters.iteration_exhaustion != fixtures.record_count +
           admissibility_cases ||
     counters.endpoint_fallback < 2 * admissibility_cases ||
     counters.invalid_state != 0)
    ghl_error("closure counters do not describe the completed fixture evaluations\n");
  ghl_m1_reset_closure_counters();
  ghl_m1_get_closure_counters(&counters);
  if(counters.ordinary_convergence || counters.endpoint_fallback ||
     counters.iteration_exhaustion || counters.invalid_state ||
     counters.downstream_repair || counters.residual_rejection)
    ghl_error("closure counter reset did not clear completed evaluations\n");
  ghl_info("THC_M1 pointwise fixtures: %zu agreement pairs (%zu sensitivity, %zu input-invariant controls); %zu distinct zero-flux policy pairs passed local admissibility checks\n",
           sensitive+controls, sensitive, controls, admissibility_cases);
  m1_thcm1_fixture_free(&fixtures);
}

static void check_shared_closure_boundary_contracts(void) {
  ghl_m1_parameters params;
  if(ghl_m1_initialize(0.5, 1.0e-12, 1.0, 1.0e-6, 1.0e-10,
                        100, 1.0e-10, &params) != ghl_success)
    ghl_error("shared closure boundary parameter initialization failed\n");

  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  ghl_primitive_quantities prims = {0};
  prims.u0 = 1.0;
  const ghl_m1_rad_state static_zero_flux = {.E = 1.0, .F = {0.0, 0.0, 0.0}};
  ghl_m1_closure closure;
  ghl_error_codes_t error = ghl_m1_compute_closure_with_primitives(
      &params, &metric, &prims, &static_zero_flux, &closure);
  if(error != ghl_success || closure.solve_status != ghl_m1_closure_solve_converged ||
     closure.root_iterations != 0 || !closure.four_point_compatibility)
    ghl_error("static zero-flux closure did not take its endpoint contract\n");
  const ghl_m1_neutrino_state static_zero_flux_state = {
      .N = 1.0, .E = static_zero_flux.E,
      .F = {static_zero_flux.F[0], static_zero_flux.F[1], static_zero_flux.F[2]}};
  check_closure(&params, &metric, &static_zero_flux_state, &closure, false, -1);
  long double static_trace = 0.0L;
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      static_trace += (long double)metric.gammaDD[i][j] * closure.P[i][j];
  if(fabsl(static_trace - static_zero_flux.E) >
     128.0L * DBL_EPSILON * fmaxl(fabsl(static_trace), static_zero_flux.E))
    ghl_error("static zero-flux closure failed its pressure trace invariant\n");

  /* With zero Eulerian flux but a moving fluid, the primary covariant tensor
   * can fail its observable trace check.  The documented fallback publishes
   * an Eulerian Minerbo tensor and exposes the exceptional status. */
  prims.vU[0] = 0.4;
  prims.u0 = 1.0 / sqrt(1.0 - prims.vU[0] * prims.vU[0]);
  error = ghl_m1_compute_closure_with_primitives(
      &params, &metric, &prims, &static_zero_flux, &closure);
  if(error != ghl_success ||
     closure.solve_status != ghl_m1_closure_solve_endpoint_fallback ||
     closure.four_point_compatibility || closure.root_iterations != 0 ||
     !isfinite(closure.xi) || closure.xi < 0.0 || closure.xi > 1.0)
    ghl_error("moving zero-flux closure did not publish its admissibility fallback\n");
  const ghl_m1_neutrino_state moving_zero_flux_state = {
      .N = 1.0, .E = static_zero_flux.E,
      .F = {static_zero_flux.F[0], static_zero_flux.F[1], static_zero_flux.F[2]}};
  check_closure(&params, &metric, &moving_zero_flux_state, &closure, false, -1);
  long double moving_trace = 0.0L;
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      moving_trace += (long double)metric.gammaDD[i][j] * closure.P[i][j];
  if(fabsl(moving_trace - static_zero_flux.E) >
     128.0L * DBL_EPSILON * fmaxl(fabsl(moving_trace), static_zero_flux.E))
    ghl_error("moving zero-flux fallback failed its pressure trace invariant\n");

  /* A positive finite residual tolerance is part of the public control
   * contract.  Setting it to the smallest representable positive value makes
   * a finite nonzero solve residual observable without changing production
   * tolerances or adding a test-only seam. */
  ghl_m1_parameters residual_params = params;
  const double true_min = nextafter(0.0, 1.0);
  if(ghl_m1_set_closure_residual_tolerance(true_min, &residual_params) != ghl_success)
    ghl_error("could not set the public residual tolerance boundary\n");
  const ghl_m1_rad_state nonzero_flux = {.E = 1.0, .F = {0.37, 0.0, 0.0}};
  error = ghl_m1_compute_closure_with_primitives(
      &residual_params, &metric, &(ghl_primitive_quantities){0},
      &nonzero_flux, &closure);
  if(error != ghl_error_m1_closure_residual_too_large)
    ghl_error("closure residual acceptance boundary was not enforced\n");

  /* A one-iteration public solver configuration must publish the explicit
   * iteration-exhausted status when residual acceptance is otherwise open. */
  ghl_m1_parameters exhausted_params = params;
  if(ghl_m1_set_closure_solver_controls(true_min, 1, &exhausted_params) != ghl_success ||
     ghl_m1_set_closure_residual_tolerance(1.0, &exhausted_params) != ghl_success)
    ghl_error("could not set the public closure iteration controls\n");
  error = ghl_m1_compute_closure_with_primitives(
      &exhausted_params, &metric, &(ghl_primitive_quantities){0},
      &nonzero_flux, &closure);
  if(error != ghl_success ||
     closure.solve_status != ghl_m1_closure_solve_iteration_exhausted ||
     closure.root_iterations != 1 || !isfinite(closure.root_residual))
    ghl_error("closure iteration exhaustion was not published\n");

  const ghl_m1_closure exhausted_sentinel = closure;
  if(ghl_m1_set_closure_residual_tolerance(true_min, &exhausted_params) != ghl_success)
    ghl_error("could not tighten exhausted-solve residual acceptance\n");
  ghl_m1_reset_closure_counters();
  error = ghl_m1_compute_closure_with_primitives(
      &exhausted_params, &metric, &(ghl_primitive_quantities){0},
      &nonzero_flux, &closure);
  ghl_m1_closure_counters rejected_counters;
  ghl_m1_get_closure_counters(&rejected_counters);
  if(error != ghl_error_m1_closure_residual_too_large ||
     memcmp(&closure, &exhausted_sentinel, sizeof(closure)) != 0 ||
     rejected_counters.iteration_exhaustion != 1 ||
     rejected_counters.residual_rejection != 1 ||
     rejected_counters.invalid_state != 0)
    ghl_error("rejected exhausted closure lost transactional/counter contract\n");

  /* This state selects the range-scaled residual path (E > sqrt(DBL_MAX))
   * while remaining finite and inside the realizability cone. */
  const double large_energy = nextafter(sqrt(DBL_MAX), INFINITY);
  const ghl_m1_rad_state scaled_state = {
      .E = large_energy, .F = {0.25 * large_energy, 0.0, 0.0}};
  const ghl_m1_rad_state unit_scale_state = {
      .E = 1.0, .F = {0.25, 0.0, 0.0}};
  ghl_m1_closure unit_scale_closure;
  error = ghl_m1_compute_closure_with_primitives(
      &params, &metric, &(ghl_primitive_quantities){0},
      &unit_scale_state, &unit_scale_closure);
  if(error != ghl_success)
    ghl_error("unit-scale closure oracle failed\n");
  error = ghl_m1_compute_closure_with_primitives(
      &params, &metric, &(ghl_primitive_quantities){0},
      &scaled_state, &closure);
  if(error != ghl_success || !isfinite(closure.chi) || !isfinite(closure.xi) ||
     !isfinite(closure.root_residual))
    ghl_error("range-scaled closure state was not published\n");
  /* The Brent implementation accepts an endpoint when its bracket midpoint
   * is no larger than tau (ghl_m1_closure.c:721-726). Each returned endpoint
   * is therefore at most 2*tau from the exact root; comparing two independent
   * homogeneous evaluations gives a 4*tau xi bound. For this static,
   * axis-aligned state, the published Minerbo polynomial is
   * chi(xi) (ghl_m1_closure.c:116-119), with
   * dchi/dxi=(12*xi-6*xi^2+24*xi^3)/15 <= 2 on [0,1]. Thus the propagated
   * chi bound is 8*tau. The corresponding normalized pressure entries are
   * Pxx/E=chi and Pyy/E=Pzz/E=(1-chi)/2, so 8*tau is also a conservative
   * component bound. Add the existing unit-scale invariant slack, rather than
   * inventing a decimal allowance, for the final floating-point operations.
   */
  const double root_endpoint_bound = 4.0 * params.closure_root_tolerance;
  const double arithmetic_bound = m1_invariant_slack(1.0);
  const double homogeneous_bound = 2.0 * root_endpoint_bound
                                 + arithmetic_bound;
  if(!m1_nearly_equal(closure.chi, unit_scale_closure.chi,
                      homogeneous_bound, homogeneous_bound) ||
     !m1_nearly_equal(closure.xi, unit_scale_closure.xi,
                      homogeneous_bound, homogeneous_bound))
    ghl_error("range-scaled closure changed its homogeneous chi/xi result\n");
  long double scaled_trace = 0.0L;
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      scaled_trace += (long double)metric.gammaDD[i][j] * closure.P[i][j];
      if(!m1_nearly_equal(closure.P[i][j] / large_energy,
                          unit_scale_closure.P[i][j],
                          homogeneous_bound, homogeneous_bound))
        ghl_error("range-scaled closure changed its normalized pressure tensor\n");
    }
  }
  if(fabsl(scaled_trace - scaled_state.E) >
     128.0L * DBL_EPSILON * fmaxl(fabsl(scaled_trace), scaled_state.E))
    ghl_error("range-scaled closure failed its pressure trace invariant\n");
  const ghl_m1_neutrino_state scaled_neutrino_state = {
      .N = 1.0, .E = scaled_state.E,
      .F = {scaled_state.F[0], scaled_state.F[1], scaled_state.F[2]}};
  check_closure(&params, &metric, &scaled_neutrino_state,
                &closure, true, -1);
}

static void check_low_lapse_shift_closure(void) {
  ghl_m1_parameters params;
  require_condition(ghl_m1_initialize(
      1.0e-12, 1.0e-12, 1.0e-12, 1.0e-6, 1.0e-12, 80, 1.0e-10,
      &params) == ghl_success,
      "low-lapse closure parameters", -1, -1);

  const double lapse_values[] = {1.0, 1.0e-2, 1.0e-3, 1.0e-4};
  const ghl_m1_rad_state rad_state = {
      .E = 1.0, .F = {0.3, 0.2, 0.1}};
  ghl_m1_closure reference = {0};
  for(size_t case_index = 0;
      case_index < sizeof(lapse_values) / sizeof(lapse_values[0]);
      ++case_index) {
    ghl_metric_quantities metric;
    m1_setup_flat_metric(&metric);
    metric.lapse = lapse_values[case_index];
    metric.lapseinv = 1.0 / metric.lapse;
    metric.lapseinv2 = metric.lapseinv * metric.lapseinv;
    metric.betaU[0] = 0.5;
    metric.betaU[1] = 0.2;

    ghl_primitive_quantities prims = {0};
    prims.u0 = metric.lapseinv;
    prims.vU[0] = -metric.betaU[0];
    prims.vU[1] = -metric.betaU[1];
    ghl_m1_closure closure;
    const ghl_error_codes_t error = ghl_m1_compute_closure_with_primitives(
        &params, &metric, &prims, &rad_state, &closure);
    require_condition(error == ghl_success,
                      "low-lapse shifted closure failed", (int)case_index, -1);
    const ghl_m1_neutrino_state neutrino_state = {
        .N = 1.0, .E = rad_state.E,
        .F = {rad_state.F[0], rad_state.F[1], rad_state.F[2]}};
    check_closure(&params, &metric, &neutrino_state, &closure, true,
                  (int)case_index);
    if(case_index == 0) {
      reference = closure;
    } else {
      require_condition(m1_nearly_equal(closure.chi, reference.chi,
                                        1.0e-10, 1.0e-12) &&
                        m1_nearly_equal(closure.xi, reference.xi,
                                        1.0e-10, 1.0e-12),
                        "low-lapse closure changed with lapse or shift",
                        (int)case_index, -1);
      for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
          require_condition(m1_nearly_equal(closure.P[i][j], reference.P[i][j],
                                             1.0e-10, 1.0e-12),
                            "low-lapse pressure changed with lapse or shift",
                            (int)case_index, -1);
    }
  }
}

static void check_near_zero_flux_closure(void) {
  ghl_m1_parameters params;
  require_condition(ghl_m1_initialize(1.0e-12, 1.0e-20, 1.0e-12, 1.0e-6,
      1.0e-12, 100, 1.0e-10, &params) == ghl_success,
      "near-zero flux parameters", -1, -1);
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);

  /* Reproduced cases: E/F^2 overflows at 1e-155 and 1e-160, while F^2
   * underflows at 1e-170. The safe nonzero reference has the same direction
   * and negligible flux, without either intermediate range failure. */
  const double fluxes[] = {1.0e-155, 1.0e-160, 1.0e-170};
  const double velocities[] = {0.0, 0.3};
  const double tolerance = m1_invariant_slack(1.0);
  for(size_t v = 0; v < sizeof(velocities) / sizeof(velocities[0]); ++v) {
    ghl_primitive_quantities prims = {0};
    prims.vU[0] = velocities[v];
    prims.u0 = 1.0 / sqrt(1.0 - prims.vU[0] * prims.vU[0]);
    const ghl_m1_rad_state reference_state = {
        .E = 1.0, .F = {1.0e-150, 0.0, 0.0}};
    ghl_m1_closure reference;
    require_condition(ghl_m1_compute_closure_with_primitives(
        &params, &metric, &prims, &reference_state, &reference) == ghl_success,
        "near-zero flux reference closure", -1, -1);

    for(size_t f = 0; f < sizeof(fluxes) / sizeof(fluxes[0]); ++f) {
      ghl_m1_rad_state state = {.E = 1.0, .F = {fluxes[f], 0.0, 0.0}};
      require_condition(ghl_m1_realizability_repair(
          &params, &metric, &state) == ghl_success,
          "near-zero flux realizability", (int)f, (int)v);
      ghl_m1_closure closure;
      require_condition(ghl_m1_compute_closure_with_primitives(
          &params, &metric, &prims, &state, &closure) == ghl_success,
          "near-zero flux closure", (int)f, (int)v);
      require_condition(closure.four_point_compatibility &&
          closure.solve_status == ghl_m1_closure_solve_converged,
          "nonzero flux incorrectly used zero-flux fallback", (int)f, (int)v);

      double trace = 0.0;
      for(int i = 0; i < 3; ++i) {
        trace += closure.P[i][i];
        /* The axis-aligned pressure is diagonal, so nonnegative diagonal
         * entries also establish positive semidefiniteness. */
        require_condition(closure.P[i][i] >= 0.0,
            "near-zero flux negative pressure", (int)f, (int)v);
        for(int j = 0; j < 3; ++j) {
          const double expected = velocities[v] == 0.0
              ? (i == j ? state.E / 3.0 : 0.0) : reference.P[i][j];
          require_condition(isfinite(closure.P[i][j]) &&
              fabs(closure.P[i][j] - expected) <= tolerance &&
              closure.P[i][j] == closure.P[j][i] &&
              (i == j || closure.P[i][j] == 0.0),
              "near-zero flux pressure changed its directional limit", (int)f, (int)v);
        }
      }
      require_condition(fabs(trace - state.E) <=
          128.0 * DBL_EPSILON * fmax(fabs(trace), state.E),
          "near-zero flux pressure trace", (int)f, (int)v);
    }
  }
}

static void check_nonzero_flux_admissibility_fallback(void) {
  ghl_m1_parameters params;
  require_condition(ghl_m1_initialize(1.0e-8, 1.0e-30, 1.0, 1.0e-6,
      1.0e-10, 100, 1.0e-10, &params) == ghl_success,
      "nonzero fallback initialization", -1, -1);
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  ghl_primitive_quantities prims = {.vU = {0.0, 0.8, 0.0}, .u0 = 5.0 / 3.0};
  /* Transverse motion makes this primary covariant pressure non-PSD. The
   * documented fallback uses Eulerian f=3/10, whose Minerbo polynomial is
   * exactly chi=27673/75000. The three energy scales exercise ordinary,
   * low-energy, and normalized comoving admissibility arithmetic for the same
   * physical tensor. */
  const double energies[] = {1.0, 1.0e-20, nextafter(sqrt(DBL_MAX), INFINITY)};
  const double expected_chi = 27673.0 / 75000.0;
  const double tolerance = m1_invariant_slack(1.0);
  for(size_t scale = 0; scale < sizeof(energies) / sizeof(energies[0]); ++scale) {
    const ghl_m1_rad_state state = {
        .E = energies[scale], .F = {0.3 * energies[scale], 0.0, 0.0}};
    ghl_m1_closure closure;
    require_condition(ghl_m1_compute_closure_with_primitives(
        &params, &metric, &prims, &state, &closure) == ghl_success,
        "nonzero flux fallback", -1, -1);
    require_condition(closure.solve_status == ghl_m1_closure_solve_endpoint_fallback &&
        !closure.four_point_compatibility && closure.root_iterations == 0,
        "nonzero flux did not select the admissibility fallback", -1, -1);
    require_condition(fabs(closure.chi - expected_chi) <= tolerance,
        "fallback changed the analytic Minerbo factor", -1, -1);
    ghl_m1_closure_failure_stage_t failure_stage;
    int validation_reason;
    ghl_m1_get_last_closure_failure_stage(&failure_stage);
    ghl_m1_get_last_closure_validation_reason(&validation_reason);
    require_condition(failure_stage == ghl_m1_closure_failure_none &&
        validation_reason == GHL_M1_CLOSURE_VALIDATION_PSD,
        "admissibility fallback lost its PSD diagnostic", -1, -1);
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        const double expected = i != j ? 0.0 :
            (i == 0 ? expected_chi : (1.0 - expected_chi) / 2.0);
        require_condition(isfinite(closure.P[i][j]) &&
            fabs(closure.P[i][j] / state.E - expected) <= tolerance,
            "fallback changed the analytic pressure tensor", -1, -1);
      }
    }
  }
}

static void check_closure_arithmetic_boundaries(void) {
  ghl_m1_parameters params;
  require_condition(ghl_m1_initialize(0.1, 1.0e-12, 1.0, 1.0e-6,
      1.0e-10, 100, 1.0e-10, &params) == ghl_success,
      "closure arithmetic initializer", -1, -1);

  const ghl_m1_closure sentinel = {
      .chi = 7.0 / 11.0, .xi = 5.0 / 13.0,
      .root_residual = 3.0 / 17.0, .root_iterations = 19,
      .solve_status = ghl_m1_closure_solve_iteration_exhausted,
      .four_point_compatibility = true,
      .P = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}}};
  ghl_m1_closure closure;
  ghl_m1_closure_failure_stage_t failure_stage;
  int validation_reason;
  ghl_m1_closure_counters counters;

  /* A finite shift and coordinate velocity can cancel to an Eulerian-rest
   * state, while beta^2 overflows the intermediate four-metric. The public
   * closure must reject the workspace before publishing partial tensors. */
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  metric.betaU[0] = DBL_MAX;
  ghl_primitive_quantities prims = {
      .vU = {-DBL_MAX, 0.0, 0.0}, .u0 = 1.0};
  const ghl_m1_rad_state finite_state = {
      .E = 1.0, .F = {0.0, 0.0, 0.0}};
  closure = sentinel;
  ghl_m1_reset_closure_counters();
  require_condition(ghl_m1_compute_closure_with_primitives(
      &params, &metric, &prims, &finite_state, &closure) ==
      ghl_error_m1_invalid_state,
      "overflowing finite shift was accepted", -1, -1);
  require_condition(memcmp(&closure, &sentinel, sizeof(closure)) == 0,
      "workspace overflow published a partial closure", -1, -1);
  ghl_m1_get_last_closure_failure_stage(&failure_stage);
  ghl_m1_get_last_closure_validation_reason(&validation_reason);
  ghl_m1_get_closure_counters(&counters);
  require_condition(failure_stage == ghl_m1_closure_failure_workspace &&
      validation_reason == 0 && counters.invalid_state == 1 &&
      counters.ordinary_convergence == 0 && counters.endpoint_fallback == 0 &&
      counters.iteration_exhaustion == 0 && counters.residual_rejection == 0,
      "workspace overflow diagnostic accounting failed", -1, -1);

  /* The scaled norm result itself can exceed double range even though every
   * metric/state operand is finite. Realizability must reject that result and
   * leave the caller's output object untouched. */
  metric = (ghl_metric_quantities){0};
  metric.lapse = metric.lapseinv = metric.lapseinv2 = 1.0;
  metric.detgamma = 1.0e-300;
  metric.sqrt_detgamma = 1.0e-150;
  for(int i = 0; i < 3; ++i) {
    metric.gammaDD[i][i] = 1.0e-100;
    metric.gammaUU[i][i] = 1.0e100;
  }
  prims = (ghl_primitive_quantities){.u0 = 1.0};
  const ghl_m1_rad_state overflowing_flux = {
      .E = 1.0, .F = {DBL_MAX, 0.0, 0.0}};
  closure = sentinel;
  ghl_m1_reset_closure_counters();
  require_condition(ghl_m1_compute_closure_with_primitives(
      &params, &metric, &prims, &overflowing_flux, &closure) ==
      ghl_error_m1_invalid_state,
      "overflowing realizability ratio was accepted", -1, -1);
  require_condition(memcmp(&closure, &sentinel, sizeof(closure)) == 0,
      "realizability overflow published a partial closure", -1, -1);
  ghl_m1_get_last_closure_failure_stage(&failure_stage);
  ghl_m1_get_last_closure_validation_reason(&validation_reason);
  ghl_m1_get_closure_counters(&counters);
  require_condition(failure_stage == ghl_m1_closure_failure_none &&
      validation_reason == 0 && counters.invalid_state == 1 &&
      counters.ordinary_convergence == 0 && counters.endpoint_fallback == 0 &&
      counters.iteration_exhaustion == 0 && counters.residual_rejection == 0,
      "realizability overflow diagnostic accounting failed", -1, -1);

  /* Keep the state realizable, but make the Eulerian fallback pressure
   * unrepresentable through the inverse spatial metric. This reaches the
   * fallback's public failure path after the primary PSD rejection. */
  prims = (ghl_primitive_quantities){.vU = {0.0, 0.8, 0.0},
      .u0 = 5.0 / 3.0};
  const ghl_m1_rad_state overflowing_fallback = {
      .E = DBL_MAX, .F = {0.0, 0.0, 0.0}};
  closure = sentinel;
  ghl_m1_reset_closure_counters();
  require_condition(ghl_m1_compute_closure_with_primitives(
      &params, &metric, &prims, &overflowing_fallback, &closure) ==
      ghl_error_m1_invalid_state,
      "overflowing fallback pressure was accepted", -1, -1);
  require_condition(memcmp(&closure, &sentinel, sizeof(closure)) == 0,
      "fallback overflow published a partial closure", -1, -1);
  ghl_m1_get_last_closure_failure_stage(&failure_stage);
  ghl_m1_get_last_closure_validation_reason(&validation_reason);
  ghl_m1_get_closure_counters(&counters);
  require_condition(failure_stage == ghl_m1_closure_failure_tensor_validation &&
      validation_reason == GHL_M1_CLOSURE_VALIDATION_PSD &&
      counters.invalid_state == 1 && counters.endpoint_fallback == 0 &&
      counters.iteration_exhaustion == 0 && counters.residual_rejection == 0,
      "fallback overflow diagnostic accounting failed", -1, -1);

  /* Finite, realizable inputs can still exceed the representable range of
   * individual closure operations. Exercise representative arithmetic
   * failures through the public API and require transactional output and
   * useful diagnostics. */
  const struct {
    double energy, lapse, shift, velocity;
  } arithmetic_cases[] = {
      {DBL_MAX, 1.0, 2.0, 0.0},
      {DBL_MAX, 2.0, 0.0, 0.0},
      {sqrt(DBL_MAX), 1.0e100, 0.0, 0.0},
      {sqrt(DBL_MAX), 1.0, 0.0, 0.8},
      /* Near-light Eulerian motion makes the comoving-energy arithmetic
       * overflow while all input fields remain finite and realizable. */
      {1.0e100, 1.0, 0.0, 0.999999999},
      /* The same motion at the scaled-energy boundary exercises the
       * long-double comoving-flux norm rejection. */
      {sqrt(DBL_MAX), 1.0, 0.0, 0.999999999}};
  for(size_t i = 0; i < sizeof(arithmetic_cases)/sizeof(arithmetic_cases[0]); ++i) {
    ghl_initialize_metric(arithmetic_cases[i].lapse,
        arithmetic_cases[i].shift, 0.0, 0.0,
        1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
    prims = (ghl_primitive_quantities){
        .vU = {arithmetic_cases[i].lapse * arithmetic_cases[i].velocity
               - arithmetic_cases[i].shift, 0.0, 0.0},
        .u0 = 1.0 / (arithmetic_cases[i].lapse
                     * sqrt(1.0 - SQR(arithmetic_cases[i].velocity)))};
    const ghl_m1_rad_state state = {
        .E = arithmetic_cases[i].energy,
        .F = {0.3 * arithmetic_cases[i].energy, 0.0, 0.0}};
    closure = sentinel;
    ghl_m1_reset_closure_counters();
    require_condition(ghl_m1_compute_closure_with_primitives(
        &params, &metric, &prims, &state, &closure) == ghl_error_m1_invalid_state,
        "closure arithmetic failure was accepted", (int)i, -1);
    require_condition(memcmp(&closure, &sentinel, sizeof(closure)) == 0,
        "closure arithmetic failure published partial output", (int)i, -1);
    ghl_m1_get_last_closure_failure_stage(&failure_stage);
    ghl_m1_get_last_closure_validation_reason(&validation_reason);
    ghl_m1_get_closure_counters(&counters);
    /* FP contraction can change which intermediate first loses range.
     * Require a recorded arithmetic failure without prescribing the
     * compiler's evaluation order. */
    require_condition(failure_stage >= ghl_m1_closure_failure_workspace &&
        failure_stage <= ghl_m1_closure_failure_residual &&
        validation_reason == 0 && counters.invalid_state == 1 &&
        counters.ordinary_convergence == 0 && counters.endpoint_fallback == 0 &&
        counters.iteration_exhaustion == 0 && counters.residual_rejection == 0,
        "closure arithmetic failure diagnostics", (int)i, -1);
  }

  /* A near-light-speed fluid in a small-lapse coordinate system reaches
   * the unbracketed endpoint policy, whose residual must still be accepted
   * before any pressure is published. Depending on contraction, arithmetic
   * can instead become invalid before that endpoint is reached. */
  const double speed = nextafter(1.0, 0.0);
  ghl_initialize_metric(1.0e-100, 0.0, 0.0, 0.0,
      1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  prims = (ghl_primitive_quantities){.vU = {1.0e-100 * speed, 0.0, 0.0},
      .u0 = 1.0 / sqrt(1.0 - speed * speed) / 1.0e-100};
  const ghl_m1_rad_state endpoint_state = {.E = 1.0, .F = {0.3, 0.0, 0.0}};
  closure = sentinel;
  ghl_m1_reset_closure_counters();
  const ghl_error_codes_t endpoint_error = ghl_m1_compute_closure_with_primitives(
      &params, &metric, &prims, &endpoint_state, &closure);
  require_condition(endpoint_error == ghl_error_m1_closure_residual_too_large ||
      endpoint_error == ghl_error_m1_invalid_state,
      "unbracketed endpoint bypassed residual gate", -1, -1);
  ghl_m1_get_closure_counters(&counters);
  require_condition(memcmp(&closure, &sentinel, sizeof(closure)) == 0 &&
      counters.endpoint_fallback ==
          (endpoint_error == ghl_error_m1_closure_residual_too_large) &&
      counters.residual_rejection ==
          (endpoint_error == ghl_error_m1_closure_residual_too_large) &&
      counters.invalid_state == (endpoint_error == ghl_error_m1_invalid_state) &&
      counters.ordinary_convergence == 0 && counters.iteration_exhaustion == 0,
      "rejected endpoint output/counter contract", -1, -1);

  /* In a sheared SPD metric, raising the primary tensor loses enough
   * precision to violate its pressure trace even though the tensor is
   * positive semidefinite. Finite flux must not take the zero-flux trace
   * fallback, and the invalid candidate must remain unpublished. */
  const double shear = 1.0 - 0x1p-10;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0,
      1.0, shear, 0.0, 1.0, 0.0, 1.0, &metric);
  prims = (ghl_primitive_quantities){.u0 = 1.0};
  const ghl_m1_rad_state trace_state = {
      .E = 1.0, .F = {0.3 * sqrt(1.0 - shear * shear), 0.0, 0.0}};
  closure = sentinel;
  ghl_m1_reset_closure_counters();
  require_condition(ghl_m1_compute_closure_with_primitives(
      &params, &metric, &prims, &trace_state, &closure) == ghl_error_m1_invalid_state,
      "inaccurate finite-flux pressure trace was accepted", -1, -1);
  ghl_m1_get_last_closure_failure_stage(&failure_stage);
  ghl_m1_get_last_closure_validation_reason(&validation_reason);
  ghl_m1_get_closure_counters(&counters);
  require_condition(memcmp(&closure, &sentinel, sizeof(closure)) == 0 &&
      failure_stage == ghl_m1_closure_failure_tensor_validation &&
      validation_reason == GHL_M1_CLOSURE_VALIDATION_TRACE &&
      counters.invalid_state == 1 && counters.endpoint_fallback == 0,
      "finite-flux trace rejection output/diagnostic contract", -1, -1);
}

static void check_supplied_closure_psd_range(void) {
  ghl_m1_parameters params;
  require_condition(ghl_m1_initialize(0.1, 1.0e-310, 1.0, 1.0e-6,
      1.0e-10, 100, 1.0e-10, &params) == ghl_success,
      "PSD range initializer", -1, -1);
  const double scales[] = {1.0e-100, 1.0};
  const double energies[] = {1.0e-300, 1.0};
  for(size_t i = 0; i < sizeof(scales)/sizeof(scales[0]); ++i) {
    ghl_metric_quantities metric;
    ghl_initialize_metric(1.0, 0.0, 0.0, 0.0,
        scales[i], 0.0, 0.0, scales[i], 0.0, scales[i], &metric);
    const ghl_m1_rad_state state = {.E = energies[i]};
    ghl_m1_closure closure = {0};
    for(int j = 0; j < 3; ++j)
      closure.P[j][j] = energies[i] / scales[i] / 3.0;
    /* First, lowering a trace-consistent tensor underflows to zero. Second,
     * a finite, symmetric supplied tensor overflows during symmetrization
     * in the PSD basis. Neither may publish stress energy. */
    if(i == 1) closure.P[0][1] = closure.P[1][0] = DBL_MAX;
    const ghl_stress_energy sentinel = {.T4 = {{17.0}}};
    ghl_stress_energy stress = sentinel;
    require_condition(ghl_m1_compute_stress_energy(
        &params, &metric, &state, &closure, &stress) == ghl_error_m1_invalid_state &&
        memcmp(&stress, &sentinel, sizeof(stress)) == 0,
        "PSD range failure published stress energy", (int)i, -1);
    int reason;
    ghl_m1_get_last_closure_validation_reason(&reason);
    require_condition(reason == (i == 0 ? GHL_M1_CLOSURE_VALIDATION_PSD :
                                         GHL_M1_CLOSURE_VALIDATION_NONFINITE),
        "PSD range failure diagnostic", (int)i, -1);
  }

  /* The private tensor validator also owns its Cholesky rejection. Test
   * that local boundary directly: public operators reject this indefinite
   * metric earlier, before reaching pressure validation. */
  const ghl_metric_quantities indefinite = {
      .gammaDD = {{1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 0.0, 1.0}}};
  const ghl_m1_closure pressure = {
      .P = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
  require_condition(ghl_m1_validate_closure_tensor_psd(&indefinite, &pressure)
      == ghl_error_m1_invalid_state,
      "PSD helper accepted an indefinite metric", -1, -1);
  int reason;
  ghl_m1_get_last_closure_validation_reason(&reason);
  require_condition(reason == GHL_M1_CLOSURE_VALIDATION_PSD,
      "PSD helper Cholesky rejection diagnostic", -1, -1);

  /* A nonzero antisymmetric tensor has a zero symmetric part. The helper
   * must reject that degenerate eigenproblem instead of normalizing by
   * zero; the public validator separately rejects its asymmetry. */
  ghl_metric_quantities flat;
  m1_setup_flat_metric(&flat);
  const ghl_m1_closure antisymmetric = {
      .P = {{0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
  require_condition(ghl_m1_validate_closure_tensor_psd(&flat, &antisymmetric)
      == ghl_error_m1_invalid_state,
      "PSD helper accepted a zero symmetric tensor", -1, -1);
  ghl_m1_get_last_closure_validation_reason(&reason);
  require_condition(reason == GHL_M1_CLOSURE_VALIDATION_PSD,
      "PSD helper zero eigen-scale rejection diagnostic", -1, -1);
}

static void check_inline_wrapper_boundaries(void) {
  ghl_m1_parameters params;
  require_condition(ghl_m1_initialize(0.5, 1.0e-12, 1.0, 1.0e-6,
      1.0e-10, 100, 1.0e-10, &params) == ghl_success,
      "inline boundary parameters", -1, -1);
  double floored_energy = -1.0;
  bool floor_applied = false;
  require_condition(ghl_m1_apply_energy_floor(&params, 0.0, &floored_energy,
      &floor_applied) == ghl_success && floor_applied &&
      floored_energy == params.E_floor,
      "energy floor flag on repair", -1, -1);
  require_condition(ghl_m1_apply_energy_floor(&params, 1.0, &floored_energy,
      &floor_applied) == ghl_success && !floor_applied && floored_energy == 1.0,
      "energy floor flag on unchanged input", -1, -1);
  ghl_m1_neutrino_parameters nu;
  m1_neutrino_seeded_default_parameters(&nu);
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  ghl_metric_quantities derivative = {0};
  ghl_extrinsic_curvature curvature = {0};
  ghl_primitive_quantities prims = {.u0 = 1.0};
  ghl_m1_neutrino_state state = {.N = 1.0, .E = 1.0};
  ghl_m1_closure closure = {.chi = 1.0/3.0,
      .P = {{1.0/3.0, 0.0, 0.0}, {0.0, 1.0/3.0, 0.0}, {0.0, 0.0, 1.0/3.0}}};
  ghl_m1_comoving comoving;
  ghl_stress_energy stress;
  ghl_m1_sources sources;
  ghl_m1_diagnostics diagnostics;
  const ghl_m1_neutrino_state *volatile no_state = NULL;
#define NULL_WRAPPER(call) require_condition((call) == ghl_error_m1_null_pointer, \
    #call, -1, -1)
  NULL_WRAPPER(ghl_m1_compute_neutrino_closure(&params, &metric, &prims, no_state, &closure));
  NULL_WRAPPER(ghl_m1_compute_neutrino_comoving_moments(
      &params, &metric, &prims, no_state, &closure, &comoving));
  NULL_WRAPPER(ghl_m1_compute_neutrino_stress_energy(
      &params, &metric, no_state, &closure, &stress));
  NULL_WRAPPER(ghl_m1_compute_neutrino_geometry_sources(&params, &metric,
      &derivative, &derivative, &derivative, &curvature, no_state, &closure, &sources));
  NULL_WRAPPER(ghl_m1_compute_neutrino_diagnostics(&params, &metric, no_state, &closure, &diagnostics));
  ghl_m1_neutrino_rates rates = {.species = ghl_m1_neutrino_nux,
      .mean_energy = 1.0, .n_eq = 1.0, .J_eq = 1.0};
  double E = -3.0, N = -4.0, F[3] = {-5.0, -6.0, -7.0};
#define RHS(s, p, r, ep, fp, np, include) ghl_m1_compute_neutrino_explicit_rhs_sources( \
    &params, &nu, &metric, &derivative, &derivative, &derivative, &curvature, \
    (p), (s), &closure, (include), (r), (ep), (fp), (np))
  double *volatile no_output = NULL;
  const ghl_primitive_quantities *volatile no_prims = NULL;
  const ghl_m1_neutrino_rates *volatile no_rates = NULL;
  NULL_WRAPPER(RHS(no_state, &prims, &rates, &E, F, &N, true));
  NULL_WRAPPER(RHS(&state, &prims, &rates, no_output, F, &N, true));
  NULL_WRAPPER(RHS(&state, &prims, &rates, &E, no_output, &N, true));
  NULL_WRAPPER(RHS(&state, no_prims, &rates, &E, F, &N, true));
  NULL_WRAPPER(RHS(&state, &prims, no_rates, &E, F, &N, true));
#undef NULL_WRAPPER
  require_condition(RHS(&state, NULL, NULL, &E, F, NULL, false) == ghl_success &&
      E == 0.0 && F[0] == 0.0 && F[1] == 0.0 && F[2] == 0.0,
      "geometry-only RHS without optional number output", -1, -1);
  metric.lapse = -1.0;
  require_condition(RHS(&state, &prims, &rates, &E, F, &N, true) ==
      ghl_error_m1_invalid_metric, "RHS geometry failure", -1, -1);
  metric.lapse = 1.0;
  rates.mean_energy = 0.0;
  require_condition(RHS(&state, &prims, &rates, &E, F, &N, true) ==
      ghl_error_m1_microphysics_failure, "RHS interaction failure", -1, -1);
  /* Isolate each publication component. Undensitized interaction sources
   * remain finite; lapse weighting alone overflows E, N, or momentum. */
  for(int component = 0; component < 3; ++component) {
    metric.lapse = component == 2 ? 100.0 : 2.0;
    metric.lapseinv = 1.0 / metric.lapse;
    metric.lapseinv2 = metric.lapseinv * metric.lapseinv;
    prims.u0 = metric.lapseinv;
    state.F[0] = component == 2 ? 0.1 : 0.0;
    rates = (ghl_m1_neutrino_rates){.species = ghl_m1_neutrino_nux,
        .mean_energy = 1.0, .n_eq = DBL_MAX, .J_eq = DBL_MAX};
    if(component == 0) {
      rates.kappa_a_E = rates.kappa_tr = 1.0;
      rates.eta_E = DBL_MAX;
    } else if(component == 1) {
      rates.kappa_a_N = 1.0;
      rates.eta_N = DBL_MAX;
    } else {
      rates.kappa_s = rates.kappa_tr = DBL_MAX;
    }
    E = -3.0; N = -4.0;
    F[0] = -5.0; F[1] = -6.0; F[2] = -7.0;
    require_condition(ghl_m1_validate_neutrino_rates(&rates, NULL) == ghl_success,
        "RHS overflow rates must be valid", -1, -1);
    double undensitized_N;
    require_condition(ghl_m1_compute_neutrino_interaction_sources(&params, &nu,
        &metric, &prims, &state, &rates, &sources, &undensitized_N) == ghl_success,
        "RHS overflow must pass undensitized interaction evaluation", -1, -1);
    require_condition(RHS(&state, &prims, &rates, &E, F, &N, true) ==
        ghl_error_m1_invalid_state, "RHS overflow publication", -1, -1);
    require_condition(E == -3.0 && N == -4.0 && F[0] == -5.0 &&
        F[1] == -6.0 && F[2] == -7.0, "RHS error changed output", -1, -1);
  }
#undef RHS
}

int main(int argc, char **argv) {
  const char *fixture_dir = "Unit_Tests/data/m1_thcm1";
  if(argc == 3 && strcmp(argv[1], "--fixture-dir") == 0) fixture_dir = argv[2];
  else if(argc != 1) ghl_error("Usage: %s [--fixture-dir PATH]\n", argv[0]);
  check_shared_closure_boundary_contracts();
  check_low_lapse_shift_closure();
  check_near_zero_flux_closure();
  check_nonzero_flux_admissibility_fallback();
  check_closure_arithmetic_boundaries();
  check_supplied_closure_psd_range();
  check_inline_wrapper_boundaries();
  check_pointwise_fixtures(fixture_dir);
  const char stress_energy_filename[] = "/stress_energy.m1";
  const size_t stress_energy_path_size =
      strlen(fixture_dir) + sizeof(stress_energy_filename);
  char *stress_energy_path = malloc(stress_energy_path_size);
  if(stress_energy_path == NULL) {
    ghl_error("stress-energy fixture path allocation failed\n");
    return EXIT_FAILURE;
  }
  snprintf(stress_energy_path, stress_energy_path_size, "%s%s", fixture_dir,
           stress_energy_filename);
  char stress_energy_error[M1_THCM1_FIXTURE_TEXT_MAX];
  size_t stress_energy_records = 0;
  if(!m1_thcm1_run_stress_energy_fixtures(
         stress_energy_path, &stress_energy_records, stress_energy_error,
         sizeof(stress_energy_error)))
    ghl_error("%s\n", stress_energy_error);
  free(stress_energy_path);
  ghl_info("THC_M1 stress-energy fixtures: %zu agreement pairs\n",
           stress_energy_records);
  const char source_filename[] = "/m1_thcm1_instantaneous_sources.m1";
  const size_t source_path_size = strlen(fixture_dir) + sizeof(source_filename);
  char *source_path = malloc(source_path_size);
  if(source_path == NULL) {
    ghl_error("source fixture path allocation failed\n");
    return EXIT_FAILURE;
  }
  snprintf(source_path, source_path_size, "%s%s", fixture_dir, source_filename);
  char source_error[M1_THCM1_FIXTURE_TEXT_MAX];
  size_t source_records = 0;
  if(!m1_thcm1_run_instantaneous_source_fixtures(source_path, &source_records,
        source_error, sizeof(source_error)))
    ghl_error("%s\n", source_error);
  free(source_path);
  ghl_info("THC_M1 instantaneous source fixtures: %zu agreement pairs\n", source_records);

  ghl_m1_parameters m1_params = {0};
  ghl_error_codes_t error = ghl_m1_initialize(
      1.0e-10, 1.0e-12, 1.0e-8, 1.0e-6, 1.0e-12,
      20, 1.0e-10, &m1_params);
  if(error != ghl_success)
    ghl_error("ghl_m1_initialize failed with code %d\n", (int)error);

  ghl_m1_neutrino_parameters nu_params;
  m1_neutrino_seeded_default_parameters(&nu_params);

  m1_neutrino_seeded_rng rng = {
      .state = M1_NEUTRINO_SEEDED_PRNG_SEED};
  unsigned long long cases_run = 0;
  unsigned long long species_runs = 0;

  for(int case_index = 0;
      case_index < M1_NEUTRINO_SEEDED_CASE_COUNT;
      ++case_index) {
    m1_neutrino_seeded_case test_case;
    m1_neutrino_seeded_make_case(&rng, case_index, &test_case);
    const bool nonzero_flux =
        m1_neutrino_seeded_flux_factor(&test_case.metric, &test_case.state) > 0.0;

    for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
      const ghl_error_codes_t rates_error =
          ghl_m1_validate_neutrino_rates(&test_case.rates[species], NULL);
      require_condition(rates_error == ghl_success,
                        "generated rates failed validation", case_index, species);
      check_generated_state(&m1_params, &nu_params, &test_case.metric,
                            &test_case.state, case_index, species);
      check_generated_state(&m1_params, &nu_params, &test_case.metric,
                            &test_case.perturbed_state, case_index, species);
      ++species_runs;
    }

    /* One shared E/F closure is intentionally reused for every species. */
    ghl_m1_closure closure;
    error = ghl_m1_compute_neutrino_closure(
        &m1_params, &test_case.metric, &test_case.prims,
        &test_case.state, &closure);
    require_condition(error == ghl_success, "baseline closure failed",
                      case_index, -1);
    check_closure(&m1_params, &test_case.metric, &test_case.state,
                  &closure, nonzero_flux, case_index);

    ghl_m1_closure_decomposition_diagnostic decomposition_diagnostic = {0};
    const ghl_m1_rad_state diagnostic_rad_state =
        ghl_m1_neutrino_project_rad_state(&test_case.state);
    error = ghl_m1_compute_closure_decomposition_diagnostic(
        &m1_params, &test_case.metric, &test_case.prims,
        &diagnostic_rad_state, &decomposition_diagnostic);
    require_condition(error == ghl_success,
                      "closure decomposition diagnostic failed",
                      case_index, -1);
    check_closure_decomposition_diagnostic(
        &test_case.metric, &test_case.state, &closure,
        &decomposition_diagnostic, case_index);

    ghl_m1_closure_failure_stage_t failure_stage = ghl_m1_closure_failure_none;
    int validation_reason = -1;
    ghl_m1_get_last_closure_failure_stage(&failure_stage);
    ghl_m1_get_last_closure_validation_reason(&validation_reason);
    const bool closure_fallback =
        closure.solve_status == ghl_m1_closure_solve_endpoint_fallback ||
        closure.solve_status == ghl_m1_closure_solve_iteration_exhausted;
    require_condition(closure_fallback ||
                      (failure_stage == ghl_m1_closure_failure_none &&
                       validation_reason == 0),
                      "successful converged closure left failure diagnostics set",
                      case_index, -1);
    ghl_m1_get_last_closure_failure_stage(NULL);
    ghl_m1_get_last_closure_validation_reason(NULL);

    ghl_m1_closure perturbed_closure;
    error = ghl_m1_compute_neutrino_closure(
        &m1_params, &test_case.metric, &test_case.prims,
        &test_case.perturbed_state, &perturbed_closure);
    require_condition(error == ghl_success, "perturbed closure failed",
                      case_index, -1);
    check_closure(&m1_params, &test_case.metric,
                  &test_case.perturbed_state, &perturbed_closure,
                  nonzero_flux, case_index);

    ghl_m1_comoving comoving;
    error = ghl_m1_compute_neutrino_comoving_moments(
        &m1_params, &test_case.metric, &test_case.prims,
        &test_case.state, &closure, &comoving);
    require_condition(error == ghl_success, "baseline comoving moments failed",
                      case_index, -1);
    check_comoving(&comoving, case_index);

    ghl_m1_comoving perturbed_comoving;
    error = ghl_m1_compute_neutrino_comoving_moments(
        &m1_params, &test_case.metric, &test_case.prims,
        &test_case.perturbed_state, &perturbed_closure, &perturbed_comoving);
    require_condition(error == ghl_success,
                      "perturbed comoving moments failed", case_index, -1);
    check_comoving(&perturbed_comoving, case_index);

    ghl_m1_diagnostics diagnostics;
    error = ghl_m1_compute_neutrino_diagnostics(
        &m1_params, &test_case.metric, &test_case.state, &closure, &diagnostics);
    require_condition(error == ghl_success, "diagnostics failed", case_index, -1);
    check_diagnostics(&m1_params, &diagnostics, case_index);

    ghl_stress_energy stress_energy;
    error = ghl_m1_compute_neutrino_stress_energy(
        &m1_params, &test_case.metric, &test_case.state, &closure,
        &stress_energy);
    require_condition(error == ghl_success, "baseline stress-energy failed",
                      case_index, -1);
    check_stress_energy(&stress_energy, case_index);

    ghl_stress_energy perturbed_stress_energy;
    error = ghl_m1_compute_neutrino_stress_energy(
        &m1_params, &test_case.metric, &test_case.perturbed_state,
        &perturbed_closure, &perturbed_stress_energy);
    require_condition(error == ghl_success,
                      "perturbed stress-energy failed", case_index, -1);
    check_stress_energy(&perturbed_stress_energy, case_index);
    require_condition(perturbed_stress_energy.T4[0][0] != stress_energy.T4[0][0],
                      "radiation perturbation produced no stress response",
                      case_index, -1);

    /* Geometry and explicit RHS are exercised with finite generated metric
     * derivatives. The no-interaction RHS must reproduce geometry exactly. */
    ghl_m1_sources geometry_sources;
    error = ghl_m1_compute_neutrino_geometry_sources(
        &m1_params, &test_case.metric, &test_case.metric_derivs[0],
        &test_case.metric_derivs[1], &test_case.metric_derivs[2],
        &test_case.curv, &test_case.state, &closure, &geometry_sources);
    require_condition(error == ghl_success, "geometry sources failed",
                      case_index, -1);
    require_finite_sources(&geometry_sources, "geometry source is nonfinite",
                           case_index, -1);

    double rhs_E = NAN;
    double rhs_F[3] = {NAN, NAN, NAN};
    double rhs_N = NAN;
    error = ghl_m1_compute_neutrino_explicit_rhs_sources(
        &m1_params, &nu_params, &test_case.metric,
        &test_case.metric_derivs[0], &test_case.metric_derivs[1],
        &test_case.metric_derivs[2], &test_case.curv, &test_case.prims,
        &test_case.state, &closure, false, NULL, &rhs_E, rhs_F, &rhs_N);
    require_condition(error == ghl_success, "geometry-only explicit RHS failed",
                      case_index, -1);
    require_condition(m1_nearly_equal(rhs_E, geometry_sources.S_E, 1.0e-13, 1.0e-14),
                      "geometry-only RHS energy mismatch", case_index, -1);
    for(int i = 0; i < 3; ++i)
      require_condition(m1_nearly_equal(rhs_F[i], geometry_sources.S[i],
                                        1.0e-13, 1.0e-14),
                        "geometry-only RHS momentum mismatch", case_index, -1);

    double baseline_number_flux[3] = {0.0, 0.0, 0.0};
    double baseline_transport_velocity[3] = {0.0, 0.0, 0.0};
    double perturbed_number_flux[3] = {0.0, 0.0, 0.0};
    double perturbed_transport_velocity[3] = {0.0, 0.0, 0.0};
    const int transport_reference_species = ghl_m1_neutrino_nue;

    /* The same transport result must be independent of the species-specific
     * frozen rate bundle, and the supplied-closure path must agree with the
     * fresh-closure path. */
    for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
      double fresh_flux[3], fresh_velocity[3];
      double null_rate_flux[3], null_rate_velocity[3];
      double supplied_flux[3], supplied_velocity[3];
      error = ghl_m1_compute_neutrino_number_flux(
          &m1_params, &nu_params, &test_case.metric, &test_case.prims,
          &test_case.state, fresh_flux, fresh_velocity);
      require_condition(error == ghl_success, "fresh number flux failed",
                        case_index, species);
      error = ghl_m1_compute_neutrino_number_flux(
          &m1_params, &nu_params, &test_case.metric, &test_case.prims,
          &test_case.state, null_rate_flux, null_rate_velocity);
      require_condition(error == ghl_success, "NULL-rate number flux failed",
                        case_index, species);
      error = ghl_m1_compute_neutrino_number_flux_from_closure(
          &m1_params, &nu_params, &test_case.metric, &test_case.prims,
          &test_case.state, &closure, supplied_flux, supplied_velocity);
      require_condition(error == ghl_success,
                        "supplied-closure number flux failed",
                        case_index, species);
      for(int i = 0; i < 3; ++i) {
        require_finite_value(fresh_flux[i], "number flux is nonfinite",
                             case_index, species);
        require_finite_value(fresh_velocity[i],
                             "transport velocity is nonfinite",
                             case_index, species);
        require_condition(m1_nearly_equal(fresh_flux[i],
                                           test_case.state.N * fresh_velocity[i],
                                           2.0e-12, 1.0e-14),
                          "number-current identity failed", case_index, species);
        require_condition(m1_nearly_equal(fresh_flux[i], null_rate_flux[i],
                                           2.0e-12, 1.0e-14) &&
                          m1_nearly_equal(fresh_velocity[i], null_rate_velocity[i],
                                          2.0e-12, 1.0e-14),
                          "repeated number transport calls disagree",
                          case_index, species);
        require_condition(m1_nearly_equal(fresh_flux[i], supplied_flux[i],
                                           2.0e-12, 1.0e-14) &&
                          m1_nearly_equal(fresh_velocity[i], supplied_velocity[i],
                                          2.0e-12, 1.0e-14),
                          "fresh and supplied-closure currents disagree",
                          case_index, species);
      }
      const double velocity_norm = m1_neutrino_seeded_metric_norm(
          test_case.metric.gammaDD, fresh_velocity);
      require_condition(velocity_norm <= 1.0 + 1.0e-10,
                        "number transport velocity is superluminal",
                        case_index, species);

      if(species == transport_reference_species)
        for(int i = 0; i < 3; ++i) {
          baseline_number_flux[i] = fresh_flux[i];
          baseline_transport_velocity[i] = fresh_velocity[i];
        }

      double perturbed_flux_from_closure[3];
      double perturbed_velocity_from_closure[3];
      error = ghl_m1_compute_neutrino_number_flux_from_closure(
          &m1_params, &nu_params, &test_case.metric, &test_case.prims,
          &test_case.perturbed_state, &perturbed_closure,
          perturbed_flux_from_closure, perturbed_velocity_from_closure);
      require_condition(error == ghl_success,
                        "perturbed supplied-closure number flux failed",
                        case_index, species);
      for(int i = 0; i < 3; ++i) {
        require_finite_value(perturbed_flux_from_closure[i],
                             "perturbed number flux is nonfinite",
                             case_index, species);
        require_finite_value(perturbed_velocity_from_closure[i],
                             "perturbed transport velocity is nonfinite",
                             case_index, species);
        if(species == transport_reference_species) {
          perturbed_number_flux[i] = perturbed_flux_from_closure[i];
          perturbed_transport_velocity[i] = perturbed_velocity_from_closure[i];
        }
      }
    }

    /* Coordinate physical number fluxes must apply the lapse/shift conversion
     * to the contravariant spatial number flux without densitization. */
    for(int direction = 0; direction < 3; ++direction) {
      double physical_number_flux = NAN;
      error = ghl_m1_compute_neutrino_physical_number_flux_from_closure(
          &m1_params, &nu_params, &test_case.metric, &test_case.prims,
          &test_case.state, &closure, (ghl_m1_direction_t)direction,
          &physical_number_flux);
      require_condition(error == ghl_success,
                        "physical number flux failed", case_index, -1);
      const double expected = test_case.metric.lapse * baseline_number_flux[direction]
                            - test_case.metric.betaU[direction] * test_case.state.N;
      require_condition(m1_nearly_equal(physical_number_flux, expected,
                                        2.0e-12, 1.0e-14),
                        "physical number-flux identity failed", case_index, -1);
    }

    const ghl_m1_direction_t direction =
        (ghl_m1_direction_t)(case_index % 3);
    double physical_E_L, physical_E_R;
    double physical_F_L[3], physical_F_R[3];
    const ghl_m1_rad_state rad_state_L =
        ghl_m1_neutrino_project_rad_state(&test_case.state);
    const ghl_m1_rad_state rad_state_R =
        ghl_m1_neutrino_project_rad_state(&test_case.perturbed_state);
    error = ghl_m1_compute_physical_flux(
        &test_case.metric, direction, &rad_state_L, &closure,
        &physical_E_L, physical_F_L);
    require_condition(error == ghl_success, "left physical flux failed",
                      case_index, -1);
    error = ghl_m1_compute_physical_flux(
        &test_case.metric, direction, &rad_state_R, &perturbed_closure,
        &physical_E_R, physical_F_R);
    require_condition(error == ghl_success, "right physical flux failed",
                      case_index, -1);
    const double physical_number_flux_L =
        test_case.metric.lapse * baseline_number_flux[direction]
        - test_case.metric.betaU[direction] * test_case.state.N;
    const double physical_number_flux_R =
        test_case.metric.lapse * perturbed_number_flux[direction]
        - test_case.metric.betaU[direction] * test_case.perturbed_state.N;
    const double speed = m1_neutrino_seeded_uniform(&rng, 0.2, 1.1);
    double flux_tilde_N = NAN;
    double flux_tilde_E = NAN;
    double flux_tilde_F[3] = {NAN, NAN, NAN};
    error = ghl_m1_compute_neutrino_rusanov_flux(
        &m1_params, &nu_params, &test_case.metric, direction,
        &test_case.state, &test_case.perturbed_state, &closure,
        &perturbed_closure, baseline_number_flux, perturbed_number_flux,
        baseline_transport_velocity, perturbed_transport_velocity, speed,
        &flux_tilde_N, &flux_tilde_E, flux_tilde_F);
    require_condition(error == ghl_success, "neutrino Rusanov flux failed",
                      case_index, -1);

    const double state_L[5] = {test_case.state.N, test_case.state.E,
                               test_case.state.F[0], test_case.state.F[1],
                               test_case.state.F[2]};
    const double state_R[5] = {test_case.perturbed_state.N,
                               test_case.perturbed_state.E,
                               test_case.perturbed_state.F[0],
                               test_case.perturbed_state.F[1],
                               test_case.perturbed_state.F[2]};
    const double physical_flux_L[5] = {physical_number_flux_L, physical_E_L,
                                       physical_F_L[0], physical_F_L[1],
                                       physical_F_L[2]};
    const double physical_flux_R[5] = {physical_number_flux_R, physical_E_R,
                                       physical_F_R[0], physical_F_R[1],
                                       physical_F_R[2]};
    const double computed_flux[5] = {flux_tilde_N, flux_tilde_E,
                                     flux_tilde_F[0], flux_tilde_F[1],
                                     flux_tilde_F[2]};
    for(int component = 0; component < 5; ++component) {
      const double expected = test_case.metric.sqrt_detgamma *
          (0.5 * (physical_flux_L[component] + physical_flux_R[component])
           - 0.5 * speed * (state_R[component] - state_L[component]));
      require_finite_value(computed_flux[component],
                           "Rusanov output is nonfinite", case_index, -1);
      require_condition(m1_nearly_equal(computed_flux[component], expected,
                                        2.0e-12, 1.0e-13),
                        "Rusanov component identity failed", case_index, -1);
    }
    double scalar_number_flux = NAN;
    error = ghl_m1_compute_number_rusanov_flux(
        test_case.state.N, test_case.perturbed_state.N,
        physical_number_flux_L, physical_number_flux_R, speed,
        &scalar_number_flux);
    require_condition(error == ghl_success &&
                      m1_nearly_equal(flux_tilde_N,
                                      test_case.metric.sqrt_detgamma *
                                          scalar_number_flux,
                                      2.0e-12, 1.0e-13),
                      "scalar number Rusanov identity failed", case_index, -1);
    require_condition(flux_tilde_N != 0.0 || test_case.state.N !=
                      test_case.perturbed_state.N,
                      "baseline/perturbed number response was not exercised",
                      case_index, -1);

    /* The wavespeed cap and all three coordinate directions are exercised for
     * each species' positive transport opacity. */
    for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
      for(int wave_direction = 0; wave_direction < 3; ++wave_direction) {
        double s_minus = NAN, s_plus = NAN;
        error = ghl_m1_compute_neutrino_wavespeeds(
            &test_case.metric, (ghl_m1_direction_t)wave_direction,
            &s_minus, &s_plus);
        require_condition(error == ghl_success, "neutrino wavespeed failed",
                          case_index, species);
        require_condition(isfinite(s_minus) && isfinite(s_plus) &&
                          s_minus <= 0.0 && s_plus >= 0.0 &&
                          s_minus <= s_plus,
                          "neutrino wavespeed contract failed",
                          case_index, species);
      }
    }

    /* Frozen-rate source, matter-coupling, backward-Euler, and lepton paths. */
    for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
      /* The single-species source API accepts only a complete one-species
       * bundle. Electron-flavor aggregate number rates with a distinct
       * charged-current subset require the paired source API; project the
       * generated bundle to its valid single-species form for these local
       * operator equivalence checks. */
      ghl_m1_neutrino_rates single_species_rates = test_case.rates[species];
      if(species != ghl_m1_neutrino_nux) {
        single_species_rates.kappa_a_N =
            single_species_rates.kappa_a_N_cc;
        single_species_rates.eta_N = single_species_rates.eta_N_cc;
      }
      ghl_m1_sources source_from_closure = {0};
      double N_source_from_closure = NAN;
      error = ghl_m1_compute_neutrino_interaction_sources_from_closure(
          &m1_params, &nu_params, &test_case.metric, &test_case.prims,
          &test_case.state, &closure, &single_species_rates,
          &source_from_closure, &N_source_from_closure);
      require_condition(error == ghl_success,
                        "supplied-closure interaction source failed",
                        case_index, species);
      require_finite_sources(&source_from_closure,
                             "interaction source is nonfinite",
                             case_index, species);
      require_finite_value(N_source_from_closure,
                           "number source is nonfinite", case_index, species);

      ghl_m1_sources source_fresh = {0};
      double N_source_fresh = NAN;
      error = ghl_m1_compute_neutrino_interaction_sources(
          &m1_params, &nu_params, &test_case.metric, &test_case.prims,
          &test_case.state, &single_species_rates,
          &source_fresh, &N_source_fresh);
      require_condition(error == ghl_success, "fresh interaction source failed",
                        case_index, species);
      require_condition(m1_nearly_equal(source_from_closure.S_E,
                                         source_fresh.S_E, 2.0e-10, 1.0e-12) &&
                        m1_nearly_equal(N_source_from_closure, N_source_fresh,
                                        2.0e-10, 1.0e-12),
                        "fresh and supplied source paths disagree",
                        case_index, species);
      for(int i = 0; i < 3; ++i)
        require_condition(m1_nearly_equal(source_from_closure.S[i],
                                           source_fresh.S[i], 2.0e-10, 1.0e-12),
                          "fresh and supplied source momentum disagree",
                          case_index, species);

      double matter_tau = NAN;
      double matter_momentum[3] = {NAN, NAN, NAN};
      error = ghl_m1_compute_neutrino_matter_coupling_sources(
          &test_case.metric, &source_from_closure, &matter_tau,
          matter_momentum);
      require_condition(error == ghl_success,
                        "matter-coupling source failed", case_index, species);
      const double alpha_sqrt_gamma =
          test_case.metric.lapse * test_case.metric.sqrt_detgamma;
      require_condition(m1_nearly_equal(matter_tau,
                                         -alpha_sqrt_gamma * source_from_closure.S_E,
                                         2.0e-12, 1.0e-13),
                        "matter energy coupling identity failed",
                        case_index, species);
      for(int i = 0; i < 3; ++i)
        require_condition(m1_nearly_equal(matter_momentum[i],
                                           -alpha_sqrt_gamma * source_from_closure.S[i],
                                           2.0e-12, 1.0e-13),
                          "matter momentum coupling identity failed",
                          case_index, species);

      double explicit_rhs_E = NAN;
      double explicit_rhs_F[3] = {NAN, NAN, NAN};
      double explicit_rhs_N = NAN;
      error = ghl_m1_compute_neutrino_explicit_rhs_sources(
          &m1_params, &nu_params, &test_case.metric,
          &test_case.metric_derivs[0], &test_case.metric_derivs[1],
          &test_case.metric_derivs[2], &test_case.curv, &test_case.prims,
          &test_case.state, &closure, true, &single_species_rates,
          &explicit_rhs_E, explicit_rhs_F, &explicit_rhs_N);
      require_condition(error == ghl_success,
                        "interaction explicit RHS failed", case_index, species);
      require_finite_value(explicit_rhs_E, "explicit RHS energy is nonfinite",
                           case_index, species);
      require_finite_value(explicit_rhs_N, "explicit RHS number is nonfinite",
                           case_index, species);
      require_condition(m1_nearly_equal(
                            explicit_rhs_E,
                            geometry_sources.S_E + alpha_sqrt_gamma *
                                source_from_closure.S_E,
                            2.0e-10, 1.0e-12) &&
                        m1_nearly_equal(explicit_rhs_N,
                                        alpha_sqrt_gamma * N_source_from_closure,
                                        2.0e-10, 1.0e-12),
                        "explicit RHS source composition failed",
                        case_index, species);
      for(int i = 0; i < 3; ++i) {
        require_finite_value(explicit_rhs_F[i],
                             "explicit RHS momentum is nonfinite",
                             case_index, species);
        require_condition(m1_nearly_equal(
                              explicit_rhs_F[i],
                              geometry_sources.S[i] + alpha_sqrt_gamma *
                                  source_from_closure.S[i],
                              2.0e-10, 1.0e-12),
                          "explicit RHS momentum composition failed",
                          case_index, species);
      }

      const double gamma_N = compute_gamma_n(
          &test_case.metric, &test_case.prims, &comoving);
      require_condition(isfinite(gamma_N) && gamma_N > 64.0 * DBL_EPSILON,
                        "generated number-current normalization is invalid",
                        case_index, species);
      const double dt_alpha = 0.05 * test_case.metric.lapse;
      double N_out = NAN;
      error = ghl_m1_update_neutrino_number_backward_euler(
          &nu_params, &single_species_rates, dt_alpha, gamma_N,
          test_case.state.N, &N_out);
      require_condition(error == ghl_success && isfinite(N_out) && N_out >= 0.0,
                        "backward-Euler number update failed",
                        case_index, species);
      const double expected_N =
          (test_case.state.N + dt_alpha * single_species_rates.eta_N) /
          (1.0 + dt_alpha * single_species_rates.kappa_a_N / gamma_N);
      require_condition(m1_nearly_equal(N_out, expected_N, 2.0e-12, 1.0e-13),
                        "backward-Euler number identity failed",
                        case_index, species);

      const double dL_rad_cc = 0.015 * test_case.rates[species].lepton_weight;
      const double baryon_density = 1.7;
      double dYe_matter = NAN;
      error = ghl_m1_compute_neutrino_lepton_increment(
          &test_case.rates[species], dL_rad_cc, baryon_density, &dYe_matter);
      require_condition(error == ghl_success &&
                        m1_nearly_equal(dYe_matter, -dL_rad_cc / baryon_density,
                                        2.0e-12, 1.0e-14),
                        "lepton increment identity failed", case_index, species);
    }

    check_transactional_contracts(&m1_params, &nu_params, &test_case,
                                  &closure, case_index);
    ++cases_run;
  }

  ghl_info("unit_test_m1_neutrino_seeded_invariants: %llu cases, "
           "%llu species cases, PRNG %s seed 0x%016llx: all tests passed\n",
           cases_run, species_runs, M1_NEUTRINO_SEEDED_PRNG_VERSION,
           (unsigned long long)M1_NEUTRINO_SEEDED_PRNG_SEED);
  return 0;
}
