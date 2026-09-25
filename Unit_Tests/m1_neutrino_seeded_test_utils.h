#ifndef UNIT_TESTS_M1_NEUTRINO_SEEDED_TEST_UTILS_H_
#define UNIT_TESTS_M1_NEUTRINO_SEEDED_TEST_UTILS_H_

/*
 * Local data construction for unit_test_m1_neutrino_seeded_invariants.c.
 *
 * The generator is deliberately versioned and deterministic.  It only emits
 * states that are admissible by construction: spatial metrics are generated
 * as L L^T, fluid velocities are scaled with the generated metric, and
 * radiation fluxes are scaled to a selected M1 cone fraction.  This helper is
 * a local unit-test data source; it is not a THC_M1 reference or equivalence
 * fixture.
 */

#include <stdint.h>

#include "m1_test_utils.h"

#define M1_NEUTRINO_SEEDED_PRNG_VERSION "splitmix64-v1"
#define M1_NEUTRINO_SEEDED_PRNG_SEED    UINT64_C(0x7d4f2b91a6c3e805)
#define M1_NEUTRINO_SEEDED_CASE_COUNT   256

typedef struct {
  uint64_t state;
} m1_neutrino_seeded_rng;

typedef struct {
  ghl_metric_quantities metric;
  ghl_primitive_quantities prims;
  ghl_m1_neutrino_state state;
  ghl_m1_neutrino_state perturbed_state;
  ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
  ghl_metric_quantities metric_derivs[3];
  ghl_extrinsic_curvature curv;
} m1_neutrino_seeded_case;

static inline uint64_t
m1_neutrino_seeded_next_u64(m1_neutrino_seeded_rng *restrict rng) {

  /* SplitMix64-v1: the state advance and finalizer are part of the replay
   * contract. Unsigned overflow is defined by C99 for uint64_t. */
  uint64_t z = (rng->state += UINT64_C(0x9e3779b97f4a7c15));
  z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
  return z ^ (z >> 31);
}

static inline double m1_neutrino_seeded_unit(m1_neutrino_seeded_rng *restrict rng) {

  return (double)(m1_neutrino_seeded_next_u64(rng) >> 11) * 0x1.0p-53;
}

static inline double m1_neutrino_seeded_uniform(
      m1_neutrino_seeded_rng *restrict rng,
      const double lower,
      const double upper) {

  return lower + (upper - lower) * m1_neutrino_seeded_unit(rng);
}

static inline double m1_neutrino_seeded_metric_norm(
      const double metric_tensor[3][3],
      const double vector[3]) {

  long double norm_squared = 0.0L;
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      norm_squared += (long double)vector[i] * metric_tensor[i][j] * vector[j];
    }
  }
  return sqrt((double)norm_squared);
}

static inline void m1_neutrino_seeded_make_metric(
      m1_neutrino_seeded_rng *restrict rng,
      const bool off_diagonal,
      ghl_metric_quantities *restrict metric) {

  *metric = (ghl_metric_quantities){ 0 };
  metric->lapse = m1_neutrino_seeded_uniform(rng, 0.78, 1.12);
  metric->lapseinv = 1.0 / metric->lapse;
  metric->lapseinv2 = metric->lapseinv * metric->lapseinv;
  for(int i = 0; i < 3; ++i) {
    metric->betaU[i] = m1_neutrino_seeded_uniform(rng, -0.12, 0.12);
  }

  if(!off_diagonal) {
    m1_setup_flat_metric(metric);
    metric->lapse = m1_neutrino_seeded_uniform(rng, 0.78, 1.12);
    metric->lapseinv = 1.0 / metric->lapse;
    metric->lapseinv2 = metric->lapseinv * metric->lapseinv;
    for(int i = 0; i < 3; ++i) {
      metric->betaU[i] = m1_neutrino_seeded_uniform(rng, -0.12, 0.12);
    }
    return;
  }

  /* A lower-triangular Cholesky factor gives an SPD metric without needing
   * rejection sampling. The inverse is formed analytically from the same
   * symmetric tensor so gammaDD/gammaUU remain a coherent pair. */
  const double L00 = m1_neutrino_seeded_uniform(rng, 0.88, 1.12);
  const double L10 = m1_neutrino_seeded_uniform(rng, -0.18, 0.18);
  const double L11 = m1_neutrino_seeded_uniform(rng, 0.88, 1.12);
  const double L20 = m1_neutrino_seeded_uniform(rng, -0.18, 0.18);
  const double L21 = m1_neutrino_seeded_uniform(rng, -0.18, 0.18);
  const double L22 = m1_neutrino_seeded_uniform(rng, 0.88, 1.12);
  const double L[3][3] = { { L00, 0.0, 0.0 }, { L10, L11, 0.0 }, { L20, L21, L22 } };

  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      double value = 0.0;
      for(int k = 0; k < 3; ++k) {
        value += L[i][k] * L[j][k];
      }
      metric->gammaDD[i][j] = value;
    }
  }

  const double a = metric->gammaDD[0][0];
  const double b = metric->gammaDD[0][1];
  const double c = metric->gammaDD[0][2];
  const double d = metric->gammaDD[1][1];
  const double e = metric->gammaDD[1][2];
  const double f = metric->gammaDD[2][2];
  const double determinant
        = a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d);

  metric->gammaUU[0][0] = (d * f - e * e) / determinant;
  metric->gammaUU[0][1] = (c * e - b * f) / determinant;
  metric->gammaUU[0][2] = (b * e - c * d) / determinant;
  metric->gammaUU[1][0] = metric->gammaUU[0][1];
  metric->gammaUU[1][1] = (a * f - c * c) / determinant;
  metric->gammaUU[1][2] = (b * c - a * e) / determinant;
  metric->gammaUU[2][0] = metric->gammaUU[0][2];
  metric->gammaUU[2][1] = metric->gammaUU[1][2];
  metric->gammaUU[2][2] = (a * d - b * b) / determinant;
  metric->detgamma = determinant;
  metric->sqrt_detgamma = sqrt(determinant);
}

static inline void m1_neutrino_seeded_make_direction(
      m1_neutrino_seeded_rng *restrict rng,
      const double metric_tensor[3][3],
      double direction[3]) {

  for(int i = 0; i < 3; ++i) {
    direction[i] = m1_neutrino_seeded_uniform(rng, -1.0, 1.0);
  }
  double norm = m1_neutrino_seeded_metric_norm(metric_tensor, direction);
  if(!(norm > 0.0)) {
    direction[0] = 1.0;
    direction[1] = 0.0;
    direction[2] = 0.0;
    norm = m1_neutrino_seeded_metric_norm(metric_tensor, direction);
  }
  for(int i = 0; i < 3; ++i) {
    direction[i] /= norm;
  }
}

static inline void m1_neutrino_seeded_make_primitives(
      m1_neutrino_seeded_rng *restrict rng,
      const bool moving,
      const ghl_metric_quantities *restrict metric,
      ghl_primitive_quantities *restrict prims) {

  *prims = (ghl_primitive_quantities){ 0 };
  prims->rho = m1_neutrino_seeded_uniform(rng, 0.4, 2.5);
  prims->eps = m1_neutrino_seeded_uniform(rng, 0.05, 0.6);
  prims->press = prims->rho * prims->eps * m1_neutrino_seeded_uniform(rng, 0.8, 1.4);
  prims->Y_e = m1_neutrino_seeded_uniform(rng, 0.05, 0.45);
  prims->temperature = m1_neutrino_seeded_uniform(rng, 0.2, 4.0);
  prims->entropy = m1_neutrino_seeded_uniform(rng, 0.05, 1.0);

  double eulerian_velocity[3] = { 0.0, 0.0, 0.0 };
  if(moving) {
    double direction[3];
    m1_neutrino_seeded_make_direction(rng, metric->gammaDD, direction);
    const double speed = m1_neutrino_seeded_uniform(rng, 0.05, 0.45);
    for(int i = 0; i < 3; ++i) {
      eulerian_velocity[i] = speed * direction[i];
    }
  }

  const double velocity_norm
        = m1_neutrino_seeded_metric_norm(metric->gammaDD, eulerian_velocity);
  const double W = 1.0 / sqrt((1.0 - velocity_norm) * (1.0 + velocity_norm));
  prims->u0 = W / metric->lapse;
  for(int i = 0; i < 3; ++i) {
    prims->vU[i] = metric->lapse * eulerian_velocity[i] - metric->betaU[i];
  }
}

static inline void m1_neutrino_seeded_make_state(
      m1_neutrino_seeded_rng *restrict rng,
      const bool anisotropic,
      const ghl_metric_quantities *restrict metric,
      ghl_m1_neutrino_state *restrict state) {

  *state = (ghl_m1_neutrino_state){ 0 };
  state->N = m1_neutrino_seeded_uniform(rng, 0.02, 2.0);
  state->E = m1_neutrino_seeded_uniform(rng, 0.05, 5.0);

  if(anisotropic) {
    double direction[3];
    m1_neutrino_seeded_make_direction(rng, metric->gammaUU, direction);
    const double flux_factor = m1_neutrino_seeded_uniform(rng, 0.08, 0.88);
    for(int i = 0; i < 3; ++i) {
      state->F[i] = state->E * flux_factor * direction[i];
    }
  }
}

static inline void m1_neutrino_seeded_make_perturbed_state(
      m1_neutrino_seeded_rng *restrict rng,
      const ghl_m1_neutrino_state *restrict state,
      ghl_m1_neutrino_state *restrict perturbed_state) {

  /* Scale E and F together so the perturbation remains in the same M1 cone;
   * perturb N independently to exercise the number-current path. */
  *perturbed_state = *state;
  const double ef_scale = m1_neutrino_seeded_uniform(rng, 0.85, 1.15);
  const double n_scale = m1_neutrino_seeded_uniform(rng, 0.8, 1.2);
  perturbed_state->E *= ef_scale;
  perturbed_state->N *= n_scale;
  for(int i = 0; i < 3; ++i) {
    perturbed_state->F[i] *= ef_scale;
  }
}

static inline void m1_neutrino_seeded_make_rates(
      m1_neutrino_seeded_rng *restrict rng,
      const ghl_m1_neutrino_species_t species,
      ghl_m1_neutrino_rates *restrict rates) {

  *rates = (ghl_m1_neutrino_rates){ 0 };
  rates->species = species;
  rates->n_eq = m1_neutrino_seeded_uniform(rng, 0.04, 1.6);
  rates->mean_energy = m1_neutrino_seeded_uniform(rng, 0.5, 18.0);
  rates->J_eq = rates->n_eq * rates->mean_energy;
  rates->kappa_a_N = m1_neutrino_seeded_uniform(rng, 0.01, 0.35);
  rates->kappa_a_E = m1_neutrino_seeded_uniform(rng, 0.01, 0.35);
  rates->kappa_s = m1_neutrino_seeded_uniform(rng, 0.01, 0.35);
  rates->kappa_tr = rates->kappa_a_E + rates->kappa_s;
  rates->eta_N = rates->kappa_a_N * rates->n_eq;
  rates->eta_E = rates->kappa_a_E * rates->J_eq;

  switch(species) {
    case ghl_m1_neutrino_nue:
      rates->lepton_weight = 1.0;
      rates->kappa_a_N_cc
            = rates->kappa_a_N * m1_neutrino_seeded_uniform(rng, 0.15, 0.85);
      rates->eta_N_cc = rates->kappa_a_N_cc * rates->n_eq;
      break;
    case ghl_m1_neutrino_anue:
      rates->lepton_weight = -1.0;
      rates->kappa_a_N_cc
            = rates->kappa_a_N * m1_neutrino_seeded_uniform(rng, 0.15, 0.85);
      rates->eta_N_cc = rates->kappa_a_N_cc * rates->n_eq;
      break;
    case ghl_m1_neutrino_nux:
      rates->lepton_weight = 0.0;
      rates->kappa_a_N_cc = 0.0;
      rates->eta_N_cc = 0.0;
      break;
    default:
      rates->lepton_weight = 0.0;
      rates->kappa_a_N_cc = 0.0;
      rates->eta_N_cc = 0.0;
      break;
  }
}

static inline void m1_neutrino_seeded_make_geometry_data(
      m1_neutrino_seeded_rng *restrict rng,
      ghl_metric_quantities metric_derivs[3],
      ghl_extrinsic_curvature *restrict curv) {

  for(int direction = 0; direction < 3; ++direction) {
    metric_derivs[direction] = (ghl_metric_quantities){ 0 };
    metric_derivs[direction].lapse = m1_neutrino_seeded_uniform(rng, -0.15, 0.15);
    for(int j = 0; j < 3; ++j) {
      metric_derivs[direction].betaU[j] = m1_neutrino_seeded_uniform(rng, -0.15, 0.15);
      for(int k = 0; k < 3; ++k) {
        metric_derivs[direction].gammaDD[j][k]
              = m1_neutrino_seeded_uniform(rng, -0.15, 0.15);
      }
    }
  }

  *curv = (ghl_extrinsic_curvature){ 0 };
  for(int i = 0; i < 3; ++i) {
    for(int j = i; j < 3; ++j) {
      const double value = m1_neutrino_seeded_uniform(rng, -0.15, 0.15);
      curv->K[i][j] = value;
      curv->K[j][i] = value;
    }
  }
}

static inline void m1_neutrino_seeded_make_case(
      m1_neutrino_seeded_rng *restrict rng,
      const int case_index,
      m1_neutrino_seeded_case *restrict test_case) {

  const bool moving = (case_index & 1) != 0;
  const bool off_diagonal = ((case_index >> 1) & 1) != 0;
  const bool anisotropic = ((case_index >> 2) & 1) != 0;

  m1_neutrino_seeded_make_metric(rng, off_diagonal, &test_case->metric);
  m1_neutrino_seeded_make_primitives(rng, moving, &test_case->metric, &test_case->prims);
  m1_neutrino_seeded_make_state(rng, anisotropic, &test_case->metric, &test_case->state);
  m1_neutrino_seeded_make_perturbed_state(
        rng, &test_case->state, &test_case->perturbed_state);
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    m1_neutrino_seeded_make_rates(
          rng, (ghl_m1_neutrino_species_t)species, &test_case->rates[species]);
  }
  m1_neutrino_seeded_make_geometry_data(rng, test_case->metric_derivs, &test_case->curv);
}

static inline double m1_neutrino_seeded_flux_factor(
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_neutrino_state *restrict state) {

  const double norm = m1_neutrino_seeded_metric_norm(metric->gammaUU, state->F);
  return norm / state->E;
}

static inline void
m1_neutrino_seeded_default_parameters(ghl_m1_neutrino_parameters *restrict nu_params) {

  *nu_params = (ghl_m1_neutrino_parameters){ 0 };
  nu_params->N_floor = 1.0e-12;
}

#endif // UNIT_TESTS_M1_NEUTRINO_SEEDED_TEST_UTILS_H_
