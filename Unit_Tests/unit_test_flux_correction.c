#include <float.h>

#include "ghl_unit_tests.h"

/*
 * Checks the two routines that support first-order flux correction (FOFC) of
 * Lemaster & Stone (2009, ApJ 691, 1092) as used in AthenaK by Fields et al. (2025,
 * ApJS 276, 35, Sec. 3.3): ghl_assess_candidate_state, the test of an estimated
 * state, and ghl_accumulate_face_transfer, the conservative update of Eq. 10 there.
 * No downloaded fixtures are needed.
 */

static ghl_parameters params;
static ghl_eos_parameters eos;
static ghl_metric_quantities metric;
static ghl_ADM_aux_quantities metric_aux;

static void
setup(const ghl_con2prim_id_t main_routine,
      const double max_Lorentz_factor,
      const double rho_min) {

  const ghl_con2prim_id_t backups[3]
        = { ghl_con2prim_id_None, ghl_con2prim_id_None, ghl_con2prim_id_None };
  ghl_initialize_params(
        main_routine, backups, false, false, false, 1e100, max_Lorentz_factor, 0.0,
        &params);

  const double rho_ppoly[1] = { 0.0 };
  const double Gamma_ppoly[1] = { 2.0 };
  ghl_error_codes_t error = ghl_initialize_hybrid_eos_functions_and_params(
        1e-12, rho_min, 1e3, 1, rho_ppoly, Gamma_ppoly, 100.0, 2.0, &eos);
  ghl_abort_if_error(error);

  ghl_initialize_metric(
        1.1, 0.01, -0.02, 0.005, 1.3, 0.02, 0.0, 1.2, 0.01, 1.4, &metric);
  ghl_compute_ADM_auxiliaries(&metric, &metric_aux);
}

// The pressure is press_factor times the cold pressure, so a factor of one is a cold
// state and a factor of three is warm enough that no limit acts.
static void make_state(
      const double rho,
      const double press_factor,
      const double vU[3],
      ghl_primitive_quantities *restrict prims,
      ghl_conservative_quantities *restrict cons) {

  double P_cold, eps_cold;
  ghl_hybrid_compute_P_cold_and_eps_cold(&eos, rho, &P_cold, &eps_cold);
  const double press = press_factor * P_cold;
  const double eps = eps_cold + (press - P_cold) / ((eos.Gamma_th - 1.0) * rho);
  ghl_initialize_primitives(
        rho, press, eps, vU[0], vU[1], vU[2], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, prims);

  bool speed_limited = false;
  ghl_abort_if_error(
        ghl_limit_v_and_compute_u0(&params, &metric, prims, &speed_limited));
  ghl_compute_conservs(&metric, &metric_aux, prims, cons);
}

static bool assess(
      const ghl_conservative_quantities *restrict cons,
      const ghl_primitive_quantities *restrict prims_guess) {

  bool flagged = false;
  ghl_abort_if_error(ghl_assess_candidate_state(
        &params, &eos, &metric, &metric_aux, cons, prims_guess, &flagged));
  return flagged;
}

static void expect_error(
      const ghl_parameters *restrict params_in,
      const ghl_eos_parameters *restrict eos_in,
      const ghl_conservative_quantities *restrict cons,
      const ghl_primitive_quantities *restrict prims_guess,
      const ghl_error_codes_t expected) {

  bool flagged = false;
  if(ghl_assess_candidate_state(
           params_in, eos_in, &metric, &metric_aux, cons, prims_guess, &flagged)
           != expected
     || flagged) {
    ghl_error("A configuration error (code %d) was not returned unchanged\n", expected);
  }
}

static void test_assess_candidate_state(void) {
  ghl_primitive_quantities prims;
  ghl_conservative_quantities cons;
  const double moving[3] = { 0.1, 0.05, -0.02 };

  setup(ghl_con2prim_id_Noble2D, 10.0, 1e-12);
  make_state(1e-3, 3.0, moving, &prims, &cons);
  if(assess(&cons, &prims)) {
    ghl_error("A healthy candidate was flagged\n");
  }

  // The inputs are copied, not changed.
  const ghl_primitive_quantities prims_before = prims;
  const ghl_conservative_quantities cons_before = cons;
  assess(&cons, &prims);
  if(memcmp(&prims, &prims_before, sizeof prims)
     || memcmp(&cons, &cons_before, sizeof cons)) {
    ghl_error("ghl_assess_candidate_state changed its inputs\n");
  }

  // Candidates that fail recovery, or that recovery cannot reproduce.
  ghl_conservative_quantities bad = cons;
  bad.SD[0] *= 1e3;
  if(!assess(&bad, &prims)) {
    ghl_error("A candidate with an unphysically large momentum was not flagged\n");
  }
  bad = cons;
  bad.rho = -1.0;
  if(!assess(&bad, &prims)) {
    ghl_error("A candidate with negative density was not flagged\n");
  }
  for(int field = 0; field < 5; field++) {
    bad = cons;
    double *const fields[5] = { &bad.rho, &bad.tau, &bad.SD[0], &bad.SD[1], &bad.SD[2] };
    *fields[field] = field % 2 ? NAN : INFINITY;
    if(!assess(&bad, &prims)) {
      ghl_error("A candidate that is not finite in field %d was not flagged\n", field);
    }
  }

  // A primitive limit is seen through closure: the density floor raises the density by
  // one part in a million, about a hundred times the closure bound, and no diagnostic
  // records it.
  setup(ghl_con2prim_id_Noble2D, 10.0, 1.000001e-3);
  if(!assess(&cons, &prims)) {
    ghl_error("A candidate raised by the density floor was not flagged\n");
  }

  // Each of these changes the candidate by less than the closure tolerance, so only
  // the diagnostic that records it can flag the candidate.
  // 1. The conservative limiter raises tau to tau_atm.
  setup(ghl_con2prim_id_Noble2D, 10.0, 1e-12);
  // The fluid is at rest relative to the normal observers, v^i = -beta^i.
  const double at_rest[3] = { -metric.betaU[0], -metric.betaU[1], -metric.betaU[2] };
  make_state(eos.rho_atm, 1.0, at_rest, &prims, &cons);
  cons.tau = 0.999 * eos.tau_atm;
  if(!assess(&cons, &prims)) {
    ghl_error("A candidate slightly below tau_atm was not flagged\n");
  }
  // 2. The velocity limit acts on a speed barely above the cap.
  make_state(1e-3, 3.0, moving, &prims, &cons);
  setup(ghl_con2prim_id_Noble2D, (1.0 - 1e-10) * metric.lapse * prims.u0, 1e-12);
  if(!assess(&cons, &prims)) {
    ghl_error("A candidate barely above the Lorentz factor cap was not flagged\n");
  }
  // 3. Font1D replaces the energy equation by the cold EOS, which a cold state
  // satisfies.
  setup(ghl_con2prim_id_Font1D, 10.0, 1e-12);
  make_state(1e-3, 1.0, moving, &prims, &cons);
  if(!assess(&cons, &prims)) {
    ghl_error("A candidate recovered by Font1D was not flagged\n");
  }

  // A configuration error is returned and leaves the flag unchanged.
  setup(ghl_con2prim_id_Newman1D, 10.0, 1e-12);
  expect_error(&params, &eos, &cons, &prims, ghl_error_invalid_c2p_key);

  setup(ghl_con2prim_id_Noble2D, 10.0, 1e-12);
  ghl_eos_parameters unknown_eos = eos;
  unknown_eos.eos_type = (ghl_eos_t)99;
  expect_error(&params, &unknown_eos, &cons, &prims, ghl_error_unknown_eos_type);

  // Noble1D_entropy2 needs Gamma_th to equal the one-piece Gamma.
  setup(ghl_con2prim_id_Noble1D_entropy2, 10.0, 1e-12);
  ghl_eos_parameters mismatched_eos = eos;
  mismatched_eos.Gamma_th = 2.5;
  expect_error(&params, &mismatched_eos, &cons, &prims, ghl_error_invalid_eos_type);
}

// A tabulated EOS needs HDF5. Without it the configuration error is returned; with it
// the repository's coarse sample table is read, so run from the repository root.
static void test_assess_tabulated(void) {
  const ghl_con2prim_id_t backups[3]
        = { ghl_con2prim_id_None, ghl_con2prim_id_None, ghl_con2prim_id_None };
  ghl_parameters tabulated_params;
  ghl_initialize_params(
        ghl_con2prim_id_Palenzuela1D, backups, false, true, false, 1e100, 10.0, 0.0,
        &tabulated_params);
  ghl_eos_parameters tabulated_eos = { 0 };
  ghl_primitive_quantities prims;
  ghl_conservative_quantities cons;

#ifdef GHL_DISABLE_HDF5
  tabulated_eos.eos_type = ghl_eos_tabulated;
  ghl_initialize_primitives(
        1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 1.0, &prims);
  ghl_initialize_conservatives(1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.3, &cons);
  expect_error(
        &tabulated_params, &tabulated_eos, &cons, &prims, ghl_error_used_disabled_hdf5);
#else
  bool flagged = false;
  ghl_abort_if_error(ghl_initialize_tabulated_eos_functions_and_params(
        "Unit_Tests/sample_table/"
        "Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817_simple.h5",
        1e-14, -1, -1, 0.3, -1, -1, 0.1, -1, -1, &tabulated_eos));

  const double rho = sqrt(tabulated_eos.rho_min * tabulated_eos.rho_max);
  const double T = sqrt(tabulated_eos.T_min * tabulated_eos.T_max);
  double press, eps;
  ghl_abort_if_error(
        ghl_tabulated_compute_P_eps_from_T(&tabulated_eos, rho, 0.3, T, &press, &eps));
  ghl_initialize_primitives(
        rho, press, eps, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, T, &prims);
  bool speed_limited = false;
  ghl_abort_if_error(
        ghl_limit_v_and_compute_u0(&tabulated_params, &metric, &prims, &speed_limited));
  ghl_compute_conservs(&metric, &metric_aux, &prims, &cons);

  ghl_abort_if_error(ghl_assess_candidate_state(
        &tabulated_params, &tabulated_eos, &metric, &metric_aux, &cons, &prims,
        &flagged));
  if(flagged) {
    ghl_error("A healthy tabulated candidate was flagged\n");
  }

  // Palenzuela takes Y_e from the guess, so only the Y_e closure comparison sees this.
  cons.Y_e *= 1.01;
  ghl_abort_if_error(ghl_assess_candidate_state(
        &tabulated_params, &tabulated_eos, &metric, &metric_aux, &cons, &prims,
        &flagged));
  if(!flagged) {
    ghl_error("A tabulated candidate with inconsistent Y_e was not flagged\n");
  }
  ghl_tabulated_free_memory(&tabulated_eos);
#endif
}

static void test_accumulate_face_transfer(void) {
  // Powers of two keep every update exact, so the sum of the two cells is conserved.
  ghl_conservative_quantities transfer, left, right;
  ghl_initialize_conservatives(0.25, -0.5, 1.5, -2.5, 0.125, 0.0625, 0.75, &transfer);
  ghl_initialize_conservatives(1.0, 2.0, 3.0, 4.0, 5.0, 8.0, 0.5, &left);
  ghl_initialize_conservatives(-1.0, 0.5, 0.25, 0.75, -3.0, 2.0, 4.0, &right);
  const ghl_conservative_quantities left_before = left;
  const ghl_conservative_quantities right_before = right;

  ghl_abort_if_error(ghl_accumulate_face_transfer(&transfer, &left, &right));
  if(left.rho != left_before.rho - transfer.rho
     || right.rho != right_before.rho + transfer.rho
     || left.tau != left_before.tau - transfer.tau
     || right.tau != right_before.tau + transfer.tau
     || left.Y_e != left_before.Y_e - transfer.Y_e
     || right.Y_e != right_before.Y_e + transfer.Y_e
     || left.entropy != left_before.entropy - transfer.entropy
     || right.entropy != right_before.entropy + transfer.entropy
     || left.SD[0] != left_before.SD[0] - transfer.SD[0]
     || right.SD[0] != right_before.SD[0] + transfer.SD[0]
     || left.SD[1] != left_before.SD[1] - transfer.SD[1]
     || right.SD[1] != right_before.SD[1] + transfer.SD[1]
     || left.SD[2] != left_before.SD[2] - transfer.SD[2]
     || right.SD[2] != right_before.SD[2] + transfer.SD[2]) {
    ghl_error("ghl_accumulate_face_transfer did not apply the transfer\n");
  }

  // A transfer that is not finite in any field changes nothing.
  for(int field = 0; field < 7; field++) {
    ghl_conservative_quantities bad = transfer;
    double *const fields[7] = { &bad.rho,   &bad.tau,   &bad.Y_e,  &bad.entropy,
                                &bad.SD[0], &bad.SD[1], &bad.SD[2] };
    *fields[field] = field % 2 ? NAN : INFINITY;
    left = left_before;
    right = right_before;
    if(ghl_accumulate_face_transfer(&bad, &left, &right)
             != ghl_error_invalid_face_transfer
       || memcmp(&left, &left_before, sizeof left)
       || memcmp(&right, &right_before, sizeof right)) {
      ghl_error("A non-finite face transfer was not rejected cleanly\n");
    }
  }
}

int main(void) {
  test_assess_candidate_state();
  test_assess_tabulated();
  test_accumulate_face_transfer();
  ghl_info("Flux correction routines test has passed!\n");
  return 0;
}
