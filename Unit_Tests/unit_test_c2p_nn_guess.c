#ifdef __APPLE__
#define _DARWIN_C_SOURCE
#else
#define _XOPEN_SOURCE 700
#endif

#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ghl_con2prim.h"

#ifndef GHL_DISABLE_HDF5
#include <errno.h>
#include <hdf5.h>
#include <unistd.h>
#endif

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#define CHECK(cond, ...)                                                       \
  do {                                                                        \
    if(!(cond)) {                                                             \
      fprintf(stderr, __VA_ARGS__);                                           \
      fprintf(stderr, "\n");                                                  \
      exit(1);                                                                \
    }                                                                         \
  } while(0)

#define CHECK_ERROR(got, expected)                                             \
  CHECK((got) == (expected), "expected error %d but got %d",                  \
        (int)(expected), (int)(got))

static float invrng(const float lo, const float hi) {
  return 1.0f / (hi - lo);
}

static void set_valid_model(
      ghl_c2p_nn_model *restrict model,
      int *restrict x_kind,
      float *restrict x_lo,
      float *restrict x_hi,
      float *restrict x_invrng,
      int *restrict out_kind,
      float *restrict out_lo,
      float *restrict out_hi,
      float *restrict out_invrng,
      float *restrict W_in,
      float *restrict b_in,
      float *restrict W_hid,
      float *restrict b_hid,
      float *restrict W_out,
      float *restrict b_out) {
  memset(model, 0, sizeof(*model));
  model->in_dim = 4;
  model->hidden_dim = 2;
  model->n_hidden = 2;
  model->out_dim = 2;
  model->q_idx = 0;
  model->s_idx = 2;
  model->x_eps = 1e-30f;
  model->y_eps = 1e-6f;
  model->dx_eps = 1e-5f;
  model->x_kind = x_kind;
  model->x_lo = x_lo;
  model->x_hi = x_hi;
  model->x_invrng = x_invrng;
  model->out_kind = out_kind;
  model->out_lo = out_lo;
  model->out_hi = out_hi;
  model->out_invrng = out_invrng;
  model->W_in = W_in;
  model->b_in = b_in;
  model->W_hid = W_hid;
  model->b_hid = b_hid;
  model->W_out = W_out;
  model->b_out = b_out;

  x_kind[0] = 0;
  x_kind[1] = 1;
  x_kind[2] = 0;
  x_kind[3] = 1;
  x_lo[0] = 0.0f;
  x_lo[1] = -2.0f;
  x_lo[2] = 0.0f;
  x_lo[3] = -1.0f;
  x_hi[0] = 10.0f;
  x_hi[1] = 2.0f;
  x_hi[2] = 10.0f;
  x_hi[3] = 1.0f;
  for(int i = 0; i < 4; ++i) {
    x_invrng[i] = invrng(x_lo[i], x_hi[i]);
  }

  out_kind[0] = 0;
  out_kind[1] = 1;
  out_lo[0] = 0.0f;
  out_hi[0] = 1.0f;
  out_invrng[0] = 1.0f;
  out_lo[1] = 2.0f;
  out_hi[1] = 6.0f;
  out_invrng[1] = invrng(out_lo[1], out_hi[1]);

  for(int i = 0; i < 8; ++i) {
    W_in[i] = 0.0f;
    W_out[i] = 0.0f;
  }
  b_in[0] = 0.25f;
  b_in[1] = -0.5f;
  for(int i = 0; i < 4; ++i) {
    W_hid[i] = 0.0f;
  }
  b_hid[0] = 0.125f;
  b_hid[1] = -0.25f;
  b_out[0] = 0.0f;
  b_out[1] = -1.0986122886681098f;
}

static ghl_c2p_nn_model valid_stack_model(void) {
  static int x_kind[4];
  static float x_lo[4], x_hi[4], x_invrng[4];
  static int out_kind[2];
  static float out_lo[2], out_hi[2], out_invrng[2];
  static float W_in[8], b_in[2], W_hid[4], b_hid[2], W_out[8], b_out[2];
  ghl_c2p_nn_model model;
  set_valid_model(&model, x_kind, x_lo, x_hi, x_invrng,
                  out_kind, out_lo, out_hi, out_invrng,
                  W_in, b_in, W_hid, b_hid, W_out, b_out);
  return model;
}

static void test_validate_model(void) {
  ghl_c2p_nn_model model = valid_stack_model();
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model), ghl_success);
  CHECK_ERROR(ghl_c2p_nn_validate_model(NULL), ghl_error_nn_c2p_model_is_null);

  model = valid_stack_model();
  model.in_dim = 3;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_invalid_dimensions);

  model = valid_stack_model();
  model.q_idx = model.in_dim;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_invalid_input_index);

  model = valid_stack_model();
  model.q_idx = 1;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_invalid_input_index);

  model = valid_stack_model();
  model.s_idx = 0;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_invalid_input_index);

  model = valid_stack_model();
  model.y_eps = 0.5f;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_invalid_number);

  model = valid_stack_model();
  model.x_lo = NULL;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_missing_array);

  model = valid_stack_model();
  model.x_kind[1] = 7;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_invalid_kind);

  model = valid_stack_model();
  model.x_invrng[0] = 0.0f;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_invalid_number);

  model = valid_stack_model();
  model.out_kind[0] = 1;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_invalid_kind);

  model = valid_stack_model();
  model.b_hid[1] = NAN;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_invalid_number);
}

static void test_guess_model(void) {
  ghl_c2p_nn_model model = valid_stack_model();
  ghl_nn_c2p_input_t input = { 2.0f, 0.25f, 0.5f, 0.1f };
  ghl_nn_c2p_guess_t guess = ghl_c2p_nn_guess(&model, input);
  CHECK(fabsf(guess.x - 4.0f) < 1e-6f,
        "x-bounded guess mismatch: %.9g", guess.x);

  model = valid_stack_model();
  model.out_kind[1] = 2;
  model.out_lo[1] = logf(2.0f);
  model.out_hi[1] = logf(8.0f);
  model.out_invrng[1] = invrng(model.out_lo[1], model.out_hi[1]);
  model.b_out[1] = 0.0f;
  guess = ghl_c2p_nn_guess(&model, input);
  CHECK(fabsf(guess.x - 4.0f) < 1e-6f,
        "x guess with extra log-linear output mismatch: %.9g", guess.x);

  model = valid_stack_model();
  input.r = NAN;
  guess = ghl_c2p_nn_guess(&model, input);
  CHECK(fabsf(guess.x - 4.0f) < 1e-6f,
        "finite fallback x mismatch: %.9g", guess.x);

  model = valid_stack_model();
  model.x_kind[1] = 9;
  input.r = 0.25f;
  guess = ghl_c2p_nn_guess(&model, input);
  CHECK(fabsf(guess.x - 4.0f) < 1e-6f,
        "invalid transform fallback x mismatch: %.9g", guess.x);

  guess = ghl_c2p_nn_guess(NULL, input);
  CHECK(guess.x == 0.0f, "NULL model fallback failed");

  model = valid_stack_model();
  input.q = FLT_MAX;
  input.s = 0.0f;
  input.r = 0.0f;
  input.t = 0.0f;
  guess = ghl_c2p_nn_guess(&model, input);
  CHECK(isfinite(guess.x) && guess.x == 0.0f,
        "overflowing midpoint did not return the finite sentinel");

  model = valid_stack_model();
  model.W_in[0] = 0.2f;
  model.W_in[1] = 0.4f;
  model.W_in[2] = 0.6f;
  model.W_in[3] = 0.8f;
  model.b_in[0] = 0.0f;
  model.W_hid[0] = 1.0f;
  model.b_hid[0] = 0.0f;
  model.W_out[0] = 1.0f;
  model.b_out[0] = 0.0f;
  input = (ghl_nn_c2p_input_t){ 2.0f, 0.25f, 0.5f, 0.1f };
  const float r_scaled = (log10f(input.r) + 2.0f) * 0.25f;
  const float t_scaled = (log10f(input.t) + 1.0f) * 0.5f;
  const float hidden = 0.2f * input.q * 0.1f
                     + 0.4f * r_scaled
                     + 0.6f * input.s * 0.1f
                     + 0.8f * t_scaled;
  const float expected_x = 1.0f + input.q - input.s
                         + (1.0f + input.q) / (1.0f + expf(-hidden));
  guess = ghl_c2p_nn_guess(&model, input);
  CHECK(fabsf(guess.x - expected_x) < 1e-6f,
        "fixed {q,r,s,t} feature ordering mismatch: %.9g vs %.9g",
        guess.x, expected_x);
}

static void fake_enforce_bounds(
      const ghl_eos_parameters *restrict eos,
      double *restrict rho,
      double *restrict Y_e,
      double *restrict eps) {
  (void)eos;
  (void)rho;
  (void)Y_e;
  (void)eps;
}

static ghl_error_codes_t fake_compute_P_S_T(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double eps,
      double *restrict press,
      double *restrict entropy,
      double *restrict temperature) {
  (void)eos;
  *press = rho * (1.0 + eps);
  *entropy = Y_e + eps;
  *temperature = 1.0 + eps;
  return ghl_success;
}

static ghl_error_codes_t fake_compute_P_S_T_nonfinite(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double eps,
      double *restrict press,
      double *restrict entropy,
      double *restrict temperature) {
  (void)eos;
  (void)rho;
  (void)Y_e;
  (void)eps;
  *press = NAN;
  *entropy = 0.0;
  *temperature = 1.0;
  return ghl_success;
}

static ghl_error_codes_t fake_compute_P_S_T_error(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double eps,
      double *restrict press,
      double *restrict entropy,
      double *restrict temperature) {
  (void)eos;
  (void)rho;
  (void)Y_e;
  (void)eps;
  (void)press;
  (void)entropy;
  (void)temperature;
  return ghl_error_table_bisection;
}

static void check_atmosphere_guess(
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric,
      const double BU[3],
      const ghl_primitive_quantities *restrict prims) {
  CHECK(prims->rho == eos->rho_atm && prims->press == eos->press_atm
        && prims->eps == eos->eps_atm && prims->entropy == eos->entropy_atm
        && prims->Y_e == eos->Y_e_atm && prims->temperature == eos->T_atm,
        "tabulated guess did not return atmosphere thermodynamics");
  CHECK(prims->vU[0] == -metric->betaU[0]
        && prims->vU[1] == -metric->betaU[1]
        && prims->vU[2] == -metric->betaU[2]
        && prims->u0 == metric->lapseinv,
        "tabulated guess did not return atmosphere velocity");
  CHECK(prims->BU[0] == BU[0] && prims->BU[1] == BU[1] && prims->BU[2] == BU[2],
        "tabulated guess did not preserve magnetic fields");
}

static void test_public_primitive_guess_helper(void) {
  const ghl_con2prim_id_t backups[3] = {
    ghl_con2prim_id_None, ghl_con2prim_id_None, ghl_con2prim_id_None
  };
  ghl_parameters params;
  ghl_initialize_params(
        ghl_con2prim_id_None, backups, false, false, false, 1e100, 20.0, 0.0,
        &params);

  ghl_metric_quantities metric;
  ghl_initialize_metric(
        0.8, 0.03, -0.02, 0.01,
        1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);

  ghl_eos_parameters eos = { 0 };
  eos.eos_type = ghl_eos_tabulated;
  eos.rho_atm = 0.01;
  eos.press_atm = 0.02;
  eos.eps_atm = 0.03;
  eos.entropy_atm = 0.04;
  eos.Y_e_atm = 0.05;
  eos.T_atm = 0.06;
  eos.T_max = 10.0;

  ghl_c2p_nn_model model = valid_stack_model();
  eos.c2p_nn = &model;
  ghl_tabulated_enforce_bounds_rho_Ye_eps = fake_enforce_bounds;
  ghl_tabulated_compute_P_S_T_from_eps = fake_compute_P_S_T;

  ghl_conservative_quantities cons;
  ghl_initialize_conservatives(1.0, 1.0, 0.05, -0.03, 0.02, 0.0, 0.2, &cons);
  ghl_primitive_quantities prims;
  ghl_initialize_primitives(
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.1, -0.2, 0.3, 0.0, 0.0, 0.0, &prims);
  const double BU[3] = { prims.BU[0], prims.BU[1], prims.BU[2] };
  ghl_c2p_nn_guess_primitives(&params, &eos, &metric, &cons, &prims);
  CHECK(isfinite(prims.rho) && prims.rho > 0.0 && prims.rho != eos.rho_atm
        && isfinite(prims.press) && isfinite(prims.u0),
        "valid public NN primitive guess did not use the completion path");

  ghl_guess_primitives(&params, &eos, &metric, &cons, &prims);
  CHECK(isfinite(prims.rho) && prims.rho > 0.0 && prims.rho != eos.rho_atm
        && isfinite(prims.press) && isfinite(prims.u0),
        "valid default tabulated guess did not use the completion path");

  ghl_initialize_conservatives(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, &cons);
  ghl_initialize_primitives(
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        sqrt(1.5), 0.0, 0.0, 0.0, 0.0, 0.0, &prims);
  ghl_tabulated_primitive_guess_aux zero_spanning_aux;
  ghl_tabulated_compute_primitive_guess_auxiliaries(
        &metric, &cons, &prims, &zero_spanning_aux);
  CHECK(1.0 + zero_spanning_aux.q - zero_spanning_aux.s < 0.0
        && 2.0 + 2.0 * zero_spanning_aux.q - zero_spanning_aux.s > 0.0,
        "NN test state does not span zero in its admissible x bracket");
  const double zero_spanning_BU[3] = { prims.BU[0], prims.BU[1], prims.BU[2] };
  ghl_c2p_nn_guess_primitives(&params, &eos, &metric, &cons, &prims);
  check_atmosphere_guess(&eos, &metric, zero_spanning_BU, &prims);

  ghl_initialize_conservatives(1.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.2, &cons);
  ghl_initialize_primitives(
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        sqrt(2.0), 0.0, 0.0, 0.0, 0.0, 0.0, &prims);
  const double negative_BU[3] = { prims.BU[0], prims.BU[1], prims.BU[2] };
  ghl_c2p_nn_guess_primitives(&params, &eos, &metric, &cons, &prims);
  CHECK(prims.rho != eos.rho_atm && isfinite(prims.rho)
        && isfinite(prims.press) && isfinite(prims.u0),
        "usable negative NN x was rejected by sign alone");

  ghl_tabulated_primitive_guess_aux negative_aux = { 0 };
  ghl_initialize_conservatives(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, &cons);
  ghl_tabulated_primitive_guess_from_x(
        &params, &eos, &metric, &cons, &negative_aux, -2.0, &prims);
  CHECK(prims.rho != eos.rho_atm && isfinite(prims.rho)
        && isfinite(prims.press) && isfinite(prims.u0),
        "usable negative x was rejected by sign alone");

  ghl_tabulated_primitive_guess_from_x(
        &params, &eos, &metric, &cons, &negative_aux, 0.0, &prims);
  check_atmosphere_guess(&eos, &metric, negative_BU, &prims);

  ghl_tabulated_enforce_bounds_rho_Ye_eps = NULL;
  ghl_tabulated_primitive_guess_from_x(
        &params, &eos, &metric, &cons, &negative_aux, -2.0, &prims);
  check_atmosphere_guess(&eos, &metric, negative_BU, &prims);
  ghl_tabulated_enforce_bounds_rho_Ye_eps = fake_enforce_bounds;

  ghl_tabulated_compute_P_S_T_from_eps = NULL;
  ghl_tabulated_primitive_guess_from_x(
        &params, &eos, &metric, &cons, &negative_aux, -2.0, &prims);
  check_atmosphere_guess(&eos, &metric, negative_BU, &prims);

  ghl_tabulated_compute_P_S_T_from_eps = fake_compute_P_S_T_error;
  ghl_tabulated_primitive_guess_from_x(
        &params, &eos, &metric, &cons, &negative_aux, -2.0, &prims);
  check_atmosphere_guess(&eos, &metric, negative_BU, &prims);
  ghl_tabulated_compute_P_S_T_from_eps = fake_compute_P_S_T;

  prims.BU[0] = BU[0];
  prims.BU[1] = BU[1];
  prims.BU[2] = BU[2];
  eos.c2p_nn = NULL;
  ghl_c2p_nn_guess_primitives(&params, &eos, &metric, &cons, &prims);
  check_atmosphere_guess(&eos, &metric, BU, &prims);

  cons.rho = 0.0;
  ghl_guess_primitives(&params, &eos, &metric, &cons, &prims);
  check_atmosphere_guess(&eos, &metric, BU, &prims);

  ghl_tabulated_primitive_guess_aux aux = { 0 };
  cons.rho = 1.0;
  aux.q = aux.r = aux.s = aux.t = 1.0;
  aux.B_squared = 1.0;
  aux.BdotS = 1.0;
  aux.SU[0] = aux.SU[1] = aux.SU[2] = 1.0;
  ghl_tabulated_primitive_guess_from_x(
        &params, &eos, &metric, &cons, &aux, 1e200, &prims);
  check_atmosphere_guess(&eos, &metric, BU, &prims);

  ghl_tabulated_primitive_guess_from_x(
        &params, &eos, &metric, &cons, &aux, DBL_MIN, &prims);
  check_atmosphere_guess(&eos, &metric, BU, &prims);

  aux.r = DBL_MAX;
  aux.q = aux.s = aux.t = 0.0;
  aux.B_squared = aux.BdotS = 0.0;
  aux.SU[0] = aux.SU[1] = aux.SU[2] = 0.0;
  ghl_tabulated_primitive_guess_from_x(
        &params, &eos, &metric, &cons, &aux, pow(DBL_MIN, 0.25), &prims);
  check_atmosphere_guess(&eos, &metric, BU, &prims);

  aux.q = aux.r = aux.s = aux.t = 1.0;
  aux.B_squared = aux.BdotS = 1.0;
  aux.SU[0] = aux.SU[1] = aux.SU[2] = 1.0;
  cons.Y_e = DBL_MAX;
  cons.rho = DBL_MIN;
  ghl_tabulated_primitive_guess_from_x(
        &params, &eos, &metric, &cons, &aux, 2.0, &prims);
  check_atmosphere_guess(&eos, &metric, BU, &prims);
  cons.rho = 1.0;
  cons.Y_e = 0.2;

  ghl_tabulated_compute_P_S_T_from_eps = fake_compute_P_S_T_nonfinite;
  ghl_tabulated_primitive_guess_from_x(
        &params, &eos, &metric, &cons, &aux, 2.0, &prims);
  check_atmosphere_guess(&eos, &metric, BU, &prims);

  ghl_tabulated_compute_P_S_T_from_eps = fake_compute_P_S_T;
  aux.q = aux.r = aux.s = aux.t = 0.0;
  aux.B_squared = aux.BdotS = 0.0;
  aux.SU[0] = DBL_MAX;
  aux.SU[1] = aux.SU[2] = 0.0;
  ghl_tabulated_primitive_guess_from_x(
        &params, &eos, &metric, &cons, &aux, 1e-70, &prims);
  check_atmosphere_guess(&eos, &metric, BU, &prims);

  aux.SU[0] = 1e200;
  ghl_tabulated_primitive_guess_from_x(
        &params, &eos, &metric, &cons, &aux, 1.0, &prims);
  check_atmosphere_guess(&eos, &metric, BU, &prims);
}

#ifdef GHL_DISABLE_HDF5
static void test_disabled_direct_tabulated_solvers(void) {
  ghl_parameters params = { 0 };
  ghl_eos_parameters eos = { 0 };
  ghl_metric_quantities metric = { 0 };
  ghl_ADM_aux_quantities metric_aux = { 0 };
  ghl_conservative_quantities cons = { 0 };
  ghl_primitive_quantities prims = { .rho = 1.0, .press = 2.0 };
  ghl_con2prim_diagnostics diagnostics = { .tau_fix = true, .n_iter = 7 };
  const ghl_primitive_quantities prims_before = prims;
  const ghl_con2prim_diagnostics diagnostics_before = diagnostics;

#define CHECK_DISABLED_SOLVER(name)                                          \
  CHECK_ERROR(name(&params, &eos, &metric, &metric_aux, &cons,              \
                   &prims, &diagnostics), ghl_error_used_disabled_hdf5)
  CHECK_DISABLED_SOLVER(ghl_tabulated_Noble2D);
  CHECK_DISABLED_SOLVER(ghl_tabulated_Palenzuela1D_energy);
  CHECK_DISABLED_SOLVER(ghl_tabulated_Palenzuela1D_entropy);
  CHECK_DISABLED_SOLVER(ghl_tabulated_Newman1D_energy);
  CHECK_DISABLED_SOLVER(ghl_tabulated_Newman1D_entropy);
#undef CHECK_DISABLED_SOLVER

  CHECK(memcmp(&prims, &prims_before, sizeof(prims)) == 0,
        "disabled direct solver mutated primitives");
  CHECK(memcmp(&diagnostics, &diagnostics_before, sizeof(diagnostics)) == 0,
        "disabled direct solver mutated diagnostics");
}
#endif

#ifndef GHL_DISABLE_HDF5
static char nn_test_directory[PATH_MAX];
static char nn_original_directory[PATH_MAX];
static bool nn_test_directory_active;

static const char *const nn_test_files[] = {
  "unit_test_c2p_nn_preserved.h5",
  "unit_test_c2p_nn_root.h5",
  "unit_test_c2p_nn_embedded.h5",
  "unit_test_c2p_nn_legacy.h5",
  "unit_test_c2p_nn_missing_scalar.h5",
  "unit_test_c2p_nn_scalar_bad_rank.h5",
  "unit_test_c2p_nn_scalar_bad_type.h5",
  "unit_test_c2p_nn_invalid_dims.h5",
  "unit_test_c2p_nn_missing_array.h5",
  "unit_test_c2p_nn_array_bad_rank.h5",
  "unit_test_c2p_nn_array_bad_size.h5",
  "unit_test_c2p_nn_array_bad_type.h5",
  "unit_test_c2p_nn_missing_out_scaling.h5",
  "unit_test_c2p_nn_validation_failure.h5"
};

static int cleanup_hdf5_test_directory(void) {
  if(!nn_test_directory_active) {
    return 0;
  }
  if(chdir(nn_original_directory) != 0) {
    fprintf(stderr, "failed to leave NN test directory %s: %s\n",
            nn_test_directory, strerror(errno));
    return 1;
  }
  int failed = 0;
  for(size_t i = 0; i < sizeof(nn_test_files)/sizeof(nn_test_files[0]); ++i) {
    char path[PATH_MAX];
    const int written = snprintf(
          path, sizeof(path), "%s/%s", nn_test_directory, nn_test_files[i]);
    if(written < 0 || (size_t)written >= sizeof(path)) {
      fprintf(stderr, "failed to form cleanup path for NN test file %s\n",
              nn_test_files[i]);
      failed = 1;
      continue;
    }
    if(unlink(path) != 0 && errno != ENOENT) {
      fprintf(stderr, "failed to remove NN test file %s: %s\n", path,
              strerror(errno));
      failed = 1;
    }
  }

  if(rmdir(nn_test_directory) != 0) {
    fprintf(stderr, "failed to remove NN test directory %s: %s\n",
            nn_test_directory, strerror(errno));
    failed = 1;
  }
  else {
    nn_test_directory_active = false;
  }
  return failed;
}

static void cleanup_hdf5_test_directory_at_exit(void) {
  (void)cleanup_hdf5_test_directory();
}

static void setup_hdf5_test_directory(void) {
  CHECK(getcwd(nn_original_directory, sizeof(nn_original_directory)) != NULL,
        "failed to record original test directory: %s", strerror(errno));
  const char *tmpdir = getenv("TMPDIR");
  if(tmpdir == NULL || tmpdir[0] == '\0') {
    tmpdir = "/tmp";
  }
  const int written = snprintf(
        nn_test_directory, sizeof(nn_test_directory), "%s%sunit_test_c2p_nn_XXXXXX",
        tmpdir, tmpdir[strlen(tmpdir) - 1] == '/' ? "" : "/");
  CHECK(written > 0 && (size_t)written < sizeof(nn_test_directory),
        "NN test directory template is too long");
  if(mkdtemp(nn_test_directory) == NULL) {
    fprintf(stderr, "failed to create NN test directory: %s\n", strerror(errno));
    exit(1);
  }
  nn_test_directory_active = true;
  if(chdir(nn_test_directory) != 0) {
    fprintf(stderr, "failed to enter NN test directory %s: %s\n",
            nn_test_directory, strerror(errno));
    if(rmdir(nn_test_directory) != 0) {
      fprintf(stderr, "failed to remove unused NN test directory %s: %s\n",
              nn_test_directory, strerror(errno));
    }
    exit(1);
  }
  if(atexit(cleanup_hdf5_test_directory_at_exit) != 0) {
    fprintf(stderr, "failed to register NN test cleanup\n");
    if(chdir(nn_original_directory) != 0) {
      fprintf(stderr, "failed to leave NN test directory %s: %s\n",
              nn_test_directory, strerror(errno));
    }
    else if(rmdir(nn_test_directory) != 0) {
      fprintf(stderr, "failed to remove unused NN test directory %s: %s\n",
              nn_test_directory, strerror(errno));
    }
    exit(1);
  }
}

static void create_group_checked(hid_t file_id, const char *name) {
  hid_t group_id = H5Gcreate2(file_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  CHECK(group_id >= 0, "failed to create HDF5 group %s", name);
  CHECK(H5Gclose(group_id) >= 0, "failed to close HDF5 group %s", name);
}

static void write_scalar_i32(hid_t file_id, const char *name, int value) {
  hid_t space_id = H5Screate(H5S_SCALAR);
  CHECK(space_id >= 0, "failed to create scalar dataspace");
  hid_t dataset_id = H5Dcreate2(file_id, name, H5T_NATIVE_INT, space_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  CHECK(dataset_id >= 0, "failed to create HDF5 dataset %s", name);
  CHECK(H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, &value) >= 0,
        "failed to write HDF5 dataset %s", name);
  H5Dclose(dataset_id);
  H5Sclose(space_id);
}

static void write_scalar_f32(hid_t file_id, const char *name, float value) {
  hid_t space_id = H5Screate(H5S_SCALAR);
  CHECK(space_id >= 0, "failed to create scalar dataspace");
  hid_t dataset_id = H5Dcreate2(file_id, name, H5T_NATIVE_FLOAT, space_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  CHECK(dataset_id >= 0, "failed to create HDF5 dataset %s", name);
  CHECK(H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, &value) >= 0,
        "failed to write HDF5 dataset %s", name);
  H5Dclose(dataset_id);
  H5Sclose(space_id);
}

static void write_scalar_string(hid_t file_id, const char *name) {
  const char value[] = "bad";
  hid_t string_type = H5Tcopy(H5T_C_S1);
  CHECK(string_type >= 0, "failed to create string datatype");
  CHECK(H5Tset_size(string_type, sizeof(value)) >= 0,
        "failed to set string datatype size");
  hid_t space_id = H5Screate(H5S_SCALAR);
  CHECK(space_id >= 0, "failed to create scalar dataspace");
  hid_t dataset_id = H5Dcreate2(file_id, name, string_type, space_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  CHECK(dataset_id >= 0, "failed to create HDF5 dataset %s", name);
  CHECK(H5Dwrite(dataset_id, string_type, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, value) >= 0,
        "failed to write HDF5 dataset %s", name);
  H5Dclose(dataset_id);
  H5Sclose(space_id);
  H5Tclose(string_type);
}

static void write_array_i32(
      hid_t file_id,
      const char *name,
      int rank,
      const hsize_t *restrict dims,
      const int *restrict values) {
  hid_t space_id = H5Screate_simple(rank, dims, NULL);
  CHECK(space_id >= 0, "failed to create array dataspace");
  hid_t dataset_id = H5Dcreate2(file_id, name, H5T_NATIVE_INT, space_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  CHECK(dataset_id >= 0, "failed to create HDF5 dataset %s", name);
  CHECK(H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, values) >= 0,
        "failed to write HDF5 dataset %s", name);
  H5Dclose(dataset_id);
  H5Sclose(space_id);
}

static void write_array_f32(
      hid_t file_id,
      const char *name,
      int rank,
      const hsize_t *restrict dims,
      const float *restrict values) {
  hid_t space_id = H5Screate_simple(rank, dims, NULL);
  CHECK(space_id >= 0, "failed to create array dataspace");
  hid_t dataset_id = H5Dcreate2(file_id, name, H5T_NATIVE_FLOAT, space_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  CHECK(dataset_id >= 0, "failed to create HDF5 dataset %s", name);
  CHECK(H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, values) >= 0,
        "failed to write HDF5 dataset %s", name);
  H5Dclose(dataset_id);
  H5Sclose(space_id);
}

static void write_array_string(
      hid_t file_id,
      const char *name,
      int rank,
      const hsize_t *restrict dims) {
  const char values[4][4] = { "bad", "bad", "bad", "bad" };
  hid_t string_type = H5Tcopy(H5T_C_S1);
  CHECK(string_type >= 0, "failed to create string datatype");
  CHECK(H5Tset_size(string_type, sizeof(values[0])) >= 0,
        "failed to set string datatype size");
  hid_t space_id = H5Screate_simple(rank, dims, NULL);
  CHECK(space_id >= 0, "failed to create array dataspace");
  hid_t dataset_id = H5Dcreate2(file_id, name, string_type, space_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  CHECK(dataset_id >= 0, "failed to create HDF5 dataset %s", name);
  CHECK(H5Dwrite(dataset_id, string_type, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, values) >= 0,
        "failed to write HDF5 dataset %s", name);
  H5Dclose(dataset_id);
  H5Sclose(space_id);
  H5Tclose(string_type);
}

static void prefixed_name(
      const char *restrict prefix,
      const char *restrict name,
      char *restrict out,
      size_t out_size) {
  const int written = (prefix == NULL || prefix[0] == '\0')
                    ? snprintf(out, out_size, "%s", name)
                    : snprintf(out, out_size, "%s/%s", prefix, name);
  CHECK(written > 0 && (size_t)written < out_size,
        "HDF5 path buffer too small");
}

static void create_nn_groups(hid_t file_id, const char *prefix) {
  const char *groups[] = {
    "dims", "meta", "scaling", "layers"
  };
  if(prefix != NULL && prefix[0] != '\0') {
    create_group_checked(file_id, prefix);
  }
  for(size_t i = 0; i < sizeof(groups) / sizeof(groups[0]); ++i) {
    char name[256];
    prefixed_name(prefix, groups[i], name, sizeof(name));
    create_group_checked(file_id, name);
  }
}

static void write_valid_hdf5_model(
      hid_t file_id,
      const char *prefix,
      const bool include_out_scaling,
      const int out_dim) {
  create_nn_groups(file_id, prefix);

  char name[256];
  prefixed_name(prefix, "dims/in_dim", name, sizeof(name));
  write_scalar_i32(file_id, name, 4);
  prefixed_name(prefix, "dims/hidden_dim", name, sizeof(name));
  write_scalar_i32(file_id, name, 2);
  prefixed_name(prefix, "dims/n_hidden", name, sizeof(name));
  write_scalar_i32(file_id, name, 2);
  prefixed_name(prefix, "dims/out_dim", name, sizeof(name));
  write_scalar_i32(file_id, name, out_dim);
  prefixed_name(prefix, "meta/q_idx", name, sizeof(name));
  write_scalar_i32(file_id, name, 0);
  prefixed_name(prefix, "meta/s_idx", name, sizeof(name));
  write_scalar_i32(file_id, name, 2);
  prefixed_name(prefix, "scaling/x_eps", name, sizeof(name));
  write_scalar_f32(file_id, name, 1e-30f);
  prefixed_name(prefix, "meta/y_eps", name, sizeof(name));
  write_scalar_f32(file_id, name, 1e-6f);
  prefixed_name(prefix, "meta/dx_eps", name, sizeof(name));
  write_scalar_f32(file_id, name, 1e-5f);

  const hsize_t dims_in[1] = { 4 };
  const hsize_t dims_hidden[1] = { 2 };
  const hsize_t dims_out[1] = { (hsize_t)out_dim };
  const hsize_t dims_w_in[2] = { 2, 4 };
  const hsize_t dims_b_hid[2] = { 1, 2 };
  const hsize_t dims_w_hid[3] = { 1, 2, 2 };
  const hsize_t dims_w_out[2] = { (hsize_t)out_dim, 2 };
  const int x_kind[4] = { 0, 1, 0, 1 };
  const float x_lo[4] = { 0.0f, -2.0f, 0.0f, -1.0f };
  const float x_hi[4] = { 10.0f, 2.0f, 10.0f, 1.0f };
  const float x_invrng[4] = { 0.1f, 0.25f, 0.1f, 0.5f };
  const int out_kind[2] = { 0, 1 };
  const float out_lo[2] = { 0.0f, 2.0f };
  const float out_hi[2] = { 1.0f, 6.0f };
  const float out_invrng[2] = { 1.0f, 0.25f };
  const float W_in[8] = { 0 };
  const float b_in[2] = { 0.25f, -0.5f };
  const float W_hid[4] = { 0 };
  const float b_hid[2] = { 0.125f, -0.25f };
  const float W_out[4] = { 0 };
  const float b_out[2] = { 0.0f, -1.0986122886681098f };

  prefixed_name(prefix, "scaling/x_kind", name, sizeof(name));
  write_array_i32(file_id, name, 1, dims_in, x_kind);
  prefixed_name(prefix, "scaling/x_lo", name, sizeof(name));
  write_array_f32(file_id, name, 1, dims_in, x_lo);
  prefixed_name(prefix, "scaling/x_hi", name, sizeof(name));
  write_array_f32(file_id, name, 1, dims_in, x_hi);
  prefixed_name(prefix, "scaling/x_invrng", name, sizeof(name));
  write_array_f32(file_id, name, 1, dims_in, x_invrng);
  if(include_out_scaling) {
    prefixed_name(prefix, "scaling/out_kind", name, sizeof(name));
    write_array_i32(file_id, name, 1, dims_out, out_kind);
    prefixed_name(prefix, "scaling/out_lo", name, sizeof(name));
    write_array_f32(file_id, name, 1, dims_out, out_lo);
    prefixed_name(prefix, "scaling/out_hi", name, sizeof(name));
    write_array_f32(file_id, name, 1, dims_out, out_hi);
    prefixed_name(prefix, "scaling/out_invrng", name, sizeof(name));
    write_array_f32(file_id, name, 1, dims_out, out_invrng);
  }
  prefixed_name(prefix, "layers/W_in", name, sizeof(name));
  write_array_f32(file_id, name, 2, dims_w_in, W_in);
  prefixed_name(prefix, "layers/b_in", name, sizeof(name));
  write_array_f32(file_id, name, 1, dims_hidden, b_in);
  prefixed_name(prefix, "layers/W_hid", name, sizeof(name));
  write_array_f32(file_id, name, 3, dims_w_hid, W_hid);
  prefixed_name(prefix, "layers/b_hid", name, sizeof(name));
  write_array_f32(file_id, name, 2, dims_b_hid, b_hid);
  prefixed_name(prefix, "layers/W_out", name, sizeof(name));
  write_array_f32(file_id, name, 2, dims_w_out, W_out);
  prefixed_name(prefix, "layers/b_out", name, sizeof(name));
  write_array_f32(file_id, name, 1, dims_out, b_out);
}

static void create_hdf5_file(
      const char *restrict path,
      const char *restrict prefix,
      const bool include_out_scaling,
      const int out_dim) {
  hid_t file_id = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  CHECK(file_id >= 0, "failed to create HDF5 file %s", path);
  write_valid_hdf5_model(file_id, prefix, include_out_scaling, out_dim);
  CHECK(H5Fclose(file_id) >= 0, "failed to close HDF5 file %s", path);
}

static void replace_dataset_i32(
      const char *restrict path,
      const char *restrict name,
      const int value) {
  hid_t file_id = H5Fopen(path, H5F_ACC_RDWR, H5P_DEFAULT);
  CHECK(file_id >= 0, "failed to open HDF5 file %s", path);
  CHECK(H5Ldelete(file_id, name, H5P_DEFAULT) >= 0,
        "failed to delete HDF5 dataset %s", name);
  write_scalar_i32(file_id, name, value);
  CHECK(H5Fclose(file_id) >= 0, "failed to close HDF5 file %s", path);
}

static void delete_dataset(
      const char *restrict path,
      const char *restrict name) {
  hid_t file_id = H5Fopen(path, H5F_ACC_RDWR, H5P_DEFAULT);
  CHECK(file_id >= 0, "failed to open HDF5 file %s", path);
  CHECK(H5Ldelete(file_id, name, H5P_DEFAULT) >= 0,
        "failed to delete HDF5 dataset %s", name);
  CHECK(H5Fclose(file_id) >= 0, "failed to close HDF5 file %s", path);
}

static void replace_dataset_with_scalar_string(
      const char *restrict path,
      const char *restrict name) {
  hid_t file_id = H5Fopen(path, H5F_ACC_RDWR, H5P_DEFAULT);
  CHECK(file_id >= 0, "failed to open HDF5 file %s", path);
  CHECK(H5Ldelete(file_id, name, H5P_DEFAULT) >= 0,
        "failed to delete HDF5 dataset %s", name);
  write_scalar_string(file_id, name);
  CHECK(H5Fclose(file_id) >= 0, "failed to close HDF5 file %s", path);
}

static void replace_dataset_with_i32_array(
      const char *restrict path,
      const char *restrict name,
      const int rank,
      const hsize_t *restrict dims,
      const int *restrict values) {
  hid_t file_id = H5Fopen(path, H5F_ACC_RDWR, H5P_DEFAULT);
  CHECK(file_id >= 0, "failed to open HDF5 file %s", path);
  CHECK(H5Ldelete(file_id, name, H5P_DEFAULT) >= 0,
        "failed to delete HDF5 dataset %s", name);
  write_array_i32(file_id, name, rank, dims, values);
  CHECK(H5Fclose(file_id) >= 0, "failed to close HDF5 file %s", path);
}

static void replace_dataset_with_string_array(
      const char *restrict path,
      const char *restrict name,
      const int rank,
      const hsize_t *restrict dims) {
  hid_t file_id = H5Fopen(path, H5F_ACC_RDWR, H5P_DEFAULT);
  CHECK(file_id >= 0, "failed to open HDF5 file %s", path);
  CHECK(H5Ldelete(file_id, name, H5P_DEFAULT) >= 0,
        "failed to delete HDF5 dataset %s", name);
  write_array_string(file_id, name, rank, dims);
  CHECK(H5Fclose(file_id) >= 0, "failed to close HDF5 file %s", path);
}

static void check_hdf5_load_failure_preserves_model(
      const char *restrict path,
      const ghl_error_codes_t expected_error) {
  ghl_eos_parameters eos = { 0 };
  create_hdf5_file("unit_test_c2p_nn_preserved.h5", "", true, 2);
  CHECK_ERROR(ghl_c2p_nn_load_hdf5("unit_test_c2p_nn_preserved.h5", &eos),
              ghl_success);
  ghl_c2p_nn_model *const preserved_model = eos.c2p_nn;
  CHECK(preserved_model != NULL, "initial NN HDF5 load returned NULL model");
  CHECK_ERROR(ghl_c2p_nn_load_hdf5(path, &eos), expected_error);
  CHECK(eos.c2p_nn == preserved_model,
        "failed HDF5 load should preserve existing model");
  ghl_c2p_nn_free(eos.c2p_nn);
}

static void test_hdf5_loaders(void) {
  const ghl_nn_c2p_input_t input = { 2.0f, 0.25f, 0.5f, 0.1f };
  const int out_kind_linear = 1;

  const char root_path[] = "unit_test_c2p_nn_root.h5";
  const char embedded_path[] = "unit_test_c2p_nn_embedded.h5";
  const char legacy_path[] = "unit_test_c2p_nn_legacy.h5";

  create_hdf5_file(root_path, "", true, 2);
  ghl_eos_parameters eos = { 0 };
  CHECK_ERROR(ghl_c2p_nn_load_hdf5(root_path, &eos), ghl_success);
  CHECK(eos.c2p_nn != NULL, "direct NN HDF5 load returned NULL model");
  CHECK(eos.c2p_nn->out_dim == 2, "direct HDF5 out_dim mismatch");
  CHECK(eos.c2p_nn->out_kind[1] == out_kind_linear,
        "direct HDF5 output-kind metadata mismatch");
  ghl_nn_c2p_guess_t guess = ghl_c2p_nn_guess(eos.c2p_nn, input);
  CHECK(fabsf(guess.x - 4.0f) < 1e-6f, "direct HDF5 x mismatch");
  ghl_c2p_nn_free(eos.c2p_nn);
  eos.c2p_nn = NULL;

  create_hdf5_file(embedded_path, "grhayl_nn_c2p", true, 2);
  CHECK_ERROR(ghl_c2p_nn_load_from_eos_hdf5(embedded_path, &eos),
              ghl_success);
  CHECK(eos.c2p_nn != NULL, "embedded NN HDF5 load returned NULL model");
  CHECK(eos.c2p_nn->out_dim == 2, "embedded HDF5 out_dim mismatch");
  CHECK(eos.c2p_nn->out_kind[1] == out_kind_linear,
        "embedded HDF5 output-kind metadata mismatch");
  guess = ghl_c2p_nn_guess(eos.c2p_nn, input);
  CHECK(fabsf(guess.x - 4.0f) < 1e-6f, "embedded HDF5 x mismatch");
  ghl_c2p_nn_free(eos.c2p_nn);
  eos.c2p_nn = NULL;

  create_hdf5_file(legacy_path, "", false, 1);
  CHECK_ERROR(ghl_c2p_nn_load_hdf5(legacy_path, &eos), ghl_success);
  CHECK(eos.c2p_nn != NULL, "legacy NN HDF5 load returned NULL model");
  CHECK(eos.c2p_nn->out_dim == 1, "legacy HDF5 out_dim mismatch");
  guess = ghl_c2p_nn_guess(eos.c2p_nn, input);
  CHECK(fabsf(guess.x - 4.0f) < 1e-6f, "legacy HDF5 x mismatch");
  ghl_c2p_nn_free(eos.c2p_nn);
  eos.c2p_nn = NULL;
}

static void test_hdf5_loader_error_paths(void) {
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  const hsize_t one_dim[1] = { 1 };
  const hsize_t three_dims[1] = { 3 };
  const hsize_t four_dims[1] = { 4 };
  const int one_value[1] = { 4 };
  const int three_values[3] = { 0, 1, 0 };

  const char missing_scalar[] = "unit_test_c2p_nn_missing_scalar.h5";
  create_hdf5_file(missing_scalar, "", true, 2);
  delete_dataset(missing_scalar, "dims/in_dim");
  check_hdf5_load_failure_preserves_model(
        missing_scalar, ghl_error_hdf5_dataset_could_not_open);

  const char scalar_bad_rank[] = "unit_test_c2p_nn_scalar_bad_rank.h5";
  create_hdf5_file(scalar_bad_rank, "", true, 2);
  replace_dataset_with_i32_array(
        scalar_bad_rank, "dims/in_dim", 1, one_dim, one_value);
  check_hdf5_load_failure_preserves_model(
        scalar_bad_rank, ghl_error_hdf5_dataset_invalid_ndims);

  const char scalar_bad_type[] = "unit_test_c2p_nn_scalar_bad_type.h5";
  create_hdf5_file(scalar_bad_type, "", true, 2);
  replace_dataset_with_scalar_string(scalar_bad_type, "dims/in_dim");
  check_hdf5_load_failure_preserves_model(
        scalar_bad_type, ghl_error_hdf5_dataset_could_not_read);

  const char invalid_dims[] = "unit_test_c2p_nn_invalid_dims.h5";
  create_hdf5_file(invalid_dims, "", true, 2);
  replace_dataset_i32(invalid_dims, "dims/in_dim", 9);
  check_hdf5_load_failure_preserves_model(
        invalid_dims, ghl_error_nn_c2p_invalid_dimensions);

  const char missing_array[] = "unit_test_c2p_nn_missing_array.h5";
  create_hdf5_file(missing_array, "", true, 2);
  delete_dataset(missing_array, "scaling/x_kind");
  check_hdf5_load_failure_preserves_model(
        missing_array, ghl_error_hdf5_dataset_could_not_open);

  const char array_bad_rank[] = "unit_test_c2p_nn_array_bad_rank.h5";
  create_hdf5_file(array_bad_rank, "", true, 2);
  replace_dataset_i32(array_bad_rank, "scaling/x_kind", 0);
  check_hdf5_load_failure_preserves_model(
        array_bad_rank, ghl_error_hdf5_dataset_invalid_ndims);

  const char array_bad_size[] = "unit_test_c2p_nn_array_bad_size.h5";
  create_hdf5_file(array_bad_size, "", true, 2);
  replace_dataset_with_i32_array(
        array_bad_size, "scaling/x_kind", 1, three_dims, three_values);
  check_hdf5_load_failure_preserves_model(
        array_bad_size, ghl_error_hdf5_dataset_size_mismatch);

  const char array_bad_type[] = "unit_test_c2p_nn_array_bad_type.h5";
  create_hdf5_file(array_bad_type, "", true, 2);
  replace_dataset_with_string_array(
        array_bad_type, "scaling/x_kind", 1, four_dims);
  check_hdf5_load_failure_preserves_model(
        array_bad_type, ghl_error_hdf5_dataset_could_not_read);

  const char missing_out_scaling[] = "unit_test_c2p_nn_missing_out_scaling.h5";
  create_hdf5_file(missing_out_scaling, "", false, 2);
  check_hdf5_load_failure_preserves_model(
        missing_out_scaling, ghl_error_hdf5_dataset_could_not_open);

  const char validation_failure[] = "unit_test_c2p_nn_validation_failure.h5";
  create_hdf5_file(validation_failure, "", true, 2);
  replace_dataset_i32(validation_failure, "meta/q_idx", 4);
  check_hdf5_load_failure_preserves_model(
        validation_failure, ghl_error_nn_c2p_invalid_input_index);
  replace_dataset_i32(validation_failure, "meta/q_idx", 0);
  replace_dataset_i32(validation_failure, "meta/s_idx", 0);
  check_hdf5_load_failure_preserves_model(
        validation_failure, ghl_error_nn_c2p_invalid_input_index);
}
#endif

int main(void) {
  test_validate_model();
  test_guess_model();
  test_public_primitive_guess_helper();
#ifdef GHL_DISABLE_HDF5
  test_disabled_direct_tabulated_solvers();
#endif
#ifndef GHL_DISABLE_HDF5
  setup_hdf5_test_directory();
  test_hdf5_loaders();
  test_hdf5_loader_error_paths();
  CHECK(cleanup_hdf5_test_directory() == 0,
        "failed to clean NN HDF5 test artifacts");
#endif
  printf("All c2p neural-network guess tests succeeded\n");
  return 0;
}
