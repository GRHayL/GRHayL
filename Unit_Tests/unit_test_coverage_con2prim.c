#define _POSIX_C_SOURCE 200809L
#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#elif defined(__APPLE__) && !defined(_DARWIN_C_SOURCE)
#define _DARWIN_C_SOURCE
#endif

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ghl_unit_tests.h"
#include "ghl_con2prim.h"
#include "../GRHayL/Con2Prim/roots.h"
#ifndef GHL_DISABLE_HDF5
#include <hdf5.h>
#endif

/*
 * Focused line-coverage tests for Con2Prim guards that existing
 * reference-data suites never drive:
 *
 *   - Brent root-finding fast paths in check_a_b_compute_fa_fb: exact endpoint
 *     root hits and the already-tightened-bracket acceptance, all reached
 *     through the public ghl_brent entry point.
 *   - Hybrid multi-method's skip-a-None-backup early break in the
 *     ghl_con2prim_hybrid_multi_method backup loop, reached by configuring
 *     every backup slot to ghl_con2prim_id_None while the main routine fails.
 *
 * Each case asserts the expected return code and the exact resulting root or
 * diagnostics state; no internal helpers are reimplemented or stubbed.
 */

static int g_fail;
#define CHECK(c)                                                        \
  do {                                                                  \
    if(!(c)) {                                                          \
      fprintf(stderr, "FAIL %s:%d  %s\n", __FILE__, __LINE__, #c);      \
      g_fail = 1;                                                       \
    }                                                                   \
  } while(0)

static double g_root_target;
static double f_linear(
      const double x,
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_conservative_quantities *restrict cons_undens,
      fparams_struct *restrict fparams,
      ghl_primitive_quantities *restrict prims) {
  (void)params; (void)eos; (void)cons_undens; (void)fparams; (void)prims;
  return x - g_root_target;
}

static void test_brent_endpoint_roots_and_tiny_interval(void) {
  ghl_parameters params;
  ghl_eos_parameters eos;
  ghl_conservative_quantities cons;
  ghl_primitive_quantities prims;
  fparams_struct fparams;
  memset(&params, 0, sizeof(params));
  memset(&eos, 0, sizeof(eos));
  memset(&cons, 0, sizeof(cons));
  memset(&prims, 0, sizeof(prims));
  memset(&fparams, 0, sizeof(fparams));

  struct { const char *label; double a, b, root_at, tol, expect_root; } cases[] = {
    /* f(a)==0 exactly: root=a with residual 0, no iterations. */
    { "endpoint-a", 0.0, 1.0, 0.0, 1e-12, 0.0 },
    /* f(b)==0 exactly: root=b with residual 0, no iterations. */
    { "endpoint-b", 0.0, 1.0, 1.0, 1e-12, 1.0 },
    /* f(a)*f(b)<0 with bracket tighter than tol: brent accepts the closest
     * endpoint b as the root, again without iterating. */
    { "tiny-interval", 0.49, 0.50, 0.495, 1.0, 0.50 },
  };
  for(unsigned k = 0; k < sizeof(cases)/sizeof(cases[0]); k++) {
    g_root_target = cases[k].root_at;
    roots_params r;
    memset(&r, 0, sizeof(r));
    r.tol = cases[k].tol;
    r.max_iters = 100;
    const ghl_error_codes_t err = ghl_brent(
          f_linear, &params, &eos, &cons, &fparams, &prims,
          cases[k].a, cases[k].b, &r);
    CHECK(err == ghl_success);
    CHECK(r.root == cases[k].expect_root);
    if(k == 2) {
      CHECK(fabs(r.residual - (cases[k].b - cases[k].root_at)) < 1e-15);
    } else {
      CHECK(r.residual == 0.0);
    }
  }
}

static void test_hybrid_multi_method_backup_none_break(void) {
  ghl_parameters params;
  ghl_eos_parameters eos;
  ghl_metric_quantities metric;
  ghl_ADM_aux_quantities metric_aux;
  ghl_conservative_quantities cons;
  ghl_primitive_quantities prims;
  ghl_con2prim_diagnostics diag;

  const double rho_ppoly[1] = { 0.0 };
  const double Gamma_ppoly[1] = { 2.0 };
  const ghl_error_codes_t setup_error = ghl_initialize_hybrid_eos_functions_and_params(
        1e-12, 1e-12, 1e300, 1, rho_ppoly, Gamma_ppoly, 1.0, 2.0, &eos);
  CHECK(setup_error == ghl_success);
  if(setup_error != ghl_success) return;
  memset(&metric_aux, 0, sizeof(metric_aux));

  const ghl_con2prim_id_t none_backups[3]
        = { ghl_con2prim_id_None, ghl_con2prim_id_None, ghl_con2prim_id_None };
  ghl_initialize_params(ghl_con2prim_id_None, none_backups,
                        false, false, false, 1e8, 60.0, 0.0, &params);
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0,
                        &metric);
  ghl_initialize_conservatives(1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, &cons);
  ghl_initialize_primitives(1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, &prims);
  ghl_initialize_diagnostics(&diag);

  const ghl_error_codes_t err = ghl_con2prim_hybrid_multi_method(
        &params, &eos, &metric, &metric_aux, &cons, &prims, &diag);
  /* main=None fails in select_method; backup loop must bail out on the first
   * None slot without ever setting any backup-flag. */
  CHECK(err == ghl_error_invalid_c2p_key);
  CHECK(diag.backup[0] == false);
  CHECK(diag.backup[1] == false);
  CHECK(diag.backup[2] == false);
}


#ifndef GHL_DISABLE_HDF5

static void write_i32_scalar(const hid_t obj, const char *const name, const int value) {
  const hid_t space = H5Screate(H5S_SCALAR);
  const hid_t dset = H5Dcreate2(obj, name, H5T_NATIVE_INT, space,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
  H5Dclose(dset);
  H5Sclose(space);
}

static void write_f32_scalar(const hid_t obj, const char *const name, const float value) {
  const hid_t space = H5Screate(H5S_SCALAR);
  const hid_t dset = H5Dcreate2(obj, name, H5T_NATIVE_FLOAT, space,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
  H5Dclose(dset);
  H5Sclose(space);
}

static void write_f32_array(const hid_t obj, const char *const name,
                            const int rank, const hsize_t *const dims,
                            const float fill) {
  hsize_t total = 1;
  for(int i = 0; i < rank; i++) total *= dims[i];
  const hid_t space = H5Screate_simple(rank, dims, NULL);
  const hid_t dset = H5Dcreate2(obj, name, H5T_NATIVE_FLOAT, space,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  float *const data = malloc((size_t)total * sizeof(*data));
  for(hsize_t i = 0; i < total; i++) data[i] = fill;
  H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  free(data);
  H5Dclose(dset);
  H5Sclose(space);
}

static void write_i32_array(const hid_t obj, const char *const name,
                            const int rank, const hsize_t *const dims,
                            const int fill) {
  hsize_t total = 1;
  for(int i = 0; i < rank; i++) total *= dims[i];
  const hid_t space = H5Screate_simple(rank, dims, NULL);
  const hid_t dset = H5Dcreate2(obj, name, H5T_NATIVE_INT, space,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  int *const data = malloc((size_t)total * sizeof(*data));
  for(hsize_t i = 0; i < total; i++) data[i] = fill;
  H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  free(data);
  H5Dclose(dset);
  H5Sclose(space);
}

/* Minimal valid embedded NN model with bounded features and zero weights. */
static void coverage_create_degenerate_nn_file(const char *const path) {
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
  const hid_t file = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  const hid_t group = H5Gcreate2(file, "grhayl_nn_c2p", H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
  const hid_t dims = H5Gcreate2(group, "dims", H5P_DEFAULT,
                                H5P_DEFAULT, H5P_DEFAULT);
  const hid_t meta = H5Gcreate2(group, "meta", H5P_DEFAULT,
                                H5P_DEFAULT, H5P_DEFAULT);
  const hid_t scaling = H5Gcreate2(group, "scaling", H5P_DEFAULT,
                                   H5P_DEFAULT, H5P_DEFAULT);
  const hid_t layers = H5Gcreate2(group, "layers", H5P_DEFAULT,
                                  H5P_DEFAULT, H5P_DEFAULT);
  const int in_dim = 4, hidden_dim = 2, n_hidden = 1, out_dim = 1;
  const hsize_t h_in_dim[1] = { in_dim };
  const hsize_t h_out_dim[1] = { out_dim };
  const hsize_t h_hidden[1] = { hidden_dim };
  const hsize_t h_w_in[2] = { hidden_dim, in_dim };
  const hsize_t h_w_out[2] = { out_dim, hidden_dim };
  write_i32_scalar(dims, "in_dim", in_dim);
  write_i32_scalar(dims, "hidden_dim", hidden_dim);
  write_i32_scalar(dims, "n_hidden", n_hidden);
  write_i32_scalar(dims, "out_dim", out_dim);
  write_i32_scalar(meta, "q_idx", 0);
  write_i32_scalar(meta, "s_idx", 2);
  write_f32_scalar(scaling, "x_eps", 1e-6f);
  write_f32_scalar(meta, "y_eps", 1e-6f);
  write_f32_scalar(meta, "dx_eps", 1e-6f);
  write_i32_array(scaling, "x_kind", 1, h_in_dim, 0);
  /* A strictly increasing feature range satisfies model validation. */
  write_f32_array(scaling, "x_lo", 1, h_in_dim, -1.0f);
  write_f32_array(scaling, "x_hi", 1, h_in_dim, 1.0f);
  write_f32_array(scaling, "x_invrng", 1, h_in_dim, 1.0f);
  write_i32_array(scaling, "out_kind", 1, h_out_dim, 0);
  write_f32_array(scaling, "out_lo", 1, h_out_dim, 0.0f);
  write_f32_array(scaling, "out_hi", 1, h_out_dim, 1.0f);
  write_f32_array(scaling, "out_invrng", 1, h_out_dim, 1.0f);
  write_f32_array(layers, "W_in", 2, h_w_in, 0.0f);
  write_f32_array(layers, "b_in", 1, h_hidden, 0.0f);
  write_f32_array(layers, "W_out", 2, h_w_out, 0.0f);
  write_f32_array(layers, "b_out", 1, h_out_dim, 0.0f);
  H5Gclose(layers);
  H5Gclose(scaling);
  H5Gclose(meta);
  H5Gclose(dims);
  H5Gclose(group);
  H5Fclose(file);
}

static void test_nn_guess_x_lo_width_overflow_and_degenerate_ranges(void) {
  ghl_eos_parameters eos = { 0 };
  char model_path[] = "/tmp/grhayl-coverage-nn.XXXXXX";
  const int fd = mkstemp(model_path);
  CHECK(fd >= 0);
  if(fd < 0) return;
  CHECK(close(fd) == 0);
  coverage_create_degenerate_nn_file(model_path);
  const ghl_error_codes_t load_err =
        ghl_c2p_nn_load_from_eos_hdf5(model_path, &eos);
  CHECK(remove(model_path) == 0);
  CHECK(load_err == ghl_success);
  CHECK(eos.c2p_nn != NULL);
  if(load_err != ghl_success || eos.c2p_nn == NULL) {
    return;
  }

  /* q = FLT_MAX saturates the bracket midpoint: width stays at FLT_MAX and
   * x_lo + 0.5f*width overflows to a non-finite value, so the midpoint guard
   * clamps the guess to the zero sentinel. */
  const ghl_nn_c2p_input_t midpoint_inf_input = { FLT_MAX, 0.25f, 0.0f, 0.1f };
  const ghl_nn_c2p_guess_t midpoint_inf_guess =
        ghl_c2p_nn_guess(eos.c2p_nn, midpoint_inf_input);
  CHECK(isfinite(midpoint_inf_guess.x) && midpoint_inf_guess.x == 0.0f);

  /* Finite operands can still overflow x_lo = 1 + q - s. */
  const ghl_nn_c2p_input_t lower_bound_inf = { FLT_MAX, 0.25f, -FLT_MAX, 0.1f };
  CHECK(ghl_c2p_nn_guess(eos.c2p_nn, lower_bound_inf).x == 0.0f);

  /* q == s == 0 keeps x_lo/width finite while NaN r/t trigger the finite
   * fallback, confirming the model is usable once the guards pass. */
  const ghl_nn_c2p_input_t nan_rt_input = { 0.0f, NAN, 0.0f, NAN };
  const ghl_nn_c2p_guess_t nan_rt_guess =
        ghl_c2p_nn_guess(eos.c2p_nn, nan_rt_input);
  CHECK(isfinite(nan_rt_guess.x));

  ghl_c2p_nn_free(eos.c2p_nn);
  eos.c2p_nn = NULL;
}

/*
 * Drives ghl_con2prim_tabulated_multi_method's neural-network retry inside
 * the backup loop (con2prim_multi_method.c:215-217 and the speed_limited
 * fold at :220), which no reference-data suite reaches: with a failed main
 * routine, a failed first backup, and enable_neural_net_c2p set, the driver
 * must copy the NN-adjusted guess back into prims, record nn_guess_used and
 * backup[0], and retry the same backup key once more.
 *
 * The NN guess at line 179 runs through the production primitive-guess
 * helper, whose contract requires a valid tabulated EOS, so this test boots
 * the checked-in sample table with NN support disabled, then attaches the
 * same synthetic model the loader produces and flips the enable flag. Every
 * configured key is Font1D, which the tabulated dispatcher rejects with
 * ghl_error_invalid_c2p_key before any EOS table access; the retry decision
 * and bookkeeping are asserted exactly.
 */
static void test_tabulated_multi_method_nn_backup_retry(void) {
  static const char table_path[]
        = "Unit_Tests/sample_table/"
          "Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817_simple.h5";
  ghl_eos_parameters eos = { 0 };
  ghl_initialize_eos_functions(ghl_eos_tabulated);
  const ghl_error_codes_t table_err = ghl_initialize_tabulated_eos(
        table_path, ghl_eos_table_stellarcollapse,
        false /* clean_sound_speed */, false /* enable_neural_net_c2p */,
        1e-12 /* rho_atm */, 1e-12 /* rho_min */, 1e300 /* rho_max */,
        0.5 /* Y_e_atm */, 0.05 /* Y_e_min */, 0.5 /* Y_e_max */,
        1e-2 /* T_atm */, 1e-2 /* T_min */, 1e2 /* T_max */, &eos);
  CHECK(table_err == ghl_success);
  if(table_err != ghl_success) return;

  /* The sample table carries no NN weights, so enable_neural_net_c2p had to
   * stay off during initialization; attach the synthetic model directly and
   * enable the retry flag exactly as a model-bearing table would. */
  char model_path[] = "/tmp/grhayl-coverage-nn-retry.XXXXXX";
  const int fd = mkstemp(model_path);
  CHECK(fd >= 0 && close(fd) == 0);
  if(fd < 0) {
    ghl_tabulated_free_memory(&eos);
    return;
  }
  coverage_create_degenerate_nn_file(model_path);
  const ghl_error_codes_t load_err =
        ghl_c2p_nn_load_from_eos_hdf5(model_path, &eos);
  CHECK(remove(model_path) == 0);
  CHECK(load_err == ghl_success);
  CHECK(eos.c2p_nn != NULL);
  if(load_err != ghl_success || eos.c2p_nn == NULL) {
    ghl_tabulated_free_memory(&eos);
    return;
  }
  eos.enable_neural_net_c2p = true;

  ghl_parameters params;
  ghl_metric_quantities metric;
  ghl_ADM_aux_quantities metric_aux;
  ghl_conservative_quantities cons;
  ghl_primitive_quantities prims;
  ghl_con2prim_diagnostics diag;

  /* main = backup[0] = Font1D: a hybrid-only key the tabulated dispatcher
   * rejects with ghl_error_invalid_c2p_key before any EOS access. Both
   * attempts fail, driving the NN retry inside the backup loop; the remaining
   * backup slots stay None so the loop stops after the first. */
  const ghl_con2prim_id_t backups[3]
        = { ghl_con2prim_id_Font1D, ghl_con2prim_id_None, ghl_con2prim_id_None };
  ghl_initialize_params(ghl_con2prim_id_Font1D, backups,
                        false, false, false, 1e8, 60.0, 0.0, &params);
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0,
                        &metric);
  memset(&metric_aux, 0, sizeof(metric_aux));
  ghl_initialize_conservatives(1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, &cons);
  ghl_initialize_primitives(1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, &prims);
  ghl_initialize_diagnostics(&diag);

  const ghl_error_codes_t err = ghl_con2prim_tabulated_multi_method(
        &params, &eos, &metric, &metric_aux, &cons, &prims, &diag);

  /* Dispatcher rejects every configured key: final error is the last backup's
   * ghl_error_invalid_c2p_key, with the retry bookkeeping recorded. */
  CHECK(err == ghl_error_invalid_c2p_key);
  CHECK(diag.nn_guess_used == true);
  CHECK(diag.backup[0] == true);
  CHECK(diag.backup[1] == false);
  CHECK(diag.backup[2] == false);

  ghl_c2p_nn_free(eos.c2p_nn);
  eos.c2p_nn = NULL;
  ghl_tabulated_free_memory(&eos);
}
#endif

int main(void) {
  g_fail = 0;
  test_brent_endpoint_roots_and_tiny_interval();
  test_hybrid_multi_method_backup_none_break();
#ifndef GHL_DISABLE_HDF5
  test_nn_guess_x_lo_width_overflow_and_degenerate_ranges();
  test_tabulated_multi_method_nn_backup_retry();
#endif
  if(g_fail) return 1;
  puts("unit_test_coverage_con2prim: all passed");
  return 0;
}
