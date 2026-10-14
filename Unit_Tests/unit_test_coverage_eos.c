/* Targeted coverage tests for EOS/Hybrid and EOS/Tabulated production paths
 * not exercised by the existing suite. Fully self-contained: HDF5 fixtures
 * are generated in C (the same libhdf5 API the production loader uses) in a
 * private per-PID temporary directory and removed before main() returns.
 * No Python or network dependency.
 *
 * Allocation-fault injection (Linux only) uses an RTLD_NEXT malloc/calloc
 * interposer.  Every request is first filtered by caller ownership: only
 * allocations whose immediate caller resolves (dladdr) to libghl are
 * counted or failed, so libhdf5 traffic is physically excluded even when
 * the armed window spans a read that calls it.  The size key then selects
 * the target allocation deterministically.  Dry runs count matches through
 * a pipe; injected runs fork() and child exit(0) so libc atexit handlers
 * (__gcov_exit) flush coverage counters.  The child detects leftover
 * pointers after failure and reports them as leaks.
 *
 */
#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif
#if defined(__APPLE__) && !defined(_DARWIN_C_SOURCE)
#define _DARWIN_C_SOURCE
#endif

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "ghl_unit_tests.h"

static void expect_code(
      const ghl_error_codes_t got,
      const ghl_error_codes_t want,
      const int key,
      const char *what) {
  if(got != want) {
    fprintf(
          stderr, "Test %d failed: %s returned %d, expected %d\n", key, what, got, want);
    exit(1);
  }
}
static void expect_true(const bool cond, const int key, const char *msg) {
  if(!cond) {
    fprintf(stderr, "Test %d failed: %s\n", key, msg);
    exit(1);
  }
}

typedef void (*child_work_fn)(void *);

/* Reap pid; waitpid is only retried on EINTR, anything else is fatal. */
static int reap_child(const int key, const pid_t pid) {
  int st = 0;
  for(;;) {
    const pid_t r = waitpid(pid, &st, 0);
    if(r == pid) {
      break;
    }
    if(r < 0 && errno == EINTR) {
      continue;
    }
    fprintf(
          stderr, "Test %d failed: waitpid(%ld): %s\n", key, (long)pid, strerror(errno));
    exit(1);
  }
  if(WIFSIGNALED(st)) {
    fprintf(stderr, "Test %d failed: child killed by signal %d\n", key, WTERMSIG(st));
    exit(1);
  }
  return WEXITSTATUS(st);
}

static void
run_in_child(const int key, child_work_fn work, void *arg, const int expected_exit) {
  fflush(NULL);
  const pid_t pid = fork();
  if(pid < 0) {
    fprintf(stderr, "Test %d failed: fork(): %s\n", key, strerror(errno));
    exit(1);
  }
  if(pid == 0) {
    work(arg);
    fflush(NULL);
    exit(0); /* libc atexit handlers (gcov counters) run */
  }
  const int got = reap_child(key, pid);
  if(got != expected_exit) {
    fprintf(
          stderr, "Test %d failed: child exit=%d expected=%d\n", key, got,
          expected_exit);
    exit(1);
  }
}

/* ------------------------------------------------------------------ */
/* EOS/Hybrid                                                          */
/* ------------------------------------------------------------------ */
static void test_hybrid_find_index_from_P(const int key0) {
  const int neos = 4;
  const double rho_b_min = 1e-12, rho_b_max = 1e300;
  const double rho_ppoly[3] = { 2.44034e+07, 3.78358e+11, 2.62780e+12 };
  const double Gamma_ppoly[4] = { 1.58425, 1.28733, 0.62223, 1.35692 };
  const double k_ppoly0 = 6.80110e-9, Gamma_th = 2.0;

  ghl_eos_parameters eos = { 0 };
  expect_code(
        ghl_initialize_hybrid_eos_functions_and_params(
              rho_b_min, rho_b_min, rho_b_max, neos, rho_ppoly, Gamma_ppoly, k_ppoly0,
              Gamma_th, &eos),
        ghl_success, key0, "ghl_initialize_hybrid_eos_functions_and_params");

  ghl_eos_parameters one = { 0 };
  one.neos = 1;
  expect_true(
        NRPyEOS_find_polytropic_index_from_P(&one, 1.0) == 0, key0,
        "neos==1 must force index 0");

  expect_true(
        NRPyEOS_find_polytropic_index_from_P(&eos, 0.0) == 0, key0,
        "P=0 must map to polytrope 0");
  expect_true(
        NRPyEOS_find_polytropic_index_from_P(&eos, eos.p_ppoly[0]) == 1, key0,
        "P==p_ppoly[0] must advance to polytrope 1");
  expect_true(
        NRPyEOS_find_polytropic_index_from_P(&eos, eos.p_ppoly[1]) == 2, key0,
        "P==p_ppoly[1] must advance to polytrope 2");
  expect_true(
        NRPyEOS_find_polytropic_index_from_P(&eos, eos.p_ppoly[2]) == 3, key0,
        "P==p_ppoly[2] must advance to polytrope 3");
  expect_true(
        NRPyEOS_find_polytropic_index_from_P(&eos, 1e300) == 3, key0,
        "P above the last transition must clamp to the last index");
}

static void test_hybrid_rho_cold_from_P_cold(const int key0) {
  const int neos = 2;
  const double rho_b_min = 1e-12, rho_b_max = 1e300;
  const double rho_ppoly[1] = { 1.234 };
  const double Gamma_ppoly[2] = { 2.0, 3.0 };
  const double k_ppoly0 = 1.0, Gamma_th = 2.0;
  ghl_eos_parameters eos = { 0 };
  expect_code(
        ghl_initialize_hybrid_eos_functions_and_params(
              rho_b_min, rho_b_min, rho_b_max, neos, rho_ppoly, Gamma_ppoly, k_ppoly0,
              Gamma_th, &eos),
        ghl_success, key0, "ghl_initialize_hybrid_eos_functions_and_params");

  const double r1 = NRPyEOS_hybrid_compute_rho_cold_from_P_cold(&eos, k_ppoly0);
  expect_true(isfinite(r1) && r1 > 0.0, key0, "rho must be finite and positive");
  expect_true(
        relative_error(r1, 1.0) < 1e-13, key0,
        "rho_cold(P=K_ppoly) must be exactly 1.0");

  const double r2 = NRPyEOS_hybrid_compute_rho_cold_from_P_cold(&eos, eos.p_ppoly[0]);
  expect_true(
        relative_error(r2, rho_ppoly[0]) < 1e-13, key0,
        "rho_cold at p_ppoly[0] must reproduce rho_ppoly[0]");

  const double rf = NRPyEOS_hybrid_compute_rho_cold_from_P_cold(&eos, 1e-300);
  expect_true(rf >= eos.rho_min, key0, "rho_cold must respect the rho_min floor");
}

/* ------------------------------------------------------------------ */
/* EOS/Tabulated                                                       */
/* ------------------------------------------------------------------ */
#ifndef GHL_DISABLE_HDF5

#include "../GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.h"
#include "../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.h"
#include <hdf5.h>

static char fixture_dir[PATH_MAX];
static pid_t fixture_owner;
static void fixture_teardown(void);

/* Every HDF5 handle and write status is checked: silently truncated or
 * unwritten synthetic fixtures would produce misleading read-path results. */
static void w_scalar_int(hid_t f, const char *name, int v, const int key) {
  hsize_t d[1] = { 1 };
  const hid_t s = H5Screate_simple(1, d, NULL);
  const hid_t g = (s >= 0) ? H5Dcreate(
                                   f, name, H5T_NATIVE_INT, s, H5P_DEFAULT, H5P_DEFAULT,
                                   H5P_DEFAULT)
                           : -1;
  herr_t wr = -1;
  if(g >= 0) {
    wr = H5Dwrite(g, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
  }
  if(g >= 0) {
    expect_true(H5Dclose(g) >= 0, key, "H5Dclose(int) failed");
  }
  if(s >= 0) {
    expect_true(H5Sclose(s) >= 0, key, "H5Sclose(int) failed");
  }
  expect_true(wr >= 0, key, "H5Dwrite(int) failed");
}
static void w_scalar_dbl(hid_t f, const char *name, double v, const int key) {
  hsize_t d[1] = { 1 };
  const hid_t s = H5Screate_simple(1, d, NULL);
  const hid_t g = (s >= 0) ? H5Dcreate(
                                   f, name, H5T_NATIVE_DOUBLE, s, H5P_DEFAULT,
                                   H5P_DEFAULT, H5P_DEFAULT)
                           : -1;
  herr_t wr = -1;
  if(g >= 0) {
    wr = H5Dwrite(g, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
  }
  if(g >= 0) {
    expect_true(H5Dclose(g) >= 0, key, "H5Dclose(scalar) failed");
  }
  if(s >= 0) {
    expect_true(H5Sclose(s) >= 0, key, "H5Sclose(scalar) failed");
  }
  expect_true(wr >= 0, key, "H5Dwrite(scalar) failed");
}
static void w_1d(hid_t f, const char *name, const double *v, hsize_t n, const int key) {
  const hid_t s = H5Screate_simple(1, &n, NULL);
  const hid_t g = (s >= 0) ? H5Dcreate(
                                   f, name, H5T_NATIVE_DOUBLE, s, H5P_DEFAULT,
                                   H5P_DEFAULT, H5P_DEFAULT)
                           : -1;
  herr_t wr = -1;
  if(g >= 0) {
    wr = H5Dwrite(g, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v);
  }
  if(g >= 0) {
    expect_true(H5Dclose(g) >= 0, key, "H5Dclose(1d) failed");
  }
  if(s >= 0) {
    expect_true(H5Sclose(s) >= 0, key, "H5Sclose(1d) failed");
  }
  expect_true(wr >= 0, key, "H5Dwrite(1d) failed");
}
static void
w_3d(hid_t f,
     const char *name,
     const double *v,
     hsize_t ny,
     hsize_t nt,
     hsize_t nr,
     const int key) {
  hsize_t d[3] = { ny, nt, nr };
  const hid_t s = H5Screate_simple(3, d, NULL);
  const hid_t g = (s >= 0) ? H5Dcreate(
                                   f, name, H5T_NATIVE_DOUBLE, s, H5P_DEFAULT,
                                   H5P_DEFAULT, H5P_DEFAULT)
                           : -1;
  herr_t wr = -1;
  if(g >= 0) {
    wr = H5Dwrite(g, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v);
  }
  if(g >= 0) {
    expect_true(H5Dclose(g) >= 0, key, "H5Dclose(3d) failed");
  }
  if(s >= 0) {
    expect_true(H5Sclose(s) >= 0, key, "H5Sclose(3d) failed");
  }
  expect_true(wr >= 0, key, "H5Dwrite(3d) failed");
}

enum fx_kind {
  FX_VALID,
  FX_BADSIZE,
  FX_BADAXIS,
  FX_BADH,
  FX_CS2_NAN,
  FX_CS2_SUPER,
  FX_CS2_NEG
};

static void write_table(const char *path, const enum fx_kind kind, const int key) {
  const int n = 4, npts = n * n * n;
  double logrho[4], logtemp[4], ye[4], data[19][64];

  const double l10 = log(10.0);
  for(int i = 0; i < n; i++) {
    const double t = 1.0 + i;
    logrho[i] = (t - log(CGS_TO_CODE_DENSITY)) / l10;
    logtemp[i] = t / l10;
    ye[i] = 0.1 + 0.1 * i;
  }

  const double conv = CGS_TO_CODE_LENGTH * CGS_TO_CODE_LENGTH
                      / (CGS_TO_CODE_TIME * CGS_TO_CODE_TIME);
  for(int i = 0; i < npts; i++) {
    const double rho = 1.0 + (i % n);
    const double T = 1.0 + ((i / n) % n);
    const double Ye = 0.1 + 0.1 * ((i / (n * n)) % n);
    data[0][i] = (rho + Ye + T - log(CGS_TO_CODE_PRESSURE)) / l10;
    data[1][i] = (rho + Ye - T - log(CGS_TO_CODE_ENERGY)) / l10;
    data[2][i] = rho - Ye + T;
    data[3][i] = rho - Ye - T;
    data[4][i] = (-rho + Ye + T) / conv;
    data[5][i] = (-rho + Ye - T) / CGS_TO_CODE_ENERGY;
    data[6][i] = (-rho - Ye + T) / (CGS_TO_CODE_PRESSURE / CGS_TO_CODE_DENSITY);
    data[7][i] = (-rho - Ye - T) / (CGS_TO_CODE_PRESSURE / CGS_TO_CODE_ENERGY);
    data[8][i] = 2 * rho + Ye + T;
    data[9][i] = 2 * rho + Ye - T;
    data[10][i] = 2 * rho - Ye + T;
    data[11][i] = 2 * rho - Ye - T;
    data[12][i] = -2 * rho + Ye + T;
    data[13][i] = -2 * rho + Ye - T;
    data[14][i] = -2 * rho - Ye + T;
    data[15][i] = -2 * rho - Ye - T;
    data[16][i] = rho + 3 * Ye - T;
    data[17][i] = rho - 3 * Ye + T;
    data[18][i] = rho - 3 * Ye - T;
  }

  if(kind == FX_CS2_NAN) {
    for(int i = 0; i < npts; i++) {
      data[4][i] = NAN;
    }
  }
  if(kind == FX_CS2_SUPER) {
    for(int i = 0; i < npts; i++) {
      data[4][i] = 1.5 / conv;
    }
  }
  if(kind == FX_CS2_NEG) {
    for(int i = 0; i < npts; i++) {
      data[4][i] = -0.5 / conv;
    }
  }
  if(kind == FX_BADH) {
    for(int i = 0; i < npts; i++) {
      data[0][i] = NAN;
      data[1][i] = NAN;
    }
  }

  const int lw = (kind == FX_BADSIZE) ? 3 : 4;
  if(kind == FX_BADAXIS) {
    logrho[2] = logrho[1];
  }

  const hid_t f = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  expect_true(f >= 0, key, "H5Fcreate failed");
  w_scalar_int(f, "have_rel_cs2", 1, key);
  w_scalar_int(f, "pointsrho", lw, key);
  w_scalar_int(f, "pointstemp", n, key);
  w_scalar_int(f, "pointsye", n, key);
  w_scalar_dbl(f, "energy_shift", 123e-45 / CGS_TO_CODE_ENERGY, key);
  w_1d(f, "logrho", logrho, (hsize_t)lw, key);
  w_1d(f, "logtemp", logtemp, (hsize_t)n, key);
  w_1d(f, "ye", ye, (hsize_t)n, key);
  static const char *q[19]
        = { "logpress", "logenergy", "entropy", "munu", "cs2",  "dedt", "dpdrhoe",
            "dpderho",  "muhat",     "mu_e",    "mu_p", "mu_n", "Xa",   "Xh",
            "Xn",       "Xp",        "Abar",    "Zbar", "gamma" };
  for(int k = 0; k < 19; k++) {
    w_3d(f, q[k], data[k], n, n, n, key);
  }
  expect_true(H5Fclose(f) >= 0, key, "H5Fclose failed");
}

static void fixture_setup(void) {
  char tmpl[PATH_MAX];
  snprintf(tmpl, sizeof(tmpl), "/tmp/grhayl_cov_eos_XXXXXX");
  expect_true(mkdtemp(tmpl) != NULL, 60, "mkdtemp failed");
  fixture_owner = getpid();
  strncpy(fixture_dir, tmpl, sizeof(fixture_dir) - 1);
  fixture_dir[sizeof(fixture_dir) - 1] = '\0';
  expect_true(atexit(fixture_teardown) == 0, 60, "register fixture cleanup");
  char p[PATH_MAX];
#define W(k, n)                                         \
  do {                                                  \
    snprintf(p, sizeof(p), "%s/%s.h5", fixture_dir, n); \
    write_table(p, FX_##k, 60);                         \
  } while(0)
  W(VALID, "valid");
  W(BADSIZE, "badsize");
  W(BADAXIS, "badaxis");
  W(BADH, "badh");
  W(CS2_NAN, "cs2_nan");
  W(CS2_SUPER, "cs2_super");
  W(CS2_NEG, "cs2_neg");
#undef W
}

static void fixture_teardown(void) {
  if(fixture_dir[0] == '\0' || getpid() != fixture_owner) {
    return;
  }
  static const char *files[]
        = { "valid", "badsize", "badaxis", "badh", "cs2_nan", "cs2_super", "cs2_neg" };
  char p[PATH_MAX];
  for(size_t i = 0; i < sizeof(files) / sizeof(files[0]); i++) {
    snprintf(p, sizeof(p), "%s/%s.h5", fixture_dir, files[i]);
    unlink(p);
  }
  rmdir(fixture_dir);
  fixture_dir[0] = '\0';
}

static ghl_eos_parameters make_eos(const bool clean) {
  ghl_eos_parameters eos = { 0 };
  eos.eos_type = ghl_eos_tabulated;
  eos.table_type = ghl_eos_table_stellarcollapse;
  eos.clean_sound_speed = clean;
  return eos;
}
static void path_of(char *out, const size_t n, const char *name) {
  snprintf(out, n, "%s/%s.h5", fixture_dir, name);
}

static void test_read_table_then_adjust_cs2(const int key0) {
  char p[PATH_MAX];

  path_of(p, sizeof(p), "cs2_nan");
  ghl_eos_parameters eos = make_eos(true);
  expect_code(
        NRPyEOS_read_table_set_EOS_params(p, &eos), ghl_success, key0,
        "read valid NaN cs2 table");
  for(int i = 0; i < 64; i++) {
    expect_true(
          eos.table_all[NRPyEOS_cs2_key + NRPyEOS_ntablekeys * i] == 1.0 - 1e-15, key0,
          "NaN cs2 must be cleaned to 1.0-1e-15");
  }
  NRPyEOS_free_memory(&eos);

  path_of(p, sizeof(p), "cs2_super");
  eos = make_eos(true);
  expect_code(
        NRPyEOS_read_table_set_EOS_params(p, &eos), ghl_success, key0,
        "read valid superluminal cs2 table");
  for(int i = 0; i < 64; i++) {
    expect_true(
          eos.table_all[NRPyEOS_cs2_key + NRPyEOS_ntablekeys * i] == 1.0, key0,
          "superluminal cs2 must be cleaned to exactly 1.0");
  }
  NRPyEOS_free_memory(&eos);

  path_of(p, sizeof(p), "cs2_neg");
  eos = make_eos(true);
  expect_code(
        NRPyEOS_read_table_set_EOS_params(p, &eos), ghl_success, key0,
        "read valid negative cs2 table");
  for(int i = 0; i < 64; i++) {
    expect_true(
          eos.table_all[NRPyEOS_cs2_key + NRPyEOS_ntablekeys * i] == 0.0, key0,
          "negative cs2 must be cleaned to exactly 0.0");
  }
  NRPyEOS_free_memory(&eos);
}

static void test_read_table_error_paths(const int key0) {
  char p[PATH_MAX];
  path_of(p, sizeof(p), "valid");

  expect_code(
        NRPyEOS_read_table_set_EOS_params(p, NULL), ghl_error_eos_struct_is_null, key0,
        "NULL eos struct must report eos_struct_is_null");

  ghl_eos_parameters wrong_family = { 0 };
  wrong_family.eos_type = ghl_eos_hybrid;
  wrong_family.table_type = ghl_eos_table_stellarcollapse;
  expect_code(
        NRPyEOS_read_table_set_EOS_params(p, &wrong_family), ghl_error_invalid_eos_type,
        key0, "hybrid eos_type must be rejected");

  ghl_eos_parameters wrong_type = { 0 };
  wrong_type.eos_type = ghl_eos_tabulated;
  wrong_type.table_type = ghl_eos_table_types;
  expect_code(
        NRPyEOS_read_table_set_EOS_params(p, &wrong_type),
        ghl_error_invalid_eos_table_type, key0, "sentinel table_type must be rejected");

  path_of(p, sizeof(p), "badsize");
  ghl_eos_parameters eos = make_eos(true);
  const ghl_error_codes_t badsize_rc = NRPyEOS_read_table_set_EOS_params(p, &eos);
  expect_true(
        badsize_rc == ghl_error_hdf5_dataset_size_mismatch, key0,
        "short logrho axis must be rejected as size mismatch");
  if(eos.table_logrho != NULL) {
    NRPyEOS_free_memory(&eos);
  }

  path_of(p, sizeof(p), "badaxis");
  eos = make_eos(true);
  expect_code(
        NRPyEOS_read_table_set_EOS_params(p, &eos), ghl_error_invalid_eos_table, key0,
        "non-monotone logrho axis must be rejected");
  if(eos.table_logrho != NULL) {
    NRPyEOS_free_memory(&eos);
  }

  path_of(p, sizeof(p), "badh");
  eos = make_eos(true);
  expect_code(
        NRPyEOS_read_table_set_EOS_params(p, &eos), ghl_error_invalid_eos_table, key0,
        "NaN enthalpy table must be rejected");
  if(eos.table_logrho != NULL) {
    NRPyEOS_free_memory(&eos);
  }
}

static void test_tabulated_beta_equilibrium_paths(const int key0) {
  char p[PATH_MAX];
  path_of(p, sizeof(p), "valid");
  ghl_eos_parameters eos = make_eos(true);
  expect_code(
        NRPyEOS_read_table_set_EOS_params(p, &eos), ghl_success, key0,
        "read valid table for beq tests");
  /* read_table fills the table but not the tabulated ghl_* function
   * pointers; those are installed by NRPyEOS_initialize_tabulated_functions()
   * (also invoked by the public ghl_initialize_tabulated_eos_functions_and_params). */
  NRPyEOS_initialize_tabulated_functions();

  expect_code(
        NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T(eos.table_T_min * 0.5, &eos),
        ghl_error_table_min_T, key0, "T below the table must report table_min_T");
  expect_code(
        NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T(eos.table_T_max * 2.0, &eos),
        ghl_error_table_max_T, key0, "T above the table must report table_max_T");

  expect_code(
        NRPyEOS_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(
              0.5 * (eos.table_T_min + eos.table_T_max), &eos),
        ghl_success, key0, "mid-table beta-equilibrium call must succeed");
  expect_true(eos.Ye_of_lr != NULL, key0, "Ye_of_lr must be allocated");
  expect_true(eos.lp_of_lr != NULL, key0, "lp_of_lr must be allocated");

  /* log(rho) outside the beta-equilibrium table must be rejected ... */
  double P_out;
  expect_code(
        ghl_tabulated_compute_P_from_rho(&eos, eos.table_rho_max * 2.0, &P_out),
        ghl_error_root_not_bracketed, key0,
        "rho above the beq table must report root_not_bracketed");

  /* ... while log(rho) == log(table_rho_max) (the exact last grid point)
   * must clamp the left index to nx-2 and succeed. */
  double P_hi;
  expect_code(
        ghl_tabulated_compute_P_from_rho(&eos, eos.table_rho_max, &P_hi), ghl_success,
        key0, "log(rho)==last grid point must clamp to the last interval");
  expect_true(isfinite(P_hi), key0, "upper-edge P must be finite");

  const double rho_mid = 0.5 * (eos.table_rho_min + eos.table_rho_max);
  const double ye_mid = 0.5 * (eos.table_Y_e_min + eos.table_Y_e_max);
  double eps_min_at_rho, P_tmp;
  expect_code(
        ghl_tabulated_compute_P_eps_from_T(
              &eos, rho_mid, ye_mid, eos.table_T_min, &P_tmp, &eps_min_at_rho),
        ghl_success, key0, "eps at T_min probe");
  double T_floor;
  expect_code(
        ghl_tabulated_compute_T_from_eps(
              &eos, rho_mid, ye_mid, eps_min_at_rho, &T_floor),
        ghl_success, key0, "T_from_eps at the lower edge must find T_min");
  expect_true(
        T_floor >= eos.table_T_min, key0,
        "T_floor must not fall below the table minimum");

  NRPyEOS_free_memory(&eos);
}

static void test_interpolation_failures_and_enthalpy(const int key) {
  char p[PATH_MAX];
  path_of(p, sizeof(p), "valid");
  ghl_eos_parameters eos = make_eos(false);
  expect_code(
        NRPyEOS_read_table_set_EOS_params(p, &eos), ghl_success, key,
        "load interpolation fixture");
  eos.root_finding_precision = 1e-12;
  const double rho = sqrt(eos.table_rho_min * eos.table_rho_max);
  const double ye = 0.5 * (eos.table_Y_e_min + eos.table_Y_e_max);
  double a = 17.0, b = 19.0, c = 23.0, T = eos.table_T_min;
  expect_code(
        NRPyEOS_cs2_from_rho_Ye_T(&eos, eos.table_rho_max * 2, ye, T, &a),
        ghl_error_table_max_rho, key, "cs2 bounds");
  expect_true(a == 17.0, key, "cs2 failure preserves output");
  expect_code(
        NRPyEOS_eps_and_T_from_rho_Ye_P(&eos, eos.table_rho_max * 2, ye, 1.0, &a, &T),
        ghl_error_table_max_rho, key, "P inversion bounds");
  expect_true(a == 17.0 && T == eos.table_T_min, key, "P inversion preserves output");
  expect_code(
        NRPyEOS_P_eps_S_and_T_from_rho_Ye_h(
              &eos, eos.table_rho_max * 2, ye, 1.0, &a, &b, &c, &T),
        ghl_error_table_max_rho, key, "h inversion bounds");
  expect_true(
        a == 17.0 && b == 19.0 && c == 23.0 && T == eos.table_T_min, key,
        "h inversion preserves outputs");
  const int keys[] = { NRPyEOS_press_key };
  expect_code(
        NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities(
              &eos, 1, 1e-12, rho, ye, -eos.energy_shift - 1.0, NRPyEOS_eps_key, keys,
              &a, &T),
        ghl_error_table_neg_energy, key, "negative shifted energy");
  expect_true(
        a == 17.0 && T == eos.table_T_min, key, "energy failure preserves outputs");
  const double bad_T[] = { NAN, INFINITY, 0.0, -1.0 };
  for(size_t i = 0; i < sizeof(bad_T) / sizeof(bad_T[0]); i++) {
    expect_true(
          NRPyEOS_tabulated_get_index_T(&eos, bad_T[i]) == -1, key,
          "invalid temperature index");
  }

  NRPyEOS_stellarcollapse_t sc = { 0 };
  expect_code(
        NRPyEOS_stellarcollapse_to_ghl(&sc, &eos), ghl_error_invalid_eos_table, key,
        "invalid conversion dimensions");
  NRPyEOS_stellarcollapse_free_table(NULL);
  NRPyEOS_free_memory(&eos);

  double logrho = 0.0, logh = 31.0;
  double all[NRPyEOS_ntablekeys] = { 0 };
  ghl_eos_parameters one = { .N_rho = 1,
                             .N_T = 1,
                             .N_Ye = 1,
                             .table_logrho = &logrho,
                             .table_logh = &logh,
                             .table_all = all,
                             .energy_shift = 4.0 };
  all[NRPyEOS_enthalpy_key] = 29.0;
  expect_code(
        NRPyEOS_tabulate_enthalpy_checked(&one), ghl_error_invalid_eos_table, key,
        "negative enthalpy rejected");
  expect_true(
        logh == 31.0 && all[NRPyEOS_enthalpy_key] == 29.0, key,
        "invalid enthalpy preserves outputs");
  one.energy_shift = 0.0;
  NRPyEOS_tabulate_enthalpy(&one);
  expect_true(
        fabs(logh - log(3.0)) < 1e-14 && all[NRPyEOS_enthalpy_key] == logh, key,
        "legacy enthalpy wrapper stores log h");
}

/* ------------------------------------------------------------------ */
/* Linux-only OOM injection via RTLD_NEXT interposer                   */
/* ------------------------------------------------------------------ */
#if defined(__linux__)

#include <dlfcn.h>

/* Caller-qualified allocation-fault injection.
 *
 * The armed window may span code that internally calls libhdf5 (a table
 * read), so windowing alone cannot exclude libhdf5 traffic: the interposer
 * qualifies the immediate caller via dladdr and only counts or fails
 * requests whose calling address is owned by libghl.  The dry run then
 * measures the exact GRHayL occurrence count for the (size,kind) key and
 * the injected run fails one deterministic occurrence, checking the error
 * code and that every to_ghl outcome pointer is NULL afterwards. */

enum { RK_MALLOC = 1, RK_CALLOC = 2 };

static void *(*ghl_real_malloc)(size_t) = NULL;
static void *(*ghl_real_calloc)(size_t, size_t) = NULL;
static void (*ghl_real_free)(void *) = NULL;

/* Record actual ownership, rather than treating nulled outputs as proof that
 * the allocations behind them were freed. Tracker nodes use the real allocator
 * and therefore cannot themselves trigger an injected GRHayL failure. */
struct tracked_allocation {
  void *pointer;
  struct tracked_allocation *next;
};
static struct tracked_allocation *tracked_allocations;

static void track_allocation(void *pointer, void *caller_address);

static long ghl_oom_arm_size = 0;         /* byte size to match while armed */
static long ghl_oom_arm_kind = 0;         /* RK_MALLOC / RK_CALLOC          */
static long ghl_oom_fail_occurrence = -1; /* >=0: fail this occurrence      */
static long ghl_oom_seen = 0;             /* matching requests so far       */
static long ghl_oom_fired = 0;            /* nonzero iff a NULL was forged  */
static int ghl_oom_armed = 0;             /* nonzero iff the size key is set */
static __thread int ghl_oom_in_hook = 0;  /* loader-allocation recursion guard */

/* Caller-ownership seam.
 *
 * libhdf5 floods every plausible size bucket, so (size,window) keys cannot
 * isolate GRHayL.  The deterministic seam tested at the head of every
 * interposed malloc/calloc is the *immediate caller*: only requests whose
 * calling instruction address resolves (dladdr) to this process's copy of
 * libghl are eligible for counting or failure.  That is recursion- and
 * concurrency-safe because it needs no state beyond the armed size key and
 * dladdr is re-entrant.  Test-side callers live in the test executable and
 * never satisfy the filter.
 *
 * The thread-local guard sends any allocator recursion from dlsym or dladdr
 * directly to the real allocator. */

static int ghl_caller_is_libghl(void *address) {
  Dl_info info;
  if(address == NULL) {
    return 0;
  }
  if(dladdr(address, &info) == 0 || info.dli_fname == NULL) {
    return 0;
  }
  /* The shared library is linked as libghl.so -> libghl_<version>.so.
   * A static link resolves to the test executable instead, which never
   * matches these names. */
  const char *base = strrchr(info.dli_fname, '/');
  base = (base != NULL) ? base + 1 : info.dli_fname;
  return strncmp(base, "libghl", 6) == 0;
}

enum h5_fault {
  H5_FAULT_NONE,
  H5_FAULT_SPACE,
  H5_FAULT_NDIMS,
  H5_FAULT_DIMS,
  H5_FAULT_READ
};
static enum h5_fault h5_fault;
static int h5_fault_fired;

static bool fail_h5(const enum h5_fault target, void *caller) {
  if(h5_fault != target || !ghl_caller_is_libghl(caller)) {
    return false;
  }
  h5_fault_fired++;
  return true;
}

hid_t H5Dget_space(hid_t dataset) {
  typedef hid_t (*fn)(hid_t);
  fn real = (fn)dlsym(RTLD_NEXT, "H5Dget_space");
  if(real == NULL) {
    exit(42);
  }
  if(fail_h5(H5_FAULT_SPACE, __builtin_return_address(0))) {
    return -1;
  }
  return real(dataset);
}

int H5Sget_simple_extent_ndims(hid_t space) {
  typedef int (*fn)(hid_t);
  fn real = (fn)dlsym(RTLD_NEXT, "H5Sget_simple_extent_ndims");
  if(real == NULL) {
    exit(42);
  }
  if(fail_h5(H5_FAULT_NDIMS, __builtin_return_address(0))) {
    return -1;
  }
  return real(space);
}

int H5Sget_simple_extent_dims(hid_t space, hsize_t dims[], hsize_t maxdims[]) {
  typedef int (*fn)(hid_t, hsize_t *, hsize_t *);
  fn real = (fn)dlsym(RTLD_NEXT, "H5Sget_simple_extent_dims");
  if(real == NULL) {
    exit(42);
  }
  if(fail_h5(H5_FAULT_DIMS, __builtin_return_address(0))) {
    return -1;
  }
  return real(space, dims, maxdims);
}

herr_t H5Dread(
      hid_t dataset,
      hid_t dtype,
      hid_t memory_space,
      hid_t file_space,
      hid_t transfer,
      void *buffer) {
  typedef herr_t (*fn)(hid_t, hid_t, hid_t, hid_t, hid_t, void *);
  fn real = (fn)dlsym(RTLD_NEXT, "H5Dread");
  if(real == NULL) {
    exit(42);
  }
  if(fail_h5(H5_FAULT_READ, __builtin_return_address(0))) {
    return -1;
  }
  return real(dataset, dtype, memory_space, file_space, transfer, buffer);
}

static void test_hdf5_faults(const int key) {
  char p[PATH_MAX];
  path_of(p, sizeof(p), "valid");
  const hid_t file = H5Fopen(p, H5F_ACC_RDONLY, H5P_DEFAULT);
  expect_true(file >= 0, key, "open helper fixture");
  const ssize_t handles = H5Fget_obj_count(file, H5F_OBJ_ALL);
  expect_true(handles >= 0, key, "count fixture handles");
  double sentinel = 17.0;
  for(enum h5_fault fault = H5_FAULT_SPACE; fault <= H5_FAULT_READ; fault++) {
    void *output = &sentinel;
    h5_fault = fault;
    h5_fault_fired = 0;
    const ghl_error_codes_t rc
          = NRPyEOS_hdf5_read_double_dataset(file, "logrho", 4, &output);
    h5_fault = H5_FAULT_NONE;
    expect_code(
          rc,
          fault == H5_FAULT_READ ? ghl_error_hdf5_dataset_could_not_read
                                 : ghl_error_hdf5_dataset_invalid_ndims,
          key, "injected HDF5 failure");
    expect_true(h5_fault_fired == 1, key, "intended HDF5 fault must fire once");
    expect_true(
          output == &sentinel && sentinel == 17.0, key,
          "HDF5 failure must preserve output");
    expect_true(
          H5Fget_obj_count(file, H5F_OBJ_ALL) == handles, key,
          "HDF5 failure must close dataset handles");
  }
  expect_true(H5Fclose(file) >= 0, key, "close helper fixture");
}

static int oom_match(const size_t sz, const int kind, void *caller_address) {
  if(!ghl_oom_armed) {
    return 0;
  }
  if((long)sz != ghl_oom_arm_size) {
    return 0;
  }
  if(kind != ghl_oom_arm_kind) {
    return 0;
  }
  if(sz == 0) {
    return 0;
  }
  if(!ghl_caller_is_libghl(caller_address)) {
    return 0;
  }
  ghl_oom_seen++;
  if(ghl_oom_fail_occurrence >= 0 && ghl_oom_seen == ghl_oom_fail_occurrence + 1) {
    ghl_oom_fired++;
    return 1;
  }
  return 0;
}

void *malloc(size_t size) {
  if(ghl_real_malloc == NULL) {
    typedef void *(*real_fn)(size_t);
    ghl_real_malloc = (real_fn)dlsym(RTLD_NEXT, "malloc");
    if(ghl_real_malloc == NULL) {
      exit(42);
    }
  }
  if(ghl_oom_in_hook) {
    return ghl_real_malloc(size);
  }
  ghl_oom_in_hook = 1;
  const int fail = oom_match(size, RK_MALLOC, __builtin_return_address(0));
  ghl_oom_in_hook = 0;
  if(fail) {
    return NULL;
  }
  void *pointer = ghl_real_malloc(size);
  track_allocation(pointer, __builtin_return_address(0));
  return pointer;
}

void *calloc(size_t nmemb, size_t size) {
  if(ghl_real_calloc == NULL) {
    typedef void *(*real_fn)(size_t, size_t);
    ghl_real_calloc = (real_fn)dlsym(RTLD_NEXT, "calloc");
    if(ghl_real_calloc == NULL) {
      exit(42);
    }
  }
  if(ghl_oom_in_hook || (size != 0 && nmemb > SIZE_MAX / size)) {
    return ghl_real_calloc(nmemb, size);
  }
  ghl_oom_in_hook = 1;
  const int fail = oom_match(nmemb * size, RK_CALLOC, __builtin_return_address(0));
  ghl_oom_in_hook = 0;
  if(fail) {
    return NULL;
  }
  void *pointer = ghl_real_calloc(nmemb, size);
  track_allocation(pointer, __builtin_return_address(0));
  return pointer;
}

static void track_allocation(void *pointer, void *caller_address) {
  if(pointer == NULL || !ghl_oom_armed || ghl_oom_in_hook) {
    return;
  }
  ghl_oom_in_hook = 1;
  if(ghl_caller_is_libghl(caller_address)) {
    struct tracked_allocation *node = ghl_real_malloc(sizeof(*node));
    if(node == NULL) {
      exit(42);
    }
    node->pointer = pointer;
    node->next = tracked_allocations;
    tracked_allocations = node;
  }
  ghl_oom_in_hook = 0;
}

void free(void *pointer) {
  if(ghl_real_free == NULL) {
    ghl_real_free = (void (*)(void *))dlsym(RTLD_NEXT, "free");
    if(ghl_real_free == NULL) {
      exit(42);
    }
  }
  struct tracked_allocation **link = &tracked_allocations;
  while(*link != NULL) {
    struct tracked_allocation *node = *link;
    if(node->pointer == pointer) {
      *link = node->next;
      ghl_real_free(node);
      break;
    }
    link = &node->next;
  }
  ghl_real_free(pointer);
}

static void oom_arm_begin(const long size, const int kind, const long occurrence) {
  ghl_oom_arm_size = size;
  ghl_oom_arm_kind = kind;
  ghl_oom_fail_occurrence = occurrence;
  ghl_oom_seen = 0;
  ghl_oom_fired = 0;
  ghl_oom_armed = 1;
}
static void oom_arm_end(void) { ghl_oom_armed = 0; }

static void init_tab_eos(ghl_eos_parameters *eos) {
  memset(eos, 0, sizeof(*eos));
  eos->eos_type = ghl_eos_tabulated;
  eos->table_type = ghl_eos_table_stellarcollapse;
  eos->clean_sound_speed = false; /* keep reads fast and silent */
}

struct oom_read_run {
  const char *path;
  long size;
  int kind;
  long occurrence;
  ghl_error_codes_t expect_rc;
};

static void oom_read_table_child(void *vp) {
  const struct oom_read_run oc = *(const struct oom_read_run *)vp;
  ghl_eos_parameters eos;
  init_tab_eos(&eos);
  oom_arm_begin(oc.size, oc.kind, oc.occurrence);
  const ghl_error_codes_t rc = NRPyEOS_read_table_set_EOS_params(oc.path, &eos);
  oom_arm_end();
  if(rc != oc.expect_rc) {
    exit(1);
  }
  if(ghl_oom_fired != 1) {
    exit(2);
  }
  /* Whatever failed, the six to_ghl outcome pointers must be NULL: a
   * partial success that left one behind would betray a leak. */
  if(eos.table_logrho != NULL || eos.table_logT != NULL || eos.table_Y_e != NULL
     || eos.table_all != NULL || eos.table_eps != NULL || eos.table_logh != NULL) {
    exit(3);
  }
  if(tracked_allocations != NULL) {
    exit(5);
  }
  exit(0);
}

struct oom_beq_run {
  const char *path;
  long size;
  long occurrence;
  ghl_error_codes_t expect_rc;
  int phase; /* 0=Ye_of_lr 1=munu_of_Ye 2=lp/le/lh triple */
};

static void oom_ye_beq_child(void *vp) {
  const struct oom_beq_run oc = *(const struct oom_beq_run *)vp;
  ghl_eos_parameters eos;
  init_tab_eos(&eos);
  const ghl_error_codes_t rerr = NRPyEOS_read_table_set_EOS_params(oc.path, &eos);
  if(rerr != ghl_success) {
    exit(4);
  }
  oom_arm_begin(oc.size, RK_MALLOC, oc.occurrence);
  const double T = 0.5 * (eos.table_T_min + eos.table_T_max);
  const ghl_error_codes_t rc
        = NRPyEOS_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(T, &eos);
  oom_arm_end();
  if(rc != oc.expect_rc) {
    exit(1);
  }
  if(ghl_oom_fired != 1) {
    exit(2);
  }
  /* Output preservation: munu_of_Ye scratch failure must not publish the
   * lp/le/lh triple; a triple failure must keep Ye_of_lr/table_all and
   * leave the triple unpublished. */
  if(eos.lp_of_lr != NULL || eos.le_of_lr != NULL || eos.lh_of_lr != NULL) {
    exit(3);
  }
  if(oc.phase == 0 && eos.Ye_of_lr != NULL) {
    exit(3);
  }
  if(oc.phase == 2 && (eos.Ye_of_lr == NULL || eos.table_all == NULL)) {
    exit(3);
  }
  NRPyEOS_free_memory(&eos);
  if(tracked_allocations != NULL) {
    exit(5);
  }
  exit(0);
}

/* Dry-run: child counts matching allocations WITHOUT ever failing one
 * (occurrence < 0) and reports the exact count through a pipe. */
static long count_matching_allocations(
      const int key,
      const char *path,
      const long size,
      const int kind,
      const int is_beq) {
  int fds[2];
  if(pipe(fds) != 0) {
    fprintf(stderr, "Test %d failed: pipe(): %s\n", key, strerror(errno));
    exit(1);
  }
  fflush(NULL);
  const pid_t pid = fork();
  if(pid < 0) {
    fprintf(stderr, "Test %d failed: fork(): %s\n", key, strerror(errno));
    exit(1);
  }
  if(pid == 0) {
    close(fds[0]);
    oom_arm_begin(size, kind, -1); /* record only */
    ghl_eos_parameters eos;
    init_tab_eos(&eos);
    if(!is_beq) {
      const ghl_error_codes_t rc = NRPyEOS_read_table_set_EOS_params(path, &eos);
      oom_arm_end();
      if(rc != ghl_success) {
        exit(11);
      }
      NRPyEOS_free_memory(&eos);
    }
    else {
      const ghl_error_codes_t rc = NRPyEOS_read_table_set_EOS_params(path, &eos);
      oom_arm_end();
      if(rc != ghl_success) {
        exit(11);
      }
      oom_arm_begin(size, kind, -1);
      const ghl_error_codes_t brc
            = NRPyEOS_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(
                  0.5 * (eos.table_T_min + eos.table_T_max), &eos);
      oom_arm_end();
      if(brc != ghl_success) {
        exit(12);
      }
      NRPyEOS_free_memory(&eos);
    }
    const long k = ghl_oom_seen;
    ssize_t off = 0;
    while(off < (ssize_t)sizeof(k)) {
      const ssize_t w = write(fds[1], (const char *)&k + off, sizeof(k) - off);
      if(w < 0 && errno == EINTR) {
        continue;
      }
      if(w <= 0) {
        exit(13);
      }
      off += w;
    }
    exit(0);
  }
  close(fds[1]);
  long k = -1;
  ssize_t off = 0;
  while(off < (ssize_t)sizeof(k)) {
    const ssize_t r = read(fds[0], (char *)&k + off, sizeof(k) - off);
    if(r < 0 && errno == EINTR) {
      continue;
    }
    if(r == 0) {
      break;
    }
    if(r < 0) {
      fprintf(stderr, "Test %d failed: read(pipe): %s\n", key, strerror(errno));
      exit(1);
    }
    off += r;
  }
  close(fds[0]);
  const int rc = reap_child(key, pid);
  if(rc != 0 || off != (ssize_t)sizeof(k) || k < 0) {
    fprintf(
          stderr, "Test %d failed: dry-run child exit=%d bytes=%zd count=%ld\n", key, rc,
          off, k);
    exit(1);
  }
  return k;
}

static void test_oom_paths(const int key0) {
  char p[PATH_MAX];
  path_of(p, sizeof(p), "valid");

  const size_t n = 4;
  const size_t npts = 64;
  const size_t nk = NRPyEOS_ntablekeys;
  const long sz_vec = (long)(sizeof(double) * n);         /*  32: 4-double axis */
  const long sz_pts = (long)(sizeof(double) * npts);      /* 512: 64-double 3D data */
  const long sz_all = (long)(sizeof(double) * npts * nk); /* 9728: table_all */
  const long sz_sc = (long)sizeof(NRPyEOS_stellarcollapse_t);

  /* --- GRHayL-owned allocations --------------------------------------------
   * Caller qualification already excludes libhdf5, so the dry-run counts
   * below are exact GRHayL occurrence numbers for the whole read:
   *   - 3 axis reads + 3 to_ghl axis copies at 32 bytes;
   *   - 19 dataset reads + 2 to_ghl (eps, logh) at 512 bytes;
   *   - 1 table_all at npoints*ntablekeys bytes;
   *   - 1 stellarcollapse table struct calloc.
   * Fail each to_ghl and stellarcollapse occurrence once and require
   * ghl_error_out_of_memory plus transactional cleanup to NULL. */
  const long n_g_all = count_matching_allocations(key0, p, sz_all, RK_MALLOC, 0);
  expect_true(
        n_g_all == 1, key0,
        "exactly one npoints*ntablekeys table_all malloc expected per read");

  const long n_sc = count_matching_allocations(key0, p, sz_sc, RK_CALLOC, 0);
  expect_true(
        n_sc == 1, key0,
        "exactly one stellarcollapse table struct calloc expected per read");

  const long n_read_32 = count_matching_allocations(key0, p, sz_vec, RK_MALLOC, 0);
  const long n_read_512 = count_matching_allocations(key0, p, sz_pts, RK_MALLOC, 0);
  expect_true(
        n_read_32 == 6, key0,
        "exactly six 32-byte GRHayL axis mallocs expected per read");
  expect_true(
        n_read_512 == 21, key0, "exactly 21 512-byte GRHayL mallocs expected per read");

  struct oom_read_run oc_read = { p, sz_all, RK_MALLOC, 0, ghl_error_out_of_memory };
  run_in_child(key0, oom_read_table_child, &oc_read, 0);

  oc_read.size = sz_sc;
  oc_read.kind = RK_CALLOC;
  run_in_child(key0, oom_read_table_child, &oc_read, 0);

  /* The to_ghl axis allocations are the final three 32-byte occurrences
   * (after the GRHayL_H5 axis reads), and the to_ghl eps/logh allocations
   * are the final two 512-byte occurrences.  Fail each one; the child
   * requires ghl_error_out_of_memory and every table pointer NULL. */
  for(long occ = 0; occ < n_read_32; occ++) {
    oc_read.size = sz_vec;
    oc_read.kind = RK_MALLOC;
    oc_read.occurrence = occ;
    run_in_child(key0, oom_read_table_child, &oc_read, 0);
  }
  for(long occ = 0; occ < n_read_512; occ++) {
    oc_read.size = sz_pts;
    oc_read.kind = RK_MALLOC;
    oc_read.occurrence = occ;
    run_in_child(key0, oom_read_table_child, &oc_read, 0);
  }

  const long scalar_sizes[] = { sizeof(int), sizeof(double) };
  const long scalar_counts[] = { 4, 1 };
  for(size_t i = 0; i < sizeof(scalar_sizes) / sizeof(scalar_sizes[0]); i++) {
    const long count
          = count_matching_allocations(key0, p, scalar_sizes[i], RK_MALLOC, 0);
    expect_true(count == scalar_counts[i], key0, "scalar allocation count");
    for(long occurrence = 0; occurrence < count; occurrence++) {
      oc_read.size = scalar_sizes[i];
      oc_read.kind = RK_MALLOC;
      oc_read.occurrence = occurrence;
      run_in_child(key0, oom_read_table_child, &oc_read, 0);
    }
  }

  /* --- Beta-equilibrium allocations ----------------------------------------
   * Caller qualification makes every 32-byte malloc inside the beq call
   * GRHayL-owned: Ye_of_lr (1), munu_of_Ye scratch (1), lp/le/lh_of_lr
   * triple (3).  Fail each occurrence once and check output preservation. */
  const long n_beq = count_matching_allocations(key0, p, sz_vec, RK_MALLOC, 1);
  expect_true(n_beq == 5, key0, "exactly five 32-byte beq mallocs expected");

  for(long occ = 0; occ < n_beq; occ++) {
    const int phase = (occ == 0) ? 0 : (occ == 1) ? 1 : 2;
    struct oom_beq_run oc = { p, sz_vec, occ, ghl_error_out_of_memory, phase };
    run_in_child(key0, oom_ye_beq_child, &oc, 0);
  }
}
#endif /* __linux__ */
#endif /* !GHL_DISABLE_HDF5 */

/* ------------------------------------------------------------------ */
int main(void) {
  test_hybrid_find_index_from_P(1);
  test_hybrid_rho_cold_from_P_cold(2);

#ifndef GHL_DISABLE_HDF5
  fixture_setup();
  test_read_table_then_adjust_cs2(3);
  test_read_table_error_paths(4);
  test_tabulated_beta_equilibrium_paths(5);
  test_interpolation_failures_and_enthalpy(8);
#if defined(__linux__)
  test_oom_paths(6);
  test_hdf5_faults(7);
#endif
  fixture_teardown();
#endif

  printf("(GRHayL) unit_test_coverage_eos: EOS boundary tests passed\n");
  return 0;
}
