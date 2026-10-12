#define _POSIX_C_SOURCE 200809L
#include "ghl_unit_tests.h"
#include <errno.h>
#include <float.h>
#include <signal.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

/* Flat-space analytic tensors: a moving unmagnetized fluid, and a static
 * magnetized fluid. The reference uses no production magnetic helper. */
static void check_stress_energy(void) {
  ghl_metric_quantities metric;
  ghl_ADM_aux_quantities aux;
  ghl_initialize_metric(1, 0, 0, 0, 1, 0, 0, 1, 0, 1, &metric);
  ghl_compute_ADM_auxiliaries(&metric, &aux);
  const ghl_primitive_quantities cases[] = {
    { .rho = 2.0, .press = 0.5, .eps = 0.25,
      .u0 = 1.25, .vU = { 0.6, 0.0, 0.0 } },
    { .rho = 2.0, .press = 0.5, .eps = 0.25,
      .u0 = 1.0, .BU = { 0.25, -0.5, 0.75 } }
  };
  for(size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); ++c) {
    const ghl_primitive_quantities *p = &cases[c];
    const double u[4] = { p->u0, p->u0 * p->vU[0], 0.0, 0.0 };
    const double b[4] = { 0.0, p->BU[0], p->BU[1], p->BU[2] };
    const double b2 = b[1] * b[1] + b[2] * b[2] + b[3] * b[3];
    const double inertia = p->rho * (1.0 + p->eps) + p->press + b2;
    ghl_stress_energy tensor;
    ghl_compute_TUPmunu(&metric, &aux, p, &tensor);
    for(int i = 0; i < 4; ++i) {
      for(int j = i; j < 4; ++j) {
        const double gij = i != j ? 0.0 : (i == 0 ? -1.0 : 1.0);
        const double expected = inertia * u[i] * u[j]
                                + (p->press + 0.5 * b2) * gij - b[i] * b[j];
        /* These binary-fraction inputs and references are exactly representable. */
        if(tensor.T4[i][j] != expected) {
          ghl_error("Analytic TUPmunu mismatch, case %zu component %d,%d\n", c, i, j);
        }
      }
    }
  }
}

static void check_simple_initializer_rejection(void) {
  const double valid[] = { 1.0, 0.0, 10.0, 0.5, 0.0, 10.0, 2.0 };
  const double nonfinite[] = { NAN, INFINITY, -INFINITY };
  for(size_t field = 0; field < sizeof(valid) / sizeof(valid[0]); ++field) {
    for(size_t c = 0; c < sizeof(nonfinite) / sizeof(nonfinite[0]); ++c) {
      double args[7];
      memcpy(args, valid, sizeof(args));
      args[field] = nonfinite[c];
      ghl_eos_parameters eos;
      memset(&eos, 0xa5, sizeof(eos));
      unsigned char original[sizeof(eos)];
      memcpy(original, &eos, sizeof(eos));
      const ghl_error_codes_t expected = field == 0 ? ghl_error_invalid_rho_atm
            : field == 3 ? ghl_error_invalid_press_atm : ghl_error_invalid_eos_parameters;
      const ghl_error_codes_t error = ghl_initialize_simple_eos(
            args[0], args[1], args[2], args[3], args[4], args[5], args[6], &eos);
      if(error != expected || memcmp(original, &eos, sizeof(eos)) != 0) {
        ghl_error("Simple EOS rejection failed for argument %zu, case %zu\n", field, c);
      }
    }
  }
}

static void check_tabulated_nonfinite_parameters(void) {
  const double valid[9] = { 1, 0, 10, 0.1, 0, 1, 1, 0, 10 };
  for(size_t c = 0; c < sizeof(valid) / sizeof(valid[0]); ++c) {
    double args[9];
    memcpy(args, valid, sizeof(args));
    args[c] = NAN;
    ghl_eos_parameters eos = { 0 };
    const ghl_error_codes_t error = ghl_initialize_tabulated_eos_functions_and_params(
          "not-read-for-nonfinite-input.h5",
          args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], &eos);
#ifdef GHL_DISABLE_HDF5
    const ghl_error_codes_t expected = ghl_error_used_disabled_hdf5;
#else
    const ghl_error_codes_t expected = ghl_error_invalid_eos_parameters;
#endif
    if(error != expected) {
      ghl_error("Tabulated EOS nonfinite argument %zu returned %d\n", c, error);
    }
  }
}

static void check_hybrid_transition_overflow(void) {
  const double rho[1] = { 8.0 };
  const double gamma[2] = { 2.0, 3.0 };
  ghl_eos_parameters eos;
  memset(&eos, 0xa5, sizeof(eos));
  unsigned char original[sizeof(eos)];
  memcpy(original, &eos, sizeof(eos));
  /* K and the integrated-energy coefficients remain finite, but the
   * transition pressure K*rho^2 = 6.4e308 exceeds binary64 range. */
  const ghl_error_codes_t error = ghl_initialize_hybrid_eos_functions_and_params(
        1.0, 0.0, 10.0, 2, rho, gamma, 1.0e307, 2.0, &eos);
  if(error != ghl_error_invalid_eos_parameters
     || memcmp(original, &eos, sizeof(eos)) != 0) {
    ghl_error("Hybrid pressure-transition overflow was not rejected transactionally\n");
  }
}

static void check_fatal_diagnostics(const bool abort_case) {
  int descriptors[2];
  if(pipe(descriptors) != 0) ghl_error("Fatal-test pipe failed\n");
  fflush(NULL);
  const pid_t child = fork();
  if(child < 0) ghl_error("Fatal-test fork failed\n");
  if(child == 0) {
    close(descriptors[0]);
    if(dup2(descriptors[1], STDERR_FILENO) < 0) _exit(127);
    close(descriptors[1]);
    const struct rlimit no_core = { 0, 0 };
    if(setrlimit(RLIMIT_CORE, &no_core) != 0) _exit(127);
    if(abort_case) {
      ghl_abort("coverage abort diagnostic\n");
    }
    else {
      ghl_abort_if_error(ghl_error_nn_c2p_invalid_number);
    }
    _exit(0);
  }
  close(descriptors[1]);
  char output[2048];
  size_t length = 0;
  while(length + 1 < sizeof(output)) {
    const ssize_t count = read(descriptors[0], output + length, sizeof(output) - length - 1);
    if(count > 0) length += (size_t)count;
    else if(count < 0 && errno == EINTR) continue;
    else if(count < 0) ghl_error("Fatal-test diagnostic read failed\n");
    else break;
  }
  output[length] = '\0';
  close(descriptors[0]);
  int status;
  pid_t waited;
  do { waited = waitpid(child, &status, 0); } while(waited < 0 && errno == EINTR);
  if(waited != child) ghl_error("Fatal-test wait failed\n");
  const bool status_ok = abort_case ? WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT
        : WIFEXITED(status) && WEXITSTATUS(status) == ghl_error_nn_c2p_invalid_number;
  const char *needle = abort_case ? "coverage abort diagnostic" : "non-finite or invalid numeric data";
  if(!status_ok || strstr(output, needle) == NULL) {
    ghl_error("Fatal error lost its status or diagnostic\n");
  }
}

int main(void) {
  check_stress_energy();
  check_simple_initializer_rejection();
  check_tabulated_nonfinite_parameters();
  check_hybrid_transition_overflow();
  check_fatal_diagnostics(false);
  check_fatal_diagnostics(true);
  ghl_abort_if_error(ghl_success);
  /* The legacy switch returned for unrecognized codes; preserve that contract. */
  ghl_abort_if_error((ghl_error_codes_t)-1);
  ghl_abort_if_error((ghl_error_codes_t)(ghl_error_m1_source_double_application + 1));
  ghl_info("Core coverage boundary tests passed.\n");
  return 0;
}
