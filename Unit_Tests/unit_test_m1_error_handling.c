#define _POSIX_C_SOURCE 200809L
#include "ghl_m1.h"
#include "m1_test_utils.h"
#include "m1_thcm1_fixture_utils.h"
#include "m1_thcm1_transport_fixture.h"

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

/*
 * Error-path coverage for shared M1 validation and fatal diagnostics. The
 * test checks error-code mappings, fixture-parser rejection, and transactional
 * output preservation across shared and neutrino operations.
 */

typedef struct {
  ghl_error_codes_t code;
  const char *diagnostic;
} m1_error_case;

static const m1_error_case m1_error_cases[] = {
  { ghl_error_m1_null_pointer, "M1 routine received a NULL pointer" },
  { ghl_error_m1_invalid_epsilon_c,
    "M1 epsilon_c must be finite and strictly between zero and one" },
  { ghl_error_m1_invalid_E_floor, "M1 E_floor must be finite and strictly positive" },
  { ghl_error_m1_invalid_zeta_min, "M1 zeta_min must be finite and strictly positive" },
  { ghl_error_m1_invalid_fd_epsilon_rel,
    "M1 relative finite-difference epsilon must be finite and strictly positive" },
  { ghl_error_m1_invalid_fd_epsilon_abs,
    "M1 absolute finite-difference epsilon must be finite and strictly positive" },
  { ghl_error_m1_invalid_newton_max_iterations,
    "M1 Newton maximum iteration count must be strictly positive" },
  { ghl_error_m1_invalid_newton_tolerance,
    "M1 Newton relative tolerance must be finite and strictly positive" },
  { ghl_error_m1_invalid_metric, "M1 operation received an invalid metric" },
  { ghl_error_m1_invalid_implicit_jacobian,
    "M1 implicit solve produced an invalid Jacobian" },
  { ghl_error_m1_invalid_state, "M1 operation received an invalid radiation state" },
  { ghl_error_m1_implicit_admissibility,
    "M1 implicit trial state is outside the admissible domain" },
  { ghl_error_m1_implicit_solve_failure, "M1 implicit Newton solve failed to converge" },
  { ghl_error_m1_implicit_terminal_fallback,
    "M1 implicit solve exhausted its terminal fallback policy" },
  { ghl_error_m1_con2prim_failure,
    "M1 source update failed during conservative-to-primitive recovery" },
  { ghl_error_m1_microphysics_failure,
    "M1 neutrino microphysics rate evaluation failed" },
  { ghl_error_m1_invalid_newton_absolute_tolerance,
    "M1 Newton absolute tolerance must be finite and strictly positive" },
  { ghl_error_m1_invalid_closure_tolerance,
    "M1 closure tolerance must be finite and strictly positive" },
  { ghl_error_m1_invalid_closure_max_iterations,
    "M1 closure maximum iteration count must be strictly positive" },
  { ghl_error_m1_invalid_repair_policy,
    "M1 received an unsupported realizability repair policy" },
  { ghl_error_m1_closure_residual_too_large,
    "M1 closure residual exceeded the configured acceptance tolerance" },
  { ghl_error_flux_source_invalid_input,
    "Flux/source operation received invalid input" },
  { ghl_error_m1_incompatible_transport_policy,
    "M1 transport policy is incompatible with the requested operation" },
  { ghl_error_m1_source_double_application,
    "M1 source update was applied more than once to the same state" }
};

static void fail_test(const char *message) {
  fprintf(stderr, "unit_test_m1_error_handling: %s\n", message);
  exit(EXIT_FAILURE);
}

static size_t
read_child_output(const int output_fd, char *output, const size_t output_size) {
  size_t length = 0;
  while(length + 1 < output_size) {
    const ssize_t count = read(output_fd, output + length, output_size - length - 1);
    if(count > 0) {
      length += (size_t)count;
      continue;
    }
    if(count < 0 && errno == EINTR) {
      continue;
    }
    if(count < 0) {
      fail_test("reading child diagnostic failed");
    }
    break;
  }
  output[length] = '\0';
  return length;
}

static void check_fatal_error_case(const m1_error_case *restrict error_case) {
  int output_pipe[2];
  if(pipe(output_pipe) != 0) {
    fail_test("could not create child diagnostic pipe");
  }

  const pid_t child = fork();
  if(child < 0) {
    fail_test("could not fork fatal-helper child");
  }

  if(child == 0) {
    close(output_pipe[0]);
    if(dup2(output_pipe[1], STDERR_FILENO) < 0) {
      _exit(127);
    }
    close(output_pipe[1]);
    ghl_abort_if_error(error_case->code);
    _exit(0);
  }

  close(output_pipe[1]);
  char output[1024];
  read_child_output(output_pipe[0], output, sizeof(output));
  close(output_pipe[0]);

  int status = 0;
  if(waitpid(child, &status, 0) < 0) {
    fail_test("could not wait for fatal-helper child");
  }
  if(!WIFEXITED(status) || WEXITSTATUS(status) != (int)error_case->code) {
    fail_test("fatal-helper child did not exit with its error code");
  }
  if(strstr(output, error_case->diagnostic) == NULL) {
    fail_test("fatal-helper diagnostic did not identify the error");
  }
}

static void check_success_returns(void) { ghl_abort_if_error(ghl_success); }

static void check_failed_initializer_stops(void) {
  int output_pipe[2];
  if(pipe(output_pipe) != 0) {
    fail_test("could not create initializer diagnostic pipe");
  }

  const pid_t child = fork();
  if(child < 0) {
    fail_test("could not fork initializer child");
  }

  if(child == 0) {
    close(output_pipe[0]);
    if(dup2(output_pipe[1], STDOUT_FILENO) < 0
       || dup2(output_pipe[1], STDERR_FILENO) < 0) {
      _exit(127);
    }
    close(output_pipe[1]);

    ghl_m1_parameters params = { 0 };
    const ghl_error_codes_t error = ghl_m1_initialize(
          0.0, 1.0e-12, 1.0e-8, 1.0e-6, 1.0e-8, 20, 1.0e-10, &params);
    if(error != ghl_error_m1_invalid_epsilon_c) {
      _exit(126);
    }
    ghl_abort_if_error(error);

    const char continuation[] = "continued\n";
    if(write(STDOUT_FILENO, continuation, sizeof(continuation) - 1)
       != (ssize_t)(sizeof(continuation) - 1)) {
      _exit(127);
    }
    _exit(0);
  }

  close(output_pipe[1]);
  char output[1024];
  read_child_output(output_pipe[0], output, sizeof(output));
  close(output_pipe[0]);

  int status = 0;
  if(waitpid(child, &status, 0) < 0) {
    fail_test("could not wait for initializer child");
  }
  if(!WIFEXITED(status) || WEXITSTATUS(status) != ghl_error_m1_invalid_epsilon_c) {
    fail_test("fatal initializer path did not terminate with its error code");
  }
  if(strstr(output, "continued") != NULL) {
    fail_test("failed initializer reached its continuation");
  }
  if(strstr(output, "M1 epsilon_c must be finite and strictly between zero and one")
     == NULL) {
    fail_test("failed initializer did not emit its diagnostic");
  }
}

/* Synthetic records test the reader, not the physical reference data. */
typedef struct {
  const char *baseline_input;
  const char *perturbed_input;
  const char *baseline_output;
  const char *perturbed_output;
  const char *normalization;
  const char *operation;
  const char *end;
  int version;
  int records;
  int available;
  int reference_status;
} fixture_reader_case;

/* One ID, the message text, and at most one decimal digit per size_t bit. */
#define M1_FIXTURE_ERROR_SIZE                                                 \
  (sizeof("synthetic") + sizeof(size_t) * CHAR_BIT                            \
   + sizeof("invalid fixture record  at index  fixture  component  response " \
            "comparison failed: actual= reference= normalization=")           \
   + 3 * (sizeof("-1.2345678901234567e+") + sizeof(int) * CHAR_BIT))

static int load_synthetic_fixture(
      const fixture_reader_case *test_case,
      m1_thcm1_fixture_collection *collection) {
  char path[] = "/tmp/grhayl-m1-reader-XXXXXX";
  const int fd = mkstemp(path);
  if(fd < 0) {
    fail_test("could not create temporary fixture");
  }
  FILE *file = fdopen(fd, "w");
  if(file == NULL) {
    close(fd);
    unlink(path);
    fail_test("could not open temporary fixture stream");
  }
  fprintf(
        file,
        "M1_THCM1_FIXTURE %d\noperation %s\n"
        "policy A1_A2\nrecord_count %d\n",
        test_case->version, test_case->operation, test_case->records);
  for(int i = 0; i < test_case->records; ++i) {
    /* Repeated IDs deliberately exercise duplicate-record rejection. */
    fprintf(
          file,
          "record_begin\ncase_id synthetic\npair_id synthetic-pair\n"
          "origin directed\nseed_id none\nfamily radiation\n"
          "perturbation radiation\nsensitivity_start 0\nsensitivity_count 1\n"
          "baseline_input %s\nperturbed_input %s\n"
          "baseline_output %s\nperturbed_output %s\nnormalization %s\n"
          "baseline_available %d\nperturbed_available 1\n"
          "baseline_status 0\nperturbed_status 0\n"
          "baseline_reference_status %d\nperturbed_reference_status 0\n"
          "record_end\n",
          test_case->baseline_input, test_case->perturbed_input,
          test_case->baseline_output, test_case->perturbed_output,
          test_case->normalization, test_case->available, test_case->reference_status);
  }
  fputs(test_case->end, file);
  if(fclose(file) != 0) {
    unlink(path);
    fail_test("could not write temporary fixture");
  }
  char error[M1_FIXTURE_ERROR_SIZE];
  const int loaded = m1_thcm1_fixture_load(
        path, "synthetic", 1, 1, collection, error, sizeof(error));
  if(unlink(path) != 0) {
    fail_test("could not remove temporary fixture");
  }
  return loaded;
}

static void check_fixture_reader_and_comparator(void) {
  const fixture_reader_case valid = { .baseline_input = "1 1",
                                      .perturbed_input = "1 2",
                                      .baseline_output = "1 3",
                                      .perturbed_output = "1 4",
                                      .normalization = "1 1",
                                      .operation = "synthetic",
                                      .end = "end\n",
                                      .version = 1,
                                      .records = 1,
                                      .available = 1,
                                      .reference_status = 0 };
  m1_thcm1_fixture_collection collection = { 0 };
  if(!load_synthetic_fixture(&valid, &collection)) {
    fail_test("valid synthetic fixture was rejected");
  }
  m1_thcm1_fixture_record *record = &collection.records[0];
  double baseline[] = { 3.0 }, perturbed[] = { 4.0 };
  const double normalization[] = { 1.0 };
  char error[M1_FIXTURE_ERROR_SIZE];
  m1_thcm1_fixture_comparison_report report;
  if(!m1_thcm1_fixture_compare_paired(
           record, normalization, baseline, perturbed, &report, error, sizeof(error))) {
    fail_test("exact synthetic pair was rejected");
  }

  /* The offline replay compares one current baseline result with the
   * trusted/perturbed THC envelope; it must not require a current perturbed
   * endpoint. */
  double current_baseline = 3.0 + 1.0e-13;
  if(!m1_thcm1_fixture_compare_envelope(
           record, "pointwise_a1_a2_v1", normalization, &current_baseline, &report,
           error, sizeof(error))) {
    fail_test("valid one-endpoint envelope was rejected");
  }
  /* This discrepancy is inside the stored THC response magnitude (1.0), but
   * outside the base A1/A2 tolerance.  The retained response must participate
   * in acceptance. */
  current_baseline = 3.5;
  if(!m1_thcm1_fixture_compare_envelope(
           record, "pointwise_a1_a2_v1", normalization, &current_baseline, &report,
           error, sizeof(error))) {
    fail_test("baseline discrepancy inside response envelope was rejected");
  }
  current_baseline = 10.0;
  if(m1_thcm1_fixture_compare_envelope(
           record, "pointwise_a1_a2_v1", normalization, &current_baseline, &report,
           error, sizeof(error))) {
    fail_test("baseline discrepancy beyond response envelope was accepted");
  }
  current_baseline = NAN;
  if(m1_thcm1_fixture_compare_envelope(
           record, "pointwise_a1_a2_v1", normalization, &current_baseline, &report,
           error, sizeof(error))) {
    fail_test("nonfinite one-endpoint value was accepted");
  }
  current_baseline = 3.0;
  if(m1_thcm1_fixture_compare_envelope(
           record, "unsupported_policy", normalization, &current_baseline, &report,
           error, sizeof(error))) {
    fail_test("unsupported envelope policy was accepted");
  }
  double control_baseline_input[] = { 1.0, 1.0 };
  double control_perturbed_input[] = { 1.0, 1.0 };
  double control_perturbed_output[] = { 3.0 };
  const char control_perturbation[] = "input_invariant_control:synthetic";
  m1_thcm1_fixture_record control = *record;
  control.baseline_input = control_baseline_input;
  control.perturbed_input = control_perturbed_input;
  control.perturbed_output = control_perturbed_output;
  control.perturbation = (char *)control_perturbation;
  control.sensitivity_start = 0;
  control.sensitivity_count = 0;
  current_baseline = 3.0 + 1.0e-13;
  if(!m1_thcm1_fixture_compare_envelope(
           &control, "pointwise_a1_a2_v1", normalization, &current_baseline, &report,
           error, sizeof(error))) {
    fail_test("zero-response control within policy was rejected");
  }
  double control_nonzero_perturbed_output[] = { 3.5 };
  control.perturbed_output = control_nonzero_perturbed_output;
  if(m1_thcm1_fixture_compare_envelope(
           &control, "pointwise_a1_a2_v1", normalization, &current_baseline, &report,
           error, sizeof(error))) {
    fail_test("nonzero-response control was accepted");
  }
  control.perturbed_output = control_perturbed_output;
  current_baseline = 3.0 + 1.0e-8;
  if(m1_thcm1_fixture_compare_envelope(
           &control, "pointwise_a1_a2_v1", normalization, &current_baseline, &report,
           error, sizeof(error))) {
    fail_test("zero-response control discrepancy was accepted");
  }

  m1_thcm1_fixture_record missing_pair = *record;
  missing_pair.pair_id = NULL;
  current_baseline = 3.0;
  if(m1_thcm1_fixture_compare_envelope(
           &missing_pair, "pointwise_a1_a2_v1", normalization, &current_baseline,
           &report, error, sizeof(error))) {
    fail_test("missing pair metadata was accepted");
  }
  const char mismatched_pair_id[] = "other-pair";
  m1_thcm1_fixture_record mismatched_pair = *record;
  mismatched_pair.pair_id = (char *)mismatched_pair_id;
  if(m1_thcm1_fixture_compare_envelope(
           &mismatched_pair, "pointwise_a1_a2_v1", normalization, &current_baseline,
           &report, error, sizeof(error))) {
    fail_test("mismatched pair metadata was accepted");
  }

  const double saved_nonfinite_output = record->perturbed_output[0];
  record->perturbed_output[0] = NAN;
  if(m1_thcm1_fixture_compare_envelope(
           record, "pointwise_a1_a2_v1", normalization, &current_baseline, &report,
           error, sizeof(error))) {
    fail_test("nonfinite response-envelope output was accepted");
  }
  record->perturbed_output[0] = saved_nonfinite_output;

  const double saved_perturbed_output = record->perturbed_output[0];
  record->perturbed_output[0] = record->baseline_output[0];
  current_baseline = 3.5;
  if(m1_thcm1_fixture_compare_envelope(
           record, "pointwise_a1_a2_v1", normalization, &current_baseline, &report,
           error, sizeof(error))) {
    fail_test("nonzero baseline discrepancy with zero response was accepted");
  }
  current_baseline = 3.0 + 1.0e-13;
  if(!m1_thcm1_fixture_compare_envelope(
           record, "pointwise_a1_a2_v1", normalization, &current_baseline, &report,
           error, sizeof(error))) {
    fail_test("zero-response envelope within policy was rejected");
  }
  current_baseline = 3.0 + 1.0e-8;
  if(m1_thcm1_fixture_compare_envelope(
           record, "pointwise_a1_a2_v1", normalization, &current_baseline, &report,
           error, sizeof(error))) {
    fail_test("zero-response envelope discrepancy was accepted");
  }
  record->perturbed_output[0] = saved_perturbed_output;

  /* These offsets exceed the reviewed A1/A2 bounds for this unit-scale pair. */
  const double state_bound = 2.0e-12 + 2.0e-10 * 4.0;
  baseline[0] = 3.0 + 2.0 * state_bound;
  if(m1_thcm1_fixture_compare_paired(
           record, normalization, baseline, perturbed, &report, error, sizeof(error))) {
    fail_test("baseline-only discrepancy was accepted");
  }
  baseline[0] = 3.0;
  perturbed[0] = 4.0 + 2.0 * state_bound;
  if(m1_thcm1_fixture_compare_paired(
           record, normalization, baseline, perturbed, &report, error, sizeof(error))) {
    fail_test("perturbed-only discrepancy was accepted");
  }
  baseline[0] = 3.0 + 2.0 * state_bound;
  if(m1_thcm1_fixture_compare_paired(
           record, normalization, baseline, perturbed, &report, error, sizeof(error))) {
    fail_test("common offset hidden from the response was accepted");
  }

  /* Each state error is within its bound; their difference exceeds the
   * response bound (2e-12 + 2e-10 for a unit response). */
  baseline[0] = 3.0 - state_bound / 2.0;
  perturbed[0] = 4.0 + state_bound / 2.0;
  if(m1_thcm1_fixture_compare_paired(
           record, normalization, baseline, perturbed, &report, error, sizeof(error))
     || strcmp(report.gate, "response") != 0) {
    fail_test("response-specific discrepancy was not detected");
  }
  baseline[0] = NAN;
  perturbed[0] = 4.0;
  if(m1_thcm1_fixture_compare_paired(
           record, normalization, baseline, perturbed, &report, error, sizeof(error))) {
    fail_test("nonfinite computed value was accepted");
  }

  baseline[0] = 3.0;
  perturbed[0] = 4.0;
  if(!m1_thcm1_transport_compare_paired(
           record, normalization, baseline, perturbed, &report, error, sizeof(error))) {
    fail_test("exact transport pair was rejected");
  }
  baseline[0] -= 2.0e-12 * 2.0;
  perturbed[0] += 2.0e-12 * 2.0;
  if(!m1_thcm1_transport_compare_paired(
           record, normalization, baseline, perturbed, &report, error, sizeof(error))) {
    fail_test("transport response rejected propagated endpoint uncertainty");
  }
  perturbed[0] += 4.0 * 2.0e-12 * 4.0;
  if(m1_thcm1_transport_compare_paired(
           record, normalization, baseline, perturbed, &report, error, sizeof(error))) {
    fail_test("transport perturbed-endpoint discrepancy was accepted");
  }
  /* The retained denominator floor is 1e-300: this difference is half its
   * 2e-12 bound, even though the expected states are exactly zero. */
  record->baseline_output[0] = record->perturbed_output[0] = 0.0;
  baseline[0] = perturbed[0] = 1.0e-300 * 2.0e-12 / 2.0;
  if(!m1_thcm1_transport_compare_paired(
           record, normalization, baseline, perturbed, &report, error, sizeof(error))) {
    fail_test("transport denominator-floor case was rejected");
  }
  baseline[0] *= 4.0;
  if(m1_thcm1_transport_compare_paired(
           record, normalization, baseline, perturbed, &report, error, sizeof(error))) {
    fail_test("transport denominator-floor discrepancy was accepted");
  }
  double span_baseline[] = { 1.0, 2.0, 3.0 };
  double span_perturbed[] = { 2.0, 2.0, 99.0 };
  m1_thcm1_fixture_record invalid_span = *record;
  invalid_span.input_count = 3;
  invalid_span.baseline_input = span_baseline;
  invalid_span.perturbed_input = span_perturbed;
  invalid_span.sensitivity_start = 0;
  invalid_span.sensitivity_count = 1;
  if(m1_thcm1_fixture_record_valid(&invalid_span, 3, 1)) {
    fail_test("changes outside the declared sensitivity span were accepted");
  }
  m1_thcm1_fixture_free(&collection);

  fixture_reader_case invalid = valid;
#define REJECT_FIXTURE(field, value, message)           \
  do {                                                  \
    invalid = valid;                                    \
    invalid.field = value;                              \
    if(load_synthetic_fixture(&invalid, &collection)) { \
      m1_thcm1_fixture_free(&collection);               \
      fail_test(message);                               \
    }                                                   \
  } while(0)
  REJECT_FIXTURE(version, 0, "wrong fixture version was accepted");
  REJECT_FIXTURE(operation, "wrong", "wrong fixture operation was accepted");
  REJECT_FIXTURE(records, 0, "empty fixture was accepted");
  REJECT_FIXTURE(records, 2, "duplicate fixture IDs were accepted");
  REJECT_FIXTURE(end, "", "truncated fixture was accepted");
  REJECT_FIXTURE(end, "end trailing", "trailing fixture data was accepted");
  REJECT_FIXTURE(baseline_input, "2 1 9", "mismatched input lengths were accepted");
  REJECT_FIXTURE(baseline_output, "2 3 9", "mismatched output lengths were accepted");
  REJECT_FIXTURE(baseline_output, "1 nan", "NaN reference was accepted");
  REJECT_FIXTURE(perturbed_output, "1 inf", "infinite reference was accepted");
  REJECT_FIXTURE(perturbed_input, "1 1", "erased sensitivity input was accepted");
  REJECT_FIXTURE(normalization, "1 0", "zero normalization was accepted");
  REJECT_FIXTURE(available, 0, "unavailable reference was accepted");
  REJECT_FIXTURE(reference_status, 1, "failed reference status was accepted");
#undef REJECT_FIXTURE
  if(m1_thcm1_fixture_load(NULL, "synthetic", 1, 1, &collection, error, sizeof(error))) {
    fail_test("null fixture path was accepted");
  }
}

static ghl_error_codes_t initialize_parameter_vector(
      const double values[7],
      const int iterations,
      ghl_m1_parameters *params) {
  return ghl_m1_initialize_with_newton_tolerances(
        values[0], values[1], values[2], values[3], values[4], iterations, values[5],
        values[6], params);
}

static void check_parameter_validation(void) {
  const double valid[] = { 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
  const ghl_error_codes_t invalid_codes[]
        = { ghl_error_m1_invalid_epsilon_c,
            ghl_error_m1_invalid_E_floor,
            ghl_error_m1_invalid_zeta_min,
            ghl_error_m1_invalid_fd_epsilon_rel,
            ghl_error_m1_invalid_fd_epsilon_abs,
            ghl_error_m1_invalid_newton_tolerance,
            ghl_error_m1_invalid_newton_absolute_tolerance };
  const double invalid[] = { 0.0, -1.0, NAN, INFINITY };
  ghl_m1_parameters params = { 0 };
  if(initialize_parameter_vector(valid, 1, &params) != ghl_success) {
    fail_test("valid initialization failed");
  }
  unsigned char before[sizeof(params)];
  memcpy(before, &params, sizeof(params));
#define REQUIRE_PARAMETER_ERROR(call, expected)                                    \
  do {                                                                             \
    if((call) != (expected) || memcmp(before, &params, sizeof(params)) != 0)       \
      fail_test("parameter rejection returned wrong status or modified settings"); \
  } while(0)
  for(size_t field = 0; field < sizeof(valid) / sizeof(valid[0]); ++field) {
    for(size_t bad = 0; bad < sizeof(invalid) / sizeof(invalid[0]); ++bad) {
      double values[sizeof(valid) / sizeof(valid[0])];
      memcpy(values, valid, sizeof(values));
      values[field] = invalid[bad];
      REQUIRE_PARAMETER_ERROR(
            initialize_parameter_vector(values, 1, &params), invalid_codes[field]);
    }
  }
  double epsilon_endpoint[sizeof(valid) / sizeof(valid[0])];
  memcpy(epsilon_endpoint, valid, sizeof(valid));
  epsilon_endpoint[0] = 1.0;
  REQUIRE_PARAMETER_ERROR(
        initialize_parameter_vector(epsilon_endpoint, 1, &params),
        ghl_error_m1_invalid_epsilon_c);
  REQUIRE_PARAMETER_ERROR(
        initialize_parameter_vector(valid, 0, &params),
        ghl_error_m1_invalid_newton_max_iterations);
  REQUIRE_PARAMETER_ERROR(
        initialize_parameter_vector(valid, 1, NULL), ghl_error_m1_null_pointer);
  REQUIRE_PARAMETER_ERROR(
        ghl_m1_set_newton_tolerances(1.0, 1.0, NULL), ghl_error_m1_null_pointer);
  REQUIRE_PARAMETER_ERROR(
        ghl_m1_set_closure_solver_controls(1.0, 1, NULL), ghl_error_m1_null_pointer);
  REQUIRE_PARAMETER_ERROR(
        ghl_m1_set_closure_residual_tolerance(1.0, NULL), ghl_error_m1_null_pointer);
  for(size_t bad = 0; bad < sizeof(invalid) / sizeof(invalid[0]); ++bad) {
    REQUIRE_PARAMETER_ERROR(
          ghl_m1_set_newton_tolerances(invalid[bad], 1.0, &params),
          ghl_error_m1_invalid_newton_tolerance);
    REQUIRE_PARAMETER_ERROR(
          ghl_m1_set_newton_tolerances(1.0, invalid[bad], &params),
          ghl_error_m1_invalid_newton_absolute_tolerance);
    REQUIRE_PARAMETER_ERROR(
          ghl_m1_set_closure_solver_controls(invalid[bad], 1, &params),
          ghl_error_m1_invalid_closure_tolerance);
    REQUIRE_PARAMETER_ERROR(
          ghl_m1_set_closure_residual_tolerance(invalid[bad], &params),
          ghl_error_m1_invalid_closure_tolerance);
  }
  REQUIRE_PARAMETER_ERROR(
        ghl_m1_set_closure_solver_controls(nextafter(1.0, INFINITY), 1, &params),
        ghl_error_m1_invalid_closure_tolerance);
  REQUIRE_PARAMETER_ERROR(
        ghl_m1_set_closure_solver_controls(1.0, 0, &params),
        ghl_error_m1_invalid_closure_max_iterations);
#undef REQUIRE_PARAMETER_ERROR
  if(ghl_m1_set_closure_solver_controls(1.0, 1, &params) != ghl_success
     || params.closure_root_tolerance != 1.0 || params.closure_root_max_iterations != 1
     || ghl_m1_set_closure_residual_tolerance(1.0, &params) != ghl_success
     || params.closure_root_residual_tolerance != 1.0
     || ghl_m1_set_newton_tolerances(0.5, 1.0, &params) != ghl_success
     || params.newton_tolerance != 0.5 || params.newton_absolute_tolerance != 1.0) {
    fail_test("valid parameter setters did not publish settings");
  }
}

static void require_error_code(
      const ghl_error_codes_t actual,
      const ghl_error_codes_t expected,
      const char *message) {
  if(actual != expected) {
    fail_test(message);
  }
}

static void check_shared_null_and_state_contracts(void) {
  ghl_m1_parameters params;
  if(ghl_m1_initialize(0.5, 1.0e-12, 1.0, 1.0e-6, 1.0e-10, 100, 1.0e-10, &params)
     != ghl_success) {
    fail_test("shared-contract parameter initialization failed");
  }

  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  ghl_primitive_quantities prims = { 0 };
  prims.u0 = 1.0;
  const ghl_m1_rad_state rad = { .E = 1.0, .F = { 0.1, 0.0, 0.0 } };
  ghl_m1_closure closure = { 0 };
  if(ghl_m1_compute_closure_with_primitives(&params, &metric, &prims, &rad, &closure)
     != ghl_success) {
    fail_test("shared-contract reference closure failed");
  }

  ghl_m1_comoving comoving = { 0 };
  ghl_m1_sources sources = { 0 };
  ghl_m1_sources interaction = { .S_E = 1.0, .S = { 2.0, 3.0, 4.0 } };
  ghl_stress_energy stress = { 0 };
  ghl_m1_diagnostics diagnostics = { 0 };
  ghl_m1_closure_decomposition_diagnostic decomposition = { 0 };
  ghl_metric_quantities derivatives[3] = { { 0 } };
  ghl_extrinsic_curvature curv = { 0 };
  double E_out = -1.0;
  bool floor_applied = false;

  require_error_code(
        ghl_m1_apply_energy_floor(NULL, 1.0, &E_out, &floor_applied),
        ghl_error_m1_null_pointer, "energy-floor NULL parameter was accepted");
  require_error_code(
        ghl_m1_apply_energy_floor(&params, 1.0, NULL, &floor_applied),
        ghl_error_m1_null_pointer, "energy-floor NULL output was accepted");
  if(ghl_m1_apply_energy_floor(&params, 1.0, &E_out, NULL) != ghl_success
     || E_out != 1.0) {
    fail_test("energy-floor optional diagnostics NULL path failed");
  }
  require_error_code(
        ghl_m1_apply_energy_floor(&params, NAN, &E_out, &floor_applied),
        ghl_error_m1_invalid_state, "nonfinite energy-floor input was accepted");

  ghl_m1_rad_state repair_state = rad;
  require_error_code(
        ghl_m1_realizability_repair(NULL, &metric, &repair_state),
        ghl_error_m1_null_pointer, "repair NULL parameters were accepted");
  require_error_code(
        ghl_m1_realizability_repair(&params, NULL, &repair_state),
        ghl_error_m1_null_pointer, "repair NULL metric was accepted");
  require_error_code(
        ghl_m1_realizability_repair(&params, &metric, NULL), ghl_error_m1_null_pointer,
        "repair NULL state was accepted");
  repair_state.E = NAN;
  require_error_code(
        ghl_m1_realizability_repair(&params, &metric, &repair_state),
        ghl_error_m1_invalid_state, "repair nonfinite energy was accepted");

  require_error_code(
        ghl_m1_compute_closure_minerbo(NULL, &metric, &prims, &rad, &closure),
        ghl_error_m1_null_pointer, "closure NULL parameters were accepted");
  require_error_code(
        ghl_m1_compute_closure_minerbo(&params, NULL, &prims, &rad, &closure),
        ghl_error_m1_null_pointer, "closure NULL metric was accepted");
  require_error_code(
        ghl_m1_compute_closure_minerbo(&params, &metric, NULL, &rad, &closure),
        ghl_error_m1_null_pointer, "closure NULL primitives were accepted");
  require_error_code(
        ghl_m1_compute_closure_minerbo(&params, &metric, &prims, NULL, &closure),
        ghl_error_m1_null_pointer, "closure NULL state was accepted");
  require_error_code(
        ghl_m1_compute_closure_minerbo(&params, &metric, &prims, &rad, NULL),
        ghl_error_m1_null_pointer, "closure NULL output was accepted");

  require_error_code(
        ghl_m1_compute_closure_decomposition_diagnostic(
              NULL, &metric, &prims, &rad, &decomposition),
        ghl_error_m1_null_pointer, "decomposition NULL parameters were accepted");
  require_error_code(
        ghl_m1_compute_closure_decomposition_diagnostic(
              &params, NULL, &prims, &rad, &decomposition),
        ghl_error_m1_null_pointer, "decomposition NULL metric was accepted");
  require_error_code(
        ghl_m1_compute_closure_decomposition_diagnostic(
              &params, &metric, NULL, &rad, &decomposition),
        ghl_error_m1_null_pointer, "decomposition NULL primitives were accepted");
  require_error_code(
        ghl_m1_compute_closure_decomposition_diagnostic(
              &params, &metric, &prims, NULL, &decomposition),
        ghl_error_m1_null_pointer, "decomposition NULL state was accepted");
  require_error_code(
        ghl_m1_compute_closure_decomposition_diagnostic(
              &params, &metric, &prims, &rad, NULL),
        ghl_error_m1_null_pointer, "decomposition NULL output was accepted");

  require_error_code(
        ghl_m1_compute_comoving_moments(
              NULL, &metric, &prims, &rad, &closure, &comoving),
        ghl_error_m1_null_pointer, "comoving NULL parameters were accepted");
  require_error_code(
        ghl_m1_compute_comoving_moments(
              &params, NULL, &prims, &rad, &closure, &comoving),
        ghl_error_m1_null_pointer, "comoving NULL metric was accepted");
  require_error_code(
        ghl_m1_compute_comoving_moments(
              &params, &metric, NULL, &rad, &closure, &comoving),
        ghl_error_m1_null_pointer, "comoving NULL primitives were accepted");
  require_error_code(
        ghl_m1_compute_comoving_moments(
              &params, &metric, &prims, NULL, &closure, &comoving),
        ghl_error_m1_null_pointer, "comoving NULL state was accepted");
  require_error_code(
        ghl_m1_compute_comoving_moments(&params, &metric, &prims, &rad, NULL, &comoving),
        ghl_error_m1_null_pointer, "comoving NULL closure was accepted");
  require_error_code(
        ghl_m1_compute_comoving_moments(&params, &metric, &prims, &rad, &closure, NULL),
        ghl_error_m1_null_pointer, "comoving NULL output was accepted");

  require_error_code(
        ghl_m1_compute_geometry_sources(
              NULL, &metric, derivatives, derivatives + 1, derivatives + 2, &curv, &rad,
              &closure, &sources),
        ghl_error_m1_null_pointer, "geometry NULL parameters were accepted");
  require_error_code(
        ghl_m1_compute_geometry_sources(
              &params, NULL, derivatives, derivatives + 1, derivatives + 2, &curv, &rad,
              &closure, &sources),
        ghl_error_m1_null_pointer, "geometry NULL metric was accepted");
  require_error_code(
        ghl_m1_compute_geometry_sources(
              &params, &metric, NULL, derivatives + 1, derivatives + 2, &curv, &rad,
              &closure, &sources),
        ghl_error_m1_null_pointer, "geometry NULL x derivatives were accepted");
  require_error_code(
        ghl_m1_compute_geometry_sources(
              &params, &metric, derivatives, NULL, derivatives + 2, &curv, &rad,
              &closure, &sources),
        ghl_error_m1_null_pointer, "geometry NULL y derivatives were accepted");
  require_error_code(
        ghl_m1_compute_geometry_sources(
              &params, &metric, derivatives, derivatives + 1, NULL, &curv, &rad,
              &closure, &sources),
        ghl_error_m1_null_pointer, "geometry NULL z derivatives were accepted");
  require_error_code(
        ghl_m1_compute_geometry_sources(
              &params, &metric, derivatives, derivatives + 1, derivatives + 2, NULL,
              &rad, &closure, &sources),
        ghl_error_m1_null_pointer, "geometry NULL curvature was accepted");
  require_error_code(
        ghl_m1_compute_geometry_sources(
              &params, &metric, derivatives, derivatives + 1, derivatives + 2, &curv,
              NULL, &closure, &sources),
        ghl_error_m1_null_pointer, "geometry NULL state was accepted");
  require_error_code(
        ghl_m1_compute_geometry_sources(
              &params, &metric, derivatives, derivatives + 1, derivatives + 2, &curv,
              &rad, NULL, &sources),
        ghl_error_m1_null_pointer, "geometry NULL closure was accepted");
  require_error_code(
        ghl_m1_compute_geometry_sources(
              &params, &metric, derivatives, derivatives + 1, derivatives + 2, &curv,
              &rad, &closure, NULL),
        ghl_error_m1_null_pointer, "geometry NULL output was accepted");

  double source_tilde_tau = -1.0;
  double source_tilde_S[3] = { -2.0, -3.0, -4.0 };
  require_error_code(
        ghl_m1_compute_matter_coupling_sources(
              NULL, &interaction, &source_tilde_tau, source_tilde_S),
        ghl_error_m1_null_pointer, "matter-coupling NULL metric was accepted");
  require_error_code(
        ghl_m1_compute_matter_coupling_sources(
              &metric, NULL, &source_tilde_tau, source_tilde_S),
        ghl_error_m1_null_pointer, "matter-coupling NULL input was accepted");
  require_error_code(
        ghl_m1_compute_matter_coupling_sources(
              &metric, &interaction, NULL, source_tilde_S),
        ghl_error_m1_null_pointer, "matter-coupling NULL energy output was accepted");
  require_error_code(
        ghl_m1_compute_matter_coupling_sources(
              &metric, &interaction, &source_tilde_tau, NULL),
        ghl_error_m1_null_pointer, "matter-coupling NULL momentum output was accepted");
  require_error_code(
        ghl_m1_compute_matter_coupling_sources(
              &metric, &interaction, &source_tilde_tau, source_tilde_S),
        ghl_success, "valid matter coupling failed");
  if(source_tilde_tau != -1.0 || source_tilde_S[0] != -2.0 || source_tilde_S[1] != -3.0
     || source_tilde_S[2] != -4.0) {
    fail_test("valid matter coupling changed its numerical result");
  }

  ghl_metric_quantities bad_matter_metric = metric;
  bad_matter_metric.lapse = NAN;
  source_tilde_tau = 15.0;
  source_tilde_S[0] = 16.0;
  source_tilde_S[1] = 17.0;
  source_tilde_S[2] = 18.0;
  require_error_code(
        ghl_m1_compute_matter_coupling_sources(
              &bad_matter_metric, &interaction, &source_tilde_tau, source_tilde_S),
        ghl_error_m1_invalid_metric, "matter coupling accepted an invalid metric");
  if(source_tilde_tau != 15.0 || source_tilde_S[0] != 16.0 || source_tilde_S[1] != 17.0
     || source_tilde_S[2] != 18.0) {
    fail_test("invalid metric changed matter-coupling outputs");
  }
  ghl_m1_sources bad_interaction = interaction;
  bad_interaction.S_E = NAN;
  source_tilde_tau = 19.0;
  source_tilde_S[0] = 20.0;
  source_tilde_S[1] = 21.0;
  source_tilde_S[2] = 22.0;
  require_error_code(
        ghl_m1_compute_matter_coupling_sources(
              &metric, &bad_interaction, &source_tilde_tau, source_tilde_S),
        ghl_error_m1_invalid_state,
        "matter coupling accepted a nonfinite energy source");
  if(source_tilde_tau != 19.0 || source_tilde_S[0] != 20.0 || source_tilde_S[1] != 21.0
     || source_tilde_S[2] != 22.0) {
    fail_test("nonfinite energy source changed matter-coupling outputs");
  }
  for(int component = 0; component < 3; ++component) {
    bad_interaction = interaction;
    bad_interaction.S[component] = NAN;
    source_tilde_tau = 11.0;
    source_tilde_S[0] = 12.0;
    source_tilde_S[1] = 13.0;
    source_tilde_S[2] = 14.0;
    require_error_code(
          ghl_m1_compute_matter_coupling_sources(
                &metric, &bad_interaction, &source_tilde_tau, source_tilde_S),
          ghl_error_m1_invalid_state,
          "matter coupling accepted a nonfinite momentum source");
    if(source_tilde_tau != 11.0 || source_tilde_S[0] != 12.0 || source_tilde_S[1] != 13.0
       || source_tilde_S[2] != 14.0) {
      fail_test("nonfinite momentum source changed matter-coupling outputs");
    }
  }

  require_error_code(
        ghl_m1_compute_stress_energy(NULL, &metric, &rad, &closure, &stress),
        ghl_error_m1_null_pointer, "stress NULL parameters were accepted");
  require_error_code(
        ghl_m1_compute_stress_energy(&params, NULL, &rad, &closure, &stress),
        ghl_error_m1_null_pointer, "stress NULL metric was accepted");
  require_error_code(
        ghl_m1_compute_stress_energy(&params, &metric, NULL, &closure, &stress),
        ghl_error_m1_null_pointer, "stress NULL state was accepted");
  require_error_code(
        ghl_m1_compute_stress_energy(&params, &metric, &rad, NULL, &stress),
        ghl_error_m1_null_pointer, "stress NULL closure was accepted");
  require_error_code(
        ghl_m1_compute_stress_energy(&params, &metric, &rad, &closure, NULL),
        ghl_error_m1_null_pointer, "stress NULL output was accepted");

  require_error_code(
        ghl_m1_compute_diagnostics(NULL, &metric, &rad, &closure, &diagnostics),
        ghl_error_m1_null_pointer, "diagnostics NULL parameters were accepted");
  require_error_code(
        ghl_m1_compute_diagnostics(&params, NULL, &rad, &closure, &diagnostics),
        ghl_error_m1_null_pointer, "diagnostics NULL metric was accepted");
  require_error_code(
        ghl_m1_compute_diagnostics(&params, &metric, NULL, &closure, &diagnostics),
        ghl_error_m1_null_pointer, "diagnostics NULL state was accepted");
  require_error_code(
        ghl_m1_compute_diagnostics(&params, &metric, &rad, NULL, &diagnostics),
        ghl_error_m1_null_pointer, "diagnostics NULL closure was accepted");
  require_error_code(
        ghl_m1_compute_diagnostics(&params, &metric, &rad, &closure, NULL),
        ghl_error_m1_null_pointer, "diagnostics NULL output was accepted");

  double delta_l = -1.0;
  require_error_code(
        ghl_m1_compute_face_normal_delta_l(NULL, ghl_m1_dirn0, 1.0, &delta_l),
        ghl_error_m1_null_pointer, "face-normal NULL metric was accepted");
  require_error_code(
        ghl_m1_compute_face_normal_delta_l(&metric, ghl_m1_dirn0, 1.0, NULL),
        ghl_error_m1_null_pointer, "face-normal NULL output was accepted");
  require_error_code(
        ghl_m1_compute_face_normal_delta_l(
              &metric, (ghl_m1_direction_t)99, 1.0, &delta_l),
        ghl_error_m1_invalid_state, "face-normal invalid direction was accepted");
  require_error_code(
        ghl_m1_compute_face_normal_delta_l(&metric, ghl_m1_dirn0, 0.0, &delta_l),
        ghl_error_m1_invalid_state, "face-normal nonpositive spacing was accepted");

  require_error_code(
        ghl_m1_compute_harmonic_diffusion_coefficient(1.0, 1.0, NULL),
        ghl_error_m1_null_pointer, "harmonic NULL output was accepted");
  require_error_code(
        ghl_m1_compute_harmonic_diffusion_coefficient(NAN, 1.0, &delta_l),
        ghl_error_m1_invalid_state, "harmonic nonfinite left opacity was accepted");
  require_error_code(
        ghl_m1_compute_harmonic_diffusion_coefficient(1.0, 0.0, &delta_l),
        ghl_error_m1_invalid_state, "harmonic nonpositive right opacity was accepted");

  /* The neutrino wrappers validate the state before dereferencing any local
   * projection or forwarding the call.  One state-NULL call covers every
   * wrapper instance without inventing a private test seam. */
  require_error_code(
        ghl_m1_compute_neutrino_Jthick(
              &params, &metric, &prims, NULL, &E_out, &floor_applied),
        ghl_error_m1_null_pointer, "neutrino Jthick NULL state was accepted");
  require_error_code(
        ghl_m1_compute_neutrino_closure(&params, &metric, &prims, NULL, &closure),
        ghl_error_m1_null_pointer, "neutrino closure NULL state was accepted");
  require_error_code(
        ghl_m1_compute_neutrino_comoving_moments(
              &params, &metric, &prims, NULL, &closure, &comoving),
        ghl_error_m1_null_pointer, "neutrino moments NULL state was accepted");
  require_error_code(
        ghl_m1_compute_neutrino_stress_energy(&params, &metric, NULL, &closure, &stress),
        ghl_error_m1_null_pointer, "neutrino stress NULL state was accepted");
  require_error_code(
        ghl_m1_compute_neutrino_geometry_sources(
              &params, &metric, derivatives, derivatives + 1, derivatives + 2, &curv,
              NULL, &closure, &sources),
        ghl_error_m1_null_pointer, "neutrino geometry NULL state was accepted");
  require_error_code(
        ghl_m1_compute_neutrino_diagnostics(
              &params, &metric, NULL, &closure, &diagnostics),
        ghl_error_m1_null_pointer, "neutrino diagnostics NULL state was accepted");
}

static void check_shared_numeric_contracts(void) {
  ghl_m1_parameters params;
  if(ghl_m1_initialize(0.5, 1.0e-12, 1.0, 1.0e-6, 1.0e-10, 100, 1.0e-10, &params)
     != ghl_success) {
    fail_test("numeric-contract parameter initialization failed");
  }

  /* The legacy initializer clamps a product that underflows to zero to the
   * documented DBL_MIN absolute tolerance.  All supplied inputs remain
   * finite and strictly positive, so this exercises the public boundary. */
  const double true_min = nextafter(0.0, 1.0);
  ghl_m1_parameters tiny_params;
  if(ghl_m1_initialize(
           true_min, true_min, true_min, true_min, true_min, 1, true_min, &tiny_params)
           != ghl_success
     || tiny_params.newton_absolute_tolerance != DBL_MIN) {
    fail_test("initializer did not clamp an underflowed tolerance product");
  }

  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  ghl_primitive_quantities prims = { 0 };
  prims.u0 = 1.0;
  const ghl_m1_rad_state rad = { .E = 1.0, .F = { 0.1, 0.0, 0.0 } };
  ghl_m1_closure closure = { 0 };
  if(ghl_m1_compute_closure_with_primitives(&params, &metric, &prims, &rad, &closure)
     != ghl_success) {
    fail_test("numeric-contract reference closure failed");
  }

  /* Each mutation isolates one metric predicate while keeping the public
   * face-normal operation's other inputs valid. */
  ghl_metric_quantities bad_metric;
  double delta_l = 17.0;
#define REQUIRE_BAD_METRIC(mutation, message)                                           \
  do {                                                                                  \
    bad_metric = metric;                                                                \
    mutation;                                                                           \
    delta_l = 17.0;                                                                     \
    require_error_code(                                                                 \
          ghl_m1_compute_face_normal_delta_l(&bad_metric, ghl_m1_dirn0, 1.0, &delta_l), \
          ghl_error_m1_invalid_metric, message);                                        \
    if(delta_l != 17.0)                                                                 \
      fail_test("face-normal published after bad metric");                              \
  } while(0)
  REQUIRE_BAD_METRIC(bad_metric.lapse = NAN, "nonfinite metric lapse was accepted");
  REQUIRE_BAD_METRIC(bad_metric.lapse = 0.0, "nonpositive metric lapse was accepted");
  REQUIRE_BAD_METRIC(
        bad_metric.detgamma = NAN, "nonfinite metric determinant was accepted");
  REQUIRE_BAD_METRIC(
        bad_metric.detgamma = 0.0, "nonpositive metric determinant was accepted");
  REQUIRE_BAD_METRIC(
        bad_metric.sqrt_detgamma = NAN,
        "nonfinite metric square-root determinant was accepted");
  REQUIRE_BAD_METRIC(
        bad_metric.sqrt_detgamma = 0.0,
        "nonpositive metric square-root determinant was accepted");
  REQUIRE_BAD_METRIC(bad_metric.betaU[0] = NAN, "nonfinite metric shift was accepted");
  REQUIRE_BAD_METRIC(
        bad_metric.gammaDD[0][0] = NAN, "nonfinite covariant metric was accepted");
  REQUIRE_BAD_METRIC(
        bad_metric.gammaDD[0][0] = 0.0, "non-SPD covariant metric was accepted");
  REQUIRE_BAD_METRIC(
        bad_metric.gammaDD[0][1] = 0.25, "nonsymmetric covariant metric was accepted");
  REQUIRE_BAD_METRIC(
        bad_metric.gammaUU[0][0] = NAN, "nonfinite inverse metric was accepted");
  REQUIRE_BAD_METRIC(
        bad_metric.gammaUU[0][0] = 0.0, "non-SPD inverse metric was accepted");
  REQUIRE_BAD_METRIC(
        bad_metric.gammaUU[0][1] = 0.25, "nonsymmetric inverse metric was accepted");
  REQUIRE_BAD_METRIC(
        bad_metric.detgamma = 2.0, "inconsistent metric determinant was accepted");
  REQUIRE_BAD_METRIC(
        bad_metric.sqrt_detgamma = 2.0, "inconsistent metric square root was accepted");
  REQUIRE_BAD_METRIC(
        bad_metric.gammaUU[0][0] = 1.01, "inconsistent inverse metric was accepted");
#undef REQUIRE_BAD_METRIC

  /* validate_parameters is reached through the public repair operation. */
  ghl_m1_rad_state repair_state;
  ghl_m1_parameters bad_params = params;
  bad_params.E_floor = NAN;
  repair_state = rad;
  require_error_code(
        ghl_m1_realizability_repair(&bad_params, &metric, &repair_state),
        ghl_error_m1_invalid_E_floor, "repair accepted a nonfinite E_floor");
  bad_params = params;
  bad_params.repair_policy = (ghl_m1_repair_policy_t)0;
  repair_state = rad;
  require_error_code(
        ghl_m1_realizability_repair(&bad_params, &metric, &repair_state),
        ghl_error_m1_invalid_repair_policy,
        "repair accepted an unsupported repair policy");
  bad_params = params;
  bad_params.epsilon_c = NAN;
  repair_state = rad;
  require_error_code(
        ghl_m1_realizability_repair(&bad_params, &metric, &repair_state),
        ghl_error_m1_invalid_epsilon_c, "repair accepted a nonfinite epsilon_c");
  bad_params = params;
  bad_params.one_minus_epsilon_c_sq = NAN;
  repair_state = rad;
  require_error_code(
        ghl_m1_realizability_repair(&bad_params, &metric, &repair_state),
        ghl_error_m1_invalid_epsilon_c, "repair accepted a nonfinite cone limit");
  bad_params = params;
  bad_params.one_minus_epsilon_c_sq = 0.25;
  repair_state = rad;
  require_error_code(
        ghl_m1_realizability_repair(&bad_params, &metric, &repair_state),
        ghl_error_m1_invalid_epsilon_c, "repair accepted an inconsistent cone limit");
  bad_params = params;
  bad_params.closure_root_tolerance = NAN;
  repair_state = rad;
  require_error_code(
        ghl_m1_realizability_repair(&bad_params, &metric, &repair_state),
        ghl_error_m1_invalid_closure_tolerance,
        "repair accepted a nonfinite closure tolerance");
  bad_params = params;
  bad_params.closure_root_max_iterations = 0;
  repair_state = rad;
  require_error_code(
        ghl_m1_realizability_repair(&bad_params, &metric, &repair_state),
        ghl_error_m1_invalid_closure_max_iterations,
        "repair accepted zero closure iterations");
  bad_params = params;
  bad_params.closure_root_residual_tolerance = NAN;
  repair_state = rad;
  require_error_code(
        ghl_m1_realizability_repair(&bad_params, &metric, &repair_state),
        ghl_error_m1_invalid_closure_tolerance,
        "repair accepted a nonfinite residual tolerance");

  ghl_m1_rad_state bad_rad = rad;
  bad_rad.E = NAN;
  require_error_code(
        ghl_m1_compute_closure_with_primitives(
              &params, &metric, &prims, &bad_rad, &closure),
        ghl_error_m1_invalid_state, "closure accepted nonfinite energy");
  bad_rad = rad;
  bad_rad.F[0] = NAN;
  require_error_code(
        ghl_m1_compute_closure_with_primitives(
              &params, &metric, &prims, &bad_rad, &closure),
        ghl_error_m1_invalid_state, "closure accepted nonfinite flux");
  bad_rad = rad;
  bad_rad.F[0] = 0.9;
  require_error_code(
        ghl_m1_compute_closure_with_primitives(
              &params, &metric, &prims, &bad_rad, &closure),
        ghl_error_m1_invalid_state, "closure accepted a flux outside its cone");

  ghl_primitive_quantities bad_prims = prims;
  bad_prims.vU[0] = NAN;
  require_error_code(
        ghl_m1_compute_closure_with_primitives(
              &params, &metric, &bad_prims, &rad, &closure),
        ghl_error_m1_invalid_state, "closure accepted nonfinite velocity");
  bad_prims = prims;
  bad_prims.vU[0] = 2.0;
  require_error_code(
        ghl_m1_compute_closure_with_primitives(
              &params, &metric, &bad_prims, &rad, &closure),
        ghl_error_u0_singular, "closure accepted a superluminal velocity");

  /* Supplied-closure validation is shared by stress and diagnostics. */
  ghl_stress_energy stress = { 0 };
  ghl_m1_diagnostics diagnostics = { 0 };
  ghl_m1_closure bad_closure = closure;
  bad_closure.P[0][0] = NAN;
  require_error_code(
        ghl_m1_compute_stress_energy(&params, &metric, &rad, &bad_closure, &stress),
        ghl_error_m1_invalid_state, "stress accepted a nonfinite pressure tensor");
  bad_closure = closure;
  bad_closure.P[0][0] += 0.1;
  require_error_code(
        ghl_m1_compute_stress_energy(&params, &metric, &rad, &bad_closure, &stress),
        ghl_error_m1_invalid_state, "stress accepted a trace-inconsistent tensor");
  bad_closure = closure;
  bad_closure.P[0][1] += 0.1;
  require_error_code(
        ghl_m1_compute_stress_energy(&params, &metric, &rad, &bad_closure, &stress),
        ghl_error_m1_invalid_state, "stress accepted an asymmetric tensor");
  bad_closure = closure;
  bad_closure.P[0][0] = 1.0;
  bad_closure.P[1][1] = 1.0;
  bad_closure.P[2][2] = -1.0;
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      if(i != j) {
        bad_closure.P[i][j] = 0.0;
      }
    }
  }
  require_error_code(
        ghl_m1_compute_stress_energy(&params, &metric, &rad, &bad_closure, &stress),
        ghl_error_m1_invalid_state, "stress accepted a non-PSD tensor");

  bad_closure = closure;
  bad_closure.xi = NAN;
  require_error_code(
        ghl_m1_compute_diagnostics(&params, &metric, &rad, &bad_closure, &diagnostics),
        ghl_error_m1_invalid_state, "diagnostics accepted nonfinite closure xi");
  bad_closure = closure;
  bad_closure.solve_status = ghl_m1_closure_solve_invalid;
  require_error_code(
        ghl_m1_compute_diagnostics(&params, &metric, &rad, &bad_closure, &diagnostics),
        ghl_error_m1_invalid_state, "diagnostics accepted an invalid solve status");

  /* Finite source inputs may still overflow the densitized publication. */
  ghl_metric_quantities large_volume = metric;
  for(int i = 0; i < 3; ++i) {
    large_volume.gammaDD[i][i] = 4.0;
    large_volume.gammaUU[i][i] = 0.25;
  }
  large_volume.detgamma = 64.0;
  large_volume.sqrt_detgamma = 8.0;
  ghl_m1_sources interaction = { .S_E = DBL_MAX, .S = { 0.0, 0.0, 0.0 } };
  double source_tilde_tau = 21.0;
  double source_tilde_S[3] = { 22.0, 23.0, 24.0 };
  require_error_code(
        ghl_m1_compute_matter_coupling_sources(
              &large_volume, &interaction, &source_tilde_tau, source_tilde_S),
        ghl_error_m1_invalid_state,
        "matter coupling published an infinite energy source");
  if(source_tilde_tau != 21.0 || source_tilde_S[0] != 22.0 || source_tilde_S[1] != 23.0
     || source_tilde_S[2] != 24.0) {
    fail_test("energy overflow changed matter-coupling outputs");
  }

  ghl_metric_quantities derivatives[3] = { { 0 } };
  ghl_extrinsic_curvature curv = { 0 };
  ghl_m1_sources geometry_sources = { 0 };
  ghl_m1_closure large_volume_closure = { 0 };
  if(ghl_m1_compute_closure_with_primitives(
           &params, &large_volume, &prims, &rad, &large_volume_closure)
     != ghl_success) {
    fail_test("large-volume reference closure failed");
  }
  derivatives[1].lapse = DBL_MAX;
  geometry_sources.S_E = 71.0;
  geometry_sources.S[0] = 72.0;
  geometry_sources.S[1] = 73.0;
  geometry_sources.S[2] = 74.0;
  require_error_code(
        ghl_m1_compute_geometry_sources(
              &params, &large_volume, derivatives, derivatives + 1, derivatives + 2,
              &curv, &rad, &large_volume_closure, &geometry_sources),
        ghl_error_m1_invalid_state, "geometry published an infinite momentum source");
  if(geometry_sources.S_E != 71.0 || geometry_sources.S[0] != 72.0
     || geometry_sources.S[1] != 73.0 || geometry_sources.S[2] != 74.0) {
    fail_test("geometry published output after an overflow rejection");
  }

  derivatives[0].lapse = NAN;
  require_error_code(
        ghl_m1_compute_geometry_sources(
              &params, &large_volume, derivatives, derivatives + 1, derivatives + 2,
              &curv, &rad, &large_volume_closure, &geometry_sources),
        ghl_error_m1_invalid_metric, "geometry accepted a nonfinite lapse derivative");
  derivatives[0].lapse = 0.0;
  derivatives[0].betaU[0] = NAN;
  require_error_code(
        ghl_m1_compute_geometry_sources(
              &params, &large_volume, derivatives, derivatives + 1, derivatives + 2,
              &curv, &rad, &large_volume_closure, &geometry_sources),
        ghl_error_m1_invalid_metric, "geometry accepted a nonfinite shift derivative");
  derivatives[0].betaU[0] = 0.0;
  derivatives[0].gammaDD[0][0] = NAN;
  require_error_code(
        ghl_m1_compute_geometry_sources(
              &params, &large_volume, derivatives, derivatives + 1, derivatives + 2,
              &curv, &rad, &large_volume_closure, &geometry_sources),
        ghl_error_m1_invalid_metric, "geometry accepted a nonfinite metric derivative");
  derivatives[0].gammaDD[0][0] = 0.0;
  curv.K[0][0] = NAN;
  require_error_code(
        ghl_m1_compute_geometry_sources(
              &params, &large_volume, derivatives, derivatives + 1, derivatives + 2,
              &curv, &rad, &large_volume_closure, &geometry_sources),
        ghl_error_m1_invalid_metric, "geometry accepted nonfinite curvature");

  ghl_metric_quantities large_shift = metric;
  large_shift.betaU[0] = DBL_MAX;
  require_error_code(
        ghl_m1_compute_stress_energy(&params, &large_shift, &rad, &closure, &stress),
        ghl_error_m1_invalid_state, "stress published an infinite tensor entry");

  ghl_metric_quantities tiny_lapse = metric;
  tiny_lapse.lapse = true_min;
  require_error_code(
        ghl_m1_compute_stress_energy(&params, &tiny_lapse, &rad, &closure, &stress),
        ghl_error_m1_invalid_metric, "stress accepted an overflowing lapse inverse");

  /* A finite state may have a flux norm above one before canonical repair.
   * The narrow but coherent metric makes original_flux_factor * E overflow,
   * exercising the diagnostic saturation without changing the state contract. */
  ghl_metric_quantities narrow_metric = metric;
  for(int i = 0; i < 3; ++i) {
    narrow_metric.gammaDD[i][i] = 0.25;
    narrow_metric.gammaUU[i][i] = 4.0;
  }
  narrow_metric.detgamma = 0.015625;
  narrow_metric.sqrt_detgamma = 0.125;
  ghl_m1_rad_state extreme_repair = { .E = DBL_MAX, .F = { DBL_MAX, 0.0, 0.0 } };
  require_error_code(
        ghl_m1_realizability_repair(&params, &narrow_metric, &extreme_repair),
        ghl_success, "extreme finite state was not repaired");
  if(extreme_repair.E != DBL_MAX || !isfinite(extreme_repair.F[0])) {
    fail_test("extreme repair did not publish a finite repaired state");
  }

  /* Positive finite opacities whose reciprocal overflows are rejected before
   * publication; the output sentinel verifies the transactional boundary. */
  double diffusion = 23.0;
  require_error_code(
        ghl_m1_compute_harmonic_diffusion_coefficient(true_min, 1.0, &diffusion),
        ghl_error_m1_invalid_state,
        "harmonic coefficient accepted an overflowing reciprocal");
  if(diffusion != 23.0) {
    fail_test("harmonic coefficient changed output after rejection");
  }
}

/* Validate the same public parameter contract at each caller boundary. This
 * exercises the independently compiled header validators in those callers;
 * testing only the initializer cannot establish their rejection behavior. */
static void check_parameter_contract_at_callers(
      const ghl_m1_parameters *params,
      const ghl_error_codes_t expected) {
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  ghl_primitive_quantities prims = { .u0 = 1.0 };
  const ghl_m1_rad_state rad = { .E = 1.0 };
  ghl_m1_rad_state repaired = rad;
  const ghl_m1_closure closure = {
    .P = { { 1.0 / 3.0, 0.0, 0.0 }, { 0.0, 1.0 / 3.0, 0.0 }, { 0.0, 0.0, 1.0 / 3.0 } },
    .chi = 1.0 / 3.0
  };
  ghl_m1_closure closure_out = { 0 };
  ghl_m1_comoving comoving;
  ghl_m1_closure_decomposition_diagnostic decomposition;
  ghl_m1_diagnostics diagnostics;
  ghl_stress_energy stress;
  ghl_m1_sources sources;
  const ghl_metric_quantities derivative = { 0 };
  const ghl_extrinsic_curvature curvature = { 0 };
  double J = 0.0;
  bool J_valid = false;
  const double zero[3] = { 0 };
  double flux_N, flux_E, flux_F[3], velocity[3];
  const ghl_m1_neutrino_parameters nu = { .N_floor = 1.0e-12 };
  const ghl_m1_neutrino_state ns = { .N = 1.0, .E = 1.0 };
  ghl_m1_neutrino_state ns_repaired = ns;
#define CHECK_CALL(call) require_error_code((call), expected, #call)
  CHECK_CALL(ghl_m1_compute_closure_with_primitives(
        params, &metric, &prims, &rad, &closure_out));
  CHECK_CALL(ghl_m1_compute_closure_decomposition_diagnostic(
        params, &metric, &prims, &rad, &decomposition));
  CHECK_CALL(ghl_m1_compute_comoving_moments(
        params, &metric, &prims, &rad, &closure, &comoving));
  CHECK_CALL(ghl_m1_compute_diagnostics(params, &metric, &rad, &closure, &diagnostics));
  CHECK_CALL(ghl_m1_compute_stress_energy(params, &metric, &rad, &closure, &stress));
  CHECK_CALL(ghl_m1_compute_geometry_sources(
        params, &metric, &derivative, &derivative, &derivative, &curvature, &rad,
        &closure, &sources));
  CHECK_CALL(ghl_m1_compute_Jthick(params, &metric, &prims, &rad, &J, &J_valid));
  CHECK_CALL(ghl_m1_realizability_repair(params, &metric, &repaired));
  CHECK_CALL(ghl_m1_repair_neutrino_state(params, &nu, &metric, &ns_repaired, NULL));
  CHECK_CALL(ghl_m1_compute_neutrino_number_flux(
        params, &nu, &metric, &prims, &ns, flux_F, velocity));
  CHECK_CALL(ghl_m1_compute_neutrino_rusanov_flux(
        params, &nu, &metric, ghl_m1_dirn0, &ns, &ns, &closure, &closure, zero, zero,
        zero, zero, 1.0, &flux_N, &flux_E, flux_F));
#undef CHECK_CALL
  if(expected != ghl_success
     && (memcmp(&repaired, &rad, sizeof(rad)) != 0
         || memcmp(&ns_repaired, &ns, sizeof(ns)) != 0)) {
    fail_test("parameter rejection changed a repair state");
  }
}

static void check_caller_parameter_matrix(void) {
  ghl_m1_parameters valid;
  require_error_code(
        ghl_m1_initialize(0.5, 1.0e-12, 1.0, 1.0e-6, 1.0e-10, 100, 1.0e-10, &valid),
        ghl_success, "caller matrix initialize");
  check_parameter_contract_at_callers(&valid, ghl_success);
#define BAD_PARAMETER(field, value, expected)            \
  do {                                                   \
    ghl_m1_parameters bad = valid;                       \
    bad.field = (value);                                 \
    check_parameter_contract_at_callers(&bad, expected); \
  } while(0)
  BAD_PARAMETER(E_floor, NAN, ghl_error_m1_invalid_E_floor);
  BAD_PARAMETER(E_floor, 0.0, ghl_error_m1_invalid_E_floor);
  BAD_PARAMETER(
        repair_policy, (ghl_m1_repair_policy_t)-1, ghl_error_m1_invalid_repair_policy);
  BAD_PARAMETER(epsilon_c, NAN, ghl_error_m1_invalid_epsilon_c);
  BAD_PARAMETER(epsilon_c, 0.0, ghl_error_m1_invalid_epsilon_c);
  BAD_PARAMETER(epsilon_c, 1.0, ghl_error_m1_invalid_epsilon_c);
  BAD_PARAMETER(one_minus_epsilon_c_sq, NAN, ghl_error_m1_invalid_epsilon_c);
  BAD_PARAMETER(one_minus_epsilon_c_sq, 0.25, ghl_error_m1_invalid_epsilon_c);
  BAD_PARAMETER(one_minus_epsilon_c_sq, 0.75, ghl_error_m1_invalid_epsilon_c);
  BAD_PARAMETER(closure_root_tolerance, NAN, ghl_error_m1_invalid_closure_tolerance);
  BAD_PARAMETER(closure_root_tolerance, 0.0, ghl_error_m1_invalid_closure_tolerance);
  BAD_PARAMETER(
        closure_root_tolerance, nextafter(1.0, INFINITY),
        ghl_error_m1_invalid_closure_tolerance);
  BAD_PARAMETER(
        closure_root_max_iterations, 0, ghl_error_m1_invalid_closure_max_iterations);
  BAD_PARAMETER(
        closure_root_residual_tolerance, NAN, ghl_error_m1_invalid_closure_tolerance);
  BAD_PARAMETER(
        closure_root_residual_tolerance, 0.0, ghl_error_m1_invalid_closure_tolerance);
#undef BAD_PARAMETER
}

static void check_state_contract_at_callers(
      const ghl_m1_parameters *params,
      const ghl_metric_quantities *metric,
      const ghl_m1_rad_state *rad,
      const ghl_m1_closure *closure,
      const bool test_closure_generation,
      const ghl_error_codes_t expected) {
  const ghl_primitive_quantities prims = { .u0 = 1.0 };
  ghl_m1_closure closure_out = { 0 };
  ghl_m1_comoving comoving;
  ghl_m1_closure_decomposition_diagnostic decomposition;
  ghl_m1_diagnostics diagnostics;
  ghl_stress_energy stress;
  ghl_m1_sources sources;
  const ghl_metric_quantities derivative = { 0 };
  const ghl_extrinsic_curvature curvature = { 0 };
  double J;
  bool J_valid;
#define CHECK_STATE(call) require_error_code((call), expected, #call)
  if(test_closure_generation) {
    CHECK_STATE(ghl_m1_compute_closure_with_primitives(
          params, metric, &prims, rad, &closure_out));
    CHECK_STATE(ghl_m1_compute_closure_decomposition_diagnostic(
          params, metric, &prims, rad, &decomposition));
    CHECK_STATE(ghl_m1_compute_Jthick(params, metric, &prims, rad, &J, &J_valid));
  }
  CHECK_STATE(ghl_m1_compute_comoving_moments(
        params, metric, &prims, rad, closure, &comoving));
  CHECK_STATE(ghl_m1_compute_diagnostics(params, metric, rad, closure, &diagnostics));
  CHECK_STATE(ghl_m1_compute_stress_energy(params, metric, rad, closure, &stress));
  CHECK_STATE(ghl_m1_compute_geometry_sources(
        params, metric, &derivative, &derivative, &derivative, &curvature, rad, closure,
        &sources));
#undef CHECK_STATE
}

static void check_caller_state_matrix(void) {
  ghl_m1_parameters params;
  require_error_code(
        ghl_m1_initialize(0.5, 1.0e-12, 1.0, 1.0e-6, 1.0e-10, 100, 1.0e-10, &params),
        ghl_success, "state matrix initialize");
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  const ghl_m1_rad_state rad = { .E = 1.0 };
  const ghl_m1_closure closure = {
    .P = { { 1.0 / 3.0, 0.0, 0.0 }, { 0.0, 1.0 / 3.0, 0.0 }, { 0.0, 0.0, 1.0 / 3.0 } },
    .chi = 1.0 / 3.0
  };
  check_state_contract_at_callers(&params, &metric, &rad, &closure, true, ghl_success);
  metric.lapse = -1.0;
  check_state_contract_at_callers(
        &params, &metric, &rad, &closure, true, ghl_error_m1_invalid_metric);
  m1_setup_flat_metric(&metric);
  ghl_m1_rad_state bad = rad;
  bad.E = NAN;
  check_state_contract_at_callers(
        &params, &metric, &bad, &closure, true, ghl_error_m1_invalid_state);
  bad.E = 0.0;
  check_state_contract_at_callers(
        &params, &metric, &bad, &closure, true, ghl_error_m1_invalid_state);
  for(int axis = 0; axis < 3; ++axis) {
    bad = rad;
    bad.F[axis] = NAN;
    check_state_contract_at_callers(
          &params, &metric, &bad, &closure, true, ghl_error_m1_invalid_state);
    bad.F[axis] = 2.0;
    check_state_contract_at_callers(
          &params, &metric, &bad, &closure, true, ghl_error_m1_invalid_state);
    ghl_m1_closure pressure = closure;
    pressure.P[axis][axis] = NAN;
    check_state_contract_at_callers(
          &params, &metric, &rad, &pressure, false, ghl_error_m1_invalid_state);
  }
  ghl_m1_closure pressure = closure;
  pressure.P[0][1] = 0.1;
  check_state_contract_at_callers(
        &params, &metric, &rad, &pressure, false, ghl_error_m1_invalid_state);
  pressure = closure;
  pressure.P[0][0] = -1.0;
  pressure.P[1][1] = pressure.P[2][2] = 1.0;
  check_state_contract_at_callers(
        &params, &metric, &rad, &pressure, false, ghl_error_m1_invalid_state);
}

static void check_remaining_range_and_axis_failures(void) {
  const double tiny = nextafter(0.0, 1.0);
  /* Each pair reaches a distinct reciprocal, sum, or final-product failure;
   * the public inputs themselves are finite and positive. */
  const double opacity_pairs[][2] = { { 1.0, tiny },
                                      { 0.5 / DBL_MAX, 0.5 / DBL_MAX },
                                      { DBL_MAX, DBL_MAX },
                                      { 1.0e-200, 1.0e-200 },
                                      { 1.0e200, 1.0e200 } };
  double value;
  for(size_t i = 0; i < sizeof(opacity_pairs) / sizeof(opacity_pairs[0]); ++i) {
    require_error_code(
          ghl_m1_compute_harmonic_diffusion_coefficient(
                opacity_pairs[i][0], opacity_pairs[i][1], &value),
          ghl_error_m1_invalid_state, "harmonic unrepresentable result accepted");
  }
  require_error_code(
        ghl_m1_compute_harmonic_diffusion_coefficient(1.0, NAN, &value),
        ghl_error_m1_invalid_state, "nonfinite right opacity");

  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  require_error_code(
        ghl_m1_compute_face_normal_delta_l(&metric, ghl_m1_dirn0, NAN, &value),
        ghl_error_m1_invalid_state, "nonfinite coordinate extent");
  require_error_code(
        ghl_m1_compute_face_normal_delta_l(&metric, ghl_m1_dirn0, 0.0, &value),
        ghl_error_m1_invalid_state, "zero coordinate extent");
  ghl_metric_quantities bad_metric = metric;
  memset(bad_metric.gammaDD, 0, sizeof(bad_metric.gammaDD));
  require_error_code(
        ghl_m1_compute_face_normal_delta_l(&bad_metric, ghl_m1_dirn0, 1.0, &value),
        ghl_error_m1_invalid_metric,
        "zero matrix with inconsistent positive determinant accepted");
  bad_metric = metric;
  bad_metric.gammaDD[0][0] = tiny;
  bad_metric.gammaDD[0][1] = bad_metric.gammaDD[1][0] = 1.0;
  require_error_code(
        ghl_m1_compute_face_normal_delta_l(&bad_metric, ghl_m1_dirn0, 1.0, &value),
        ghl_error_m1_invalid_metric,
        "finite matrix with overflowing Cholesky pivot accepted");
  /* Unit determinant, mutually inverse diagonal tensors. Coordinate extent
   * overflow/underflow is independent of metric admissibility. */
  for(int case_index = 0; case_index < 2; ++case_index) {
    const double factor = case_index == 0 ? 1.0e200 : 1.0e-200;
    metric.gammaDD[0][0] = factor;
    metric.gammaUU[0][0] = 1.0 / factor;
    metric.gammaDD[1][1] = metric.gammaDD[2][2] = 1.0 / sqrt(factor);
    metric.gammaUU[1][1] = metric.gammaUU[2][2] = sqrt(factor);
    require_error_code(
          ghl_m1_compute_face_normal_delta_l(
                &metric, ghl_m1_dirn0, case_index == 0 ? DBL_MAX : tiny, &value),
          ghl_error_m1_invalid_state, "unrepresentable proper extent accepted");
  }

  m1_setup_flat_metric(&metric);
  ghl_m1_parameters params;
  require_error_code(
        ghl_m1_initialize(0.5, 1.0e-12, 1.0, 1.0e-6, 1.0e-10, 100, 1.0e-10, &params),
        ghl_success, "axis checks initialize");
  /* Isolate each momentum component in the downstream-repair predicate. */
  for(int axis = 0; axis < 3; ++axis) {
    ghl_m1_rad_state repaired = { .E = 1.0 };
    repaired.F[axis] = 2.0;
    require_error_code(
          ghl_m1_realizability_repair(&params, &metric, &repaired), ghl_success,
          "single-axis super-cone repair failed");
    if(repaired.E != 1.0 || repaired.F[axis] != 0.25) {
      ghl_error("Single-axis repair violated the squared-ratio rule\n");
    }
    for(int other = 0; other < 3; ++other) {
      if(other != axis && repaired.F[other] != 0.0) {
        ghl_error("Single-axis repair changed a zero component\n");
      }
    }
  }
  const ghl_m1_rad_state rad = { .E = 1.0 };
  const ghl_m1_closure closure = {
    .P = { { 1.0 / 3.0, 0.0, 0.0 }, { 0.0, 1.0 / 3.0, 0.0 }, { 0.0, 0.0, 1.0 / 3.0 } },
    .chi = 1.0 / 3.0
  };
  const ghl_extrinsic_curvature curvature = { 0 };
  const ghl_m1_sources original = { .S_E = -3.0, .S = { -2.0, -1.0, 1.0 } };
  for(int axis = 0; axis < 3; ++axis) {
    ghl_metric_quantities derivatives[3] = { { 0 } };
    derivatives[axis].lapse = NAN;
    ghl_m1_sources sources = original;
    require_error_code(
          ghl_m1_compute_geometry_sources(
                &params, &metric, &derivatives[0], &derivatives[1], &derivatives[2],
                &curvature, &rad, &closure, &sources),
          ghl_error_m1_invalid_metric, "nonfinite derivative axis accepted");
    if(memcmp(&sources, &original, sizeof(sources)) != 0) {
      fail_test("derivative rejection changed sources");
    }
  }
  metric.lapse = 2.0;
  metric.lapseinv = 0.5;
  metric.lapseinv2 = 0.25;
  ghl_m1_sources interaction = { .S_E = 1.0, .S = { DBL_MAX, 0.0, 0.0 } };
  double matter_energy = 31.0;
  double momentum[3] = { 32.0, 33.0, 34.0 };
  require_error_code(
        ghl_m1_compute_matter_coupling_sources(
              &metric, &interaction, &matter_energy, momentum),
        ghl_error_m1_invalid_state, "unrepresentable matter momentum accepted");
  if(matter_energy != 31.0 || momentum[0] != 32.0 || momentum[1] != 33.0
     || momentum[2] != 34.0) {
    fail_test("momentum overflow changed matter-coupling outputs");
  }
}

static void check_remaining_shared_rejections(void) {
  ghl_m1_parameters params;
  require_error_code(
        ghl_m1_initialize(0.5, 1.0e-12, 1.0, 1.0e-6, 1.0e-10, 100, 1.0e-10, &params),
        ghl_success, "shared rejection initializer failed");
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  ghl_primitive_quantities prims = { .u0 = 1.0 };
  const ghl_m1_rad_state rad = { .E = 1.0, .F = { 0.25, 0.0, 0.0 } };
  ghl_m1_closure closure = { 0 };
  require_error_code(
        ghl_m1_compute_closure_with_primitives(&params, &metric, &prims, &rad, &closure),
        ghl_success, "shared rejection reference closure failed");

  /* Change just the scalar metadata: the pressure remains admissible, so
   * these calls reach the documented diagnostics metadata checks. Each
   * rejected call must leave the caller's diagnostic packet untouched. */
  const ghl_m1_diagnostics sentinel = { .closure_xi = -17.0, .r = -19.0 };
  ghl_m1_diagnostics diagnostics;
#define REJECT_CLOSURE_METADATA(field, value)                                     \
  do {                                                                            \
    ghl_m1_closure bad = closure;                                                 \
    bad.field = (value);                                                          \
    diagnostics = sentinel;                                                       \
    require_error_code(                                                           \
          ghl_m1_compute_diagnostics(&params, &metric, &rad, &bad, &diagnostics), \
          ghl_error_m1_invalid_state, "invalid closure " #field " accepted");     \
    if(memcmp(&diagnostics, &sentinel, sizeof(diagnostics)) != 0)                 \
      fail_test("invalid closure metadata changed diagnostics");                  \
  } while(0)
  REJECT_CLOSURE_METADATA(xi, nextafter(0.0, -INFINITY));
  REJECT_CLOSURE_METADATA(xi, nextafter(1.0, INFINITY));
  REJECT_CLOSURE_METADATA(chi, NAN);
  REJECT_CLOSURE_METADATA(chi, nextafter(1.0 / 3.0, -INFINITY));
  REJECT_CLOSURE_METADATA(chi, nextafter(1.0, INFINITY));
  REJECT_CLOSURE_METADATA(root_residual, NAN);
  REJECT_CLOSURE_METADATA(root_residual, nextafter(0.0, -INFINITY));
  REJECT_CLOSURE_METADATA(root_iterations, -1);
#undef REJECT_CLOSURE_METADATA

  ghl_metric_quantities bad_metric = metric;
  bad_metric.gammaDD[0][0] = -1.0;
  require_error_code(
        ghl_m1_compute_diagnostics(&params, &bad_metric, &rad, &closure, &diagnostics),
        ghl_error_m1_invalid_metric, "diagnostics accepted non-SPD metric");
  ghl_m1_rad_state bad_state = rad;
  bad_state.E = 0.0;
  require_error_code(
        ghl_m1_compute_diagnostics(&params, &metric, &bad_state, &closure, &diagnostics),
        ghl_error_m1_invalid_state, "diagnostics accepted subfloor state");
  ghl_m1_closure bad_pressure = closure;
  bad_pressure.P[0][0] += 1.0;
  require_error_code(
        ghl_m1_compute_diagnostics(&params, &metric, &rad, &bad_pressure, &diagnostics),
        ghl_error_m1_invalid_state, "diagnostics accepted incorrect pressure trace");

  for(int component = 0; component < 3; ++component) {
    bad_state = rad;
    bad_state.F[component] = NAN;
    const ghl_m1_rad_state original = bad_state;
    require_error_code(
          ghl_m1_realizability_repair(&params, &metric, &bad_state),
          ghl_error_m1_invalid_state, "repair accepted nonfinite flux component");
    if(memcmp(&bad_state, &original, sizeof(bad_state)) != 0) {
      fail_test("failed repair changed state");
    }
  }

  /* Finite curvature with alpha=2 produces an unrepresentable energy source:
   * P:K=E*DBL_MAX by the pressure trace identity, before lapse multiplication.
   * This checks energy rejection separately from existing momentum cases. */
  metric.lapse = 2.0;
  metric.lapseinv = 0.5;
  metric.lapseinv2 = 0.25;
  prims.u0 = 0.5;
  require_error_code(
        ghl_m1_compute_closure_with_primitives(&params, &metric, &prims, &rad, &closure),
        ghl_success, "energy overflow reference closure failed");
  ghl_extrinsic_curvature curvature = { 0 };
  for(int i = 0; i < 3; ++i) {
    curvature.K[i][i] = DBL_MAX;
  }
  ghl_metric_quantities derivative = { 0 };
  const ghl_m1_sources sources_before = { .S_E = -17.0, .S = { -1.0, -2.0, -3.0 } };
  ghl_m1_sources sources = sources_before;
  require_error_code(
        ghl_m1_compute_geometry_sources(
              &params, &metric, &derivative, &derivative, &derivative, &curvature, &rad,
              &closure, &sources),
        ghl_error_m1_invalid_state, "geometry accepted overflowing energy source");
  if(memcmp(&sources, &sources_before, sizeof(sources)) != 0) {
    fail_test("energy source rejection changed outputs");
  }
}

static void check_comoving_arithmetic_rejections(void) {
  ghl_m1_parameters params;
  require_error_code(
        ghl_m1_initialize(0.1, 1.0e-12, 1.0, 1.0e-6, 1.0e-10, 100, 1.0e-10, &params),
        ghl_success, "comoving range initializer failed");
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  const ghl_m1_comoving sentinel = {
    .J = -17.0, .HU = { -1.0, -2.0, -3.0 }, .HD = { -4.0, -5.0, -6.0 }, .Hn = -19.0
  };

  /* At v=4/5, W=5/3. For isotropic P and F=0, J=91E/27
   * and H^x=-400E/81. E=DBL_MAX/2 overflows J; E=DBL_MAX/4
   * keeps J finite but overflows H^x. The third case has finite negative
   * J: transverse pressure and F^x=0.9 give J=(25/9)(1-1.44).
   * All input tensors pass the public pressure/realizability validation. */
  const double energies[] = { DBL_MAX / 2.0, DBL_MAX / 4.0, 1.0 };
  for(size_t k = 0; k < sizeof(energies) / sizeof(energies[0]); ++k) {
    ghl_m1_rad_state rad = { .E = energies[k] };
    ghl_m1_closure closure = { .chi = 1.0 / 3.0 };
    if(k == 2) {
      rad.F[0] = 0.9;
      closure.P[1][1] = rad.E;
    }
    else {
      for(int i = 0; i < 3; ++i) {
        closure.P[i][i] = rad.E / 3.0;
      }
    }
    ghl_m1_diagnostics diagnostics;
    require_error_code(
          ghl_m1_compute_diagnostics(&params, &metric, &rad, &closure, &diagnostics),
          ghl_success, "comoving range case did not pass public tensor validation");
    ghl_primitive_quantities prims = { .u0 = 1.0 };
    ghl_m1_comoving result = sentinel;
    require_error_code(
          ghl_m1_compute_comoving_moments(
                &params, &metric, &prims, &rad, &closure, &result),
          ghl_success, "stationary comoving range control failed");
    if(result.J != rad.E || result.Hn != 0.0) {
      fail_test("stationary comoving range control changed energy");
    }
    for(int i = 0; i < 3; ++i) {
      if(result.HU[i] != rad.F[i] || result.HD[i] != rad.F[i]) {
        fail_test("stationary comoving range control changed flux");
      }
    }
    prims.vU[0] = 0.8;
    result = sentinel;
    require_error_code(
          ghl_m1_compute_comoving_moments(
                &params, &metric, &prims, &rad, &closure, &result),
          ghl_error_m1_invalid_state, "invalid comoving arithmetic accepted");
    if(memcmp(&result, &sentinel, sizeof(result)) != 0) {
      fail_test("comoving arithmetic failure published partial moments");
    }

    if(k == 2) {
      prims.vU[0] = 1.0;
      require_error_code(
            ghl_m1_compute_comoving_moments(
                  &params, &metric, &prims, &rad, &closure, &result),
            ghl_error_u0_singular, "luminal comoving velocity accepted");
      if(memcmp(&result, &sentinel, sizeof(result)) != 0) {
        fail_test("comoving velocity failure published moments");
      }
    }
  }
}

static void check_scalar_boundary_rejections(void) {
  ghl_m1_parameters params;
  require_error_code(
        ghl_m1_initialize(0.1, 1.0e-12, 1.0, 1.0e-6, 1.0e-10, 100, 1.0e-10, &params),
        ghl_success, "scalar boundary initializer failed");
  const double invalid_floors[] = { NAN, 0.0 };
  for(size_t k = 0; k < sizeof(invalid_floors) / sizeof(invalid_floors[0]); ++k) {
    params.E_floor = invalid_floors[k];
    double energy = -17.0;
    bool applied = true;
    require_error_code(
          ghl_m1_apply_energy_floor(&params, 1.0, &energy, &applied),
          ghl_error_m1_invalid_E_floor, "invalid direct energy floor accepted");
    if(energy != -17.0 || !applied) {
      fail_test("invalid direct energy floor published outputs");
    }
  }
  double flux = -17.0;
  require_error_code(
        ghl_m1_compute_number_rusanov_flux(1.0, NAN, 0.0, 0.0, 1.0, &flux),
        ghl_error_m1_invalid_state, "nonfinite right number state accepted");
  if(flux != -17.0) {
    fail_test("invalid number Rusanov state published flux");
  }

  ghl_metric_quantities metric;
  ghl_initialize_metric(
        nextafter(0.0, 1.0), 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  /* The raw-speed API validates lapse and spatial geometry; its expression
   * reads no cached inverse lapse. Half the smallest positive double rounds
   * to zero, although the input lapse and spatial metric are admissible. */
  double delta;
  require_error_code(
        ghl_m1_compute_face_normal_delta_l(&metric, ghl_m1_dirn0, 1.0, &delta),
        ghl_success, "speed underflow metric rejected before arithmetic");
  if(delta != 2.0) {
    fail_test("speed underflow metric control mismatch");
  }
  double minus = -17.0, plus = -19.0;
  require_error_code(
        ghl_m1_compute_raw_lightcone_speeds(&metric, ghl_m1_dirn0, &minus, &plus),
        ghl_error_m1_invalid_metric, "underflowed lightcone scale accepted");
  if(minus != -17.0 || plus != -19.0) {
    fail_test("underflowed lightcone scale published speeds");
  }
}

static void check_closure_shift_overflow(void) {
  ghl_m1_parameters params;
  require_error_code(
        ghl_m1_initialize(0.1, 1.0e-12, 1.0, 1.0e-6, 1.0e-10, 100, 1.0e-10, &params),
        ghl_success, "closure shift initializer failed");
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  metric.betaU[0] = 4.0;
  /* Coordinate velocity -beta leaves Eulerian velocity exactly zero. */
  ghl_primitive_quantities prims = { .u0 = 1.0, .vU = { -4.0, 0.0, 0.0 } };
  ghl_m1_rad_state rad = { .E = 1.0, .F = { 0.5, 0.0, 0.0 } };
  ghl_m1_closure closure = { 0 };
  require_error_code(
        ghl_m1_compute_closure_with_primitives(&params, &metric, &prims, &rad, &closure),
        ghl_success, "finite closure shift control failed");
  const ghl_m1_closure sentinel = closure;
  rad.E = DBL_MAX;
  rad.F[0] = DBL_MAX / 2.0;
  ghl_m1_reset_closure_counters();
  /* F_mu n^mu=0 requires F_0=beta^i F_i=2*DBL_MAX. The
   * admissible spatial state therefore fails in workspace construction. */
  require_error_code(
        ghl_m1_compute_closure_with_primitives(&params, &metric, &prims, &rad, &closure),
        ghl_error_m1_invalid_state, "overflowing temporal flux accepted");
  if(!m1_closure_identical(&closure, &sentinel)) {
    fail_test("workspace failure published closure");
  }
  ghl_m1_closure_failure_stage_t stage;
  ghl_m1_get_last_closure_failure_stage(&stage);
  ghl_m1_closure_counters counters;
  ghl_m1_get_closure_counters(&counters);
  if(stage != ghl_m1_closure_failure_workspace || counters.invalid_state != 1
     || counters.ordinary_convergence || counters.endpoint_fallback
     || counters.iteration_exhaustion || counters.downstream_repair
     || counters.residual_rejection) {
    fail_test("workspace overflow diagnostic accounting mismatch");
  }
}

static void check_comoving_energy_failure(void) {
  ghl_m1_parameters params;
  /* DBL_MIN is accepted by the public epsilon_c contract, but
   * 1 - epsilon_c rounds to one.  Thus the unit-flux state remains accepted
   * by the realizability check and can reach the comoving-energy guard. */
  require_error_code(
        ghl_m1_initialize(
              DBL_MIN, DBL_MIN, 1.0e-8, 1.0e-6, 1.0e-12, 20, 1.0e-10, &params),
        ghl_success, "comoving-energy boundary initializer failed");

  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  ghl_primitive_quantities prims
        = { .u0 = 1.0, .vU = { nextafter(1.0, 0.0), 0.0, 0.0 } };
  const ghl_m1_rad_state rad = { .E = 1.0, .F = { 1.0, 0.0, 0.0 } };
  ghl_m1_closure closure = { .xi = -17.0, .chi = -19.0 };

  ghl_m1_reset_closure_counters();
  require_error_code(
        ghl_m1_compute_closure_with_primitives(&params, &metric, &prims, &rad, &closure),
        ghl_error_m1_invalid_state, "comoving-energy boundary was accepted");

  ghl_m1_closure_failure_stage_t stage;
  ghl_m1_get_last_closure_failure_stage(&stage);
  ghl_m1_closure_counters counters;
  ghl_m1_get_closure_counters(&counters);
  if(stage != ghl_m1_closure_failure_comoving_energy || counters.invalid_state != 1
     || counters.ordinary_convergence || counters.endpoint_fallback
     || counters.iteration_exhaustion || counters.downstream_repair
     || counters.residual_rejection) {
    fail_test("comoving-energy failure diagnostic accounting mismatch");
  }
}

int main(void) {
  check_success_returns();

  for(size_t index = 0; index < sizeof(m1_error_cases) / sizeof(m1_error_cases[0]);
      ++index) {
    check_fatal_error_case(&m1_error_cases[index]);
  }

  check_failed_initializer_stops();
  check_fixture_reader_and_comparator();
  check_parameter_validation();
  check_shared_null_and_state_contracts();
  check_shared_numeric_contracts();
  check_caller_parameter_matrix();
  check_caller_state_matrix();
  check_remaining_range_and_axis_failures();
  check_remaining_shared_rejections();
  check_comoving_arithmetic_rejections();
  check_scalar_boundary_rejections();
  check_closure_shift_overflow();
  check_comoving_energy_failure();

  ghl_info(
        "unit_test_m1_error_handling: fatal mappings, fixture rejection, and "
        "shared M1 validation checks passed\n");
  return 0;
}
