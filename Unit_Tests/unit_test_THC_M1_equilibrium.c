#include <stdio.h>

#include "ghl_unit_tests_et.h"

#define VALIDATION_THRESHOLD (1e-10)

int main(int argc, char **argv) {

  if(argc != 3) {
    printf("Usage: %s <eos table path> <input data file>\n", argv[0]);
    return 1;
  }

  const char *input_filename = argv[2];
  GHL_TEST_INITIALIZE(input_filename);

  int n_points = -1;
  GHL_TEST_LOG_READ(n_points);
  if(n_points < 0) {
    ghl_error("Failed to read n_points from file\n");
  }

  int n_neutrinos = -1;
  GHL_TEST_LOG_READ(n_neutrinos);
  if(n_neutrinos < 0) {
    ghl_error("Failed to read n_neutrinos from file\n");
  }

  for(int k = 0; k < n_points; k++) {
    for(int j = 0; j < n_points; j++) {
      for(int i = 0; i < n_points; i++) {
        GHL_TEST_LOG_READ_METRIC_AND_PRIMS(adm_metric, aux_metric, curv, prims);

        for(int n = 0; n < n_neutrinos; n++) {
          int ig = -1;
          GHL_TEST_LOG_READ(ig);

          // Compute nudens_0 and nudens_1
          double nudens_0 = 0, nudens_1 = 0;

          // Validate nudens_0 and nudens_1
          GHL_TEST_LOG_READ_AND_VALIDATE(nudens_0, VALIDATION_THRESHOLD);
          GHL_TEST_LOG_READ_AND_VALIDATE(nudens_1, VALIDATION_THRESHOLD);

          // Compute radiation stuff
          double rJ = 0, rE = 0, rnnu = 0, rN = 0;

          // Validate radiation stuff
          GHL_TEST_LOG_READ_AND_VALIDATE(rJ, VALIDATION_THRESHOLD);
          GHL_TEST_LOG_READ_AND_VALIDATE(rE, VALIDATION_THRESHOLD);
          GHL_TEST_LOG_READ_AND_VALIDATE(rnnu, VALIDATION_THRESHOLD);
          GHL_TEST_LOG_READ_AND_VALIDATE(rN, VALIDATION_THRESHOLD);

          // Read magic number, as a sanity check
          GHL_TEST_LOG_READ_AND_CHECK_MAGIC_NUMBER;
        }
      }
    }
  }

  GHL_TEST_FINALIZE(input_filename);
  return 0;
}
