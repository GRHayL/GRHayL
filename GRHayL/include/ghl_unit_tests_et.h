#ifndef GHL_UNIT_TESTS_ET_H
#define GHL_UNIT_TESTS_ET_H

#include "ghl_unit_tests.h"

#include <stdio.h>
#include <stdlib.h>

#define GHL_TEST_MAGIC_NUMBER (4224)

#define GHL_PRINTF(...)  \
  printf("(GHL TEST) "); \
  printf(__VA_ARGS__);

// Initialize unit testing for the Einstein Toolkit.
#define GHL_TEST_INITIALIZE(ghl_test_filename)                \
  FILE *ghl_test_fp_ = fopen(ghl_test_filename, "rb");        \
  if(ghl_test_fp_ == NULL) {                                  \
    ghl_error("Failed to open file %s\n", ghl_test_filename); \
  }                                                           \
  GHL_PRINTF("Succesfully opened file %s\n", ghl_test_filename);

// Close log file. This must appear at the end of every tested function.
#define GHL_TEST_FINALIZE(ghl_test_filename) \
  fclose(ghl_test_fp_);                      \
  GHL_PRINTF("Successfully closed file %s\n", ghl_test_filename);

// Generic macro to read "anything" from file. Uses typeof instead of decltype.
#define GHL_TEST_LOG_READ(var)                                     \
  if(fread(&var, sizeof(__typeof__(var)), 1, ghl_test_fp_) != 1) { \
    GHL_PRINTF("Failed to read '%s' from log file\n", #var);       \
    exit(1);                                                       \
  }

// Read a magic number and perform sanity check
#define GHL_TEST_LOG_READ_AND_CHECK_MAGIC_NUMBER          \
  {                                                       \
    int ghl_test_magic_number_ = 0;                       \
    GHL_TEST_LOG_READ(ghl_test_magic_number_);            \
    if(ghl_test_magic_number_ != GHL_TEST_MAGIC_NUMBER) { \
      GHL_PRINTF("Magic number sanity check failed\n");   \
      exit(2);                                            \
    }                                                     \
  }

// Read metric and primitives from log file
#define GHL_TEST_LOG_READ_METRIC_AND_PRIMS(adm_, aux_, curv_, prims_) \
  GHL_TEST_LOG_READ(adm_);                                            \
  GHL_TEST_LOG_READ(aux_);                                            \
  GHL_TEST_LOG_READ(curv_);                                           \
  GHL_TEST_LOG_READ(prims_);

// Read variable from file with same type as result. Check if relative error is smaller
// than threshold.
#define GHL_TEST_LOG_READ_AND_VALIDATE(result, threshold)                        \
  {                                                                              \
    __typeof__(result) trusted;                                                  \
    GHL_TEST_LOG_READ(trusted);                                                  \
    if(relative_error(result, trusted) > threshold) {                            \
      GHL_PRINTF(                                                                \
            "Relative error threshold (>%g) -- %s: %g %g\n", threshold, #result, \
            result, trusted);                                                    \
    }                                                                            \
  }

#endif // GHL_UNIT_TESTS_ET_H
