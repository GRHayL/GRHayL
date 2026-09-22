#ifndef GHL_TEST_HELPERS_H_
#define GHL_TEST_HELPERS_H_

#include <stddef.h>
#include <stdio.h>

size_t ghl_test_read_grid_size(
      FILE *restrict infile,
      const char *restrict filename,
      int minimum,
      int *restrict dirlength);

#endif // GHL_TEST_HELPERS_H_
