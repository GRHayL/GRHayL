#include "ghl_test_helpers.h"
#include "ghl.h"

#include <limits.h>
#include <stdint.h>

size_t ghl_test_read_grid_size(
      FILE *restrict infile,
      const char *restrict filename,
      const int minimum,
      int *restrict dirlength) {
  if(fread(dirlength, sizeof(*dirlength), 1, infile) != 1) {
    ghl_error("Could not read the grid size from %s.\n", filename);
  }
  if(*dirlength < minimum) {
    ghl_error(
          "%s requires a grid size of at least %d; got %d.\n", filename, minimum,
          *dirlength);
  }

  const size_t n = (size_t)*dirlength;
  if(n > SIZE_MAX / n || n * n > SIZE_MAX / n) {
    ghl_error("The grid size in %s is too large.\n", filename);
  }
  const size_t arraylength = n * n * n;
  if(arraylength > INT_MAX || arraylength > SIZE_MAX / sizeof(double)
     || arraylength > SIZE_MAX / 25) {
    ghl_error("The grid size in %s is too large.\n", filename);
  }
  return arraylength;
}
