#include <stdint.h>

#undef UINT64_C

#include "ghl_radiation.h"

int main(void) {
  if(!robust_isfinite(1.0) || robust_isfinite(INFINITY) || !robust_isnan(NAN)) {
    ghl_error("NRPyLeakage portable classifier fallback failed\n");
  }
  return 0;
}
