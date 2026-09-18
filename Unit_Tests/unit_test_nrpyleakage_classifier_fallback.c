#include <stdint.h>

#define GHL_NRPYLEAKAGE_FORCE_PORTABLE_CLASSIFIERS 1

#include "ghl_radiation.h"

int main(void) {
  const double finite_values[] = {
    0.0,
    -0.0,
    1.0,
    DBL_MAX,
    -DBL_MAX,
    ldexp(1.0, DBL_MIN_EXP - DBL_MANT_DIG),
    -ldexp(1.0, DBL_MIN_EXP - DBL_MANT_DIG),
  };
  for(size_t i = 0; i < sizeof(finite_values) / sizeof(finite_values[0]); i++) {
    if(!robust_isfinite(finite_values[i]) || robust_isnan(finite_values[i])) {
      ghl_error("Portable classifier rejected finite value %zu\n", i);
    }
  }

  const double infinite_values[] = { INFINITY, -INFINITY };
  for(size_t i = 0; i < sizeof(infinite_values) / sizeof(infinite_values[0]); i++) {
    if(robust_isfinite(infinite_values[i]) || robust_isnan(infinite_values[i])) {
      ghl_error("Portable classifier misclassified infinity %zu\n", i);
    }
  }

  const double nan_values[] = { NAN, -NAN };
  for(size_t i = 0; i < sizeof(nan_values) / sizeof(nan_values[0]); i++) {
    if(robust_isfinite(nan_values[i]) || !robust_isnan(nan_values[i])) {
      ghl_error("Portable classifier rejected NaN %zu\n", i);
    }
  }
  return 0;
}
