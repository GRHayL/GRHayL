#include "ghl_flux_source.h"

#include <math.h>

ghl_error_codes_t ghl_calculate_Rusanov_flux(
      const double *restrict state_L,
      const double *restrict state_R,
      const double *restrict physical_flux_L,
      const double *restrict physical_flux_R,
      const int component_count,
      const double speed,
      double *restrict flux) {
  if(state_L == NULL || state_R == NULL || physical_flux_L == NULL ||
     physical_flux_R == NULL || flux == NULL || component_count <= 0 ||
     !isfinite(speed) || speed < 0.0)
    return ghl_error_flux_source_invalid_input;

  /* Validate and evaluate every candidate before publishing any component. */
  for(int component = 0; component < component_count; ++component) {
    if(!isfinite(state_L[component]) || !isfinite(state_R[component]) ||
       !isfinite(physical_flux_L[component]) ||
       !isfinite(physical_flux_R[component]))
      return ghl_error_flux_source_invalid_input;

    const double candidate =
        0.5 * (physical_flux_L[component] + physical_flux_R[component])
        - 0.5 * speed * (state_R[component] - state_L[component]);
    if(!isfinite(candidate))
      return ghl_error_flux_source_invalid_input;
  }

  for(int component = 0; component < component_count; ++component) {
    flux[component] =
        0.5 * (physical_flux_L[component] + physical_flux_R[component])
        - 0.5 * speed * (state_R[component] - state_L[component]);
  }
  return ghl_success;
}
