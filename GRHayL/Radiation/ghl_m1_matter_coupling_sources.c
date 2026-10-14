#include "ghl_m1.h"
#include "ghl_m1_utils.h"

ghl_error_codes_t ghl_m1_compute_matter_coupling_sources(
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_sources *restrict interaction_sources,
      double *restrict source_tildetau,
      double source_tildeS[3]) {
  if(metric == NULL || interaction_sources == NULL || source_tildetau == NULL
     || source_tildeS == NULL) {
    return ghl_error_m1_null_pointer;
  }

  if(!isfinite(metric->lapse) || metric->lapse <= 0.0 || !isfinite(metric->sqrt_detgamma)
     || metric->sqrt_detgamma <= 0.0) {
    return ghl_error_m1_invalid_metric;
  }
  const double alpha_sqrt_detgamma = metric->lapse * metric->sqrt_detgamma;

  if(!isfinite(interaction_sources->S_E)) {
    return ghl_error_m1_invalid_state;
  }

  const double source_tildetau_candidate
        = -alpha_sqrt_detgamma * interaction_sources->S_E;
  if(!isfinite(source_tildetau_candidate)) {
    return ghl_error_m1_invalid_state;
  }

  double source_tildeS_candidate[3];
  for(int i = 0; i < 3; i++) {
    if(!isfinite(interaction_sources->S[i])) {
      return ghl_error_m1_invalid_state;
    }
    source_tildeS_candidate[i] = -alpha_sqrt_detgamma * interaction_sources->S[i];
    if(!isfinite(source_tildeS_candidate[i])) {
      return ghl_error_m1_invalid_state;
    }
  }

  *source_tildetau = source_tildetau_candidate;
  for(int i = 0; i < 3; i++) {
    source_tildeS[i] = source_tildeS_candidate[i];
  }

  return ghl_success;
}
