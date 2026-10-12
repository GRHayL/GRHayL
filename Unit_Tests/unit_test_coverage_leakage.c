#include "ghl_unit_tests.h"
#include "ghl_nrpyleakage_nucleon_blocking.h"
#include "ghl_nrpyleakage_rate_helpers.h"

/* Error returns must not turn invalid inputs into published beta moments. */
static void check_beta_input_rejection(void) {
  nrpyl_beta_moments emitted = { -7.0, -9.0 };
  nrpyl_beta_moments absorbed = emitted;
  if(nrpyl_compute_beta_emission_moments(INFINITY, 0.0, 0.0, 1, 0.0, &emitted)
           != ghl_error_nrpyleakage_blocking
     || nrpyl_compute_beta_absorption_moments(INFINITY, 0.0, 0.0, -1, 0.0, &absorbed)
              != ghl_error_nrpyleakage_blocking
     || emitted.number != -7.0 || emitted.energy != -9.0
     || absorbed.number != -7.0 || absorbed.energy != -9.0) {
    ghl_error("Infinite beta temperature was not rejected transactionally\n");
  }
  const nrpyl_beta_moments invalid[] = {
    { -1.0, 1.0 }, { 1.0, -1.0 }, { 1.0, INFINITY }, { NAN, 1.0 }
  };
  for(size_t c = 0; c < sizeof(invalid) / sizeof(invalid[0]); ++c) {
    if(nrpyl_validate_beta_moments(&invalid[c]) != ghl_error_nrpyleakage_blocking) {
      ghl_error("Invalid beta moment %zu accepted\n", c);
    }
  }
}

static void check_population_rejection(void) {
  nrpyl_nucleon_population population = { 0.0 };
  if(nrpyl_compute_population(NAN, &population) != ghl_error_nrpyleakage_blocking
     || nrpyl_compute_population(-INFINITY, &population)
              != ghl_error_nrpyleakage_blocking) {
    ghl_error("Nonfinite logarithmic population accepted\n");
  }
  const double cases[][4] = {
    { INFINITY, 0.5, 0.1, 0.2 }, { 0.5, INFINITY, 0.1, 0.2 },
    { 0.5, 0.5, INFINITY, 0.2 }, { 0.5, 0.5, 0.1, INFINITY },
    { -0.1, 0.5, 0.1, 0.2 }, { 0.8, 0.5, 0.1, 0.2 },
    { 0.5, -0.1, 0.1, 0.2 }, { 0.5, 0.8, 0.1, 0.2 },
    { 0.5, 0.5, -0.1, 0.2 }, { 0.5, 0.5, 0.1, -0.2 }
  };
  for(size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); ++c) {
    double forward = cases[c][2], reverse = cases[c][3];
    if(nrpyl_validate_blocking_outputs(
             0.7, 0.6, cases[c][0], cases[c][1], &forward, &reverse)
             != ghl_error_nrpyleakage_blocking
       || forward != cases[c][2] || reverse != cases[c][3]) {
      ghl_error("Invalid blocking output %zu was accepted or modified\n", c);
    }
  }
}

int main(void) {
  check_beta_input_rejection();
  check_population_rejection();
  ghl_info("Leakage coverage boundary tests passed.\n");
  return 0;
}
