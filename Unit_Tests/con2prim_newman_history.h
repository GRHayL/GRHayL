#ifndef CON2PRIM_NEWMAN_HISTORY_H_
#define CON2PRIM_NEWMAN_HISTORY_H_

#include "ghl_con2prim.h"

static unsigned Newman_history_calls;

static ghl_error_codes_t Newman_history_pressure(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double eps_or_entropy,
      double *restrict press,
      double *restrict temperature) {

  (void)eos;
  (void)rho;
  (void)Y_e;
  (void)eps_or_entropy;
  ++Newman_history_calls;
  // Each attempt starts at zero: exact integer increments keep AtR == 1,
  // preventing Aitken resets and convergence through the final iteration.
  *press += 1.0;
  *temperature = 1.0;
  return ghl_success;
}

// Call after initializing tabulated EOS function pointers; no table is needed.
// This reaches AtP[maxsteps]; detecting an undersized history array requires
// ASan/UBSan. The status and call-count checks alone do not detect that overflow.
static void test_Newman_history_bounds(void) {
  const ghl_parameters params = {
    .max_Lorentz_factor = 2.0,
    .inv_sq_max_Lorentz_factor = 0.25,
  };
  const ghl_eos_parameters eos = {
    .eos_type = ghl_eos_tabulated,
    .rho_min = 0.5,
    .rho_max = 2.0,
    .Y_e_min = 0.0,
    .Y_e_max = 1.0,
    .T_min = 1.0,
    .T_max = 2.0,
    .eps_min = 0.0,
    .eps_max = 2.0,
    .entropy_min = 0.0,
    .entropy_max = 2.0,
  };
  const ghl_conservative_quantities cons = {
    .rho = 1.0,
    .tau = 1.0,
    .Y_e = 0.1,
    .entropy = 1.0,
  };
  ghl_metric_quantities metric;
  ghl_initialize_metric(1, 0, 0, 0, 1, 0, 0, 1, 0, 1, &metric);
  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(&metric, &metric_aux);

  ghl_error_codes_t (*saved_energy)(
        const ghl_eos_parameters *, double, double, double, double *, double *)
        = ghl_tabulated_compute_P_T_from_eps;
  ghl_error_codes_t (*saved_entropy)(
        const ghl_eos_parameters *, double, double, double, double *, double *)
        = ghl_tabulated_compute_P_T_from_S;
  ghl_tabulated_compute_P_T_from_eps = Newman_history_pressure;
  ghl_tabulated_compute_P_T_from_S = Newman_history_pressure;

  ghl_error_codes_t errors[2];
  ghl_con2prim_diagnostics diagnostics[2];
  unsigned calls[2];
  for(int variant = 0; variant < 2; ++variant) {
    ghl_primitive_quantities prims = { .temperature = 1.5 };
    ghl_initialize_diagnostics(&diagnostics[variant]);
    Newman_history_calls = 0;
    if(variant == 0) {
      errors[variant] = ghl_tabulated_Newman1D_energy(
            &params, &eos, &metric, &metric_aux, &cons, &prims, &diagnostics[variant]);
    }
    else {
      errors[variant] = ghl_tabulated_Newman1D_entropy(
            &params, &eos, &metric, &metric_aux, &cons, &prims, &diagnostics[variant]);
    }
    calls[variant] = Newman_history_calls;
  }

  ghl_tabulated_compute_P_T_from_eps = saved_energy;
  ghl_tabulated_compute_P_T_from_S = saved_entropy;
  for(int variant = 0; variant < 2; ++variant) {
    // maxsteps=300 in each Newman helper, plus one T_min retry: 2*300 calls.
    // This checks that attempt schedule; failure-path n_iter is unspecified.
    if(errors[variant] != ghl_error_c2p_max_iter || calls[variant] != 600) {
      ghl_error(
            "Newman %s history: expected max-iteration error and "
            "600 EOS calls; got error=%d and %u calls\n",
            variant == 0 ? "energy" : "entropy", (int)errors[variant], calls[variant]);
    }
  }
}

#endif // CON2PRIM_NEWMAN_HISTORY_H_
