#include "ghl_con2prim.h"

/*
 * Function     : ghl_con2prim_hybrid_select_method()
 * Description  : Calls the Con2Prim routine designated by c2p_key
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_con2prim_hybrid_select_method
*/
GHL_DEVICE
ghl_error_codes_t ghl_con2prim_hybrid_select_method(
      const ghl_con2prim_method_t c2p_key,
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  switch(c2p_key) {
    // Noble routines (see https://arxiv.org/pdf/astro-ph/0512420.pdf)
    case Noble2D:
      return ghl_hybrid_Noble2D(params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);
    case Noble1D:
      return ghl_hybrid_Noble1D(params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);
    // Font routine (see https://arxiv.org/abs/gr-qc/9811015)
    case Font1D:
      return ghl_hybrid_Font1D(params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);
    // Palenzuela1D routine (see https://arxiv.org/pdf/1712.07538.pdf)
    case Palenzuela1D:
      return ghl_hybrid_Palenzuela1D_energy(params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);
    // Entropy routines (see https://arxiv.org/abs/2208.14487)
    case Noble1D_entropy:
      return ghl_hybrid_Noble1D_entropy(params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);
    //case Noble1D_entropy2:
    //  return ghl_hybrid_Noble1D_entropy2(params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);
    case Palenzuela1D_entropy:
      return ghl_hybrid_Palenzuela1D_entropy(params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);
    default:
      return ghl_error_invalid_c2p_key;
  }
}

/*
 * Function     : ghl_con2prim_tabulated_select_method()
 * Description  : Calls the Con2Prim routine designated by c2p_key
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_con2prim_tabulated_select_method
*/
GHL_DEVICE
ghl_error_codes_t ghl_con2prim_tabulated_select_method(
      const ghl_con2prim_method_t c2p_key,
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  switch(c2p_key) {
    // Palenzuela1D routine (see https://arxiv.org/pdf/1712.07538.pdf)
    case Palenzuela1D:
      return ghl_tabulated_Palenzuela1D_energy(params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);
    // Newman 1D routine (see https://escholarship.org/uc/item/0s53f84b)
    case Newman1D:
      return ghl_tabulated_Newman1D_energy(params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);
    // Entropy routines (see https://arxiv.org/abs/2208.14487)
    case Palenzuela1D_entropy:
      return ghl_tabulated_Palenzuela1D_entropy(params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);
    case Newman1D_entropy:
      return ghl_tabulated_Newman1D_entropy(params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);
    default:
      return ghl_error_invalid_c2p_key;
  }
}

/*
 * Function     : ghl_con2prim_tabulated_select_method()
 * Description  : Calls the Con2Prim routine designated by c2p_key
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_con2prim_tabulated_select_method
*/
GHL_DEVICE
ghl_error_codes_t ghl_con2prim_hybrid_multi_method(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  if(params->calc_prim_guess)
    ghl_guess_primitives(eos, ADM_metric, cons_undens, prims);

  // Store primitive guesses (used if con2prim fails)
  const ghl_primitive_quantities prims_guess = *prims;

  ghl_error_codes_t error = ghl_con2prim_hybrid_select_method(params->main_routine, params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);

  if(error && params->backup_routine[0] != None) {
    // Backup 1 triggered
    diagnostics->backup[0] = true;
    // Reset guesses
    *prims = prims_guess;
    // Backup routine #1
    error = ghl_con2prim_hybrid_select_method(params->backup_routine[0], params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);

    if(error && params->backup_routine[1] != None) {
      // Backup 2 triggered
      diagnostics->backup[1] = true;
      // Reset guesses
      *prims = prims_guess;
      // Backup routine #2
      error = ghl_con2prim_hybrid_select_method(params->backup_routine[1], params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);

      if(error && params->backup_routine[2] != None) {
        // Backup 3 triggered
        diagnostics->backup[2] = true;
        // Reset guesses
        *prims = prims_guess;
        // Backup routine #3
        error = ghl_con2prim_hybrid_select_method(params->backup_routine[2], params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);
      }
    }
  }
  return error;
}

/*
 * Function     : ghl_con2prim_tabulated_select_method()
 * Description  : Calls the Con2Prim routine designated by c2p_key
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_con2prim_tabulated_select_method
*/
GHL_DEVICE
ghl_error_codes_t ghl_con2prim_tabulated_multi_method(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  if(params->calc_prim_guess)
    ghl_guess_primitives(eos, ADM_metric, cons_undens, prims);

  // Store primitive guesses (used if con2prim fails)
  const ghl_primitive_quantities prims_guess = *prims;

  ghl_error_codes_t error = ghl_con2prim_tabulated_select_method(params->main_routine, params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);

  if(error && params->backup_routine[0] != None) {
    // Backup 1 triggered
    diagnostics->backup[0] = true;
    // Reset guesses
    *prims = prims_guess;
    // Backup routine #1
    error = ghl_con2prim_tabulated_select_method(params->backup_routine[0], params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);

    if(error && params->backup_routine[1] != None) {
      // Backup 2 triggered
      diagnostics->backup[1] = true;
      // Reset guesses
      *prims = prims_guess;
      // Backup routine #2
      error = ghl_con2prim_tabulated_select_method(params->backup_routine[1], params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);

      if(error && params->backup_routine[2] != None) {
        // Backup 3 triggered
        diagnostics->backup[2] = true;
        // Reset guesses
        *prims = prims_guess;
        // Backup routine #3
        error = ghl_con2prim_tabulated_select_method(params->backup_routine[2], params, eos, ADM_metric, metric_aux, cons_undens, prims, diagnostics);
      }
    }
  }
  return error;
}
