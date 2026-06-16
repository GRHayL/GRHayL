#include "ghl_con2prim.h"

/*
 * Function     : ghl_con2prim_hybrid_select_method()
 * Description  : Calls the Con2Prim routine designated by c2p_key
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_con2prim_hybrid_select_method
*/

ghl_error_codes_t ghl_con2prim_hybrid_select_method(
      const ghl_con2prim_id_t c2p_key,
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  switch(c2p_key) {
    // Noble routines (see https://arxiv.org/pdf/astro-ph/0512420.pdf)
    case ghl_con2prim_id_Noble2D:
      return ghl_hybrid_Noble2D(params, eos, metric_adm, metric_aux, cons_undens, prims, diagnostics);
    case ghl_con2prim_id_Noble1D:
      return ghl_hybrid_Noble1D(params, eos, metric_adm, metric_aux, cons_undens, prims, diagnostics); 
    // Font routine (see https://arxiv.org/abs/gr-qc/9811015)
    case ghl_con2prim_id_Font1D:
      return ghl_hybrid_Font1D(params, eos, metric_adm, metric_aux, cons_undens, prims, diagnostics);
    // Palenzuela1D routine (see https://arxiv.org/pdf/1712.07538.pdf)
    case ghl_con2prim_id_Palenzuela1D:
      return ghl_hybrid_Palenzuela1D_energy(params, eos, metric_adm, metric_aux, cons_undens, prims, diagnostics);
    // Entropy routines (see https://arxiv.org/abs/2208.14487)
    case ghl_con2prim_id_Noble1D_entropy:
      return ghl_hybrid_Noble1D_entropy(params, eos, metric_adm, metric_aux, cons_undens, prims, diagnostics); 
    case ghl_con2prim_id_Palenzuela1D_entropy:
      return ghl_hybrid_Palenzuela1D_entropy(params, eos, metric_adm, metric_aux, cons_undens, prims, diagnostics);
    default:
      return ghl_error_invalid_c2p_key;
  }
}

/*
 * Function     : ghl_con2prim_tabulated_select_method()
 * Description  : Calls the Con2Prim routine designated by c2p_key
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_con2prim_tabulated_select_method
*/

ghl_error_codes_t ghl_con2prim_tabulated_select_method(
      const ghl_con2prim_id_t c2p_key,
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

#ifdef GHL_DISABLE_HDF5
  return ghl_error_used_disabled_hdf5;
#else
  switch(c2p_key) {
    // Noble routine (see https://arxiv.org/pdf/astro-ph/0512420.pdf)
    case ghl_con2prim_id_Noble2D:
      return ghl_tabulated_Noble2D(params, eos, metric_adm, metric_aux, cons_undens, prims, diagnostics);
    // Palenzuela1D routine (see https://arxiv.org/pdf/1712.07538.pdf)
    case ghl_con2prim_id_Palenzuela1D:
      return ghl_tabulated_Palenzuela1D_energy(params, eos, metric_adm, metric_aux, cons_undens, prims, diagnostics);
    // Newman 1D routine (see https://escholarship.org/uc/item/0s53f84b)
    case ghl_con2prim_id_Newman1D:
      return ghl_tabulated_Newman1D_energy(params, eos, metric_adm, metric_aux, cons_undens, prims, diagnostics);
    // Entropy routines (see https://arxiv.org/abs/2208.14487)
    case ghl_con2prim_id_Palenzuela1D_entropy:
      return ghl_tabulated_Palenzuela1D_entropy(params, eos, metric_adm, metric_aux, cons_undens, prims, diagnostics);
    case ghl_con2prim_id_Newman1D_entropy:
      return ghl_tabulated_Newman1D_entropy(params, eos, metric_adm, metric_aux, cons_undens, prims, diagnostics);
    default:
      return ghl_error_invalid_c2p_key;
  }
#endif
}

/*
 * Function     : ghl_con2prim_tabulated_select_method()
 * Description  : Calls the Con2Prim routine designated by c2p_key
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_con2prim_tabulated_select_method
*/

ghl_error_codes_t ghl_con2prim_hybrid_multi_method(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  if(params->calc_prim_guess)
    ghl_guess_primitives(eos, metric_adm, cons_undens, prims);

  // Store primitive guesses (used if con2prim fails)
  const ghl_primitive_quantities prims_guess = *prims;

  ghl_error_codes_t error;
  error = ghl_con2prim_hybrid_select_method(params->main_routine,
                                            params, eos, metric_adm, metric_aux,
                                            cons_undens, prims, diagnostics);

  // Note(Leo): this updated backup strategy works for any number of backup
  //            routines, cleaning up the logic, removing duplicated code, and
  //            minimizing the chances of a bug. The variable `n_backups` can
  //            be set from `params`, but since this requires an API change
  //            I have not implemented it yet. This comment and the variable
  //            should be removed in a future PR.
  const size_t n_backups = 3;
  for(size_t n = 0; (n < n_backups) && (error != ghl_success); n++) {
    // FIXME(Leo): once the comment above is addressed, this if statement can
    //             be removed, as the backup routines will be set at startup
    //             and we won't need to check if they are None anymore.
    if(params->backup_routine[n] == ghl_con2prim_id_None) {
      break;
    }
    // Backup triggered
    diagnostics->backup[n] = true;
    // Reset guesses
    *prims = prims_guess;
    // Backup routine
    error = ghl_con2prim_hybrid_select_method(params->backup_routine[n],
                                              params, eos, metric_adm, metric_aux,
                                              cons_undens, prims, diagnostics);
  }
  return error;
}

/*
 * Function     : ghl_con2prim_tabulated_select_method()
 * Description  : Calls the Con2Prim routine designated by c2p_key
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_con2prim_tabulated_select_method
*/

ghl_error_codes_t ghl_con2prim_tabulated_multi_method(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

#ifdef GHL_DISABLE_HDF5
  return ghl_error_used_disabled_hdf5;
#else
  if(params->calc_prim_guess)
    ghl_guess_primitives(eos, metric_adm, cons_undens, prims);

  // Store primitive guesses (used if con2prim fails)
  const ghl_primitive_quantities prims_guess = *prims;

  ghl_error_codes_t error;
  error = ghl_con2prim_tabulated_select_method(params->main_routine,
                                               params, eos, metric_adm, metric_aux,
                                               cons_undens, prims, diagnostics);

  // Note(Leo): this updated backup strategy works for any number of backup
  //            routines, cleaning up the logic, removing duplicated code, and
  //            minimizing the chances of a bug. The variable `n_backups` can
  //            be set from `params`, but since this requires an API change
  //            I have not implemented it yet. This comment and the variable
  //            should be removed in a future PR.
  const size_t n_backups = 3;
  for(size_t n = 0; (n < n_backups) && (error != ghl_success); n++) {
    // FIXME(Leo): once the comment above is addressed, this if statement can
    //             be removed, as the backup routines will be set at startup
    //             and we won't need to check if they are None anymore.
    if(params->backup_routine[n] == ghl_con2prim_id_None) {
      break;
    }
    // Backup triggered
    diagnostics->backup[n] = true;
    // Reset guesses
    *prims = prims_guess;
    // Backup routine
    error = ghl_con2prim_tabulated_select_method(params->backup_routine[n],
                                                 params, eos, metric_adm, metric_aux,
                                                 cons_undens, prims, diagnostics);
  }
  return error;
#endif
}
