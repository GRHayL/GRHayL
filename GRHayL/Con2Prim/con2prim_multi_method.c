#include "con2prim.h"

int ghl_con2prim_select_method(
      const con2prim_method_t c2p_key,
      const ghl_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics ) {

  switch( c2p_key ) {
    // Noble2D routine (see https://arxiv.org/pdf/astro-ph/0512420.pdf)
    case Noble2D:
      return ghl_hybrid_Noble2D(params, eos, ADM_metric, metric_aux, cons, prims, diagnostics);
    // "Font fix" routine (see )
    case FontFix:
      return ghl_hybrid_Font_fix(params, eos, ADM_metric, metric_aux, cons, prims, diagnostics);
    // Palenzuela1D routine (see https://arxiv.org/pdf/1712.07538.pdf)
    case Palenzuela1D:
      return ghl_tabulated_Palenzuela1D_energy(params, eos, ADM_metric, metric_aux, cons, prims, diagnostics);
    // Palenzuela1D routine with entropy (see https://arxiv.org/pdf/2208.14487.pdf)
    case Palenzuela1D_entropy:
      return ghl_tabulated_Palenzuela1D_entropy(params, eos, ADM_metric, metric_aux, cons, prims, diagnostics);
    // Newman 1D routine (see https://escholarship.org/content/qt0s53f84b/qt0s53f84b.pdf)
    case Newman1D:
      return ghl_tabulated_Newman1D_energy(params, eos, ADM_metric, metric_aux, cons, prims, diagnostics);
    // Newman1D routine with entropy (see https://arxiv.org/pdf/2208.14487.pdf)
    case Newman1D_entropy:
      return ghl_tabulated_Newman1D_entropy(params, eos, ADM_metric, metric_aux, cons, prims, diagnostics);
    default:
      ghl_Error(100, "Unknown c2p key (%d).", c2p_key);
      return 100;
  }
}

int ghl_con2prim_multi_method(
      const ghl_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics ) {

  if(params->calc_prim_guess)
    ghl_guess_primitives(eos, ADM_metric, metric_aux, cons, prims);

  // Store primitive guesses (used if con2prim fails)
  const primitive_quantities prims_guess = *prims;

  int failed = ghl_con2prim_select_method(params->main_routine, params, eos, ADM_metric, metric_aux, cons, prims, diagnostics);

  if( failed && params->backup_routine[0] != None ) {
    // Backup 1 triggered
    diagnostics->backup[0] = 1;
    // Reset guesses
    *prims = prims_guess;
    // Backup routine #1
    failed = ghl_con2prim_select_method(params->backup_routine[0], params, eos, ADM_metric, metric_aux, cons, prims, diagnostics);

    if( failed && params->backup_routine[1] != None ) {
      // Backup 2 triggered
      diagnostics->backup[1] = 1;
      // Reset guesses
      *prims = prims_guess;
      // Backup routine #2
      failed = ghl_con2prim_select_method(params->backup_routine[1], params, eos, ADM_metric, metric_aux, cons, prims, diagnostics);

      if( failed && params->backup_routine[2] != None ) {
        // Backup 3 triggered
        diagnostics->backup[2] = 1;
        // Reset guesses
        *prims = prims_guess;
        // Backup routine #3
        failed = ghl_con2prim_select_method(params->backup_routine[2], params, eos, ADM_metric, metric_aux, cons, prims, diagnostics);
      }
    }
  }
  return failed;
}
