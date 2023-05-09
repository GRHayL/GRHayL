#include "con2prim.h"

int grhayl_con2prim_select_method(
      const con2prim_method_t c2p_key,
      const GRHayL_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics ) {

  switch( c2p_key ) {
    // Noble2D routine (see https://arxiv.org/pdf/astro-ph/0512420.pdf)
    case Noble2D:
      return Hybrid_Noble2D(params, eos, metric, cons, prims, diagnostics);
    // Palenzuela1D routine (see https://arxiv.org/pdf/1712.07538.pdf)
    case Palenzuela1D:
      return Tabulated_Palenzuela1D_energy(params, eos, metric, cons, prims, diagnostics);
    // Palenzuela1D routine with entropy (see https://arxiv.org/pdf/2208.14487.pdf)
    case Palenzuela1D_entropy:
      return Tabulated_Palenzuela1D_entropy(params, eos, metric, cons, prims, diagnostics);
    // Newman 1D routine (see https://escholarship.org/content/qt0s53f84b/qt0s53f84b.pdf)
    case Newman1D:
      return Tabulated_Newman1D_energy(params, eos, metric, cons, prims, diagnostics);
    // Newman1D routine with entropy (see https://arxiv.org/pdf/2208.14487.pdf)
    case Newman1D_entropy:
      return Tabulated_Newman1D_entropy(params, eos, metric, cons, prims, diagnostics);
    default:
      grhayl_Error(100, "Unknown c2p key (%d).", c2p_key);
      return 100;
  }
}

int grhayl_con2prim_multi_method(
      const GRHayL_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics ) {

  if(params->calc_prim_guess)
    guess_primitives(eos, metric, cons, prims);

  // Store primitive guesses (used if con2prim fails)
  const primitive_quantities prims_guess = *prims;

  int failed = grhayl_con2prim_select_method(params->main_routine, params, eos, metric, cons, prims, diagnostics);

  if( failed && params->backup_routine[0] != None ) {
    // Backup 1 triggered
    diagnostics->backup[0] = 1;
    // Reset guesses
    *prims = prims_guess;
    // Backup routine #1
    failed = grhayl_con2prim_select_method(params->backup_routine[0], params, eos, metric, cons, prims, diagnostics);

    if( failed && params->backup_routine[1] != None ) {
      // Backup 2 triggered
      diagnostics->backup[1] = 1;
      // Reset guesses
      *prims = prims_guess;
      // Backup routine #2
      failed = grhayl_con2prim_select_method(params->backup_routine[1], params, eos, metric, cons, prims, diagnostics);

      if( failed && params->backup_routine[2] != None ) {
        // Backup 3 triggered
        diagnostics->backup[2] = 1;
        // Reset guesses
        *prims = prims_guess;
        // Backup routine #3
        failed = grhayl_con2prim_select_method(params->backup_routine[2], params, eos, metric, cons, prims, diagnostics);
      }
    }
  }
  return failed;
}
