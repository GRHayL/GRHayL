#include "con2prim.h"

int Tabulated_Select_Method(
      const GRHayL_parameters *restrict params,
      const eos_parameters *restrict eos,
      const con2prim_method_t c2p_key,
      const metric_quantities *restrict metric,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics ) {

  switch( c2p_key ) {
    // Palenzuela1D routine (see https://arxiv.org/pdf/1712.07538.pdf)
    case Palenzuela1D:
      return Tabulated_Palenzuela1D_energy(params, eos, metric, cons, prims, diagnostics);
    case Palenzuela1D_entropy:
      return Tabulated_Palenzuela1D_entropy(params, eos, metric, cons, prims, diagnostics);
    case Newman1D:
      return Tabulated_Newman1D_energy(params, eos, metric, cons, prims, diagnostics);
    case Newman1D_entropy:
      return Tabulated_Newman1D_entropy(params, eos, metric, cons, prims, diagnostics);
    default:
      grhayl_Error(100, "Unknown c2p key in Tabulated_Select_Method (%d).", c2p_key);
      return 100;
  }
}

int Tabulated_Multi_Method(
      const GRHayL_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics ) {

  primitive_quantities prims_guess = *prims;
  if(params->calc_prim_guess)
    guess_primitives(eos, metric, prims, cons, &prims_guess);

  int check = Tabulated_Select_Method(params, eos, params->main_routine, metric, cons, prims, diagnostics);

  if( (check != 0) && (params->backup_routine[0] != None) ) {
    // Backup 1 triggered
    diagnostics->backup[0] = 1;
    // Reset guesses
    *prims = prims_guess;
    // Backup routine #1
    check = Tabulated_Select_Method(params, eos, params->backup_routine[0], metric, cons, prims, diagnostics);

    if( (check != 0) && (params->backup_routine[1] != None) ) {
      // Backup 2 triggered
      diagnostics->backup[1] = 1;
      // Reset guesses
      *prims = prims_guess;
      // Backup routine #2
      check = Tabulated_Select_Method(params, eos, params->backup_routine[1], metric, cons, prims, diagnostics);

      if( (check != 0) && (params->backup_routine[2] != None) ) {
        // Backup 3 triggered
        diagnostics->backup[2] = 1;
        // Reset guesses
        *prims = prims_guess;
        // Backup routine #3
        check = Tabulated_Select_Method(params, eos, params->backup_routine[2], metric, cons, prims, diagnostics);
      }
    }
  }
  return check;
}
