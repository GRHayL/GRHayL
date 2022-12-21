#include <stdio.h>
#include <stdlib.h>
#include "con2prim.h"

int Hybrid_Select_Method(const GRHayL_parameters *restrict params,
                         const eos_parameters *restrict eos,
                         const int c2p_key,
                         const metric_quantities *restrict metric,
                         const conservative_quantities *restrict cons,
                         primitive_quantities *restrict prims,
                         con2prim_diagnostics *restrict diagnostics );

int Hybrid_Multi_Method( const GRHayL_parameters *restrict params,
                         const eos_parameters *restrict eos,
                         const metric_quantities *restrict metric,
                         const conservative_quantities *restrict cons,
                         const primitive_quantities *restrict prims,
                         primitive_quantities *restrict prims_guess,
                         con2prim_diagnostics *restrict diagnostics ) {

  // Compute primitive initial guesses if they are not provided
  if(params->calc_prim_guess) guess_primitives(eos, metric, prims, cons, prims_guess);
  int check = Hybrid_Select_Method(params, eos, params->main_routine, metric, cons, prims_guess, diagnostics);

  if( (check != 0) && (params->backup_routine[0] != None) ) {
    // Backup 1 triggered
    diagnostics->backup[0] = 1;
    // Recompute guesses
    if(params->calc_prim_guess) guess_primitives(eos, metric, prims, cons, prims_guess);
    // Backup routine #1
    check = Hybrid_Select_Method(params, eos, params->backup_routine[0], metric, cons, prims_guess, diagnostics);

    if( (check != 0) && (params->backup_routine[1] != None) ) {
      // Backup 2 triggered
      diagnostics->backup[1] = 1;
      // Recompute guesses
      if(params->calc_prim_guess) guess_primitives(eos, metric, prims, cons, prims_guess);
      // Backup routine #2
      check = Hybrid_Select_Method(params, eos, params->backup_routine[1], metric, cons, prims_guess, diagnostics);

      if( (check != 0) && (params->backup_routine[2] != None) ) {
        // Backup 3 triggered
        diagnostics->backup[2] = 1;
        // Recompute guesses
        if(params->calc_prim_guess) guess_primitives(eos, metric, prims, cons, prims_guess);
        // Backup routine #3
        check = Hybrid_Select_Method(params, eos, params->backup_routine[2], metric, cons, prims_guess, diagnostics);
      }
    }
  }
  return check;
}

int Hybrid_Select_Method(const GRHayL_parameters *restrict params,
                         const eos_parameters *restrict eos,
                         const int c2p_key,
                         const metric_quantities *restrict metric,
                         const conservative_quantities *restrict cons,
                         primitive_quantities *restrict prims,
                         con2prim_diagnostics *restrict diagnostics ) {

  switch( c2p_key ) {

    // Noble2D routine (see https://arxiv.org/pdf/astro-ph/0512420.pdf)
    case Noble2D:
      return( Hybrid_Noble2D(params, eos, metric, cons, prims, diagnostics) );
      break;

//    // Noble1D routine (see https://arxiv.org/pdf/astro-ph/0512420.pdf)
//    case Noble1D:
//      return( con2prim_Noble1D(eos,metric,cons,prim,stats) );
//      break;
//
//    // Noble1D entropy routine (see https://arxiv.org/pdf/0808.3140.pdf)
//    case Noble1D_entropy:
//      return( con2prim_Noble1D_entropy(eos,metric,cons,prim,stats) );
//      break;
//
//    // Noble1D entropy 2 routine (see https://arxiv.org/pdf/0808.3140.pdf)
//    case Noble1D_entropy2:
//      return( con2prim_Noble1D_entropy2(eos,metric,cons,prim,stats) );
//      break;
//
//    // Cerda-Duran et al. 2D routine (see https://arxiv.org/pdf/0804.4572.pdf)
//    case CerdaDuran2D:
//      return( con2prim_CerdaDuran2D(eos,metric,cons,prim,stats) );
//
//    // Cerda-Duran et al. 3D routine (see https://arxiv.org/pdf/0804.4572.pdf)
//    case CerdaDuran3D:
//      return( con2prim_CerdaDuran3D(eos,metric,cons,prim,stats) );
//
//    // Palenzuela 1D routine (see https://arxiv.org/pdf/1712.07538.pdf)
//    case Palenzuela1D:
//      return( con2prim_Palenzuela1D(eos,metric,cons,prim,stats) );
//      break;
//
//    // Newman 1D routine (see https://escholarship.org/content/qt0s53f84b/qt0s53f84b.pdf)
//    case Newman1D:
//      return( con2prim_Newman1D(eos,metric,cons,prim,stats) );
//      break;

    default:
      printf("Unknown c2p key in Hybrid_Select_Method (%d). ABORTING!", c2p_key);
      exit(100);
      return 100;
      break;
  }
}
