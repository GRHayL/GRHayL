#include "../../con2prim_gem.h"

int C2P_Select_Hybrid_Method( const GRHayL_parameters *restrict params,
                              const eos_parameters *restrict eos, const int c2p_key,
                              const metric_quantities *restrict metric,
                              const conservative_quantities *restrict cons,
                              primitive_quantities *restrict prims,
                              con2prim_diagnostics *restrict diagnostics ) {

  switch( c2p_key ) {

    // Noble2D routine (see https://arxiv.org/pdf/astro-ph/0512420.pdf)
    case Noble2D:
      return( C2P_Hybrid_Noble2D(params, eos, metric, cons, prims, diagnostics) );
      break;

    case OldNoble2D:
      return( C2P_Hybrid_OldNoble2D(eos,metric,cons,prims,diagnostics) );
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
//TODO: errors      CCTK_VERROR("Unknown c2p key in con2prim_select (%d). ABORTING!",c2p_key);
      return 100;
      break;

  }
}
