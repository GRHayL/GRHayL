#include "con2prim_header.h"

int con2prim( const GRMHD_parameters *restrict params, const eos_parameters *restrict eos,
              const metric_quantities *restrict metric, conservative_quantities *restrict cons,
              primitive_quantities *restrict prims, con2prim_diagnostics *restrict diagnostics ) {

  // declare some variables for HARM.
  conservative_quantities cons_undens;
  primitive_quantities prims_guess;

  /*
    -- Driver for new prim. var. solver.  The driver just translates
    between the two sets of definitions for U and P.  The user may
    wish to alter the translation as they see fit.


    //        /     rho u^t     \                           //
    //    U = | T^t_t + rho u^t | * sqrt(-det(g_{\mu\nu}))  //
    //        |     T^t_i       |                           //
    //        \      B^i        /                           //
    //                                                      //
    //        /     rho     \                               //
    //    P = |     uu      |                               //
    //        | \tilde{u}^i |                               //
    //        \     B^i     /                               //

    (above equations have been fixed by Yuk Tung & Zach)
  */

  // U[NPR]    = conserved variables (current values on input/output);
  // g4dn[NDIM][NDIM] = covariant form of the 4-metric ;
  // g4up[NDIM][NDIM] = contravariant form of the 4-metric ;
  // gdet             = sqrt( - determinant of the 4-metric) ;
  // prim[NPR] = primitive variables (guess on input, calculated values on
  //                     output if there are no problems);

  // U[1]   =
  // U[2-4] =  stildei + rhostar

  // Other ideas for setting the gamma speed limit
  //double GAMMA_SPEED_LIMIT = 100.0;
  //if(METRIC_LAP_PSI4[PSI6]>Psi6threshold) GAMMA_SPEED_LIMIT=500.0;
  //if(METRIC_LAP_PSI4[PSI6]>Psi6threshold) GAMMA_SPEED_LIMIT=100.0;

  //FIXME: Only works if poisoning is turned on. Otherwise will access unknown memory. This trick alone speeds up the whole code (Cowling) by 2%.
  //int startguess=0;
  //if(robust_isnan(PRIMS[VX])) startguess=1;
  int startguess=1;

  double u0L=1.0;

  for(int which_guess=startguess;which_guess<3;which_guess++) {

    // Set the conserved variables required by the con2prim routine
    undensitize( eos, params->main_routine, metric, prims, cons, &cons_undens );

    // Set primitive guesses
    guess_primitives( eos, params->main_routine, which_guess, metric, prims, cons, &prims_guess );

    /************* Conservative-to-primitive recovery ************/

    int check = con2prim_select(eos, params->main_routine, metric, &cons_undens, &prims_guess, diagnostics);

    if( (check != 0) && (params->backup_routine[0] != None) ) {
      // Backup 1 triggered
      diagnostics->backup[0] = 1;
      // Recompute guesses
      undensitize( eos,params->backup_routine[0], metric, prims, cons, &cons_undens );
      guess_primitives( eos,params->backup_routine[0], which_guess, metric, prims, cons, &prims_guess );
      // Backup routine #1
      check = con2prim_select(eos, params->backup_routine[0], metric, &cons_undens, &prims_guess, diagnostics);

      if( (check != 0) && (params->backup_routine[1] != None) ) {
        // Backup 2 triggered
        diagnostics->backup[1] = 1;
        // Recompute guesses
        undensitize( eos,params->backup_routine[1], metric, prims, cons, &cons_undens );
        guess_primitives( eos,params->backup_routine[1], which_guess, metric, prims, cons, &prims_guess );
        // Backup routine #2
        check = con2prim_select(eos, params->backup_routine[1], metric, &cons_undens, &prims_guess, diagnostics);

        if( (check != 0) && (params->backup_routine[2] != None) ) {
          // Backup 3 triggered
          diagnostics->backup[2] = 1;
          // Recompute guesses
          undensitize( eos,params->backup_routine[2], metric, prims, cons, &cons_undens );
          guess_primitives( eos,params->backup_routine[2], which_guess, metric, prims, cons, &prims_guess );
          // Backup routine #3
          check = con2prim_select(eos, params->backup_routine[2], metric, &cons_undens, &prims_guess, diagnostics);
        }
      }
    }
    /*************************************************************/

    if(check!=0) {
      check = font_fix(eos, metric, &cons_undens, prims, &prims_guess, &diagnostics, &u0L);
      diagnostics->font_fixes+=1;
    }

    if(check==0) {
//TODO: set up errors
      // Check for NAN!
//      if( robust_isnan(prim[RHO]*prim[TEMP]*prim[YE]*prim[PRESS]*prim[EPS]*prim[ENT]*utx_new*uty_new*utz_new*u0L) ) {
//        CCTK_VINFO("***********************************************************");
//        CCTK_VINFO("NAN found in function %s (file: %s)",__func__,__FILE__);
//        CCTK_VINFO("Input IllinoisGRMHD conserved variables:");
//        CCTK_VINFO("rho_*, Ye_*, ~tau, ~S_{i}: %e %e %e %e %e %e",CONSERVS[RHOSTAR],CONSERVS[YESTAR],CONSERVS[TAUENERGY],CONSERVS[STILDEX],CONSERVS[STILDEY],CONSERVS[STILDEZ]);
//        CCTK_VINFO("Input con2prim conserved variables:");
//        CCTK_VINFO("D, DYe, tau, S_{i}: %e %e %e %e %e %e",cons[RHO],cons[YE],cons[TAU],cons[S1_cov],cons[S2_cov],cons[S3_cov]);
//        CCTK_VINFO("Output primitive variables:");
//        CCTK_VINFO("rho, T, Ye: %e %e %e",prim[RHO],prim[TEMP],prim[YE]);
//        CCTK_VINFO("P, eps, S : %e %e %e",prim[PRESS],prim[EPS],prim[ENT]);
//        CCTK_VINFO("u^{mu}    : %e %e %e %e",u0L,utx_new,uty_new,utz_new);
//        CCTK_VINFO("***********************************************************");
//      }

        *prims = prims_guess;

	return 0;
//      }
    } else {
      //If we didn't find a root, then try again with a different guess.
    }
  }
  return 1;
}
