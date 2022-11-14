// Thorn      : IllinoisGRMHD
// File       : con2prim_set_cons_and_prim_from_CONSERVS_and_PRIMS.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This provides functions which 1. convert IllinoisGRMHD's set
//              of conservative variables into the appropriate variables
//              required by the C2P routines and 2. set appropriate primitive
//              guesses.

#include "cctk.h"
#include "con2prim_header.h"
#include "EOS_hybrid_header.h"

void guess_primitives( const eos_parameters *restrict eos,
                       const int c2p_key, const int which_guess,
                       const metric_quantities *restrict metric,
                       const primitive_quantities *restrict prims,
                       const conservative_quantities *restrict cons,
                       primitive_quantities *restrict prims_guess ) {

  *prims_guess = *prims;
  // First the Noble et al. con2prim guesses, which are the more complicated ones
  if( (c2p_key == Noble2D         ) ||
      (c2p_key == Noble1D         ) ||
      (c2p_key == Noble1D_entropy ) ||
      (c2p_key == Noble1D_entropy2) ||
      (c2p_key == CerdaDuran2D    ) ||
      (c2p_key == CerdaDuran3D    ) ) {
    int  polytropic_index = 0;
    double K_ppoly_tab      = 0.0;
    double Gamma_ppoly_tab  = 0.0;
    double rho_b_oldL       = prims->rho;
    double P_oldL           = prims->press;

    if(which_guess==1) {
      //Use a different initial guess:
      rho_b_oldL = cons->rho/metric->psi6;

      /**********************************
       * Piecewise Polytropic EOS Patch *
       *  Finding Gamma_ppoly_tab and K_ppoly_tab *
       **********************************/
      /* Here we use our newly implemented
       * find_polytropic_K_and_Gamma() function
       * to determine the relevant polytropic
       * Gamma and K parameters to be used
       * within this function.
       */
      polytropic_index = find_polytropic_K_and_Gamma_index(eos,rho_b_oldL);
      K_ppoly_tab     = eos->K_ppoly_tab[polytropic_index];
      Gamma_ppoly_tab = eos->Gamma_ppoly_tab[polytropic_index];

      // After that, we compute P_cold
      P_oldL = K_ppoly_tab*pow(rho_b_oldL,Gamma_ppoly_tab);
    }

    if(which_guess==2) {
      //Use atmosphere as initial guess:
      rho_b_oldL = 100.0*eos->rho_atm;

      /**********************************
       * Piecewise Polytropic EOS Patch *
       *  Finding Gamma_ppoly_tab and K_ppoly_tab *
       **********************************/
      /* Here we use our newly implemented
       * find_polytropic_K_and_Gamma() function
       * to determine the relevant polytropic
       * Gamma and K parameters to be used
       * within this function.
       */
      polytropic_index = find_polytropic_K_and_Gamma_index(eos,rho_b_oldL);
      K_ppoly_tab     = eos->K_ppoly_tab[polytropic_index];
      Gamma_ppoly_tab = eos->Gamma_ppoly_tab[polytropic_index];

      // After that, we compute P_cold
      P_oldL = K_ppoly_tab*pow(rho_b_oldL,Gamma_ppoly_tab);
    }

    prims_guess->rho = rho_b_oldL;
    prims_guess->press = P_oldL;
  }

  if( eos->eos_type == 1 ) {
    // This one is very simple! The only guess required is the temperature
    if( which_guess == 1 ) {
      prims_guess->temp = eos->temp_atm;
    }
    else {
      prims_guess->temp = eos->temp_max;
    }
    prims_guess->Y_e = cons->Y_e/cons->rho;
  }

}
