#include "con2prim_header.h"
#include "EOS/Hybrid/EOS_hybrid_header.h"
#include <stdio.h>

/* Function    : guess_primitives()
 * Authors     : Leo Werneck and Samuel Cupp
 * Description : This function computes initial guesses for the primitives
                 for use in the Con2Prim solver.
 * Dependencies: 
 *
 * Inputs      : eos            - an initialized eos_parameters struct
 *                                with data for the EOS of the simulation
 *             : which_guess    - an integer which selects the initial guess
 *                                to be used for the primitives
 *             : metric         - an initialized metric_quantities struct
 *                                with data for the gridpoint of interest
 *             : cons           - an initialized conservative_quantities
 *                                struct with data for the gridpoint of
 *                                interest
 *
 * Outputs     : prims_guess    - primitive_quantities struct filled with
 *                                an initial primitives guess for the
 *                                Con2Prim solver
 */

//TODO: consider passing in cons_undens instead. I don't think we need densitized.
void guess_primitives( const eos_parameters *restrict eos,
                       const int which_guess,
                       const metric_quantities *restrict metric,
                       const primitive_quantities *restrict prims,
                       const conservative_quantities *restrict cons,
                       primitive_quantities *restrict prims_guess ) {

  *prims_guess = *prims;

  if(which_guess==1) {
    //Use a different initial guess:
    prims_guess->rho = cons->rho/metric->psi6;
    prims_guess->temp = eos->temp_atm;
  } else if(which_guess==2) {
    //Use atmosphere as initial guess:
    prims_guess->rho = 100.0*eos->rho_atm;
    prims_guess->temp = eos->temp_max;
  } else {
    printf("WARNING: guess_primitives was passed an unknown guess type of %d.", which_guess);
  }

  // TODO: Hybrid only?
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
  int polytropic_index = find_polytropic_K_and_Gamma_index(eos, prims_guess->rho);
  double K_ppoly_tab      = eos->K_ppoly_tab[polytropic_index];
  double Gamma_ppoly_tab  = eos->Gamma_ppoly_tab[polytropic_index];

  // After that, we compute P_cold
  prims_guess->press = K_ppoly_tab*pow(prims_guess->rho, Gamma_ppoly_tab);
  prims_guess->Y_e = cons->Y_e/cons->rho;
}
