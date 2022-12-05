#include "EOS_hybrid_header.h"

/* Function    : setup_K_ppoly_and_eps_integ_consts()
 * Authors     : Leo Werneck
 * Description : For a given set of EOS inputs, determine
 *               values of K_ppoly that will result in a
 *               everywhere continuous P_cold function.
 *
 * Inputs      : eos                  - an initialized eos_parameters struct
 *                                      with data for the EOS of the simulation
 *
 * Outputs     : eos->K_ppoly     - fully populated array of K_ppoly
 *                                      to be used in each polytropic EOS.
 *             : eos->eps_integ_const - fully populated array of C_{j}'s,
 *                                      used to compute eps_cold for
 *                                      a piecewise polytropic EOS.
 */
void setup_K_ppoly_and_eps_integ_consts(eos_parameters *restrict eos) {

  /* When neos = 1, we will only need the value K_ppoly[0] and eps_integ_const[0].
   * Since our only polytropic EOS is given by
   *  -----------------------------------
   * | P_{0} = K_{0} * rho ^ (Gamma_{0}) | ,
   *  -----------------------------------
   * and, therefore,
   *  ---------------------------------------------------------------
   * | eps_{0} = K_{0} * rho ^ (Gamma_{0}-1) / (Gamma_{0}-1) + C_{0} | ,
   *  ---------------------------------------------------------------
   * we only need to set up K_{0} := K_ppoly[0] and C_{0} := eps_integ_const[0].
   * K_{0} is a user input, so we need to do nothing. C_{0}, on the other hand,
   * is fixed by demanding that eps(rho) -> 0 as rho -> 0. Thus, C_{0} = 0.
   */
  eos->eps_integ_const[0] = 0.0;
  if(eos->neos==1) return;

  /********************
   * Setting up K_{j} *
   ********************/
  /* When neos > 1, we have the following structure
   *
   *           /      K_{0} * rho^(Gamma_{0})      ,                 rho <  rho_{0}
   *           |      K_{1} * rho^(Gamma_{1})      ,      rho_{0} <= rho <  rho_{1}
   *           |              ...
   * P(rho) = <       K_{j} * rho^(Gamma_{j})      ,    rho_{j-1} <= rho <  rho_{j}
   *           |              ...
   *           | K_{neos-2} * rho^(Gamma_{neos-2}) , rho_{neos-3} <= rho <  rho_{neos-2}
   *           \ K_{neos-1} * rho^(Gamma_{neos-1}) ,                 rho >= rho_{neos-2}
   *
   * Imposing that P(rho) be everywhere continuous, we have
   *  -------------------------------------------------
   * | K_{j} = K_{j-1} * rho^(Gamma_{j-1} - Gamma_{j}) |
   *  -------------------------------------------------
   */
  for(int j=1; j<eos->neos; j++) {
    // Set a useful auxiliary variable to keep things more compact:
    // First, (Gamma_{j-1} - Gamma_{j}):
    double Gamma_diff = eos->Gamma_ppoly[j-1] - eos->Gamma_ppoly[j];

    // Implement the boxed equation above, using our auxiliary variable:
    eos->K_ppoly[j] = eos->K_ppoly[j-1] * pow(eos->rho_ppoly[j-1],Gamma_diff);
  }

  /********************
   * Setting up C_{j} *
   ********************/
  /* When neos > 1, we have the following structure (let neos->N):
   *
   *             /      K_{0}*rho^(Gamma_{0}-1)/(Gamma_{0}-1)  + C_{0},                rho <  rho_{0}
   *             |      K_{1}*rho^(Gamma_{1}-1)/(Gamma_{1}-1)  + C_{1},     rho_{0} <= rho <  rho_{1}
   *             |                     ...
   * eps(rho) = <       K_{j}*rho^(Gamma_{j}-1)/(Gamma_{j}-1)  + C_{j},   rho_{j-1} <= rho <  rho_{j}
   *             |                     ...
   *             | K_{N-2}*rho^(Gamma_{N-2}-1)/(Gamma_{N-2}-1) + C_{N-2}, rho_{N-3} <= rho <  rho_{N-2}
   *             \ K_{N-1}*rho^(Gamma_{N-1}-1)/(Gamma_{N-1}-1) + C_{N-1},              rho >= rho_{N-2}
   *
   * Imposing that eps_{cold}(rho) be everywhere continuous, we have
   *  ---------------------------------------------------------------
   * | C_{j} = C_{j-1}                                               |
   * |       + ( K_{j-1}*rho_{j-1}^(Gamma_{j-1}-1) )/(Gamma_{j-1}-1) |
   * |       - ( K_{j+0}*rho_{j-1}^(Gamma_{j+0}-1) )/(Gamma_{j+0}-1) |
   *  ---------------------------------------------------------------
   */
  for(int j=1; j<eos->neos; j++) {
    // Set a few useful auxiliary variables to keep things more compact:
    // First, (Gamma_{j-1}-1):
    double Gammajm1m1 = eos->Gamma_ppoly[j-1] - 1.0;

    // Then, (Gamma_{j+0}-1):
    double Gammajp0m1 = eos->Gamma_ppoly[j+0] - 1.0;

    // Next, ( K_{j-1}*rho_{j-1}^(Gamma_{j-1}-1) )/(Gamma_{j-1}-1):
    double aux_epsm1  = eos->K_ppoly[j-1]*pow(eos->rho_ppoly[j-1],Gammajm1m1)/Gammajm1m1;

    // Finally, ( K_{j+0}*rho_{j+0}^(Gamma_{j+0}-1) )/(Gamma_{j+0}-1):
    double aux_epsp0  = eos->K_ppoly[j+0]*pow(eos->rho_ppoly[j-1],Gammajp0m1)/Gammajp0m1;

    // Implement the boxed equation above, using our auxiliary variables:
    eos->eps_integ_const[j] = eos->eps_integ_const[j-1] + aux_epsm1 - aux_epsp0;
  }
}

