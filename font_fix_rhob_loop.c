#include "con2prim_header.h"
#include "EOS_hybrid_header.h"

/* Function    : font_fix_rhob_loop()
 * Authors     : Leo Werneck
 * Description : Determines rhob using the font fix prescription
 * Dependencies: find_polytropic_K_and_Gamma_index()
 *             : compute_P_cold__eps_cold()
 * Reference   : Etienne et al. (2011) [https://arxiv.org/pdf/1112.0568.pdf]
 *
 * Inputs      : maxits          - maximum number of iterations allowed
 *             : tol             - font fix tolerance
 *             : W               - See eq. (A26)
 *             : Sf2             - S_{fluid}^{2}, see eq. (A24)
 *             : Psim6           - This is equal to sqrt(\gamma)
 *             : sdots           - \tilde{S}_{\mu}\tilde{S}^{\mu}
 *             : BbardotS2       - (\bar{B}^{\mu}S_{\mu})^{2},
 *             : B2bar           - \bar{B}^{2}, see eq. (A28)
 *             : cons            - Array of conservative variables
 *             : eos             - Struct of EOS parameters
 *             : rhob_in         - Initial value of rhob
 *             : rhob_out        - Output variable
 *
 * Outputs     : rhob_out        - Updated value of rhob
 *             : return value: 0 - Font fix worked
 *             : return value: 1 - Font fix failed
 */
int font_fix_rhob_loop( const eos_parameters *restrict eos,
                        const int maxits, const double tol,
                        const double W_in, const double Sf2_in,
                        const double Psim6, const double sdots,
                        const double BbardotS2, const double B2bar,
                        const conservative_quantities *restrict cons,
                        const double rhob_in, double *restrict rhob_out_ptr ) {

  double rhob_out = *rhob_out_ptr;

  /* Declare basic variables */
  bool fontcheck=true;
  int itcount = 0, j0, j1;
  double W0, Sf20, rhob0, rhob1, h, P_cold, eps_cold;
  double W   = W_in;
  double Sf2 = Sf2_in;

  //////////////////////
  // OUTER LOOP START //
  //////////////////////
  while(fontcheck && itcount < maxits) {

    /* Set variables to their input values */
    itcount++;
    W0    = W;
    Sf20  = Sf2;
    rhob1 = rhob_in;

    /* Based on rhob_in (i.e. rhob1), determine the
     * polytropic index j1
     */
    j1 = find_polytropic_K_and_Gamma_index(eos,rhob1);

    //////////////////////
    // INNER LOOP START //
    //////////////////////
    do {

      /* Set rhob0/j0 to be equal to the rhob/j used
       * in the previous iteration, i.e. rhob1/j1.
       */
      rhob0 = rhob1;
      j0    = j1;

      /* Compute h using h_cold and our polytropic EOS
       * .------------------------------------------.
       * | h = h_cold = 1 + eps_cold + P_cold/rhob. |
       * .------------------------------------------.
       */
      compute_P_cold__eps_cold(eos,rhob0, &P_cold, &eps_cold);
      h = 1.0 + eps_cold + P_cold/rhob0;

      /* Update rhob using eq. (A62) in Etienne et al. (2011)
       *          https://arxiv.org/pdf/1112.0568.pdf
       * .---------------------------------------------------------------------------.
       * | rhob = rho_star * Psi^{-6} / sqrt( 1 + S_fluid^{2}/( (rho_star*h)^{2} ) ) |
       * .---------------------------------------------------------------------------.
       */
      rhob1 = cons->rho*Psim6/sqrt(1.0+Sf20/SQR(cons->rho*h));

      /* Update j1 */
      j1 = find_polytropic_K_and_Gamma_index(eos,rhob1);

    }  while( fabs(rhob1-rhob0) > rhob1*tol || j1 != j0);
    //////////////////////
    //  INNER LOOP END  //
    //////////////////////

    /* Output the last value of rhob */
    rhob_out = rhob1;

    /* Perform physical checks on the variables
     * and output the last value of h obtained
     */
    compute_P_cold__eps_cold(eos,rhob_out, &P_cold, &eps_cold);
    h = 1.0 + eps_cold + P_cold/rhob_out;

    /* Set W based on eq. (A60) in Etienne et al. (2011)
     *       https://arxiv.org/pdf/1112.0568.pdf
     * .-------------------------------------------------------.
     * | W = psi^{-6} * sqrt( S_fluid^{2} + (rho_star*h)^{2} ) |
     * .-------------------------------------------------------.
     */
    W = sqrt( Sf20 + SQR(cons->rho*h))*Psim6;

    /* Then update S_{fluid}^{2} using eq. (A61) in Etienne et al. (2011)
     *           https://arxiv.org/pdf/1112.0568.pdf
     * .---------------------------------------------------------------------------.
     * | S_fluid^{2} = ( W^{2}*S^{2} + (B.S)^2*(B^{2} + 2W) )/( ( W + B^{2} )^{2} )|
     * .---------------------------------------------------------------------------.
     */
    Sf2 = (SQR(W)*sdots + BbardotS2*(B2bar + 2.0*W))/SQR(W+B2bar);

    if ( fabs(W-W0) < W*tol && fabs(Sf20-Sf2) < Sf2*tol) fontcheck=false;

  }
  //////////////////////
  //  OUTER LOOP END  //
  //////////////////////

  /* If the code converged before the max
   * number of iterations were exceeded,
   * return 0, otherwise return 1.
   */
  if(fontcheck || itcount >= maxits) {
    *rhob_out_ptr = rhob_out;
    return 1;
  }
  else {
    *rhob_out_ptr = rhob_out;
    return 0;
  }
}
