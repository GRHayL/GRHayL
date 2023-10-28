#include "con2prim.h"

int ghl_hybrid_Font1D_loop(
      const ghl_eos_parameters *restrict eos,
      const int maxits, const double tol,
      const double W_in, const double Sf2_in,
      const double Psim6, const double sdots,
      const double BdotS2, const double B2,
      const ghl_conservative_quantities *restrict cons,
      const double rhob_in, double *restrict rhob_out_ptr) {

  /* Declare basic variables */
  bool Fontcheck=true;
  int itcount = 0, j0, j1;
  double W0, Sf20, rhob0, rhob1, h, P_cold, eps_cold;
  double W   = W_in;
  double Sf2 = Sf2_in;
  double rhob_out=0;

  //////////////////////
  // OUTER LOOP START //
  //////////////////////
  while(Fontcheck && itcount < maxits) {

    /* Set variables to their input values */
    itcount++;
    W0    = W;
    Sf20  = Sf2;
    rhob1 = rhob_in;

    /* Based on rhob_in (i.e. rhob1), determine the
     * polytropic index j1
     */
    j1 = ghl_hybrid_find_polytropic_index(eos,rhob1);

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
      ghl_hybrid_compute_P_cold_and_eps_cold(eos, rhob0, &P_cold, &eps_cold);
      h = 1.0 + eps_cold + P_cold/rhob0;

      /* Update rhob using eq. (A62) in Etienne et al. (2011)
       *          https://arxiv.org/pdf/1112.0568.pdf
       * .---------------------------------------------------------------------------.
       * | rhob = rho_star * Psi^{-6} / sqrt( 1 + S_fluid^{2}/( (rho_star*h)^{2} ) ) |
       * .---------------------------------------------------------------------------.
       */
      rhob1 = cons->rho*Psim6/sqrt(1.0+Sf20/SQR(cons->rho*h));

      /* Update j1 */
      j1 = ghl_hybrid_find_polytropic_index(eos,rhob1);

    }  while(fabs(rhob1-rhob0) > rhob1*tol || j1 != j0);
    //////////////////////
    //  INNER LOOP END  //
    //////////////////////

    /* Output the last value of rhob */
    rhob_out = rhob1;

    /* Perform physical checks on the variables
     * and output the last value of h obtained
     */
    ghl_hybrid_compute_P_cold_and_eps_cold(eos, rhob_out, &P_cold, &eps_cold);
    h = 1.0 + eps_cold + P_cold/rhob_out;

    /* Set W based on eq. (A60) in Etienne et al. (2011)
     *       https://arxiv.org/pdf/1112.0568.pdf
     * .-------------------------------------------------------.
     * | W = psi^{-6} * sqrt( S_fluid^{2} + (rho_star*h)^{2} ) |
     * .-------------------------------------------------------.
     */
    W = sqrt(Sf20 + SQR(cons->rho*h))*Psim6;

    /* Then update S_{fluid}^{2} using eq. (A61) in Etienne et al. (2011)
     *           https://arxiv.org/pdf/1112.0568.pdf
     * .---------------------------------------------------------------------------.
     * | S_fluid^{2} = ( W^{2}*S^{2} + (B.S)^2*(B^{2} + 2W) )/( ( W + B^{2} )^{2} )|
     * .---------------------------------------------------------------------------.
     */
    Sf2 = (SQR(W)*sdots + BdotS2*(B2 + 2.0*W))/SQR(W+B2);

    if (fabs(W-W0) < W*tol && fabs(Sf20-Sf2) < Sf2*tol) Fontcheck=false;

  }
  //////////////////////
  //  OUTER LOOP END  //
  //////////////////////

  /* If the code converged before the max
   * number of iterations were exceeded,
   * return 0, otherwise return 1.
   */
  if(Fontcheck || itcount >= maxits) {
    *rhob_out_ptr = rhob_out;
    return 1;
  }
  else {
    *rhob_out_ptr = rhob_out;
    return 0;
  }
}
