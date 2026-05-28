#include "ghl_con2prim.h"

/**
 * @ingroup c2p_internal
 * @brief Computes the density using the Font1D method
 *
 * @details
 * This routine is the core algorithm used by @ref ghl_hybrid_Font1D . It
 * solves for the primitive density by assuming that there are now thermal
 * contributions to the solution.
 *
 * @param[in] eos:          pointer to ghl_eos_parameters struct
 *
 * @param[in] maxits:       maximum number of iterations
 *
 * @param[in] tol:          maximum allowed error for convergence
 *
 * @param[in] W_in:         initial guess for quantity \f$ W \f$
 *
 * @param[in] Sf2_in:       initial guess for quantityf \f$ S_\mathrm{fluid}^2 \f$
 *
 * @param[in] Psim6:        inverse of \f$ \psi^6 = \sqrt{|\gamma|} \f$
 *
 * @param[in] sdots:        conservative quantity \f$ S^2 \f$
 *
 * @param[in] BdotS2:       conservative quantity \f$ B \cdot S\f$
 *
 * @param[in] B2:           magnetic quantity \f$ B^2 \f$
 *
 * @param[in] cons_undens:  pointer to ghl_conservative_quantities struct with
 *                          **undensitized** conservative variables
 *
 * @param[in] rhob_in:      initial guess for the density \f$ \rho \f$
 *
 * @param[out] rhob_out_ptr: returned density value
 *
 * @returns error code for any Con2Prim failures
 */
ghl_error_codes_t ghl_hybrid_Font1D_loop(
      const ghl_eos_parameters *restrict eos,
      const int maxits,
      const double tol,
      const double W_in,
      const double Sf2_in,
      const double Psim6,
      const double sdots,
      const double BdotS2,
      const double B2,
      const ghl_conservative_quantities *restrict cons_undens,
      const double rhob_in,
      double *restrict rhob_out_ptr) {

  bool Font_success=false;
  int itcount = 0, j0, j1;
  double W0, Sf20, rhob0, rhob1, h, P_cold, eps_cold;
  double W   = W_in;
  double Sf2 = Sf2_in;
  double rhob_out=0;

  /**
   * # Function Step-by-Step
   *
   * This function iterates for maxits iterations or converging to within
   * tolerance of a solution.
   */
  while(!Font_success && itcount < maxits) {

    itcount++;
    W0    = W;
    Sf20  = Sf2;
    rhob1 = rhob_in;

    /**
     * For each iteration, we find the polytropic index \f$ j \f$ using
     * @ref ghl_hybrid_find_polytropic_index . Then, we enter a subloop to find
     * the next guess for the density using eq. (A62) of \cite Etienne_2012
     *
     * \f[
     * \rho_1 = \frac{D}{\psi^6 \sqrt{1 + \left(\frac{S_\mathrm{fluid}^2}{D h}\right)^2}}
     * \f]
     *
     * where
     *
     * \f[
     * h = 1 + \epsilon_\mathrm{cold} + \frac{P_\mathrm{cold}}{\rho}
     * \f]
     *
     * The subloop continues until the density is in the same polytropic piece
     * and
     *
     * \f[
     * \left| \rho_1 - \rho_0 \right| < t \rho_1
     * \f]
     *
     * where \f$ \rho_0 \f$ is the value from the previous iteration and \f$ t \f$
     * is the tolerance.
     */
    j1 = ghl_hybrid_find_polytropic_index(eos,rhob1);

    do {
      rhob0 = rhob1;
      j0    = j1;

      ghl_hybrid_compute_P_cold_and_eps_cold(eos, rhob0, &P_cold, &eps_cold);
      h = 1.0 + eps_cold + P_cold/rhob0;

      rhob1 = cons_undens->rho*Psim6/sqrt(1.0+Sf20/SQR(cons_undens->rho*h));

      j1 = ghl_hybrid_find_polytropic_index(eos,rhob1);

    }  while(fabs(rhob1-rhob0) > rhob1*tol || j1 != j0);

    rhob_out = rhob1;

    /**
     * Once we find a sufficiently converged density, we compute the new
     * \f$ W \f$ and \f$ S_\mathrm{fluid}^2 \f$ using eqs. A60 and A61 from
     * \cite Etienne_2012, respectively:
     *
     * \f[
     * \begin{aligned}
     * W &= \frac{\sqrt{S_\mathrm{fluid}^4 + \left( D h \right)^2}}{\psi^6} \\
     * S_\mathrm{fluid}^2 &= \frac{W^2 S^2
     *                           + \left( B \cdot S \right)^2 \left( B^2 + 2W \right)}
     *                           {\left( W + B^2 \right)^2}
     * \end{aligned}
     * \f]
     *
     * If these are within the tolerance bounds
     *
     * \f[
     * \begin{aligned}
     * \left| W_1 - W_0 \right| &< t W_1 \\
     * \left| \left(S^2_\mathrm{fluid}\right)_1 - \left(S^2_\mathrm{fluid}\right)_0 \right|
     *     &< t \left(S^2_\mathrm{fluid}\right)_1
     * \end{aligned}
     * \f]
     *
     * then the loop ends and the function returns the density value. As before,
     * \f$ t \f$ is the input tolerance and the (0,1) subscripts represent the
     * previous and current iterations, respectively.
     */
    ghl_hybrid_compute_P_cold_and_eps_cold(eos, rhob_out, &P_cold, &eps_cold);
    h = 1.0 + eps_cold + P_cold/rhob_out;

    W = sqrt(Sf20 + SQR(cons_undens->rho*h))*Psim6;
    Sf2 = (SQR(W)*sdots + BdotS2*(B2 + 2.0*W))/SQR(W+B2);

    if (fabs(W-W0) < W*tol && fabs(Sf20-Sf2) < Sf2*tol) Font_success=true;
  }

  if(!Font_success || itcount >= maxits) {
    *rhob_out_ptr = rhob_out;
    return ghl_error_c2p_max_iter;
  } else {
    *rhob_out_ptr = rhob_out;
    return ghl_success;
  }
}
