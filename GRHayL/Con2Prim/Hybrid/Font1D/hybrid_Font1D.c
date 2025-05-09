#include "ghl_con2prim.h"

ghl_error_codes_t ghl_hybrid_Font1D_loop(
      const ghl_eos_parameters *restrict eos,
      const int maxits, const double tol, const double W_in,
      const double Sf2_in, const double Psim6, const double sdots,
      const double BdotS2, const double B2,
      const ghl_conservative_quantities *restrict cons_undens,
      const double rhob_in, double *restrict rhob_out_ptr);

/**
 * @ingroup hyb_c2p
 * @brief Computes primitive variables using the Font1D method
 *
 * @details
 * This function uses the method described in Appendix A.4 of \cite Etienne_IGM
 * to recover the primitive variables from the conservative variables. This
 * method comes from \cite Font_2000, which outlines the idea of the method in
 * the text preceeding Eq. 84. This method serves as a backup for very low
 * densities where other hybrid Con2Prim routines might recover a negative
 * pressure. Font1D assumes that \f$ P=P_\mathrm{cold} \f$ and solves for
 * \f$ \rho \f$ using an iterative method. This method is guaranteed to
 * succees so long as the conservatives are physically reasonable: 
 * \f$ \rho_* > 0 \f$ and \f$ \tilde{S}_i \in (-\infty, \infty) \f$.
 * The return value gives information on the success or failure of the
 * recovery attempt.
 *
 * @param[in] params: pointer to ghl_parameters struct
 *
 * @param[in] eos:    pointer to ghl_eos_parameters struct
 *
 * @param[in] ADM_metric:   pointer to ghl_metric_quantities struct with ADM metric
 *
 * @param[in] metric_aux:   pointer to ghl_ADM_aux_quantities struct
 *
 * @param[in] cons_undens:  pointer to ghl_conservative_quantities struct with
 *                          **undensitized** conservative variables
 *
 * @param[in,out] prims:    pointer to ghl_primitive_quantities struct;
 *                          input is initial guess for iterative solver;
 *                          output is the primitives consistent with the
 *                          input conservatives
 *
 * @param[out] diagnostics: pointer to ghl_con2prim_diagnostics struct; returns
 *                          with several Con2Prim solver diagnostics
 *
 * @returns error code for any Con2Prim failures
 */
ghl_error_codes_t ghl_hybrid_Font1D(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  double utU[3];

  /**
   * # Function Step-by-Step
   *
   * We start by computing several intermediate quantities dependent on
   * \f$ B^i \f$ and \f$ S^i \f$ using @ref ghl_compute_vec2_from_vec3D().
   * Additionally, if \f$ S^2 \approx 0 \f$ we can immediately assume the
   * velocity is zero and skip the iterative solve.
   */
  const double sdots = ghl_compute_vec2_from_vec3D(ADM_metric->gammaUU, cons_undens->SD);

  const double B2 = ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, prims->BU);

  double BdotS, BdotS2, hatBdotS;
  if(B2 < 1e-150) {
    BdotS = BdotS2 = hatBdotS = 0.0;
  } else {
    const double B_mag = sqrt(B2);
    BdotS = prims->BU[0]*cons_undens->SD[0] + prims->BU[1]*cons_undens->SD[1] + prims->BU[2]*cons_undens->SD[2];
    BdotS2 = BdotS*BdotS;
    hatBdotS = BdotS/B_mag;
  }

  double Psim6 = 1.0/ADM_metric->sqrt_detgamma;

  double rhob;
  if (sdots<1.e-300) {
    utU[0] = 0.0;
    utU[1] = 0.0;
    utU[2] = 0.0;
    prims->u0    = ADM_metric->lapseinv;
    prims->vU[0] = -ADM_metric->betaU[0];
    prims->vU[1] = -ADM_metric->betaU[1];
    prims->vU[2] = -ADM_metric->betaU[2];
  } else {
    /**
     * Now we apply core @ref ghl_hybrid_Font1D_loop method with several
     * attempts. As this method is intended as a final resort backup, it tries
     * 5 times, with the max iterations scaling from 300 to 500 and the
     * tolerance cutoff from \f$ 10^{-15} \f$ to \f$ 2.56 \times 10^{-13} \f$.
     * Each attempt also inherits the converged result from the previous attempt
     * to improve the likelihood of success.
     *
     * We also compute several initial guess quantities outside the loop so
     * we can reuse the values.
     */
    double W0    = sqrt( SQR(hatBdotS) + SQR(cons_undens->rho) ) * Psim6;
    double Sf2  = (SQR(W0)*sdots + BdotS2*(B2 + 2.0*W0))/SQR(W0+B2);
    double rhob0 = cons_undens->rho*Psim6/sqrt(1.0+Sf2/SQR(cons_undens->rho));

    const int maxits = 300;
    double tol = 1.e-15;
    int error;
    const int Font1D_attempts = 5;
    for(int n=0; n<Font1D_attempts; n++) {
      const int loop_maxits = maxits + n*50; // From 300 to 500 for 5 iterations
      const double loop_tol = tol*pow(4,n); // tolerance multipliers are {0,4,16,64,256}
      error = ghl_hybrid_Font1D_loop(
            eos, loop_maxits, loop_tol, W0, Sf2, Psim6,
            sdots, BdotS2, B2, cons_undens, rhob0, &rhob);
      rhob0 = rhob;
      if(!error) break;
    }

    if(error)
      return error;

    /**
     * If @ref ghl_hybrid_Font1D_loop succeeds, then we compute the other
     * primitive variables. First we compute \f$ P_\mathrm{cold} \f$,
     * \f$ \epsilon_\mathrm{cold} \f$, and
     *
     * \f[
     * h_\mathrm{cold} = 1 + \epsilon_\mathrm{cold} + \frac{P_\mathrm{cold}}{\rho}
     * \f]
     */
    double P_cold, eps_cold;
    ghl_hybrid_compute_P_cold_and_eps_cold(eos, rhob, &P_cold, &eps_cold);
    double h = 1.0 + eps_cold + P_cold/rhob;

    /**
     * We then compute \f$ \gamma_v \f$ using equation (A19) in \cite Etienne_2012 :
     *
     * \f[
     * \gamma_v = \frac{D}{\psi^6 \rho}
     * \f]
     */
    double gammav = cons_undens->rho*Psim6/rhob;

    /**
     * Finally, we compute and speed-limit the velocity
     * \f[
     * \tilde{u}^i = f_2 \left( S^i + f_1 B^i \right)
     * \f]
     *
     * where
     *
     * \f[
     * f_1 = \frac{\sqrt{|\gamma|} B \cdot S}{\gamma_v D h}
     * \f]
     *
     * and
     *
     * \f[
     * f_2 = \left[ D h +  \frac{\sqrt{|\gamma|} B^2}{\gamma_v} \right]^{-1}
     * \f]
     */
    double rhosh = cons_undens->rho*h;
    double fac1 = ADM_metric->sqrt_detgamma*BdotS/(gammav*rhosh);
    double fac2 = 1.0/(rhosh + ADM_metric->sqrt_detgamma*B2/gammav);

    double SU[3];
    ghl_raise_lower_vector_3D(ADM_metric->gammaUU, cons_undens->SD, SU);

    utU[0] = fac2*(SU[0] + fac1*prims->BU[0]);
    utU[1] = fac2*(SU[1] + fac1*prims->BU[1]);
    utU[2] = fac2*(SU[2] + fac1*prims->BU[2]);
    diagnostics->speed_limited = ghl_limit_utilde_and_compute_v(params, ADM_metric, utU, prims);
  }

  /**
   * The Font fix only sets the velocities. We set the remaining primitives using
   * @ref ghl_hybrid_compute_P_cold, @ref ghl_hybrid_compute_entropy_function,
   * and
   *
   * \f[
   * \begin{aligned}
   * \rho &= \frac{D}{\alpha u^0 \sqrt{|\gamma|}} \\
   * \epsilon &= \frac{P}{\rho (\Gamma - 1)}
   * \end{aligned}
   * \f]
   */
  prims->rho = cons_undens->rho/(ADM_metric->lapse*prims->u0*ADM_metric->sqrt_detgamma);

  double K_ppoly, Gamma_ppoly;
  ghl_hybrid_get_K_and_Gamma(eos, prims->rho, &K_ppoly, &Gamma_ppoly);

  ghl_hybrid_compute_P_cold(eos, prims->rho, &prims->press);

  prims->eps = prims->press/(prims->rho*(Gamma_ppoly-1.0));
  if(params->evolve_entropy)
    prims->entropy = ghl_hybrid_compute_entropy_function(eos, prims->rho, prims->press);

  diagnostics->which_routine = Font1D;
  return ghl_success;
}
