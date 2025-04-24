#include "ghl_con2prim.h"

/**
 * @ingroup Con2Prim
 * @brief Enforces constraining inequalities on \f$ \tilde{\tau} \f$ and \f$ \tilde{S} \f$.
 *
 * @details
 * This function applies limits to \f$ \tilde{\tau} \f$ and \f$ \tilde{S_i} \f$
 * based on the work by \cite Faber_2007 and \cite Etienne_2012.
 * ghl_con2prim_diagnostics::tau_fix and ghl_con2prim_diagnostics::Stilde_fix
 * record whether any limits were applied.
 *
 * The only primitive that needs to be provided is ghl_primitives_quantities::BU.
 * This is used before the Con2Prim solver, so \f$ B^i \f$ is the only available
 * primitive.
 *
 * ## Function Step-by-Step
 *
 * ## Computing Intermediate Quantities
 *
 * We first compute \f$ \tilde{S}^2 \f$ and \f$B^2\f$ using
 * @ref ghl_compute_vec2_from_vec3D(). Then, we compute \f$ B\cdot S \f$,
 * \f$ \hat{B}\cdot S \f$ and some other intermediate quantities:
 *
 * \f[
 * \Omega \equiv \mathrm{\texttt{BdotS}} = B^i \tilde{S}_i
 * \f]
 *
 * \f[
 * \hat{\Omega} \equiv \mathrm{\texttt{hatBdotS}} = \frac{\Omega}{|B|}
 * \f]
 *
 * \f[
 * W_m \equiv \mathrm{\texttt{Wm}} = \frac{\sqrt{\hat{\Omega}^2 + \rho^2}}{\sqrt{\gamma}}
 * \f]
 *
 * \f[
 * S_m^2 \equiv \mathrm{\texttt{Sm2}}
 *   = \frac{W_m^2 \tilde{S}^2 + \Omega^2 (B^2 + 2 W_m)}{W_m + B^2}
 * \f]
 *
 * \f[
 * W_\mathrm{min} \equiv \mathrm{\texttt{Wmin}}
 *   = \frac{\sqrt{S_m^2 + \rho^2}}{\sqrt{\gamma}}
 * \f]
 *
 * \f[
 * \mathrm{\texttt{half_psi6_B2}} = \frac{\sqrt{\gamma}B^2}{2}
 * \f]
 *
 * \f[
 * \tilde{\tau}_3 \equiv \mathrm{\texttt{tau_fluid_term3}}
 *   = \frac{B^2 \tilde{S}^2
 *   - \Omega^2}{2 \sqrt{\gamma}\left( W_\mathrm{min} + B^2\right)^2}
 * \f]
 *
 * To prevent possible [floating-point underflow](https://en.wikipedia.org/wiki/Arithmetic_underflow)
 * or division by zero errors, we check if \f$ B^2 < 1e-150 \f$. For the case
 * of very small \f$ B \f$, we compute the previous values using expressions
 * with \f$ B=0 \f$.
 *
 * ## Modifying Energy variable
 *
 * To start, we enforce that \f$ \tilde{\tau} > \tilde{\tau}_\mathrm{atm} \f$.
 * Then, we apply the limit to \f$ \tilde{\tau} \f$ with the conditional
 *
 * \f[
 * \tilde{\tau} < \frac{\psi^6 B^2}{2}
 * \f]
 *
 * If this condition is met, then we reset \f$ \tilde{\tau} \f$ to
 *
 * \f[
 * \tilde{\tau} = \tilde{\tau}_\mathrm{atm} + \frac{\psi^6 B^2}{2}
 * \f]
 *
 * ## Modifying Momentum variable
 *
 * The logic for modifying \f$ \tilde{S}_i \f$ changes depending on the size of
 * the magnetic field:
 *
 * If \f$ B^2 < \frac{P_\mathrm{atm}}{10^{-32}} \f$ and
 * \f$ \tilde{S}^2 > \tilde{\tau}(\tilde{\tau} + 2\tilde{\tilde{D}}) \f$
 * then
 *
 * \f[
 * \tilde{S}_i = \tilde{S}_i\sqrt{\frac{\tilde{\tau}(\tilde{\tau}
 *                                                   + 2\tilde{D})}{\tilde{S}^2}}
 * \f]
 *
 * ***
 *
 * Otherwise, If \f$ \sqrt{\gamma} > \psi^6_\mathrm{threshold} \f$ then
 *
 * \f[
 * \tilde{\tau}_\mathrm{min} = \tilde{\tau} + \frac{\sqrt{\gamma}B^2}{2} + \tilde{\tau}_3
 * \f]
 *
 * If \f$ \tilde{\tau}_\mathrm{min} < 1.001 \tilde{\tau}_\mathrm{atm} \f$ then
 *
 * \f[
 * \tilde{\tau}_\mathrm{min} = 1.001 \tilde{\tau}_\mathrm{atm}
 * \f]
 *
 * \f[
 * \tilde{\tau} = \tilde{\tau}_\mathrm{min} + \frac{\psi^6 B^2}{2} + \tilde{\tau}_3
 * \f]
 *
 * Finally, if \f$ \tilde{S}^2 > \tilde{\tau}_\mathrm{min}
 *                               (\tilde{\tau}_\mathrm{min} + 2\tilde{\tilde{D}}) \f$
 * then
 *
 * \f[
 * \tilde{S}_i = \tilde{S}_i\sqrt{\frac{\tilde{\tau}_\mathrm{min}(\tilde{\tau}_\mathrm{min} + 2\tilde{D})}{\tilde{S}^2}}
 * \f]
 */
void ghl_apply_conservative_limits(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_primitive_quantities *restrict prims,
      ghl_conservative_quantities *restrict cons,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  //First, prepare for the tau and stilde fixes:
  const double sdots = ghl_compute_vec2_from_vec3D(ADM_metric->gammaUU, cons->SD);

  const double B2 = ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, prims->BU);

  double BdotS, hatBdotS, half_psi6_B2, tau_fluid_term3;
  if(B2 < 1e-150) {
    BdotS = hatBdotS = half_psi6_B2 = tau_fluid_term3 = 0.0;
  } else {
    const double Bmag = sqrt(B2);
    BdotS = prims->BU[0]*cons->SD[0] + prims->BU[1]*cons->SD[1] + prims->BU[2]*cons->SD[2];
    hatBdotS = BdotS/Bmag;
    const double Wm = sqrt(SQR(hatBdotS) + SQR(cons->rho))/ADM_metric->sqrt_detgamma;
    const double Sm2 = (SQR(Wm)*sdots + SQR(BdotS)*(B2+2.0*Wm))/SQR(Wm+B2);
    const double Wmin = sqrt(Sm2 + SQR(cons->rho))/ADM_metric->sqrt_detgamma;
    half_psi6_B2 = 0.5*ADM_metric->sqrt_detgamma*B2;
    tau_fluid_term3 = (B2*sdots - SQR(BdotS))*0.5/(ADM_metric->sqrt_detgamma*SQR(Wmin+B2));
  }

  cons->tau = MAX(cons->tau, eos->tau_atm);

  //tau fix, applicable when B==0 and B!=0:
  if(cons->tau < half_psi6_B2) {
    cons->tau = eos->tau_atm + half_psi6_B2;
    diagnostics->tau_fix = true;
  }

  //Apply Stilde fix when B==0.
  if(B2 < eos->press_atm*1e-32) {
    const double rhot = 0.999999*cons->tau*(cons->tau+2.0*cons->rho);

    if(sdots > rhot) {
      const double rfactm1 = sqrt(rhot/sdots);
      cons->SD[0]*=rfactm1;
      cons->SD[1]*=rfactm1;
      cons->SD[2]*=rfactm1;
      diagnostics->Stilde_fix = true;
    }
  //Apply new Stilde fix.
  } else if(ADM_metric->sqrt_detgamma>params->psi6threshold) {
    double tau_fluid_min = cons->tau - half_psi6_B2 - tau_fluid_term3;
    if (tau_fluid_min < eos->tau_atm*1.001) {
      tau_fluid_min = eos->tau_atm*1.001;
      cons->tau = tau_fluid_min + half_psi6_B2 + tau_fluid_term3;
      diagnostics->tau_fix = true;
    }

    const double rhot = 0.999999*tau_fluid_min*(tau_fluid_min+2.0*cons->rho);

    if(sdots > rhot) {
      const double rfactm1 = sqrt(rhot/sdots);
      cons->SD[0]*=rfactm1;
      cons->SD[1]*=rfactm1;
      cons->SD[2]*=rfactm1;
      diagnostics->Stilde_fix = true;
    }
  }
}
