#include "con2prim.h"

/*
 * Function     : ghl_apply_conservative_limits()
 * Description  : Checks whether \tilde{tau} and \tilde{S} obey the inequalities
 *                constraining their possible values. If they are outside
 *                these bounds, the function enforces the inequalities
 *                and updates the diagnostic data.
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_apply_conservative_limits
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
