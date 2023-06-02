#include "con2prim.h"

/* Function    : grhayl_apply_inequality_fixes()
 * Description : Checks whether \tilde{tau} \tilde{S} obey the inequalities
 *               constraining their possible values. If they are outside
 *               these bounds, the function enforces the inequalities
 *               and updates the diagnostic data.
 *
 * Inputs      : params         - grhayl_parameters struct with parameters
 *                                for the simulation
 *             : eos            - eos_parameters struct with data for the
 *                                EOS of the simulation
 *             : metric         - metric_quantities struct with data for
 *                                the gridpoint of interest
 *             : prims          - primitive_quantities struct with data
 *                                for the gridpoint of interest
 *
 * Outputs     : cons           - returns potentially different values for cons->tau and cons->SD
 *             : diagnostics    - tracks whether any inequalities were violated
 *
 */

void grhayl_apply_inequality_fixes(
      const grhayl_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const primitive_quantities *restrict prims,
      conservative_quantities *restrict cons,
      con2prim_diagnostics *restrict diagnostics) {

  //First, prepare for the tau and stilde fixes:
  const double sdots = grhayl_compute_vec2_from_vecD(ADM_metric->gammaUU, cons->SD);

  const double Bbar[3] = {prims->BU[0]*ONE_OVER_SQRT_4PI, prims->BU[1]*ONE_OVER_SQRT_4PI, prims->BU[2]*ONE_OVER_SQRT_4PI};

  const double Bbar2 = grhayl_compute_vec2_from_vecU(ADM_metric->gammaDD, Bbar);

  double BbardotS, hatBbardotS;
  double Wm, Sm2, Wmin, half_psi6_Bbar2, tau_fluid_term3;
  if(Bbar2 < 1e-150) {
    BbardotS = hatBbardotS = half_psi6_Bbar2 = tau_fluid_term3 = 0.0;
    Wm = cons->rho/metric_aux->psi6;
    Sm2 = (Wm*sdots + 2.0)/Wm;
    Wmin = sqrt(Sm2 + SQR(cons->rho))/metric_aux->psi6;
  } else {
    const double Bbar_mag = sqrt(Bbar2);
    BbardotS = Bbar[0]*cons->SD[0] + Bbar[1]*cons->SD[1] + Bbar[2]*cons->SD[2];
    hatBbardotS = BbardotS/Bbar_mag;
    Wm = sqrt(SQR(hatBbardotS)+ SQR(cons->rho))/metric_aux->psi6;
    Sm2 = (SQR(Wm)*sdots + SQR(BbardotS)*(Bbar2+2.0*Wm))/SQR(Wm+Bbar2);
    Wmin = sqrt(Sm2 + SQR(cons->rho))/metric_aux->psi6;
    half_psi6_Bbar2 = 0.5*metric_aux->psi6*Bbar2;
    tau_fluid_term3 = (Bbar2*sdots - SQR(BbardotS))*0.5/(metric_aux->psi6*SQR(Wmin+Bbar2));
  }

  //tau fix, applicable when B==0 and B!=0:
  if(cons->tau < half_psi6_Bbar2) {
    cons->tau = eos->tau_atm + half_psi6_Bbar2;
    diagnostics->failure_checker+=1000000;
  }

  double tau_fluid_min = cons->tau - half_psi6_Bbar2 - tau_fluid_term3;

  //Apply Stilde fix when B==0.
  if(Bbar2 < eos->press_atm*1e-32) {
    const double rhot = 0.999999*cons->tau*(cons->tau+2.0*cons->rho);

    if(sdots > rhot) {
      const double rfactm1 = sqrt(rhot/sdots);
      cons->SD[0]*=rfactm1;
      cons->SD[1]*=rfactm1;
      cons->SD[2]*=rfactm1;
      diagnostics->failure_checker+=10000000;
    }
  //Apply new Stilde fix.
  } else if(metric_aux->psi6>params->psi6threshold) {
    if (tau_fluid_min < eos->tau_atm*1.001) {
      tau_fluid_min = eos->tau_atm*1.001;
      cons->tau = tau_fluid_min + half_psi6_Bbar2 + tau_fluid_term3;
    }

    const double rhot = 0.999999*tau_fluid_min*(tau_fluid_min+2.0*cons->rho);

    if(sdots > rhot) {
      const double rfactm1 = sqrt(rhot/sdots);
      cons->SD[0]*=rfactm1;
      cons->SD[1]*=rfactm1;
      cons->SD[2]*=rfactm1;
      diagnostics->failure_checker+=100000000;
    }
  }
}
