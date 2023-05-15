#include "con2prim.h"

/* Function    : apply_inequality_fixes()
 * Description : Checks whether \tilde{tau} \tilde{S} obey the inequalities
 *               constraining their possible values. If they are outside
 *               these bounds, the function enforces the inequalities
 *               and updates the diagnostic data.
 *
 * Inputs      : params         - GRHayL_parameters struct with parameters
 *                                for the simulation
 *             : eos            - eos_parameters struct with data for the
 *                                EOS of the simulation
 *             : metric         - metric_quantities struct with data for
 *                                the gridpoint of interest
 *             : prims          - primitive_quantities struct with data
 *                                for the gridpoint of interest
 *
 * Outputs     : cons           - returns potentially different values for cons->tau and cons->S_*
 *             : diagnostics    - tracks whether any inequalities were violated
 *
 */

void apply_inequality_fixes(
      const GRHayL_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      const primitive_quantities *restrict prims,
      conservative_quantities *restrict cons,
      con2prim_diagnostics *restrict diagnostics) {

  //First, prepare for the tau and stilde fixes:
  const double S[3] = {cons->S_x, cons->S_y, cons->S_z};
  const double sdots = compute_vec2_from_vcov(metric, S);

  const double Bbar[3] = {prims->Bx*ONE_OVER_SQRT_4PI, prims->By*ONE_OVER_SQRT_4PI, prims->Bz*ONE_OVER_SQRT_4PI};

  const double Bbar2 = compute_vec2_from_vcon(metric, Bbar);

  double BbardotS, hatBbardotS;
  double Wm, Sm2, Wmin, half_psi6_Bbar2, tau_fluid_term3;
  if(Bbar2 < 1e-150) {
    BbardotS = hatBbardotS = half_psi6_Bbar2 = tau_fluid_term3 = 0.0;
    Wm = cons->rho/metric->psi6;
    Sm2 = (Wm*sdots + 2.0)/Wm;
    Wmin = sqrt(Sm2 + SQR(cons->rho))/metric->psi6;
  } else {
    const double Bbar_mag = sqrt(Bbar2);
    BbardotS = Bbar[0]*cons->S_x + Bbar[1]*cons->S_y + Bbar[2]*cons->S_z;
    hatBbardotS = BbardotS/Bbar_mag;
    Wm = sqrt(SQR(hatBbardotS)+ SQR(cons->rho))/metric->psi6;
    Sm2 = (SQR(Wm)*sdots + SQR(BbardotS)*(Bbar2+2.0*Wm))/SQR(Wm+Bbar2);
    Wmin = sqrt(Sm2 + SQR(cons->rho))/metric->psi6;
    half_psi6_Bbar2 = 0.5*metric->psi6*Bbar2;
    tau_fluid_term3 = (Bbar2*sdots - SQR(BbardotS))*0.5/(metric->psi6*SQR(Wmin+Bbar2));
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
      cons->S_x*=rfactm1;
      cons->S_y*=rfactm1;
      cons->S_z*=rfactm1;
      diagnostics->failure_checker+=10000000;
    }
  //Apply new Stilde fix.
  } else if(metric->psi6>params->psi6threshold) {
    if (tau_fluid_min < eos->tau_atm*1.001) {
      tau_fluid_min = eos->tau_atm*1.001;
      cons->tau = tau_fluid_min + half_psi6_Bbar2 + tau_fluid_term3;
    }

    const double rhot = 0.999999*tau_fluid_min*(tau_fluid_min+2.0*cons->rho);

    if(sdots > rhot) {
      const double rfactm1 = sqrt(rhot/sdots);
      cons->S_x*=rfactm1;
      cons->S_y*=rfactm1;
      cons->S_z*=rfactm1;
      diagnostics->failure_checker+=100000000;
    }
  }
}
