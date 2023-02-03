#include "con2prim.h"

/* Function    : apply_inequality_fixes()
 * Description : This function checks whether tau tilde and S tilde
 *               obey the inequalities constraining their possible values.
 *               If they are outside these bounds, the function enforces
 *               the inequalities and updates the diagnostic data.
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
 * Outputs     : cons           - conservative_quantities struct with potentially
 *                                different values for cons->tau and cons->S_*
 *             : diagnostics    - con2prim_diagnostics struct which tracks whether
 *                                any of the inequality conditions were violated
 *                                and required the conservatives to be adjusted
 */

void apply_inequality_fixes(
      const GRHayL_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      const primitive_quantities *restrict prims,
      conservative_quantities *restrict cons,
      con2prim_diagnostics *restrict diagnostics) {

  //First, prepare for the tau and stilde fixes:

  const double Bxbar = prims->Bx*ONE_OVER_SQRT_4PI;
  const double Bybar = prims->By*ONE_OVER_SQRT_4PI;
  const double Bzbar = prims->Bz*ONE_OVER_SQRT_4PI;

  const double Bbar_x = metric->adm_gxx*Bxbar + metric->adm_gxy*Bybar + metric->adm_gxz*Bzbar;
  const double Bbar_y = metric->adm_gxy*Bxbar + metric->adm_gyy*Bybar + metric->adm_gyz*Bzbar;
  const double Bbar_z = metric->adm_gxz*Bxbar + metric->adm_gyz*Bybar + metric->adm_gzz*Bzbar;
  const double Bbar2 = Bxbar*Bbar_x + Bybar*Bbar_y + Bzbar*Bbar_z;
  double Bbar = sqrt(Bbar2);

  const double check_B_small = fabs(Bxbar)+fabs(Bybar)+fabs(Bzbar);
  if (check_B_small>0 && check_B_small<1.e-150) {
    // need to compute Bbar specially to prevent floating-point underflow
    double Bmax = fabs(Bxbar);
    if (Bmax < fabs(Bybar)) Bmax = fabs(Bybar);
    if (Bmax < fabs(Bzbar)) Bmax = fabs(Bzbar);
    const double Bxtmp = Bxbar/Bmax, Bytemp=Bybar/Bmax, Bztemp=Bzbar/Bmax;
    const double B_xtemp = Bbar_x/Bmax, B_ytemp=Bbar_y/Bmax, B_ztemp=Bbar_z/Bmax;
    Bbar = sqrt(Bxtmp*B_xtemp + Bytemp*B_ytemp + Bztemp*B_ztemp)*Bmax;
  }

  double BbardotS = Bxbar*cons->S_x + Bybar*cons->S_y + Bzbar*cons->S_z;
  double hatBbardotS = BbardotS/Bbar;
  if (Bbar<1.e-300) hatBbardotS = 0.0;

  // Limit hatBbardotS
  //double max_gammav = 100.0;
  //double rhob_max = CONSERVS[RHOSTAR]/METRIC_LAP_PSI4[PSI6];
  //double hmax = 1.0 + 2.0*rhob_max;
  //double abs_hatBbardotS_max = sqrt(SQR(max_gammav)-1.0)*CONSERVS[RHOSTAR]*hmax;
  //if (fabs(hatBbardotS) > abs_hatBbardotS_max) {
  //   double fac_reduce = abs_hatBbardotS_max/fabs(hatBbardotS);
  //   double hatBbardotS_max = hatBbardotS*fac_reduce;
  //   double Bbar_inv = 1.0/Bbar;
  //   double hat_Bbar_x = Bbar_x*Bbar_inv;
  //   double hat_Bbar_y = Bbar_y*Bbar_inv;
  //   double hat_Bbar_z = Bbar_z*Bbar_inv;
  //   double sub_fact = hatBbardotS_max - hatBbardotS;
  //   CONSERVS[STILDEX] += sub_fact*hat_Bbar_x;
  //   CONSERVS[STILDEY] += sub_fact*hat_Bbar_y;
  //   CONSERVS[STILDEZ] += sub_fact*hat_Bbar_z;
  //   hatBbardotS = hatBbardotS_max;
  //   BbardotS *= fac_reduce;
  //   CONSERVS[STILDEX] = CONSERVS[STILDEX]; CONSERVS[STILDEY] = CONSERVS[STILDEY]; CONSERVS[STILDEZ] = CONSERVS[STILDEZ];
  //}

  const double sdots= metric->adm_gupxx*SQR(cons->S_x)+metric->adm_gupyy*SQR(cons->S_y)+metric->adm_gupzz*SQR(cons->S_z)+2.0*
    (metric->adm_gupxy*cons->S_x*cons->S_y+metric->adm_gupxz*cons->S_x*cons->S_z+metric->adm_gupyz*cons->S_y*cons->S_z);

  const double Wm = sqrt(SQR(hatBbardotS)+ SQR(cons->rho))/metric->psi6;
  const double Sm2 = (SQR(Wm)*sdots + SQR(BbardotS)*(Bbar2+2.0*Wm))/SQR(Wm+Bbar2);
  const double Wmin = sqrt(Sm2 + SQR(cons->rho))/metric->psi6;
  const double sdots_fluid_max = sdots;

  //tau fix, applicable when B==0 and B!=0:
  if(cons->tau < 0.5*metric->psi6*Bbar2) {
    cons->tau = eos->tau_atm+0.5*metric->psi6*Bbar2;
    diagnostics->failure_checker+=1000000;
  }

  double tau_fluid_min = cons->tau - 0.5*metric->psi6*Bbar2 - (Bbar2*sdots - SQR(BbardotS))*0.5/(metric->psi6*SQR(Wmin+Bbar2));

  //Apply Stilde fix when B==0.
  if(check_B_small*check_B_small < eos->press_atm*1e-32) {
    const double rhot = cons->tau*(cons->tau+2.0*cons->rho);
    const double safetyfactor = 0.999999;
    //if(METRIC_LAP_PSI4[PSI6]>Psi6threshold) safetyfactor=0.99;

    if(sdots > safetyfactor*rhot) {
      const double rfactm1 = sqrt((safetyfactor*rhot)/sdots);
      cons->S_x*=rfactm1;
      cons->S_y*=rfactm1;
      cons->S_z*=rfactm1;
      diagnostics->failure_checker+=10000000;
    }
  //Apply new Stilde fix.
  } else if(metric->psi6>params->psi6threshold) {
    if (tau_fluid_min < eos->tau_atm*1.001) {
      tau_fluid_min = eos->tau_atm*1.001;
      cons->tau = tau_fluid_min + 0.5*metric->psi6*Bbar2 + (Bbar2*sdots - SQR(BbardotS))*0.5/(metric->psi6*SQR(Wmin+Bbar2));
    }

    const double LHS = tau_fluid_min*(tau_fluid_min+2.0*cons->rho);
    const double RHS = sdots_fluid_max;

    const double safetyfactor = 0.999999;
    if(safetyfactor*LHS < RHS) {
      const double rfactm1 = sqrt((safetyfactor*LHS)/RHS);
      cons->S_x*=rfactm1;
      cons->S_y*=rfactm1;
      cons->S_z*=rfactm1;
      diagnostics->failure_checker+=100000000;
    }
  }
  return;
}
