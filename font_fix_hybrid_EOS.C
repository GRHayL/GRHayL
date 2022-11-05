#include "con2prim_header.h"
#include "EOS_Hybrid.h"

/**********************************
 * Piecewise Polytropic EOS Patch *
 *    Font fix: function call     *
 **********************************/
int font_fix_hybrid_EOS( const eos_parameters *restrict eos,
                         const metric_quantities *restrict metric,
                         const conservative_quantities *restrict cons_undens,
                         const primitive_quantities *restrict prims,
                         CCTK_REAL &u_x, CCTK_REAL &u_y, CCTK_REAL &u_z ) {

  CCTK_REAL Bxbar = prims->Bx*ONE_OVER_SQRT_4PI;
  CCTK_REAL Bybar = prims->By*ONE_OVER_SQRT_4PI;
  CCTK_REAL Bzbar = prims->Bz*ONE_OVER_SQRT_4PI;
  CCTK_REAL Bbar_x = metric->gxx*Bxbar + metric->gxy*Bybar + metric->gxz*Bzbar;
  CCTK_REAL Bbar_y = metric->gxy*Bxbar + metric->gyy*Bybar + metric->gyz*Bzbar;
  CCTK_REAL Bbar_z = metric->gxz*Bxbar + metric->gyz*Bybar + metric->gzz*Bzbar;
  CCTK_REAL B2bar = Bxbar*Bbar_x + Bybar*Bbar_y + Bzbar*Bbar_z;
  CCTK_REAL Bbar = sqrt(B2bar);

  CCTK_REAL check_B_small = fabs(Bxbar)+fabs(Bybar)+fabs(Bzbar);
  if (check_B_small>0 && check_B_small<1.e-150) {
    // need to compute B2bar specially to prevent floating-point underflow
    CCTK_REAL Bmax = fabs(Bxbar);
    if (Bmax < fabs(Bybar)) Bmax=fabs(Bybar);
    if (Bmax < fabs(Bzbar)) Bmax=fabs(Bzbar);
    CCTK_REAL Bxtmp=Bxbar/Bmax, Bytemp=Bybar/Bmax, Bztemp=Bzbar/Bmax;
    CCTK_REAL B_xtemp=Bbar_x/Bmax, B_ytemp=Bbar_y/Bmax, B_ztemp=Bbar_z/Bmax;
    Bbar = sqrt(Bxtmp*B_xtemp + Bytemp*B_ytemp + Bztemp*B_ztemp)*Bmax;
  }
  CCTK_REAL BbardotS = Bxbar*cons_undens->S_x + Bybar*cons_undens->S_y + Bzbar*cons_undens->S_z;
  CCTK_REAL BbardotS2 = BbardotS*BbardotS;
  CCTK_REAL hatBbardotS = BbardotS/Bbar;
  if (Bbar<1.e-300) hatBbardotS = 0.0;
  CCTK_REAL Psim6 = 1.0/metric->psi6;

  // Limit hatBbardotS
  //CCTK_REAL max_gammav = 100.0;
  //CCTK_REAL rhob_max = CONSERVS[RHOSTAR]*Psim6;
  //CCTK_REAL hmax = 1.0 + gam_gamm1_kpoly*pow(rhob_max,gam1);
  //CCTK_REAL abs_hatBbardotS_max = sqrt(SQR(max_gammav)-1.0)*CONSERVS[RHOSTAR]*hmax;
  //if (fabs(hatBbardotS) > abs_hatBbardotS_max) {
  //   CCTK_REAL fac_reduce = abs_hatBbardotS_max/fabs(hatBbardotS);
  //   CCTK_REAL hatBbardotS_max = hatBbardotS*fac_reduce;
  //   CCTK_REAL Bbar_inv = 1.0/Bbar;
  //   CCTK_REAL hat_Bbar_x = Bbar_x*Bbar_inv;
  //   CCTK_REAL hat_Bbar_y = Bbar_y*Bbar_inv;
  //   CCTK_REAL hat_Bbar_z = Bbar_z*Bbar_inv;
  //   CCTK_REAL sub_fact = hatBbardotS_max - hatBbardotS;
  //   cons->S_x += sub_fact*hat_Bbar_x;
  //   cons->S_y += sub_fact*hat_Bbar_y;
  //   cons->S_z += sub_fact*hat_Bbar_z;
  //   hatBbardotS = hatBbardotS_max;
  //   BbardotS *= fac_reduce;
  //   BbardotS2 = BbardotS*BbardotS;
  //}

  CCTK_REAL sdots = metric->adm_gupxx*SQR(cons_undens->S_x) + metric->adm_gupyy*SQR(cons_undens->S_y) + metric->adm_gupzz*SQR(cons_undens->S_z)
    + 2.0*( metric->adm_gupxy*cons_undens->S_x*cons_undens->S_y + metric->adm_gupxz*cons_undens->S_x*cons_undens->S_z
            + metric->adm_gupyz*cons_undens->S_y*cons_undens->S_z);


  CCTK_REAL rhob;
  if (sdots<1.e-300) {
    rhob = cons_undens->rho*Psim6;
    u_x=0.0; u_y=0.0; u_z=0.0;
    return 0;
  }
  /* This test has some problem.
     if (fabs(BbardotS2 - sdots*B2bar) > 1e-8) {
     CCTK_VInfo(CCTK_THORNSTRING,"(Bbar dot S)^2, Bbar^2 * sdotS, %e %e",SQR(BbardotS),sdots*B2bar);
     CCTK_VInfo(CCTK_THORNSTRING,"Cauchy-Schwartz inequality is violated!");
     }
  */


  // Initial guess for W, S_fluid and rhob
  CCTK_REAL W0    = sqrt( SQR(hatBbardotS) + SQR(cons_undens->rho) ) * Psim6;
  CCTK_REAL Sf20  = (SQR(W0)*sdots + BbardotS2*(B2bar + 2.0*W0))/SQR(W0+B2bar);
  CCTK_REAL rhob0 = cons_undens->rho*Psim6/sqrt(1.0+Sf20/SQR(cons_undens->rho));


  //****************************************************************
  //                          FONT FIX
  // Impose Font fix when HARM primitives solver fails to find
  //   acceptable set of primitives.
  //****************************************************************

  /* Set the maximum number of iterations */
  int maxits = 500;

  /* Set the allowed tolerance */
  CCTK_REAL tol = 1.e-15;

  /* Declare basic variables */
  int font_fix_status;

    /**********************
   * FONT FIX MAIN LOOP *
   **********************
   * Perform the font fix routine until convergence
   * is obtained and the algorithm returns with no
   * error. Every time the Font fix fails, increase
   * the tolerance by a factor of 10.
   */
  int font_fix_attempts = 5;
  CCTK_REAL font_fix_tol_factor = 10.0;
  for(int n=0; n<font_fix_attempts; n++) {

    tol *= pow(font_fix_tol_factor,n);
    font_fix_status = font_fix__rhob_loop(maxits,tol, W0,Sf20,Psim6,sdots,BbardotS2,B2bar, cons_undens, eos, rhob0,rhob);
    rhob0 = rhob;
    if(font_fix_status==0) break;

  }


  //**************************************************************************************************************

  /* Font fix works! */
  /* First compute P_cold, eps_cold, then h = h_cold */
  CCTK_REAL P_cold, eps_cold;
  compute_P_cold__eps_cold(eos, rhob, &P_cold, &eps_cold);
  CCTK_REAL h = 1.0 + eps_cold + P_cold/rhob;

  /* Then compute gamma_v using equation (A19) in
   * Etienne et al. (2011) [https://arxiv.org/pdf/1112.0568.pdf]
   * .-----------------------------------------.
   * | gamma_v = psi^{-6} * (rho_star / rho_b) |
   * .-----------------------------------------.
   */
  CCTK_REAL gammav = cons_undens->rho*Psim6/rhob;

  /* Finally, compute u_{i} */
  CCTK_REAL rhosh = cons_undens->rho*h;
  CCTK_REAL fac1 = metric->psi6*BbardotS/(gammav*rhosh);
  CCTK_REAL fac2 = 1.0/(rhosh + metric->psi6*B2bar/gammav);
  u_x = fac2*(cons_undens->S_x + fac1*Bbar_x);
  u_y = fac2*(cons_undens->S_y + fac1*Bbar_y);
  u_z = fac2*(cons_undens->S_z + fac1*Bbar_z);

  return 0;
}
