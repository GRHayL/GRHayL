#include "con2prim.h"
// HEADER_TODO: remove need for EOS header
#include "EOS_hybrid.h"

/**********************************
 * Piecewise Polytropic EOS Patch *
 *    Font fix: function call     *
 **********************************/
int font_fix_hybrid_EOS( const eos_parameters *restrict eos,
                         const metric_quantities *restrict metric,
                         const conservative_quantities *restrict cons,
                         const primitive_quantities *restrict prims,
                         double *restrict u_x_ptr, double *restrict u_y_ptr, double *restrict u_z_ptr ) {

  double Bxbar = prims->Bx*ONE_OVER_SQRT_4PI;
  double Bybar = prims->By*ONE_OVER_SQRT_4PI;
  double Bzbar = prims->Bz*ONE_OVER_SQRT_4PI;
  double Bbar_x = metric->adm_gxx*Bxbar + metric->adm_gxy*Bybar + metric->adm_gxz*Bzbar;
  double Bbar_y = metric->adm_gxy*Bxbar + metric->adm_gyy*Bybar + metric->adm_gyz*Bzbar;
  double Bbar_z = metric->adm_gxz*Bxbar + metric->adm_gyz*Bybar + metric->adm_gzz*Bzbar;
  double B2bar = Bxbar*Bbar_x + Bybar*Bbar_y + Bzbar*Bbar_z;
  double Bbar = sqrt(B2bar);

  double check_B_small = fabs(Bxbar)+fabs(Bybar)+fabs(Bzbar);
  if (check_B_small>0 && check_B_small<1.e-150) {
    // need to compute B2bar specially to prevent floating-point underflow
    double Bmax = fabs(Bxbar);
    if (Bmax < fabs(Bybar)) Bmax=fabs(Bybar);
    if (Bmax < fabs(Bzbar)) Bmax=fabs(Bzbar);
    double Bxtmp=Bxbar/Bmax, Bytemp=Bybar/Bmax, Bztemp=Bzbar/Bmax;
    double B_xtemp=Bbar_x/Bmax, B_ytemp=Bbar_y/Bmax, B_ztemp=Bbar_z/Bmax;
    Bbar = sqrt(Bxtmp*B_xtemp + Bytemp*B_ytemp + Bztemp*B_ztemp)*Bmax;
  }
  double BbardotS = Bxbar*cons->S_x + Bybar*cons->S_y + Bzbar*cons->S_z;
  double BbardotS2 = BbardotS*BbardotS;
  double hatBbardotS = BbardotS/Bbar;
  if (Bbar<1.e-300) hatBbardotS = 0.0;
  double Psim6 = 1.0/metric->psi6;

  // Limit hatBbardotS
  //double max_gammav = 100.0;
  //double rhob_max = CONSERVS[RHOSTAR]*Psim6;
  //double hmax = 1.0 + gam_gamm1_kpoly*pow(rhob_max,gam1);
  //double abs_hatBbardotS_max = sqrt(SQR(max_gammav)-1.0)*CONSERVS[RHOSTAR]*hmax;
  //if (fabs(hatBbardotS) > abs_hatBbardotS_max) {
  //   double fac_reduce = abs_hatBbardotS_max/fabs(hatBbardotS);
  //   double hatBbardotS_max = hatBbardotS*fac_reduce;
  //   double Bbar_inv = 1.0/Bbar;
  //   double hat_Bbar_x = Bbar_x*Bbar_inv;
  //   double hat_Bbar_y = Bbar_y*Bbar_inv;
  //   double hat_Bbar_z = Bbar_z*Bbar_inv;
  //   double sub_fact = hatBbardotS_max - hatBbardotS;
  //   cons->S_x += sub_fact*hat_Bbar_x;
  //   cons->S_y += sub_fact*hat_Bbar_y;
  //   cons->S_z += sub_fact*hat_Bbar_z;
  //   hatBbardotS = hatBbardotS_max;
  //   BbardotS *= fac_reduce;
  //   BbardotS2 = BbardotS*BbardotS;
  //}

  double sdots = metric->adm_gupxx*SQR(cons->S_x) + metric->adm_gupyy*SQR(cons->S_y) + metric->adm_gupzz*SQR(cons->S_z)
    + 2.0*( metric->adm_gupxy*cons->S_x*cons->S_y + metric->adm_gupxz*cons->S_x*cons->S_z
            + metric->adm_gupyz*cons->S_y*cons->S_z);


  double rhob;
  if (sdots<1.e-300) {
    rhob = cons->rho*Psim6;
    *u_x_ptr=0.0; *u_y_ptr=0.0; *u_z_ptr=0.0;
    return 0;
  }
  /* This test has some problem.
     if (fabs(BbardotS2 - sdots*B2bar) > 1e-8) {
     CCTK_VINFO("(Bbar dot S)^2, Bbar^2 * sdotS, %e %e",SQR(BbardotS),sdots*B2bar);
     CCTK_VINFO("Cauchy-Schwartz inequality is violated!");
     }
  */


  // Initial guess for W, S_fluid and rhob
  double W0    = sqrt( SQR(hatBbardotS) + SQR(cons->rho) ) * Psim6;
  double Sf20  = (SQR(W0)*sdots + BbardotS2*(B2bar + 2.0*W0))/SQR(W0+B2bar);
  double rhob0 = cons->rho*Psim6/sqrt(1.0+Sf20/SQR(cons->rho));


  //****************************************************************
  //                          FONT FIX
  // Impose Font fix when HARM primitives solver fails to find
  //   acceptable set of primitives.
  //****************************************************************

  /* Set the maximum number of iterations */
  int maxits = 500;

  /* Set the allowed tolerance */
  double tol = 1.e-15;

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
  double font_fix_tol_factor = 10.0;
  for(int n=0; n<font_fix_attempts; n++) {

    tol *= pow(font_fix_tol_factor,n);
    font_fix_status = font_fix_rhob_loop(eos, maxits, tol, W0, Sf20, Psim6, sdots,
                                         BbardotS2, B2bar, cons, rhob0, &rhob);
    rhob0 = rhob;
    if(font_fix_status==0) break;

  }


  //**************************************************************************************************************

  /* Font fix works! */
  /* First compute P_cold, eps_cold, then h = h_cold */
  double P_cold, eps_cold;
  (*eos->hybrid_compute_P_cold_and_eps_cold)(eos, rhob, &P_cold, &eps_cold);
  double h = 1.0 + eps_cold + P_cold/rhob;

  /* Then compute gamma_v using equation (A19) in
   * Etienne et al. (2011) [https://arxiv.org/pdf/1112.0568.pdf]
   * .-----------------------------------------.
   * | gamma_v = psi^{-6} * (rho_star / rho_b) |
   * .-----------------------------------------.
   */
  double gammav = cons->rho*Psim6/rhob;

  /* Finally, compute u_{i} */
  double rhosh = cons->rho*h;
  double fac1 = metric->psi6*BbardotS/(gammav*rhosh);
  double fac2 = 1.0/(rhosh + metric->psi6*B2bar/gammav);
  *u_x_ptr = fac2*(cons->S_x + fac1*Bbar_x);
  *u_y_ptr = fac2*(cons->S_y + fac1*Bbar_y);
  *u_z_ptr = fac2*(cons->S_z + fac1*Bbar_z);

  return 0;
}
