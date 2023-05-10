#include "../../harm_u2p_util.h"

/**********************************
 * Piecewise Polytropic EOS Patch *
 *    Font fix: function call     *
 **********************************/
int Hybrid_Font_Fix(
      const GRHayL_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics) {

  diagnostics->failure_checker+=10000;
  diagnostics->font_fix=1;
  double u_x, u_y, u_z;

  const double Bxbar = prims->Bx*ONE_OVER_SQRT_4PI;
  const double Bybar = prims->By*ONE_OVER_SQRT_4PI;
  const double Bzbar = prims->Bz*ONE_OVER_SQRT_4PI;
  const double Bbar_x = metric->adm_gxx*Bxbar + metric->adm_gxy*Bybar + metric->adm_gxz*Bzbar;
  const double Bbar_y = metric->adm_gxy*Bxbar + metric->adm_gyy*Bybar + metric->adm_gyz*Bzbar;
  const double Bbar_z = metric->adm_gxz*Bxbar + metric->adm_gyz*Bybar + metric->adm_gzz*Bzbar;
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

  double sdots = metric->adm_gupxx*SQR(cons->S_x) + metric->adm_gupyy*SQR(cons->S_y) + metric->adm_gupzz*SQR(cons->S_z)
    + 2.0*( metric->adm_gupxy*cons->S_x*cons->S_y + metric->adm_gupxz*cons->S_x*cons->S_z
            + metric->adm_gupyz*cons->S_y*cons->S_z);

  double rhob;
  if (sdots<1.e-300) {
    u_x=0.0;
    u_y=0.0;
    u_z=0.0;
  } else {
    // Initial guess for W, S_fluid and rhob
    double W0    = sqrt( SQR(hatBbardotS) + SQR(cons->rho) ) * Psim6;
    double Sf20  = (SQR(W0)*sdots + BbardotS2*(B2bar + 2.0*W0))/SQR(W0+B2bar);
    double rhob0 = cons->rho*Psim6/sqrt(1.0+Sf20/SQR(cons->rho));

    //****************************************************************
    //                          FONT FIX
    // Impose Font fix when HARM primitives solver fails to find
    //   acceptable set of primitives.
    //****************************************************************

    const int maxits = 300;
    double tol = 1.e-15;
    int font_fix_status;
    const int font_fix_attempts = 5;
    for(int n=0; n<font_fix_attempts; n++) {
      const int loop_maxits = maxits + n*50; // From 300 to 500 for 5 iterations
      const double loop_tol = tol*pow(4,n); // tolerance multipliers are {0,4,16,64,256}
      font_fix_status = hybrid_font_fix_loop(eos, loop_maxits, loop_tol, W0, Sf20, Psim6, sdots,
                                           BbardotS2, B2bar, cons, rhob0, &rhob);
      rhob0 = rhob;
      if(font_fix_status==0) break;
    }

    //**************************************************************************************************************

    if(font_fix_status==1)
      return 1;

    /* Font fix works! */
    /* First compute P_cold, eps_cold, then h = h_cold */
    double P_cold, eps_cold;
    eos->hybrid_compute_P_cold_and_eps_cold(eos, rhob, &P_cold, &eps_cold);
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
    u_x = fac2*(cons->S_x + fac1*Bbar_x);
    u_y = fac2*(cons->S_y + fac1*Bbar_y);
    u_z = fac2*(cons->S_z + fac1*Bbar_z);
  }

  //Translate to HARM primitive now:
  double utcon1 = metric->adm_gupxx*u_x + metric->adm_gupxy*u_y + metric->adm_gupxz*u_z;
  double utcon2 = metric->adm_gupxy*u_x + metric->adm_gupyy*u_y + metric->adm_gupyz*u_z;
  double utcon3 = metric->adm_gupxz*u_x + metric->adm_gupyz*u_y + metric->adm_gupzz*u_z;

  //The Font fix only sets the velocities.  Here we set the pressure & density HARM primitives.
  limit_utilde_and_compute_v(eos, metric, &utcon1, &utcon2, &utcon3, prims, &diagnostics->speed_limited);
  prims->rho = cons->rho/(metric->lapse*prims->u0*metric->psi6);

  double K_ppoly, Gamma_ppoly;
  eos->hybrid_get_K_and_Gamma(eos, prims->rho, &K_ppoly, &Gamma_ppoly);

  // After that, we set P = P_cold
  eos->hybrid_compute_P_cold(eos, prims->rho, &prims->press);

  // and compute epsilon from rho and pressure
  prims->eps = prims->press/(prims->rho*(Gamma_ppoly-1.0));
  if( params->evolve_entropy ) eos->hybrid_compute_entropy_function(eos, prims->rho, prims->press, &prims->entropy);
  return 0;
}
