#include "../../harm_u2p_util.h"

/**********************************
 * Piecewise Polytropic EOS Patch *
 *    Font fix: function call     *
 **********************************/
int ghl_hybrid_Font_fix(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  double utU[3];

  const double sdots = ghl_compute_vec2_from_vec3D(ADM_metric->gammaUU, cons->SD);

  const double BbarU[3] = {prims->BU[0]*ONE_OVER_SQRT_4PI, prims->BU[1]*ONE_OVER_SQRT_4PI, prims->BU[2]*ONE_OVER_SQRT_4PI};
  const double Bbar2 = ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, BbarU);

  double BbardotS, BbardotS2, hatBbardotS;
  if(Bbar2 < 1e-150) {
    BbardotS = BbardotS2 = hatBbardotS = 0.0;
  } else {
    const double Bbar_mag = sqrt(Bbar2);
    BbardotS = BbarU[0]*cons->SD[0] + BbarU[1]*cons->SD[1] + BbarU[2]*cons->SD[2];
    BbardotS2 = BbardotS*BbardotS;
    hatBbardotS = BbardotS/Bbar_mag;
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
    // Initial guess for W, S_fluid and rhob
    double W0    = sqrt( SQR(hatBbardotS) + SQR(cons->rho) ) * Psim6;
    double Sf20  = (SQR(W0)*sdots + BbardotS2*(Bbar2 + 2.0*W0))/SQR(W0+Bbar2);
    double rhob0 = cons->rho*Psim6/sqrt(1.0+Sf20/SQR(cons->rho));

    //****************************************************************
    //                          FONT FIX
    //****************************************************************

    const int maxits = 300;
    double tol = 1.e-15;
    int Font_fix_status;
    const int Font_fix_attempts = 5;
    for(int n=0; n<Font_fix_attempts; n++) {
      const int loop_maxits = maxits + n*50; // From 300 to 500 for 5 iterations
      const double loop_tol = tol*pow(4,n); // tolerance multipliers are {0,4,16,64,256}
      Font_fix_status = ghl_hybrid_Font_fix_loop(eos, loop_maxits, loop_tol, W0, Sf20, Psim6, sdots,
                                           BbardotS2, Bbar2, cons, rhob0, &rhob);
      rhob0 = rhob;
      if(Font_fix_status==0) break;
    }

    //****************************************************************

    if(Font_fix_status==1)
      return 1;

    /* First compute P_cold, eps_cold, then h = h_cold */
    double P_cold, eps_cold;
    ghl_hybrid_compute_P_cold_and_eps_cold(eos, rhob, &P_cold, &eps_cold);
    double h = 1.0 + eps_cold + P_cold/rhob;

    /* Then compute gamma_v using equation (A19) in
     * Etienne et al. (2011) [https://arxiv.org/pdf/1112.0568.pdf]
     * .-----------------------------------------.
     * | gamma_v = psi^{-6} * (rho_star / rho_b) |
     * .-----------------------------------------.
     */
    double gammav = cons->rho*Psim6/rhob;

    /* Finally, compute u^{i} */
    double rhosh = cons->rho*h;
    double fac1 = ADM_metric->sqrt_detgamma*BbardotS/(gammav*rhosh);
    double fac2 = 1.0/(rhosh + ADM_metric->sqrt_detgamma*Bbar2/gammav);

    double SU[3];
    ghl_raise_lower_vector_3D(ADM_metric->gammaUU, cons->SD, SU);

    utU[0] = fac2*(SU[0] + fac1*BbarU[0]);
    utU[1] = fac2*(SU[1] + fac1*BbarU[1]);
    utU[2] = fac2*(SU[2] + fac1*BbarU[2]);
    diagnostics->speed_limited = ghl_limit_utilde_and_compute_v(eos, ADM_metric, utU, prims);
  }

  //The Font fix only sets the velocities.  Here we set the primitives.
  prims->rho = cons->rho/(ADM_metric->lapse*prims->u0*ADM_metric->sqrt_detgamma);

  double K_ppoly, Gamma_ppoly;
  ghl_hybrid_get_K_and_Gamma(eos, prims->rho, &K_ppoly, &Gamma_ppoly);

  // After that, we set P = P_cold
  ghl_hybrid_compute_P_cold(eos, prims->rho, &prims->press);

  // and compute epsilon from rho and pressure
  prims->eps = prims->press/(prims->rho*(Gamma_ppoly-1.0));
  if( params->evolve_entropy ) ghl_hybrid_compute_entropy_function(eos, prims->rho, prims->press, &prims->entropy);

  /* Font fix works! */
  diagnostics->which_routine = Font1D;
  return 0;
}
