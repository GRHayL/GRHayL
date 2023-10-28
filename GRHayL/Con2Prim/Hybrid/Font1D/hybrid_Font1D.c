#include "con2prim.h"

int ghl_hybrid_Font1D_loop(
      const ghl_eos_parameters *restrict eos,
      const int maxits, const double tol, const double W_in,
      const double Sf2_in, const double Psim6, const double sdots,
      const double BdotS2, const double B2,
      const ghl_conservative_quantities *restrict cons,
      const double rhob_in, double *restrict rhob_out_ptr);

/*
 * Function      : ghl_hybrid_Font1D()
 * Description   : Determines rhob using the Font et al prescription
 * Documentation : https://github.com/GRHayL/GRHayL/wiki/ghl_hybrid_Font1D
*/
int ghl_hybrid_Font1D(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  double utU[3];

  const double sdots = ghl_compute_vec2_from_vec3D(ADM_metric->gammaUU, cons->SD);

  const double B2 = ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, prims->BU);

  double BdotS, BdotS2, hatBdotS;
  if(B2 < 1e-150) {
    BdotS = BdotS2 = hatBdotS = 0.0;
  } else {
    const double B_mag = sqrt(B2);
    BdotS = prims->BU[0]*cons->SD[0] + prims->BU[1]*cons->SD[1] + prims->BU[2]*cons->SD[2];
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
    // Initial guess for W, S_fluid and rhob
    double W0    = sqrt( SQR(hatBdotS) + SQR(cons->rho) ) * Psim6;
    double Sf20  = (SQR(W0)*sdots + BdotS2*(B2 + 2.0*W0))/SQR(W0+B2);
    double rhob0 = cons->rho*Psim6/sqrt(1.0+Sf20/SQR(cons->rho));

    //****************************************************************
    //                          FONT FIX
    //****************************************************************

    const int maxits = 300;
    double tol = 1.e-15;
    int failure;
    const int Font1D_attempts = 5;
    for(int n=0; n<Font1D_attempts; n++) {
      const int loop_maxits = maxits + n*50; // From 300 to 500 for 5 iterations
      const double loop_tol = tol*pow(4,n); // tolerance multipliers are {0,4,16,64,256}
      failure = ghl_hybrid_Font1D_loop(
            eos, loop_maxits, loop_tol, W0, Sf20, Psim6,
            sdots, BdotS2, B2, cons, rhob0, &rhob);
      rhob0 = rhob;
      if(!failure) break;
    }

    //****************************************************************

    if(failure)
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
    double fac1 = ADM_metric->sqrt_detgamma*BdotS/(gammav*rhosh);
    double fac2 = 1.0/(rhosh + ADM_metric->sqrt_detgamma*B2/gammav);

    double SU[3];
    ghl_raise_lower_vector_3D(ADM_metric->gammaUU, cons->SD, SU);

    utU[0] = fac2*(SU[0] + fac1*prims->BU[0]);
    utU[1] = fac2*(SU[1] + fac1*prims->BU[1]);
    utU[2] = fac2*(SU[2] + fac1*prims->BU[2]);
    diagnostics->speed_limited = ghl_limit_utilde_and_compute_v(params, ADM_metric, utU, prims);
  }

  //The Font fix only sets the velocities.  Here we set the primitives.
  prims->rho = cons->rho/(ADM_metric->lapse*prims->u0*ADM_metric->sqrt_detgamma);

  double K_ppoly, Gamma_ppoly;
  ghl_hybrid_get_K_and_Gamma(eos, prims->rho, &K_ppoly, &Gamma_ppoly);

  // After that, we set P = P_cold
  ghl_hybrid_compute_P_cold(eos, prims->rho, &prims->press);

  // and compute epsilon from rho and pressure
  prims->eps = prims->press/(prims->rho*(Gamma_ppoly-1.0));
  if(params->evolve_entropy)
    prims->entropy = ghl_hybrid_compute_entropy_function(eos, prims->rho, prims->press);

  /* Font fix works! */
  diagnostics->which_routine = Font1D;
  return 0;
}
