#include "../../utils_Palenzuela1D.h"

static ghl_error_codes_t ghl_newman_energy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const double S_squared,
      const double BdotS,
      const double B_squared,
      const double *restrict SU,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_conservative_quantities *restrict con,
      ghl_primitive_quantities *restrict prims,
      const double tol_x,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  // Set basic quantities from input
  double invD   = 1.0/con->rho;
  double xye    = con->Y_e * invD;
  double xtemp  = prims->temperature;

  //Now, begin iterative procedure to derive the primitive variables
  unsigned step=0;
  const unsigned maxsteps = 300;
  double AtP[maxsteps]; // not length 3 in case of extrap. probs
  double AtR;
  int AtStep=0;

  // Floor tau
  const double tau = MAX(con->tau, eos->tau_atm);

  // d = 0.5( S^{2}*B^{2} - (B.S)^{2} ) (eq. 5.7 in Newman & Hamlin 2014)
  double d = 0.5*(S_squared*B_squared-BdotS*BdotS);
  if(d < 1e-20) d = 0.0;

  // e = tau + D
  const double e = tau + con->rho;
  // z = rho*h*W^{2} = D*h*W; initialize to zero.
  double z=0, invW=0, W=0;

  // Useful auxiliary variable
  const double BdotSsq = BdotS * BdotS;

  // Now initialize the pressure and AtP
  double P_old, xprs = 0.0;
  AtP[0] = xprs;

  // Leo mod: compute auxiliary variables so we can find eps
  const double q = tau * invD;
  const double s = B_squared * invD;
  const double t = BdotS/(pow(con->rho,1.5));

  bool conacc = false;
  do {
    P_old = xprs;
    step++;
    //This check is required to ensure that |cos(cPhi)|<1 below

    // a = e + P + B^{2}/2 (eq. 5.7 in Newman & Hamlin 2014)
    double a = e + xprs + 0.5*B_squared;

    // Eq. (5.9) of Newman & Hamlin 2014
    if(d > 4.0*a*a*a/27.0) return 1;

    // phi = acos( sqrt(27d/4a^{3}) ) (eq. 5.10 in Newman & Hamlin 2014)
    double phi = acos( sqrt(27.0*d/(4.0*a*a*a)) );
    // Eps = a/3 - 2a/3 cos( 2phi/3 + 2pi/3 )
    double Eps = (a/3.0)*(1.0 - 2.0*cos((2.0/3.0)*(phi + M_PI)));

    // From the definition of Eps: Eps = z + B^{2} => z = Eps - B^{2}
    z = Eps-B_squared;

    // Now compute (eq. 5.2 in Newman & Hamlin 2014)
    //
    // v^{2} = ( z^{2}S^{2} + (2z + B^{2})(B.S)^{2}/(z^{2}(z+B^{2})^{2})
    const double Epssq = Eps*Eps;
    const double zsq   = z*z;
    const double vsq   = (zsq * S_squared + (z+Eps)*BdotSsq)/(zsq*Epssq);

    // Impose physical limits and compute W
    invW = MIN(MAX(sqrt(1.0-vsq), 1.0/params->max_Lorentz_factor), 1.0);
    W    = 1.0/invW;

    // Set the prims (eps is computed as in Palenzuela et al.)
    // See e.g., Eq. (44) of https://arxiv.org/pdf/1712.07538.pdf
    double x    = z * invD;
    double xrho = con->rho * invW;
    double xeps = - 1.0 + (1.0-W*W)*x*invW
                  + W*( 1.0 + q - s + 0.5*( s*invW*invW + (t*t)/(x*x) ) );
    ghl_tabulated_enforce_bounds_rho_Ye_eps(eos, &xrho, &xye, &xeps);
    ghl_tabulated_compute_P_T_from_eps(eos, xrho, xye, xeps, &xprs, &xtemp);

    AtStep++;
    AtP[AtStep]=xprs;

    if(AtStep>=2) {   //Aitken extrapolation
      AtR = (AtP[AtStep]-AtP[AtStep-1])/(AtP[AtStep-1]-AtP[AtStep-2]);
      if(AtR<1. && AtR>0.) {
        xprs=AtP[AtStep-1]+(AtP[AtStep]-AtP[AtStep-1])/(1.-AtR);
        AtStep=0;
        conacc = 1;
        AtP[0]=xprs;   //starting value for next Aitken extrapolation
      }
    }
  }
  while(fabs(xprs-P_old)>tol_x*(xprs+P_old) && step<maxsteps);

  if(step >= maxsteps)
    return ghl_error_c2p_max_iter;

  diagnostics->n_iter = step;

  if(conacc) {     //converged on an extrap. so recompute vars
    const double a     = e + xprs + 0.5*B_squared;
    const double phi   = acos(sqrt(27.0*d/(4.0*a))/a);
    const double Eps   = a/3.0*( 1.0 - 2.0*cos( (2.0/3.0)*(phi + M_PI) ) );
    z                  = Eps-B_squared;
    const double Epssq = Eps*Eps;
    const double zsq   = z*z;
    const double vsq   = (zsq * S_squared + (z+Eps)*BdotSsq)/(zsq*Epssq);
    invW               = MIN(MAX(sqrt(1.0-vsq), 1.0/params->max_Lorentz_factor), 1.0);
    W                  = 1.0/invW;
  }

  if(isnan(z*W)) return ghl_error_c2p_singular;

  // Set the primitives
  double utildeU[3] = {
    W*(SU[0] + BdotS*prims->BU[0]/z)/(z+B_squared),
    W*(SU[1] + BdotS*prims->BU[1]/z)/(z+B_squared),
    W*(SU[2] + BdotS*prims->BU[2]/z)/(z+B_squared)
  };
  prims->rho         = con->rho*invW;
  prims->Y_e         = xye;
  prims->temperature = xtemp;
  ghl_tabulated_enforce_bounds_rho_Ye_T(eos, &prims->rho, &prims->Y_e, &prims->temperature);
  diagnostics->speed_limited = ghl_limit_utilde_and_compute_v(params, ADM_metric, utildeU, prims);
  ghl_tabulated_compute_P_eps_S_from_T(eos, prims->rho, prims->Y_e, prims->temperature,
                                       &prims->press, &prims->eps, &prims->entropy);

  return ghl_success;
}

ghl_error_codes_t ghl_tabulated_Newman1D_energy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  // Step 1: Compute auxiliary quantities
  double SU[3], Bsq, Ssq, BdotS;
  compute_SU_Bsq_Ssq_BdotS(ADM_metric, cons_undens, prims,
                           SU, &Bsq, &Ssq, &BdotS);

  // Step 2: Call the Newman routine that uses the energy to recover T
  const double tol_x = 1e-15;
  diagnostics->which_routine = Newman1D;

  ghl_error_codes_t error = ghl_newman_energy(params, eos, Ssq, BdotS, Bsq, SU, ADM_metric,
				cons_undens, prims, tol_x, diagnostics);

  if(error) {
    prims->temperature = eos->T_min;
    error = ghl_newman_energy(params, eos, Ssq, BdotS, Bsq, SU, ADM_metric,
			      cons_undens, prims, tol_x, diagnostics);
  }
  return error;
}
