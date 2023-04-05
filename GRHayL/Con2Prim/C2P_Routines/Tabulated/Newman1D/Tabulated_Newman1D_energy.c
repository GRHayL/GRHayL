#include "../utils.h"

int newman_energy(
      const eos_parameters *restrict eos,
      const double S_squared,
      const double BdotS,
      const double B_squared,
      const double *restrict BU,
      const double *restrict SU,
      const conservative_quantities *restrict con,
      primitive_quantities *restrict prim,
      const double tol_x ) {

  // Set basic quantities from input
  double invD   = 1.0/con->rho;
  double xye    = con->Y_e * invD;
  double xtemp  = eos->T_atm;

  //Now, begin iterative procedure to derive the primitive variables
  unsigned step=0;
  const unsigned maxsteps = 300;
  double AtP[maxsteps]; // not length 3 in case of extrap. probs
  double AtR;
  int AtStep=0;

  // d = 0.5( S^{2}*B^{2} - (B.S)^{2} ) (eq. 5.7 in Newman & Hamlin 2014)
  const double d = MAX(0.5*(S_squared*B_squared-BdotS*BdotS),0.0);
  // e = tau + D
  const double e = con->tau + con->rho;
  // z = rho*h*W^{2} = D*h*W; initialize to zero.
  double z=0, invW=0, W=0;

  // Useful auxiliary variable
  const double BdotSsq = BdotS * BdotS;

  // Now initialize the pressure and AtP
  double P_old, xprs = 0.0;
  AtP[0] = xprs;

  // Leo mod: compute auxiliary variables so we can find eps
  const double q = con->tau * invD;
  const double s = B_squared * invD;
  const double t = BdotS/(pow(con->rho,1.5));

  bool conacc = false;
  do {
    P_old = xprs;
    step++;
    //This check is required to ensure that |cos(cPhi)|<1 below

    // a = e + P + B^{2}/2 (eq. 5.7 in Newman & Hamlin 2014)
    double a = e + xprs + 0.5*B_squared;

    // phi = acos( sqrt(27d/4a^{3}) ) (eq. 5.10 in Newman & Hamlin 2014)
    double phi = acos( sqrt(27.0*d/(4.0*a))/a );
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
    invW = MIN(MAX(sqrt(1.0-vsq), 1.0/eos->W_max), 1.0);
    W    = 1.0/invW;

    // Set the prims (eps is computed as in Palenzuela et al.)
    // See e.g., Eq. (44) of https://arxiv.org/pdf/1712.07538.pdf
    const double x    = z * invD;
    const double xrho = con->rho * invW;
    const double xeps = - 1.0 + (1.0-W*W)*x*invW
                      + W*( 1.0 + q - s + 0.5*( s*invW*invW + (t*t)/(x*x) ) );
    eos->tabulated_compute_P_T_from_eps(eos, xrho, xye, xeps, &xprs, &xtemp);

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

  if (step >= maxsteps)
    return roots_error_max_iter;

  if( conacc ) {     //converged on an extrap. so recompute vars
    const double a     = e + xprs + 0.5*B_squared;
    const double phi   = acos(sqrt(27.0*d/(4.0*a))/a);
    const double Eps   = a/3.0*( 1.0 - 2.0*cos( (2.0/3.0)*(phi + M_PI) ) );
    z                  = Eps-B_squared;
    const double Epssq = Eps*Eps;
    const double zsq   = z*z;
    const double vsq   = (zsq * S_squared + (z+Eps)*BdotSsq)/(zsq*Epssq);
    invW               = MIN(MAX(sqrt(1.0-vsq), 1.0/eos->W_max), 1.0);
    W                  = 1.0/invW;
    prim->rho          = con->rho*invW;
  }

  // Set the primitives
  prim->rho         = con->rho*invW;
  prim->Y_e         = xye;
  prim->temperature = xtemp;
  prim->vx          = W*(SU[0] + BdotS*BU[0]/z)/(z+B_squared);
  prim->vy          = W*(SU[1] + BdotS*BU[1]/z)/(z+B_squared);
  prim->vz          = W*(SU[2] + BdotS*BU[2]/z)/(z+B_squared);
  eos->tabulated_compute_P_eps_from_T( eos, prim->rho, prim->Y_e, prim->temperature,
                                       &prim->press, &prim->eps );

  return grhayl_success;
}

int Tabulated_Newman1D_energy(
      const GRHayL_parameters *restrict grhayl_params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      const conservative_quantities *restrict cons_undens,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics ) {

  // Step 1: Compute auxiliary quantities
  double BU[3], SU[3], Bsq, Ssq, BdotS;
  compute_BU_SU_Bsq_Ssq_BdotS(metric, cons_undens, prims,
                              BU, SU, &Bsq, &Ssq, &BdotS);

  // Step 2: Call the Newman routine that uses the energy to recover T
  const double tol_x = 1e-15;
  return newman_energy(eos, Ssq, BdotS, Bsq, BU, SU,
                       cons_undens, prims, tol_x);
}