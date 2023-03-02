#include "../utils.h"

int newman(
      const eos_parameters *restrict eos,
      const double S_squared,
      const double BdotS,
      const double B_squared,
      const double *restrict BU,
      const double *restrict SU,
      const conservative_quantities *restrict con,
      primitive_quantities *restrict prim,
      const double tol_x ) {

  bool conacc = false;

  // We always guess W = 1, since we don't keep
  // track of the velocities in between time
  // steps and therefore do not have a good
  // guess for them
  double invD   = 1.0/con->rho;
  double xrho   = con->rho;
  double xye    = con->Y_e * invD;
  double xtemp  = eos->T_atm;
  double xprs   = 0.0;
  double xeps   = 0.0;

  // Compute P and eps
  // eos->tabulated_compute_P_eps_from_T(eos, xrho, xye, xtemp, &xprs, &xeps);

  //Now, begin iterative procedure to derive the primitive variables
  double P_old = xprs; // -> setup pressure initial guess
  unsigned step=0;

  const unsigned maxsteps = 300;
  double AtP[maxsteps]; // not length 3 in case of extrap. probs
  double AtR;
  int AtStep=0;
  AtP[0]=xprs;

  // Leo mod: compute auxiliary variables so we can find eps
  const double q = con->tau * invD;
  const double s = B_squared * invD;
  const double t = BdotS/(pow(con->rho,1.5));

  // d = 0.5( S^{2}*B^{2} - (B.S)^{2} ) (eq. 5.7 in Newman & Hamlin 2014)
  const double d = MAX(0.5*(S_squared*B_squared-BdotS*BdotS),0.0);
  // e = tau + D
  const double e = con->tau + con->rho;
  // z = rho*h*W^{2} = D*h*W; initialize to zero.
  double z=0, W=0;

  // Useful auxiliary variable
  const double BdotSsq = BdotS * BdotS;

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
    const double zBsq   = z + B_squared;
    const double zBsqsq = zBsq*zBsq;
    const double zsq    = z*z;
    const double vsq    = (zsq * S_squared + (z+zBsq)*BdotSsq)/(zsq*zBsqsq);

    // Then compute W^{2} = 1/(1-v^{2})
    const double Wsq    = 1.0/(1.0-vsq);

    // Impose physical limits on W
    W = MIN(MAX(sqrt(Wsq), 1.0), eos->W_max);

    // Then compute rho = D/W
    double xrho = con->rho/W;

    // Initialize P, eps, and S to zero
    xprs = 0.0;
    xeps = 0.0;

    // Remember that
    //
    // z = rho * h * W^{2} = D * h * W
    //
    // and therefore:
    //
    // x = h * W = z / D
    const double x = z * invD;

    // Now use the Palenzuela formula:
    //
    // eps = - 1.0 + x(1-W^{2})/W + W( 1 + q - s + 0.5*( s/W^{2} + t^{2}/x^{2} ) )
    xeps = - 1.0 + (1.0-W*W)*x/W + W*( 1.0 + q - s + 0.5*( s/(W*W) + (t*t)/(x*x) ) );
    // Then compute P, S, and T using (rho,Ye,eps)
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
    const double a      = e + xprs + 0.5*B_squared;
    const double phi    = acos(sqrt(27.0*d/(4.0*a))/a);
    const double Eps    = a/3.0*( 1.0 - 2.0*cos( (2.0/3.0)*(phi + M_PI) ) );
    z = Eps-B_squared;
    const double zBsq   = z + B_squared;
    const double zBsqsq = zBsq*zBsq;
    const double zsq    = z*z;
    const double vsq    = (zsq * S_squared + (z+zBsq)*BdotSsq)/(zsq*zBsqsq);
    const double Wsq    = 1.0/(1.0-vsq);
    W                   = MIN(MAX(sqrt(Wsq), 1.0), eos->W_max);
    prim->rho           = con->rho/W; //    rho[s] = tildeD[s]/(sqrtDetg[s]*W[s]);
  }

  prim->rho         = con->rho/W;
  prim->Y_e         = xye;
  prim->temperature = xtemp;
  prim->vx          = W*(SU[0] + BdotS*BU[0]/z)/(z+B_squared);
  prim->vy          = W*(SU[1] + BdotS*BU[1]/z)/(z+B_squared);
  prim->vz          = W*(SU[2] + BdotS*BU[2]/z)/(z+B_squared);

  // If using the specific internal energy, remember that
  //
  // z = rho * h * W^{2} = D * h * W
  //
  // and therefore:
  //
  // x = h * W = z / D
  const double x = z * invD;

  // Now use the Palenzuela formula:
  //
  // eps = - 1.0 + x(1-W^{2})/W + W( 1 + q - s + 0.5*( s/W^{2} + t^{2}/x^{2} ) )
  // prim->eps = - 1.0 + (1.0-W*W)*x/W + W*( 1.0 + q - s + 0.5*( s/(W*W) + (t*t)/(x*x) ) );
  // Then compute P, S, and T using (rho,Ye,eps)
  // eos->tabulated_compute_P_T_from_eps( eos, prim->rho, prim->Y_e, prim->eps,
                                       // &prim->press, &prim->temperature );
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

  // Step 1: Compute S^{2} = gamma^{ij}S_{i}S_{j}
  double SD[3] = {cons_undens->S_x, cons_undens->S_y, cons_undens->S_z};
  double S_squared = compute_S_squared(metric, SD);

  // Step 2: Enforce ceiling on S^{2} (Eq. A5 of [1])
  // Step 2.1: Compute maximum allowed value for S^{2}
  const double S_squared_max = SQR(cons_undens->tau + cons_undens->rho);
  if( S_squared > S_squared_max ) {
    // Step 2.2: Rescale S_{i}
    const double rescale_factor = sqrt(0.9999*S_squared_max/S_squared);
    for(int i=0;i<3;i++)
      SD[i] *= rescale_factor;

    // Step 2.3: Recompute S^{2}
    S_squared = compute_S_squared(metric, SD);
  }

  // Step 3: Compute B^{2} = gamma_{ij}B^{i}B^{j}
  const double BU[3] = {prims->Bx * ONE_OVER_SQRT_4PI,
                        prims->By * ONE_OVER_SQRT_4PI,
                        prims->Bz * ONE_OVER_SQRT_4PI};
  const double B_squared = compute_Bsq_from_Bup(metric, BU);

  // Step 4: Compute B.S = B^{i}S_{i}
  double BdotS = 0.0;
  for(int i=0;i<3;i++) BdotS += BU[i]*SD[i];

  // Step 5: Compute S^{i}
  double SU[3];
  raise_vector_3d(metric, SD, SU);

  const double tol_x = 1e-10;
  return newman(eos, S_squared, BdotS, B_squared, BU, SU,
                cons_undens, prims, tol_x);
}
