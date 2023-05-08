#include "../harm_u2p_util.h"

/* Function    : Hybrid_Noble2D()
 * Description : Unpacks the primitive_quantities struct into the variables
                 needed by the Newton-Rapson solver provided by HARM, then
                 repacks the  primitives. This function
                 is adapted from the HARM function provided by IllinoisGRMHD. The
                 original HARM copyright is included below.

 * Inputs      : params         - GRHayL_parameters struct with parameters
 *                                for the simulation
 *             : eos            - eos_parameters struct with data for the
 *                                EOS of the simulation
 *             : metric         - metric_quantities struct with data for
 *                                the gridpoint of interest
 *             : cons           - conservative_quantities struct with data
 *                                for the gridpoint of interest
 *
 * Outputs     : prims          - returns computed primitives if Newton-Rapson
                                  method converges
 *             : diagnostics    - tracks the number of iterations for convergence
 *
 */

void NR_3D_WZT(
      const eos_parameters *restrict eos,
      const double tol_x,
      const double S_squared,
      const double BdotS,
      const double B_squared,
      const double *restrict Sup,
      const conservative_quantities *restrict cons_undens,
      primitive_quantities *restrict prims_guess,
      bool *restrict c2p_failed );

static inline double compute_W_from_cons(
      const eos_parameters *restrict eos,
      const conservative_quantities *restrict cons_undens,
      const double S_squared,
      const double B_squared,
      const double BdotS ) {

  // Use the same as the Palenzuela routine
  const double q = con[TAU]/con[DD];
  const double r = S_squared/(con[DD]*con[DD]);
  const double s = B_squared/con[DD];
  const double t = BdotS/(pow(con[DD],1.5));
  const double x = 2.0+2.0*q-s;

  double Wminus2 = 1.0 - ( x*x*r + (2*x+s)*t*t ) / ( x*x*(x+s)*(x+s) );
  Wminus2           = fmin(fmax(Wminus2,eos.inv_W_max_squared ), 1.0);
  const double W = pow(Wminus2, -0.5);
  return W;

}

/*****************************************************************************/
/************** ILLINOISGRMHD INTERFACE FOR CERDA-DURAN CON2PRIM *************/
/*****************************************************************************/
//
// This routine expects 8 conservs as input:
//
// -> D     = rho_star      / sqrt(gamma) = W * rho (W is the Lorentz factor)
// -> DYe   = Ye_Star       / sqrt(gamma)
// -> B^{i} = \tilde{B}^{i} / sqrt(gamma) / sqrt(4pi)
// -> S_{i} = \tilde{S}_{i} / sqrt(gamma)
//
// From the input quantities, we compute B_{i} and S^{i}
int Tabulated_CerdaDuran3D(
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      const conservative_quantities *restrict cons_undens,
      primitive_quantities *restrict prims_guess,
      con2prim_diagnostics *restrict diagnostics ) {

  const double Bup[3] = {prims_guess->Bx * ONE_OVER_SQRT_4PI,
                         prims_guess->By * ONE_OVER_SQRT_4PI,
                         prims_guess->Bz * ONE_OVER_SQRT_4PI};
  const double B_squared = compute_vec2_from_vcon(metric, Bup);

  const double Sdn[3] = {cons_undens->S_x,
                         cons_undens->S_y,
                         cons_undens->S_z};
  double Sup[3]; raise_vector_3d(metric, Sdn, Sup);
  const double S_squared = compute_vec2_from_vcov(metric, Sdn);


  // Enforce ceiling on S^{2} (A5 of Palenzuela et al. https://arxiv.org/pdf/1505.01607.pdf)
  double S_squared_max = SQR( con[DD] + con[TAU] );
  if( S_squared > 0.9999 * S_squared_max ) {
    // Compute rescaling factor
    double rescale_factor_must_be_less_than_one = sqrt(0.9999*S_squared_max/S_squared);
    // Rescale S_{i}
    for(int i=0;i<3;i++) SD[i] *= rescale_factor_must_be_less_than_one;
    // S_{i} has been rescaled. Recompute S^{i}.
    raise_or_lower_indices_3d( SD,gammaUU, SU );
    // Now recompute S^{2} := gamma^{ij}S^{i}S_{j}.
    S_squared = 0.0;
    for(int i=0;i<3;i++) S_squared += SU[i] * SD[i];
    // Check if the fix was successful
    if( simple_rel_err(S_squared,0.9999*S_squared_max) > 1e-12 ) CCTK_VError(VERR_DEF_PARAMS,"Incompatible values of S_squared after rescaling: %.15e %.15e\n",S_squared,0.9999*S_squared_max);
  }

  // Need to calculate for (21) and (22) in Cerda-Duran 2008
  // B * S = B^i * S_i
  double BdotS = 0.0;
  for(int i=0;i<3;i++) BdotS += Bup[i] * Sdn[i];


  bool c2p_failed = false;
  const double tol_x    = 5e-9;
  NR_3D_WZT(eos, tol_x, S_squared, BdotS, B_squared, SU, cons_undens, prims_guess, &c2p_failed);

  return c2p_failed;
}

/*****************************************************************************/
/*********************** CERDA-DURAN CON2PRIM FUNCTIONS **********************/
/*****************************************************************************/
void calc_WZT_guesses(
      const eos_parameters *restrict eos,
      const double S_squared,
      const double B_squared,
      const double BdotS,
      const conservative_quantities *restrict cons_undens,
      double *restrict xmax) {

  // First compute W following Palenzuela
  const double W = compute_W_from_cons(eos, cons_undens, S_squared, B_squared, BdotS);

  // Compute P and eps
  const double xrho  = cons_undens->rho/W;
  const double xye   = cons_undens->Y_e/cons_undens->rho;
  const double xtemp = eos.T_max; // initial guess, choose large enough
  double xprs        = 0.0;
  double xeps        = 0.0;
  WVU_EOS_P_and_eps_from_rho_Ye_T( xrho, xye, xtemp, &xprs, &xeps );

  // Compute z = rho * h * W^2
  const double h = 1.0 + xeps + xprs/xrho;
  const double z = xrho * h * W * W;

  // Now set W, z, and T
  xmax[0] = W;
  xmax[1] = z;
  xmax[2] = xtemp;

}


void calc_prim_from_x_3D_WZT(
      const eos_parameters *restrict eos,
      const double BdotS,
      const double B_squared,
      const double *restrict Sup,
      const conservative_quantities *restrict cons_undens,
      primitive_quantities *restrict prims_guess,
      double *restrict x ) {

  // Recover the primitive variables from the scalars (W,Z)
  // and conserved variables, Eq. (23)-(25) in Cerdá-Durán et al. 2008
  double W = x[0];
  double Z = x[1];
  double T = x[2];

  // Calculate press, eps etc. from (rho, temp, Ye) using EOS,
  // required for consistency
  double xrho  = cons_undens->rho/W;
  double xye   = cons_undens->Y_e/cons_undens->rho;
  double xtemp = T;
  double xeps  = 0.0;
  double xprs  = 0.0;
  eos->tabulated_compute_P_eps_from_T(eos, xrho, xye, xtemp, &xprs, &xeps);

  // Eq. (24) in Siegel et al. 2018, with S^{i} := gamma^{ij} S_{j}
  // The extra factor of W converts v^{i} to tilde(u)^{i}.
  prims_guess->vx = (Sup[0] + BdotS*prims_guess->Bx/Z)/(Z+B_squared);
  prims_guess->vy = (Sup[1] + BdotS*prims_guess->By/Z)/(Z+B_squared);
  prims_guess->vz = (Sup[2] + BdotS*prims_guess->Bz/Z)/(Z+B_squared);
  prims_guess->rho = xrho;
  prims_guess->Y_e = xye;
  prims_guess->temperature = T;
  prims_guess->press = xprs;
  prims_guess->eps = xeps;
}

void NR_step_3D_eps(
      const eos_parameters *restrict eos,
      const double S_squared,
      const double BdotS,
      const double B_squared,
      const conservative_quantities *restrict cons_undens,
      double *restrict x,
      double *restrict dx,
      double *restrict f ) {
  // Finding the roots of f(x):
  //
  // x_{n+1} = x_{n} - f(x)/J = x_{n} + dx_{n}
  //
  // where J is the Jacobian matrix J_{ij} = df_i/dx_j
  //
  // Here, compute dx = [dW, dT]
  double W = x[0];
  double Z = x[1];
  double T = x[2];

  // Need partial derivatives of specific internal energy and pressure wrt density and
  // temperature. Those need to be based on primitives computed from Newton-Raphson state
  // vector x and conservatives
  const double xrho      = con[DD]/W;
  const double xye       = con[YE]/con[DD];
  double xprs      = 0.0;
  double xeps      = 0.0;
  double xdPdrho   = 0.0;
  double xdPdT     = 0.0;
  double xdepsdrho = 0.0;
  double xdepsdT   = 0.0;
  WVU_EOS_P_eps_dPdrho_dPdT_depsdrho_and_depsdT_from_rho_Ye_T(xrho,xye,T,&xprs,&xeps,&xdPdrho,&xdPdT,&xdepsdrho,&xdepsdT);

  // Some useful auxiliary variables
  double BdotSsqr       = BdotS*BdotS;
  double z_plus_Bsq     = z + B_squared;
  double z_plus_Bsq_sqr = z_plus_Bsq*z_plus_Bsq;
  double z_sqr          = z*z;
  double W_sqr          = W*W;
  double inv_W          = 1.0/W;
  double inv_W_sqr      = 1.0/W_sqr;
  double inv_D          = 1.0/con[DD];

  // Now compute eps(W,z)
  double eps = - 1.0 + (z*inv_W - xprs*W)*inv_D;
  eps = fmax(eps,eos.eps_min);

  //---------------------------------------------------------------
  //------------ Equation (29) in Siegel et al. (2017) ------------
  //---------------------------------------------------------------
  // f0 = ( (z+B^{2})^{2} - S^{2} - (2z+B^{2})(B.S)^{2}/z^{2} )W^{2} - (z+B^{2})^{2}
  double f0 = (z_plus_Bsq_sqr - S_squared - (z+z_plus_Bsq)*BdotSsqr/z_sqr)*W_sqr - z_plus_Bsq_sqr;

  // df0/dW = 2W( (z+B^{2})^{2} - S^{2} - (2z+B^{2})(B.S)^{2}/z^{2} )
  double a  = 2.0*(z_plus_Bsq_sqr - S_squared - (z+z_plus_Bsq)*(BdotSsqr)/z_sqr)*W;

  // df0/dz = ( 2(z+B^{2}) + 2(B.S)^{2}/z^{2} + 2(B^{2})(B.S)^{2}/z^{3} )W^{2} - 2(z+B^{2})
  double b  = (2.0*z_plus_Bsq + (2.0/(z_sqr) + 2.0*B_squared/(z_sqr*z))*BdotSsqr)*W_sqr - 2.0*z_plus_Bsq;

  // df0/dT = 0
  double c  = 0.0;
  //---------------------------------------------------------------

  //---------------------------------------------------------------
  //------------ Equation (30) in Siegel et al. (2017) ------------
  //---------------------------------------------------------------
  // f1 = (tau + D - z - B^{2} + (B.S)^{2}/(2z^{2}) + P)W^{2} + B^{2}/2
  double f1 = (con[TAU] + con[DD] - z - B_squared + (BdotSsqr)/(2.0*z_sqr) + xprs)*W_sqr + 0.5*B_squared;

  // Recall that D = W rho => rho = D / W. Thus
  //
  // dP/dW = (dP/drho)(drho/dW) = -W^{2}D(dP/drho)
  //
  // df1/dW = 2(tau + D - z - B^{2} + (B.S)^{2}/(2z^{2}) + P)W + dPdW * W^{2}
  double d  = 2.0*(con[TAU] + con[DD] - z - B_squared + BdotSsqr/(2.0*z_sqr) + xprs)*W - con[DD]*xdPdrho;

  // df1/dz = ( -1 - (B.S)^{2}/z^{3} ) W^{2}
  double e  = (-1.0-BdotSsqr/(z_sqr*z))*W_sqr;

  // df1/dT = W^{2}dPdT
  double fc = W_sqr * xdPdT;

  //---------------------------------------------------------------

  //---------------------------------------------------------------
  //------------ Equation (31) in Siegel et al. (2017) ------------
  //---------------------------------------------------------------
  // Note that we use the *pressure* instead of eps
  // f2 = P(W,z) - P(rho,Ye,T)
  double f2 = eps - xeps;

  // df2/dW = deps(W,z)/dW - (deps(rho,Ye,T)/drho)*(drho/dW)
  //        = -z/(D W^{2}) - P/W + (1/W)(dP/drho) + (D/W^{2})(deps/drho)
  double g  = (con[DD]*xdepsdrho - z*inv_D)*inv_W_sqr + (xdPdrho - xprs)*inv_W;

  // df2/dz = deps(W,z)/dz = 1/(D W)
  double h  = inv_W*inv_D;

  // df2/dT = -(W/D)(dP/dT) - deps(rho,Ye,T)/dT
  double k  = -W*inv_D*xdPdT - xdepsdT;

  // Set the vector of equations
  f[0] = f0;
  f[1] = f1;
  f[2] = f2;

  // Compute the determinant
  double A    = e*k-fc*h;
  double B    = fc*g-d*k;
  double C    = d*h-e*g;
  double detJ = a*(A) + b*(B) + c*(C);

  //Compute the matrix inverse
  double Ji[3][3];
  Ji[0][0] = A/detJ;
  Ji[1][0] = B/detJ;
  Ji[2][0] = C/detJ;
  Ji[0][1] = (c*h-b*k)/detJ;
  Ji[1][1] = (a*k-c*g)/detJ;
  Ji[2][1] = (g*b-a*h)/detJ;
  Ji[0][2] = (b*fc-c*e)/detJ;
  Ji[1][2] = (c*d-a*fc)/detJ;
  Ji[2][2] = (a*e-b*d)/detJ;

  // Compute the step size
  dx[0] = -Ji[0][0]* f[0] - Ji[0][1]* f[1] - Ji[0][2]* f[2];
  dx[1] = -Ji[1][0]* f[0] - Ji[1][1]* f[1] - Ji[1][2]* f[2];
  dx[2] = -Ji[2][0]* f[0] - Ji[2][1]* f[1] - Ji[2][2]* f[2];

}

void NR_3D_WZT(
      const eos_parameters *restrict eos,
      const double tol_x,
      const double S_squared,
      const double BdotS,
      const double B_squared,
      const double *restrict Sup,
      const conservative_quantities *restrict cons_undens,
      primitive_quantities *restrict prims_guess,
      bool *restrict c2p_failed ) {

  // 2D Newton-Raphson scheme, using state vector x = (W, T) and 2D function
  // f(x) = (f1(x), f2(x)) given by Eqs. (27), (28) of Siegel et al. 2018

  // initialize Newton-Raphson state vector x  = (W, T)
  double x[3];
  for(int i=0; i<3; i++) {
    x[i] = 0.0;
  }

  int count = 0;      // count number of iterations
  double f[3];        // root finding function f
  double x_lowlim[3]; // lower limits on W, T
  double dx[3];       // displacement vector
  double x_old[3];    // old state vector
  double error[3];    // error vector abs(dx/x)

  x_lowlim[0] = 1.0;
  x_lowlim[1] = eos->rho_min;
  x_lowlim[2] = eos->T_atm; // exp( ( log(eos.T_max)+log(eos.T_min) )/2.0 );

  // set initial guess for Newton-Raphson state vector x = (W, T)
  calc_WZT_guesses(eos, S_squared, B_squared, BdotS, cons_undens, x);


  // initialize variables
  for(int i=0; i<3; i++) {
    x_old[i] = 0.0;
    dx[i]    = 0.0;
    f[i]     = 0.0;
    error[i] = 0.0;
    //check for NaNs
    if( x[i] != x[i] ) {
      *c2p_failed = true;
      return;
    }
  }

  int i_extra     = 0;
  int doing_extra = 0;

  if (EXTRA_NEWT_ITER==0) {
    i_extra=-1;
  }

  bool keep_iterating=true;
  double maxerror;

  while(keep_iterating) {

    // do Newton-Raphson step
    NR_step_3D_eps(eos, S_squared, BdotS, B_squared, cons_undens, x, dx, f);

    // Update x vector and compute error
    for(int i=0; i<3; i++) {

      x_old[i] = x[i];
      // update x and exclude unphysical regime
      x[i]     = fmax(x[i]+ dx[i], x_lowlim[i]);
      error[i] = fabs((x[i]-x_old[i])/x[i]);

      // Check for NaNs
      if( x[i] != x[i] ) {
        *c2p_failed = true;
        return;
      }

    }
    maxerror = MAX(error[0], error[1]);
    count++;

    // termination criterion
    if( (fabs(maxerror) <= tol_x) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++;

    if( ((fabs(maxerror) <= tol_x)&&(doing_extra == 0))
        || (i_extra >= EXTRA_NEWT_ITER) || (count >= (MAX_NEWT_ITER)) ) {
      keep_iterating = false;
    }

  } // END of while(keep_iterating)

  //  Check for bad untrapped divergences
  //if( (!robust_isfinite(f[0])) || (!robust_isfinite(f[1])) ) {
  //  *c2p_failed = true;
  //}

  if( fabs(maxerror) <= tol_x ){
    *c2p_failed = false;
  }
  else if( (fabs(maxerror) <= tol_x) && (fabs(maxerror) > tol_x) ){
    *c2p_failed = false;
  }
  else {
    *c2p_failed = true;
  }

  // Recover the primitive variables from the final scalars x = (W, T)
  calc_prim_from_x_3D_WZT(eos, BdotS, B_squared, Sup, cons_undens, prims_guess, x);
}
