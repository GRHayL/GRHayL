#include "../harm_u2p_util.h"

/* Function    : Hybrid_Noble2D()
 * Description : Unpacks the ghl_primitive_quantities struct into the variables
                 needed by the Newton-Rapson solver provided by HARM, then
                 repacks the  primitives. This function
                 is adapted from the HARM function provided by IllinoisGRMHD. The
                 original HARM copyright is included below.

 * Inputs      : params         - ghl_parameters struct with parameters
 *                                for the simulation
 *             : eos            - ghl_eos_parameters struct with data for the
 *                                EOS of the simulation
 *             : metric         - ghl_metric_quantities struct with data for
 *                                the gridpoint of interest
 *             : cons           - ghl_conservative_quantities struct with data
 *                                for the gridpoint of interest
 *
 * Outputs     : prims          - returns computed primitives if Newton-Rapson
                                  method converges
 *             : diagnostics    - tracks the number of iterations for convergence
 *
 */

void NR_2D_WT(
      const ghl_eos_parameters *restrict eos,
      const int safe_guess,
      const double tol_x,
      const double S_squared,
      const double BdotS,
      const double B_squared,
      const double *restrict SU,
      const ghl_conservative_quantities *restrict cons_undens,
      double W,
      ghl_primitive_quantities *restrict prims_guess,
      bool *restrict c2p_failed );

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
int Tabulated_CerdaDuran2D(
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims_guess,
      ghl_con2prim_diagnostics *restrict diagnostics ) {

  const double BU[3] = {prims_guess->BU[0] * ONE_OVER_SQRT_4PI,
                        prims_guess->BU[1] * ONE_OVER_SQRT_4PI,
                        prims_guess->BU[2] * ONE_OVER_SQRT_4PI};
  const double B_squared = ghl_compute_vec2_from_vec3D(metric->gammaDD, BU);

  const double SD[3] = {cons_undens->SD[0],
                        cons_undens->SD[1],
                        cons_undens->SD[2]};
  double SU[3]; ghl_raise_lower_vector_3D(metric->gammaUU, SD, SU);
  const double S_squared = ghl_compute_vec2_from_vec3D(metric->gammaUU, SD);

  // Need to calculate for (21) and (22) in Cerda-Duran 2008
  // B * S = B^i * S_i
  double BdotS = 0.0;
  for(int i=0;i<3;i++) BdotS += BU[i] * SD[i];

  double W = 0;

  bool c2p_failed = false;
  int safe_guess  = 0;
  const double tol_x    = 5e-9;
  NR_2D_WT(eos, safe_guess, tol_x, S_squared, BdotS, B_squared, SU, cons_undens, W, prims_guess, &c2p_failed);

  if( c2p_failed ) {
    // If failed to recover the prims, try again with safe guesses
    int safe_guess=1;
    NR_2D_WT(eos, safe_guess, tol_x, S_squared, BdotS, B_squared, SU, cons_undens, W, prims_guess, &c2p_failed);
  }

  return c2p_failed;
}

/*****************************************************************************/
/*********************** CERDA-DURAN CON2PRIM FUNCTIONS **********************/
/*****************************************************************************/
void calc_WT_max(
      const ghl_eos_parameters *restrict eos,
      const double B_squared,
      const ghl_conservative_quantities *restrict cons_undens,
      double *restrict xmax) {

  // Calculate maximum values for x = (rho, T) ("safe guess" initial values)
  // cf. Cerda-Duran et al. 2008, Eq. (39)-(42)

  double rhomax = cons_undens->rho;
  double epsmax = (cons_undens->tau - B_squared/2.0) / cons_undens->rho;

  // ensure that rhomax and epsmax are in validity range of EOS
  // Note that in setting the IllinoisGRMHD EOS parameters, we
  // already impose a safety factor on the table floors and
  // ceilings, so we don't need to use the prefactors of 95%
  if(rhomax > eos->rho_max) rhomax = eos->rho_max;
  if(epsmax > eos->eps_max) epsmax = eos->eps_max;

  // Now compute P max and T max
  double xye   = cons_undens->Y_e/cons_undens->rho;
  double xtemp = eos->T_max; // initial guess, choose large enough
  double xprs  = 0.0;
  ghl_tabulated_compute_P_T_from_eps(eos, rhomax, xye, epsmax, &xprs, &xtemp);

  // Now set W_max and T_max
  xmax[0] = 1.0e4;
  xmax[1] = xtemp;
}


void calc_prim_from_x_2D_WT(
      const ghl_eos_parameters *restrict eos,
      const double BdotS,
      const double B_squared,
      const double *restrict SU,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims_guess,
      double *restrict x ) {

  // Recover the primitive variables from the scalars (W,Z)
  // and conserved variables, Eq. (23)-(25) in Cerdá-Durán et al. 2008
  double W = x[0];
  double T = x[1];

  // Calculate press, eps etc. from (rho, temp, Ye) using EOS,
  // required for consistency
  double xrho  = cons_undens->rho/W;
  double xye   = cons_undens->Y_e/cons_undens->rho;
  double xtemp = T;
  double xeps  = 0.0;
  double xprs  = 0.0;
  ghl_tabulated_compute_P_eps_from_T(eos, xrho, xye, xtemp, &xprs, &xeps);

  const double Z = cons_undens->rho * (1.0 + xeps + xprs/xrho) * W;

  // Eq. (24) in Siegel et al. 2018, with S^{i} := gamma^{ij} S_{j}
  // The extra factor of W converts v^{i} to tilde(u)^{i}.
  prims_guess->vU[0] = (SU[0] + BdotS*prims_guess->BU[0]/Z)/(Z+B_squared);
  prims_guess->vU[1] = (SU[1] + BdotS*prims_guess->BU[1]/Z)/(Z+B_squared);
  prims_guess->vU[2] = (SU[2] + BdotS*prims_guess->BU[2]/Z)/(Z+B_squared);
  prims_guess->rho = xrho;
  prims_guess->Y_e = xye;
  prims_guess->temperature = T;
  prims_guess->press = xprs;
  prims_guess->eps = xeps;
}

void NR_step_2D_WT(
      const ghl_eos_parameters *restrict eos,
      const double S_squared,
      const double BdotS,
      const double B_squared,
      const ghl_conservative_quantities *restrict cons_undens,
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
  double T = x[1];

  // Need partial derivatives of specific internal energy and pressure wrt density and
  // temperature. Those need to be based on primitives computed from Newton-Raphson state
  // vector x and conservatives
  const double rho      = cons_undens->rho/W;
  const double ye       = cons_undens->Y_e/cons_undens->rho;
  double P        = 0.0;
  double eps      = 0.0;
  double dPdrho   = 0.0;
  double dPdT     = 0.0;
  double depsdrho = 0.0;
  double depsdT   = 0.0;
  //ghl_tabulated_compute_P_eps_dPdrho_dPdT_depsdrho_depsdT_from_T(eos, rho, ye, T, &P, &eps, &dPdrho, &dPdT, &depsdrho, &depsdT);

  // h = 1 + eps + P/rho
  const double h     = 1.0 + eps + P / rho;

  // z = rho*h*W^{2} = ( rho + eps*rho + P )*W^{2}
  const double W_sqr = W*W;
  const double z     = rho * h * W_sqr;

  // Some useful auxiliary variables
  const double BdotSsqr       = BdotS*BdotS;
  const double z_plus_Bsq     = z + B_squared;
  const double z_plus_Bsq_sqr = z_plus_Bsq*z_plus_Bsq;
  const double z_sqr          = z*z;

  // dz/dW = 2 * rho * h * W = 2 * D * h
  const double dz_dW = cons_undens->rho*h;

  // dz/dT = ( rho*deps/dT + dP/dT )*W^{2}
  const double dz_dT = (rho * depsdT + dPdT)*W_sqr;

  // f1 = ( (z+B^{2})^{2} - S^{2} - (2z+B^{2})(B.S)^{2}/z^{2} )W^{2} - (z+B^{2})^{2}
  f[0] = (z_plus_Bsq_sqr - S_squared - (z+z_plus_Bsq)*BdotSsqr/z_sqr)*W_sqr - z_plus_Bsq_sqr;

  // df1/dz = ( 2(z+B^{2}) + 2(B.S)^{2}/z^{2} + 2(B^{2})(B.S)^{2}/z^{3} )W^{2} - 2(z+B^{2})
  const double df1_dz = (2.0*z_plus_Bsq + (2.0/(z_sqr) + 2.0*B_squared/(z_sqr*z))*BdotSsqr)*W_sqr - 2.0*z_plus_Bsq;

  // Ignoring the fact that z = z(W) for now, we have:
  // df1/dW = 2W( (z+B^{2})^{2} - S^{2} - (2z+B^{2})(B.S)^{2}/z^{2} )
  const double df1_dW = 2.0*(z_plus_Bsq_sqr - S_squared - (z+z_plus_Bsq)*(BdotSsqr)/z_sqr)*W;

  // f2 = (tau + D - z - B^{2} + (B.S)^{2}/(2z^{2}) + P)W^{2} + B^{2}/2
  f[1] = (cons_undens->tau + cons_undens->rho - z - B_squared + (BdotSsqr)/(2.0*z_sqr) + P)*W_sqr + 0.5*B_squared;

  // df2/dz = ( -1 - (B.S)^{2}/z^{3} ) W^{2}
  const double df2_dz = (-1.0-BdotSsqr/(z_sqr*z))*W_sqr;

  // df2/dW = 2(tau + D - z - B^{2} + (B.S)^{2}/(2z^{2}) + P)W
  const double df2_dW = 2.0*(cons_undens->tau + cons_undens->rho - z - B_squared + BdotSsqr/(2.0*z_sqr) + P)*W;

  // df2/dP = W^{2}
  const double df2_dP = W_sqr;

  // Now compute the Jacobian matrix
  double J[2][2];
  // df1_dW
  J[0][0] = df1_dz * dz_dW + df1_dW;
  const double a = J[0][0];

  // df1_dT
  J[0][1] = df1_dz * dz_dT;
  const double b = J[0][1];

  // df2_dW
  J[1][0] = df2_dz * dz_dW + df2_dW;
  const double c = J[1][0];

  // df2_dT
  J[1][1] = df2_dz * dz_dT + df2_dP * dPdT;
  const double d = J[1][1];

  // Then the inverse Jacobian Matrix
  double Ji[2][2];
  const double detJ = a*d - b*c;
  Ji[0][0]    = +d/detJ;
  Ji[0][1]    = -b/detJ;
  Ji[1][0]    = -c/detJ;
  Ji[1][1]    = +a/detJ;

  // Compute the step size
  dx[0] = -Ji[0][0]* f[0] - Ji[0][1]* f[1];
  dx[1] = -Ji[1][0]* f[0] - Ji[1][1]* f[1];

}


void NR_2D_WT(
      const ghl_eos_parameters *restrict eos,
      const int safe_guess,
      const double tol_x,
      const double S_squared,
      const double BdotS,
      const double B_squared,
      const double *restrict SU,
      const ghl_conservative_quantities *restrict cons_undens,
      double W,
      ghl_primitive_quantities *restrict prims_guess,
      bool *restrict c2p_failed ) {

  // 2D Newton-Raphson scheme, using state vector x = (W, T) and 2D function
  // f(x) = (f1(x), f2(x)) given by Eqs. (27), (28) of Siegel et al. 2018

  // initialize Newton-Raphson state vector x  = (W, T)
  double x[2];
  for(int i=0; i<2; i++) {
    x[i] = 0.0;
  }

  int count = 0;      // count number of iterations
  double f[2];        // root finding function f
  double x_lowlim[2]; // lower limits on W, T
  double dx[2];       // displacement vector
  double x_old[2];    // old state vector
  double error[2];    // error vector abs(dx/x)

  x_lowlim[0] = 1.0;
  x_lowlim[1] = eos->T_min;

  // set initial guess for Newton-Raphson state vector x = (W, T)
  if( safe_guess==0 ) {
    x[0] = W;
    x[1] = prims_guess->temperature;
  }
  else if( safe_guess==1 ) {
    calc_WT_max(eos, B_squared, cons_undens, x);
  }

  // initialize variables
  for(int i=0; i<2; i++) {
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
    NR_step_2D_WT(eos, S_squared, BdotS, B_squared, cons_undens, x, dx, f);

    // Update x vector and compute error
    for(int i=0; i<2; i++) {

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
  calc_prim_from_x_2D_WT(eos, BdotS, B_squared, SU, cons_undens, prims_guess, x);
}
