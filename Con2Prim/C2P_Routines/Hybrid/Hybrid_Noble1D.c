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

/**********************************************************************************
  Copyright 2005 Scott C. Noble, Charles F. Gammie,
  Jonathan C. McKinney, and Luca Del Zanna


  This file is part of PVS-GRMHD.

  PVS-GRMHD is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  PVS-GRMHD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with PVS-GRMHD; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  -------------------------------------------------------------------------------
*/

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************

utoprim_1d.c:
---------------

    Uses the 1D_W method:
       -- solves for one independent variable (W) via a 1D
          Newton-Raphson method
       -- can be used (in principle) with a general equation of state.

  -- Currently returns with an error state (>0) if a negative rest-mass
      density or internal energy density is calculated.  You may want
      to change this aspect of the code so that it still calculates the
      velocity and so that you can floor the densities.  If you want to
      change this aspect of the code please comment out the "return(retval)"
      statement after "retval = 5;" statement in Utoprim_new_body();

******************************************************************************/

#define NEWT_DIM (1)

double dpdW_calc_vsq(const eos_parameters *restrict eos, double W, double vsq);

/**********************************************************************/
/******************************************************************

  Utoprim_1d():

  -- Driver for new prim. var. solver.  The driver just translates
     between the two sets of definitions for U and P.  The user may
     wish to alter the translation as they see fit.

     It assumes that on input/output:


              /  rho u^t           \
         U =  |  T^t_t   + rho u^t |  sqrt(-det(g_{\mu\nu}))
              |  T^t_\mu           |
              \   B^i              /


             /    rho        \
         P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /


     ala HARM.

   Arguments:
       U[NPR]    = conserved variables (current values on input/output);
       gcov[NDIM][NDIM] = covariant form of the metric;
       gcon[NDIM][NDIM] = contravariant form of the metric;
       gdet             = sqrt( - determinant of the metric);
       prim[NPR] = primitive variables (guess on input, calculated values on
                                        output if there are no problems);

   -- NOTE: for those using this routine for special relativistic MHD and are
            unfamiliar with metrics, merely set
              gcov = gcon = diag(-1,1,1,1)  and gdet = 1.;

******************************************************************/



/**********************************************************************/
/**********************************************************************************

  Hybrid_Noble1D():

     -- Attempt an inversion from U to prim using the initial guess prim.

     -- This is the main routine that calculates auxiliary quantities for the
        Newton-Raphson routine.

  -- assumes that
             /   rho gamma   \
         U = | alpha T^t_\mu |
             \   alpha B^i   /



             /     rho     \
      prim = |     uu      |
             | \tilde{u}^i |
             \  alpha B^i  /


return:  (i*100 + j)  where
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used)
                   or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the
                   given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence
                   (occurrence of "nan" or "+/-inf";

         j = 0 -> success
             1 -> failure: some sort of failure in Newton-Raphson;
             2 -> failure: utsq<0 w/ initial p[] guess;
             3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: v^2 > 1
             5 -> failure: rho,uu <= 0;

**********************************************************************************/

int Hybrid_Noble1D(
      const GRHayL_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      const conservative_quantities *restrict cons_undens,
      primitive_quantities *restrict prims_guess,
      con2prim_diagnostics *restrict diagnostics ) {

  double gnr_out[NEWT_DIM];

  // Contains Bsq,QdotBsq,Qsq,Qtsq,Qdotn,QdotB,D,W,W_times_S,ye
  harm_aux_vars_struct harm_aux;

  const int n = NEWT_DIM;

  // Assume ok initially:
  int retval = 0;

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  const double Bup[4] = {0.0, prims_guess->Bx * ONE_OVER_SQRT_4PI,
                              prims_guess->By * ONE_OVER_SQRT_4PI,
                              prims_guess->Bz * ONE_OVER_SQRT_4PI};

  double Bdn[4]; lower_vector(metric, Bup, Bdn);

  const double uu = - cons_undens->tau*metric->lapse
                    - (metric->lapse-1.0)*cons_undens->rho
                    + metric->betax*cons_undens->S_x
                    + metric->betay*cons_undens->S_y
                    + metric->betaz*cons_undens->S_z;

  const double Qdn[4] = {uu - cons_undens->rho,
                              cons_undens->S_x,
                              cons_undens->S_y,
                              cons_undens->S_z};

  double Qup[4]; raise_vector(metric, Qdn, Qup);

  harm_aux.Bsq = 0. ;
  for(int i=1; i<4; i++) harm_aux.Bsq += Bup[i]*Bdn[i];

  harm_aux.QdotB = 0. ;
  for(int i=0; i<4; i++) harm_aux.QdotB += Qdn[i]*Bup[i];
  harm_aux.QdotBsq = harm_aux.QdotB*harm_aux.QdotB;

  // n_{\mu}Q^{\mu} = -alpha Q^{0}, since n_{\mu} = (-alpha,0,0,0)
  harm_aux.Qdotn = -metric->lapse*Qup[0];

  harm_aux.Qsq = 0.0;
  for(int i=0; i<4; i++) harm_aux.Qsq += Qdn[i]*Qup[i] ;

  harm_aux.Qtsq = harm_aux.Qsq + harm_aux.Qdotn*harm_aux.Qdotn;

  harm_aux.D    = cons_undens->rho;

  /* calculate W from last timestep and use for guess */
  double utsq = 0.0;
  // IGM always set the velocity guesses to 0; not sure how ut^i in harm relates to v^i
  //for(int i=1; i<4; i++)
  //  for(int j=1; j<4; j++) utsq += metric->gdn[i][j]*prims[UTCON1+i-1]*prims[UTCON1+j-1];

  if( (utsq < 0.) && (fabs(utsq) < 1.0e-13) ) {
    utsq = fabs(utsq);
  }
  if(utsq < 0.0 || utsq > UTSQ_TOO_BIG) {
    retval = 2;
    return(retval);
  }

  double Wsq = 1.0 + utsq;   // Lorentz factor squared
  harm_aux.W = sqrt(Wsq); // Lorentz factor, not to be confused with Gamma

  // Always calculate rho from D and W so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq).
  double rho0 = harm_aux.D / harm_aux.W;

  // p = 0.0;
  // if( eos.is_Hybrid ) {
    const int polytropic_index = eos->hybrid_find_polytropic_index(eos, prims_guess->rho);
    const double Gamma_ppoly = eos->Gamma_ppoly[polytropic_index];
    double u = prims_guess->press/(Gamma_ppoly - 1.0);
    double p = pressure_rho0_u(eos, rho0, u);
  // } else if( eos.is_Tabulated ) {
  //   harm_aux.ye            = U[YE]/U[RHO];
  //   harm_aux.W_times_S = U[WS];
  //   harm_aux.use_entropy   = false;
  //   harm_aux.T_guess       = prim[TEMP];
  //   double xrho         = rho0;
  //   double xye          = harm_aux.ye;
  //   double xtemp        = harm_aux.T_guess;
  //   double xprs         = 0.0;
  //   double xeps         = 0.0;
  //   double xdepsdT      = 0.0;

  //   // Now compute P and eps from (rho,Ye,T). Note that
  //   // at this point we do not know W, so we do not
  //   // use the entropy in this function call.
  //   WVU_EOS_P_eps_and_depsdT_from_rho_Ye_T( xrho,xye,xtemp, &xprs,&xeps,&xdepsdT );
  //   p = xprs;
  //   u = xeps*xrho;
  //   if( xdepsdT < eos.depsdT_threshold ) harm_aux.use_entropy = true;
  // }

  double w = rho0 + u + p;
  double W_last = w*Wsq;

  // Make sure that W is large enough so that v^2 < 1 :
  int i_increase = 0;
  while( (( W_last*W_last*W_last * ( W_last + 2.*harm_aux.Bsq )
            - harm_aux.QdotBsq*(2.*W_last + harm_aux.Bsq) ) <= W_last*W_last*(harm_aux.Qtsq-harm_aux.Bsq*harm_aux.Bsq))
         && (i_increase < 10) ) {
    W_last *= 10.;
    i_increase++;
  }

  // Calculate W:
  gnr_out[0] = W_last;

  // We need a dummy variable to keep the function call in this file
  // consistent with the ones in the con2prim_Noble1D_entropy.cc and
  // con2prim_Noble1D_entropy2.cc files.
  double dummy = 0.0;
  retval = newton_raphson_1d(eos, &harm_aux, gnr_out, n, &diagnostics->n_iter, dummy, func_1d_orig);

  const double W = gnr_out[0];

  /* Problem with solver, so return denoting error before doing anything further */
  if( (retval != 0) || (W == FAIL_VAL) ) {
    retval = retval*100+1;
    return(retval);
  } else if(W <= 0. || W > W_TOO_BIG) {
    retval = 3;
    return(retval);
  }

  // Calculate v^2:
  double vsq = vsq_calc(&harm_aux, W);
//TODO: differs from Noble2D
  if( vsq >= 1. ) {
    retval = 4;
    return(retval);
  }

  // Recover the primitive variables from the scalars and conserved variables:
  const double gtmp = sqrt(1. - vsq);
  harm_aux.W = 1./gtmp;
  rho0 = harm_aux.D * gtmp;

  w = W * (1. - vsq);

  if( eos->eos_type == 0 ) {
    p = pressure_rho0_w(eos, rho0, w);
    u = w - (rho0 + p); // u = rho0 eps, w = rho0 h
  } else {
    grhayl_warn("Tabulated not implemented!");
//    double xrho  = rho0;
//    double xye   = harm_aux.ye;
//    double xtemp = harm_aux.T_guess;
//    double xent  = harm_aux.W_times_S / W;
//    double xprs  = 0.0;
//    double xuu   = 0.0;
//    double xeps  = 0.0;
//    if( harm_aux.use_entropy ) {
//      WVU_EOS_P_eps_and_T_from_rho_Ye_S( xrho,xye,xent, &xprs,&xeps,&xtemp );
//    } else {
//      xprs  = -0.5*harm_aux.Bsq/(harm_aux.W*harm_aux.W)+harm_aux.Qdotn+W+harm_aux.Bsq-0.5*harm_aux.QdotBsq/(W*W);;
//      xuu   = (W-harm_aux.D*harm_aux.W-xprs*harm_aux.W*harm_aux.W)/(harm_aux.D*harm_aux.W) * rho0;
//      xeps  = xuu/xrho;
//      WVU_EOS_P_S_and_T_from_rho_Ye_eps( xrho,xye,xeps, &xprs,&xent,&xtemp );
//    }
//
//    // Update P and T in the prim array
//    prim[RHO  ] = xrho;
//    prim[YE   ] = harm_aux.ye;
//    prim[TEMP ] = MIN(MAX(xtemp,eos.T_atm),eos.T_max);
//    prim[PRESS] = xprs;
//    prim[EPS  ] = xeps;
//    prim[ENT  ] = xent;
  }

  if( ((rho0 <= 0.) || (u <= 0.)) ) {
    // User may want to handle this case differently, e.g. do NOT return upon
    // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:
    retval = 6;
  }

const double nup[4] = {metric->lapseinv,
		      -metric->lapseinv*metric->betax,
		      -metric->lapseinv*metric->betay,
		      -metric->lapseinv*metric->betaz};

  double Qtcon[4];
  const double g_o_WBsq = harm_aux.W/(W+harm_aux.Bsq);
  const double QdB_o_W  = harm_aux.QdotB / W;

  for(int i=1; i<4; i++) Qtcon[i] = Qup[i] + nup[i] * harm_aux.Qdotn;
  double utx = g_o_WBsq * ( Qtcon[1] + QdB_o_W*Bup[1] ) ;
  double uty = g_o_WBsq * ( Qtcon[2] + QdB_o_W*Bup[2] ) ;
  double utz = g_o_WBsq * ( Qtcon[3] + QdB_o_W*Bup[3] ) ;

  prims_guess->rho = rho0;
  //Aditional tabulated code here

  double u0;
  limit_utilde_and_compute_v(eos, metric, &u0, &utx, &uty,
                                         &utz, prims_guess, diagnostics);

  if(diagnostics->vel_limited_ptcount==1)
    prims_guess->rho = cons_undens->rho/(metric->lapse*u0);

  prims_guess->press = pressure_rho0_u(eos, prims_guess->rho, u);
  prims_guess->eps = u/prims_guess->rho;
  if( params->evolve_entropy ) eos->hybrid_compute_entropy_function(eos, prims_guess->rho, prims_guess->press, &prims_guess->entropy);

  /* Done! */
  return(retval);
}
