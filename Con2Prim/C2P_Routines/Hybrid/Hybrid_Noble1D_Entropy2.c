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

/***********************************************************************************
  Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble,
                 Gabor Toth, and Luca Del Zanna

                      HARM  version 1.0   (released May 1, 2006)

  This file is part of HARM.  HARM is a program that solves hyperbolic
  partial differential equations in conservative form using high-resolution
  shock-capturing techniques.  This version of HARM has been configured to
  solve the relativistic magnetohydrodynamic equations of motion on a
  stationary black hole spacetime in Kerr-Schild coordinates to evolve
  an accretion disk model.

  You are morally obligated to cite the following two papers in his/her
  scientific literature that results from use of any part of HARM:

  [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003,
      Astrophysical Journal, 589, 444.

  [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006,
      Astrophysical Journal, 641, 626.


  Further, we strongly encourage you to obtain the latest version of
  HARM directly from our distribution website:
  http://rainman.astro.uiuc.edu/codelib/


  HARM is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  HARM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with HARM; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************

utoprim_1d_ee2.c:
---------------

  -- uses eq. (27) of Noble  et al. or the "momentum equation" and ignores
        the energy equation (29) in order to use the additional EOS, which
        is

             P = Sc rho^(GAMMA-1) / gamma

    Uses a method similiar to  1D_W method:
       -- solves for one independent variable (rho) via a 1D
          Newton-Raphson method
       -- by substituting
          W = Dc ( Dc + GAMMA Sc rho^(GAMMA-1) / (GAMMA-1) ) / rho
           into Qtsq equation, one can get one equation for
           one unknown (rho)

       -- can be used (in principle) with a general equation of state.

  -- Currently returns with an error state (>0) if a negative rest-mass
      density or internal energy density is calculated.  You may want
      to change this aspect of the code so that it still calculates the
      velocity and so that you can floor the densities.  If you want to
      change this aspect of the code please comment out the "return(retval)"
      statement after "retval = 5;" statement in Utoprim_new_body();

******************************************************************************/

#define NEWT_DIM (1)

/**********************************************************************************

  Hybrid_Noble1D_Entropy2():

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
             2 -> failure: vsq<0 w/ initial p[] guess;
             3 -> failure: rho<0 or rho > eos->rho_max
             4 -> failure: v^2 > 1
             5 -> failure: rho,uu <= 0

**********************************************************************************/

int Hybrid_Noble1D_Entropy2(
      const GRHayL_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      const conservative_quantities *restrict cons_undens,
      primitive_quantities *restrict prims_guess,
      con2prim_diagnostics *restrict diagnostics ) {

  double gnr_out[NEWT_DIM];

  // Contains Bsq,QdotBsq,Qsq,Qtsq,Qdotn,QdotB,D,W,W_times_S,ye
  harm_aux_vars_struct harm_aux;

  const int ndim = NEWT_DIM;

  // Assume ok initially:
  int retval = 0;

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  const double Bup[4] = {0.0, prims_guess->Bx * ONE_OVER_SQRT_4PI,
                              prims_guess->By * ONE_OVER_SQRT_4PI,
                              prims_guess->Bz * ONE_OVER_SQRT_4PI};

  double Bdn[4]; lower_vector(metric, Bup, Bdn);

  // W_times_S
  harm_aux.W_times_S = cons_undens->entropy;

  const double uu = - cons_undens->tau*metric->lapse
                    - metric->lapse*cons_undens->rho
                    + metric->betax*cons_undens->S_x
                    + metric->betay*cons_undens->S_y
                    + metric->betaz*cons_undens->S_z;

  const double Qdn[4] = {uu,
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
  double vsq = 0.0;
  // IGM always set the velocity guesses to 0; not sure how ut^i in harm relates to v^i
  //for(int i=1; i<4; i++)
  //  for(int j=1; j<4; j++) vsq += metric->gdn[i][j]*prims[UTCON1+i-1]*prims[UTCON1+j-1];

  if( (vsq < 0.) && (fabs(vsq) < 1.0e-13) ) {
    vsq = fabs(vsq);
  }
  if(vsq < 0.0 || vsq > UTSQ_TOO_BIG) {
    retval = 2;
    return(retval);
  }

  double Wsq = 1.0 + vsq;   // Lorentz factor squared
  harm_aux.W = sqrt(Wsq);

  // Always calculate rho from D and W so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq).
  double rho0 = harm_aux.D / harm_aux.W;
  double u = 0;
  double p = 0;
  double w = 0;

  if( eos->eos_type == 0 ) {
    const double Gamma_ppoly = eos->Gamma_ppoly[eos->hybrid_find_polytropic_index(eos, rho0)];
    const double Gm1        = Gamma_ppoly - 1.0;                            // HARM auxiliary variable
    const double rho_Gm1    = pow(rho0,Gm1);                               // HARM auxiliary variable

    /* The definition of the entropy density, S, is
     *
     * S = P / rho^(Gamma - 1)
     *
     * Thus we have
     * .-------------------------.
     * | P = rho^(Gamma - 1) * S |
     * .-------------------------.
     */
    p = cons_undens->entropy * rho_Gm1;
    u = p/Gm1;
    w = rho0 + u + p;
  } else if( eos->eos_type == 1 ) {
    grhayl_warn("No tabulated EOS support yet! Sorry!");
  }

  gnr_out[0] = rho0;
  // To be consistent with entropy variants, unused argument 0.0 is needed
  retval = newton_raphson_1d(eos, &harm_aux, ndim, 0.0, &diagnostics->n_iter, gnr_out, func_rho2);
  rho0 = gnr_out[0];

  /* Problem with solver, so return denoting error before doing anything further */
  // HARM uses error checks for Z, which doesn't make sense for this function since
  // it returns rho; thus, we have changed the errors
  if( (retval != 0) ) {
    retval = retval*100+1;
    return(retval);
  } else if(rho0 < 0 || rho0 > eos->rho_max) {
    retval = 3;
    return(retval);
  }

  // Calculate v^2 :
  if( eos->eos_type == 0 ) {
    const double Gamma_ppoly = eos->Gamma_ppoly[eos->hybrid_find_polytropic_index(eos, rho0)];
    const double Gm1        = Gamma_ppoly - 1.0;
    const double rho_Gm1    = pow(rho0,Gm1);
    p = cons_undens->entropy * rho_Gm1;
    u = p/Gm1;
    w = rho0 + u + p;
  } else if( eos->eos_type == 1 ) {
    grhayl_warn("No tabulated EOS support yet! Sorry!");
  }

  double rel_err = (harm_aux.D != 0.0) ? fabs((harm_aux.D-rho0)/harm_aux.D) : ( (rho0 != 0.0) ? fabs((harm_aux.D-rho0)/rho0) : 0.0 );
  vsq = ( rel_err > 1e-15 ) ? (harm_aux.D-rho0)*(harm_aux.D+rho0)/(rho0*rho0) : 0.0;

  if( vsq < 0. ) {
    retval = 4;
    return(retval);
  }

  // Recover the primitive variables from the scalars and conserved variables:
  Wsq        = 1.0+vsq;
  harm_aux.W = sqrt(Wsq);
  const double Z = w * Wsq;

  prims_guess->rho = rho0;


  if( ((prims_guess->rho <= 0.0) || (u <= 0.0)) ) {
    // User may want to handle this case differently, e.g. do NOT return upon
    // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:
    retval = 5;
  }

  const double nup[4] = {metric->lapseinv,
                        -metric->lapseinv*metric->betax,
                        -metric->lapseinv*metric->betay,
                        -metric->lapseinv*metric->betaz};

  double Qtcon[4];
  const double g_o_ZBsq = harm_aux.W/(Z+harm_aux.Bsq);
  const double QdB_o_Z  = harm_aux.QdotB / Z;

  for(int i=1; i<4; i++) Qtcon[i] = Qup[i] + nup[i] * harm_aux.Qdotn;
  double utx = g_o_ZBsq * ( Qtcon[1] + QdB_o_Z*Bup[1] ) ;
  double uty = g_o_ZBsq * ( Qtcon[2] + QdB_o_Z*Bup[2] ) ;
  double utz = g_o_ZBsq * ( Qtcon[3] + QdB_o_Z*Bup[3] ) ;

  //Additional tabulated code here

  double u0;
  limit_utilde_and_compute_v(eos, metric, &u0, &utx, &uty,
                                         &utz, prims_guess, diagnostics);

  if(diagnostics->vel_limited_ptcount==1)
    prims_guess->rho = cons_undens->rho/(metric->lapse*u0);

  // Since u is dependent on rho, it seems weird to not reset u if rho has changed from above
  prims_guess->press = pressure_rho0_u(eos, prims_guess->rho, u);
  prims_guess->eps = u/prims_guess->rho;
  if( params->evolve_entropy ) eos->hybrid_compute_entropy_function(eos, prims_guess->rho, prims_guess->press, &prims_guess->entropy);

  /* Done! */
  return(retval);
}
