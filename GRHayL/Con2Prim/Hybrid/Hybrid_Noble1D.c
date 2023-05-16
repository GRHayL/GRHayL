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

***********************************************************************************/

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
             2 -> failure: vsq<0 w/ initial p[] guess;
             3 -> failure: Z<0 or Z>Z_TOO_BIG
             4 -> failure: v^2 > 1
             5 -> failure: rho,uu <= 0

**********************************************************************************/

int Hybrid_Noble1D(
      const GRHayL_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      const conservative_quantities *restrict cons_undens,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics ) {

  double gnr_out[NEWT_DIM];

  // Contains Bsq,QdotBsq,Qsq,Qtsq,Qdotn,QdotB,D,W,W_times_S,ye
  harm_aux_vars_struct harm_aux;

  const int ndim = NEWT_DIM;

  // Assume ok initially:
  int retval = 0;

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  const double Bup[4] = {0.0, prims->Bx * ONE_OVER_SQRT_4PI,
                              prims->By * ONE_OVER_SQRT_4PI,
                              prims->Bz * ONE_OVER_SQRT_4PI};

  double Bdn[4]; lower_vector(metric, Bup, Bdn);

  const double uu = - cons_undens->tau*metric->lapse
                    - metric->lapse*cons_undens->rho
                    + metric->betax*cons_undens->S_x
                    + metric->betay*cons_undens->S_y
                    + metric->betaz*cons_undens->S_z;

  const double Qdn[4] = {uu,
                         cons_undens->S_x,
                         cons_undens->S_y,
                         cons_undens->S_z};

  double Qup[4]; raise_vector_4D(metric, Qdn, Qup);

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

  const double tmp_u = metric->adm_gxx * SQR(prims->vx + metric->betax) +
                                             2.0*metric->adm_gxy*(prims->vx + metric->betax)*(prims->vy + metric->betay) +
                                             2.0*metric->adm_gxz*(prims->vx + metric->betax)*(prims->vz + metric->betaz) +
                                             metric->adm_gyy * SQR(prims->vy + metric->betay) +
                                             2.0*metric->adm_gyz*(prims->vy + metric->betay)*(prims->vz + metric->betaz) +
                                             metric->adm_gzz * SQR(prims->vz + metric->betaz);

  double u0 = 1.0/sqrt(1.0-tmp_u);
  const double utilde[3] = {u0*(prims->vx + metric->betax),
                            u0*(prims->vy + metric->betay),
                            u0*(prims->vz + metric->betaz)};

  /* calculate Z from last timestep and use for guess */
  double vsq = 0.0;
  for(int i=1; i<4; i++)
    for(int j=1; j<4; j++) vsq += metric->g4dn[i][j]*-utilde[i-1]*utilde[j-1];

  if( (vsq < 0.) && (fabs(vsq) < 1.0e-13) ) {
    vsq = fabs(vsq);
  }
  if(vsq < 0.0 || vsq > UTSQ_TOO_BIG) {
    retval = 2;
    return(retval);
  }

  const double Wsq = 1.0 + vsq;   // Lorentz factor squared
  harm_aux.W = sqrt(Wsq);

  // Always calculate rho from D and W so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq).
  const double rho0 = harm_aux.D / harm_aux.W;
  double u = 0;
  double p = 0;
  double w = 0;

  if( eos->eos_type == grhayl_eos_hybrid ) {
    const double Gamma_ppoly = eos->Gamma_ppoly[eos->hybrid_find_polytropic_index(eos, prims->rho)];
    u = prims->press/(Gamma_ppoly - 1.0);
    p = pressure_rho0_u(eos, rho0, u);
    w = rho0 + u + p;
  } else if( eos->eos_type == grhayl_eos_tabulated) {
    grhayl_warn("No tabulated EOS support yet! Sorry!");
  }

  double Z_last = w*Wsq;

  // Make sure that Z is large enough so that v^2 < 1 :
  int i_increase = 0;
  while( (( Z_last*Z_last*Z_last * ( Z_last + 2.*harm_aux.Bsq )
            - harm_aux.QdotBsq*(2.*Z_last + harm_aux.Bsq) ) <= Z_last*Z_last*(harm_aux.Qtsq-harm_aux.Bsq*harm_aux.Bsq))
         && (i_increase < 10) ) {
    Z_last *= 10.;
    i_increase++;
  }

  // Calculate Z:
  gnr_out[0] = Z_last;

  // To be consistent with entropy variants, unused argument 0.0 is needed
  retval = newton_raphson_1d(eos, &harm_aux, ndim, 0.0, &diagnostics->n_iter, gnr_out, func_1d_orig);

  const double Z = gnr_out[0];

  /* Problem with solver, so return denoting error before doing anything further */
  if( (retval != 0) || (Z == FAIL_VAL) ) {
    retval = retval*100+1;
    return(retval);
  } else if(Z <= 0. || Z > Z_TOO_BIG) {
    retval = 3;
    return(retval);
  }

  // Calculate v^2:
  vsq = vsq_calc(&harm_aux, Z);
//TODO: differs from Noble2D
  if( vsq >= 1. ) {
    retval = 4;
    return(retval);
  }

  // Recover the primitive variables from the scalars and conserved variables:
  const double gtmp = sqrt(1. - vsq);
  harm_aux.W = 1.0/gtmp;
  w = Z * (1.0 - vsq);

  prims->rho = harm_aux.D * gtmp;

  if( eos->eos_type == grhayl_eos_hybrid ) {
    p = pressure_rho0_w(eos, prims->rho, w);
    u = w - (prims->rho + p); // u = rho eps, w = rho0 h
  } else if( eos->eos_type == grhayl_eos_tabulated ) {
    grhayl_warn("No tabulated EOS support yet! Sorry!");
  }

  if( ((prims->rho <= 0.0) || (u <= 0.0)) ) {
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

  limit_utilde_and_compute_v(eos, metric, &utx, &uty, &utz, prims, &diagnostics->speed_limited);

  if(diagnostics->speed_limited==1)
    prims->rho = cons_undens->rho/(metric->lapse*prims->u0);

  // Since u is dependent on rho, it seems weird to not reset u if rho has changed from above
  prims->press = pressure_rho0_u(eos, prims->rho, u);
  prims->eps = u/prims->rho;
  if( params->evolve_entropy ) eos->hybrid_compute_entropy_function(eos, prims->rho, prims->press, &prims->entropy);

  /* Done! */
  return(retval);
}
