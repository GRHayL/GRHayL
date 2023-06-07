#include "../harm_u2p_util.h"

/* Function    : Hybrid_Noble2D()
 * Description : Unpacks the primitive_quantities struct into the variables
                 needed by the Newton-Rapson solver provided by HARM, then
                 repacks the  primitives. This function
                 is adapted from the HARM function provided by IllinoisGRMHD. The
                 original HARM copyright is included below.

 * Inputs      : params         - grhayl_parameters struct with parameters
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

/**********************************************************************************

  ghl_hybrid_Noble1D():

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


return: i where
        i = 0 -> success
            1 -> initial v^2 < 0 with initial primitive guess;
            2 -> Newton-Raphson solver did not converge to a solution with the
                 given tolerances;
            3 -> Newton-Raphson procedure encountered a numerical divergence
                 (occurrence of "nan" or "+/-inf");
            4 -> Z<0 or Z>Z_TOO_BIG
            5 -> v^2 > 1 returned by the Newton-Raphson solver;
            6 -> rho <= 0 computed by returned quantities; note that this error code
                 is bypassed by the Cupp_fix parameter, known cases of this error
                 are resolved by ghl_enforce_primitive_limits_and_compute_u0()

**********************************************************************************/

int ghl_hybrid_Noble1D(
      const grhayl_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const conservative_quantities *restrict cons_undens,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics ) {

  double gnr_out[NEWT_DIM];

  // Contains Bsq,QdotBsq,Qsq,Qtsq,Qdotn,QdotB,D,W,W_times_S,ye
  harm_aux_vars_struct harm_aux;

  const int ndim = NEWT_DIM;

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  const double BbarU[3] = {prims->BU[0] * ONE_OVER_SQRT_4PI,
                           prims->BU[1] * ONE_OVER_SQRT_4PI,
                           prims->BU[2] * ONE_OVER_SQRT_4PI};
  harm_aux.Bsq = ghl_compute_vec2_from_vecU(ADM_metric->gammaDD, BbarU);

  const double uu = - cons_undens->tau*ADM_metric->lapse
                    - ADM_metric->lapse*cons_undens->rho
                    + ADM_metric->betaU[0]*cons_undens->SD[0]
                    + ADM_metric->betaU[1]*cons_undens->SD[1]
                    + ADM_metric->betaU[2]*cons_undens->SD[2];

  const double QD[4] = {uu,
                        cons_undens->SD[0],
                        cons_undens->SD[1],
                        cons_undens->SD[2]};

  double QU[4]; ghl_raise_vector_4D(metric_aux->g4UU, QD, QU);
  harm_aux.Qsq = 0.0;
  for(int i=0; i<4; i++) harm_aux.Qsq += QD[i]*QU[i] ;

  harm_aux.QdotB = 0. ;
  for(int i=0; i<3; i++) harm_aux.QdotB += QD[i+1]*BbarU[i];
  harm_aux.QdotBsq = harm_aux.QdotB*harm_aux.QdotB;

  // n_{\mu}Q^{\mu} = -alpha Q^{0}, since n_{\mu} = (-alpha,0,0,0)
  harm_aux.Qdotn = -ADM_metric->lapse*QU[0];

  harm_aux.Qtsq = harm_aux.Qsq + harm_aux.Qdotn*harm_aux.Qdotn;

  harm_aux.D    = cons_undens->rho;

  /* calculate Z from last timestep and use for guess */
  const double utU_guess[3] = {prims->vU[0] + ADM_metric->betaU[0],
                               prims->vU[1] + ADM_metric->betaU[1],
                               prims->vU[2] + ADM_metric->betaU[2]};
  const double tmp_u = ghl_compute_vec2_from_vecU(ADM_metric->gammaDD, utU_guess);

  double vsq = tmp_u/(1.0-tmp_u);

  if( (vsq < 0.) && (fabs(vsq) < 1.0e-13) ) {
    vsq = fabs(vsq);
  } else if(vsq < 0.0 || vsq > 10.0) {
    return 1;
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
    const double Gamma_ppoly = eos->Gamma_ppoly[ghl_hybrid_find_polytropic_index(eos, prims->rho)];
    u = prims->press/(Gamma_ppoly - 1.0);
    p = ghl_pressure_rho0_u(eos, rho0, u);
    w = rho0 + u + p;
  } else if(eos->eos_type == grhayl_eos_tabulated) {
    ghl_warn("No tabulated EOS support yet! Sorry!");
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
  // Noble2D has fabs(Z)
  gnr_out[0] = Z_last;

  // To be consistent with entropy variants, unused argument 0.0 is needed
  int retval = ghl_newton_raphson_1d(eos, &harm_aux, ndim, 0.0, &diagnostics->n_iter, gnr_out, ghl_func_1d_orig);

  const double Z = gnr_out[0];

  /* Problem with solver, so return denoting error before doing anything further */
  if(retval != 0) {
    return retval;
  } else if(Z <= 0. || Z > Z_TOO_BIG) {
    return 4;
  }

  // Calculate v^2:
  vsq = ghl_vsq_calc(&harm_aux, Z);
//TODO: differs from Noble2D
  if( vsq >= 1. ) {
    return 5;
  }

  // Recover the primitive variables from the scalars and conserved variables:
  const double gtmp = sqrt(1. - vsq);
  harm_aux.W = 1.0/gtmp;
  w = Z * (1.0 - vsq);

  prims->rho = harm_aux.D * gtmp;

  // Cupp Fix logic:
  // If the returned value is 5, then the Newton-Rapson method converged, but the values were so small
  // that u or rho were negative (usually u). Since the method converged, we only need to fix the values
  // using enforce_primitive_limits_and_output_u0(). There's no need to trigger a Font fix. In my experience,
  // Font Fix returns nearly the same values as this, but takes longer to run (we already did the work for
  // these results, after all!
  // Also note that we have completely eliminated u, so that check doesn't exist any longer.
  if( !params->Cupp_Fix && prims->rho <= 0.0) {
    return 6;
  }

  const double nU[4] = {ADM_metric->lapseinv,
                       -ADM_metric->lapseinv*ADM_metric->betaU[0],
                       -ADM_metric->lapseinv*ADM_metric->betaU[1],
                       -ADM_metric->lapseinv*ADM_metric->betaU[2]};

  double Qtcon[4];
  const double g_o_ZBsq = harm_aux.W/(Z+harm_aux.Bsq);
  const double QdB_o_Z  = harm_aux.QdotB / Z;

  for(int i=1; i<4; i++) Qtcon[i] = QU[i] + nU[i] * harm_aux.Qdotn;

  double utU[3] = {g_o_ZBsq * ( Qtcon[1] + QdB_o_Z*BbarU[0] ),
                   g_o_ZBsq * ( Qtcon[2] + QdB_o_Z*BbarU[1] ),
                   g_o_ZBsq * ( Qtcon[3] + QdB_o_Z*BbarU[2] )};

  //Additional tabulated code here

  ghl_limit_utilde_and_compute_v(eos, ADM_metric, utU, prims, &diagnostics->speed_limited);

  if(diagnostics->speed_limited==1)
    prims->rho = cons_undens->rho/(ADM_metric->lapse*prims->u0);

  if( eos->eos_type == grhayl_eos_hybrid ) {
    prims->press = ghl_pressure_rho0_w(eos, prims->rho, w);
    double P_cold = 0.0;
    double eps_cold = 0.0;
    ghl_hybrid_compute_P_cold_and_eps_cold(eos, prims->rho, &P_cold, &eps_cold);
    prims->eps = eps_cold + (prims->press-P_cold)/(eos->Gamma_th-1.0)/prims->rho;
    if( params->evolve_entropy ) ghl_hybrid_compute_entropy_function(eos, prims->rho, prims->press, &prims->entropy);
  } else if(eos->eos_type == grhayl_eos_tabulated) {
    ghl_warn("No tabulated EOS support yet! Sorry!");
  }

  /* Done! */
  return 0;
}
