#include "../../../utils_Noble.h"

/* Function    :  Hybrid_Noble1D_entropy()
 * Description :  Unpacks the ghl_primitive_quantities struct into the variables
 *                needed by the Newton-Rapson solver, then repacks the  primitives.
 *                This function is adapted from the HARM function provided by IllinoisGRMHD.
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_hybrid_Noble1D_entropy
*/

/***********************************************************************************
******************************* HARM License ***************************************
************************************************************************************

    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble, 
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

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

ghl_error_codes_t ghl_hybrid_Noble1D_entropy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  double gnr_out[1];

  harm_aux_vars_struct harm_aux;

  double rho0, Z_last;
  if( ghl_initialize_Noble_entropy(params, eos, ADM_metric, metric_aux,
                                   cons_undens, prims, &harm_aux, &rho0, &Z_last) )
    return 1;

  // Make sure that Z is large enough so that v^2 < 1 :
  int i_increase = 0;
  while( (( Z_last*Z_last*Z_last * ( Z_last + 2.*harm_aux.Bsq )
            - harm_aux.QdotBsq*(2.*Z_last + harm_aux.Bsq) ) <= Z_last*Z_last*(harm_aux.Qtsq-harm_aux.Bsq*harm_aux.Bsq))
         && (i_increase < 10) ) {
    Z_last *= 10.;
    i_increase++;
  }

  gnr_out[0] = Z_last;

  ghl_error_codes_t retval = ghl_general_newton_raphson(eos, &harm_aux, 1, rho0, gnr_out, ghl_validate_1D_entropy, ghl_func_Z);

  const double Z = gnr_out[0];

  /* Problem with solver, so return denoting error before doing anything further */
  if(retval != 0) {
    return retval;
  } else if(Z <= 0. || Z > 1e20) {
    return 4;
  }

  const int n_iter = harm_aux.n_iter;

  double rho_g    =  rho0;
  gnr_out[0] = rho0;

  int ntries = 0;
  while (
    (retval = ghl_general_newton_raphson(eos, &harm_aux, 1, Z, gnr_out, ghl_validate_1D_entropy, ghl_func_rho))
    && (ntries++ < 10) ) {
    rho_g *= 10.0;
    gnr_out[0] = rho_g;
  }

  if(retval != 0) {
    return retval;
  }

  // Combine count for both loops
  harm_aux.n_iter += n_iter;

  // Calculate v^2:
  rho0       = gnr_out[0];

  const double rel_err = (harm_aux.D != 0.0) ? fabs((harm_aux.D-rho0)/harm_aux.D) :
                         (     (rho0 != 0.0) ? fabs((harm_aux.D-rho0)/rho0) : 0.0);
  const double utsq = (rel_err > 1e-15) ? (harm_aux.D-rho0)*(harm_aux.D+rho0)/(rho0*rho0) : 0.0;

  if(utsq < 0.0) {
    return ghl_error_neg_vsq;
  }

  // Recover the primitive variables from the scalars and conserved variables:
  const double Wsq = 1.0+utsq;
  const double W   = sqrt(Wsq);

  prims->rho = rho0;

  if(prims->rho <= 0.0) {
    return ghl_error_neg_rho;
  }

  diagnostics->speed_limited = ghl_finalize_Noble_entropy(params, eos, ADM_metric, metric_aux, cons_undens, &harm_aux, Z, W, prims);
  if(prims->press <= 0.0) {
    return ghl_error_neg_pressure;
  }

  /* Done! */
  diagnostics->n_iter = harm_aux.n_iter;
  diagnostics->which_routine = Noble1D_entropy;
  return ghl_success;
}
