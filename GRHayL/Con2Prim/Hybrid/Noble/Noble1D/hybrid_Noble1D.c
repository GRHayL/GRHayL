#include "../../../utils_Noble.h"

/* Function    :  Hybrid_Noble1D()
 * Description :  Unpacks the ghl_primitive_quantities struct into the variables
 *                needed by the Newton-Rapson solver, then repacks the  primitives.
 *                This function is adapted from the HARM function provided by IllinoisGRMHD.
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_hybrid_Noble1D
*/

/***********************************************************************************
******************************* HARM License ***************************************
************************************************************************************

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

int ghl_hybrid_Noble1D(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  double gnr_out[1];

  harm_aux_vars_struct harm_aux;

  // Calculate Z:
  if( ghl_initialize_Noble(params, eos, ADM_metric, metric_aux, cons_undens,
                           prims, &harm_aux, &gnr_out[0]) )
    return 1;

  // To be consistent with entropy variants, unused argument 0.0 is needed
  const int retval = ghl_general_newton_raphson(eos, &harm_aux, 1, 0.0, gnr_out, ghl_validate_1D, ghl_func_1D);

  const double Z = gnr_out[0];

  /* Problem with solver, so return denoting error before doing anything further */
  if(retval != 0) {
    return retval;
  } else if(Z <= 0. || Z > 1e20) {
    return 4;
  }

  // Calculate v^2:
  double vsq = ghl_vsq_calc(&harm_aux, Z);
  if(vsq >= 1.0) {
    vsq = 1.-2.e-16;
  } else if(vsq < 0.0) {
    //v should be real!
    return 5;
  }

  // Recover the primitive variables from the scalars and conserved variables:
  ghl_finalize_Noble(params, eos, ADM_metric, metric_aux, cons_undens, &harm_aux, Z, vsq, prims);
  if(!params->ignore_negative_pressure && prims->press <= 0.0) {
    return 6;
  }

  /* Done! */
  diagnostics->n_iter = harm_aux.n_iter;
  diagnostics->which_routine = Noble1D;
  return 0;
}
