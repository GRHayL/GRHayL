#include "con2prim.h"

/* Function    : ghl_limit_utilde_and_compute_v()
 * Description : Applies speed limit to \tilde{u}^i and computes v^i and u^0
 *
 * Inputs      : eos            - ghl_eos_parameters struct with data for the
 *                                EOS of the simulation
 *             : metric         - ghl_metric_quantities struct with data for
 *                                the gridpoint of interest
 *             : utcon1_ptr     - pointer to the x component of \tilde{u}^1
 *             : utcon2_ptr     - pointer to the y component of \tilde{u}^2
 *             : utcon3_ptr     - pointer to the z component of \tilde{u}^3
 *
 * Outputs     : utcon1_ptr     - returns velocity-limited \tilde{u}^1
 *             : utcon2_ptr     - returns velocity-limited \tilde{u}^2
 *             : utcon3_ptr     - returns velocity-limited \tilde{u}^3
 *             : prims          - returns prims->v^i and prims->u0 computed from
 *                                velocity-limited \tilde{u}^i
 *             : diagnostics    - tracks if the velocity was limited
 *
 */

//Now that we have found some solution, we first limit velocity:
//FIXME: Probably want to use exactly the same velocity limiter function here as in mhdflux.C
void ghl_limit_utilde_and_compute_v(
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      double utU[3],
      ghl_primitive_quantities *restrict prims,
      int *restrict speed_limit) {

  //Velocity limiter:
  double ut2 = ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, utU);
  double au0m1 = ut2/( 1.0+sqrt(1.0+ut2) );

  (*speed_limit) = 0;
  // *** Limit velocity
  if (au0m1 > 0.9999999*(eos->W_max-1.0)) {
    double fac = sqrt((SQR(eos->W_max)-1.0)/(SQR(1.0+au0m1) - 1.0));
    utU[0] *= fac;
    utU[1] *= fac;
    utU[2] *= fac;
    ut2 = ut2 * SQR(fac);
    au0m1 = ut2/( 1.0+sqrt(1.0+ut2) );
    (*speed_limit) = 1;
  } //Finished limiting velocity

  // Calculate v^i and u^0 from \tilde{u}^i
  prims->u0 = (au0m1+1.0)*ADM_metric->lapseinv;
  prims->vU[0] = utU[0]/prims->u0 - ADM_metric->betaU[0];
  prims->vU[1] = utU[1]/prims->u0 - ADM_metric->betaU[1];
  prims->vU[2] = utU[2]/prims->u0 - ADM_metric->betaU[2];
}
