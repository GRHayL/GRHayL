#include "ghl_con2prim.h"

/* Function    : ghl_limit_utilde_and_compute_v()
 * Description : Applies speed limit to \tilde{u}^i and computes v^i and u^0
*/

//Now that we have found some solution, we first limit velocity:
//FIXME: Probably want to use exactly the same velocity limiter function here as in mhdflux.C
bool ghl_limit_utilde_and_compute_v(
      const ghl_parameters *restrict params,
      const ghl_metric_quantities *restrict ADM_metric,
      double utU[3],
      ghl_primitive_quantities *restrict prims) {

  bool speed_limited = false;
  //Velocity limiter:
  double ut2 = ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, utU);
  double au0m1 = ut2/( 1.0+sqrt(1.0+ut2) );

  // *** Limit velocity
  if (au0m1 > 0.9999999*(params->max_Lorentz_factor-1.0)) {
    double fac = sqrt((SQR(params->max_Lorentz_factor)-1.0)/(SQR(1.0+au0m1) - 1.0));
    utU[0] *= fac;
    utU[1] *= fac;
    utU[2] *= fac;
    ut2 = ut2 * SQR(fac);
    au0m1 = ut2/( 1.0+sqrt(1.0+ut2) );
    speed_limited = true;
  } //Finished limiting velocity

  // Calculate v^i and u^0 from \tilde{u}^i
  prims->u0 = (au0m1+1.0)*ADM_metric->lapseinv;
  prims->vU[0] = utU[0]/prims->u0 - ADM_metric->betaU[0];
  prims->vU[1] = utU[1]/prims->u0 - ADM_metric->betaU[1];
  prims->vU[2] = utU[2]/prims->u0 - ADM_metric->betaU[2];
  return speed_limited;
}
