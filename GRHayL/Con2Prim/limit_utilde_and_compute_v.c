#include "con2prim.h"

/* Function    : limit_utilde_and_compute_v()
 * Description : Applies speed limit to \tilde{u}^i and computes v^i and u^0
 *
 * Inputs      : eos            - eos_parameters struct with data for the
 *                                EOS of the simulation
 *             : metric         - metric_quantities struct with data for
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
void limit_utilde_and_compute_v(
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      double *restrict utcon1_ptr,
      double *restrict utcon2_ptr,
      double *restrict utcon3_ptr,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics ) {

  double utcon1 = *utcon1_ptr;
  double utcon2 = *utcon2_ptr;
  double utcon3 = *utcon3_ptr;

  //Velocity limiter:
  double gijuiuj = metric->adm_gxx*SQR(utcon1 ) +
    2.0*metric->adm_gxy*utcon1*utcon2 + 2.0*metric->adm_gxz*utcon1*utcon3 +
    metric->adm_gyy*SQR(utcon2) + 2.0*metric->adm_gyz*utcon2*utcon3 +
    metric->adm_gzz*SQR(utcon3);
  double au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
  prims->u0 = (au0m1+1.0)*metric->lapseinv;

  diagnostics->vel_limited_ptcount=0;
  // *** Limit velocity
  if (au0m1 > 0.9999999*(eos->W_max-1.0)) {
    double fac = sqrt((SQR(eos->W_max)-1.0)/(SQR(1.0+au0m1) - 1.0));
    utcon1 *= fac;
    utcon2 *= fac;
    utcon3 *= fac;
    gijuiuj = gijuiuj * SQR(fac);
    au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
    // Reset rho_b and u0
    prims->u0 = (au0m1+1.0)*metric->lapseinv;
    diagnostics->vel_limited_ptcount=1;
    diagnostics->failure_checker+=1000;
  } //Finished limiting velocity

  *utcon1_ptr = utcon1;
  *utcon2_ptr = utcon2;
  *utcon3_ptr = utcon3;

  // Calculate v^i from \tilde{u}^i
  prims->vx = utcon1/prims->u0 - metric->betax;
  prims->vy = utcon2/prims->u0 - metric->betay;
  prims->vz = utcon3/prims->u0 - metric->betaz;
}