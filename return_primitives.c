#include "con2prim_header.h"
#include "cctk.h"

/* This function fills the struct primitive_quantities with data. 
   For more information on the arguments, see the definition of the
   struct in new_header.h. */
void return_primitives(const eos_parameters *restrict eos, const primitive_quantities *restrict prims,
                      double *restrict rho, double *restrict press, double *restrict epsilon,
                      double *restrict vx, double *restrict vy, double *restrict vz,
                      double *restrict Bx, double *restrict By, double *restrict Bz,
                      double *restrict entropy, double *restrict Y_e, double *restrict temp) {
//  CCTK_VINFO("Setting rho=%f to prims->rho=%f",*rho,prims->rho);
//  CCTK_VINFO("Setting press=%f to prims->press=%f",*press,prims->press);
//  CCTK_VINFO("Setting vx=%f to prims->vx=%f",*vx,prims->vx);
//  CCTK_VINFO("Setting vy=%f to prims->vy=%f",*vy,prims->vy);
//  CCTK_VINFO("Setting vz=%f to prims->vz=%f",*vz,prims->vz);
//  CCTK_VINFO("Setting epsilon=%f to prims->eps=%f",*epsilon,prims->eps);

  *rho = prims->rho;
  *press = prims->press;
  *vx = prims->vx;
  *vy = prims->vy;
  *vz = prims->vz;
  *epsilon = prims->eps;
  *entropy = prims->entropy;
  // Tabulated EOS quantities
  if( eos->eos_type == 1 ) {
    *Y_e = prims->Y_e;
    *temp = prims->temp;
  }
}
