#include "cctk.h"
#include "IGM.h"

void compute_characteristic_speeds(
      const cGH *restrict cctkGH,
      const int flux_dir,
      const eos_parameters *restrict eos,
      const double **in_metric,
      /*const*/ double **in_prims_r,
      /*const*/ double **in_prims_l,
      double *restrict cmin,
      double *restrict cmax) {

  const double poison = 0.0/0.0;

  void (*calculate_characteristic_speed)(const primitive_quantities *restrict prims_r,
                                         const primitive_quantities *restrict prims_l,
                                         struct eos_parameters const *restrict eos,
                                         const metric_quantities *restrict metric_face,
                                         double *cmin, double *cmax);

  // Set function pointer to specific function for a given direction
  switch(flux_dir) {
    case 0:
      calculate_characteristic_speed = &calculate_characteristic_speed_dirn0;
      break;
    case 1:
      calculate_characteristic_speed = &calculate_characteristic_speed_dirn1;
      break;
    case 2:
      calculate_characteristic_speed = &calculate_characteristic_speed_dirn2;
      break;
    default:
      CCTK_VERROR("Warning: invalid flux_dir value (not 0, 1, or 2) has been passed to compute_characteristic_speeds.");
  }

  const int imin = cctkGH->cctk_nghostzones[0];
  const int jmin = cctkGH->cctk_nghostzones[1];
  const int kmin = cctkGH->cctk_nghostzones[2];
  const int imax = cctkGH->cctk_lsh[0]-(cctkGH->cctk_nghostzones[0]-1);
  const int jmax = cctkGH->cctk_lsh[1]-(cctkGH->cctk_nghostzones[1]-1);
  const int kmax = cctkGH->cctk_lsh[2]-(cctkGH->cctk_nghostzones[2]-1);

#pragma omp parallel for
  for(int k=kmin; k<kmax; k++)
    for(int j=jmin; j<jmax; j++)
      for(int i=imin; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        metric_quantities metric_face;
        primitive_quantities prims_r, prims_l;
        interpolate_to_face_and_initialize_metric(
                          cctkGH, i, j, k,
                          flux_dir, in_metric[LAPSE],
                          in_metric[BETAX], in_metric[BETAY], in_metric[BETAZ],
                          in_metric[GXX], in_metric[GXY], in_metric[GXZ],
                          in_metric[GYY], in_metric[GYZ], in_metric[GZZ],
                          &metric_face);

        initialize_primitives(in_prims_r[RHOB][index], in_prims_r[PRESSURE][index], poison,
                              in_prims_r[VX][index], in_prims_r[VY][index], in_prims_r[VZ][index],
                              in_prims_r[BX_CENTER][index], in_prims_r[BY_CENTER][index], in_prims_r[BZ_CENTER][index],
                              poison, poison, poison, // entropy, Y_e, temp
                              &prims_r);

        initialize_primitives(in_prims_l[RHOB][index], in_prims_l[PRESSURE][index], poison,
                              in_prims_l[VX][index], in_prims_l[VY][index], in_prims_l[VZ][index],
                              in_prims_l[BX_CENTER][index], in_prims_l[BY_CENTER][index], in_prims_l[BZ_CENTER][index],
                              poison, poison, poison, // entropy, Y_e, temp
                              &prims_l);

        int speed_limited = poison;
        limit_v_and_compute_u0(eos, &metric_face, &prims_r, &speed_limited);
        limit_v_and_compute_u0(eos, &metric_face, &prims_l, &speed_limited);

        calculate_characteristic_speed(&prims_r, &prims_l, eos, &metric_face, &cmin[index], &cmax[index]);
  }
}
