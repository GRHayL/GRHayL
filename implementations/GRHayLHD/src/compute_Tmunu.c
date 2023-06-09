#include "GRHayLHD.h"

void GRHayLHD_compute_Tmunu(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLHD_compute_Tmunu;
  DECLARE_CCTK_PARAMETERS;

  const double poison = 0.0/0.0;

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // Read in ADM metric quantities from gridfunctions and
        // set auxiliary and ADM metric quantities
        metric_quantities ADM_metric;
        ghl_enforce_detgtij_and_initialize_ADM_metric(
              alp[index],
              betax[index], betay[index], betaz[index],
              gxx[index], gxy[index], gxz[index],
              gyy[index], gyz[index], gzz[index],
              &ADM_metric);

        ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        // Read in primitive variables from gridfunctions
        primitive_quantities prims;
        ghl_initialize_primitives(
              rho_b[index], pressure[index], eps[index],
              vx[index], vy[index], vz[index],
              0.0, 0.0, 0.0,
              poison, poison, poison, &prims);

        prims.u0 = u0[index];

        stress_energy Tmunu;
        ghl_compute_TDNmunu(&ADM_metric, &metric_aux, &prims, &Tmunu);

        ghl_return_stress_energy(
              &Tmunu, &eTtt[index], &eTtx[index],
              &eTty[index], &eTtz[index], &eTxx[index],
              &eTxy[index], &eTxz[index], &eTyy[index],
              &eTyz[index], &eTzz[index]);
      }
    }
  }
}
