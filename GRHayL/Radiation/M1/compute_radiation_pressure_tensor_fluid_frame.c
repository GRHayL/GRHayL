#include "ghl.h"
#include "ghl_m1.h"

// This function computes Eq. (6) of Radice et al. (2022)
void ghl_radiation_compute_pressure_tensor_fluid_frame(
      const ghl_ADM_aux_quantities *adm_aux,
      const ghl_primitive_quantities *prims,
      const double J,
      ghl_radiation_pressure_tensor *K) {

  // u^{i} = v^{i} * u^{0}
  const double u4U[4] = { prims->u0, prims->u0 * prims->vU[0], prims->u0 * prims->vU[1],
                          prims->u0 * prims->vU[2] };

  double u4D[4] = { 0, 0, 0, 0 };
  for(int mu = 0; mu < 4; mu++) {
    for(int nu = 0; nu < 4; nu++) {
      u4D[mu] += adm_aux->g4DD[mu][nu] * u4U[nu];
    }
  }

  for(int mu = 0; mu < 4; mu++) {
    for(int nu = 0; nu < 4; nu++) {
      K->DD[mu][nu] = J * (adm_aux->g4DD[mu][nu] + u4D[mu] * u4D[nu]);
    }
  }
}
