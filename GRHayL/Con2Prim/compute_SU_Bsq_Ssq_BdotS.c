#include "ghl.h"

/*
 * Function : compute_BU_SU_Bsq_Ssq_BdotS
 * Author   : Leo Werneck
 *
 * Compute S^{i}, B^2, S^2, B.S = B^{i}S_{i}.
 *
 * Parameters : metric       - Metric quantities
 *            : cons_undens  - Undensitized conservatives
 *            : prims        - Input primitives (for B^{i})
 *            : SU           - Stores S^{i}.
 *            : Bsq          - Stores B^2.
 *            : Ssq          - Stores S^2.
 *            : BdotS        - Stores B.S = B^{i}S_{i}.
 *
 * Returns    : Nothing.
 */
void ghl_compute_SU_Bsq_Ssq_BdotS(
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons_undens,
      const ghl_primitive_quantities *restrict prims,
      double *restrict SU,
      double *restrict Bsq,
      double *restrict Ssq,
      double *restrict BdotS) {

  // Step 1: Compute S^{2} = gamma^{ij}S_{i}S_{j}
  double SD[3] = {cons_undens->SD[0], cons_undens->SD[1], cons_undens->SD[2]};
  double S_squared = ghl_compute_vec2_from_vec3D(metric_adm->gammaUU, SD);

  // Step 2: Enforce ceiling on S^{2} (Eq. A5 of [1])
  // Step 2.1: Compute maximum allowed value for S^{2}
  const double S_squared_max = SQR(cons_undens->tau + cons_undens->rho);
  if(S_squared > S_squared_max) {
    // Step 2.2: Rescale S_{i}
    const double rescale_factor = sqrt(0.9999*S_squared_max/S_squared);
    for(int i=0;i<3;i++)
      SD[i] *= rescale_factor;

    // Step 2.3: Recompute S^{2}
    S_squared = ghl_compute_vec2_from_vec3D(metric_adm->gammaUU, SD);
  }
  *Ssq = S_squared;

  // Step 3: Compute B^{2} = gamma_{ij}B^{i}B^{j}
  *Bsq = ghl_compute_vec2_from_vec3D(metric_adm->gammaDD, prims->BU);

  // Step 4: Compute B.S = B^{i}S_{i}
  *BdotS = prims->BU[0]*SD[0] + prims->BU[1]*SD[1] + prims->BU[2]*SD[2];

  // Step 5: Compute S^{i}
  ghl_raise_lower_vector_3D(metric_adm->gammaUU, SD, SU);
}
