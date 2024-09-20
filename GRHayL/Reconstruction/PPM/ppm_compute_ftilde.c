#include "ghl_reconstruction.h"

void ghl_compute_ftilde(
      const ghl_parameters *restrict params,
      const double pressure[6],
      const double v_flux_dirn[6],
      double ftilde[2]) {

  ftilde[0] = ghl_shock_detection_ftilde(params, pressure, v_flux_dirn);
  ftilde[1] = ghl_shock_detection_ftilde(params, &pressure[1], &v_flux_dirn[1]);
}
