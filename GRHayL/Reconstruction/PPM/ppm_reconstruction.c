#include "ghl_reconstruction.h"

void ghl_ppm_reconstruction(
      const double ftilde[2],
      const double var_data[6],
      double *restrict var_datar,
      double *restrict var_datal) {

    double tmpr, tmpl;

  ghl_ppm_compute_for_cell(ftilde[0], var_data, &tmpr, &tmpl);
  *var_datal = tmpr;

  ghl_ppm_compute_for_cell(ftilde[1], &var_data[1], &tmpr, &tmpl);
  *var_datar = tmpl;
}
