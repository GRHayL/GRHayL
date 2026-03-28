#include "ghl_reconstruction.h"

void ghl_weno7z_reconstruction(
      const double U[8],
      double *restrict Ur,
      double *restrict Ul) {

  double tmpr, tmpl;

  ghl_weno7z_reconstruction_for_cell(&U[0], &tmpr, &tmpl);
  *Ul = tmpr;

  ghl_weno7z_reconstruction_for_cell(&U[1], &tmpr, &tmpl);
  *Ur = tmpl;
}
