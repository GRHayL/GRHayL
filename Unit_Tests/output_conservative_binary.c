#include "unit_tests.h"

void output_conservative_binary(
                     const bool evolve_entropy,
                     const conservative_quantities *restrict cons_orig,
                     const conservative_quantities *restrict cons,
                     FILE *restrict outfile) {

  conservative_quantities cons_error;

  cons_error.rho = relative_error(cons->rho, cons_orig->rho);
  cons_error.tau = relative_error(cons->tau, cons_orig->tau);
  cons_error.S_x = relative_error(cons->S_x, cons_orig->S_x);
  cons_error.S_y = relative_error(cons->S_y, cons_orig->S_y);
  cons_error.S_z = relative_error(cons->S_z, cons_orig->S_z);

  fwrite(&cons_orig->rho, sizeof(double), 1, outfile);
  fwrite(&cons->rho,      sizeof(double), 1, outfile);
  fwrite(&cons_error.rho, sizeof(double), 1, outfile);
  
  fwrite(&cons_orig->tau, sizeof(double), 1, outfile);
  fwrite(&cons->tau,      sizeof(double), 1, outfile);
  fwrite(&cons_error.tau, sizeof(double), 1, outfile);

  fwrite(&cons_orig->S_x, sizeof(double), 1, outfile);
  fwrite(&cons->S_x,      sizeof(double), 1, outfile);
  fwrite(&cons_error.S_x, sizeof(double), 1, outfile);

  fwrite(&cons_orig->S_y, sizeof(double), 1, outfile);
  fwrite(&cons->S_y,      sizeof(double), 1, outfile);
  fwrite(&cons_error.S_y, sizeof(double), 1, outfile);

  fwrite(&cons_orig->S_z, sizeof(double), 1, outfile);
  fwrite(&cons->S_z,      sizeof(double), 1, outfile);
  fwrite(&cons_error.S_z, sizeof(double), 1, outfile);

  if(evolve_entropy) {
    fwrite(&cons_orig->rho, sizeof(double), 1, outfile);
    fwrite(&cons->rho,      sizeof(double), 1, outfile);
    fwrite(&cons_error.rho, sizeof(double), 1, outfile);
  }
}
