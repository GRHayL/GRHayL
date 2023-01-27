#include "unit_tests.h"

void read_conservative_binary(
                     const bool evolve_entropy,
                     double *restrict rho,
                     double *restrict tau,
                     double *restrict S_x,
                     double *restrict S_y,
                     double *restrict S_z,
                     double *restrict entropy,
                     FILE *restrict infile) {

  int key;
  key  = fread(rho, sizeof(double), 1, infile);
  key += fread(tau, sizeof(double), 1, infile);
  key += fread(S_x, sizeof(double), 1, infile);
  key += fread(S_y, sizeof(double), 1, infile);
  key += fread(S_z, sizeof(double), 1, infile);

  if(evolve_entropy)
    key += fread(entropy, sizeof(double), 1, infile);

  // Since each read only reads a single double, the key should just be a sum of every read
  // that happens. Hence, 5 variables and +1 for entropy.
  if( key != 5 + evolve_entropy)
    grhayl_error("An error has occured with reading in trusted conservative data. "
                 "Please check that comparison data "
                 "is up-to-date with current test version.\n");
}

void write_conservative_binary(
                     const bool evolve_entropy,
                     const conservative_quantities *restrict cons,
                     FILE *restrict outfile) {

  fwrite(&cons->rho, sizeof(double), 1, outfile);
  fwrite(&cons->tau, sizeof(double), 1, outfile);
  fwrite(&cons->S_x, sizeof(double), 1, outfile);
  fwrite(&cons->S_y, sizeof(double), 1, outfile);
  fwrite(&cons->S_z, sizeof(double), 1, outfile);

  if(evolve_entropy)
    fwrite(&cons->entropy, sizeof(double), 1, outfile);
}
