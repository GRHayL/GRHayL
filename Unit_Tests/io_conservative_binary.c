#include "unit_tests.h"

void read_conservative_binary(
                     const bool evolve_entropy,
                     const conservative_quantities *restrict cons_orig,
                     const conservative_quantities *restrict cons,
                     FILE *restrict infile) {

  fread(&cons_orig->rho, sizeof(double), 1, infile);
  fread(&cons->rho,      sizeof(double), 1, infile);
  
  fread(&cons_orig->tau, sizeof(double), 1, infile);
  fread(&cons->tau,      sizeof(double), 1, infile);

  fread(&cons_orig->S_x, sizeof(double), 1, infile);
  fread(&cons->S_x,      sizeof(double), 1, infile);

  fread(&cons_orig->S_y, sizeof(double), 1, infile);
  fread(&cons->S_y,      sizeof(double), 1, infile);

  fread(&cons_orig->S_z, sizeof(double), 1, infile);
  fread(&cons->S_z,      sizeof(double), 1, infile);

  if(evolve_entropy) {
    fread(&cons_orig->entropy, sizeof(double), 1, infile);
    fread(&cons->entropy,      sizeof(double), 1, infile);
  }
}

void write_conservative_binary(
                     const bool evolve_entropy,
                     const conservative_quantities *restrict cons_orig,
                     const conservative_quantities *restrict cons,
                     FILE *restrict outfile) {

  fwrite(&cons_orig->rho, sizeof(double), 1, outfile);
  fwrite(&cons->rho,      sizeof(double), 1, outfile);
  
  fwrite(&cons_orig->tau, sizeof(double), 1, outfile);
  fwrite(&cons->tau,      sizeof(double), 1, outfile);

  fwrite(&cons_orig->S_x, sizeof(double), 1, outfile);
  fwrite(&cons->S_x,      sizeof(double), 1, outfile);

  fwrite(&cons_orig->S_y, sizeof(double), 1, outfile);
  fwrite(&cons->S_y,      sizeof(double), 1, outfile);

  fwrite(&cons_orig->S_z, sizeof(double), 1, outfile);
  fwrite(&cons->S_z,      sizeof(double), 1, outfile);

  if(evolve_entropy) {
    fwrite(&cons_orig->entropy, sizeof(double), 1, outfile);
    fwrite(&cons->entropy,      sizeof(double), 1, outfile);
  }
}
