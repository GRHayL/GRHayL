#include "unit_tests.h"

void read_conservative_binary(
                     const bool evolve_entropy,
                     conservative_quantities *restrict cons_orig,
                     conservative_quantities *restrict cons,
                     FILE *restrict infile) {

  int key;
  key = fread(&cons_orig->rho, sizeof(double), 1, infile);
  key += fread(&cons->rho,      sizeof(double), 1, infile);

  key += fread(&cons_orig->tau, sizeof(double), 1, infile);
  key += fread(&cons->tau,      sizeof(double), 1, infile);

  key += fread(&cons_orig->S_x, sizeof(double), 1, infile);
  key += fread(&cons->S_x,      sizeof(double), 1, infile);

  key += fread(&cons_orig->S_y, sizeof(double), 1, infile);
  key += fread(&cons->S_y,      sizeof(double), 1, infile);

  key += fread(&cons_orig->S_z, sizeof(double), 1, infile);
  key += fread(&cons->S_z,      sizeof(double), 1, infile);

  if(evolve_entropy) {
    key += fread(&cons_orig->entropy, sizeof(double), 1, infile);
    key += fread(&cons->entropy,      sizeof(double), 1, infile);
  }

  // Since each read only reads a single double, the key should just be a sum of every read
  // that happens. Hence, 5 variables and +1 for entropy. This is then doubled, as it is
  // reading cons and cons_orig.
  const int correct_key = 2*(5 + evolve_entropy);
  if( key != correct_key) {
    printf("An error has occured with reading in trusted conservative data."
           "Please check that comparison data"
           "is up-to-date with current test version.\n");
    exit(1);
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
