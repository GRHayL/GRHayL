#include "unit_tests.h"

void read_primitive_binary(
      const int eos_type,
      const bool evolve_entropy,
      double *restrict rho,
      double *restrict press,
      double *restrict vx,
      double *restrict vy,
      double *restrict vz,
      double *restrict eps,
      double *restrict Bx,
      double *restrict By,
      double *restrict Bz,
      double *restrict entropy,
      double *restrict Y_e,
      double *restrict temperature,
      FILE *restrict infile) {

  int key;
  key  = fread(rho,   sizeof(double), 1, infile);
  key += fread(press, sizeof(double), 1, infile);
  key += fread(vx,    sizeof(double), 1, infile);
  key += fread(vy,    sizeof(double), 1, infile);
  key += fread(vz,    sizeof(double), 1, infile);
  key += fread(eps,   sizeof(double), 1, infile);
  key += fread(Bx,    sizeof(double), 1, infile);
  key += fread(By,    sizeof(double), 1, infile);
  key += fread(Bz,    sizeof(double), 1, infile);

  if(evolve_entropy)
    key += fread(entropy, sizeof(double), 1, infile);

  if(eos_type == 1) { //Tabulated
    key += fread(Y_e, sizeof(double), 1, infile);
    key += fread(temperature, sizeof(double), 1, infile);
  }

  // Since each read only reads a single double, the key should just be a sum of every read
  // that happens. Hence, 8 for main primitives, +1 for entropy, and
  // 2 for tabulated (Y_e and temperature).
  const int correct_key = 9 + evolve_entropy + (eos_type == 2)*2;
  if( key != correct_key)
    grhayl_error("An error has occured with reading in trusted primitive data. "
                 "Please check that comparison data "
                 "is up-to-date with current test version.\n");
}

void write_primitive_binary(
      const int eos_type,
      const bool evolve_entropy,
      const primitive_quantities *restrict prims,
      FILE *restrict outfile) {

  fwrite(&prims->rho, sizeof(double), 1, outfile);
  fwrite(&prims->press, sizeof(double), 1, outfile);
  fwrite(&prims->vx, sizeof(double), 1, outfile);
  fwrite(&prims->vy, sizeof(double), 1, outfile);
  fwrite(&prims->vz, sizeof(double), 1, outfile);
  fwrite(&prims->eps, sizeof(double), 1, outfile);
  fwrite(&prims->Bx, sizeof(double), 1, outfile);
  fwrite(&prims->By, sizeof(double), 1, outfile);
  fwrite(&prims->Bz, sizeof(double), 1, outfile);

  if(evolve_entropy)
    fwrite(&prims->entropy, sizeof(double), 1, outfile);

  if(eos_type == 1) { //Tabulated
    fwrite(&prims->Y_e, sizeof(double), 1, outfile);
    fwrite(&prims->temperature, sizeof(double), 1, outfile);
  }
}
