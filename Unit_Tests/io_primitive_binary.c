#include "unit_tests.h"

void read_primitive_binary(
      const int eos_type,
      const bool velocity_only,
      const bool evolve_entropy,
      const primitive_quantities *restrict prims_orig,
      const primitive_quantities *restrict prims,
      FILE *restrict infile) {

  int key;
  key = fread(&prims_orig->vx, sizeof(double), 1, infile);
  key = fread(&prims->vx,      sizeof(double), 1, infile);

  key = fread(&prims_orig->vy, sizeof(double), 1, infile);
  key = fread(&prims->vy,      sizeof(double), 1, infile);

  key = fread(&prims_orig->vz, sizeof(double), 1, infile);
  key = fread(&prims->vz,      sizeof(double), 1, infile);

  if(!velocity_only) {
    key = fread(&prims_orig->rho, sizeof(double), 1, infile);
    key = fread(&prims->rho,      sizeof(double), 1, infile);

    key = fread(&prims_orig->press, sizeof(double), 1, infile);
    key = fread(&prims->press,      sizeof(double), 1, infile);

    if(evolve_entropy) {
      key = fread(&prims_orig->entropy, sizeof(double), 1, infile);
      key = fread(&prims->entropy,      sizeof(double), 1, infile);
    }
    if(eos_type == 2) { //Tabulated
      key = fread(&prims_orig->Y_e, sizeof(double), 1, infile);
      key = fread(&prims->Y_e,      sizeof(double), 1, infile);

      key = fread(&prims_orig->temperature, sizeof(double), 1, infile);
      key = fread(&prims->temperature,      sizeof(double), 1, infile);
    }
  }
}

void write_primitive_binary(
      const int eos_type,
      const bool velocity_only,
      const bool evolve_entropy,
      const primitive_quantities *restrict prims_orig,
      const primitive_quantities *restrict prims,
      FILE *restrict outfile) {

  fwrite(&prims_orig->vx, sizeof(double), 1, outfile);
  fwrite(&prims->vx,      sizeof(double), 1, outfile);

  fwrite(&prims_orig->vy, sizeof(double), 1, outfile);
  fwrite(&prims->vy,      sizeof(double), 1, outfile);

  fwrite(&prims_orig->vz, sizeof(double), 1, outfile);
  fwrite(&prims->vz,      sizeof(double), 1, outfile);

  if(!velocity_only) {
    fwrite(&prims_orig->rho, sizeof(double), 1, outfile);
    fwrite(&prims->rho,      sizeof(double), 1, outfile);

    fwrite(&prims_orig->press, sizeof(double), 1, outfile);
    fwrite(&prims->press,      sizeof(double), 1, outfile);

    if(evolve_entropy) {
      fwrite(&prims_orig->entropy, sizeof(double), 1, outfile);
      fwrite(&prims->entropy,      sizeof(double), 1, outfile);
    }
    if(eos_type == 2) { //Tabulated
      fwrite(&prims_orig->Y_e, sizeof(double), 1, outfile);
      fwrite(&prims->Y_e,      sizeof(double), 1, outfile);

      fwrite(&prims_orig->temperature, sizeof(double), 1, outfile);
      fwrite(&prims->temperature,      sizeof(double), 1, outfile);
    }
  }
}
