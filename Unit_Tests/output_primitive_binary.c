#include "unit_tests.h"

void output_primitive_binary(
                     const int eos_type,
                     const bool velocity_only,
                     const bool evolve_entropy,
                     const primitive_quantities *restrict prims_orig, 
                     const primitive_quantities *restrict prims, 
                     FILE *restrict outfile) {

  primitive_quantities prims_error;
  prims_error.rho         = relative_error(prims->rho,         prims_orig->rho);
  prims_error.press       = relative_error(prims->press,       prims_orig->press);
  prims_error.eps         = relative_error(prims->eps,         prims_orig->eps);
  prims_error.vx          = relative_error(prims->vx,          prims_orig->vx);
  prims_error.vy          = relative_error(prims->vy,          prims_orig->vy);
  prims_error.vz          = relative_error(prims->vz,          prims_orig->vz);
  prims_error.Bx          = relative_error(prims->Bx,          prims_orig->Bx);
  prims_error.By          = relative_error(prims->By,          prims_orig->By);
  prims_error.Bz          = relative_error(prims->Bz,          prims_orig->Bz);
  prims_error.entropy     = relative_error(prims->entropy,     prims_orig->entropy);
  prims_error.Y_e         = relative_error(prims->Y_e,         prims_orig->Y_e);
  prims_error.temperature = relative_error(prims->temperature, prims_orig->temperature);

  fwrite(&prims_orig->vx, sizeof(double), 1, outfile);
  fwrite(&prims->vx,      sizeof(double), 1, outfile);
  fwrite(&prims_error.vx, sizeof(double), 1, outfile);

  fwrite(&prims_orig->vy, sizeof(double), 1, outfile);
  fwrite(&prims->vy,      sizeof(double), 1, outfile);
  fwrite(&prims_error.vy, sizeof(double), 1, outfile);

  fwrite(&prims_orig->vz, sizeof(double), 1, outfile);
  fwrite(&prims->vz,      sizeof(double), 1, outfile);
  fwrite(&prims_error.vz, sizeof(double), 1, outfile);

  if(!velocity_only) {
    fwrite(&prims_orig->rho, sizeof(double), 1, outfile);
    fwrite(&prims->rho,      sizeof(double), 1, outfile);
    fwrite(&prims_error.rho, sizeof(double), 1, outfile);
  
    fwrite(&prims_orig->press, sizeof(double), 1, outfile);
    fwrite(&prims->press,      sizeof(double), 1, outfile);
    fwrite(&prims_error.press, sizeof(double), 1, outfile);

    if(evolve_entropy) {
      fwrite(&prims_orig->entropy, sizeof(double), 1, outfile);
      fwrite(&prims->entropy,      sizeof(double), 1, outfile);
      fwrite(&prims_error.entropy, sizeof(double), 1, outfile);
    }
    if(eos_type == 2) { //Tabulated
      fwrite(&prims_orig->Y_e, sizeof(double), 1, outfile);
      fwrite(&prims->Y_e,      sizeof(double), 1, outfile);
      fwrite(&prims_error.Y_e, sizeof(double), 1, outfile);

      fwrite(&prims_orig->temperature, sizeof(double), 1, outfile);
      fwrite(&prims->temperature,      sizeof(double), 1, outfile);
      fwrite(&prims_error.temperature, sizeof(double), 1, outfile);
    }
  }
}
