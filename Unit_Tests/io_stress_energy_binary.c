#include "unit_tests.h"

void read_stress_energy_binary(
                     double *restrict Ttt,
                     double *restrict Ttx,
                     double *restrict Tty,
                     double *restrict Ttz,
                     double *restrict Txx,
                     double *restrict Txy,
                     double *restrict Txz,
                     double *restrict Tyy,
                     double *restrict Tyz,
                     double *restrict Tzz,
                     FILE *restrict infile) {

  int key;
  key  = fread(Ttt, sizeof(double), 1, infile);
  key += fread(Ttx, sizeof(double), 1, infile);
  key += fread(Tty, sizeof(double), 1, infile);
  key += fread(Ttz, sizeof(double), 1, infile);
  key += fread(Txx, sizeof(double), 1, infile);
  key += fread(Txy, sizeof(double), 1, infile);
  key += fread(Txz, sizeof(double), 1, infile);
  key += fread(Tyy, sizeof(double), 1, infile);
  key += fread(Tyz, sizeof(double), 1, infile);
  key += fread(Tzz, sizeof(double), 1, infile);

  // Since each read only reads a single double, the key should just be a sum of
  // every read that happens.
  if( key != 10)
    grhayl_error("An error has occured with reading in trusted stress energy tensor data."
                 "Please check that comparison data"
                 "is up-to-date with current test version.\n");
}

void write_stress_energy_binary(
                     const stress_energy *restrict Tmunu,
                     FILE *restrict outfile) {

  fwrite(&Tmunu->Ttt, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Ttx, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Tty, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Ttz, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Txx, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Txy, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Txz, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Tyy, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Tyz, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Tzz, sizeof(double), 1, outfile);
}
