#include "unit_tests.h"

void read_stress_energy_binary(
                     stress_energy *restrict Tmunu,
                     FILE *restrict infile) {

  int key;
  key  = fread(&Tmunu->Ttt,      sizeof(double), 1, infile);
  key += fread(&Tmunu->Ttx,      sizeof(double), 1, infile);
  key += fread(&Tmunu->Tty,      sizeof(double), 1, infile);
  key += fread(&Tmunu->Ttz,      sizeof(double), 1, infile);
  key += fread(&Tmunu->Txx,      sizeof(double), 1, infile);
  key += fread(&Tmunu->Txy,      sizeof(double), 1, infile);
  key += fread(&Tmunu->Txz,      sizeof(double), 1, infile);
  key += fread(&Tmunu->Tyy,      sizeof(double), 1, infile);
  key += fread(&Tmunu->Tyz,      sizeof(double), 1, infile);
  key += fread(&Tmunu->Tzz,      sizeof(double), 1, infile);

  // Since each read only reads a single double, the key should just be a sum of
  // every read that happens.
  if( key != 10) {
    printf("An error has occured with reading in trusted stress energy tensor data."
           "Please check that comparison data"
           "is up-to-date with current test version.\n");
    exit(1);
  }
}

void write_stress_energy_binary(
                     const stress_energy *restrict Tmunu,
                     FILE *restrict outfile) {

  fwrite(&Tmunu->Ttt,      sizeof(double), 1, outfile);
  fwrite(&Tmunu->Ttx,      sizeof(double), 1, outfile);
  fwrite(&Tmunu->Tty,      sizeof(double), 1, outfile);
  fwrite(&Tmunu->Ttz,      sizeof(double), 1, outfile);
  fwrite(&Tmunu->Txx,      sizeof(double), 1, outfile);
  fwrite(&Tmunu->Txy,      sizeof(double), 1, outfile);
  fwrite(&Tmunu->Txz,      sizeof(double), 1, outfile);
  fwrite(&Tmunu->Tyy,      sizeof(double), 1, outfile);
  fwrite(&Tmunu->Tyz,      sizeof(double), 1, outfile);
  fwrite(&Tmunu->Tzz,      sizeof(double), 1, outfile);
}
