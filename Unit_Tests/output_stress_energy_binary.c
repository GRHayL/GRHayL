#include "unit_tests.h"

void output_stress_energy_binary(
                     const stress_energy *restrict Tmunu_orig,
                     const stress_energy *restrict Tmunu,
                     FILE *restrict outfile) {

  stress_energy Tmunu_error;

  Tmunu_error.Ttt = relative_error(Tmunu->Ttt, Tmunu_orig->Ttt);
  Tmunu_error.Ttx = relative_error(Tmunu->Ttx, Tmunu_orig->Ttx);
  Tmunu_error.Tty = relative_error(Tmunu->Tty, Tmunu_orig->Tty);
  Tmunu_error.Ttz = relative_error(Tmunu->Ttz, Tmunu_orig->Ttz);
  Tmunu_error.Txx = relative_error(Tmunu->Txx, Tmunu_orig->Txx);
  Tmunu_error.Txy = relative_error(Tmunu->Txy, Tmunu_orig->Txy);
  Tmunu_error.Txz = relative_error(Tmunu->Txz, Tmunu_orig->Txz);
  Tmunu_error.Tyy = relative_error(Tmunu->Tyy, Tmunu_orig->Tyy);
  Tmunu_error.Tyz = relative_error(Tmunu->Tyz, Tmunu_orig->Tyz);
  Tmunu_error.Tzz = relative_error(Tmunu->Tzz, Tmunu_orig->Tzz);

  fwrite(&Tmunu_orig->Ttt, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Ttt,      sizeof(double), 1, outfile);
  fwrite(&Tmunu_error.Ttt, sizeof(double), 1, outfile);

  fwrite(&Tmunu_orig->Ttx, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Ttx,      sizeof(double), 1, outfile);
  fwrite(&Tmunu_error.Ttx, sizeof(double), 1, outfile);

  fwrite(&Tmunu_orig->Tty, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Tty,      sizeof(double), 1, outfile);
  fwrite(&Tmunu_error.Tty, sizeof(double), 1, outfile);

  fwrite(&Tmunu_orig->Ttz, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Ttz,      sizeof(double), 1, outfile);
  fwrite(&Tmunu_error.Ttz, sizeof(double), 1, outfile);

  fwrite(&Tmunu_orig->Txx, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Txx,      sizeof(double), 1, outfile);
  fwrite(&Tmunu_error.Txx, sizeof(double), 1, outfile);

  fwrite(&Tmunu_orig->Txy, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Txy,      sizeof(double), 1, outfile);
  fwrite(&Tmunu_error.Txy, sizeof(double), 1, outfile);

  fwrite(&Tmunu_orig->Txz, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Txz,      sizeof(double), 1, outfile);
  fwrite(&Tmunu_error.Txz, sizeof(double), 1, outfile);

  fwrite(&Tmunu_orig->Tyy, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Tyy,      sizeof(double), 1, outfile);
  fwrite(&Tmunu_error.Tyy, sizeof(double), 1, outfile);

  fwrite(&Tmunu_orig->Tyz, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Tyz,      sizeof(double), 1, outfile);
  fwrite(&Tmunu_error.Tyz, sizeof(double), 1, outfile);

  fwrite(&Tmunu_orig->Tzz, sizeof(double), 1, outfile);
  fwrite(&Tmunu->Tzz,      sizeof(double), 1, outfile);
  fwrite(&Tmunu_error.Tzz, sizeof(double), 1, outfile);
}
