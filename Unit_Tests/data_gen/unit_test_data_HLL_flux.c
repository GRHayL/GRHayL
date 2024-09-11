#include "ghl_unit_tests.h"

int main(int argc, char **argv) {

  const int dirlength = 20;
  const int arraylength = dirlength*dirlength*dirlength;

  double *phi_bssn = (double*) malloc(sizeof(double)*arraylength);

  double *cmin[3];
  double *cmax[3];
  for(int i=0; i<3; i++) {
    cmin[i] = (double*) malloc(sizeof(double)*arraylength);
    cmax[i] = (double*) malloc(sizeof(double)*arraylength);
  }

  double *vrr[3];
  double *vrl[3];
  double *vlr[3];
  double *vll[3];
  for(int i=0; i<3; i++) {
    vrr[i] = (double*) malloc(sizeof(double)*arraylength);
    vrl[i] = (double*) malloc(sizeof(double)*arraylength);
    vlr[i] = (double*) malloc(sizeof(double)*arraylength);
    vll[i] = (double*) malloc(sizeof(double)*arraylength);
  }

  double *Br[3];
  double *Bl[3];
  for(int i=0; i<3; i++) {
    Br[i] = (double*) malloc(sizeof(double)*arraylength);
    Bl[i] = (double*) malloc(sizeof(double)*arraylength);
  }

  double *A_rhs[3];
  for(int i=0; i<3; i++) {
    A_rhs[i] = (double*) malloc(sizeof(double)*arraylength);
  }

  // Initialize random data. Note that for this test,
  // we needn't worry too much with physical reasonableness.
  // Discontinuities or v>1 might give physically strange
  // output, but as long as cmin+cmax != 0 the function will
  // produce numbers.
  for(int k=1; k<dirlength; k++)
    for(int j=1; j<dirlength; j++)
      for(int i=1; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);

        phi_bssn[index] = randf(0.0,5.0);

        for(int coord=0; coord<3; coord++) {
          Br[coord][index]  = randf(-1.0,1.0);
          Bl[coord][index]  = randf(-1.0,1.0);
          vrr[coord][index] = randf(-1.0,1.0);
          vrl[coord][index] = randf(-1.0,1.0);
          vlr[coord][index] = randf(-1.0,1.0);
          vll[coord][index] = randf(-1.0,1.0);

          bool zero = true;
          while(zero) {
            cmin[coord][index] = randf(-10.0,10.0);
            cmax[coord][index] = randf(-10.0,10.0);
            zero = cmin[coord][index]+cmax[coord][index] == 0.0;
          }
          A_rhs[coord][index] = 0.0;
        }
  }

  FILE* outfile;
  outfile = fopen_with_check("HLL_flux_input.bin", "wb");

  fwrite(&dirlength, sizeof(int), 1, outfile);
  fwrite(phi_bssn, sizeof(double), arraylength, outfile);

  for(int coord=0; coord<3; coord++) {
    fwrite(Br[coord], sizeof(double), arraylength, outfile);
    fwrite(Bl[coord], sizeof(double), arraylength, outfile);

    fwrite(vrr[coord], sizeof(double), arraylength, outfile);
    fwrite(vrl[coord], sizeof(double), arraylength, outfile);
    fwrite(vlr[coord], sizeof(double), arraylength, outfile);
    fwrite(vll[coord], sizeof(double), arraylength, outfile);

    fwrite(cmin[coord], sizeof(double), arraylength, outfile);
    fwrite(cmax[coord], sizeof(double), arraylength, outfile);
  }
  fclose(outfile);

  for(int A_dir=1; A_dir<4; A_dir++) {
    int dir1 = A_dir%3, dir2 = (A_dir+1)%3;
    ghl_test_compute_A_flux_with_B(
              dirlength, A_dir, phi_bssn,
              cmin[dir1], cmax[dir1], cmin[dir2], cmax[dir2],
              vrr[dir1], vrl[dir1], vlr[dir1], vll[dir1],
              vrr[dir2], vrl[dir2], vlr[dir2], vll[dir2],
              Br[dir1], Bl[dir1], Br[dir2], Bl[dir2],
              A_rhs[A_dir-1]);
  }

  outfile = fopen_with_check("HLL_flux_with_B_output.bin","wb");
  fwrite(A_rhs[0], sizeof(double), arraylength, outfile);
  fwrite(A_rhs[1], sizeof(double), arraylength, outfile);
  fwrite(A_rhs[2], sizeof(double), arraylength, outfile);
  fclose(outfile);

  for(int A_dir=1; A_dir<4; A_dir++) {
    int dir1 = A_dir%3, dir2 = (A_dir+1)%3;
    ghl_test_compute_A_flux_with_Btilde(
              dirlength, A_dir,
              cmin[dir1], cmax[dir1], cmin[dir2], cmax[dir2],
              vrr[dir1], vrl[dir1], vlr[dir1], vll[dir1],
              vrr[dir2], vrl[dir2], vlr[dir2], vll[dir2],
              Br[dir1], Bl[dir1], Br[dir2], Bl[dir2],
              A_rhs[A_dir-1]);
  }

  outfile = fopen_with_check("HLL_flux_with_Btilde_output.bin","wb");
  fwrite(A_rhs[0], sizeof(double), arraylength, outfile);
  fwrite(A_rhs[1], sizeof(double), arraylength, outfile);
  fwrite(A_rhs[2], sizeof(double), arraylength, outfile);
  fclose(outfile);

  for(int index=0; index<arraylength; index++) {
    phi_bssn[index] *= (1.0 + randf(-1,1)*1.0e-14);

    for(int coord=0; coord<3; coord++) {
      Br[coord][index] *= (1.0 + randf(-1,1)*1.0e-14);
      Bl[coord][index] *= (1.0 + randf(-1,1)*1.0e-14);

      vrr[coord][index] *= (1.0 + randf(-1,1)*1.0e-14);
      vrl[coord][index] *= (1.0 + randf(-1,1)*1.0e-14);
      vlr[coord][index] *= (1.0 + randf(-1,1)*1.0e-14);
      vll[coord][index] *= (1.0 + randf(-1,1)*1.0e-14);

      cmin[coord][index] *= (1.0 + randf(-1,1)*1.0e-14);
      cmax[coord][index] *= (1.0 + randf(-1,1)*1.0e-14);
      cmin[coord][index] *= (1.0 + randf(-1,1)*1.0e-14);
      cmax[coord][index] *= (1.0 + randf(-1,1)*1.0e-14);
    }
  }

  for(int A_dir=1; A_dir<4; A_dir++) {
    int dir1 = A_dir%3, dir2 = (A_dir+1)%3;
    ghl_test_compute_A_flux_with_B(
              dirlength, A_dir, phi_bssn,
              cmin[dir1], cmax[dir1], cmin[dir2], cmax[dir2],
              vrr[dir1], vrl[dir1], vlr[dir1], vll[dir1],
              vrr[dir2], vrl[dir2], vlr[dir2], vll[dir2],
              Br[dir1], Bl[dir1], Br[dir2], Bl[dir2],
              A_rhs[A_dir-1]);
  }

  outfile = fopen_with_check("HLL_flux_with_B_output_pert.bin","wb");
  fwrite(A_rhs[0], sizeof(double), arraylength, outfile);
  fwrite(A_rhs[1], sizeof(double), arraylength, outfile);
  fwrite(A_rhs[2], sizeof(double), arraylength, outfile);
  fclose(outfile);

  for(int A_dir=1; A_dir<4; A_dir++) {
    int dir1 = A_dir%3, dir2 = (A_dir+1)%3;
    ghl_test_compute_A_flux_with_Btilde(
              dirlength, A_dir,
              cmin[dir1], cmax[dir1], cmin[dir2], cmax[dir2],
              vrr[dir1], vrl[dir1], vlr[dir1], vll[dir1],
              vrr[dir2], vrl[dir2], vlr[dir2], vll[dir2],
              Br[dir1], Bl[dir1], Br[dir2], Bl[dir2],
              A_rhs[A_dir-1]);
  }

  outfile = fopen_with_check("HLL_flux_with_Btilde_output_pert.bin","wb");
  fwrite(A_rhs[0], sizeof(double), arraylength, outfile);
  fwrite(A_rhs[1], sizeof(double), arraylength, outfile);
  fwrite(A_rhs[2], sizeof(double), arraylength, outfile);
  fclose(outfile);
}
