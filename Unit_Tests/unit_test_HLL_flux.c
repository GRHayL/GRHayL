#include "ghl_unit_tests.h"

int main(int argc, char **argv) {

  FILE* infile = fopen_with_check("HLL_flux_input.bin", "rb");

  int dirlength;
  int key = fread(&dirlength, sizeof(int), 1, infile);
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
  double *A_trusted[3];
  double *A_pert[3];
  for(int i=0; i<3; i++) {
    A_rhs[i] = (double*) malloc(sizeof(double)*arraylength);
    A_trusted[i] = (double*) malloc(sizeof(double)*arraylength);
    A_pert[i] = (double*) malloc(sizeof(double)*arraylength);
  }

  key = fread(phi_bssn, sizeof(double), arraylength, infile);

  for(int coord=0; coord<3; coord++) {
    key += fread(Br[coord], sizeof(double), arraylength, infile);
    key += fread(Bl[coord], sizeof(double), arraylength, infile);

    key += fread(vrr[coord], sizeof(double), arraylength, infile);
    key += fread(vrl[coord], sizeof(double), arraylength, infile);
    key += fread(vlr[coord], sizeof(double), arraylength, infile);
    key += fread(vll[coord], sizeof(double), arraylength, infile);

    key += fread(cmin[coord], sizeof(double), arraylength, infile);
    key += fread(cmax[coord], sizeof(double), arraylength, infile);
  }
  fclose(infile);
  if(key != arraylength*(8*3 + 1))
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  for(int k=1; k<dirlength; k++) {
    for(int j=1; j<dirlength; j++) {
      for(int i=1; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);
        A_rhs[0][index] = 0.0;
        A_rhs[1][index] = 0.0;
        A_rhs[2][index] = 0.0;
      }
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

  infile = fopen_with_check("HLL_flux_with_B_output.bin","rb");

  key = 0;
  for(int coord=0; coord<3; coord++)
    key += fread(A_trusted[coord], sizeof(double), arraylength, infile);
  fclose(infile);
  if(key != arraylength*3)
    ghl_error("An error has occured with reading in trusted data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  infile = fopen_with_check("HLL_flux_with_B_output_pert.bin","rb");
  key = 0;
  for(int coord=0; coord<3; coord++)
    key += fread(A_pert[coord], sizeof(double), arraylength, infile);
  fclose(infile);
  if(key != arraylength*3)
    ghl_error("An error has occured with reading in perturbed data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  for(int k=2; k<dirlength-2; k++) {
    for(int j=2; j<dirlength-2; j++) {
      for(int i=2; i<dirlength-2; i++) {
        const int index = indexf(dirlength,i,j,k);
        if( ghl_pert_test_fail(A_trusted[0][index], A_rhs[0][index], A_pert[0][index]) )
          ghl_error("Test unit_test_HLL_flux_with_B has failed for variable Ax_rhs.\n"
                       "  trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n",
                       A_trusted[0][index], A_rhs[0][index], A_pert[0][index],
                       relative_error(A_trusted[0][index], A_rhs[0][index]), relative_error(A_trusted[0][index], A_pert[0][index]));
        if( ghl_pert_test_fail(A_trusted[1][index], A_rhs[1][index], A_pert[1][index]) )
          ghl_error("Test unit_test_HLL_flux_with_B has failed for variable Ay_rhs.\n"
                       "  trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n",
                       A_trusted[1][index], A_rhs[1][index], A_pert[1][index],
                       relative_error(A_trusted[1][index], A_rhs[1][index]), relative_error(A_trusted[1][index], A_pert[1][index]));
        if( ghl_pert_test_fail(A_trusted[2][index], A_rhs[2][index], A_pert[2][index]) )
          ghl_error("Test unit_test_HLL_flux_with_B has failed for variable Az_rhs.\n"
                       "  trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n",
                       A_trusted[2][index], A_rhs[2][index], A_pert[2][index],
                       relative_error(A_trusted[2][index], A_rhs[2][index]), relative_error(A_trusted[2][index], A_pert[2][index]));
      }
    }
  }

  for(int k=1; k<dirlength; k++) {
    for(int j=1; j<dirlength; j++) {
      for(int i=1; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);
        A_rhs[0][index] = 0.0;
        A_rhs[1][index] = 0.0;
        A_rhs[2][index] = 0.0;
      }
    }
  }

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

  infile = fopen_with_check("HLL_flux_with_Btilde_output.bin","rb");

  key = 0;
  for(int coord=0; coord<3; coord++)
    key += fread(A_trusted[coord], sizeof(double), arraylength, infile);
  fclose(infile);
  if(key != arraylength*3)
    ghl_error("An error has occured with reading in trusted data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  infile = fopen_with_check("HLL_flux_with_Btilde_output_pert.bin","rb");
  key = 0;
  for(int coord=0; coord<3; coord++)
    key += fread(A_pert[coord], sizeof(double), arraylength, infile);
  fclose(infile);
  if(key != arraylength*3)
    ghl_error("An error has occured with reading in perturbed data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  for(int k=2; k<dirlength-2; k++) {
    for(int j=2; j<dirlength-2; j++) {
      for(int i=2; i<dirlength-2; i++) {
        const int index = indexf(dirlength,i,j,k);
        if( ghl_pert_test_fail(A_trusted[0][index], A_rhs[0][index], A_pert[0][index]) )
          ghl_error("Test unit_test_HLL_flux_with_Btilde has failed for variable Ax_rhs.\n"
                       "  trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n",
                       A_trusted[0][index], A_rhs[0][index], A_pert[0][index],
                       relative_error(A_trusted[0][index], A_rhs[0][index]), relative_error(A_trusted[0][index], A_pert[0][index]));
        if( ghl_pert_test_fail(A_trusted[1][index], A_rhs[1][index], A_pert[1][index]) )
          ghl_error("Test unit_test_HLL_flux_with_Btilde has failed for variable Ay_rhs.\n"
                       "  trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n",
                       A_trusted[1][index], A_rhs[1][index], A_pert[1][index],
                       relative_error(A_trusted[1][index], A_rhs[1][index]), relative_error(A_trusted[1][index], A_pert[1][index]));
        if( ghl_pert_test_fail(A_trusted[2][index], A_rhs[2][index], A_pert[2][index]) )
          ghl_error("Test unit_test_HLL_flux_with_Btilde has failed for variable Az_rhs.\n"
                       "  trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n",
                       A_trusted[2][index], A_rhs[2][index], A_pert[2][index],
                       relative_error(A_trusted[2][index], A_rhs[2][index]), relative_error(A_trusted[2][index], A_pert[2][index]));
      }
    }
  }

  ghl_info("HLL magnetic flux test has passed!\n");
  free(phi_bssn);

  for(int i=0; i<3; i++) {
    free(cmin[i]);
    free(cmax[i]);
    free(vrr[i]);
    free(vrl[i]);
    free(vlr[i]);
    free(vll[i]);
    free(Br[i]);
    free(Bl[i]);
    free(A_rhs[i]);
    free(A_trusted[i]);
    free(A_pert[i]);
  }
}
