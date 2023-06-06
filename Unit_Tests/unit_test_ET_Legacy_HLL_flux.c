#include "unit_tests.h"
#define IPH(METRICm1,METRICp0,METRICp1,METRICp2) (-0.0625*((METRICm1) + (METRICp2)) + 0.5625*((METRICp0) + (METRICp1)))

// each point requires phi from -1 to +2
// c{1,2}{min,max}, v{1,2}{r,l}{r,l} B{1,2}{r,l} from 0 to +1
void A_rhs_dir(const int dirlength, const int A_dir,
               const double *restrict phi_bssn,
               const double *restrict cmin_1, const double *restrict cmax_1,
               const double *restrict cmin_2, const double *restrict cmax_2,
               const double *restrict v1rr, const double *restrict v1rl,
               const double *restrict v1lr, const double *restrict v1ll,
               const double *restrict v2rr, const double *restrict v2rl,
               const double *restrict v2lr, const double *restrict v2ll,
               const double *restrict B1r, const double *restrict B1l,
               const double *restrict B2r, const double *restrict B2l,
               double *restrict A_rhs);

int main(int argc, char **argv) {

  FILE* infile = fopen("ET_Legacy_HLL_flux_input.bin", "rb");
  check_file_was_successfully_open(infile, "ET_Legacy_HLL_flux_input.bin");

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
  
  for(int k=1; k<dirlength; k++)
    for(int j=1; j<dirlength; j++)
      for(int i=1; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);
        A_rhs[0][index] = 0.0;
        A_rhs[1][index] = 0.0;
        A_rhs[2][index] = 0.0;
  }

  for(int A_dir=1; A_dir<4; A_dir++) {
  int dir1 = A_dir%3, dir2 = (A_dir+1)%3;
  A_rhs_dir(dirlength, A_dir, phi_bssn,
            cmin[dir1], cmax[dir1], cmin[dir2], cmax[dir2],
            vrr[dir1], vrl[dir1], vlr[dir1], vll[dir1],
            vrr[dir2], vrl[dir2], vlr[dir2], vll[dir2],
            Br[dir1], Bl[dir1], Br[dir2], Bl[dir2],
            A_rhs[A_dir-1]);
  }

  infile = fopen("ET_Legacy_HLL_flux_output.bin","rb");
  check_file_was_successfully_open(infile, "ET_Legacy_HLL_flux_output.bin");
  FILE *inpert = fopen("ET_Legacy_HLL_flux_output_pert.bin","rb");
  check_file_was_successfully_open(infile, "ET_Legacy_HLL_flux_output_pert.bin");

  key = 0;
  for(int coord=0; coord<3; coord++)
    key += fread(A_trusted[coord], sizeof(double), arraylength, infile);
  fclose(infile);
  if(key != arraylength*3)
    ghl_error("An error has occured with reading in trusted data. Please check that data\n"
                 "is up-to-date with current test version.\n");
  key = 0;
  for(int coord=0; coord<3; coord++)
    key += fread(A_pert[coord], sizeof(double), arraylength, inpert);
  fclose(inpert);
  if(key != arraylength*3)
    ghl_error("An error has occured with reading in perturbed data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  for(int k=2; k<dirlength-2; k++)
    for(int j=2; j<dirlength-2; j++)
      for(int i=2; i<dirlength-2; i++) {
        const int index = indexf(dirlength,i,j,k);
        if( validate(A_trusted[0][index], A_rhs[0][index], A_pert[0][index]) )
          ghl_error("Test unit_test_HLL_flux has failed for variable Ax_rhs.\n"
                       "  trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n",
                       A_trusted[0][index], A_rhs[0][index], A_pert[0][index],
                       relative_error(A_trusted[0][index], A_rhs[0][index]), relative_error(A_trusted[0][index], A_pert[0][index]));
        if( validate(A_trusted[1][index], A_rhs[1][index], A_pert[1][index]) )
          ghl_error("Test unit_test_HLL_flux has failed for variable Ay_rhs.\n"
                       "  trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n",
                       A_trusted[1][index], A_rhs[1][index], A_pert[1][index],
                       relative_error(A_trusted[1][index], A_rhs[1][index]), relative_error(A_trusted[1][index], A_pert[1][index]));
        if( validate(A_trusted[2][index], A_rhs[2][index], A_pert[2][index]) )
          ghl_error("Test unit_test_HLL_flux has failed for variable Az_rhs.\n"
                       "  trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n",
                       A_trusted[2][index], A_rhs[2][index], A_pert[2][index],
                       relative_error(A_trusted[2][index], A_rhs[2][index]), relative_error(A_trusted[2][index], A_pert[2][index]));
  }
  ghl_info("ET_Legacy HLL magnetic flux test has passed!\n");
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

void A_rhs_dir(const int dirlength,
               const int A_dir,
               const double *restrict phi_bssn,
               const double *restrict cmin_1,
               const double *restrict cmax_1,
               const double *restrict cmin_2,
               const double *restrict cmax_2,
               const double *restrict v1rr,
               const double *restrict v1rl,
               const double *restrict v1lr,
               const double *restrict v1ll,
               const double *restrict v2rr,
               const double *restrict v2rl,
               const double *restrict v2lr,
               const double *restrict v2ll,
               const double *restrict B1r,
               const double *restrict B1l,
               const double *restrict B2r,
               const double *restrict B2l,
               double *restrict A_rhs) {

  const int xdir = (A_dir==1);
  const int ydir = (A_dir==2);
  const int zdir = (A_dir==3);

  // This offsets the index by +1 in the perpendicular directions
  const int v_offset[3] = { !xdir, !ydir, !zdir };

  // This offsets the index by +1 in the permuted direction (x<-y<-z)
  const int B1_offset[3] = { ydir, zdir, xdir };

  // This offsets the index by +1 in the permuted direction (x->y->z)
  const int B2_offset[3] = { zdir, xdir, ydir };

  const int imax = dirlength-2;
  const int jmax = dirlength-2;
  const int kmax = dirlength-2;

  for(int k=2; k<kmax; k++)
    for(int j=2; j<jmax; j++)
      for(int i=2; i<imax; i++) {
        const int index    = indexf(dirlength,i,j,k);
        const int index_v  = indexf(dirlength,i+v_offset[0], j+v_offset[1], k+v_offset[2]);
        const int index_B1 = indexf(dirlength,i+B1_offset[0],j+B1_offset[1],k+B1_offset[2]);
        const int index_B2 = indexf(dirlength,i+B2_offset[0],j+B2_offset[1],k+B2_offset[2]);

        HLL_2D_vars vars;

        // This computes psi6 at the point staggered with respect to the two perpendicular
        // directions using the variable phi, which naturally lives at (i, j, k).
        // E.g. A_x needs phi at (i, j+1/2, k+1/2), so it must be interpolated to that point.
        // With the IPH macro, we first interpolate to the points
        // (i, j+1/2, k-1), (i, j+1/2, k), (i, j+1/2, k+1), (i, j+1/2, k+2) and use
        // those to compute phi at (i, j+1/2, k+1/2).
        const double psi6 =
          exp(6.0*IPH(
            IPH(phi_bssn[indexf(dirlength,i-!xdir  , j-xdir  -zdir,   k-!zdir)],
                phi_bssn[indexf(dirlength,i        , j       -zdir,   k-!zdir)],
                phi_bssn[indexf(dirlength,i+!xdir  , j+xdir  -zdir,   k-!zdir)],
                phi_bssn[indexf(dirlength,i+2*!xdir, j+2*xdir-zdir,   k-!zdir)]),

            IPH(phi_bssn[indexf(dirlength,i-!xdir,   j-xdir,          k      )],
                phi_bssn[index],
                phi_bssn[indexf(dirlength,i+!xdir,   j+xdir,          k      )],
                phi_bssn[indexf(dirlength,i+2*!xdir, j+2*xdir,        k      )]),

            IPH(phi_bssn[indexf(dirlength,i-!xdir,   j-xdir  +zdir,   k+!zdir)],
                phi_bssn[indexf(dirlength,i,         j       +zdir,   k+!zdir)],
                phi_bssn[indexf(dirlength,i+!xdir,   j+xdir  +zdir,   k+!zdir)],
                phi_bssn[indexf(dirlength,i+2*!xdir, j+2*xdir+zdir,   k+!zdir)]),

            IPH(phi_bssn[indexf(dirlength,i-!xdir,   j-xdir  +2*zdir, k+2*!zdir)],
                phi_bssn[indexf(dirlength,i,         j       +2*zdir, k+2*!zdir)],
                phi_bssn[indexf(dirlength,i+!xdir,   j+xdir  +2*zdir, k+2*!zdir)],
                phi_bssn[indexf(dirlength,i+2*!xdir, j+2*xdir+2*zdir, k+2*!zdir)])));

        vars.v1rr = v1rr[index_v];
        vars.v1rl = v1rl[index_v];
        vars.v1lr = v1lr[index_v];
        vars.v1ll = v1ll[index_v];

        vars.v2rr = v2rr[index_v];
        vars.v2rl = v2rl[index_v];
        vars.v2lr = v2lr[index_v];
        vars.v2ll = v2ll[index_v];

        vars.B1r = B1r[index_B1];
        vars.B1l = B1l[index_B1];

        vars.B2r = B2r[index_B2];
        vars.B2l = B2l[index_B2];

        vars.c1_min = cmin_1[index_B2];
        vars.c1_max = cmax_1[index_B2];
        vars.c2_min = cmin_2[index_B1];
        vars.c2_max = cmax_2[index_B1];

        A_rhs[index] = ghl_HLL_2D_flux(psi6, &vars);
  }
}
