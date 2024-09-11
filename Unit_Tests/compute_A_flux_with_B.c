#include "ghl_unit_tests.h"
#define IPH(METRICm1,METRICp0,METRICp1,METRICp2) (-0.0625*((METRICm1) + (METRICp2)) + 0.5625*((METRICp0) + (METRICp1)))

void ghl_test_compute_A_flux_with_B(
      const int dirlength,
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

  for(int k=2; k<kmax; k++) {
    for(int j=2; j<jmax; j++) {
      for(int i=2; i<imax; i++) {
        const int index    = indexf(dirlength,i,j,k);
        const int index_v  = indexf(dirlength,i+v_offset[0], j+v_offset[1], k+v_offset[2]);
        const int index_B1 = indexf(dirlength,i+B1_offset[0],j+B1_offset[1],k+B1_offset[2]);
        const int index_B2 = indexf(dirlength,i+B2_offset[0],j+B2_offset[1],k+B2_offset[2]);

        HLL_vars vars;

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

        A_rhs[index] = ghl_HLL_flux_with_B(psi6, &vars);
      }
    }
  }
}
