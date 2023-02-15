#include "cctk.h"
#include "GRHayLET.h"
#define IPH(METRICm1,METRICp0,METRICp1,METRICp2) (-0.0625*((METRICm1) + (METRICp2)) + 0.5625*((METRICp0) + (METRICp1)))

void A_no_gauge_rhs(const cGH *restrict cctkGH, const int A_dir,
               /*const*/ double **out_prims_r,
               /*const*/ double **out_prims_l,
               const double *restrict phi_bssn,
               const double *restrict cmin_1,
               const double *restrict cmax_1,
               const double *restrict cmin_2,
               const double *restrict cmax_2,
               double *restrict A_rhs) {
  const int xdir = (A_dir==1);
  const int ydir = (A_dir==2);
  const int zdir = (A_dir==3);

  // These are used to determine which components of v and
  // B are used in the computation.
  const int v1_offset = A_dir%3, v2_offset = (A_dir+1)%3;

  // This offsets the index by +1 in the perpendicular directions
  const int v_offset[3] = { !xdir, !ydir, !zdir };

  // This offsets the index by +1 in the permuted direction (x<-y<-z)
  const int B1_offset[3] = { ydir, zdir, xdir };

  // This offsets the index by +1 in the permuted direction (x->y->z)
  const int B2_offset[3] = { zdir, xdir, ydir };

  const int imin = cctkGH->cctk_nghostzones[0];
  const int imax = cctkGH->cctk_lsh[0]-cctkGH->cctk_nghostzones[0];
  const int jmin = cctkGH->cctk_nghostzones[1];
  const int jmax = cctkGH->cctk_lsh[1]-cctkGH->cctk_nghostzones[1];
  const int kmin = cctkGH->cctk_nghostzones[2];
  const int kmax = cctkGH->cctk_lsh[2]-cctkGH->cctk_nghostzones[2];

#pragma omp parallel for
  for(int k=kmin; k<kmax; k++)
    for(int j=jmin; j<jmax; j++)
      for(int i=imin; i<imax; i++) {
        const int index    = CCTK_GFINDEX3D(cctkGH,i,j,k);
        const int index_v  = CCTK_GFINDEX3D(cctkGH,i+v_offset[0], j+v_offset[1], k+v_offset[2]);
        const int index_B1 = CCTK_GFINDEX3D(cctkGH,i+B1_offset[0],j+B1_offset[1],k+B1_offset[2]);
        const int index_B2 = CCTK_GFINDEX3D(cctkGH,i+B2_offset[0],j+B2_offset[1],k+B2_offset[2]);

        A_no_gauge_vars vars;

        // This computes psi6 at the point staggered with respect to the two perpendicular
        // directions using the variable phi, which naturally lives at (i, j, k).
        // E.g. A_x needs phi at (i, j+1/2, k+1/2), so it must be interpolated to that point.
        // With the IPH macro, we first interpolate to the points
        // (i, j+1/2, k-1), (i, j+1/2, k), (i, j+1/2, k+1), (i, j+1/2, k+2) and use
        // those to compute phi at (i, j+1/2, k+1/2).
        vars.psi6 =
          exp(6.0*IPH(
            IPH(phi_bssn[CCTK_GFINDEX3D(cctkGH,i-!xdir  , j-xdir  -zdir,   k-!zdir)],
                phi_bssn[CCTK_GFINDEX3D(cctkGH,i        , j       -zdir,   k-!zdir)],
                phi_bssn[CCTK_GFINDEX3D(cctkGH,i+!xdir  , j+xdir  -zdir,   k-!zdir)],
                phi_bssn[CCTK_GFINDEX3D(cctkGH,i+2*!xdir, j+2*xdir-zdir,   k-!zdir)]),

            IPH(phi_bssn[CCTK_GFINDEX3D(cctkGH,i-!xdir,   j-xdir,          k      )],
                phi_bssn[index],
                phi_bssn[CCTK_GFINDEX3D(cctkGH,i+!xdir,   j+xdir,          k      )],
                phi_bssn[CCTK_GFINDEX3D(cctkGH,i+2*!xdir, j+2*xdir,        k      )]),

            IPH(phi_bssn[CCTK_GFINDEX3D(cctkGH,i-!xdir,   j-xdir  +zdir,   k+!zdir)],
                phi_bssn[CCTK_GFINDEX3D(cctkGH,i,         j       +zdir,   k+!zdir)],
                phi_bssn[CCTK_GFINDEX3D(cctkGH,i+!xdir,   j+xdir  +zdir,   k+!zdir)],
                phi_bssn[CCTK_GFINDEX3D(cctkGH,i+2*!xdir, j+2*xdir+zdir,   k+!zdir)]),

            IPH(phi_bssn[CCTK_GFINDEX3D(cctkGH,i-!xdir,   j-xdir  +2*zdir, k+2*!zdir)],
                phi_bssn[CCTK_GFINDEX3D(cctkGH,i,         j       +2*zdir, k+2*!zdir)],
                phi_bssn[CCTK_GFINDEX3D(cctkGH,i+!xdir,   j+xdir  +2*zdir, k+2*!zdir)],
                phi_bssn[CCTK_GFINDEX3D(cctkGH,i+2*!xdir, j+2*xdir+2*zdir, k+2*!zdir)])));

        vars.v1rr=out_prims_r[VXR+v1_offset][index_v];
        vars.v1rl=out_prims_l[VXR+v1_offset][index_v];
        vars.v1lr=out_prims_r[VXL+v1_offset][index_v];
        vars.v1ll=out_prims_l[VXL+v1_offset][index_v];

        vars.v2rr=out_prims_r[VXR+v2_offset][index_v];
        vars.v2rl=out_prims_l[VXR+v2_offset][index_v];
        vars.v2lr=out_prims_r[VXL+v2_offset][index_v];
        vars.v2ll=out_prims_l[VXL+v2_offset][index_v];

        vars.B1r=out_prims_r[BX_STAGGER+v1_offset][index_B1];
        vars.B1l=out_prims_l[BX_STAGGER+v1_offset][index_B1];

        vars.B2r=out_prims_r[BX_STAGGER+v2_offset][index_B2];
        vars.B2l=out_prims_l[BX_STAGGER+v2_offset][index_B2];

        vars.c1_min = cmin_1[index_B2];
        vars.c1_max = cmax_1[index_B2];
        vars.c2_min = cmin_2[index_B1];
        vars.c2_max = cmax_2[index_B1];

        A_rhs[index] = HLL_flux(&vars);
  }
}
