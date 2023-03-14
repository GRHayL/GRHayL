#include "cctk.h"
#include "IGM.h"
#define IPH(METRICm1,METRICp0,METRICp1,METRICp2) (-0.0625*((METRICm1) + (METRICp2)) + 0.5625*((METRICp0) + (METRICp1)))

void A_no_gauge_rhs(const cGH *restrict cctkGH, const int A_dir,
               const double **metric,
               /*const*/ double **in_prims_r,
               /*const*/ double **in_prims_l,
               const double *restrict phi_bssn,
//               const double **cmin,
//               const double **cmax,
               double *restrict A_rhs) {

  const double poison = 0.0/0.0;
  const int xdir = (A_dir==0);
  const int ydir = (A_dir==1);
  const int zdir = (A_dir==2);

  void (*calculate_characteristic_speed)(const primitive_quantities *restrict prims_r,
                                         const primitive_quantities *restrict prims_l,
                                         struct eos_parameters const *restrict eos,
                                         const metric_quantities *restrict metric_face,
                                         double *cmin, double *cmax);

  // Set function pointer to specific function for a given direction
  switch(A_dir) {
    case 0:
      calculate_characteristic_speed = &calculate_characteristic_speed_dirn0;
      break;
    case 1:
      calculate_characteristic_speed = &calculate_characteristic_speed_dirn1;
      break;
    case 2:
      calculate_characteristic_speed = &calculate_characteristic_speed_dirn2;
      break;
    default:
      CCTK_VERROR("Warning: invalid A_dir value (not 0, 1, or 2) has been passed to A_no_gauge_rhs.");
  }

  // These are used to determine which components of v and
  // B are used in the computation.
  const int v1_offset = (A_dir+1)%3, v2_offset = (A_dir-1)%3;

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

        vars.v1rr=in_prims_r[VXR+v1_offset][index_v];
        vars.v1rl=in_prims_l[VXR+v1_offset][index_v];
        vars.v1lr=in_prims_r[VXL+v1_offset][index_v];
        vars.v1ll=in_prims_l[VXL+v1_offset][index_v];

        vars.v2rr=in_prims_r[VXR+v2_offset][index_v];
        vars.v2rl=in_prims_l[VXR+v2_offset][index_v];
        vars.v2lr=in_prims_r[VXL+v2_offset][index_v];
        vars.v2ll=in_prims_l[VXL+v2_offset][index_v];

        vars.B1r=in_prims_r[BX_STAGGER+v1_offset][index_B1];
        vars.B1l=in_prims_l[BX_STAGGER+v1_offset][index_B1];

        vars.B2r=in_prims_r[BX_STAGGER+v2_offset][index_B2];
        vars.B2l=in_prims_l[BX_STAGGER+v2_offset][index_B2];

        //vars.c1_min = cmin[v1_offset][index_B2];
        //vars.c1_max = cmax[v1_offset][index_B2];
        //vars.c2_min = cmin[v2_offset][index_B1];
        //vars.c2_max = cmax[v2_offset][index_B1];

//TODO: we should probably get cmin/max out of the flux function so we don't recompute all of this
        metric_quantities metric_face;
        primitive_quantities prims_r, prims_l;
        interpolate_to_face_and_initialize_metric(
                          cctkGH, i+B2_offset[0], j+B2_offset[1], k+B2_offset[2],
                          A_dir, metric[LAPSE],
                          metric[BETAX], metric[BETAY], metric[BETAZ],
                          metric[GXX], metric[GXY], metric[GXZ],
                          metric[GYY], metric[GYZ], metric[GZZ],
                          &metric_face);

        initialize_primitives(in_prims_r[RHOB][index_B2], in_prims_r[PRESSURE][index_B2], poison,
                              in_prims_r[VX][index_B2], in_prims_r[VY][index_B2], in_prims_r[VZ][index_B2],
                              in_prims_r[BX_CENTER][index_B2], in_prims_r[BY_CENTER][index_B2], in_prims_r[BZ_CENTER][index_B2],
                              poison, poison, poison, // entropy, Y_e, temp
                              &prims_r);

        initialize_primitives(in_prims_l[RHOB][index_B2], in_prims_l[PRESSURE][index_B2], poison,
                              in_prims_l[VX][index_B2], in_prims_l[VY][index_B2], in_prims_l[VZ][index_B2],
                              in_prims_l[BX_CENTER][index_B2], in_prims_l[BY_CENTER][index_B2], in_prims_l[BZ_CENTER][index_B2],
                              poison, poison, poison, // entropy, Y_e, temp
                              &prims_l);


        calculate_characteristic_speed(&prims_r, &prims_l, grhayl_eos, &metric_face, &vars.c1_min, &vars.c1_max);

        interpolate_to_face_and_initialize_metric(
                          cctkGH, i+B1_offset[0], j+B1_offset[1], k+B1_offset[2],
                          A_dir, metric[LAPSE],
                          metric[BETAX], metric[BETAY], metric[BETAZ],
                          metric[GXX], metric[GXY], metric[GXZ],
                          metric[GYY], metric[GYZ], metric[GZZ],
                          &metric_face);

        initialize_primitives(in_prims_r[RHOB][index_B1], in_prims_r[PRESSURE][index_B1], poison,
                              in_prims_r[VX][index_B1], in_prims_r[VY][index_B1], in_prims_r[VZ][index_B1],
                              in_prims_r[BX_CENTER][index_B1], in_prims_r[BY_CENTER][index_B1], in_prims_r[BZ_CENTER][index_B1],
                              poison, poison, poison, // entropy, Y_e, temp
                              &prims_r);

        initialize_primitives(in_prims_l[RHOB][index_B1], in_prims_l[PRESSURE][index_B1], poison,
                              in_prims_l[VX][index_B1], in_prims_l[VY][index_B1], in_prims_l[VZ][index_B1],
                              in_prims_l[BX_CENTER][index_B1], in_prims_l[BY_CENTER][index_B1], in_prims_l[BZ_CENTER][index_B1],
                              poison, poison, poison, // entropy, Y_e, temp
                              &prims_l);

        calculate_characteristic_speed(&prims_r, &prims_l, grhayl_eos, &metric_face, &vars.c2_min, &vars.c2_max);

        A_rhs[index] = HLL_flux(&vars);
  }
}
