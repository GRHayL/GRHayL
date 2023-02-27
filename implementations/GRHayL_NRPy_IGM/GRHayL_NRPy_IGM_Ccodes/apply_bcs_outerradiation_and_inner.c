#include "././NRPy_basic_defines.h"
#include "././NRPy_function_prototypes.h"
/*
 * Compute 1st derivative finite-difference derivative with arbitrary upwind
 */
static inline REAL FD1_arbitrary_upwind_x0_dirn(const paramstruct *restrict params, const REAL *restrict gf,
                                                const int i0,const int i1,const int i2, const int offset) {
#include "./set_Cparameters.h"

  switch(offset) {
  case 0:
    return (-1.0/2.0*gf[IDX3S(i0-1,i1,i2)]
            +1.0/2.0*gf[IDX3S(i0+1,i1,i2)]) * invdx0;
  case 1:
    return (-3.0/2.0*gf[IDX3S(i0,i1,i2)]
            +2*gf[IDX3S(i0+1,i1,i2)]
            -1.0/2.0*gf[IDX3S(i0+2,i1,i2)]) * invdx0;
  case -1:
    return (+1.0/2.0*gf[IDX3S(i0-2,i1,i2)]
            -2*gf[IDX3S(i0-1,i1,i2)]
            +3.0/2.0*gf[IDX3S(i0,i1,i2)]) * invdx0;
  }
  return 0.0 / 0.0;  // poison output if offset computed incorrectly
}
/*
 * Compute 1st derivative finite-difference derivative with arbitrary upwind
 */
static inline REAL FD1_arbitrary_upwind_x1_dirn(const paramstruct *restrict params, const REAL *restrict gf,
                                                const int i0,const int i1,const int i2, const int offset) {
#include "./set_Cparameters.h"

  switch(offset) {
  case 0:
    return (-1.0/2.0*gf[IDX3S(i0,i1-1,i2)]
            +1.0/2.0*gf[IDX3S(i0,i1+1,i2)]) * invdx1;
  case 1:
    return (-3.0/2.0*gf[IDX3S(i0,i1,i2)]
            +2*gf[IDX3S(i0,i1+1,i2)]
            -1.0/2.0*gf[IDX3S(i0,i1+2,i2)]) * invdx1;
  case -1:
    return (+1.0/2.0*gf[IDX3S(i0,i1-2,i2)]
            -2*gf[IDX3S(i0,i1-1,i2)]
            +3.0/2.0*gf[IDX3S(i0,i1,i2)]) * invdx1;
  }
  return 0.0 / 0.0;  // poison output if offset computed incorrectly
}
/*
 * Compute 1st derivative finite-difference derivative with arbitrary upwind
 */
static inline REAL FD1_arbitrary_upwind_x2_dirn(const paramstruct *restrict params, const REAL *restrict gf,
                                                const int i0,const int i1,const int i2, const int offset) {
#include "./set_Cparameters.h"

  switch(offset) {
  case 0:
    return (-1.0/2.0*gf[IDX3S(i0,i1,i2-1)]
            +1.0/2.0*gf[IDX3S(i0,i1,i2+1)]) * invdx2;
  case 1:
    return (-3.0/2.0*gf[IDX3S(i0,i1,i2)]
            +2*gf[IDX3S(i0,i1,i2+1)]
            -1.0/2.0*gf[IDX3S(i0,i1,i2+2)]) * invdx2;
  case -1:
    return (+1.0/2.0*gf[IDX3S(i0,i1,i2-2)]
            -2*gf[IDX3S(i0,i1,i2-1)]
            +3.0/2.0*gf[IDX3S(i0,i1,i2)]) * invdx2;
  }
  return 0.0 / 0.0;  // poison output if offset computed incorrectly
}
/*
 * Compute r(xx0,xx1,xx2) and partial_r x^i.
 */
static inline void r_and_partial_xi_partial_r_derivs(const paramstruct *restrict params,const REAL xx0,const REAL xx1,const REAL xx2,
                                                     REAL *r,
                                                     REAL *partial_x0_partial_r,REAL *partial_x1_partial_r,REAL *partial_x2_partial_r) {
#include "./set_Cparameters.h"

  const double tmp_0 = ((xx2)*(xx2));
  const double tmp_3 = ((xx0)*(xx0)) + ((xx1)*(xx1));
  const double tmp_4 = tmp_0 + tmp_3;
  const double tmp_5 = sqrt(tmp_4);
  const double tmp_6 = 1/(tmp_3*sqrt(-tmp_0/tmp_4 + 1));
  const double tmp_7 = tmp_6*((xx0)*(xx0));
  const double tmp_8 = tmp_0/((tmp_4)*(tmp_4));
  const double tmp_9 = tmp_6*((xx1)*(xx1));
  const double tmp_10 = (1.0/(tmp_5));
  const double tmp_11 = pow(tmp_4, -3.0/2.0);
  const double tmp_12 = -tmp_0*tmp_11 + tmp_10;
  const double tmp_14 = (1.0/(tmp_10*tmp_12*tmp_7 + tmp_10*tmp_12*tmp_9 + tmp_7*tmp_8 + tmp_8*tmp_9));
  const double tmp_15 = tmp_12*tmp_14*tmp_6;
  *r = tmp_5;
  *partial_x0_partial_r = tmp_15*xx0;
  *partial_x1_partial_r = tmp_15*xx1;
  *partial_x2_partial_r = tmp_14*(tmp_11*tmp_7*xx2 + tmp_11*tmp_9*xx2);
}
/*
 * Compute \partial_r f
 */
static inline REAL compute_partial_r_f(const paramstruct *restrict params, REAL *restrict xx[3], const REAL *restrict gfs,
                                       const int which_gf, const int dest_i0,const int dest_i1,const int dest_i2,
                                       const int FACEi0,const int FACEi1,const int FACEi2,
                                       const REAL partial_x0_partial_r, const REAL partial_x1_partial_r, const REAL partial_x2_partial_r) {
#include "./set_Cparameters.h"

  ///////////////////////////////////////////////////////////

  // FD1_stencil_radius = radiation_BC_FD_order/2 = 1
  const int FD1_stencil_radius = 1;

  const int ntot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

  ///////////////////////////////////////////////////////////
  // Next we'll compute partial_xi f, using a maximally-centered stencil.
  //   The {i0,i1,i2}_offset parameters set the offset of the maximally-centered
  //   stencil, such that an offset=0 implies a centered stencil.

  // CHECK: Nxx_plus_2NGHOSTS0=10; FD1_stencil_radius=2. Then Nxx_plus_2NGHOSTS0-FD1_stencil_radius-1 = 7
  //  if dest_i0 = 9, we get i0_offset=7-9=-2, so the (4th order) deriv
  //  stencil is: -4,-3,-2,-1,0

  // CHECK: if FD1_stencil_radius=2 and dest_i0 = 1, we get i0_offset = FD1_stencil_radius-dest_i0 = 1,
  //  so the (4th order) deriv stencil is: -1,0,1,2,3

  // CHECK: if FD1_stencil_radius=2 and dest_i0 = 0, we get i0_offset = FD1_stencil_radius-1 = 2,
  //  so the (4th order) deriv stencil is: 0,1,2,3,4
  int i0_offset = FACEi0;  // Shift stencil away from the face we're updating.
  // Next adjust i0_offset so that FD stencil never goes out of bounds.
  if(dest_i0 < FD1_stencil_radius) i0_offset = FD1_stencil_radius-dest_i0;
  else if(dest_i0 > (Nxx_plus_2NGHOSTS0-FD1_stencil_radius-1)) i0_offset = (Nxx_plus_2NGHOSTS0-FD1_stencil_radius-1) - dest_i0;
  const REAL partial_x0_f=FD1_arbitrary_upwind_x0_dirn(params,&gfs[which_gf*ntot],dest_i0,dest_i1,dest_i2,i0_offset);
  int i1_offset = FACEi1;  // Shift stencil away from the face we're updating.
  // Next adjust i1_offset so that FD stencil never goes out of bounds.
  if(dest_i1 < FD1_stencil_radius) i1_offset = FD1_stencil_radius-dest_i1;
  else if(dest_i1 > (Nxx_plus_2NGHOSTS1-FD1_stencil_radius-1)) i1_offset = (Nxx_plus_2NGHOSTS1-FD1_stencil_radius-1) - dest_i1;
  const REAL partial_x1_f=FD1_arbitrary_upwind_x1_dirn(params,&gfs[which_gf*ntot],dest_i0,dest_i1,dest_i2,i1_offset);
  int i2_offset = FACEi2;  // Shift stencil away from the face we're updating.
  // Next adjust i2_offset so that FD stencil never goes out of bounds.
  if(dest_i2 < FD1_stencil_radius) i2_offset = FD1_stencil_radius-dest_i2;
  else if(dest_i2 > (Nxx_plus_2NGHOSTS2-FD1_stencil_radius-1)) i2_offset = (Nxx_plus_2NGHOSTS2-FD1_stencil_radius-1) - dest_i2;
  const REAL partial_x2_f=FD1_arbitrary_upwind_x2_dirn(params,&gfs[which_gf*ntot],dest_i0,dest_i1,dest_i2,i2_offset);
  return partial_x0_partial_r*partial_x0_f + partial_x1_partial_r*partial_x1_f + partial_x2_partial_r*partial_x2_f;
}

/*
 * *** Apply radiation BCs to all outer boundaries. ***
 */
static inline REAL radiation_bcs(const paramstruct *restrict params, const bc_struct *restrict bcstruct,REAL *restrict xx[3],
                                 const REAL *restrict gfs, REAL *restrict gfs_rhss,
                                 const int which_gf, const REAL gf_wavespeed, const REAL gf_f_infinity,
                                 const int dest_i0,const int dest_i1,const int dest_i2,
                                 const short FACEi0,const short FACEi1,const short FACEi2) {
#include "./set_Cparameters.h"

  // Nearest "interior" neighbor of this gridpoint, based on current face
  const int dest_i0_int=dest_i0+1*FACEi0, dest_i1_int=dest_i1+1*FACEi1, dest_i2_int=dest_i2+1*FACEi2;
  REAL r, partial_x0_partial_r,partial_x1_partial_r,partial_x2_partial_r;
  REAL r_int, partial_x0_partial_r_int,partial_x1_partial_r_int,partial_x2_partial_r_int;
  r_and_partial_xi_partial_r_derivs(params,xx[0][dest_i0],xx[1][dest_i1],xx[2][dest_i2],
                                    &r, &partial_x0_partial_r, &partial_x1_partial_r,  &partial_x2_partial_r);
  r_and_partial_xi_partial_r_derivs(params, xx[0][dest_i0_int], xx[1][dest_i1_int], xx[2][dest_i2_int],
                                    &r_int, &partial_x0_partial_r_int, &partial_x1_partial_r_int, &partial_x2_partial_r_int);
  const REAL partial_r_f     = compute_partial_r_f(params,xx,gfs, which_gf,dest_i0,    dest_i1,    dest_i2,
                                                   FACEi0,FACEi1,FACEi2,
                                                   partial_x0_partial_r    ,partial_x1_partial_r    ,partial_x2_partial_r);
  const REAL partial_r_f_int = compute_partial_r_f(params,xx,gfs, which_gf,dest_i0_int,dest_i1_int,dest_i2_int,
                                                   FACEi0,FACEi1,FACEi2,
                                                   partial_x0_partial_r_int,partial_x1_partial_r_int,partial_x2_partial_r_int);

  const int idx3 = IDX3S(dest_i0,dest_i1,dest_i2);
  const int idx3_int = IDX3S(dest_i0_int,dest_i1_int,dest_i2_int);

  const REAL partial_t_f_int = gfs_rhss[IDX4ptS(which_gf, idx3_int)];

  const REAL c = gf_wavespeed;
  const REAL f_infinity = gf_f_infinity;
  const REAL f     = gfs[IDX4ptS(which_gf, idx3)];
  const REAL f_int = gfs[IDX4ptS(which_gf, idx3_int)];
  const REAL partial_t_f_int_outgoing_wave = -c * (partial_r_f_int + (f_int - f_infinity) / r_int);

  const REAL k = r_int*r_int*r_int * (partial_t_f_int - partial_t_f_int_outgoing_wave);

  const REAL rinv = 1.0 / r;
  const REAL partial_t_f_outgoing_wave = -c * (partial_r_f + (f - f_infinity) * rinv);

  return partial_t_f_outgoing_wave + k * rinv*rinv*rinv;
}

void apply_bcs_outerradiation_and_inner(const paramstruct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict xx[3],
                                        const REAL custom_wavespeed[NUM_EVOL_GFS],
                                        const REAL custom_f_infinity[NUM_EVOL_GFS],
                                        REAL *restrict gfs, REAL *restrict rhs_gfs) {
#include "./set_Cparameters.h"


  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;

  ////////////////////////////////////////////////////////
  // STEP 1 of 2: Apply BCs to pure outer boundary points.
  //              By "pure" we mean that these points are
  //              on the outer boundary and not also on
  //              an inner boundary.
  //              Here we fill in the innermost ghost zone
  //              layer first and move outward. At each
  //              layer, we fill in +/- x0 faces first,
  //              then +/- x1 faces, finally +/- x2 faces,
  //              filling in the edges as we go.
  // Spawn N OpenMP threads, either across all cores, or according to e.g., taskset.
#pragma omp parallel
  {
    for(int which_gz=0;which_gz<NGHOSTS;which_gz++) for(int dirn=0;dirn<3;dirn++) {
        // This option results in about 1.6% slower runtime for SW curvilinear at 64x24x24 on 8-core Ryzen 9 4900HS
        //#pragma omp for collapse(2)
        //for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) for(int idx2d=0;idx2d<bc_info->num_pure_outer_boundary_points[which_gz][dirn];idx2d++) {
        //  {
        // Don't spawn a thread if there are no boundary points to fill; results in a nice little speedup.
        if(bc_info->num_pure_outer_boundary_points[which_gz][dirn] > 0) {
#pragma omp for  // threads have been spawned; here we distribute across them
          for(int idx2d=0;idx2d<bc_info->num_pure_outer_boundary_points[which_gz][dirn];idx2d++) {
            const short i0 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i0;
            const short i1 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i1;
            const short i2 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i2;
            const short FACEX0 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX0;
            const short FACEX1 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX1;
            const short FACEX2 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX2;
            const int idx3 = IDX3S(i0,i1,i2);
            for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
              // *** Apply radiation BCs to all outer boundary points. ***
              rhs_gfs[IDX4ptS(which_gf, idx3)] = radiation_bcs(params, bcstruct, xx, gfs, rhs_gfs, which_gf,
                                                               custom_wavespeed[which_gf], custom_f_infinity[which_gf],
                                                               i0,i1,i2, FACEX0,FACEX1,FACEX2);
            }
          }
        }
      }
  }

  ///////////////////////////////////////////////////////
  // STEP 2 of 2: Apply BCs to inner boundary points.
  //              These map to either the grid interior
  //              ("pure inner") or to pure outer boundary
  //              points ("inner maps to outer"). Those
  //              that map to outer require that outer be
  //              populated first; hence this being
  //              STEP 2 OF 2.
  apply_bcs_inner_only(params, bcstruct, rhs_gfs); // <- apply inner BCs to RHS gfs only
}
