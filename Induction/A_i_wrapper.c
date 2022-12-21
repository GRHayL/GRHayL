/* Compute the part of A_i_rhs that excludes the gauge terms. I.e., we set
 *   A_i_rhs = \partial_t A_i = \psi^{6} (v^z B^x - v^x B^z)   here.
 */

#include "cctk.h"
#include "induction.h"

static const int BX_STAGGER=8, VXR=11, VXL=14;
typedef struct gf_and_gz_struct {
  double *gf;
  int gz_lo[4],gz_hi[4];
} gf_and_gz_struct;

void A_i_wrapper( const int A_dirn,
                  CCTK_POINTER_TO_CONST void_cctkGH,
                  CCTK_POINTER_TO_CONST void_prims_r,
                  CCTK_POINTER_TO_CONST void_prims_l,
                  CCTK_POINTER_TO_CONST void_phi_interped,
                  CCTK_POINTER_TO_CONST void_cmax_1,
                  CCTK_POINTER_TO_CONST void_cmin_1,
                  CCTK_POINTER_TO_CONST void_cmax_2,
                  CCTK_POINTER_TO_CONST void_cmin_2,
                  CCTK_POINTER void_A3_rhs ) {

  const cGH *cctkGH = (const cGH *)void_cctkGH;
  int *bounds = cctkGH->cctk_lsh;
  int *ghostzones = cctkGH->cctk_nghostzones;
  const gf_and_gz_struct *prims_r = (const gf_and_gz_struct *)void_prims_r;
  const gf_and_gz_struct *prims_l = (const gf_and_gz_struct *)void_prims_l;
  const CCTK_REAL *phi_interped = (const CCTK_REAL *)void_phi_interped;
  const CCTK_REAL *cmin_1 = (const CCTK_REAL *)void_cmin_1;
  const CCTK_REAL *cmax_1 = (const CCTK_REAL *)void_cmax_1;
  const CCTK_REAL *cmin_2 = (const CCTK_REAL *)void_cmin_2;
  const CCTK_REAL *cmax_2 = (const CCTK_REAL *)void_cmax_2;
  CCTK_REAL *A3_rhs = (CCTK_REAL *)void_A3_rhs;

  // If A_dirn=1, then v1_offset=1 (v1=VY) and v2_offset=2 (v2=VZ)
  // If A_dirn=2, then v1_offset=2 (v1=VZ) and v2_offset=0 (v2=VX)
  // If A_dirn=3, then v1_offset=0 (v1=VX) and v2_offset=1 (v2=VY)
  int v1_offset = ((A_dirn-1)+1)%3,           v2_offset = ((A_dirn-1)+2)%3;

  CCTK_REAL *v1rr=prims_r[VXR+v1_offset].gf, *v2rr=prims_r[VXR+v2_offset].gf;
  CCTK_REAL *v1rl=prims_l[VXR+v1_offset].gf, *v2rl=prims_l[VXR+v2_offset].gf;
  CCTK_REAL *v1lr=prims_r[VXL+v1_offset].gf, *v2lr=prims_r[VXL+v2_offset].gf;
  CCTK_REAL *v1ll=prims_l[VXL+v1_offset].gf, *v2ll=prims_l[VXL+v2_offset].gf;

  CCTK_REAL *B1r=prims_r[BX_STAGGER+v1_offset].gf, *B1l=prims_l[BX_STAGGER+v1_offset].gf;
  CCTK_REAL *B2r=prims_r[BX_STAGGER+v2_offset].gf, *B2l=prims_l[BX_STAGGER+v2_offset].gf;

  /**** V DEPENDENCIES ****/
  /* In the case of Ax_rhs, we need v{y,z}{r,l} at (i,j+1/2,k+1/2).
   *    However, v{y,z}{r,l}{r,l} are defined at (i,j-1/2,k-1/2), so
   *    v{y,z}{r,l} at (i,j+1/2,k+1/2) is stored at v{y,z}{r,l}{r,l}(i,j+1,k+1).
   * In the case of Ay_rhs, we need v{x,z}{r,l} at (i+1/2,j,k+1/2).
   *    However, v{x,z}{r,l}{r,l} are defined at (i-1/2,j,k-1/2), so
   *    v{x,z}{r,l} at (i+1/2,j,k+1/2) is stored at v{x,z}{r,l}{r,l}(i+1,j,k+1).
   * In the case of Az_rhs, we need v{x,y}{r,l} at (i+1/2,j+1/2,k).
   *    However, v{x,y}{r,l}{r,l} are defined at (i-1/2,j-1/2,k), so
   *    v{x,y}{r,l} at (i+1/2,j+1/2,k) is stored at v{x,y}{r,l}{r,l}(i+1,j+1,k). */
  static const int vs_ijk_offset[4][3] = { {0,0,0} , {0,1,1} , {1,0,1} , {1,1,0} };

  /**** B DEPENDENCIES ****/
  /* In the case of Ax_rhs, we need B{y,z}{r,l} at (i,j+1/2,k+1/2).
   *    However, By_stagger{r,l} is defined at (i,j+1/2,k-1/2), and
   *             Bz_stagger{r,l} is defined at (i,j-1/2,k+1/2), so
   *             By_stagger{r,l} at (i,j+1/2,k+1/2) is stored at By_stagger{r,l}(i,j,k+1), and
   *             Bz_stagger{r,l} at (i,j+1/2,k+1/2) is stored at Bz_stagger{r,l}(i,j+1,k).
   * In the case of Ay_rhs, we need B{z,x}_stagger{r,l} at (i+1/2,j,k+1/2).
   *    However, Bz_stagger{r,l} is defined at (i-1/2,j,k+1/2), and
   *             Bx_stagger{r,l} is defined at (i+1/2,j,k-1/2), so
   *             Bz_stagger{r,l} at (i+1/2,j,k+1/2) is stored at Bz_stagger{r,l}(i+1,j,k), and
   *             Bx_stagger{r,l} at (i+1/2,j,k+1/2) is stored at Bx_stagger{r,l}(i,j,k+1).
   * In the case of Az_rhs, we need B{x,y}_stagger{r,l} at (i+1/2,j+1/2,k).
   *    However, Bx_stagger{r,l} is defined at (i+1/2,j-1/2,k), and
   *             By_stagger{r,l} is defined at (i-1/2,j+1/2,k), so
   *             Bx_stagger{r,l} at (i+1/2,j+1/2,k) is stored at Bx_stagger{r,l}(i,j+1,k), and
   *             By_stagger{r,l} at (i+1/2,j+1/2,k) is stored at By_stagger{r,l}(i+1,j,k).
   */
  static const int B1_ijk_offset[4][3] = { {0,0,0} , {0,0,1} , {1,0,0} , {0,1,0} };
  static const int B2_ijk_offset[4][3] = { {0,0,0} , {0,1,0} , {0,0,1} , {1,0,0} };

#pragma omp parallel for
  for(int k=ghostzones[2];k<bounds[2]-ghostzones[2];k++)
    for(int j=ghostzones[1];j<bounds[1]-ghostzones[1];j++)
      for(int i=ghostzones[0];i<bounds[0]-ghostzones[0];i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
        // The following lines set the indices appropriately. See justification in exorbitant comments above.
        int index_v =CCTK_GFINDEX3D(cctkGH,i+vs_ijk_offset[A_dirn][0],j+vs_ijk_offset[A_dirn][1],k+vs_ijk_offset[A_dirn][2]);
        int index_B1=CCTK_GFINDEX3D(cctkGH,i+B1_ijk_offset[A_dirn][0],j+B1_ijk_offset[A_dirn][1],k+B1_ijk_offset[A_dirn][2]);
        int index_B2=CCTK_GFINDEX3D(cctkGH,i+B2_ijk_offset[A_dirn][0],j+B2_ijk_offset[A_dirn][1],k+B2_ijk_offset[A_dirn][2]);

        induction_lr vars;

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

        // Stores 1/sqrt(gamma)==exp(6 phi) at (i+1/2,j+1/2,k) for Az, (i+1/2,j,k+1/2) for Ay, and (i,j+1/2,k+1/2) for Az.
        CCTK_REAL psi6_interped=exp(6.0*(phi_interped[index]));

        A3_rhs[index] = A_i_rhs_no_gauge_terms(psi6_interped, &vars);
      }
}

