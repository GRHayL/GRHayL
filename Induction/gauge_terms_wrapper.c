#include "cctk.h"
#include "induction_gem.h"

static const int SHIFTXI=0,SHIFTYI=1,SHIFTZI=2,GUPXXI=3,GUPXYI=4,GUPXZI=5,GUPYYI=6,GUPYZI=7,GUPZZI=8,
  PSII=9,LAPM1I=10,A_XI=11,A_YI=12,A_ZI=13;

void A_i_gauge_rhs(CCTK_POINTER_TO_CONST void_cctkGH,
                   CCTK_POINTER_TO_CONST void_dX,
                   CCTK_POINTER_TO_CONST void_in_vars,
                   CCTK_POINTER_TO_CONST void_phitilde,
                   CCTK_REAL Lorenz_damping_factor,
                   CCTK_REAL *restrict shiftx_interp,
                   CCTK_REAL *restrict shifty_interp,
                   CCTK_REAL *restrict shiftz_interp,
                   CCTK_REAL *restrict alpha_interp,
                   CCTK_REAL *restrict alpha_Phi_minus_betaj_A_j_interp,
                   CCTK_REAL *restrict alpha_sqrtg_Ax_interp,
                   CCTK_REAL *restrict alpha_sqrtg_Ay_interp,
                   CCTK_REAL *restrict alpha_sqrtg_Az_interp,
                   CCTK_REAL *restrict phitilde_rhs,
                   CCTK_REAL *restrict Ax_rhs,
                   CCTK_REAL *restrict Ay_rhs,
                   CCTK_REAL *restrict Az_rhs) {

  const cGH *cctkGH = (const cGH *)void_cctkGH;
  int *bounds = cctkGH->cctk_lsh;
  int *ghostzones = cctkGH->cctk_nghostzones;

  const CCTK_REAL *dX = (const CCTK_REAL *)void_dX;
  const CCTK_REAL **in_vars = (const CCTK_REAL **)void_in_vars;
  const CCTK_REAL *phitilde = (const CCTK_REAL *)void_phitilde;

  /* Compute \partial_t psi6phi = -\partial_i (  \alpha psi^6 A^i - psi6phi \beta^i)
   *    (Eq 13 of http://arxiv.org/pdf/1110.4633.pdf), using Lorenz gauge.
   * Note that the RHS consists of a shift advection term on psi6phi and
   *    a term depending on the vector potential.
   * psi6phi is defined at (i+1/2,j+1/2,k+1/2), but instead of reconstructing
   *    to compute the RHS of \partial_t psi6phi, we instead use standard
   *    interpolations.
   */

  // The stencil here is {-1,1},{-1,1},{-1,1} for x,y,z directions, respectively.
  //     Note that ALL input variables are defined at ALL gridpoints, so no
  //     worries about ghostzones.

#pragma omp parallel for
  for(int k=1;k<bounds[2]-1;k++) for(int j=1;j<bounds[1]-1;j++) for(int i=1;i<bounds[0]-1;i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k);

        induction_gauge gauge_vars;
        induction_gauge_rhs gauge_rhs_vars;

        // Read in variable at interp. stencil points from main memory, store in INTERP_VARS.
        gauge_vars.phitilde = phitilde[index];
        for(int iterz=0; iterz<2; iterz++)
        for(int itery=0; itery<2; itery++)
        for(int iterx=0; iterx<2; iterx++) {
          gauge_vars.gupxx[iterz][itery][iterx]  = in_vars[GUPXXI][CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz)];
          gauge_vars.gupxy[iterz][itery][iterx]  = in_vars[GUPXYI][CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz)];
          gauge_vars.gupxz[iterz][itery][iterx]  = in_vars[GUPXZI][CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz)];
          gauge_vars.gupyy[iterz][itery][iterx]  = in_vars[GUPYYI][CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz)];
          gauge_vars.gupyz[iterz][itery][iterx]  = in_vars[GUPYZI][CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz)];
          gauge_vars.gupzz[iterz][itery][iterx]  = in_vars[GUPZZI][CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz)];
          gauge_vars.lapm1[iterz][itery][iterx]  = in_vars[LAPM1I][CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz)];
          gauge_vars.psi[iterz][itery][iterx]    = in_vars[PSII][CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz)];
          gauge_vars.shiftx[iterz][itery][iterx] = in_vars[SHIFTXI][CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz)];
          gauge_vars.shifty[iterz][itery][iterx] = in_vars[SHIFTYI][CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz)];
          gauge_vars.shiftz[iterz][itery][iterx] = in_vars[SHIFTZI][CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz)];
          //A_x needs a stencil s.t. interp_limits={0,1,-1,1,-1,1}:
          gauge_vars.A_x[iterz+1][itery+1][iterx] = in_vars[A_XI][CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz)];
          //A_y needs a stencil s.t. interp_limits={-1,1,0,1,-1,1}:
          gauge_vars.A_y[iterz+1][itery][iterx+1] = in_vars[A_XI][CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz)];
          //A_z needs a stencil s.t. interp_limits={-1,1,-1,1,0,1}:
          gauge_vars.A_z[iterz][itery+1][iterx+1] = in_vars[A_ZI][CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz)];
        }
        for(int iter=0; iter<2; iter++) {
          gauge_vars.A_x[0][0][iter] = in_vars[A_XI][CCTK_GFINDEX3D(cctkGH,i+iter,j,k)];
          gauge_vars.A_y[0][iter][0] = in_vars[A_YI][CCTK_GFINDEX3D(cctkGH,i,j+iter,k)];
          gauge_vars.A_z[iter][0][0] = in_vars[A_ZI][CCTK_GFINDEX3D(cctkGH,i,j,k+iter)];
        }

        interpolate_for_A_i_rhs(&gauge_vars, &gauge_rhs_vars);

        alpha_interp[index] = gauge_rhs_vars.alpha_interp;
        alpha_sqrtg_Ax_interp[index] = gauge_rhs_vars.alpha_sqrtg_Ax_interp[0];
        alpha_sqrtg_Az_interp[index] = gauge_rhs_vars.alpha_sqrtg_Ay_interp[0];
        alpha_sqrtg_Ay_interp[index] = gauge_rhs_vars.alpha_sqrtg_Az_interp[0];
        alpha_Phi_minus_betaj_A_j_interp[index] = gauge_rhs_vars.alpha_Phi_minus_betaj_A_j_interp[0];
        shiftx_interp[index] = gauge_rhs_vars.shiftx_interp[0];
        shifty_interp[index] = gauge_rhs_vars.shifty_interp[0];
        shiftz_interp[index] = gauge_rhs_vars.shiftz_interp[0];
      }

#pragma omp parallel for
  for(int k=ghostzones[2];k<bounds[2]-ghostzones[2];k++) for(int j=ghostzones[1];j<bounds[1]-ghostzones[1];j++) for(int i=ghostzones[0];i<bounds[0]-ghostzones[0];i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        induction_gauge_rhs gauge_rhs_vars;

        gauge_rhs_vars.alpha_interp = alpha_interp[index];;

        gauge_rhs_vars.alpha_Phi_minus_betaj_A_j_interp[0] = alpha_Phi_minus_betaj_A_j_interp[index];
        gauge_rhs_vars.alpha_Phi_minus_betaj_A_j_interp[1] = alpha_Phi_minus_betaj_A_j_interp[CCTK_GFINDEX3D(cctkGH,i-1,j,k)];
        gauge_rhs_vars.alpha_Phi_minus_betaj_A_j_interp[2] = alpha_Phi_minus_betaj_A_j_interp[CCTK_GFINDEX3D(cctkGH,i,j-1,k)];
        gauge_rhs_vars.alpha_Phi_minus_betaj_A_j_interp[3] = alpha_Phi_minus_betaj_A_j_interp[CCTK_GFINDEX3D(cctkGH,i,j,k-1)];

        gauge_rhs_vars.alpha_sqrtg_Ax_interp[0] = alpha_sqrtg_Ax_interp[index];
        gauge_rhs_vars.alpha_sqrtg_Ay_interp[0] = alpha_sqrtg_Ay_interp[index];
        gauge_rhs_vars.alpha_sqrtg_Az_interp[0] = alpha_sqrtg_Az_interp[index];
        gauge_rhs_vars.alpha_sqrtg_Ax_interp[1] = alpha_sqrtg_Ax_interp[CCTK_GFINDEX3D(cctkGH,i+1,j,k)];
        gauge_rhs_vars.alpha_sqrtg_Ay_interp[1] = alpha_sqrtg_Ay_interp[CCTK_GFINDEX3D(cctkGH,i,j+1,k)];
        gauge_rhs_vars.alpha_sqrtg_Az_interp[1] = alpha_sqrtg_Az_interp[CCTK_GFINDEX3D(cctkGH,i,j,k+1)];

        for(int iter=-2; iter<3; iter++) {
          int indexx = CCTK_GFINDEX3D(cctkGH,i+iter,j,     k     );
          int indexy = CCTK_GFINDEX3D(cctkGH,i,     j+iter,k     );
          int indexz = CCTK_GFINDEX3D(cctkGH,i,     j,     k+iter);
          gauge_rhs_vars.shiftx_interp[iter+2] = shiftx_interp[indexx];
          gauge_rhs_vars.shifty_interp[iter+2] = shifty_interp[indexy];
          gauge_rhs_vars.shiftz_interp[iter+2] = shiftz_interp[indexz];
          gauge_rhs_vars.phitildex[iter+2] = phitilde[indexx];
          gauge_rhs_vars.phitildey[iter+2] = phitilde[indexy];
          gauge_rhs_vars.phitildez[iter+2] = phitilde[indexz];
        }
        calculate_phitilde_and_A_i_rhs(dX, Lorenz_damping_factor, &gauge_rhs_vars);
        phitilde_rhs[index] = gauge_rhs_vars.phitilde_rhs;
      }
}
