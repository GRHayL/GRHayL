#include "cctk.h"
#include "IGM.h"

void phitilde_and_A_gauge_rhs(const cGH *cctkGH,
                   const double *restrict dX,
                   const double *restrict gupxx,
                   const double *restrict gupxy,
                   const double *restrict gupxz,
                   const double *restrict gupyy,
                   const double *restrict gupyz,
                   const double *restrict gupzz,
                   const double *restrict psi,
                   const double *restrict lapse,
                   const double *restrict betax,
                   const double *restrict betay,
                   const double *restrict betaz,
                   const double *restrict Ax,
                   const double *restrict Ay,
                   const double *restrict Az,
                   const double *restrict phitilde,
                   double Lorenz_damping_factor,
                   double *restrict shiftx_interp,
                   double *restrict shifty_interp,
                   double *restrict shiftz_interp,
                   double *restrict alpha_interp,
                   double *restrict alpha_Phi_minus_betaj_A_j_interp,
                   double *restrict alpha_sqrtg_Ax_interp,
                   double *restrict alpha_sqrtg_Ay_interp,
                   double *restrict alpha_sqrtg_Az_interp,
                   double *restrict phitilde_rhs,
                   double *restrict Ax_rhs,
                   double *restrict Ay_rhs,
                   double *restrict Az_rhs) {

  /* Compute \partial_t psi6phi = -\partial_i (  \alpha psi^6 A^i - psi6phi \beta^i)
   *    (Eq 13 of http://arxiv.org/pdf/1110.4633.pdf), using Lorenz gauge.
   * Note that the RHS consists of a shift advection term on psi6phi and
   *    a term depending on the vector potential.
   * psi6phi is defined at (i+1/2,j+1/2,k+1/2), but instead of reconstructing
   *    to compute the RHS of \partial_t psi6phi, we instead use standard
   *    interpolations.
   */

  // We declare these values to be over the interior so the setting of ghostzone points is
  // more transparent.
  const int imin = cctkGH->cctk_nghostzones[0];
  const int jmin = cctkGH->cctk_nghostzones[1];
  const int kmin = cctkGH->cctk_nghostzones[2];
  const int imax = cctkGH->cctk_lsh[0]-cctkGH->cctk_nghostzones[0];
  const int jmax = cctkGH->cctk_lsh[1]-cctkGH->cctk_nghostzones[1];
  const int kmax = cctkGH->cctk_lsh[2]-cctkGH->cctk_nghostzones[2];

  // The RHS loop requires 2 ghostzones, so this loop must set those values, hence this loop
  // goes into the ghostzones. This loop requires a stencil of {-1,1},{-1,1},{-1,1},
  // which uses the last of the 3 ghostzones required by the simulation.
  // Note that ALL input variables are defined at ALL gridpoints, so no
  // worries about ghostzones.
#pragma omp parallel for
  for(int k=1; k<cctkGH->cctk_lsh[2]-1; k++)
    for(int j=1; j<cctkGH->cctk_lsh[1]-1; j++)
      for(int i=1; i<cctkGH->cctk_lsh[0]-1; i++) {
//  for(int k=kmin-2; k<kmax+2; k++)
//    for(int j=jmin-2; j<jmax+2; j++)
//      for(int i=imin-2; i<imax+2; i++) {
        const int index=CCTK_GFINDEX3D(cctkGH,i,j,k);

        // First compute \partial_j \alpha \sqrt{\gamma} A^j (RHS of \partial_i psi6phi)
        // FIXME: Would be much cheaper & easier to unstagger A_i, raise, then interpolate A^i.
        //        However, we keep it this way to be completely compatible with the original
        //        Illinois GRMHD thorn, called mhd_evolve.
        //
        //Step 1) j=x: Need to raise A_i, but to do that, we must have all variables at the same gridpoints:
        // The goal is to compute \partial_j (\alpha \sqrt{\gamma} A^j) at (i+1/2,j+1/2,k+1/2)
        //    We do this by first interpolating (RHS1x) = (\alpha \sqrt{\gamma} A^x) at
        //    (i,j+1/2,k+1/2)and (i+1,j+1/2,k+1/2), then taking \partial_x (RHS1x) =
        //    [ RHS1x(i+1,j+1/2,k+1/2) - RHS1x(i,j+1/2,k+1/2) ]/dX.
        // First bring gup's, psi, and alpha to (i,j+1/2,k+1/2):
        A_gauge_vars gauge_vars;
        A_gauge_rhs_vars gauge_rhs_vars;

        // Read in variable at interpolation stencil points from main memory.
        gauge_vars.phitilde = phitilde[index];
        for(int iterz=0; iterz<2; iterz++)
          for(int itery=0; itery<2; itery++)
            for(int iterx=0; iterx<2; iterx++) {
              const int ind = CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz);
              gauge_vars.gupxx[iterz][itery][iterx]  = gupxx[ind];
              gauge_vars.gupxy[iterz][itery][iterx]  = gupxy[ind];
              gauge_vars.gupxz[iterz][itery][iterx]  = gupxz[ind];
              gauge_vars.gupyy[iterz][itery][iterx]  = gupyy[ind];
              gauge_vars.gupyz[iterz][itery][iterx]  = gupyz[ind];
              gauge_vars.gupzz[iterz][itery][iterx]  = gupzz[ind];
              gauge_vars.lapse[iterz][itery][iterx]  = lapse[ind];
              gauge_vars.psi[iterz][itery][iterx]    = psi[ind];
              gauge_vars.shiftx[iterz][itery][iterx] = betax[ind];
              gauge_vars.shifty[iterz][itery][iterx] = betay[ind];
              gauge_vars.shiftz[iterz][itery][iterx] = betaz[ind];
        }
        // A_x needs a stencil s.t. interp_limits={ 0,1,-1,1,-1,1}.
        // A_y needs a stencil s.t. interp_limits={-1,1, 0,1,-1,1}.
        // A_z needs a stencil s.t. interp_limits={-1,1,-1,1, 0,1}.
        // We could fill only the needed elements, but it is cleaner
        // to fill in the whole 3x3x3 array.
        // TODO: the old code explicitly only filled in the necessary
        // elements. If we want to remove ~15 memcopies, do that here.
        for(int iterz=-1; iterz<2; iterz++)
          for(int itery=-1; itery<2; itery++)
            for(int iterx=-1; iterx<2; iterx++) {
              const int ind = CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz);
              gauge_vars.A_x[iterz+1][itery+1][iterx+1] = Ax[ind];
              gauge_vars.A_y[iterz+1][itery+1][iterx+1] = Ay[ind];
              gauge_vars.A_z[iterz+1][itery+1][iterx+1] = Az[ind];
        }
// This code should only copy the needed data that isn't copied in the loop for other variables, but it is untested.
//        for(int iter2=0; iter2<2; iter2++)
//        for(int iter1=0; iter1<2; iter1++) {
//          gauge_vars.A_x[iter2+1][0][iter1+1] = in_vars[A_XI][CCTK_GFINDEX3D(cctkGH, i+iter1,     j-1, k+iter2)]; // { (0,1),    -1, (0,1)}
//          gauge_vars.A_x[0][iter2+1][iter1+1] = in_vars[A_XI][CCTK_GFINDEX3D(cctkGH, i+iter1, j+iter2, k-1    )]; // { (0,1), (0,1), -1}
//          gauge_vars.A_y[iter2+1][iter1+1][0] = in_vars[A_YI][CCTK_GFINDEX3D(cctkGH,     i-1, j+iter1, k+iter2)]; // { -1,    (0,1), (0,1)}
//          gauge_vars.A_y[0][iter1+1][iter2+1] = in_vars[A_YI][CCTK_GFINDEX3D(cctkGH, i+iter2, j+iter1, k-1    )]; // { (0,1), (0,1), -1}
//          gauge_vars.A_z[iter1+1][iter2+1][0] = in_vars[A_ZI][CCTK_GFINDEX3D(cctkGH,     i-1, j+iter2, k+iter1)]; // { -1,    (0,1), (0,1)}
//          gauge_vars.A_z[iter1+1][0][iter2+1] = in_vars[A_ZI][CCTK_GFINDEX3D(cctkGH, i+iter2,     j-1, k+iter1)]; // { (0,1),    -1, (0,1)}
//        }

        interpolate_for_A_gauge_rhs(&gauge_vars, &gauge_rhs_vars);

        alpha_interp[index] = gauge_rhs_vars.alpha_interp;
        alpha_sqrtg_Ax_interp[index] = gauge_rhs_vars.alpha_sqrtg_Ax_interp[0];
        alpha_sqrtg_Ay_interp[index] = gauge_rhs_vars.alpha_sqrtg_Ay_interp[0];
        alpha_sqrtg_Az_interp[index] = gauge_rhs_vars.alpha_sqrtg_Az_interp[0];
        alpha_Phi_minus_betaj_A_j_interp[index] = gauge_rhs_vars.alpha_Phi_minus_betaj_A_j_interp[0];
        shiftx_interp[index] = gauge_rhs_vars.shiftx_interp[0];
        shifty_interp[index] = gauge_rhs_vars.shifty_interp[0];
        shiftz_interp[index] = gauge_rhs_vars.shiftz_interp[0];
      }

  const double dxinv[3] = {1.0/dX[0], 1.0/dX[1], 1.0/dX[2]};

#pragma omp parallel for
  for(int k=cctkGH->cctk_nghostzones[2]; k<cctkGH->cctk_lsh[2]-cctkGH->cctk_nghostzones[2]; k++)
    for(int j=cctkGH->cctk_nghostzones[1]; j<cctkGH->cctk_lsh[1]-cctkGH->cctk_nghostzones[1]; j++)
      for(int i=cctkGH->cctk_nghostzones[0]; i<cctkGH->cctk_lsh[0]-cctkGH->cctk_nghostzones[0]; i++) {
//  for(int k=kmin; k<kmax; k++)
//    for(int j=jmin; j<jmax; j++)
//      for(int i=imin; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // \partial_t A_i = [reconstructed stuff] + [gauge stuff],
        //    where [gauge stuff] = -\partial_i (\alpha \Phi - \beta^j A_j)
        A_gauge_rhs_vars gauge_rhs_vars;
    
        gauge_rhs_vars.alpha_interp = alpha_interp[index];
    
        gauge_rhs_vars.dxi[0] = dxinv[0];
        gauge_rhs_vars.dxi[1] = dxinv[1];
        gauge_rhs_vars.dxi[2] = dxinv[2];
    
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
          const int indexx = CCTK_GFINDEX3D(cctkGH,i+iter,j,     k     );
          const int indexy = CCTK_GFINDEX3D(cctkGH,i,     j+iter,k     );
          const int indexz = CCTK_GFINDEX3D(cctkGH,i,     j,     k+iter);
          gauge_rhs_vars.shiftx_interp[iter+2] = shiftx_interp[indexx];
          gauge_rhs_vars.shifty_interp[iter+2] = shifty_interp[indexy];
          gauge_rhs_vars.shiftz_interp[iter+2] = shiftz_interp[indexz];
          gauge_rhs_vars.phitildex[iter+2] = phitilde[indexx];
          gauge_rhs_vars.phitildey[iter+2] = phitilde[indexy];
          gauge_rhs_vars.phitildez[iter+2] = phitilde[indexz];
        }
        calculate_phitilde_and_A_gauge_rhs(Lorenz_damping_factor, &gauge_rhs_vars);
    
        phitilde_rhs[index] = gauge_rhs_vars.phitilde_rhs;
        Ax_rhs[index] += gauge_rhs_vars.A_x_gauge_rhs;
        Ay_rhs[index] += gauge_rhs_vars.A_y_gauge_rhs;
        Az_rhs[index] += gauge_rhs_vars.A_z_gauge_rhs;
  }
}
