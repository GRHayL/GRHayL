#include "induction.h"

static int MINUS1=0, PLUS0=1, PLUS1=2;

static inline double avg3(const double f[3][3][3],const int imin,const int imax, const int jmin,const int jmax, const int kmin,const int kmax);
static inline double avg2(const double f[2][2][2],const int imin,const int imax, const int jmin,const int jmax, const int kmin,const int kmax);

/* Function    : interpolate_for_A_gauge_rhs()
 * Description : computes several interpolated quantities for computing the RHS
 *               for tilde{phi} and the gauge contributions to A_i; these are
 *               used in calculate_phitilde_and_A_gauge_rhs() function
 *
 * Inputs      : gauge_vars      - A_gauge_vars struct containing stencils for 
 *                                 variables to compute interpolated values
 *
 * Outputs     : gauge_rhs_vars  - A_gauge_rhs_vars struct containing interpolated
 *                                 values for use in
 *
 */
void interpolate_for_A_gauge_rhs(
            const A_gauge_vars *restrict gauge_vars,
            A_gauge_rhs_vars *restrict gauge_rhs_vars ) {
  /* Compute \partial_t psi6phi = -\partial_i (  \alpha psi^6 A^i - psi6phi \beta^i)
   *    (Eq 13 of http://arxiv.org/pdf/1110.4633.pdf), using Lorenz gauge.
   * Note that the RHS consists of a shift advection term on psi6phi and
   *    a term depending on the vector potential.
   * psi6phi is defined at (i+1/2,j+1/2,k+1/2), but instead of reconstructing
   *    to compute the RHS of \partial_t psi6phi, we instead use standard
   *    interpolations.
   */

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

  // First \alpha at (i+1/2,j+1/2,k+1/2). Will come in handy when computing damping term later.
  gauge_rhs_vars->alpha_interp = avg2(gauge_vars->lapse , 0,1, 0,1, 0,1);

  // Next, do A^X TERM (interpolate to (i,j+1/2,k+1/2) )
  // \alpha \sqrt{\gamma} A^x = \alpha psi^6 A^x (RHS of \partial_i psi6phi)
  // Note that gupij is \tilde{\gamma}^{ij}, so we need to multiply by \psi^{-4}.
  const double gupxx_jk = avg2(gauge_vars->gupxx, 0,0, 0,1, 0,1);
  const double gupxy_jk = avg2(gauge_vars->gupxy, 0,0, 0,1, 0,1);
  const double gupxz_jk = avg2(gauge_vars->gupxz, 0,0, 0,1, 0,1);

  double lapse_psi2[2][2][2];
  double lapse_over_psi6[2][2][2];
  for(int kk=0;kk<=1;kk++) for(int jj=0;jj<=1;jj++) for(int ii=0;ii<=1;ii++) {
        const double Psi2 = gauge_vars->psi[kk][jj][ii]*gauge_vars->psi[kk][jj][ii];
        const double alpha = gauge_vars->lapse[kk][jj][ii];
        lapse_psi2[kk][jj][ii]=alpha*Psi2;
        lapse_over_psi6[kk][jj][ii]=alpha/(Psi2*Psi2*Psi2);
      }

  const double lapse_Psi2_jk = avg2(lapse_psi2, 0,0, 0,1, 0,1);

  const double A_x_jk   = gauge_vars->A_x[1][1][1]; // @ (i,j+1/2,k+1/2)
  const double A_y_jk   = avg3(gauge_vars->A_y,MINUS1,PLUS0, PLUS0,PLUS1, PLUS0,PLUS0); // @ (i+1/2,j,k+1/2)
  const double A_z_jk   = avg3(gauge_vars->A_z,MINUS1,PLUS0, PLUS0,PLUS0, PLUS0,PLUS1); // @ (i+1/2,j+1/2,k)

  gauge_rhs_vars->alpha_sqrtg_Ax_interp[0] = lapse_Psi2_jk*
    ( gupxx_jk*A_x_jk + gupxy_jk*A_y_jk + gupxz_jk*A_z_jk );

  // DO A^Y TERM (interpolate to (i+1/2,j,k+1/2) )
  // \alpha \sqrt{\gamma} A^y = \alpha psi^6 A^y (RHS of \partial_i psi6phi)
  // Note that gupij is \tilde{\gamma}^{ij}, so we need to multiply by \psi^{-4}.
  const double gupxy_ik = avg2(gauge_vars->gupxy, 0,1, 0,0, 0,1);
  const double gupyy_ik = avg2(gauge_vars->gupyy, 0,1, 0,0, 0,1);
  const double gupyz_ik = avg2(gauge_vars->gupyz, 0,1, 0,0, 0,1);

  const double lapse_Psi2_ik = avg2(lapse_psi2, 0,1, 0,0, 0,1);

  const double A_x_ik   = avg3(gauge_vars->A_x, PLUS0,PLUS1,MINUS1,PLUS0, PLUS0,PLUS0); // @ (i,j+1/2,k+1/2)
  const double A_y_ik   = gauge_vars->A_y[1][1][1]; // @ (i+1/2,j,k+1/2)
  const double A_z_ik   = avg3(gauge_vars->A_z, PLUS0,PLUS0,MINUS1,PLUS0, PLUS0,PLUS1); // @ (i+1/2,j+1/2,k)

  gauge_rhs_vars->alpha_sqrtg_Ay_interp[0] = lapse_Psi2_ik*
    ( gupxy_ik*A_x_ik + gupyy_ik*A_y_ik + gupyz_ik*A_z_ik );

  // DO A^Z TERM (interpolate to (i+1/2,j+1/2,k) )
  // \alpha \sqrt{\gamma} A^z = \alpha psi^6 A^z (RHS of \partial_i psi6phi)
  // Note that gupij is \tilde{\gamma}^{ij}, so we need to multiply by \psi^{-4}.
  const double gupxz_ij = avg2(gauge_vars->gupxz, 0,1, 0,1, 0,0);
  const double gupyz_ij = avg2(gauge_vars->gupyz, 0,1, 0,1, 0,0);
  const double gupzz_ij = avg2(gauge_vars->gupzz, 0,1, 0,1, 0,0);

  const double lapse_Psi2_ij = avg2(lapse_psi2, 0,1, 0,1, 0,0);

  const double A_x_ij   = avg3(gauge_vars->A_x, PLUS0,PLUS1, PLUS0,PLUS0,MINUS1,PLUS0); // @ (i,j+1/2,k+1/2)
  const double A_y_ij   = avg3(gauge_vars->A_y, PLUS0,PLUS0, PLUS0,PLUS1,MINUS1,PLUS0); // @ (i+1/2,j,k+1/2)
  const double A_z_ij   = gauge_vars->A_z[1][1][1]; // @ (i+1/2,j+1/2,k)

  gauge_rhs_vars->alpha_sqrtg_Az_interp[0] = lapse_Psi2_ij*
    ( gupxz_ij*A_x_ij + gupyz_ij*A_y_ij + gupzz_ij*A_z_ij );


  // Next set \alpha \Phi - \beta^j A_j at (i+1/2,j+1/2,k+1/2):
  //   We add a "L" suffix to shifti_ijk to denote "LOCAL", as we set
  //      shifti_ijk[] gridfunction below.
  gauge_rhs_vars->shiftx_interp[0] = avg2(gauge_vars->shiftx, 0,1, 0,1, 0,1);
  gauge_rhs_vars->shifty_interp[0] = avg2(gauge_vars->shifty, 0,1, 0,1, 0,1);
  gauge_rhs_vars->shiftz_interp[0] = avg2(gauge_vars->shiftz, 0,1, 0,1, 0,1);
  const double  lapse_over_Psi6_ijk = avg2(lapse_over_psi6, 0,1, 0,1, 0,1);
  const double  A_x_ijk = avg3(gauge_vars->A_x, PLUS0,PLUS1, PLUS0,PLUS0, PLUS0,PLUS0); // @ (i,j+1/2,k+1/2)
  const double  A_y_ijk = avg3(gauge_vars->A_y, PLUS0,PLUS0, PLUS0,PLUS1, PLUS0,PLUS0); // @ (i+1/2,j,k+1/2)
  const double  A_z_ijk = avg3(gauge_vars->A_z, PLUS0,PLUS0, PLUS0,PLUS0, PLUS0,PLUS1); // @ (i+1/2,j+1/2,k)

  gauge_rhs_vars->alpha_Phi_minus_betaj_A_j_interp[0] = gauge_vars->phitilde*lapse_over_Psi6_ijk
                                                      - (gauge_rhs_vars->shiftx_interp[0]*A_x_ijk
                                                         + gauge_rhs_vars->shifty_interp[0]*A_y_ijk
                                                         + gauge_rhs_vars->shiftz_interp[0]*A_z_ijk);
}

static inline double avg3(const double f[3][3][3], const int imin, const int imax, const int jmin, const int jmax, const int kmin, const int kmax) {
  double retval=0.0,num_in_sum=0.0;
  for(int kk=kmin;kk<=kmax;kk++) for(int jj=jmin;jj<=jmax;jj++) for(int ii=imin;ii<=imax;ii++) {
        retval+=f[kk][jj][ii]; num_in_sum++;
      }
  return retval/num_in_sum;
}

static inline double avg2(const double f[2][2][2], const int imin, const int imax, const int jmin, const int jmax, const int kmin, const int kmax) {
  double retval=0.0,num_in_sum=0.0;
  for(int kk=kmin;kk<=kmax;kk++) for(int jj=jmin;jj<=jmax;jj++) for(int ii=imin;ii<=imax;ii++) {
        retval+=f[kk][jj][ii]; num_in_sum++;
      }
  return retval/num_in_sum;
}
