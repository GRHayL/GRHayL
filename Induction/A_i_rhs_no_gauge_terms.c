#include "induction.h"

/* Compute the part of A_i_rhs that excludes the gauge terms. I.e., we set
 *   A_i_rhs = \partial_t A_i = \psi^{6} (v^z B^x - v^x B^z)   here.
 */
double A_i_rhs_no_gauge_terms( const double psi6,
                             induction_lr *restrict vars ) {

  double A3_rhs_rr = psi6*(vars->v1rr*vars->B2r - vars->v2rr*vars->B1r);
  double A3_rhs_rl = psi6*(vars->v1rl*vars->B2r - vars->v2rl*vars->B1l);
  double A3_rhs_lr = psi6*(vars->v1lr*vars->B2l - vars->v2lr*vars->B1r);
  double A3_rhs_ll = psi6*(vars->v1ll*vars->B2l - vars->v2ll*vars->B1l);

  // All variables for the A_i_rhs computation are now at the appropriate staggered point,
  //   so it's time to compute the HLL flux!

  double B1tilder_minus_B1tildel = psi6*( vars->B1r - vars->B1l );
  double B2tilder_minus_B2tildel = psi6*( vars->B2r - vars->B2l );

  /*---------------------------
   * Implement 2D HLL flux
   * [see Del Zanna, Bucciantini & Londrillo A&A 400, 397 (2003), Eq. (44)]
   *
   * Note that cmax/cmin (\alpha^{\pm}  as defined in that paper) is at a slightly DIFFERENT
   * point than that described in the Del Zanna et al paper (e.g., (i+1/2,j,k) instead of
   * (i+1/2,j+1/2,k) for F3).  Yuk Tung Liu discussed this point with M. Shibata,
   * who found that the effect is negligible.
   ---------------------------*/

  double c_array[4] = {vars->c1_min*vars->c2_min, vars->c1_min*vars->c2_max, vars->c1_max*vars->c2_min, vars->c1_max*vars->c2_max};
  double c_sums[2] = {vars->c1_min+vars->c1_max, vars->c2_min+vars->c2_max};
  double c_squares[2] = {vars->c1_min*vars->c1_max, vars->c2_min*vars->c2_max};

  return (c_array[3]*A3_rhs_ll + c_array[2]*A3_rhs_lr
          + c_array[1]*A3_rhs_rl + c_array[0]*A3_rhs_rr)
         /( c_sums[0]*c_sums[1] )
         - c_squares[0]*(B2tilder_minus_B2tildel)/c_sums[0]
         + c_squares[1]*(B1tilder_minus_B1tildel)/c_sums[1];
}
