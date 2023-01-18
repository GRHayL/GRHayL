#include "induction.h"

/* Compute the part of A_i_rhs that excludes the gauge terms. I.e., we set
 *   A_i_rhs = \partial_t A_i = \psi^{6} (v^z B^x - v^x B^z)   here.
 */
double A_no_gauge_rhs_stencil(const A_no_gauge_vars *restrict vars) {

  const double A3_rhs_rr = vars->psi6*(vars->v1rr*vars->B2r - vars->v2rr*vars->B1r);
  const double A3_rhs_rl = vars->psi6*(vars->v1rl*vars->B2r - vars->v2rl*vars->B1l);
  const double A3_rhs_lr = vars->psi6*(vars->v1lr*vars->B2l - vars->v2lr*vars->B1r);
  const double A3_rhs_ll = vars->psi6*(vars->v1ll*vars->B2l - vars->v2ll*vars->B1l);

  // All variables for the A_i_rhs computation are now at the appropriate staggered point,
  //   so it's time to compute the HLL flux!

  const double B1tilder_minus_B1tildel = vars->psi6*( vars->B1r - vars->B1l );
  const double B2tilder_minus_B2tildel = vars->psi6*( vars->B2r - vars->B2l );

  /*---------------------------
   * Implement 2D HLL flux
   * [see Del Zanna, Bucciantini & Londrillo A&A 400, 397 (2003), Eq. (44)]
   *
   * Note that cmax/cmin (\alpha^{\pm}  as defined in that paper) is at a slightly DIFFERENT
   * point than that described in the Del Zanna et al paper (e.g., (i+1/2,j,k) instead of
   * (i+1/2,j+1/2,k) for F3).  Yuk Tung Liu discussed this point with M. Shibata,
   * who found that the effect is negligible.
   ---------------------------*/

  const double c1_sum = vars->c1_min+vars->c1_max;
  const double c2_sum = vars->c2_min+vars->c2_max;

  return (  vars->c2_max*vars->c1_max*A3_rhs_ll
          + vars->c2_min*vars->c1_max*A3_rhs_lr
          + vars->c2_max*vars->c1_min*A3_rhs_rl
          + vars->c2_min*vars->c1_min*A3_rhs_rr)
         /( c2_sum*c1_sum )
         - vars->c1_min*vars->c1_max*(B2tilder_minus_B2tildel)/c1_sum
         + vars->c2_min*vars->c2_max*(B1tilder_minus_B1tildel)/c2_sum;
}
