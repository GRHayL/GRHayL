#include "ghl_induction.h"

/**
 * @ingroup mag_flux
 * @brief Compute RHS for \f$ A_i \f$, excluding gauge contributions; e.g.
 * @sp10 \f$ A_y^\mathrm{rhs} = \partial_t A_y = \psi^{6} (v^z B^x - v^x B^z) \f$
 *
 * @details
 * This function computes the right-hand side of the induction equation for a
 * single direction of \f$ A_i \f$ using the undensitized magnetic field
 * \f$ B^i \f$. The specific direction \f$ i \f$ is based on the data in @ref ghl_HLL_vars.
 *
 * As this nearly is identical to @ref ghl_HLL_flux_with_Btilde, we only note that
 * this function differs by requiring \f$ \psi^6 \f$ to compute the returned
 * value. The final equation is just multiplied by this factor.
 *
 * @param[in] psi6 spacetime quantity \f$ \psi^6 = \sqrt{|\gamma|} \f$ at the
 *                  gridpoint of \f$ A_i \f$
 *
 * @param[in] vars pointer to a ghl_HLL_vars struct
 *
 * @returns flux contribution to \f$ A_i^\mathrm{rhs} \f$
 */
double ghl_HLL_flux_with_B(
      const double psi6,
      const ghl_HLL_vars *restrict vars) {

  const double c1_sum = vars->c1_min+vars->c1_max;
  const double c2_sum = vars->c2_min+vars->c2_max;

  /*
    To compute A_i_rhs, we use the HLL flux from Eq. 44 of
      Del Zanna, Bucciantini & Londrillo A&A 400, 397 (2003)
    To explain the terms, let's consider the flux RHS for the z component. Then,
    1->x, 2->y, and 3->z. We first compute the B field contributions
    Bxterm = \psi^6 * cy_min * cy_max (Bx^R - Bx^L)/ (cy_min + cy_max)
    Byterm = \psi^6 * cx_min * cx_max (By^R - By^L)/ (cx_min + cx_max)
  */
  const double B1term = vars->c2_min * vars->c2_max * (vars->B1r - vars->B1l) / c2_sum;
  const double B2term = vars->c1_min * vars->c1_max * (vars->B2r - vars->B2l) / c1_sum;

  /*
     Additionally, we need the electric field contribution
     E_z = -(v_x B_y - v_y B_x)
     That means we need to compute
     E^LL, E^LR, E^RL, and E^RR.
  */
  const double A3_rhs_rr = vars->v1rr*vars->B2r - vars->v2rr*vars->B1r;
  const double A3_rhs_rl = vars->v1rl*vars->B2r - vars->v2rl*vars->B1l;
  const double A3_rhs_lr = vars->v1lr*vars->B2l - vars->v2lr*vars->B1r;
  const double A3_rhs_ll = vars->v1ll*vars->B2l - vars->v2ll*vars->B1l;

  return psi6*(B1term - B2term
              + ( vars->c2_max*vars->c1_max*A3_rhs_ll
                + vars->c2_min*vars->c1_max*A3_rhs_lr
                + vars->c2_max*vars->c1_min*A3_rhs_rl
                + vars->c2_min*vars->c1_min*A3_rhs_rr)
                /(c2_sum*c1_sum));
}
