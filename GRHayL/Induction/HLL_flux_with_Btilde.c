#include "ghl_induction.h"

/**
 * @ingroup mag_flux
 * @brief Compute RHS for \f$ A_i \f$, excluding gauge contributions; e.g.
 * @sp10 \f$ A_y^\mathrm{rhs} = \partial_t A_y = v^z \tilde{B}^x - v^x \tilde{B}^z \f$
 *
 * @details
 * This function computes the right-hand side of the induction equation for a
 * single direction of \f$ A_i \f$ using the densitized magnetic field
 * \f$ \tilde{B}^i \f$.
 *
 * The specific direction \f$ i \f$ is based on the data in @ref ghl_HLL_vars.
 * Since filling this struct is less trivial than for most functions, we
 * give a specific example assuming that we are calculating the \f$ A_z \f$ term.
 * In this case,
 *
 * \f[
 * v_1 = v_x \\
 * v_2 = v_y \\
 * B_1 = \tilde{B}_\mathrm{stagger}^x \\
 * B_2 = \tilde{B}_\mathrm{stagger}^y
 * \f]
 *
 * Note that the velocity and magnetic fields are always vectors (\f$ B^i \f$),
 * and we simply use lowered integers like \f$ B_2 \f$ to avoid confusing labels
 * with exponents. Now, in this example, we see that all the "1" variables
 * are the \f$ x \f$ components, and the "2" variables are the \f$ y \f$
 * components. Since the velocities are cell-centered---which we can refer
 * to as having index \f$ (i, j, k) \f$---they need to be computed at the
 * correct index for \f$ A_i \f$ (staggered in perpendicular directions).
 * Generally, the hydrodynamic variables cannot simply be interpolated to
 * due the possible presence of shocks. As such, these quantities are
 * reconstructed to faces. Thus, we have a left and right value for a
 * given point. The velocities have to be reconstructed twice, leading
 * to four values.
 *
 * | \f$ A_i^\mathrm{RHS} \f$ Direction | \f$ v \f$ Reconstructed in | \f$ B_1 \f$ reconstructed in | \f$ B_2 \f$ reconstructed in |
 * |:----------------------------------:|:--------------------------:|:----------------------------:|:---------------------------:|
 * | x                                  | y, z                       | z                            | y                           |
 * | y                                  | z, x                       | x                            | z                           |
 * | z                                  | x, y                       | y                            | x                           |
 *
 * The evolution equation for \f$ A_i \f$ is given by
 *
 * \f[
 * \frac{d}{dt} A_i = - E_i
 * \f]
 *
 * where these quantities must be evaluated at the staggered gridpoints. Recall
 * from the @ref Induction page that \f$ A_i \f$ is the vector potential
 * associated with \f$ \tilde{B}^i \f$. Due to this, the equation naturally
 * uses \f$ \tilde{B}^i \f$. Additionally, we know that e.g.
 *
 * \f[
 * A^{LR}_{3,\mathrm{rhs}} = v^1_{LR} B^2_L - v^2_{LR} B^1_R
 * \f]
 *
 * To evaluate this derivative, we use the HLL electric-field flux from
 * equation 44 of \cite DelZanna_2003, negated since
 * \f$ A^\mathrm{rhs}_3 = -E_3 \f$:
 *
 * \f[
 * A^\mathrm{HLL}_{3,\mathrm{rhs}} =
 * \frac{c^+_2 c^-_2}{c_2^\mathrm{sum}} \left( B^R_1 - B^L_1 \right)
 * - \frac{c^+_1 c^-_1}{c_1^\mathrm{sum}} \left( B^R_2 - B^L_2 \right)
 * + \frac{1}{c_1^\mathrm{sum} c_2^\mathrm{sum}}\left(
 * c^+_1 c^+_2 A^{LL}_{3,\mathrm{rhs}}
 * + c^+_1 c^-_2 A^{LR}_{3,\mathrm{rhs}}
 * + c^-_1 c^+_2 A^{RL}_{3,\mathrm{rhs}}
 * + c^-_1 c^-_2 A^{RR}_{3,\mathrm{rhs}} \right)
 * \f]
 *
 * where \f$ c_i^\pm \f$ represents the characteristic speeds returned by
 * @ref ghl_calculate_characteristic_speed_dirn0,
 * @ref ghl_calculate_characteristic_speed_dirn1, and
 * @ref ghl_calculate_characteristic_speed_dirn2, and \f$ c_i^\mathrm{sum} \f$
 * is given by
 *
 * \f[
 * c_i^\mathrm{sum} \equiv c^+_i + c^-_i
 * \f]
 *
 * The number subscripts represent spatial coordinates \f$ {x,y,z} \f$, which for the
 * characteristic speeds correspond to the direction of the reconstruction from which
 * these are taken. For example, \f$ c_2 \f$ is computed using the primitives reconstructed
 * in the \f$ y \f$ direction.
 *
 * The \f$ R \f$ and \f$ L \f$ superscripts refer to the right/left reconstructed values for
 * an edge. As discussed before, the velocities have to be reconstructed twice, and so we have
 * to list the superscript for each reconstruction. This leads to the notation in the
 * equation above and defines our prescription for computing the
 * flux contributions to the RHS of \f$ A_i \f$.
 *
 * @param[in] vars pointer to a ghl_HLL_vars struct
 *
 * @returns flux contribution to \f$ A_i^\mathrm{rhs} \f$
 */
double ghl_HLL_flux_with_Btilde(
      const ghl_HLL_vars *restrict vars) {

  const double c1_sum = vars->c1_min+vars->c1_max;
  const double c2_sum = vars->c2_min+vars->c2_max;

  /*
    To compute A_i_rhs, we use the HLL flux from Eq. 44 of
      Del Zanna, Bucciantini & Londrillo A&A 400, 397 (2003)
    To explain the terms, let's consider the flux RHS for the z component. Then,
    1->x, 2->y, and 3->z. We first compute the B field contributions
    Bxterm = cy_min * cy_max (Bx^R - Bx^L)/ (cy_min + cy_max)
    Byterm = cx_min * cx_max (By^R - By^L)/ (cx_min + cx_max)
  */
  const double B1term = vars->c2_min * vars->c2_max * (vars->B1r - vars->B1l) / c2_sum;
  const double B2term = vars->c1_min * vars->c1_max * (vars->B2r - vars->B2l) / c1_sum;

  /*
     Additionally, we need the electric field contribution
     E_z = -(v_x B_y - v_y B_x)
     For the 2D flux, that means we need to compute
     E^LL, E^LR, E^RL, and E^RR.
  */
  const double A3_rhs_rr = vars->v1rr*vars->B2r - vars->v2rr*vars->B1r;
  const double A3_rhs_rl = vars->v1rl*vars->B2r - vars->v2rl*vars->B1l;
  const double A3_rhs_lr = vars->v1lr*vars->B2l - vars->v2lr*vars->B1r;
  const double A3_rhs_ll = vars->v1ll*vars->B2l - vars->v2ll*vars->B1l;

  return B1term - B2term
         + ( vars->c2_max*vars->c1_max*A3_rhs_ll
           + vars->c2_min*vars->c1_max*A3_rhs_lr
           + vars->c2_max*vars->c1_min*A3_rhs_rl
           + vars->c2_min*vars->c1_min*A3_rhs_rr)
           /(c2_sum*c1_sum);
}
