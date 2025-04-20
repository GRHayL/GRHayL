#ifndef GHL_INDUCTION_H_
#define GHL_INDUCTION_H_

#include "ghl.h"

/**
 * @ingroup mag_flux
 * @struct ghl_HLL_vars
 * @brief Stores all quantities needed to compute the magnetic flux term.
 *
 * @todo
 * References @ref ghl_calculate_characteristic_speed instead of specific flavors
 *
 * @details
 * This struct is used by the @ref Induction to pass needed quantities
 * into the @ref mag_flux functions that compute the \f$ A_k^\mathrm{RHS} \f$
 * term. As these functions are performing a cross-product, the variables
 * are not immediately obvious. For a given component \f$ A_k \f$, the input
 * components are
 *
 * - \f$ k=z \f$: \f$ 1=x \f$, \f$ 2=y \f$
 * - \f$ k=y \f$: \f$ 1=z \f$, \f$ 2=x \f$
 * - \f$ k=x \f$: \f$ 1=y \f$, \f$ 2=z \f$
 *
 * The numbering in the variable names follows the
 * [common convention](https://en.wikipedia.org/wiki/Cross_product)
 * where e.g. the z component of \f$ a \times b \f$ is
 *
 * \f[
 * a_1 b_2 - a_2 b_1
 * \f]
 *
 * Of course, this notation is inconvenient for vectors as it can be confused
 * with exponents. We use e.g. B1 instead of \f$ B^1 \f$ for clarity. The
 * variable reconstruction procedure to get the correct l/r suffixes is
 * described in more detail in @ref ghl_HLL_flux_with_Btilde .
 */
typedef struct ghl_HLL_vars {
  /** Value of B1 on the right side */
  double B1r;
  /** Value of B1 on the left side */
  double B1l;
  /** Value of B2 on the right side */
  double B2r;
  /** Value of B2 on the left side */
  double B2l;
  /** Minimal characteristic speed as computed by @ref
      ghl_calculate_characteristic_speed in the '1' direction */
  double c1_min;
  /** Maximal characteristic speed as computed by @ref
      ghl_calculate_characteristic_speed in the '1' direction */
  double c1_max;
  /** Minimal characteristic speed as computed by @ref
      ghl_calculate_characteristic_speed in the '2' direction */
  double c2_min;
  /** Maximal characteristic speed as computed by @ref
      ghl_calculate_characteristic_speed in the '2' direction */
  double c2_max;
  /** Value of v1 on the right side of both reconstructions */
  double v1rr;
  /** Value of v1 on the right (left) side of the first (second) reconstruction */
  double v1rl;
  /** Value of v1 on the left (right) side of the first (second) reconstruction */
  double v1lr;
  /** Value of v1 on the left side of both reconstructions */
  double v1ll;
  /** Value of v2 on the right side of both reconstructions */
  double v2rr;
  /** Value of v2 on the right (left) side of the first (second) reconstruction */
  double v2rl;
  /** Value of v2 on the left (right) side of the first (second) reconstruction */
  double v2lr;
  /** Value of v2 on the left side of both reconstructions */
  double v2ll;
} ghl_HLL_vars;

/**
 * @ingroup mag_gauge
 * @struct ghl_induction_interp_vars
 * @brief Stores all interpolated quantities for induction functions.
 *
 * @details
 * This struct is used by the @ref Induction interpolators to store
 * the interpolated quantities, which can then be unpacked by the
 * user as desired.
 */
typedef struct ghl_induction_interp_vars {
  /** Interpolated lapse \f$ \alpha \f$ */
  double alpha;
  /** Interpolated shift \f$ \beta^i \f$ */
  double betai[3];
  /** Interpolated quantity \f$ \alpha \Phi - \beta^j A_j \f$ */
  double alpha_Phi_minus_betaj_A_j;
  /** Interpolated quantity \f$ \sqrt{g} A^i \f$ */
  double sqrtg_Ai[3];
} ghl_induction_interp_vars;

//------------------ Functions ---------------------

#ifdef __cplusplus
extern "C" {
#endif

double ghl_HLL_flux_with_B(
      const double psi6,
      const ghl_HLL_vars *restrict vars);

double ghl_HLL_flux_with_Btilde(
      const ghl_HLL_vars *restrict vars);

void ghl_interpolate_with_cell_centered_ADM(
      const ghl_metric_quantities metric_stencil[2][2][2],
      const double Ax_stencil[3][3][3],
      const double Ay_stencil[3][3][3],
      const double Az_stencil[3][3][3],
      const double phitilde,
      ghl_induction_interp_vars *restrict interp_vars);

void ghl_interpolate_with_cell_centered_BSSN(
      const ghl_metric_quantities metric[2][2][2],
      const double psi_stencil[2][2][2],
      const double Ax_stencil[3][3][3],
      const double Ay_stencil[3][3][3],
      const double Az_stencil[3][3][3],
      const double phitilde,
      ghl_induction_interp_vars *restrict interp_vars);

void ghl_interpolate_with_vertex_centered_ADM(
      const ghl_metric_quantities metric_stencil[2][2][2],
      const double Ax_stencil[3][3][3],
      const double Ay_stencil[3][3][3],
      const double Az_stencil[3][3][3],
      const double phitilde,
      ghl_induction_interp_vars *restrict interp_vars);

double ghl_calculate_phitilde_rhs(
      const double dxi[3],
      const double Lorenz_damping_factor,
      const double alpha_interp,
      const double shiftx_interp[5],
      const double shifty_interp[5],
      const double shiftz_interp[5],
      const double sqrtg_Ai_stencil[3][2],
      const double phitilde_stencil[3][5]);

#ifdef __cplusplus
}
#endif

//--------------------------------------------------

#endif // GHL_INDUCTION_H_
