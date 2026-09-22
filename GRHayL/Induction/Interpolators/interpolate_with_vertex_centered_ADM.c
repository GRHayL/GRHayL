#include "ghl_induction.h"
#include "ghl_induction_helpers.h"

/**
 * @ingroup mag_gauge
 * @brief Interpolate induction gauge quantities using a vertex-centered ADM
 *        metric.
 *
 * @details
 * C array order is `[z][y][x]`. The caller places the vertex at the
 * \f$ \tilde{\Phi} \f$ point in `metric_stencil[1][1][1]` and its backward
 * x, y, and z neighbors in `[1][1][0]`, `[1][0][1]`, and `[0][1][1]`.
 * These pairs collocate \f$\alpha\sqrt{\gamma}\gamma^{ij}\f$ with the
 * corresponding staggered \f$A_i\f$ values. The other four metric entries
 * are accepted as part of the fixed public array shape but are not read.
 *
 * Each used metric element must contain initialized `lapse`, `betaU`,
 * `gammaUU`, and `sqrt_detgamma` fields, as produced by
 * @ref ghl_initialize_metric.
 * The vector-potential arrays use the same `[z][y][x]` order and the
 * 3-point ranges documented below.
 *
 * This function assigns only `sqrtg_Ai` and
 * `alpha_Phi_minus_betaj_A_j`. It deliberately leaves `alpha` and `betai`
 * unchanged because lapse and shift already live at the fully staggered
 * vertex. It performs no bounds or centering checks; packing is the caller's
 * responsibility.
 *
 * @param[in] metric_stencil Vertex-centered ADM stencil spanning one vertex
 *                           backward through the current vertex on each axis,
 *                           with the current vertex at `[1][1][1]`.
 * @param[in] Ax_stencil \f$A_x\f$ from
 *                       \f$(i-1,j-\frac12,k-\frac12)\f$ through
 *                       \f$(i+1,j+\frac32,k+\frac32)\f$.
 * @param[in] Ay_stencil \f$A_y\f$ from
 *                       \f$(i-\frac12,j-1,k-\frac12)\f$ through
 *                       \f$(i+\frac32,j+1,k+\frac32)\f$.
 * @param[in] Az_stencil \f$A_z\f$ from
 *                       \f$(i-\frac12,j-\frac12,k-1)\f$ through
 *                       \f$(i+\frac32,j+\frac32,k+1)\f$.
 * @param[in] phitilde \f$\tilde{\Phi}\f$ at the current vertex.
 * @param[out] interp_vars Output whose `sqrtg_Ai` and
 *                         `alpha_Phi_minus_betaj_A_j` fields are assigned.
 */
void ghl_interpolate_with_vertex_centered_ADM(
      const ghl_metric_quantities metric_stencil[2][2][2],
      const double Ax_stencil[3][3][3],
      const double Ay_stencil[3][3][3],
      const double Az_stencil[3][3][3],
      const double phitilde,
      ghl_induction_interp_vars *restrict interp_vars) {
  /*
     We need to interpolate several quantities to several different points depending on
     the quantities we're computing. The staggered gridpoints for these variables are
       phitilde: (i+1/2, j+1/2, k+1/2)
       A_x:      (i,     j+1/2, k+1/2)
       A_y:      (i+1/2, j,     k+1/2)
       A_z:      (i+1/2, j+1/2, k    )
     For metric quantities, we use ghl_ADM_vertex_interp(), which computes most of the
     needed quantities. It interpolates the metric to 3 different points:
       gammaUU[0][i] is at A_x's location
       gammaUU[1][i] is at A_y's location
       gammaUU[2][i] is at A_z's location
     Note that we actually store detg*gammaUU to reduce the memory usage.
     Similarly, the function ghl_A_i_avg() interpolates A_i to these points, storing the
     interpolated data in the 4 arrays A_to_phitilde, A_to_Ax, A_to_Ay, and A_to_Az.

     These two averaging loops are split because the stencils are of different sizes. The
     metric has an even stencil, and the A_i have an odd stencil.
  */
  double gammaUU_interp[3][3];
  ghl_ADM_vertex_interp(metric_stencil, gammaUU_interp);

  double A_to_phitilde[3], A_to_Ax[3], A_to_Ay[3], A_to_Az[3];
  ghl_A_i_avg(Ax_stencil, Ay_stencil, Az_stencil, A_to_phitilde, A_to_Ax, A_to_Ay, A_to_Az);

  // Compute \partial_t phitilde = -\partial_i (  \alpha psi^6 A^i - phitilde \beta^i)
  // A^x term (interpolated to (i, j+1/2, k+1/2) )
  // \sqrt{-g} A^x = \alpha \sqrt{\gamma} A^x (RHS of \partial_i psi6phi)
  interp_vars->sqrtg_Ai[0] = gammaUU_interp[0][0]*A_to_Ax[0]
                           + gammaUU_interp[0][1]*A_to_Ax[1]
                           + gammaUU_interp[0][2]*A_to_Ax[2];

  // A^y term (interpolated to (i+1/2, j, k+1/2) )
  // \sqrt{-g} A^y = \alpha \sqrt{\gamma} A^y (RHS of \partial_i psi6phi)
  interp_vars->sqrtg_Ai[1] = gammaUU_interp[1][0]*A_to_Ay[0]
                           + gammaUU_interp[1][1]*A_to_Ay[1]
                           + gammaUU_interp[1][2]*A_to_Ay[2];

  // A^z term (interpolated to (i+1/2, j+1/2, k) )
  // \sqrt{-g} A^z = \alpha \sqrt{\gamma} A^z (RHS of \partial_i psi6phi)
  interp_vars->sqrtg_Ai[2] = gammaUU_interp[2][0]*A_to_Az[0]
                           + gammaUU_interp[2][1]*A_to_Az[1]
                           + gammaUU_interp[2][2]*A_to_Az[2];

  const ghl_metric_quantities *restrict metric_at_phitilde
        = &metric_stencil[1][1][1];

  // Next set \alpha \Phi - \beta^j A_j at (i+1/2,j+1/2,k+1/2)
  // \alpha \Phi = \alpha \tilde{\Phi} / psi^6
  //             = \alpha \tilde{\Phi} / \sqrt{\gamma}
  interp_vars->alpha_Phi_minus_betaj_A_j
        = phitilde * metric_at_phitilde->lapse / metric_at_phitilde->sqrt_detgamma
          - (metric_at_phitilde->betaU[0] * A_to_phitilde[0]
             + metric_at_phitilde->betaU[1] * A_to_phitilde[1]
             + metric_at_phitilde->betaU[2] * A_to_phitilde[2]);
}
