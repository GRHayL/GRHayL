#include "ghl_induction.h"
#include "ghl_induction_helpers.h"

/**
 * @ingroup mag_gauge_internal
 * @brief Interpolates vector potential to all staggered grid points
 *
 * @details
 * \f$ A_i \f$ is staggered in the two directions perpendicular to \f$ i \f$
 * while \f$ \tilde{\Phi} \f$ is staggered in all directions. Thus,
 * we need to interpolate each component to the three other grid locations.
 * They are interpolated through averaging the two adjacent points. For example,
 * The function will average the values at the points \f$ i \f$ and
 * \f$ i+1 \f$ to get the value at \f$ i+\frac{1}{2} \f$.
 *
 * @param[in] Ax_stencil:     3D stencil array containing \f$ A_x \f$ from
 *                            \f$ (i-1, j-\frac{1}{2}, k-\frac{1}{2}) \f$ to
 *                            \f$ (i+1, j+\frac{3}{2}, k+\frac{3}{2}) \f$
 *
 * @param[in] Ay_stencil:     3D stencil array containing \f$ A_y \f$ from
 *                            \f$ (i-\frac{1}{2}, j-1, k-\frac{1}{2}) \f$ to
 *                            \f$ (i+\frac{3}{2}, j+1, k+\frac{3}{2}) \f$
 *
 * @param[in] Az_stencil:     3D stencil array containing \f$ A_z \f$ from
 *                            \f$ (i-\frac{1}{2}, j-\frac{1}{2}, k-1) \f$ to
 *                            \f$ (i+\frac{3}{2}, j+\frac{3}{2}, k+1) \f$
 *
 * @param[out] A_to_phitilde: \f$ A_i \f$ interpolated to
 *                            \f$ (i+\frac{1}{2}, j+\frac{1}{2}, k+\frac{1}{2}) \f$
 *
 * @param[out] A_to_Ax:       \f$ A_i \f$ interpolated to
 *                            \f$ (i, j+\frac{1}{2}, k+\frac{1}{2}) \f$
 *
 * @param[out] A_to_Ay:       \f$ A_i \f$ interpolated to
 *                            \f$ (i+\frac{1}{2}, j, k+\frac{1}{2}) \f$
 *
 * @param[out] A_to_Az:       \f$ A_i \f$ interpolated to
 *                            \f$ (i+\frac{1}{2}, j+\frac{1}{2}, k) \f$
 *
 * @returns void
 */
void ghl_A_i_avg(
      const double Ax_stencil[3][3][3],
      const double Ay_stencil[3][3][3],
      const double Az_stencil[3][3][3],
      double A_to_phitilde[3],
      double A_to_Ax[3],
      double A_to_Ay[3],
      double A_to_Az[3]) {

  A_to_phitilde[0] = 0.0;
  A_to_phitilde[1] = 0.0;
  A_to_phitilde[2] = 0.0;
  A_to_Ax[1]       = 0.0;
  A_to_Ax[2]       = 0.0;
  A_to_Ay[0]       = 0.0;
  A_to_Ay[2]       = 0.0;
  A_to_Az[0]       = 0.0;
  A_to_Az[1]       = 0.0;

  // A_i is already at the location of A_i :)
  A_to_Ax[0] = Ax_stencil[1][1][1];
  A_to_Ay[1] = Ay_stencil[1][1][1];
  A_to_Az[2] = Az_stencil[1][1][1];

  /*
     The stenciling sum here is reasonably simple. The index corresponding to
     the component needs to be interpolated to \f$ +\frac{1}{2} \f$ (the kk
     loop). Thus, the kk index is always in the same index as the component
     of \f$ A_i \f$ being interpolated. For the jj loop, the index corresponding
     to the interpolated quantity (i.e. x index for A_to_Ax) needs to be
     interpolated back by \f$ -\frac{1}{2} \f$. Hence, the kk loop is in the
     range [1,size) and the jj loop is in the range [0,size-1).
  */
  for(int kk=1; kk<3; kk++) {
    A_to_phitilde[0] += Ax_stencil[1][1][kk];
    A_to_phitilde[1] += Ay_stencil[1][kk][1];
    A_to_phitilde[2] += Az_stencil[kk][1][1];
    for(int jj=0; jj<2; jj++) {
      A_to_Ax[1] += Ay_stencil[1][kk][jj];
      A_to_Ax[2] += Az_stencil[kk][1][jj];

      A_to_Ay[0] += Ax_stencil[1][jj][kk];
      A_to_Ay[2] += Az_stencil[kk][jj][1];

      A_to_Az[0] += Ax_stencil[jj][1][kk];
      A_to_Az[1] += Ay_stencil[jj][kk][1];
    }
  }
  A_to_phitilde[0] /= 2.0;
  A_to_phitilde[1] /= 2.0;
  A_to_phitilde[2] /= 2.0;
  A_to_Ax[1]       /= 4.0;
  A_to_Ax[2]       /= 4.0;
  A_to_Ay[0]       /= 4.0;
  A_to_Ay[2]       /= 4.0;
  A_to_Az[0]       /= 4.0;
  A_to_Az[1]       /= 4.0;
}

/**
 * @ingroup mag_gauge_internal
 * @brief Interpolates cell-centered BSSN quantities to staggered grid points
 *
 * @details
 * This function interpolates the cell-centered BSSN quantities to the
 * following grid points via averaging:
 *
 * - \f$ \alpha \f$, \f$ \beta^i \f$, and \f$ \frac{\alpha}{\psi^6} \f$ to
 *   - \f$ (i+\frac{1}{2}, j+\frac{1}{2}, k+\frac{1}{2}) \f$
 *
 * - \f$ \gamma^{ij} \f$ and \f$ \frac{\alpha}{\psi^2} \f$ to
 *   - \f$ (i,             j+\frac{1}{2}, k+\frac{1}{2}) \f$
 *   - \f$ (i+\frac{1}{2}, j,             k+\frac{1}{2}) \f$
 *   - \f$ (i+\frac{1}{2}, j+\frac{1}{2}, k) \f$
 *
 * Note that the last quantity is returned in ghl_metric_quantities::gammaUU
 * with the first index corresponding to the non-staggered component.
 *
 * @param[in] metric_stencil: 3D stencil array of ghl_metric_quantities from
 *                            \f$ (i, j, k) \f$ to \f$ (i+1, j+1, k+1) \f$
 *
 * @param[in] psi_stencil:    3D stencil array of \f$ \psi^6 \f$ from
 *                            \f$ (i, j, k) \f$ to \f$ (i+1, j+1, k+1) \f$
 *
 * @param[out] metric_interp: ghl_metric_quantities containing interpolated
 *                            quantities \f$ \alpha \f$, \f$ \beta^i \f$,
 *                            and \f$ \gamma^{ij} \f$
 *
 * @param[out] lapse_psi2_interp: 1D array containing interpolated values
 *                                for \f$ \frac{\alpha}{\psi^2} \f$
 *
 * @param[out] lapse_over_psi6_interp: interpolated value of \f$ \frac{\alpha}{\psi^6} \f$
 *
 * @returns void
 */
void ghl_BSSN_cell_interp(
      const ghl_metric_quantities metric_stencil[2][2][2],
      const double psi_stencil[2][2][2],
      ghl_metric_quantities *restrict metric_interp,
      double lapse_psi2_interp[3],
      double *restrict lapse_over_psi6_interp) {

  metric_interp->lapse         = 0.0;
  metric_interp->betaU[0]      = 0.0;
  metric_interp->betaU[1]      = 0.0;
  metric_interp->betaU[2]      = 0.0;
  metric_interp->gammaUU[0][0] = 0.0;
  metric_interp->gammaUU[0][1] = 0.0;
  metric_interp->gammaUU[0][2] = 0.0;
  metric_interp->gammaUU[1][0] = 0.0;
  metric_interp->gammaUU[1][1] = 0.0;
  metric_interp->gammaUU[1][2] = 0.0;
  metric_interp->gammaUU[2][0] = 0.0;
  metric_interp->gammaUU[2][1] = 0.0;
  metric_interp->gammaUU[2][2] = 0.0;
  lapse_psi2_interp[0]         = 0.0;
  lapse_psi2_interp[1]         = 0.0;
  lapse_psi2_interp[2]         = 0.0;
  *lapse_over_psi6_interp      = 0.0;

  double lapse_psi2[2][2][2];
  double lapse_over_psi6[2][2][2];
  for(int kk=0; kk<2; kk++) {
    for(int jj=0; jj<2; jj++) {
      for(int ii=0; ii<2; ii++) {
        const double psi2 = psi_stencil[kk][jj][ii]*psi_stencil[kk][jj][ii];
        lapse_psi2[kk][jj][ii] = metric_stencil[kk][jj][ii].lapse*psi2;
        lapse_over_psi6[kk][jj][ii] = metric_stencil[kk][jj][ii].lapse/(psi2*psi2*psi2);
      }
    }
  }

  for(int kk=0; kk<2; kk++) {
    for(int jj=0; jj<2; jj++) {
      // Interpolate xx, xy, xz from ccc to cvv centering for A_x
      metric_interp->gammaUU[0][0] += metric_stencil[kk][jj][0].gammaUU[0][0];
      metric_interp->gammaUU[0][1] += metric_stencil[kk][jj][0].gammaUU[0][1];
      metric_interp->gammaUU[0][2] += metric_stencil[kk][jj][0].gammaUU[0][2];
      lapse_psi2_interp[0]         += lapse_psi2[kk][jj][0];

      // Interpolate yx, yy, yz from ccc to vcv centering for A_y
      metric_interp->gammaUU[1][0] += metric_stencil[kk][0][jj].gammaUU[0][1];
      metric_interp->gammaUU[1][1] += metric_stencil[kk][0][jj].gammaUU[1][1];
      metric_interp->gammaUU[1][2] += metric_stencil[kk][0][jj].gammaUU[1][2];
      lapse_psi2_interp[1]         += lapse_psi2[kk][0][jj];

      // Interpolate zx, zy, zz from ccc to vvc centering for A_z
      metric_interp->gammaUU[2][0] += metric_stencil[0][kk][jj].gammaUU[0][2];
      metric_interp->gammaUU[2][1] += metric_stencil[0][kk][jj].gammaUU[1][2];
      metric_interp->gammaUU[2][2] += metric_stencil[0][kk][jj].gammaUU[2][2];
      lapse_psi2_interp[2]         += lapse_psi2[0][kk][jj];

      for(int ii=0; ii<2; ii++) {
        metric_interp->lapse    += metric_stencil[kk][jj][ii].lapse;
        metric_interp->betaU[0] += metric_stencil[kk][jj][ii].betaU[0];
        metric_interp->betaU[1] += metric_stencil[kk][jj][ii].betaU[1];
        metric_interp->betaU[2] += metric_stencil[kk][jj][ii].betaU[2];
        *lapse_over_psi6_interp += lapse_over_psi6[kk][jj][ii];
      }
    }
  }
  metric_interp->lapse    /= 8.0;
  metric_interp->betaU[0] /= 8.0;
  metric_interp->betaU[1] /= 8.0;
  metric_interp->betaU[2] /= 8.0;
  *lapse_over_psi6_interp /= 8.0;

  metric_interp->gammaUU[0][0] /= 4.0;
  metric_interp->gammaUU[0][1] /= 4.0;
  metric_interp->gammaUU[0][2] /= 4.0;
  metric_interp->gammaUU[1][0] /= 4.0;
  metric_interp->gammaUU[1][1] /= 4.0;
  metric_interp->gammaUU[1][2] /= 4.0;
  metric_interp->gammaUU[2][0] /= 4.0;
  metric_interp->gammaUU[2][1] /= 4.0;
  metric_interp->gammaUU[2][2] /= 4.0;
  lapse_psi2_interp[0]         /= 4.0;
  lapse_psi2_interp[1]         /= 4.0;
  lapse_psi2_interp[2]         /= 4.0;
}

/**
 * @ingroup mag_gauge_internal
 * @brief Interpolates cell-centered ADM quantities to staggered grid points
 *
 * @details
 * This function interpolates the cell-centered ADM quantities to the
 * following grid points via averaging:
 *
 * - \f$ \alpha \f$, \f$ \beta^i \f$, and \f$ \frac{\alpha}{\psi^6} \f$ to
 *   - \f$ (i+\frac{1}{2}, j+\frac{1}{2}, k+\frac{1}{2}) \f$
 *
 * - \f$ \alpha\sqrt{\gamma}\gamma^{ij} \f$ to
 *   - \f$ (i,             j+\frac{1}{2}, k+\frac{1}{2}) \f$
 *   - \f$ (i+\frac{1}{2}, j,             k+\frac{1}{2}) \f$
 *   - \f$ (i+\frac{1}{2}, j+\frac{1}{2}, k) \f$
 *
 * Note that the last quantity is returned in ghl_metric_quantities::gammaUU
 * with the first index corresponding to the non-staggered component.
 *
 * @param[in] metric_stencil: 3D stencil array of ghl_metric_quantities from
 *                            \f$ (i, j, k) \f$ to \f$ (i+1, j+1, k+1) \f$
 *
 * @param[out] metric_interp: ghl_metric_quantities containing interpolated
 *                            quantities \f$ \alpha \f$, \f$ \beta^i \f$,
 *                            and \f$ \alpha\sqrt{\gamma}\gamma^{ij} \f$
 *
 * @param[out] lapse_over_psi6_interp: interpolated value of \f$ \frac{\alpha}{\psi^6} \f$
 *
 * @returns void
 */
void ghl_ADM_cell_interp(
      const ghl_metric_quantities metric_stencil[2][2][2],
      ghl_metric_quantities *restrict metric_interp,
      double *restrict lapse_over_psi6_interp) {

  metric_interp->lapse         = 0.0;
  metric_interp->betaU[0]      = 0.0;
  metric_interp->betaU[1]      = 0.0;
  metric_interp->betaU[2]      = 0.0;
  metric_interp->gammaUU[0][0] = 0.0;
  metric_interp->gammaUU[0][1] = 0.0;
  metric_interp->gammaUU[0][2] = 0.0;
  metric_interp->gammaUU[1][0] = 0.0;
  metric_interp->gammaUU[1][1] = 0.0;
  metric_interp->gammaUU[1][2] = 0.0;
  metric_interp->gammaUU[2][0] = 0.0;
  metric_interp->gammaUU[2][1] = 0.0;
  metric_interp->gammaUU[2][2] = 0.0;
  *lapse_over_psi6_interp      = 0.0;

  for(int kk=0; kk<2; kk++) {
    for(int jj=0; jj<2; jj++) {
      // Interpolate xx, xy, xz from ccc to cvv centering for A_x
      const double detgx = metric_stencil[kk][jj][0].lapse*metric_stencil[kk][jj][0].sqrt_detgamma;
      metric_interp->gammaUU[0][0] += detgx*metric_stencil[kk][jj][0].gammaUU[0][0];
      metric_interp->gammaUU[0][1] += detgx*metric_stencil[kk][jj][0].gammaUU[0][1];
      metric_interp->gammaUU[0][2] += detgx*metric_stencil[kk][jj][0].gammaUU[0][2];

      // Interpolate yx, yy, yz from ccc to vcv centering for A_y
      const double detgy = metric_stencil[kk][0][jj].lapse*metric_stencil[kk][0][jj].sqrt_detgamma;
      metric_interp->gammaUU[1][0] += detgy*metric_stencil[kk][0][jj].gammaUU[0][1];
      metric_interp->gammaUU[1][1] += detgy*metric_stencil[kk][0][jj].gammaUU[1][1];
      metric_interp->gammaUU[1][2] += detgy*metric_stencil[kk][0][jj].gammaUU[1][2];

      // Interpolate zx, zy, zz from ccc to vvc centering for A_z
      const double detgz = metric_stencil[0][kk][jj].lapse*metric_stencil[0][kk][jj].sqrt_detgamma;
      metric_interp->gammaUU[2][0] += detgz*metric_stencil[0][kk][jj].gammaUU[0][2];
      metric_interp->gammaUU[2][1] += detgz*metric_stencil[0][kk][jj].gammaUU[1][2];
      metric_interp->gammaUU[2][2] += detgz*metric_stencil[0][kk][jj].gammaUU[2][2];

      for(int ii=0; ii<2; ii++) {
        metric_interp->lapse    += metric_stencil[kk][jj][ii].lapse;
        metric_interp->betaU[0] += metric_stencil[kk][jj][ii].betaU[0];
        metric_interp->betaU[1] += metric_stencil[kk][jj][ii].betaU[1];
        metric_interp->betaU[2] += metric_stencil[kk][jj][ii].betaU[2];
        *lapse_over_psi6_interp += metric_stencil[kk][jj][ii].lapse/metric_stencil[kk][jj][ii].sqrt_detgamma;
      }
    }
  }
  metric_interp->lapse    /= 8.0;
  metric_interp->betaU[0] /= 8.0;
  metric_interp->betaU[1] /= 8.0;
  metric_interp->betaU[2] /= 8.0;
  *lapse_over_psi6_interp /= 8.0;

  metric_interp->gammaUU[0][0] /= 4.0;
  metric_interp->gammaUU[0][1] /= 4.0;
  metric_interp->gammaUU[0][2] /= 4.0;
  metric_interp->gammaUU[1][0] /= 4.0;
  metric_interp->gammaUU[1][1] /= 4.0;
  metric_interp->gammaUU[1][2] /= 4.0;
  metric_interp->gammaUU[2][0] /= 4.0;
  metric_interp->gammaUU[2][1] /= 4.0;
  metric_interp->gammaUU[2][2] /= 4.0;
}

/**
 * @ingroup mag_gauge_internal
 * @brief Interpolates vertex-centered ADM quantities to staggered grid points
 *
 * @details
 * Since the spacetime quantities are already fully staggered, this function
 * only interpolates
 *
 * - \f$ \alpha\sqrt{\gamma}\gamma^{ij} \f$ to
 *   - \f$ (i,             j+\frac{1}{2}, k+\frac{1}{2}) \f$
 *   - \f$ (i+\frac{1}{2}, j,             k+\frac{1}{2}) \f$
 *   - \f$ (i+\frac{1}{2}, j+\frac{1}{2}, k) \f$
 *
 * Note that this quantity is returned in ghl_metric_quantities::gammaUU
 * with the first index corresponding to the non-staggered component.
 *
 * @param[in] metric_stencil: 3D stencil array of ghl_metric_quantities from
 *                            \f$ (i-\frac{1}{2}, j-\frac{1}{2}, k-\frac{1}{2}) \f$ to
 *                            \f$ (i+\frac{3}{2}, j+\frac{3}{2}, k+\frac{3}{2}) \f$
 *
 * @param[out] gammaUU_interp: 2D array containing interpolated
 *                            quantity \f$ \alpha\sqrt{\gamma}\gamma^{ij} \f$
 *
 * @returns void
 */
void ghl_ADM_vertex_interp(
      const ghl_metric_quantities metric_stencil[2][2][2],
      double gammaUU_interp[3][3]) {

  gammaUU_interp[0][0] = 0.0;
  gammaUU_interp[0][1] = 0.0;
  gammaUU_interp[0][2] = 0.0;
  gammaUU_interp[1][0] = 0.0;
  gammaUU_interp[1][1] = 0.0;
  gammaUU_interp[1][2] = 0.0;
  gammaUU_interp[2][0] = 0.0;
  gammaUU_interp[2][1] = 0.0;
  gammaUU_interp[2][2] = 0.0;

  for(int ii=0; ii<2; ii++) {
    // Interpolate xx, xy, xz from vvv to cvv centering for A_x
    const double detgx = metric_stencil[0][0][ii].lapse*metric_stencil[0][0][ii].sqrt_detgamma;
    gammaUU_interp[0][0] += detgx*metric_stencil[0][0][ii].gammaUU[0][0];
    gammaUU_interp[0][1] += detgx*metric_stencil[0][0][ii].gammaUU[0][1];
    gammaUU_interp[0][2] += detgx*metric_stencil[0][0][ii].gammaUU[0][2];

    // Interpolate yx, yy, yz from vvv to vcv centering for A_y
    const double detgy = metric_stencil[0][ii][0].lapse*metric_stencil[0][ii][0].sqrt_detgamma;
    gammaUU_interp[1][0] += detgy*metric_stencil[0][ii][0].gammaUU[0][1];
    gammaUU_interp[1][1] += detgy*metric_stencil[0][ii][0].gammaUU[1][1];
    gammaUU_interp[1][2] += detgy*metric_stencil[0][ii][0].gammaUU[1][2];

    // Interpolate zx, zy, zz from vvv to vvc centering for A_z
    const double detgz = metric_stencil[ii][0][0].lapse*metric_stencil[ii][0][0].sqrt_detgamma;
    gammaUU_interp[2][0] += detgz*metric_stencil[ii][0][0].gammaUU[0][2];
    gammaUU_interp[2][1] += detgz*metric_stencil[ii][0][0].gammaUU[1][2];
    gammaUU_interp[2][2] += detgz*metric_stencil[ii][0][0].gammaUU[2][2];
  }
  gammaUU_interp[0][0] /= 2.0;
  gammaUU_interp[0][1] /= 2.0;
  gammaUU_interp[0][2] /= 2.0;
  gammaUU_interp[1][0] /= 2.0;
  gammaUU_interp[1][1] /= 2.0;
  gammaUU_interp[1][2] /= 2.0;
  gammaUU_interp[2][0] /= 2.0;
  gammaUU_interp[2][1] /= 2.0;
  gammaUU_interp[2][2] /= 2.0;
}
