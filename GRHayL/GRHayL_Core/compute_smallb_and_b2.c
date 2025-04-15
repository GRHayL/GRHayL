#include "ghl_con2prim.h"

/**
 * @ingroup GRHayL_Core
 * @brief Computes magnetic quantities \f$ b^0 \f$, \f$ b^i \f$, and \f$ b^2 \f$.
 *
 * @details
 * This function computes the quantity \f$ b^\mu \f$ and its square using equations from
 * [Duez _et al_](https://arxiv.org/abs/astro-ph/0503420v2). Note that GRHayL's definition
 * of \f$ B \f$ differs by a factor of \f$ \frac{1}{\sqrt{4\pi}} \f$ with respect to that
 * paper. The 0th component is computed using Eqs. 23 and 31:
 *
 * \f[
 * b^t = \frac{u_i B^i}{\alpha}
 * \f]
 *
 * The spatial components are computed using Eq. 24:
 *
 * \f[
 * \begin{align}
 * b^i &= B^i_u \\
 *     &= \frac{1}{u^0} \left( B^i/\alpha + B^0_u u^i \right) \\
 *     &= \frac{1}{u^0} \left( B^i/\alpha + b^t u^i \right) \\
 *     &= \frac{1}{\alpha u^0} \left( B^i + b^t u^i \right) \\
 *     &= \frac{1}{\alpha} \left( B^i/u^0 + b^t u^i/u^0 \right) \\
 *     &= \frac{1}{\alpha} \left( B^i/u^0 + b^t v^i \right) \\
 *     &= \frac{B^i}{\alpha u^0} + b^t v^i
 * \end{align}
 * \f]
 *
 * Then, we compute \f$ b^2 \f$ with
 *
 * \f[
 * \begin{align}
 * b^2 &= g_{\mu\nu} \\
 *     &= g_{tt}\left(b^t\right)^2 + g_{xx}\left(b^x\right)^2 + g_{yy}\left(b^y\right)^2 + g_{zz}\left(b^z\right)^2 + 2\left(g_{ti} b^t b^i + g_{xy} b^x b^y + g_{xz} b^x b^z + g_{yz} b^y b^z\right) \\
 *     &= \left( -\alpha^2 + \gamma_{i j} \beta^i \beta^j \right) \left(b^t\right)^2 + \gamma_{i j} b^i b^j + 2 g_{t i} b^t b^i \\
 *     &= -\left( \alpha b^t \right)^2 + \gamma_{i j} \beta^i \beta^j \left(b^t\right)^2 + \gamma_{i j} b^i b^j + 2 b^t (\gamma_{i j}\beta^j) b^i \\
 *     &= -\left( \alpha b^t \right)^2 + \gamma_{i j} \left( \beta^i \beta^j \left(b^t\right)^2 + b^i b^j + 2 b^t \beta^j b^i \right) \\
 *     &= -\left( \alpha b^t \right)^2 + \gamma_{i j} \left(b^i + b^t\beta^i\right) \left(b^j + b^t\beta^j\right)
 * \end{align}
 * \f]
 *
 * @param[in] ADM_metric: pointer to a ghl_metric_quantities struct containing the ADM metric
 *
 * @param[in] prims:      pointer to a ghl_primitive_quantities struct. Note that
 *                        ghl_primitive_quantities::u0 **must** be valid for this function.
 *
 * @param[in] uDN:        lowered 4-velocity \f$ u_\mu \f$
 *
 * @param[out] smallb:    returned 4-vector \f$ b^\mu \f$
 *
 * @param[out] smallb2:   returned scalar \f$ b^2 \f$
 *
 * @returns void
 */
void ghl_compute_smallb_and_b2(
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_primitive_quantities *restrict prims,
      const double uDN[4],
      double smallb[4],
      double *restrict smallb2) {

  const double one_over_W = ADM_metric->lapseinv/prims->u0;

  /*
     Compute b^t (Eqs. 23 and 31 in http://arxiv.org/pdf/astro-ph/0503420.pdf):
       b^t = u_i B^i/alpha sqrt(4 pi)
     Note that our B = B/sqrt(4pi) of the paper's B. Then,
       b^t = u_i B^i/alpha
  */
  smallb[0] = (uDN[1]*prims->BU[0] + uDN[2]*prims->BU[1] + uDN[3]*prims->BU[2])*ADM_metric->lapseinv;

  /*
     Eq. 24 in http://arxiv.org/pdf/astro-ph/0503420.pdf
       b^i = B^i_u / sqrt(4 pi)
       b^i = ( B^i/alpha + B^0_u u^i ) / ( u^0 sqrt(4 pi) )
       b^i = ( B^i/alpha +  sqrt(4 pi) b^t u^i ) / ( u^0 sqrt(4 pi) )
       b^i = ( B^i +  alpha sqrt(4 pi) b^t u^i ) / ( alpha u^0 sqrt(4 pi) )
       b^i = ( B^i/u^0 +  alpha sqrt(4 pi) b^t u^i/u^0 ) / ( alpha sqrt(4 pi) )
       b^i = ( B^i/u^0 +  alpha sqrt(4 pi) b^t v^i ) / ( alpha sqrt(4 pi) )

     Note that our B = B/sqrt(4pi) of the paper's B. Then,
       b^i = ( B^i/u^0 +  alpha b^t v^i ) / alpha
       b^i = B^i/W + b^t v^i
  */
  smallb[1] = prims->BU[0]*one_over_W + prims->vU[0]*smallb[0];
  smallb[2] = prims->BU[1]*one_over_W + prims->vU[1]*smallb[0];
  smallb[3] = prims->BU[2]*one_over_W + prims->vU[2]*smallb[0];

  // b^2 = g_{\mu \nu} b^{\mu} b^{\nu}
  //     = gtt bt^2 + gxx bx^2 + gyy by^2 + gzz bz^2 + 2 (gtx bt bx + gty bt by + gtz bt bz + gxy bx by + gxz bx bz + gyz by bz)
  //     = (-al^2 + gamma_{ij} betai betaj) bt^2 + b^i b^j gamma_{ij} + 2 g_{t i} b^t b^i
  //     = - (alpha b^t)^2 + (b^t)^2 gamma_{ij} beta^i beta^j + b^i b^j gamma_{ij} + 2 b^t g_{t i} b^i
  //     = - (alpha b^t)^2 + (b^t)^2 gamma_{ij} beta^i beta^j + b^i b^j gamma_{ij} + 2 b^t (gamma_{ij} beta^j) b^i
  //     = - (alpha b^t)^2 + gamma_{ij} ((b^t)^2 beta^i beta^j + b^i b^j + 2 b^t beta^j b^i)
  //     = - (alpha b^t)^2 + gamma_{ij} ((b^t)^2 beta^i beta^j + 2 b^t beta^j b^i + b^i b^j)
  //     = - (alpha b^t)^2 + gamma_{ij} (b^i + b^t beta^i) (b^j + b^t beta^j)
  const double bi_plus_bt_betai[3] = {smallb[1] + smallb[0]*ADM_metric->betaU[0],
                                      smallb[2] + smallb[0]*ADM_metric->betaU[1],
                                      smallb[3] + smallb[0]*ADM_metric->betaU[2]};

  *smallb2 = -SQR(ADM_metric->lapse*smallb[0]) + ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, bi_plus_bt_betai);
}
