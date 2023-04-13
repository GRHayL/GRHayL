#include "con2prim.h"

/* Function    : compute_smallb_and_b2()
 * Description : Computes magnetic quantities b^0, b^i and b^2 (see Eqs. 23
                 and 24 in https://arxiv.org/abs/astro-ph/0503420).
 *
 * Inputs      : params         - GRHayL_parameters struct with parameters
 *                                for the simulation
 *             : metric         - metric_quantities struct with data for
 *                                the gridpoint of interest
 *             : prims          - primitive_quantities struct with data
 *                                for the gridpoint of interest
 *             : uDN[4]         - lowered 4-velocity u_\mu
 *
 * Outputs     : smallb         - returns computed b^\mu array
 *             : smallb2        - returns computed b^2
 *
 */

void compute_smallb_and_b2(
      const metric_quantities *restrict metric,
      const primitive_quantities *restrict prims,
      const double uDN[4],
      double *restrict smallb,
      double *restrict smallb2) {

  const double ONE_OVER_LAPSE_SQRT_4PI = metric->lapseinv*ONE_OVER_SQRT_4PI;
  const double ONE_OVER_U0 = 1.0/prims->u0;

  // Eqs. 23 and 31 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
  //   Compute alpha sqrt(4 pi) b^t = u_i B^i
  const double alpha_sqrt_4pi_bt = uDN[1]*prims->Bx + uDN[2]*prims->By + uDN[3]*prims->Bz;

  // Eq. 24 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
  // b^i = B^i_u / sqrt(4 pi)
  // b^i = ( B^i/alpha + B^0_u u^i ) / ( u^0 sqrt(4 pi) )
  // b^i = ( B^i/alpha +  sqrt(4 pi) b^t u^i ) / ( u^0 sqrt(4 pi) )
  // b^i = ( B^i +  alpha sqrt(4 pi) b^t u^i ) / ( alpha u^0 sqrt(4 pi) )
  // b^i = ( B^i/u^0 +  alpha sqrt(4 pi) b^t u^i/u^0 ) / ( alpha sqrt(4 pi) )
  // b^i = ( B^i/u^0 +  alpha sqrt(4 pi) b^t v^i ) / ( alpha sqrt(4 pi) )
  smallb[1] = (prims->Bx*ONE_OVER_U0 + prims->vx*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
  smallb[2] = (prims->By*ONE_OVER_U0 + prims->vy*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
  smallb[3] = (prims->Bz*ONE_OVER_U0 + prims->vz*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
  // Eq. 23 in http://arxiv.org/pdf/astro-ph/0503420.pdf, with alpha sqrt (4 pi) b^2 = u_i B^i already computed above
  smallb[0] = alpha_sqrt_4pi_bt * ONE_OVER_LAPSE_SQRT_4PI;

  // b^2 = g_{\mu \nu} b^{\mu} b^{\nu}
  //     = gtt bt^2 + gxx bx^2 + gyy by^2 + gzz bz^2 + 2 (gtx bt bx + gty bt by + gtz bt bz + gxy bx by + gxz bx bz + gyz by bz)
  //     = (-al^2 + gamma_{ij} betai betaj) bt^2 + b^i b^j gamma_{ij} + 2 g_{t i} b^t b^i
  //     = - (alpha b^t)^2 + (b^t)^2 gamma_{ij} beta^i beta^j + b^i b^j gamma_{ij} + 2 b^t g_{t i} b^i
  //     = - (alpha b^t)^2 + (b^t)^2 gamma_{ij} beta^i beta^j + b^i b^j gamma_{ij} + 2 b^t (gamma_{ij} beta^j) b^i
  //     = - (alpha b^t)^2 + gamma_{ij} ((b^t)^2 beta^i beta^j + b^i b^j + 2 b^t beta^j b^i)
  //     = - (alpha b^t)^2 + gamma_{ij} ((b^t)^2 beta^i beta^j + 2 b^t beta^j b^i + b^i b^j)
  //     = - (alpha b^t)^2 + gamma_{ij} (b^i + b^t beta^i) (b^j + b^t beta^j)
  const double bx_plus_shiftx_bt = smallb[1]+metric->betax*smallb[0];
  const double by_plus_shifty_bt = smallb[2]+metric->betay*smallb[0];
  const double bz_plus_shiftz_bt = smallb[3]+metric->betaz*smallb[0];
  *smallb2 = -SQR(metric->lapse*smallb[0]) +
    metric->adm_gxx*SQR(bx_plus_shiftx_bt) + metric->adm_gyy*SQR(by_plus_shifty_bt) + metric->adm_gzz*SQR(bz_plus_shiftz_bt) +
       2.0*( metric->adm_gxy*(bx_plus_shiftx_bt)*(by_plus_shifty_bt) +
             metric->adm_gxz*(bx_plus_shiftx_bt)*(bz_plus_shiftz_bt) +
             metric->adm_gyz*(by_plus_shifty_bt)*(bz_plus_shiftz_bt) );
}
