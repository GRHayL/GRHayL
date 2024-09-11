#include "ghl_con2prim.h"

/*
 * Function     : ghl_compute_smallb_and_b2()
 * Description  : Computes magnetic quantities b^0, b^i, and b^2 (see Eqs. 23
 *                and 24 in https://arxiv.org/abs/astro-ph/0503420).
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_compute_smallb_and_b2
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
