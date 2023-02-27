#include "././NRPy_basic_defines.h"
/*
 * Compute Cartesian coordinates given local grid coordinate (xx0,xx1,xx2),   accounting for the origin of this grid being possibly offcenter.
 */
void xx_to_Cart(const paramstruct *restrict params, REAL *restrict xx[3],const int i0,const int i1,const int i2, REAL xCart[3]) {
#include "./set_Cparameters.h"


    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];
      /*
   *  Original SymPy expressions:
   *  "[xCart[0] = Cart_originx + xx0,
   *    xCart[1] = Cart_originy + xx1,
   *    xCart[2] = Cart_originz + xx2]"
   */
  {
    xCart[0] = Cart_originx + xx0;
    xCart[1] = Cart_originy + xx1;
    xCart[2] = Cart_originz + xx2;
  }
}
