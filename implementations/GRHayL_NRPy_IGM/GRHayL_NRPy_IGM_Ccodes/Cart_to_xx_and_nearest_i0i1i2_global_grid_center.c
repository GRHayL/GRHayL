#include "././NRPy_basic_defines.h"
/*
 * Given Cartesian point (x,y,z), this function outputs the corresponding
 *   (xx0,xx1,xx2) and the "closest" (i0,i1,i2) for the given grid
 */
void Cart_to_xx_and_nearest_i0i1i2_global_grid_center(const paramstruct *restrict params, const REAL xCart[3], REAL xx[3], int Cart_to_i0i1i2[3]) {
#include "./set_Cparameters.h"

  const REAL Cartx = xCart[0];
  const REAL Carty = xCart[1];
  const REAL Cartz = xCart[2];

  /*
   *  Original SymPy expressions:
   *  "[xx[0] = Cartx,
   *    xx[1] = Carty,
   *    xx[2] = Cartz]"
   */
  xx[0] = Cartx;
  xx[1] = Carty;
  xx[2] = Cartz;

  // Then find the nearest index (i0,i1,i2) on underlying grid to (x,y,z)
  // Recall that:
  // xx[0][j] = xxmin[0] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*params->dxx0; // Cell-centered grid.
  //   --> j = (int) ( (xx[0][j] - xxmin[0]) / params->dxx0 + (1.0/2.0) + NGHOSTS )
  Cart_to_i0i1i2[0] = (int)( ( xx[0] - (xmin) ) / params->dxx0 + (1.0/2.0) + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
  Cart_to_i0i1i2[1] = (int)( ( xx[1] - (ymin) ) / params->dxx1 + (1.0/2.0) + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
  Cart_to_i0i1i2[2] = (int)( ( xx[2] - (zmin) ) / params->dxx2 + (1.0/2.0) + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
}
