#include "././NRPy_basic_defines.h"
/*
 * Find the CFL-constrained timestep
 */
REAL find_timestep(const paramstruct *restrict params, REAL *restrict xx[3], const REAL CFL_FACTOR) {
#include "./set_Cparameters.h"
  REAL dsmin = 1e38; // Start with a crazy high value... close to the largest number in single precision.
  for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++) {
    const REAL xx2 = xx[2][i2];
    for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++) {
      const REAL xx1 = xx[1][i1];
      for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++) {
        const REAL xx0 = xx[0][i0];
          REAL ds_dirn0, ds_dirn1, ds_dirn2;
            /*
             *  Original SymPy expressions:
             *  "[ds_dirn0 = dxx0,
             *    ds_dirn1 = dxx1,
             *    ds_dirn2 = dxx2]"
             */
            {
              ds_dirn0 = dxx0;
              ds_dirn1 = dxx1;
              ds_dirn2 = dxx2;
            }
        
        #ifndef MIN
        #define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
        #endif
          // Set dsmin = MIN(dsmin, ds_dirn0, ds_dirn1, ds_dirn2):
          dsmin = MIN(dsmin, MIN(fabs(ds_dirn0), MIN(fabs(ds_dirn1), fabs(ds_dirn2))));
        
      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
return dsmin*CFL_FACTOR/1.0;
}
