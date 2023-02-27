#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Generate 2D Loop Advection initial data
 */
void Loop_Advection(const paramstruct *restrict params,REAL *restrict xx[3], REAL *restrict evol_gfs, REAL *restrict auxevol_gfs) {
#include "./set_Cparameters.h"

  #pragma omp parallel for
  for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
    const REAL xx2 = xx[2][i2];
    for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
      const REAL xx1 = xx[1][i1];
      for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
        const REAL xx0 = xx[0][i0];
        {
          /*
           * NRPy+ Finite Difference Code Generation, Step 1 of 1: Evaluate SymPy expressions and write to main memory:
           */
          const double FDPart3_1 = sqrt(((xx2)*(xx2)) + ((0.5*dxx0 + xx0)*(0.5*dxx0 + xx0)) + ((0.5*dxx1 + xx1)*(0.5*dxx1 + xx1)));
          const double FDPart3_2 = sqrt(((xx0)*(xx0)) + ((xx1)*(xx1)) + ((xx2)*(xx2)));
          const double FDPart3_3 = ((1.0/2.0)*FDPart3_2 - 1.0/2.0*fabs(FDPart3_2 - 3.0/10.0) - 3.0/20.0)/(FDPart3_2 - TINYDOUBLE - 3.0/10.0);
          evol_gfs[IDX4S(AD0GF, i0, i1, i2)] = 0;
          evol_gfs[IDX4S(AD1GF, i0, i1, i2)] = 0;
          evol_gfs[IDX4S(AD2GF, i0, i1, i2)] = -1152921504606847.0/2305843009213693952.0*FDPart3_1 + (1.0/2.0)*fabs(3458764513820541.0/11529215046068469760.0 - 1152921504606847.0/1152921504606846976.0*FDPart3_1) + 3458764513820541.0/23058430092136939520.0;
          auxevol_gfs[IDX4S(VU0GF, i0, i1, i2)] = (1.0/12.0)*FDPart3_3;
          auxevol_gfs[IDX4S(VU1GF, i0, i1, i2)] = (1.0/24.0)*FDPart3_3;
          auxevol_gfs[IDX4S(VU2GF, i0, i1, i2)] = (1.0/24.0)*FDPart3_3;
          auxevol_gfs[IDX4S(RHOBGF, i0, i1, i2)] = 1;
          auxevol_gfs[IDX4S(PGF, i0, i1, i2)] = 3;
          evol_gfs[IDX4S(PSI6PHIGF, i0, i1, i2)] = 0;
        }
      } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
    } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
  } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
}
