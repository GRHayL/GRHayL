#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Generate Balsara4 1D initial data
 */
void Balsara4(const paramstruct *restrict params,REAL *restrict xx[3], REAL *restrict evol_gfs, REAL *restrict auxevol_gfs) {
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
          const double FDPart3_0 = (1.0/(-TINYDOUBLE + xx0));
          const double FDPart3_1 = (1.0/2.0)*xx0;
          const double FDPart3_2 = FDPart3_1 - 1.0/2.0*fabs(xx0);
          const double FDPart3_3 = TINYDOUBLE + xx0;
          const double FDPart3_4 = FDPart3_1 + (1.0/2.0)*TINYDOUBLE;
          const double FDPart3_5 = (FDPart3_4 + (1.0/2.0)*fabs(FDPart3_3))/FDPart3_3;
          const double FDPart3_6 = 0.5*dxx0;
          const double FDPart3_7 = FDPart3_6 + xx0;
          const double FDPart3_8 = 0.25*dxx0;
          evol_gfs[IDX4S(AD0GF, i0, i1, i2)] = (0.5*dxx2 + xx2)*(7*FDPart3_0*FDPart3_2 - 7*FDPart3_5);
          evol_gfs[IDX4S(AD1GF, i0, i1, i2)] = FDPart3_7*(7*(FDPart3_1 + FDPart3_8 - 1.0/2.0*fabs(FDPart3_7))/(FDPart3_7 - TINYDOUBLE) - 7*(FDPart3_4 + FDPart3_8 + (1.0/2.0)*fabs(FDPart3_3 + FDPart3_6))/(FDPart3_3 + FDPart3_6));
          evol_gfs[IDX4S(AD2GF, i0, i1, i2)] = 5.0*dxx1 + 10*xx1;
          auxevol_gfs[IDX4S(VU0GF, i0, i1, i2)] = (999.0/1000.0)*FDPart3_0*FDPart3_2 - 999.0/1000.0*FDPart3_5;
          auxevol_gfs[IDX4S(VU1GF, i0, i1, i2)] = 0;
          auxevol_gfs[IDX4S(VU2GF, i0, i1, i2)] = 0;
          auxevol_gfs[IDX4S(RHOBGF, i0, i1, i2)] = 1;
          auxevol_gfs[IDX4S(PGF, i0, i1, i2)] = 1.0/10.0;
          evol_gfs[IDX4S(PSI6PHIGF, i0, i1, i2)] = 0;
        }
      } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
    } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
  } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
}
