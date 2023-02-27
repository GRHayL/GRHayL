#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Primitives to Conservatives Routine
 */
void prims_to_cons(const paramstruct *restrict params, REAL *restrict xx[3], const REAL *restrict auxevol_gfs, REAL *restrict evol_gfs) {
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
           * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
           */
          const double alpha = auxevol_gfs[IDX4S(ALPHAGF, i0,i1,i2)];
          const double betaU0 = auxevol_gfs[IDX4S(BETAU0GF, i0,i1,i2)];
          const double betaU1 = auxevol_gfs[IDX4S(BETAU1GF, i0,i1,i2)];
          const double betaU2 = auxevol_gfs[IDX4S(BETAU2GF, i0,i1,i2)];
          const double gammaDD00 = auxevol_gfs[IDX4S(GAMMADD00GF, i0,i1,i2)];
          const double gammaDD01 = auxevol_gfs[IDX4S(GAMMADD01GF, i0,i1,i2)];
          const double gammaDD02 = auxevol_gfs[IDX4S(GAMMADD02GF, i0,i1,i2)];
          const double gammaDD11 = auxevol_gfs[IDX4S(GAMMADD11GF, i0,i1,i2)];
          const double gammaDD12 = auxevol_gfs[IDX4S(GAMMADD12GF, i0,i1,i2)];
          const double gammaDD22 = auxevol_gfs[IDX4S(GAMMADD22GF, i0,i1,i2)];
          const double u4Ut = auxevol_gfs[IDX4S(U4UTGF, i0,i1,i2)];
          const double BU0 = auxevol_gfs[IDX4S(BU0GF, i0,i1,i2)];
          const double BU1 = auxevol_gfs[IDX4S(BU1GF, i0,i1,i2)];
          const double BU2 = auxevol_gfs[IDX4S(BU2GF, i0,i1,i2)];
          const double rhob = auxevol_gfs[IDX4S(RHOBGF, i0,i1,i2)];
          const double P = auxevol_gfs[IDX4S(PGF, i0,i1,i2)];
          const double h = auxevol_gfs[IDX4S(HGF, i0,i1,i2)];
          /*
           * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
           */
          const double FDPart3_0 = betaU0*gammaDD00 + betaU1*gammaDD01 + betaU2*gammaDD02;
          const double FDPart3_2 = betaU0*gammaDD01 + betaU1*gammaDD11 + betaU2*gammaDD12;
          const double FDPart3_4 = betaU0*gammaDD02 + betaU1*gammaDD12 + betaU2*gammaDD22;
          const double FDPart3_6 = BU0*FDPart3_0 + BU1*FDPart3_2 + BU2*FDPart3_4;
          const double FDPart3_7 = ((alpha)*(alpha));
          const double FDPart3_8 = (1.0/(FDPart3_7));
          const double FDPart3_9 = FDPart3_8/((sqrt4pi)*(sqrt4pi));
          const double FDPart3_10 = FDPart3_6*FDPart3_9;
          const double FDPart3_11 = ((u4Ut)*(u4Ut));
          const double FDPart3_12 = FDPart3_9/FDPart3_11;
          const double FDPart3_14 = BU0*BU1*FDPart3_12;
          const double FDPart3_15 = BU0*BU2*FDPart3_12*gammaDD02;
          const double FDPart3_16 = BU1*BU2*FDPart3_12*gammaDD12;
          const double FDPart3_17 = ((BU0)*(BU0))*FDPart3_12*gammaDD00;
          const double FDPart3_18 = ((BU1)*(BU1))*FDPart3_12*gammaDD11;
          const double FDPart3_19 = ((BU2)*(BU2))*FDPart3_12*gammaDD22;
          const double FDPart3_20 = BU0*FDPart3_0*FDPart3_10;
          const double FDPart3_21 = BU1*FDPart3_10*FDPart3_2;
          const double FDPart3_22 = BU2*FDPart3_10*FDPart3_4;
          const double FDPart3_23 = FDPart3_11*((FDPart3_6)*(FDPart3_6))*FDPart3_9;
          const double FDPart3_24 = FDPart3_23*(FDPart3_0*betaU0 + FDPart3_2*betaU1 + FDPart3_4*betaU2 - FDPart3_7);
          const double FDPart3_25 = FDPart3_14*gammaDD01 + FDPart3_15 + FDPart3_16 + (1.0/2.0)*FDPart3_17 + (1.0/2.0)*FDPart3_18 + (1.0/2.0)*FDPart3_19 + FDPart3_20 + FDPart3_21 + FDPart3_22 + (1.0/2.0)*FDPart3_24 + P;
          const double FDPart3_26 = -BU0*FDPart3_10 + FDPart3_25*FDPart3_8*betaU0;
          const double FDPart3_27 = -BU1*FDPart3_10 + FDPart3_25*FDPart3_8*betaU1;
          const double FDPart3_28 = -BU2*FDPart3_10 + FDPart3_25*FDPart3_8*betaU2;
          const double FDPart3_29 = FDPart3_11*(2*FDPart3_14*gammaDD01 + 2*FDPart3_15 + 2*FDPart3_16 + FDPart3_17 + FDPart3_18 + FDPart3_19 + 2*FDPart3_20 + 2*FDPart3_21 + 2*FDPart3_22 + FDPart3_24 + h*rhob) - FDPart3_23 - FDPart3_25*FDPart3_8;
          const double FDPart3_30 = sqrt(gammaDD00*gammaDD11*gammaDD22 - gammaDD00*((gammaDD12)*(gammaDD12)) - ((gammaDD01)*(gammaDD01))*gammaDD22 + 2*gammaDD01*gammaDD02*gammaDD12 - ((gammaDD02)*(gammaDD02))*gammaDD11);
          const double FDPart3_31 = FDPart3_30*alpha;
          const double FDPart3_32 = FDPart3_31*rhob*u4Ut;
          evol_gfs[IDX4S(STILDED0GF, i0, i1, i2)] = FDPart3_31*(FDPart3_0*FDPart3_29 + FDPart3_26*gammaDD00 + FDPart3_27*gammaDD01 + FDPart3_28*gammaDD02);
          evol_gfs[IDX4S(STILDED1GF, i0, i1, i2)] = FDPart3_31*(FDPart3_2*FDPart3_29 + FDPart3_26*gammaDD01 + FDPart3_27*gammaDD11 + FDPart3_28*gammaDD12);
          evol_gfs[IDX4S(STILDED2GF, i0, i1, i2)] = FDPart3_31*(FDPart3_26*gammaDD02 + FDPart3_27*gammaDD12 + FDPart3_28*gammaDD22 + FDPart3_29*FDPart3_4);
          evol_gfs[IDX4S(RHO_STARGF, i0, i1, i2)] = FDPart3_32;
          evol_gfs[IDX4S(TAU_TILDEGF, i0, i1, i2)] = FDPart3_29*FDPart3_30*FDPart3_7 - FDPart3_32;
        }
      } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
    } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
  } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
}
