#include "./NRPy_basic_defines.h"
/*
 * Calculate psi3
 */
static inline void calculate_psi3(const paramstruct *restrict params, const int i0, const int i1, const int i2, const REAL *restrict auxevol_gfs, REAL *psi3) {
#include "./set_Cparameters.h"

{
  /*
   * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
   */
  const double gammaDD00 = auxevol_gfs[IDX4S(GAMMADD00GF, i0,i1,i2)];
  const double gammaDD01 = auxevol_gfs[IDX4S(GAMMADD01GF, i0,i1,i2)];
  const double gammaDD02 = auxevol_gfs[IDX4S(GAMMADD02GF, i0,i1,i2)];
  const double gammaDD11 = auxevol_gfs[IDX4S(GAMMADD11GF, i0,i1,i2)];
  const double gammaDD12 = auxevol_gfs[IDX4S(GAMMADD12GF, i0,i1,i2)];
  const double gammaDD22 = auxevol_gfs[IDX4S(GAMMADD22GF, i0,i1,i2)];
  /*
   * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
   */
  *psi3 = pow(gammaDD00*gammaDD11*gammaDD22 - gammaDD00*((gammaDD12)*(gammaDD12)) - ((gammaDD01)*(gammaDD01))*gammaDD22 + 2*gammaDD01*gammaDD02*gammaDD12 - ((gammaDD02)*(gammaDD02))*gammaDD11, 1.0/4.0);
}}

/*
 * Calculate magnetic field components from the vector potential
 */
void compute_B_and_Bstagger_from_A(const paramstruct *restrict params, const REAL *evol_gfs, REAL *auxevol_gfs) {
#include "./set_Cparameters.h"

  #pragma omp parallel for
  for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
    for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
      for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
        
        
        const int i = i0;
        const int j = i1;
        const int k = i2;
        
        int shiftedim1 = (i-1)*(i!=0); // This way, i=0 yields shiftedim1=0 and shiftedi=1, used below for our COPY boundary condition.
        int shiftedi   = shiftedim1+1;
        
        int shiftedjm1 = (j-1)*(j!=0);
        int shiftedj   = shiftedjm1+1;
        
        int shiftedkm1 = (k-1)*(k!=0);
        int shiftedk   = shiftedkm1+1;
        
        int index,indexim1,indexjm1,indexkm1;
        
        //int actual_index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        int actual_index = IDX3S(i,j,k);
        
        /**************/
        /* Bx_stagger */
        /**************/
        
        index    = IDX3S(i,shiftedj,shiftedk);
        indexjm1 = IDX3S(i,shiftedjm1,shiftedk);
        indexkm1 = IDX3S(i,shiftedj,shiftedkm1);
        
        // Set Bx_stagger = \partial_y A_z - partial_z A_y
        // "Grid" Ax(i,j,k) is actually Ax(i,j+1/2,k+1/2)
        // "Grid" Ay(i,j,k) is actually Ay(i+1/2,j,k+1/2)
        // "Grid" Az(i,j,k) is actually Az(i+1/2,j+1/2,k)
        // Therefore, the 2nd order derivative \partial_z A_y at (i+1/2,j,k) is:
        //          ["Grid" Ay(i,j,k) - "Grid" Ay(i,j,k-1)]/dZ
        //Bx_stagger[actual_index] = (Az[index]-Az[indexjm1])*dyi - (Ay[index]-Ay[indexkm1])*dzi;
        
        REAL dAz_dy = (evol_gfs[IDX4ptS(AD2GF, index)] - evol_gfs[IDX4ptS(AD2GF, indexjm1)])/dxx1;
        REAL dAy_dz = (evol_gfs[IDX4ptS(AD1GF, index)] - evol_gfs[IDX4ptS(AD1GF, indexkm1)])/dxx2;
        auxevol_gfs[IDX4ptS(BSTAGGERU0GF, actual_index)] = dAz_dy - dAy_dz;
        
        // Now multiply Bx and Bx_stagger by 1/sqrt(gamma(i+1/2,j,k)]) = 1/sqrt(1/2 [gamma + gamma_ip1]) = exp(-6 x 1/2 [phi + phi_ip1] )
        // int imax_minus_i = (cctk_lsh[0]-1)-i;
        // int indexip1jk = CCTK_GFINDEX3D(cctkGH,i + ( (imax_minus_i > 0) - (0 > imax_minus_i) ),j,k);
        // CCTK_REAL Psi_ip1 = psi_bssn[indexip1jk];
        // Bx_stagger[actual_index] *= Psim3/(Psi_ip1*Psi_ip1*Psi_ip1);
        
        int imax_minus_i = (Nxx_plus_2NGHOSTS0-1)-i;
        REAL Psi3, Psi3_ip1;
        calculate_psi3(params, i, j, k, auxevol_gfs, &Psi3);
        calculate_psi3(params, i + ( (imax_minus_i > 0) - (0 > imax_minus_i) ), j, k, auxevol_gfs, &Psi3_ip1);
        auxevol_gfs[IDX4ptS(BSTAGGERU0GF, actual_index)] /= (Psi3*Psi3_ip1);
        
        
        /**************/
        /* By_stagger */
        /**************/
        
        index    = IDX3S(shiftedi,j,shiftedk);
        indexim1 = IDX3S(shiftedim1,j,shiftedk);
        indexkm1 = IDX3S(shiftedi,j,shiftedkm1);
        // Set By_stagger = \partial_z A_x - \partial_x A_z
        //By_stagger[actual_index] = (Ax[index]-Ax[indexkm1])*dzi - (Az[index]-Az[indexim1])*dxi;
        
        REAL dAx_dz = (evol_gfs[IDX4ptS(AD0GF, index)] - evol_gfs[IDX4ptS(AD0GF, indexkm1)])/dxx2;
        REAL dAz_dx = (evol_gfs[IDX4ptS(AD2GF, index)] - evol_gfs[IDX4ptS(AD2GF, indexim1)])/dxx0;
        auxevol_gfs[IDX4ptS(BSTAGGERU1GF, actual_index)] = dAx_dz - dAz_dx;
        
        
        // Now multiply By and By_stagger by 1/sqrt(gamma(i,j+1/2,k)]) = 1/sqrt(1/2 [gamma + gamma_jp1]) = exp(-6 x 1/2 [phi + phi_jp1] )
        //int jmax_minus_j = (cctk_lsh[1]-1)-j;
        //int indexijp1k = CCTK_GFINDEX3D(cctkGH,i,j + ( (jmax_minus_j > 0) - (0 > jmax_minus_j) ),k);
        //CCTK_REAL Psi_jp1 = psi_bssn[indexijp1k];
        //By_stagger[actual_index] *= Psim3/(Psi_jp1*Psi_jp1*Psi_jp1);
        
        int jmax_minus_j = (Nxx_plus_2NGHOSTS1-1)-j;
        REAL Psi3_jp1;
        calculate_psi3(params, i,j + ( (jmax_minus_j > 0) - (0 > jmax_minus_j) ),k , auxevol_gfs, &Psi3_jp1);
        auxevol_gfs[IDX4ptS(BSTAGGERU1GF, actual_index)] /= (Psi3*Psi3_jp1);
        
        
        /**************/
        /* Bz_stagger */
        /**************/
        
        index    = IDX3S(shiftedi,shiftedj,k);
        indexim1 = IDX3S(shiftedim1,shiftedj,k);
        indexjm1 = IDX3S(shiftedi,shiftedjm1,k);
        // Set Bz_stagger = \partial_x A_y - \partial_y A_x
        //Bz_stagger[actual_index] = (Ay[index]-Ay[indexim1])*dxi - (Ax[index]-Ax[indexjm1])*dyi;
        
        REAL dAy_dx = (evol_gfs[IDX4ptS(AD1GF, index)] - evol_gfs[IDX4ptS(AD1GF, indexim1)])/dxx0;
        REAL dAx_dy = (evol_gfs[IDX4ptS(AD0GF, index)] - evol_gfs[IDX4ptS(AD0GF, indexjm1)])/dxx1;
        auxevol_gfs[IDX4ptS(BSTAGGERU2GF, actual_index)] = dAy_dx - dAx_dy;
        
        
        // Now multiply Bz_stagger by 1/sqrt(gamma(i,j,k+1/2)]) = 1/sqrt(1/2 [gamma + gamma_kp1]) = exp(-6 x 1/2 [phi + phi_kp1] )
        // int kmax_minus_k = (cctk_lsh[2]-1)-k;
        // int indexijkp1 = CCTK_GFINDEX3D(cctkGH,i,j,k + ( (kmax_minus_k > 0) - (0 > kmax_minus_k) ));
        // CCTK_REAL Psi_kp1 = psi_bssn[indexijkp1];
        // Bz_stagger[actual_index] *= Psim3/(Psi_kp1*Psi_kp1*Psi_kp1);
        
        int kmax_minus_k = (Nxx_plus_2NGHOSTS2-1)-k;
        REAL Psi3_kp1;
        calculate_psi3(params, i,j,k + ( (kmax_minus_k > 0) - (0 > kmax_minus_k) ), auxevol_gfs, &Psi3_kp1);
        auxevol_gfs[IDX4ptS(BSTAGGERU2GF, actual_index)] /= (Psi3*Psi3_kp1);
        
        
      } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
    } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
  } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)


LOOP_REGION(0, Nxx_plus_2NGHOSTS0,
            0, Nxx_plus_2NGHOSTS1,
            0, Nxx_plus_2NGHOSTS2) {
            
    const int i = i0;
    const int j = i1;
    const int k = i2;

    // Look Mom, no if() statements!
    int shiftedim1 = (i-1)*(i!=0); // This way, i=0 yields shiftedim1=0 and shiftedi=1, used below for our COPY boundary condition.
    int shiftedi   = shiftedim1+1;

    int shiftedjm1 = (j-1)*(j!=0);
    int shiftedj   = shiftedjm1+1;

    int shiftedkm1 = (k-1)*(k!=0);
    int shiftedk   = shiftedkm1+1;

    int index,indexim1,indexjm1,indexkm1;

    int actual_index = IDX3S(i,j,k);

    // For the lower boundaries, the following applies a "copy"
    //    boundary condition on Bi and Bi_stagger where needed.
    //    E.g., Bx(imin,j,k) = Bx(imin+1,j,k)
    //    We find the copy BC works better than extrapolation.
    /******/
    /* Bx */
    /******/
    index = IDX3S(shiftedi,j,k);
    indexim1 = IDX3S(shiftedim1,j,k);
    // Set Bx = 0.5 ( Bx_stagger + Bx_stagger_im1 )
    // "Grid" Bx_stagger(i,j,k) is actually Bx_stagger(i+1/2,j,k)
    //Bx[actual_index] = 0.5 * ( Bx_stagger[index] + Bx_stagger[indexim1] );
    auxevol_gfs[IDX4ptS(BU0GF, actual_index)] = 0.5*(auxevol_gfs[IDX4ptS(BSTAGGERU0GF, index)] + 
                                                     auxevol_gfs[IDX4ptS(BSTAGGERU0GF, indexim1)]);

    /******/
    /* By */
    /******/
    index = IDX3S(i,shiftedj,k);
    indexjm1 = IDX3S(i,shiftedjm1,k);
    // Set By = 0.5 ( By_stagger + By_stagger_im1 )
    // "Grid" By_stagger(i,j,k) is actually By_stagger(i,j+1/2,k)
    //By[actual_index] = 0.5 * ( By_stagger[index] + By_stagger[indexjm1] );
    auxevol_gfs[IDX4ptS(BU1GF, actual_index)] = 0.5*(auxevol_gfs[IDX4ptS(BSTAGGERU1GF, index)] + 
                                                     auxevol_gfs[IDX4ptS(BSTAGGERU1GF, indexjm1)]);

    /******/
    /* Bz */
    /******/
    index = IDX3S(i,j,shiftedk);
    indexkm1 = IDX3S(i,j,shiftedkm1);
    // Set Bz = 0.5 ( Bz_stagger + Bz_stagger_im1 )
    // "Grid" Bz_stagger(i,j,k) is actually Bz_stagger(i,j+1/2,k)
    //Bz[actual_index] = 0.5 * ( Bz_stagger[index] + Bz_stagger[indexkm1] );
    auxevol_gfs[IDX4ptS(BU2GF, actual_index)] = 0.5*(auxevol_gfs[IDX4ptS(BSTAGGERU2GF, index)] + 
                                                     auxevol_gfs[IDX4ptS(BSTAGGERU2GF, indexkm1)]);       
                                                     
  }

}
