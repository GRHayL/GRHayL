#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "GRHayLib.h"


#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

#define CHECK_PARAMETER(par) if(par==-1) CCTK_VERROR("Please set %s::%s in your parfile",CCTK_THORNSTRING,#par);

/*
 *
 * (c) 2021 Leo Werneck
 *
 * This is the thorn's driver function, responsible
 * for setting the initial data to that of an isotropic
 * gas of constant density, temperature, and electron
 * fraction in Minkowski space.
 */
void IsotropicGasID(CCTK_ARGUMENTS) {


  // Step 1: Get access to gridfunctions and parameters
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Step 2: Check correct usage
  CHECK_PARAMETER(IsotropicGasID_rho);
  CHECK_PARAMETER(IsotropicGasID_Y_e);
  CHECK_PARAMETER(IsotropicGasID_temperature);

  CCTK_INFO("Beginning initial data");

  // Step 3: Compute local pressure and epsilon
  CCTK_REAL IsotropicGasID_press, IsotropicGasID_eps, IsotropicGasID_entropy;
  ghl_tabulated_compute_P_eps_S_from_T(ghl_eos,
                                               IsotropicGasID_rho,
                                               IsotropicGasID_Y_e,
                                               IsotropicGasID_temperature,
                                               &IsotropicGasID_press,
                                               &IsotropicGasID_eps,
                                               &IsotropicGasID_entropy);

  // Step 4: Loop over the grid and set the ID
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) {
    for(int j=0;j<cctk_lsh[1];j++) {
      for(int i=0;i<cctk_lsh[0];i++) {

        const int idx = CCTK_GFINDEX3D(cctkGH,i,j,k);

        Y_e[idx] = IsotropicGasID_Y_e;
        alp[idx] = 1;
        betax[idx] = 0;
        betay[idx] = 0;
        betaz[idx] = 0;
        eps[idx] = IsotropicGasID_eps;
        gxx[idx] = 1;
        gxy[idx] = 0;
        gxz[idx] = 0;
        gyy[idx] = 1;
        gyz[idx] = 0;
        gzz[idx] = 1;
        kxx[idx] = 0;
        kxy[idx] = 0;
        kxz[idx] = 0;
        kyy[idx] = 0;
        kyz[idx] = 0;
        kzz[idx] = 0;
        press[idx] = IsotropicGasID_press;
        rho[idx] = IsotropicGasID_rho;
        temperature[idx] = IsotropicGasID_temperature;
        velx[idx] = 0;
        vely[idx] = 0;
        velz[idx] = 0;
        entropy[idx] = IsotropicGasID_entropy;

      }
    }
  }

  CCTK_INFO("All done!");
}
