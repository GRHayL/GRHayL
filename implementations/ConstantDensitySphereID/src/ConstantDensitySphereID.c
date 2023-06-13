#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "GRHayLib.h"


#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

#define CHECK_PARAMETER(par) if(par==-1) CCTK_VError(__LINE__,__FILE__,CCTK_THORNSTRING,"Please set %s::%s in your parfile",CCTK_THORNSTRING,#par);

/*
 *
 * (c) 2021 Leo Werneck
 *
 * This is the thorn's driver function, responsible
 * for setting the initial data to that of an constant
 * density sphere in Minkowski space.
 */
void ConstantDensitySphereID(CCTK_ARGUMENTS) {


  // Step 1: Get access to gridfunctions and parameters
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Step 2: Check correct usage
  CHECK_PARAMETER(ConstantDensitySphereID_sphere_radius);
  CHECK_PARAMETER(ConstantDensitySphereID_rho_interior);
  CHECK_PARAMETER(ConstantDensitySphereID_Y_e_interior);
  CHECK_PARAMETER(ConstantDensitySphereID_T_interior);
  CHECK_PARAMETER(ConstantDensitySphereID_rho_exterior);
  CHECK_PARAMETER(ConstantDensitySphereID_Y_e_exterior);
  CHECK_PARAMETER(ConstantDensitySphereID_T_exterior);

  CCTK_INFO("Beginning initial data");

  // Step 3: Compute hydro quantities inside and outside the sphere
  // Step 3.a: Sphere interior
  const CCTK_REAL rho_interior = ConstantDensitySphereID_rho_interior;
  const CCTK_REAL Y_e_interior = ConstantDensitySphereID_Y_e_interior;
  const CCTK_REAL T_interior   = ConstantDensitySphereID_T_interior;
  CCTK_REAL P_interior, eps_interior, S_interior;
  ghl_tabulated_compute_P_eps_S_from_T( ghl_eos,
                                                rho_interior,
                                                Y_e_interior,
                                                T_interior,
                                                &P_interior,
                                                &eps_interior,
                                                &S_interior );

  // Step 3.b: Sphere exterior
  const CCTK_REAL rho_exterior = ConstantDensitySphereID_rho_exterior;
  const CCTK_REAL Y_e_exterior = ConstantDensitySphereID_Y_e_exterior;
  const CCTK_REAL T_exterior   = ConstantDensitySphereID_T_exterior;
  CCTK_REAL P_exterior, eps_exterior, S_exterior;
  ghl_tabulated_compute_P_eps_S_from_T( ghl_eos,
                                                rho_exterior,
                                                Y_e_exterior,
                                                T_exterior,
                                                &P_exterior,
                                                &eps_exterior,
                                                &S_exterior );

  // Step 4: Loop over the grid and set the ID
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) {
    for(int j=0;j<cctk_lsh[1];j++) {
      for(int i=0;i<cctk_lsh[0];i++) {

        // Step 4.a: Get gridpoint index
        const int idx = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // Step 4.b: Compute current radius
        const CCTK_REAL xL = x[idx];
        const CCTK_REAL yL = y[idx];
        const CCTK_REAL zL = z[idx];
        const CCTK_REAL rL = sqrt( xL*xL + yL*yL + zL*zL );

        // Step 4.c: Initialize the ADMBase gridfunctions
        alp[idx] = 1;
        betax[idx] = 0;
        betay[idx] = 0;
        betaz[idx] = 0;
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

        // Step 4.d: Initialize the HydroBase gridfunctions
        velx[idx] = 0;
        vely[idx] = 0;
        velz[idx] = 0;
        if( rL > ConstantDensitySphereID_sphere_radius ) {
          // Step 4.d.i: Outside the sphere
          rho        [idx] = rho_exterior;
          Y_e        [idx] = Y_e_exterior;
          temperature[idx] = T_exterior;
          press      [idx] = P_exterior;
          eps        [idx] = eps_exterior;
          entropy    [idx] = S_exterior;
        }
        else {
          // Step 4.d.ii: Inside the sphere
          rho        [idx] = rho_interior;
          Y_e        [idx] = Y_e_interior;
          temperature[idx] = T_interior;
          press      [idx] = P_interior;
          eps        [idx] = eps_interior;
          entropy    [idx] = S_interior;
        }
      }
    }
  }

  CCTK_INFO("All done!");
}
