#include "././NRPy_basic_defines.h"
/*
 * Override default values for Nxx{0,1,2}, Nxx_plus_2NGHOSTS{0,1,2}, dxx{0,1,2}, and invdx{0,1,2}; and set xx[3][]
 */
void set_Nxx_dxx_invdx_params__and__xx(const int EigenCoord, const int Nxx[3],paramstruct *restrict params, REAL *restrict xx[3]) {


  // Set CoordSystemName
  snprintf(params->CoordSystemName, 100, "Cartesian");

  // Override parameter defaults with values based on command line arguments and NGHOSTS.
  params->Nxx0 = Nxx[0];
  params->Nxx1 = Nxx[1];
  params->Nxx2 = Nxx[2];

  params->Nxx_plus_2NGHOSTS0 = Nxx[0] + 2*NGHOSTS;
  params->Nxx_plus_2NGHOSTS1 = Nxx[1] + 2*NGHOSTS;
  params->Nxx_plus_2NGHOSTS2 = Nxx[2] + 2*NGHOSTS;
  // Now that params->Nxx_plus_2NGHOSTS* has been set, and below we need e.g., Nxx_plus_2NGHOSTS*, we include set_Cparameters.h here:
#include "./set_Cparameters.h"
  // Step 0d: Set up space and time coordinates
  // Step 0d.i: Declare Delta x^i=dxx{0,1,2} and invdxx{0,1,2}, as well as xxmin[3] and xxmax[3]:
  REAL xxmin[3],xxmax[3];
  if(EigenCoord == 0) {
    xxmin[0] = xmin;
    xxmax[0] = xmax;
    xxmin[1] = ymin;
    xxmax[1] = ymax;
    xxmin[2] = zmin;
    xxmax[2] = zmax;
  } else { // if (EigenCoord == 1)
    xxmin[0] = xmin;
    xxmax[0] = xmax;
    xxmin[1] = ymin;
    xxmax[1] = ymax;
    xxmin[2] = zmin;
    xxmax[2] = zmax;
  }
  // Step 0d.iii: Set params.dxx{0,1,2}, params.invdx{0,1,2}, and uniform coordinate grids xx[3][]
  params->dxx0 = (xxmax[0] - xxmin[0]) / ((REAL)Nxx[0]);
  params->invdx0 = 1.0/params->dxx0;
  xx[0] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0);
#pragma omp parallel for
  for(int j=0;j<Nxx_plus_2NGHOSTS0;j++)
    xx[0][j] = xxmin[0] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*params->dxx0; // Cell-centered grid.

  params->dxx1 = (xxmax[1] - xxmin[1]) / ((REAL)Nxx[1]);
  params->invdx1 = 1.0/params->dxx1;
  xx[1] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS1);
#pragma omp parallel for
  for(int j=0;j<Nxx_plus_2NGHOSTS1;j++)
    xx[1][j] = xxmin[1] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*params->dxx1; // Cell-centered grid.

  params->dxx2 = (xxmax[2] - xxmin[2]) / ((REAL)Nxx[2]);
  params->invdx2 = 1.0/params->dxx2;
  xx[2] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS2);
#pragma omp parallel for
  for(int j=0;j<Nxx_plus_2NGHOSTS2;j++)
    xx[2][j] = xxmin[2] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*params->dxx2; // Cell-centered grid.
}
