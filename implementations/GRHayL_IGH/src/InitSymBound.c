/*
  Set the symmetries for the GRHayL_IGH variables
*/

#include "IGH.h"
#include "Symmetry.h"

void GRHayL_IGH_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_GRHayL_IGH_InitSymBound;
  DECLARE_CCTK_PARAMETERS;


  if( sizeof(CCTK_REAL) < 8 ) CCTK_VERROR("Error: GRHayL_IGH assumes that CCTK_REAL is a double precision number. Setting otherwise will likely cause havoc with the conserv_to_prims solver.");

  if ((cctk_nghostzones[0]<3 || cctk_nghostzones[1]<3 || cctk_nghostzones[2]<3))
    CCTK_VERROR("ERROR: The version of PPM in this thorn requires 3 ghostzones. You only have (%d,%d,%d) ghostzones!",cctk_nghostzones[0],cctk_nghostzones[1],cctk_nghostzones[2]);

  if(cctk_iteration==0) {
    CCTK_VINFO("Setting Symmetry = %s... at iteration = %d",Symmetry,cctk_iteration);

    int sym[3];

    if(CCTK_EQUALS(Symmetry,"none")) {
      /* FIRST SET NO SYMMETRY OPTION */
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      SetCartSymGN(cctkGH,sym,"GRHayL_IGH::grmhd_conservatives");
      SetCartSymGN(cctkGH,sym,"GRHayL_IGH::grmhd_primitives");
    } else if(CCTK_EQUALS(Symmetry,"equatorial")) {
      /* THEN SET EQUATORIAL SYMMETRY OPTION */
      // Set default to no symmetry, which is correct for scalars and most vectors:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      SetCartSymGN(cctkGH,sym,"GRHayL_IGH::grmhd_conservatives");
      SetCartSymGN(cctkGH,sym,"GRHayL_IGH::grmhd_primitives");

      sym[2] = -1;
      SetCartSymVN(cctkGH, sym,"GRHayL_IGH::Stilde_z");
      SetCartSymVN(cctkGH, sym,"GRHayL_IGH::vz");
    } else {
      CCTK_VERROR("GRHayL_IGH_initsymbound: Should not be here; picked an impossible symmetry.");
    }
  }
}
