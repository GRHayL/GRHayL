/*
  Set the symmetries for the GRHayLMHD variables
*/

#include "GRHayLMHD.h"
#include "Symmetry.h"

void GRHayLMHD_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_GRHayLMHD_InitSymBound;
  DECLARE_CCTK_PARAMETERS;

  if( ( CCTK_EQUALS(Matter_BC,"frozen") && !CCTK_EQUALS(EM_BC,"frozen") ) ||
      ( !CCTK_EQUALS(Matter_BC,"frozen") && CCTK_EQUALS(EM_BC,"frozen") ) )
    CCTK_VERROR("If Matter_BC or EM_BC is set to FROZEN, BOTH must be set to frozen!");

  if( sizeof(CCTK_REAL) < 8 ) CCTK_VERROR("Error: GRHayLMHD assumes that CCTK_REAL is a double precision number. Setting otherwise will likely cause havoc with the conserv_to_prims solver.");

  if ((cctk_nghostzones[0]<3 || cctk_nghostzones[1]<3 || cctk_nghostzones[2]<3))
    CCTK_VERROR("ERROR: The version of PPM in this thorn requires 3 ghostzones. You only have (%d,%d,%d) ghostzones!",cctk_nghostzones[0],cctk_nghostzones[1],cctk_nghostzones[2]);

  if(cctk_iteration==0) {
    CCTK_VINFO("Setting Symmetry = %s... at iteration = %d",Symmetry,cctk_iteration);

    int sym[3];

    if(CCTK_EQUALS(Symmetry,"none")) {
      /* FIRST SET NO SYMMETRY OPTION */
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      SetCartSymGN(cctkGH,sym,"GRHayLMHD::grmhd_conservatives");
      SetCartSymGN(cctkGH,sym,"GRHayLMHD::em_Ax");
      SetCartSymGN(cctkGH,sym,"GRHayLMHD::em_Ay");
      SetCartSymGN(cctkGH,sym,"GRHayLMHD::em_Az");
      SetCartSymGN(cctkGH,sym,"GRHayLMHD::em_psi6phi");
      SetCartSymGN(cctkGH,sym,"GRHayLMHD::grmhd_primitives_allbutBi");
    } else if(CCTK_EQUALS(Symmetry,"equatorial")) {
      /* THEN SET EQUATORIAL SYMMETRY OPTION */
      // Set default to no symmetry, which is correct for scalars and most vectors:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      SetCartSymGN(cctkGH,sym,"GRHayLMHD::grmhd_conservatives");
      // Don't worry about the wrong sym values since A_{\mu} is staggered
      // and we're going to impose the symmetry separately
      SetCartSymGN(cctkGH,sym,"GRHayLMHD::em_Ax");
      SetCartSymGN(cctkGH,sym,"GRHayLMHD::em_Ay");
      SetCartSymGN(cctkGH,sym,"GRHayLMHD::em_Az");
      SetCartSymGN(cctkGH,sym,"GRHayLMHD::em_psi6phi");

      SetCartSymGN(cctkGH,sym,"GRHayLMHD::grmhd_primitives_allbutBi");

      // Then set unstaggered B field variables
      sym[2] = -Sym_Bz;
      SetCartSymVN(cctkGH, sym,"GRHayLMHD::Bx_center");
      SetCartSymVN(cctkGH, sym,"GRHayLMHD::By_center");

      sym[2] = Sym_Bz;
      SetCartSymVN(cctkGH, sym,"GRHayLMHD::Bz_center");

      sym[2] = -1;
      SetCartSymVN(cctkGH, sym,"GRHayLMHD::Stilde_z");
      SetCartSymVN(cctkGH, sym,"GRHayLMHD::vz");
    } else {
      CCTK_VERROR("GRHayLMHD_initsymbound: Should not be here; picked an impossible symmetry.");
    }
  }
}
