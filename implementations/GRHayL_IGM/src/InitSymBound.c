/*
  Set the symmetries for the GRHayL_IGM variables
*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

void GRHayL_IGM_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_GRHayL_IGM_InitSymBound;
  DECLARE_CCTK_PARAMETERS;

  if( ( CCTK_EQUALS(Matter_BC,"frozen") && !CCTK_EQUALS(EM_BC,"frozen") ) ||
      ( !CCTK_EQUALS(Matter_BC,"frozen") && CCTK_EQUALS(EM_BC,"frozen") ) )
    CCTK_VERROR("If Matter_BC or EM_BC is set to FROZEN, BOTH must be set to frozen!");

  if ((cctk_nghostzones[0]<3 || cctk_nghostzones[1]<3 || cctk_nghostzones[2]<3))
    CCTK_VERROR("ERROR: The version of PPM in this thorn requires 3 ghostzones. You only have (%d,%d,%d) ghostzones!",cctk_nghostzones[0],cctk_nghostzones[1],cctk_nghostzones[2]);

  if(cctk_iteration==0) {
    CCTK_VINFO("Setting Symmetry = %s... at iteration = %d",Symmetry,cctk_iteration);

    int sym[3];

    if(CCTK_EQUALS(Symmetry,"none")) {
      /* FIRST SET NO SYMMETRY OPTION */
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      SetCartSymGN(cctkGH,sym,"GRHayL_IGM::grmhd_conservatives");
      SetCartSymGN(cctkGH,sym,"GRHayL_IGM::em_Ax");
      SetCartSymGN(cctkGH,sym,"GRHayL_IGM::em_Ay");
      SetCartSymGN(cctkGH,sym,"GRHayL_IGM::em_Az");
      SetCartSymGN(cctkGH,sym,"GRHayL_IGM::em_psi6phi");
      SetCartSymGN(cctkGH,sym,"GRHayL_IGM::grmhd_primitives_allbutBi");
    } else if(CCTK_EQUALS(Symmetry,"equatorial")) {
      /* THEN SET EQUATORIAL SYMMETRY OPTION */
      // Set default to no symmetry, which is correct for scalars and most vectors:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      SetCartSymGN(cctkGH,sym,"GRHayL_IGM::grmhd_conservatives");
      // Don't worry about the wrong sym values since A_{\mu} is staggered
      // and we're going to impose the symmetry separately
      SetCartSymGN(cctkGH,sym,"GRHayL_IGM::em_Ax");
      SetCartSymGN(cctkGH,sym,"GRHayL_IGM::em_Ay");
      SetCartSymGN(cctkGH,sym,"GRHayL_IGM::em_Az");
      SetCartSymGN(cctkGH,sym,"GRHayL_IGM::em_psi6phi");

      SetCartSymGN(cctkGH,sym,"GRHayL_IGM::grmhd_primitives_allbutBi");

      // Then set unstaggered B field variables
      sym[2] = -Sym_Bz;
      SetCartSymVN(cctkGH, sym,"GRHayL_IGM::Bx_center");
      SetCartSymVN(cctkGH, sym,"GRHayL_IGM::By_center");
      sym[2] = Sym_Bz;
      SetCartSymVN(cctkGH, sym,"GRHayL_IGM::Bz_center");

      sym[2] = -1;
      SetCartSymVN(cctkGH, sym,"GRHayL_IGM::Stilde_z");
      SetCartSymVN(cctkGH, sym,"GRHayL_IGM::vz");
    } else {
      CCTK_VERROR("GRHayL_IGM_initsymbound: Should not be here; picked an impossible symmetry.");
    }
  }
}
