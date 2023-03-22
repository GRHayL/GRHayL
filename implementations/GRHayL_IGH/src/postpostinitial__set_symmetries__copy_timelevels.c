//-------------------------------------------------
// Stuff to run right after initial data is set up
//-------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"
#include "IGH.h"

void GRHayL_IGH_PostPostInitial_Set_Symmetries__Copy_Timelevels(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayL_IGH_PostPostInitial_Set_Symmetries__Copy_Timelevels;
  DECLARE_CCTK_PARAMETERS;

  //For emfields, we assume that you've set Bx, By, Bz (the UN-tilded B^i's)
  // or Ax, Ay, Az (if using constrained transport scheme of Del Zanna)

  if(CCTK_EQUALS(Symmetry,"equatorial")) {
    // SET SYMMETRY GHOSTZONES ON ALL CONSERVATIVE AND PRIMIIVE VARIABLES!
    int ierr;
    ierr=CartSymGN(cctkGH,"GRHayL_IGH::grmhd_conservatives"); if(ierr!=0) CCTK_VERROR("GRHayL_IGH error code #1874109358120048. Grep it in the source code");
    ierr=CartSymGN(cctkGH,"GRHayL_IGH::grmhd_primitives"); if(ierr!=0) CCTK_VERROR("GRHayL_IGH error code #1874109358120049. Grep it in the source code");
  }


  //------------------------------------------------------------------
  // FILL _p AND _p_p TIMELEVELS. Probably don't need to do this if
  // Carpet::init_fill_timelevels=yes  and
  // MoL::initial_data_is_crap = yes
  // NOTE: We don't fill metric data here.
  // FIXME: Do we really need this?

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++)
    for(int j=0; j<cctk_lsh[1]; j++)
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        rho_star_p[index] = rho_star[index];
        tau_p[index] = tau[index];
        Stildex_p[index] = Stildex[index];
        Stildey_p[index] = Stildey[index];
        Stildez_p[index] = Stildez[index];

        rho_star_p_p[index] = rho_star[index];
        tau_p_p[index] = tau[index];
        Stildex_p_p[index] = Stildex[index];
        Stildey_p_p[index] = Stildey[index];
        Stildez_p_p[index] = Stildez[index];
  }
}
