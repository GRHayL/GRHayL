//-------------------------------------------------
// Stuff to run right after initial data is set up
//-------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"
#include "IGM.h"

void GRHayL_IGM_PostPostInitial_Set_Symmetries__Copy_Timelevels(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayL_IGM_PostPostInitial_Set_Symmetries__Copy_Timelevels;
  DECLARE_CCTK_PARAMETERS;

  //For emfields, we assume that you've set Bx, By, Bz (the UN-tilded B^i's)
  // or Ax, Ay, Az (if using constrained transport scheme of Del Zanna)

  if(CCTK_EQUALS(Symmetry,"equatorial")) {
    // SET SYMMETRY GHOSTZONES ON ALL CONSERVATIVE AND PRIMIIVE VARIABLES!
    int ierr;
    ierr=CartSymGN(cctkGH,"GRHayL_IGM::grmhd_conservatives"); if(ierr!=0) CCTK_VERROR("GRHayL_IGM error code #1874109358120048. Grep it in the source code");
    ierr=CartSymGN(cctkGH,"GRHayL_IGM::grmhd_primitives_allbutBi"); if(ierr!=0) CCTK_VERROR("GRHayL_IGM error code #1874109358120049. Grep it in the source code");

    // Finish up by setting symmetry ghostzones on Bx, By, Bz, and their staggered variants.
    const double gridfunc_syms_Bx[3] = { -1,  1, -Sym_Bz};
    GRHayL_IGM_set_symmetry_gzs_staggered(cctkGH, x, y, z, Bx_center,  gridfunc_syms_Bx, 0, 0, 0);
    GRHayL_IGM_set_symmetry_gzs_staggered(cctkGH, x, y, z, Bx_stagger, gridfunc_syms_Bx, 1, 0, 0);
    const double gridfunc_syms_By[3] = {  1, -1, -Sym_Bz};
    GRHayL_IGM_set_symmetry_gzs_staggered(cctkGH, x, y, z, By_center,  gridfunc_syms_Bx, 0, 0, 0);
    GRHayL_IGM_set_symmetry_gzs_staggered(cctkGH, x, y, z, By_stagger, gridfunc_syms_By, 0, 1, 0);
    const double gridfunc_syms_Bz[3] = {  1,  1,  Sym_Bz};
    GRHayL_IGM_set_symmetry_gzs_staggered(cctkGH, x, y, z, Bz_center,  gridfunc_syms_Bz, 0, 0, 0);
    GRHayL_IGM_set_symmetry_gzs_staggered(cctkGH, x, y, z, Bz_stagger, gridfunc_syms_Bz, 0, 0, 1);

    const double gridfunc_syms_phitilde[3] = { 1,  1,  1};
    GRHayL_IGM_set_symmetry_gzs_staggered(cctkGH, x, y, z, phitilde,    gridfunc_syms_phitilde,1,1,1);
    const double gridfunc_syms_Ax[3]      = {-1,  1,  Sym_Bz};
    GRHayL_IGM_set_symmetry_gzs_staggered(cctkGH, x, y, z, Ax,         gridfunc_syms_Ax, 0, 1, 1);
    const double gridfunc_syms_Ay[3]      = { 1, -1,  Sym_Bz};
    GRHayL_IGM_set_symmetry_gzs_staggered(cctkGH, x, y, z, Ay,         gridfunc_syms_Ay, 1, 0, 1);
    const double gridfunc_syms_Az[3]      = { 1,  1, -Sym_Bz};
    GRHayL_IGM_set_symmetry_gzs_staggered(cctkGH, x, y, z, Az,         gridfunc_syms_Az, 1, 1, 0);
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

        phitilde_p[index] = phitilde[index];
        Ax_p[index] = Ax[index];
        Ay_p[index] = Ay[index];
        Az_p[index] = Az[index];

        rho_star_p_p[index] = rho_star[index];
        tau_p_p[index] = tau[index];
        Stildex_p_p[index] = Stildex[index];
        Stildey_p_p[index] = Stildey[index];
        Stildez_p_p[index] = Stildez[index];

        phitilde_p_p[index] = phitilde[index];
        Ax_p_p[index] = Ax[index];
        Ay_p_p[index] = Ay[index];
        Az_p_p[index] = Az[index];
  }
}
