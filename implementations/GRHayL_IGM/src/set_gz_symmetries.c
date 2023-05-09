//-------------------------------------------------
// Stuff to run right after initial data is set up
//-------------------------------------------------

#include "IGM.h"
#include "Symmetry.h"

void GRHayL_IGM_set_gz_symmetries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayL_IGM_set_gz_symmetries;
  DECLARE_CCTK_PARAMETERS;

  //For emfields, we assume that you've set Bx, By, Bz (the UN-tilded B^i's)
  // or Ax, Ay, Az (if using constrained transport scheme of Del Zanna)

  if(CCTK_EQUALS(Symmetry,"equatorial")) {
    // Set up symmetry ghostzones on Bx, By, Bz, and their staggered variants.
    const double gridfunc_syms_Bx[3] = { -1,  1, -Sym_Bz};
    GRHayL_IGM_set_symmetry_gzs_staggered(cctkGH, x, y, z, Bx_center,  gridfunc_syms_Bx, 0, 0, 0);
    GRHayL_IGM_set_symmetry_gzs_staggered(cctkGH, x, y, z, Bx_stagger, gridfunc_syms_Bx, 1, 0, 0);
    const double gridfunc_syms_By[3] = {  1, -1, -Sym_Bz};
    GRHayL_IGM_set_symmetry_gzs_staggered(cctkGH, x, y, z, By_center,  gridfunc_syms_By, 0, 0, 0);
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
}
