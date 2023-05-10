//--------------------------------------------------------------------------
// Register with the time stepper
// (MoL thorn, found in arrangements/CactusBase/MoL)
// To understand this, read documentation in arrangements/CactusBase/MoL/doc
//--------------------------------------------------------------------------

#include "IGH.h"
#include "Symmetry.h"

void GRHayL_IGH_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_GRHayL_IGH_RegisterVars;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, var, rhs;

  //***********************************************
  // Register evolution & RHS gridfunction variables
  var = CCTK_GroupIndex("GRHayL_IGH::grmhd_conservatives");
  rhs = CCTK_GroupIndex("GRHayL_IGH::grmhd_conservatives_rhs");
  ierr += MoLRegisterEvolvedGroup(var, rhs);

  if (ierr) CCTK_ERROR("Problems registering with MoL");
  //***********************************************

  //***********************************************
  // Next register ADMBase variables needed by
  //    GRHayL_IGH as SaveAndRestore, so that
  //    they are not set to NaN at the start of
  //    each timestep (requiring that they be
  //    e.g., recomputed from BSSN variables
  //    in the BSSN solver, like Baikal or
  //    ML_BSSN)
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("admbase::lapse"));
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("admbase::shift"));
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("admbase::metric"));
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("admbase::curv"));
  if (ierr) CCTK_ERROR("Problems registering with MoLRegisterSaveAndRestoreGroup");
}
