//--------------------------------------------------------------------------
// Register with the time stepper
// (MoL thorn, found in arrangements/CactusBase/MoL)
// To understand this, read documentation in arrangements/CactusBase/MoL/doc
//--------------------------------------------------------------------------

#include "GRHayLMHD.h"
#include "Symmetry.h"

void GRHayLMHD_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_GRHayLMHD_RegisterVars;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, var, rhs;

  //***********************************************
  // Register evolution & RHS gridfunction variables

  /* Ax and Ax_rhs */
  var = CCTK_VarIndex("GRHayLMHD::Ax");
  rhs = CCTK_VarIndex("GRHayLMHD::Ax_rhs");
  ierr += MoLRegisterEvolved(var, rhs);

  /* Ay and Ay_rhs */
  var = CCTK_VarIndex("GRHayLMHD::Ay");
  rhs = CCTK_VarIndex("GRHayLMHD::Ay_rhs");
  ierr += MoLRegisterEvolved(var, rhs);

  /* Az and Az_rhs */
  var = CCTK_VarIndex("GRHayLMHD::Az");
  rhs = CCTK_VarIndex("GRHayLMHD::Az_rhs");
  ierr += MoLRegisterEvolved(var, rhs);

  /* phitilde and phitilde_rhs */
  var = CCTK_VarIndex("GRHayLMHD::phitilde");
  rhs = CCTK_VarIndex("GRHayLMHD::phitilde_rhs");
  ierr += MoLRegisterEvolved(var, rhs);

  /* ALL OTHER EVOLVED VARIABLES (rho_star,tau,Stilde_x,Stilde_y,Stilde_z) */
  var = CCTK_GroupIndex("GRHayLMHD::grmhd_conservatives");
  rhs = CCTK_GroupIndex("GRHayLMHD::grmhd_conservatives_rhs");
  ierr += MoLRegisterEvolvedGroup(var, rhs);

  if (ierr) CCTK_ERROR("Problems registering with MoL");
  //***********************************************

  //***********************************************
  // Next register ADMBase variables needed by
  //    GRHayLMHD as SaveAndRestore, so that
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
