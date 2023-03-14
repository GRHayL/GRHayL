//--------------------------------------------------------------------------
// Register with the time stepper
// (MoL thorn, found in arrangements/CactusBase/MoL)
// To understand this, read documentation in arrangements/CactusBase/MoL/doc
//--------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "Symmetry.h"

void GRHayL_IGM_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_GRHayL_IGM_RegisterVars;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, var, rhs;

  //***********************************************
  // Register evolution & RHS gridfunction variables

  /* Ax and Ax_rhs */
  var = CCTK_VarIndex("GRHayL_IGM::Ax");
  rhs = CCTK_VarIndex("GRHayL_IGM::Ax_rhs");
  ierr += MoLRegisterEvolved(var, rhs);

  /* Ay and Ay_rhs */
  var = CCTK_VarIndex("GRHayL_IGM::Ay");
  rhs = CCTK_VarIndex("GRHayL_IGM::Ay_rhs");
  ierr += MoLRegisterEvolved(var, rhs);

  /* Az and Az_rhs */
  var = CCTK_VarIndex("GRHayL_IGM::Az");
  rhs = CCTK_VarIndex("GRHayL_IGM::Az_rhs");
  ierr += MoLRegisterEvolved(var, rhs);

  /* psi6phi and psi6phi_rhs */
  var = CCTK_VarIndex("GRHayL_IGM::psi6phi");
  rhs = CCTK_VarIndex("GRHayL_IGM::psi6phi_rhs");
  ierr += MoLRegisterEvolved(var, rhs);

  /* ALL OTHER EVOLVED VARIABLES (rho_star,tau,Stilde_x,Stilde_y,Stilde_z) */
  var = CCTK_GroupIndex("GRHayL_IGM::grmhd_conservatives");
  rhs = CCTK_GroupIndex("GRHayL_IGM::grmhd_conservatives_rhs");
  ierr += MoLRegisterEvolvedGroup(var, rhs);

  if (ierr) CCTK_ERROR("Problems registering with MoL");
  //***********************************************

  //***********************************************
  // Next register ADMBase variables needed by
  //    GRHayL_IGM as SaveAndRestore, so that
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
