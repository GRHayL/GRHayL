#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

void NRPyLeakageET_InitSym(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int sym[3] = {1,1,1};
  SetCartSymGN(cctkGH,sym,"NRPyLeakageET::NRPyLeakageET_opacities");
  SetCartSymGN(cctkGH,sym,"NRPyLeakageET::NRPyLeakageET_optical_depths");
}
