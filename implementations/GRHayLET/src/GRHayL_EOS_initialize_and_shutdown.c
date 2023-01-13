#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "GRHayLET.h"

eos_parameters *grhayl_eos;

// extern eos_parameters *grhayl_eos_parameters;

#ifdef __cplusplus
extern "C"
#endif
void GRHayLET_EOS_initialize(CCTK_ARGUMENTS) {

  DECLARE_CCTK_PARAMETERS;

  grhayl_eos = (eos_parameters *)malloc(sizeof(eos_parameters));
  initialize_tabulated_functions(grhayl_eos);
  grhayl_eos->tabulated_read_table_set_EOS_params(EOS_tablepath, grhayl_eos);
}

void GRHayLET_EOS_terminate(CCTK_ARGUMENTS) {
  free(grhayl_eos->table_all);
  free(grhayl_eos->table_logrho);
  free(grhayl_eos->table_logT);
  free(grhayl_eos->table_Ye);
  free(grhayl_eos->table_eps);
  free(grhayl_eos);
}
