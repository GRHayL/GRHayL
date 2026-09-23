#include "ghl_nrpyeos_tabulated.h"
#include "ghl_con2prim.h"

/*
 * (c) 2022 Leo Werneck
 */
// In HDF5-enabled builds, the outer struct must be successful, failed-empty,
// previously cleaned, or properly zero-initialized. Cleanup is repeatable and
// leaves owned pointers null. Disabled-HDF5 initialization needs no cleanup.
void NRPyEOS_free_memory(ghl_eos_parameters *restrict eos) {
#ifdef GHL_DISABLE_HDF5
  GHL_HDF5_ERROR_IF_USED;
#else
  ghl_info("*******************************\n");
  ghl_info("Freeing up memory.\n");

  // Free memory allocated for the table
  free(eos->table_logrho);
  eos->table_logrho = NULL;
  free(eos->table_logT);
  eos->table_logT = NULL;
  free(eos->table_Y_e);
  eos->table_Y_e = NULL;
  free(eos->table_all);
  eos->table_all = NULL;
  free(eos->table_eps);
  eos->table_eps = NULL;
  free(eos->table_logh);
  eos->table_logh = NULL;
  NRPyEOS_tabulated_free_beq_quantities(eos);
  ghl_c2p_nn_free(eos->c2p_nn);
  eos->c2p_nn = NULL;

  ghl_info("All done!\n");
  ghl_info("*******************************\n");
#endif
}
