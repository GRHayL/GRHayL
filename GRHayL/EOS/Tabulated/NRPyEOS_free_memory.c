#include "nrpyeos_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_free_memory(eos_parameters *restrict eos_params) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
#else
  ghl_info("*******************************\n");
  ghl_info("Freeing up memory.\n");

  // Free memory allocated for the table
  free(eos_params->table_logrho);
  free(eos_params->table_logT);
  free(eos_params->table_Y_e);
  free(eos_params->table_all);
  free(eos_params->table_eps);

 ghl_info("All done!\n");
 ghl_info("*******************************\n");
#endif
}
