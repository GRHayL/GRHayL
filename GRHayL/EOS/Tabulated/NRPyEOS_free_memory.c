#include "nrpyeos_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_free_memory(eos_parameters *restrict eos_params) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
#else
  grhayl_info("*******************************\n");
  grhayl_info("Freeing up memory.\n");

  // Free memory allocated for the table
  free(eos_params->table_logrho);
  free(eos_params->table_logT);
  free(eos_params->table_Y_e);
  free(eos_params->table_all);
  free(eos_params->table_eps);

 grhayl_info("All done!\n");
 grhayl_info("*******************************\n");
#endif
}
