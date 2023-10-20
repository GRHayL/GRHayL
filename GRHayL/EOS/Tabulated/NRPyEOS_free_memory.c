#include "nrpyeos_tabulated.h"

/// @brief Free memory allocated for EOS table quantities in GRHayL EOS parameter struct
/// @param eos Pointer to GRHayL EOS parameters struct
void NRPyEOS_free_memory(ghl_eos_parameters *restrict eos) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
#else
  ghl_info("*******************************\n");
  ghl_info("Freeing up memory.\n");

  // Free memory allocated for the table
  free(eos->table_logrho);
  free(eos->table_logT);
  free(eos->table_Y_e);
  free(eos->table_all);
  free(eos->table_eps);
  NRPyEOS_tabulated_free_beq_quantities(eos);

  ghl_info("All done!\n");
  ghl_info("*******************************\n");
#endif
}
