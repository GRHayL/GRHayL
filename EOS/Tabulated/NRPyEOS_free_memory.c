#include "GRHayL_EOS_Tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_free_memory(eos_parameters *restrict eos_params) {


 printf("(GRHayL - EOS) *******************************\n");
 printf("(GRHayL - EOS) Freeing up memory.\n");

  // Free memory allocated for the table
  free(eos_params->table_logrho);
  free(eos_params->table_logT);
  free(eos_params->table_Ye);
  free(eos_params->table_all);
  free(eos_params->table_eps);

 printf("(GRHayL - EOS) All done!\n");
 printf("(GRHayL - EOS) *******************************\n");
}
