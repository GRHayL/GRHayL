#include "./EOS_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_free_memory(eos_parameters *restrict eos_params) {


 fprintf(stderr,"(NRPyEOS) *******************************\n");
 fprintf(stderr,"(NRPyEOS) Freeing up memory.\n");

  // Free memory allocated for the table
  free(eos_params->table_logrho);
  free(eos_params->table_logT);
  free(eos_params->table_Ye);
  free(eos_params->table_all);
  free(eos_params->table_eps);

 fprintf(stderr,"(NRPyEOS) All done!\n");
 fprintf(stderr,"(NRPyEOS) *******************************\n");
}
