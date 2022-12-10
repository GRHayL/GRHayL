#include "./NRPy_basic_defines.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_free_memory(NRPyEOS_params *restrict eos_params) {


 fprintf(stderr,"(NRPyEOS) *******************************\n");
 fprintf(stderr,"(NRPyEOS) Freeing up memory.\n");

  // Free memory allocated for the table
  free(eos_params->logrho);
  free(eos_params->logtemp);
  free(eos_params->yes);
  free(eos_params->alltables);
  free(eos_params->epstable);

 fprintf(stderr,"(NRPyEOS) All done!\n");
 fprintf(stderr,"(NRPyEOS) *******************************\n");
}
