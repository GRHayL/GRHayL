#include "GRHayL.h"
#include "NRPyEOS_Tabulated.h"

double rel_err(const double a, const double b) {
  if( a == 0 ) return fabs(1 - a/b);
  if( b == 0 ) return fabs(1 - b/a);
  return 0.0;
}

double get_table_quantity(
      const int table_key,
      const int npoints,
      const int index ) {
  const double quantity = table_key*10 + index*9/npoints;

  switch (table_key) {
  case NRPyEOS_press_key:
    return quantity * log(10.0) + log(CGS_TO_CODE_PRESSURE);
  case NRPyEOS_eps_key:
    return quantity * log(10.0) + log(CGS_TO_CODE_ENERGY);
  case NRPyEOS_cs2_key:
    return quantity * CGS_TO_CODE_LENGTH*CGS_TO_CODE_LENGTH/CGS_TO_CODE_TIME/CGS_TO_CODE_TIME;
  case NRPyEOS_depsdT_key:
    return quantity * CGS_TO_CODE_ENERGY;
  case NRPyEOS_dPdrho_key:
    return quantity * CGS_TO_CODE_PRESSURE/CGS_TO_CODE_DENSITY;
  case NRPyEOS_dPdeps_key:
    return quantity * CGS_TO_CODE_PRESSURE/CGS_TO_CODE_ENERGY;
  default:
    return quantity;
  }
}

/*
 * (c) 2022 Leo Werneck
 */
int main(int argc, char **argv) {

  if( argc != 2 )
    grhayl_error("Correct usage is %s <table path>\n", argv[0]);

  grhayl_info("Beginning readtable unit test\n");

  // Step 1: Initialize the EOS struct
  eos_parameters eos;
  initialize_tabulated_functions(&eos);
  eos.tabulated_read_table_set_EOS_params(argv[1], &eos);

  if( eos.N_rho != 7 || eos.N_T != 5 || eos.N_Ye != 3 )
    grhayl_error("Table dimension error: expected 7 x 5 x 3, but got %d x %d x %d\n",
                 eos.N_rho, eos.N_T, eos.N_Ye);

  grhayl_info("Table dimensions: %d x %d x %d\n", eos.N_rho, eos.N_T, eos.N_Ye);

  // Step 2: Begin test
  const int npoints = eos.N_rho * eos.N_T * eos.N_Ye;
  for(int index=0; index<npoints; index++) {
    for(int var_key=0; var_key<NRPyEOS_ntablekeys; var_key++) {
      const double var       = get_table_quantity(var_key, npoints, index);
      const double table_var = eos.table_all[var_key + NRPyEOS_ntablekeys*index];
      if( rel_err(var, table_var) > 1e-15 )
        grhayl_error("Relative error exceeds round-off: %.15e vs. %.15e\n", var, table_var);
    }
  }

  // Step 4: Free memory
  eos.tabulated_free_memory(&eos);

  grhayl_info("Test finished successfully.\n");

  // All done!
  return 0;
}
