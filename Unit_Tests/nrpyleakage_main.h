int main(int argc, char **argv) {

  if( argc != 3 ) {
    grhayl_info("Usage: %s <EOS table path> <test key>\n", argv[0]);
    grhayl_info("Available test keys:\n");
    grhayl_info("  0 : Generate data\n");
    grhayl_info("  1 : Run unit test\n");
    exit(1);
  }

  const char *tablepath  = argv[1];
  const int test_key     = atoi(argv[2]);

  if( test_key != 0 && test_key != 1 ) {
    grhayl_info("Available test keys:\n");
    grhayl_info("  0 : Generate data\n");
    grhayl_info("  1 : Run unit test\n");
    grhayl_error("Unsupported test key: %d\n", test_key);
  }

  const double W_max     = 10.0;
  const double rho_b_atm = 1e-12;
  const double rho_b_min = -1;
  const double rho_b_max = -1;
  const double Y_e_atm   = 0.5;
  const double Y_e_min   = -1;
  const double Y_e_max   = -1;
  const double T_atm     = 1e-2;
  const double T_min     = -1;
  const double T_max     = -1;

  eos_parameters eos;
  grhayl_initialize_tabulated_eos_functions_and_params(tablepath, W_max,
                                                rho_b_atm, rho_b_min, rho_b_max,
                                                Y_e_atm, Y_e_min, Y_e_max,
                                                T_atm, T_min, T_max, &eos);

  if( test_key ) {
    run_unit_test(&eos);
    grhayl_info("Test finished successfully\n");
  }
  else {
    generate_test_data(&eos);
    grhayl_info("Data generation finished successfully\n");
  }

  eos.tabulated_free_memory(&eos);

  return 0;
}
