#include "ghl.h"
#include "unit_tests.h"

#include <assert.h>

int main(int argc, char **argv) {

  if(argc != 3) {
    fprintf(stderr, "Usage is: ./data_gen <eos table> <input_file>\n");
    exit(1);
  }

  const char *tablepath = argv[1];
  const char *input_file = argv[2];

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const int main_routine = None;
  const int backup_routine[3] = { None, None, None };
  const bool evolve_entropy = true;
  const bool evolve_temperature = true;
  const bool calc_prims_guess = true;
  const double Psi6threshold = 1e100;
  const double W_max = 10.0;
  const double Lorenz_damping_factor = 0.0;

  const double rho_b_atm = 1e-12;
  const double rho_b_min = rho_b_atm;
  const double rho_b_max = 1e300; // IGM default
  const double Y_e_atm = 0.5;
  const double Y_e_min = 0.05;
  const double Y_e_max = Y_e_atm;
  const double T_atm = 1e-2;
  const double T_min = T_atm;
  const double T_max = 1e2;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  ghl_parameters params;
  ghl_initialize_params(
        main_routine, backup_routine, evolve_entropy, evolve_temperature,
        calc_prims_guess, Psi6threshold, W_max, Lorenz_damping_factor, &params);

  ghl_eos_parameters eos;
  ghl_initialize_tabulated_eos_functions_and_params(
        tablepath, rho_b_atm, rho_b_min, rho_b_max, Y_e_atm, Y_e_min, Y_e_max, T_atm,
        T_min, T_max, &eos);
  eos.root_finding_precision = 1e-10;

  FILE *fp = fopen_with_check(input_file, "rb");
  int nrho, nt, nye;
  fread(&nrho, sizeof(int), 1, fp);
  fread(&nt, sizeof(int), 1, fp);
  fread(&nye, sizeof(int), 1, fp);
  assert(nrho == 50);
  assert(nt == 50);
  assert(nye == 50);

  fclose(fp);

  return 0;
}
