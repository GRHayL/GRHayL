#include "ghl.h"
#include "ghl_unit_tests.h"

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
  if(fread(&nrho, sizeof(int), 1, fp) != 1) {
    fclose(fp);
    ghl_error("Failed to read nrho\n");
  }
  if(fread(&nt, sizeof(int), 1, fp) != 1) {
    fclose(fp);
    ghl_error("Failed to read nt\n");
  }
  if(fread(&nye, sizeof(int), 1, fp) != 1) {
    fclose(fp);
    ghl_error("Failed to read nye\n");
  }
  assert(nrho == 50);
  assert(nt == 50);
  assert(nye == 50);

  const ghl_con2prim_method_t methods[4]
        = { Palenzuela1D, Newman1D, Palenzuela1D_entropy, Newman1D_entropy };
  const char *methodnames[4]
        = { "Palenzuela1D", "Newman1D", "Palenzuela1D_entropy", "Newman1D_entropy" };
  int fails[4] = { 0, 0, 0, 0 };

  for(int k = 0; k < nye; k++) {
    for(int j = 0; j < nt; j++) {
      for(int i = 0; i < nrho; i++) {
        ghl_metric_quantities ADM_metric;
        if(fread(&ADM_metric, sizeof(ADM_metric), 1, fp) != 1) {
          fclose(fp);
          ghl_error("Failed to read ADM metric at index %d, %d, %d\n", i, j, k);
        }

        ghl_ADM_aux_quantities AUX_metric;
        if(fread(&AUX_metric, sizeof(AUX_metric), 1, fp) != 1) {
          fclose(fp);
          ghl_error("Failed to read AUX metric at index %d, %d, %d\n", i, j, k);
        }

        ghl_primitive_quantities prims_orig;
        if(fread(&prims_orig, sizeof(prims_orig), 1, fp) != 1) {
          fclose(fp);
          ghl_error("Failed to read prims at index %d, %d, %d\n", i, j, k);
        }

        ghl_conservative_quantities cons;
        if(fread(&cons, sizeof(cons), 1, fp) != 1) {
          fclose(fp);
          ghl_error("Failed to read cons at index %d, %d, %d\n", i, j, k);
        }

        ghl_conservative_quantities cons_undens;
        ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);

        ghl_con2prim_diagnostics diagnostics;
        ghl_initialize_diagnostics(&diagnostics);

        for(int n = 0; n < 4; n++) {
          params.main_routine = methods[n];

          ghl_primitive_quantities prims;
          if(ghl_con2prim_tabulated_multi_method(
                   &params, &eos, &ADM_metric, &AUX_metric, &cons_undens, &prims,
                   &diagnostics)) {
            fails[n]++;
          }
        }
      }
    }
  }

  fclose(fp);

  const int npts = nrho * nt * nye;
  ghl_info("Failure rates:\n");
  for(int n = 0; n < 4; n++) {
    const char *name = methodnames[n];
    const int failcount = fails[n];
    const double pct = ((double)failcount) / ((double)npts) * 100;
    ghl_info("    %-20s : %06d/%06d : %5.1lf%%\n", name, failcount, npts, pct);
  }

  return 0;
}
