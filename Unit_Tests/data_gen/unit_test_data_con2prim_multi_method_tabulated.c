#include "ghl.h"
#include "ghl_unit_tests.h"

static inline double random_perturbation(const double magnitude) {
  return randf(1, 1) * 1e-14;
}

static void ghl_perturb_conservs(ghl_conservative_quantities *restrict cons) {
  cons->rho     *= 1.0 + random_perturbation(1e-14);
  cons->tau     *= 1.0 + random_perturbation(1e-14);
  cons->SD[0]   *= 1.0 + random_perturbation(1e-14);
  cons->SD[1]   *= 1.0 + random_perturbation(1e-14);
  cons->SD[2]   *= 1.0 + random_perturbation(1e-14);
  cons->Y_e     *= 1.0 + random_perturbation(1e-14);
  cons->entropy *= 1.0 + random_perturbation(1e-14);
}

int main(int argc, char **argv) {

  if(argc != 2) {
    fprintf(stderr, "Usage is: ./data_gen <eos table>\n");
    exit(1);
  }

  const char *tablepath = argv[1];
  {
    FILE *fp = fopen(tablepath, "r");
    if(!fp) {
      fprintf(stderr, "Could not open table file %s\n", tablepath);
      exit(2);
    }
    fclose(fp);
  }

  // These variables set up the tested range of values and number of sampling points.
  // number of sampling points in density and pressure
  const int nrho = 30;
  const int nt   = 30;
  const int nye  = 30;

  const double test_rho_min = 1e-12; // Minimum input density
  const double test_rho_max = 1e-3;  // Maximum input density
  const double test_T_min   = 1e-2;  // Minimum input temperature
  const double test_T_max   = 1e+2;  // Maximum input temperature
  const double test_ye_min  = 0.075; // Minimum input electron fraction
  const double test_ye_max  = 0.500; // Maximum input electron fraction

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const int main_routine             = None;
  const int backup_routine[3]        = { None, None, None };
  const bool evolve_entropy          = true;
  const bool evolve_temperature      = true;
  const bool calc_prims_guess        = true;
  const double Psi6threshold         = 1e100;
  const double W_max                 = 10.0;
  const double Lorenz_damping_factor = 0.0;

  const double rho_b_atm = 1e-12;
  const double rho_b_min = rho_b_atm;
  const double rho_b_max = 1e300; // IGM default
  const double Y_e_atm   = 0.5;
  const double Y_e_min   = 0.05;
  const double Y_e_max   = Y_e_atm;
  const double T_atm     = 1e-2;
  const double T_min     = T_atm;
  const double T_max     = 1e2;

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

  const double lr_min = log(test_rho_min);
  const double lr_max = log(test_rho_max);
  const double dlr    = (lr_max - lr_min) / (nrho - 1);

  const double lt_min = log(test_T_min);
  const double lt_max = log(test_T_max);
  const double dlt    = (lt_max - lt_min) / (nt - 1);

  const double ye_min = test_ye_min;
  const double ye_max = test_ye_max;
  const double dye    = (ye_max - ye_min) / (nye - 1);

  const char *output_filename = "ghl_unit_test_con2prim_tabulated.bin";
  FILE *fp                    = fopen(output_filename, "wb");
  if(!fp) {
    fprintf(stderr, "Could not open output file %s\n", output_filename);
    exit(3);
  }
  ghl_info("Generating input file for unit test...\n");

  fwrite(&nrho, sizeof(int), 1, fp);
  fwrite(&nt, sizeof(int), 1, fp);
  fwrite(&nye, sizeof(int), 1, fp);
  for(int k = 0; k < nye; k++) {
    const double xye = ye_min + dye * k;
    for(int j = 0; j < nt; j++) {
      const double xtemperature = exp(lt_min + dlt * j);
      for(int i = 0; i < nrho; i++) {
        const double xrho = exp(lr_min + dlr * i);

        // Initialize primitives
        ghl_primitive_quantities prims;
        ghl_prims_with_random_velocities_and_magnetic_fields(
              &eos, xrho, xye, xtemperature, &prims);

        // Initialize a random metric
        ghl_metric_quantities ADM_metric;
        ghl_random_metric(&ADM_metric);

        // Initialize auxiliary metric quantities
        ghl_ADM_aux_quantities AUX_metric;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &AUX_metric);

        bool speed_limited;
        ghl_limit_v_and_compute_u0(&params, &ADM_metric, &prims, &speed_limited);

        // Compute conservatives
        ghl_conservative_quantities cons;
        ghl_compute_conservs(&ADM_metric, &AUX_metric, &prims, &cons);

        // Now get perturbed conservatives
        // ghl_conservative_quantities cons_perturbed = cons;
        // ghl_perturb_conservs(&cons_perturbed);

        // Metric, primitives, and conservatives to file
        fwrite(&ADM_metric, sizeof(ADM_metric), 1, fp);
        fwrite(&AUX_metric, sizeof(AUX_metric), 1, fp);
        fwrite(&prims, sizeof(prims), 1, fp);
        fwrite(&cons, sizeof(cons), 1, fp);
        // fwrite(&cons_perturbed, sizeof(cons_perturbed), 1, fp);
      }
    }
  }
  fclose(fp);
  ghl_info("Done!\n");

  return 0;
}
