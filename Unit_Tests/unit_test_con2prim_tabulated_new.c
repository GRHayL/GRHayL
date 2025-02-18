#include "ghl.h"
#include "ghl_unit_tests.h"

#include <assert.h>

extern FILE *flog;

void print_metric(FILE *fp, const ghl_metric_quantities *restrict metric) {
  if(!fp) {
    return;
  }
  fputs("Metric:\n", fp);
  fprintf(fp, "  lapse    : %.15e\n", metric->lapse);
  fprintf(fp, "  betax    : %.15e\n", metric->betaU[0]);
  fprintf(fp, "  betay    : %.15e\n", metric->betaU[1]);
  fprintf(fp, "  betaz    : %.15e\n", metric->betaU[2]);
  fprintf(fp, "  gammaDD00: %.15e\n", metric->gammaDD[0][0]);
  fprintf(fp, "  gammaDD01: %.15e\n", metric->gammaDD[0][1]);
  fprintf(fp, "  gammaDD02: %.15e\n", metric->gammaDD[0][2]);
  fprintf(fp, "  gammaDD11: %.15e\n", metric->gammaDD[1][1]);
  fprintf(fp, "  gammaDD12: %.15e\n", metric->gammaDD[1][2]);
  fprintf(fp, "  gammaDD22: %.15e\n", metric->gammaDD[2][2]);
  fprintf(fp, "  gammaUU00: %.15e\n", metric->gammaUU[0][0]);
  fprintf(fp, "  gammaUU01: %.15e\n", metric->gammaUU[0][1]);
  fprintf(fp, "  gammaUU02: %.15e\n", metric->gammaUU[0][2]);
  fprintf(fp, "  gammaUU11: %.15e\n", metric->gammaUU[1][1]);
  fprintf(fp, "  gammaUU12: %.15e\n", metric->gammaUU[1][2]);
  fprintf(fp, "  gammaUU22: %.15e\n", metric->gammaUU[2][2]);
}

void print_4metric(FILE *fp, const ghl_ADM_aux_quantities *restrict metric) {
  if(!fp) {
    return;
  }
  fputs("Four Metric:\n", fp);
  fprintf(fp, "  g4DD00: %.15e\n", metric->g4DD[0][0]);
  fprintf(fp, "  g4DD01: %.15e\n", metric->g4DD[0][1]);
  fprintf(fp, "  g4DD02: %.15e\n", metric->g4DD[0][2]);
  fprintf(fp, "  g4DD03: %.15e\n", metric->g4DD[0][3]);
  fprintf(fp, "  g4DD10: %.15e\n", metric->g4DD[1][0]);
  fprintf(fp, "  g4DD11: %.15e\n", metric->g4DD[1][1]);
  fprintf(fp, "  g4DD12: %.15e\n", metric->g4DD[1][2]);
  fprintf(fp, "  g4DD13: %.15e\n", metric->g4DD[1][3]);
  fprintf(fp, "  g4DD20: %.15e\n", metric->g4DD[2][0]);
  fprintf(fp, "  g4DD21: %.15e\n", metric->g4DD[2][1]);
  fprintf(fp, "  g4DD22: %.15e\n", metric->g4DD[2][2]);
  fprintf(fp, "  g4DD23: %.15e\n", metric->g4DD[2][3]);
  fprintf(fp, "  g4DD30: %.15e\n", metric->g4DD[3][0]);
  fprintf(fp, "  g4DD31: %.15e\n", metric->g4DD[3][1]);
  fprintf(fp, "  g4DD32: %.15e\n", metric->g4DD[3][2]);
  fprintf(fp, "  g4DD33: %.15e\n", metric->g4DD[3][3]);
  fprintf(fp, "  g4UU00: %.15e\n", metric->g4UU[0][0]);
  fprintf(fp, "  g4UU01: %.15e\n", metric->g4UU[0][1]);
  fprintf(fp, "  g4UU02: %.15e\n", metric->g4UU[0][2]);
  fprintf(fp, "  g4UU03: %.15e\n", metric->g4UU[0][3]);
  fprintf(fp, "  g4UU10: %.15e\n", metric->g4UU[1][0]);
  fprintf(fp, "  g4UU11: %.15e\n", metric->g4UU[1][1]);
  fprintf(fp, "  g4UU12: %.15e\n", metric->g4UU[1][2]);
  fprintf(fp, "  g4UU13: %.15e\n", metric->g4UU[1][3]);
  fprintf(fp, "  g4UU20: %.15e\n", metric->g4UU[2][0]);
  fprintf(fp, "  g4UU21: %.15e\n", metric->g4UU[2][1]);
  fprintf(fp, "  g4UU22: %.15e\n", metric->g4UU[2][2]);
  fprintf(fp, "  g4UU23: %.15e\n", metric->g4UU[2][3]);
  fprintf(fp, "  g4UU30: %.15e\n", metric->g4UU[3][0]);
  fprintf(fp, "  g4UU31: %.15e\n", metric->g4UU[3][1]);
  fprintf(fp, "  g4UU32: %.15e\n", metric->g4UU[3][2]);
  fprintf(fp, "  g4UU33: %.15e\n", metric->g4UU[3][3]);
}

void print_cons(FILE *fp, const ghl_conservative_quantities *restrict cons) {
  if(!fp) {
    return;
  }
  fputs("Conservatives:\n", fp);
  fprintf(fp, "  rho    : %.15e\n", cons->rho);
  fprintf(fp, "  tau    : %.15e\n", cons->tau);
  fprintf(fp, "  SD0    : %.15e\n", cons->SD[0]);
  fprintf(fp, "  SD1    : %.15e\n", cons->SD[1]);
  fprintf(fp, "  SD2    : %.15e\n", cons->SD[2]);
  fprintf(fp, "  Y_e    : %.15e\n", cons->Y_e);
  fprintf(fp, "  entropy: %.15e\n", cons->entropy);
}

void print_prims(FILE *fp, const ghl_primitive_quantities *restrict prims) {
  if(!fp) {
    return;
  }
  fputs("Primitives:\n", fp);
  fprintf(fp, "  rho        : %.15e\n", prims->rho);
  fprintf(fp, "  temperature: %.15e\n", prims->temperature);
  fprintf(fp, "  Y_e        : %.15e\n", prims->Y_e);
  fprintf(fp, "  press      : %.15e\n", prims->press);
  fprintf(fp, "  eps        : %.15e\n", prims->eps);
  fprintf(fp, "  entropy    : %.15e\n", prims->entropy);
  fprintf(fp, "  vU[0]      : %.15e\n", prims->vU[0]);
  fprintf(fp, "  vU[1]      : %.15e\n", prims->vU[1]);
  fprintf(fp, "  vU[2]      : %.15e\n", prims->vU[2]);
  fprintf(fp, "  BU[0]      : %.15e\n", prims->BU[0]);
  fprintf(fp, "  BU[1]      : %.15e\n", prims->BU[1]);
  fprintf(fp, "  BU[2]      : %.15e\n", prims->BU[2]);
}

int main(int argc, char **argv) {

  if(argc != 3) {
    fprintf(stderr, "Usage is: ./data_gen <eos table> <input_file>\n");
    exit(1);
  }

  const char *tablepath  = argv[1];
  const char *input_file = argv[2];

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
        main_routine,
        backup_routine,
        evolve_entropy,
        evolve_temperature,
        calc_prims_guess,
        Psi6threshold,
        W_max,
        Lorenz_damping_factor,
        &params);

  ghl_eos_parameters eos;
  ghl_initialize_tabulated_eos_functions_and_params(
        tablepath,
        rho_b_atm,
        rho_b_min,
        rho_b_max,
        Y_e_atm,
        Y_e_min,
        Y_e_max,
        T_atm,
        T_min,
        T_max,
        &eos);
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
  assert(nrho == 30);
  assert(nt == 30);
  assert(nye == 30);

  const ghl_con2prim_method_t methods[4]
        = { Palenzuela1D, Newman1D, Palenzuela1D_entropy, Newman1D_entropy };
  const char *methodnames[4]
        = { "Palenzuela1D", "Newman1D", "Palenzuela1D_entropy", "Newman1D_entropy" };
  int fails[4] = { 0, 0, 0, 0 };

  // flog = fopen_with_check("log_ghl.txt", "w");
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

        ghl_con2prim_diagnostics diagnostics;
        ghl_initialize_diagnostics(&diagnostics);

        if(flog) {
          print_metric(flog, &ADM_metric);
          print_4metric(flog, &AUX_metric);
        }

        for(int n = 0; n < 4; n++) {
          params.main_routine = methods[n];

          ghl_conservative_quantities cons_undens = { 0 };
          ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);

          ghl_primitive_quantities prims = { 0 };
          ghl_guess_primitives(&eos, &ADM_metric, &cons, &prims);
          prims.BU[0] = prims_orig.BU[0];
          prims.BU[1] = prims_orig.BU[1];
          prims.BU[2] = prims_orig.BU[2];

          if(flog) {
            fprintf(flog, "%02d, %02d, %02d --- %s:\n", i, j, k, methodnames[n]);
            print_cons(flog, &cons_undens);
            print_prims(flog, &prims);
          }

          if(ghl_con2prim_tabulated_select_method(
                   params.main_routine,
                   &params,
                   &eos,
                   &ADM_metric,
                   &AUX_metric,
                   &cons_undens,
                   &prims,
                   &diagnostics)) {
            fails[n]++;
            if(flog) {
              fputs("Failed!\n", flog);
            }
          }

          if(flog) {
            fputs("-----------\n", flog);
          }
        }
      }
    }
  }

  fclose(fp);
  if(flog) {
    fclose(flog);
  }

  const int npts = nrho * nt * nye;
  ghl_info("Failure rates:\n");
  for(int n = 0; n < 4; n++) {
    const char *name    = methodnames[n];
    const int failcount = fails[n];
    const double pct    = ((double)failcount) / ((double)npts) * 100;
    ghl_info("    %-20s : %06d/%06d : %5.1lf%%\n", name, failcount, npts, pct);
  }

  return 0;
}
