#include <string.h>
#include "unit_tests.h"
#include <string.h>

void init_metric_from_string( const char *lapse_shift_string_in,
                              const char *metric_string_in,
                              ghl_metric_quantities *ADM_metric ) {

  char *token;
  char lapse_shift_string[strlen(lapse_shift_string_in)+1];
  char metric_string[strlen(metric_string_in)+1];
  strcpy(&lapse_shift_string[0], lapse_shift_string_in);
  strcpy(&metric_string     [0], metric_string_in     );

  // Lapse/shift string is: "lapse shift0 shift1 shift2"
  token = strtok(lapse_shift_string, " "); const double lapse  = strtod(token, NULL);
  token = strtok(NULL              , " "); const double betaU0 = strtod(token, NULL);
  token = strtok(NULL              , " "); const double betaU1 = strtod(token, NULL);
  token = strtok(NULL              , " "); const double betaU2 = strtod(token, NULL);

  // Metric string is: "gxx gxy gxz gyy gyz gzz"
  token = strtok(metric_string, " "); const double gxx = strtod(token, NULL);
  token = strtok(NULL         , " "); const double gxy = strtod(token, NULL);
  token = strtok(NULL         , " "); const double gxz = strtod(token, NULL);
  token = strtok(NULL         , " "); const double gyy = strtod(token, NULL);
  token = strtok(NULL         , " "); const double gyz = strtod(token, NULL);
  token = strtok(NULL         , " "); const double gzz = strtod(token, NULL);

  ghl_initialize_metric(lapse, betaU0, betaU1, betaU2,
                        gxx, gxy, gxz, gyy, gyz, gzz,
                        ADM_metric);
}

void init_conservs_from_string( const char *conservs_string_in,
                                ghl_conservative_quantities *cons ) {

  char *token;
  char conservs_string[strlen(conservs_string_in)+1];
  strcpy(&conservs_string[0], conservs_string_in);

  // Conservs string string is: "rho, tau, SD0, SD1, SD2, ent, Y_e"
  token = strtok(conservs_string, " "); const double rho = strtod(token, NULL);
  token = strtok(NULL           , " "); const double tau = strtod(token, NULL);
  token = strtok(NULL           , " "); const double SD0 = strtod(token, NULL);
  token = strtok(NULL           , " "); const double SD1 = strtod(token, NULL);
  token = strtok(NULL           , " "); const double SD2 = strtod(token, NULL);
  token = strtok(NULL           , " "); const double ent = strtod(token, NULL);
  token = strtok(NULL           , " "); const double Y_e = strtod(token, NULL);

  ghl_initialize_conservatives(rho, tau, SD0, SD1, SD2, ent, Y_e, cons);
}

void init_Bfield_from_string( const char *B_string_in,
                              ghl_primitive_quantities *prims ) {

  char *token;
  char B_string[strlen(B_string_in)+1];
  strcpy(&B_string[0], B_string_in);

  // B field string is: "BU0 BU1 BU2"
  token = strtok(B_string, " "); prims->BU[0] = strtod(token, NULL);
  token = strtok(NULL    , " "); prims->BU[1] = strtod(token, NULL);
  token = strtok(NULL    , " "); prims->BU[2] = strtod(token, NULL);
}

int main(int argc, char **argv) {

  char *tablepath = NULL;
  if( argc != 2 ) {
    ghl_info("Usage: %s <EOS table path>\n", argv[0]);
    exit(1);
  }
  else if( argc == 2 )
    tablepath = argv[1];

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const ghl_con2prim_method_t main_routine       = Newman1D;
  const ghl_con2prim_method_t backup_routines[3] = {None,None,None};
  const bool calc_prims_guess                = false;
  const bool evolve_entropy                  = false;
  const bool evolve_temperature              = true;
  const bool use_cupp_fix                    = true;
  const double Psi6threshold                 = 60; //Taken from magnetizedTOV.par
  const double Lorenz_damping_factor         = 0.0;

  const double W_max     = 10.0;
  const double rho_b_atm = 1e-12;
  const double rho_b_min = rho_b_atm;
  const double rho_b_max = 1e300;
  const double Y_e_atm   = 0.5;
  const double Y_e_min   = 0.05;
  const double Y_e_max   = Y_e_atm;
  const double T_atm     = 1e-2;
  const double T_min     = T_atm;
  const double T_max     = 1e2;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  ghl_parameters params;
  ghl_initialize_params(main_routine,
                    backup_routines,
                    evolve_entropy,
                    evolve_temperature,
                    calc_prims_guess,
                    Psi6threshold,
                    use_cupp_fix,
                    Lorenz_damping_factor,
                    &params);

  ghl_eos_parameters eos;
  ghl_initialize_tabulated_eos_functions_and_params(tablepath, W_max,
                                                    rho_b_atm, rho_b_min, rho_b_max,
                                                    Y_e_atm, Y_e_min, Y_e_max,
                                                    T_atm, T_min, T_max, &eos);
  eos.root_finding_precision=1e-10;

  const char *lapse_shift_string = "8.787479e-01 8.657114e-03 -2.066025e-02 2.384216e-04";
  const char *metric_string      = "1.315665e+00 -6.933605e-02 5.229179e-03 1.358808e+00 -8.966700e-03 1.272508e+00";
  const char *conservs_string    = "5.132187e-10 -2.636040e-05 2.826288e-06 2.442308e-06 1.950662e-06 -1.984531e+06 3.046700e-10";
  const char *B_string           = "-2.447595e-03 8.354818e-03 -4.903698e-03";

  ghl_metric_quantities ADM_metric;
  init_metric_from_string(lapse_shift_string, metric_string, &ADM_metric);

  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

  ghl_conservative_quantities cons;
  init_conservs_from_string(conservs_string, &cons);

  ghl_conservative_quantities cons_undens;
  ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);

  ghl_primitive_quantities prims;
  init_Bfield_from_string(B_string, &prims);

  ghl_con2prim_diagnostics diagnostics;
  ghl_initialize_diagnostics(&diagnostics);

  int check = ghl_con2prim_tabulated_multi_method(&params, &eos, &ADM_metric, &metric_aux, &cons_undens, &prims, &diagnostics);

  ghl_debug_print_cons(&cons);

  if( check )
    ghl_info("Con2Prim failed!\n");
  else {
    ghl_info("Con2Prim succeeded!\n");
    if( isnan(prims.rho*prims.press*prims.eps*prims.vU[0]*prims.vU[1]*prims.vU[2]*prims.entropy*prims.Y_e*prims.temperature) ) {
      ghl_info("Con2Prim succeeded, but found NAN in prims!!!\n");
      ghl_debug_print_prims(&prims);
    }
  }

  ghl_tabulated_free_memory(&eos);

  return 0;
}
