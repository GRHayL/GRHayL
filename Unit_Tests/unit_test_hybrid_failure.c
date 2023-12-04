#include "unit_tests.h"
int main(int argc, char **argv) {

  const int arraylength = 3;

  const double poison = 0.0/0.0;

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const int backup_routine[3] = {Noble2D,Noble2D,Noble2D};
  const bool evolve_entropy = false;
  const bool evolve_temperature = false;
  const bool calc_prims_guess = true;
  const double Psi6threshold = 1e100;

  const int neos = 1;
  const double W_max = 10.0;
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300;
  const double Gamma_th = 2.0;
  const double rho_ppoly[1] = {0.0};
  const double Gamma_ppoly[1] = {2.0};
  const double k_ppoly0 = 1.0;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  ghl_parameters params;
  ghl_initialize_params(
        Noble2D, backup_routine, evolve_entropy, evolve_temperature, calc_prims_guess,
        Psi6threshold, W_max, 0.0, &params);

  ghl_eos_parameters eos;
  ghl_initialize_hybrid_eos_functions_and_params(
        rho_b_min, rho_b_min, rho_b_max,
        neos, rho_ppoly, Gamma_ppoly,
        k_ppoly0, Gamma_th, &eos);

  double *lapse = (double*) malloc(sizeof(double)*arraylength);
  double *betax = (double*) malloc(sizeof(double)*arraylength);
  double *betay = (double*) malloc(sizeof(double)*arraylength);
  double *betaz = (double*) malloc(sizeof(double)*arraylength);

  double *gxx = (double*) malloc(sizeof(double)*arraylength);
  double *gxy = (double*) malloc(sizeof(double)*arraylength);
  double *gxz = (double*) malloc(sizeof(double)*arraylength);
  double *gyy = (double*) malloc(sizeof(double)*arraylength);
  double *gyz = (double*) malloc(sizeof(double)*arraylength);
  double *gzz = (double*) malloc(sizeof(double)*arraylength);

  double *Bx = (double*) malloc(sizeof(double)*arraylength);
  double *By = (double*) malloc(sizeof(double)*arraylength);
  double *Bz = (double*) malloc(sizeof(double)*arraylength);

  double *rho_star = (double*) malloc(sizeof(double)*arraylength);
  double *tau = (double*) malloc(sizeof(double)*arraylength);
  double *S_x = (double*) malloc(sizeof(double)*arraylength);
  double *S_y = (double*) malloc(sizeof(double)*arraylength);
  double *S_z = (double*) malloc(sizeof(double)*arraylength);

  /*
     Hybrid_Noble2D failures
     1) rho=0: nans as inputs to Newton-Raphson, causing error code 201
     2) rho=1e15: triggers failure to converge, causing error code 101
  */

  for(int i=0; i<arraylength; i++) {
    lapse[i] = 1.0;
    betax[i] = 0.0;
    betay[i] = 0.0;
    betaz[i] = 0.0;
    gxx[i] = 1.0;
    gxy[i] = 0.0;
    gxz[i] = 0.0;
    gyy[i] = 1.0;
    gyz[i] = 0.0;
    gzz[i] = 1.0;
    Bx[i] = By[i] = Bz[i] = 0.0;
    rho_star[i] = 1e-2;
    tau[i]   = 1e-2;
    S_x[i] = S_y[i] = S_z[i] = 1000.0*tau[i]*(tau[i] + 2.0e-2);
  }
  rho_star[1] = 1e15;
  rho_star[2] = 0.0;

  for(int i=0; i<arraylength; i++) {
    ghl_con2prim_diagnostics diagnostics;
    ghl_initialize_diagnostics(&diagnostics);
    ghl_metric_quantities ADM_metric;
    ghl_primitive_quantities prims;
    ghl_conservative_quantities cons, cons_undens;

    ghl_initialize_metric(lapse[i],
                      betax[i], betay[i], betaz[i],
                      gxx[i], gxy[i], gxz[i],
                      gyy[i], gyz[i], gzz[i],
                      &ADM_metric);

    ghl_ADM_aux_quantities metric_aux;
    ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

    ghl_initialize_primitives(
                        poison, poison, poison,
                        poison, poison, poison,
                        Bx[i], By[i], Bz[i],
                        poison, poison, poison, &prims);

    if (i==0) {
      params.calc_prim_guess = false;
      prims.rho   = cons.rho/ADM_metric.sqrt_detgamma;
      prims.u0    = 1.0;
      prims.vU[0] = 2.0;
      prims.vU[1] = 2.0;
      prims.vU[2] = 2.0;
      prims.Y_e = cons.Y_e/cons.rho;
      prims.temperature = eos.T_max;
      ghl_hybrid_compute_P_cold(&eos, prims.rho, &prims.press);
    }

    ghl_initialize_conservatives(rho_star[i], tau[i],
             S_x[i], S_y[i], S_z[i],
             poison, poison, &cons);

    ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);
    int check = ghl_con2prim_hybrid_multi_method(&params, &eos, &ADM_metric, &metric_aux, &cons_undens, &prims, &diagnostics);
    if(check != i+1)
      ghl_error("Noble2D has returned a different failure code: old %d and new %d", i+1, check);

    if(i==0) {
      // This just gets coverage for the success branches
      params.main_routine = Font1D;
      int check = ghl_con2prim_hybrid_multi_method(&params, &eos, &ADM_metric, &metric_aux, &cons_undens, &prims, &diagnostics);
      if(check != 0)
        ghl_error("Font1D has returned a different failure code: old %d and new %d", 0, check);
      params.main_routine = Noble2D;
      for (int j=0; j<3; j++) {
        params.backup_routine[2-j] = Font1D;
        int check = ghl_con2prim_hybrid_multi_method(&params, &eos, &ADM_metric, &metric_aux, &cons_undens, &prims, &diagnostics);
        if(check != 0)
          ghl_error("Font1D has returned a different failure code: old %d and new %d", 0, check);
        params.backup_routine[2-j] = Noble2D;
      }
      params.calc_prim_guess = true;
    } else if (i==3) {
      // Here, we can check the Font Fix failure condition (there's just one return value)
      params.backup_routine[0] = Font1D;
      int check = ghl_con2prim_hybrid_multi_method(&params, &eos, &ADM_metric, &metric_aux, &cons_undens, &prims, &diagnostics);
      if(check != 1)
        ghl_error("Font1D has returned a different failure code: old %d and new %d", 1, check);
      params.backup_routine[0] = None;
    }
  }
  return 0;
}
