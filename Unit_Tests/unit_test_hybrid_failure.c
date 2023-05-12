#include "unit_tests.h"
int main(int argc, char **argv) {

  const int arraylength = 4;

  const double poison = 0.0/0.0;

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const int backup_routine[3] = {Noble2D,Noble2D,Noble2D};
  const bool evolve_entropy = false;
  const bool evolve_temperature = false;
  const bool calc_prims_guess = true;
  const double Psi6threshold = 1e100;
  const bool Cupp_fix = true;

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
  GRHayL_parameters params;
  grhayl_initialize(Noble2D, backup_routine, evolve_entropy, evolve_temperature, calc_prims_guess,
                    Psi6threshold, Cupp_fix, 0.0, &params);

  eos_parameters eos;
  initialize_hybrid_eos_functions_and_params(W_max,
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
  rho_star[2] = 1e15;
  rho_star[3] = 0.0;

  for(int i=0; i<arraylength; i++) {
    // Setting backups to Noble2D allows us to traverse the multi-method
    // while doing the failure checking
    if(i>0 && i<4)
      params.backup_routine[3-i] = None;

    con2prim_diagnostics diagnostics;
    initialize_diagnostics(&diagnostics);
    metric_quantities metric;
    primitive_quantities prims;
    conservative_quantities cons, cons_undens;

    initialize_metric(lapse[i],
                      gxx[i], gxy[i], gxz[i],
                      gyy[i], gyz[i], gzz[i],
                      betax[i], betay[i], betaz[i],
                      &metric);

    initialize_primitives(
                        poison, poison, poison,
                        poison, poison, poison,
                        Bx[i], By[i], Bz[i],
                        poison, poison, poison, &prims);

    if (i==0 || i==1) {
      params.calc_prim_guess = false;
      prims.rho = cons.rho/metric.psi6;
      prims.vx = 0.9;
      prims.vy = 0.9;
      prims.vz = 0.9;
      prims.Y_e = cons.Y_e/cons.rho;
      prims.temperature = eos.T_max;
      eos.hybrid_compute_P_cold(&eos, prims.rho, &prims.press);
      if(i==1) {
        prims.vx = 0.5;
        prims.vy = 0.5;
        prims.vz = 0.5;
      }
    }

    initialize_conservatives(rho_star[i], tau[i],
             S_x[i], S_y[i], S_z[i],
             poison, poison, &cons);

    undensitize_conservatives(&metric, &cons, &cons_undens);
    int check = grhayl_con2prim_multi_method(&params, &eos, &metric, &cons_undens, &prims, &diagnostics);
    if(check != i+1)
      grhayl_error("Noble2D has returned a different failure code: old %d and new %d", check, i+1);

    // All this nonsense is so we can check out all the different ways to traverse grhayl_con2prim_multi_method
    // while doing the error checking
    if (i==1) {
      params.backup_routine[2] = FontFix;
      int check = grhayl_con2prim_multi_method(&params, &eos, &metric, &cons_undens, &prims, &diagnostics);
      if(check != 0)
        grhayl_error("FontFix has returned a different failure code: old %d and new %d", check, 0);
      params.backup_routine[2] = None;

      params.calc_prim_guess = true;
    } else if (i==2) {
      params.backup_routine[1] = Noble2D;
      params.backup_routine[2] = FontFix;
      int check = grhayl_con2prim_multi_method(&params, &eos, &metric, &cons_undens, &prims, &diagnostics);
      if(check != 0)
        grhayl_error("FontFix has returned a different failure code: old %d and new %d", check, 0);
      params.backup_routine[1] = None;
      params.backup_routine[2] = None;
    } else if (i==3) {
      params.backup_routine[0] = Noble2D;
      params.backup_routine[1] = FontFix;
      int check = grhayl_con2prim_multi_method(&params, &eos, &metric, &cons_undens, &prims, &diagnostics);
      if(check != 1)
        grhayl_error("FontFix has returned a different failure code: old %d and new %d", check, 1);
      params.backup_routine[0] = None;
      params.backup_routine[1] = None;
    }
  }


  return 0;
}
