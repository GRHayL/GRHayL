// Thorn      : GRHayL
// File       : unit_test_data_primitives.c
// Author(s)  : Leo Werneck & Samuel Cupp
// Description: In this file we provide an extensive unit test of
//              the Con2Prim gem.
#include "ghl_unit_tests.h"


int main(int argc, char **argv) {

  // These variables set up the tested range of values and number of sampling points.
  // number of sampling points in density and pressure
  const int npoints = 80;
  const int arraylength = npoints*npoints;

  double test_rho_min = 1e-12; //Minimum input density
  double test_rho_max = 1e-3; //Maximum input density

  double dummy, dummy1, dummy2, dummy3;
  double poison = 1e200;
  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  int backup_routine[3] = {None,None,None};
  bool calc_prims_guess = true;
  double Psi6threshold = 1e100; //Taken from magnetizedTOV.par
  const bool evolve_entropy = false;
  const bool evolve_temperature = false;
  const double Lorenz_damping_factor = 0;
  const double W_max = 10.0;

  int neos = 1;
  double rho_b_min = 1e-12;
  double rho_b_max = 1e300; //IGM default
  double Gamma_th = 2.0; //Taken from magnetizedTOV.par
  double rho_ppoly[1] = {0.0};
  double Gamma_ppoly[1] = {2.0};
  double k_ppoly0 = 1.0;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  ghl_parameters params;
  ghl_initialize_params(
        Noble2D, backup_routine, evolve_entropy, evolve_temperature, calc_prims_guess,
        Psi6threshold, W_max, Lorenz_damping_factor, &params);

  ghl_eos_parameters eos;
  ghl_initialize_hybrid_eos_functions_and_params(
        rho_b_min, rho_b_min, rho_b_max,
        neos, rho_ppoly, Gamma_ppoly,
        k_ppoly0, Gamma_th, &eos);

  // Compute the density step size
  const double lrmin        = log(test_rho_min);
  const double lrmax        = log(test_rho_max);
  const double dlr          = (lrmax - lrmin)/(npoints-1);

  // Allocate memory for the metrics
  double *gxx = (double*) malloc(sizeof(double)*arraylength);
  double *gxy = (double*) malloc(sizeof(double)*arraylength);
  double *gxz = (double*) malloc(sizeof(double)*arraylength);
  double *gyy = (double*) malloc(sizeof(double)*arraylength);
  double *gyz = (double*) malloc(sizeof(double)*arraylength);
  double *gzz = (double*) malloc(sizeof(double)*arraylength);

  double *lapse = (double*) malloc(sizeof(double)*arraylength);
  double *betax = (double*) malloc(sizeof(double)*arraylength);
  double *betay = (double*) malloc(sizeof(double)*arraylength);
  double *betaz = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for the primitives
  double *rho_b = (double*) malloc(sizeof(double)*arraylength);
  double *press = (double*) malloc(sizeof(double)*arraylength);
  double *vx = (double*) malloc(sizeof(double)*arraylength);
  double *vy = (double*) malloc(sizeof(double)*arraylength);
  double *vz = (double*) malloc(sizeof(double)*arraylength);
  double *Bx = (double*) malloc(sizeof(double)*arraylength);
  double *By = (double*) malloc(sizeof(double)*arraylength);
  double *Bz = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for the conservatives
  double *rho_star = (double*) malloc(sizeof(double)*arraylength);
  double *tau = (double*) malloc(sizeof(double)*arraylength);
  double *S_x = (double*) malloc(sizeof(double)*arraylength);
  double *S_y = (double*) malloc(sizeof(double)*arraylength);
  double *S_z = (double*) malloc(sizeof(double)*arraylength);

  srand(0);

  for(int j=0;j<npoints;j++) { // Density loop
    const double xrho = exp(lrmin + dlr*j);
    double P_cold   = 0.0;
    double eps_cold = 0.0;
    ghl_hybrid_compute_P_cold_and_eps_cold(&eos, xrho, &P_cold, &eps_cold);

    // Compute the pressure step size
    const double lpmin = log(P_cold);
    const double lpmax = log(P_cold*1.0e6);
    const double dlp   = (lpmax - lpmin)/(npoints-1);
    for(int i=0;i<npoints;i++) { // Pressure loop
      const int index = i + j*npoints;
      // Start by generating the metric and primitives

      rho_b[index] = xrho;
      press[index] = exp(lpmin + dlp*i);

      ghl_randomize_metric(
            &lapse[index], &betax[index], &betay[index], &betaz[index],
            &gxx[index], &gxy[index], &gxz[index],
            &gyy[index], &gyz[index], &gzz[index]);

      ghl_randomize_primitives(
            &eos, rho_b[index], press[index],
            &vx[index], &vy[index], &vz[index],
            &Bx[index], &By[index], &Bz[index]);

      ghl_con2prim_diagnostics diagnostics;
      ghl_initialize_diagnostics(&diagnostics);
      ghl_metric_quantities ADM_metric;
      ghl_primitive_quantities prims;
      ghl_conservative_quantities cons;
      ghl_stress_energy Tmunu;

      ghl_initialize_metric(
            lapse[index], betax[index], betay[index], betaz[index],
            gxx[index], gxy[index], gxz[index],
            gyy[index], gyz[index], gzz[index],
            &ADM_metric);

      ghl_ADM_aux_quantities metric_aux;
      ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

      ghl_initialize_primitives(
            rho_b[index], press[index], poison,
            vx[index], vy[index], vz[index],
            Bx[index], By[index], Bz[index],
            poison, poison, poison, // entropy, Y_e, temp
            &prims);

      // We need epsilon to compute the enthalpy in ghl_compute_conservs_and_Tmunu;
      // This normally happens in the ghl_enforce_primitive_limits_and_compute_u0 function
      prims.eps = eps_cold + (prims.press-P_cold)/(eos.Gamma_th-1.0)/prims.rho;

      ghl_error_codes_t error = ghl_limit_v_and_compute_u0(
            &params, &ADM_metric, &prims, &diagnostics.speed_limited);
      if(error)
        ghl_read_error_codes(error);

      // Compute conservatives based on these primitives
      ghl_compute_conservs_and_Tmunu(
            &ADM_metric, &metric_aux, &prims, &cons, &Tmunu);

      ghl_return_primitives(
            &prims, &rho_b[index], &press[index], &dummy,
            &vx[index], &vy[index], &vz[index],
            &Bx[index], &By[index], &Bz[index],
            &dummy1, &dummy2, &dummy3);

      ghl_return_conservatives(
            &cons, &rho_star[index], &tau[index],
            &S_x[index], &S_y[index], &S_z[index],
            &dummy1, &dummy2);
    }
  }

  // Apply inequality fix branches
  // Edge case 1) small B
  //   Subcase a) S fix 1, Bx < By (also catches edge case of Bbar < 1e-300)
  //   Subcase b) Don't fix S, Bx > By
  //   Subcase b) Don't fix S, Bx < Bz
  // Edge case 2) \psi^6 too large
  //   Subcase a) \tau fix 1 & 2, don't fix S
  //   Subcase b) Don't fix \tau, S fix 2

  ghl_metric_quantities ADM_metric;
  ghl_initialize_metric(
        1.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 0.0, 1.0,
        &ADM_metric);
  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

  // Ouput data to files and generate perturbed data
  FILE* input = fopen_with_check("ET_Legacy_primitives_input.bin", "wb");

  fwrite(&arraylength, sizeof(int), 1, input);
  fwrite(gxx,   sizeof(double), arraylength, input);
  fwrite(gxy,   sizeof(double), arraylength, input);
  fwrite(gxz,   sizeof(double), arraylength, input);
  fwrite(gyy,   sizeof(double), arraylength, input);
  fwrite(gyz,   sizeof(double), arraylength, input);
  fwrite(gzz,   sizeof(double), arraylength, input);
  fwrite(lapse, sizeof(double), arraylength, input);
  fwrite(betax, sizeof(double), arraylength, input);
  fwrite(betay, sizeof(double), arraylength, input);
  fwrite(betaz, sizeof(double), arraylength, input);

  fwrite(Bx,    sizeof(double), arraylength, input);
  fwrite(By,    sizeof(double), arraylength, input);
  fwrite(Bz,    sizeof(double), arraylength, input);

  fwrite(rho_star, sizeof(double), arraylength, input);
  fwrite(tau,      sizeof(double), arraylength, input);
  fwrite(S_x,      sizeof(double), arraylength, input);
  fwrite(S_y,      sizeof(double), arraylength, input);
  fwrite(S_z,      sizeof(double), arraylength, input);

  fclose(input);

  const double perturb = 1.0e-12;

  for(int index=0; index<arraylength; index++) { // Pressure loop
    //This is meant to simulate some round-off error that deviates from the "true" values
    Bx[index]       = (1.0 + randf(-1,1)*perturb)*Bx[index];
    By[index]       = (1.0 + randf(-1,1)*perturb)*By[index];
    Bz[index]       = (1.0 + randf(-1,1)*perturb)*Bz[index];
    rho_star[index] = (1.0 + randf(-1,1)*perturb)*rho_star[index];
    tau[index]      = (1.0 + randf(-1,1)*perturb)*tau[index];
    S_x[index]      = (1.0 + randf(-1,1)*perturb)*S_x[index];
    S_y[index]      = (1.0 + randf(-1,1)*perturb)*S_y[index];
    S_z[index]      = (1.0 + randf(-1,1)*perturb)*S_z[index];
  }

  input = fopen_with_check("ET_Legacy_primitives_input_pert.bin", "wb");

  fwrite(gxx,   sizeof(double), arraylength, input);
  fwrite(gxy,   sizeof(double), arraylength, input);
  fwrite(gxz,   sizeof(double), arraylength, input);
  fwrite(gyy,   sizeof(double), arraylength, input);
  fwrite(gyz,   sizeof(double), arraylength, input);
  fwrite(gzz,   sizeof(double), arraylength, input);
  fwrite(lapse, sizeof(double), arraylength, input);
  fwrite(betax, sizeof(double), arraylength, input);
  fwrite(betay, sizeof(double), arraylength, input);
  fwrite(betaz, sizeof(double), arraylength, input);

  fwrite(Bx,    sizeof(double), arraylength, input);
  fwrite(By,    sizeof(double), arraylength, input);
  fwrite(Bz,    sizeof(double), arraylength, input);

  fwrite(rho_star, sizeof(double), arraylength, input);
  fwrite(tau,      sizeof(double), arraylength, input);
  fwrite(S_x,      sizeof(double), arraylength, input);
  fwrite(S_y,      sizeof(double), arraylength, input);
  fwrite(S_z,      sizeof(double), arraylength, input);

  fclose(input);
}
