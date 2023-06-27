// Thorn      : GRHayL
// File       : unit_test_data_con2prim.c
// Author(s)  : Leo Werneck & Samuel Cupp
// Description: In this file we provide an extensive unit test of
//              the Con2Prim gem.
#include "unit_tests.h"


int main(int argc, char **argv) {

  FILE* input = fopen_with_check("ET_Legacy_primitives_input.bin", "rb");

  int arraylength;
  int key = fread(&arraylength, sizeof(int), 1, input);
  if( key != 1 || arraylength < 1 )
    ghl_error("An error has occured with reading the grid size. "
                 "Please check that Noble2D_initial_data.bin"
                 "is up-to-date with current test version.\n");

  const double poison = 1e300;
  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const int backup_routine[3] = {Font1D,None,None};
  const bool calc_prims_guess = true;
  const double Psi6threshold = 1e100; //Taken from magnetizedTOV.par

  const int neos = 1;
  const double W_max = 10.0; //IGM default
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300; //IGM default
  const double Gamma_th = 2.0; //Taken from magnetizedTOV.par
  const double rho_ppoly[1] = {0.0};
  const double Gamma_ppoly[1] = {2.0};
  const double k_ppoly0 = 1.0;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  ghl_parameters params;
  ghl_initialize_params(Noble2D, backup_routine, false /*evolve entropy*/, false /*evolve temperature*/, calc_prims_guess,
                    Psi6threshold, 0 /*Cupp Fix*/, 0 /*Lorenz damping factor*/, &params);

  eos_parameters eos;
  ghl_initialize_hybrid_eos_functions_and_params(W_max,
                                             rho_b_min, rho_b_min, rho_b_max,
                                             neos, rho_ppoly, Gamma_ppoly,
                                             k_ppoly0, Gamma_th, &eos);

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
  double *eps = (double*) malloc(sizeof(double)*arraylength);
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

  // Allocate memory for for trusted output
  double *rho_b_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *press_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *vx_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *vy_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *vz_trusted = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for for perturbed output
  double *rho_b_pert = (double*) malloc(sizeof(double)*arraylength);
  double *press_pert = (double*) malloc(sizeof(double)*arraylength);
  double *vx_pert = (double*) malloc(sizeof(double)*arraylength);
  double *vy_pert = (double*) malloc(sizeof(double)*arraylength);
  double *vz_pert = (double*) malloc(sizeof(double)*arraylength);

  // Read in data from file to ensure portability
  key  = fread(gxx,   sizeof(double), arraylength, input);
  key += fread(gxy,   sizeof(double), arraylength, input);
  key += fread(gxz,   sizeof(double), arraylength, input);
  key += fread(gyy,   sizeof(double), arraylength, input);
  key += fread(gyz,   sizeof(double), arraylength, input);
  key += fread(gzz,   sizeof(double), arraylength, input);
  key += fread(lapse, sizeof(double), arraylength, input);
  key += fread(betax, sizeof(double), arraylength, input);
  key += fread(betay, sizeof(double), arraylength, input);
  key += fread(betaz, sizeof(double), arraylength, input);

  key += fread(Bx, sizeof(double), arraylength, input);
  key += fread(By, sizeof(double), arraylength, input);
  key += fread(Bz, sizeof(double), arraylength, input);

  key += fread(rho_star, sizeof(double), arraylength, input);
  key += fread(tau, sizeof(double), arraylength, input);
  key += fread(S_x, sizeof(double), arraylength, input);
  key += fread(S_y, sizeof(double), arraylength, input);
  key += fread(S_z, sizeof(double), arraylength, input);

  if( key != (10+3+5)*arraylength)
    ghl_error("An error has occured with reading in initial data. "
                 "Please check that comparison data "
                 "is up-to-date with current test version.\n");

  fclose(input);

  // The output for this test is provided by IllinoisGRMHD via the ET
  // for validation with legacy code
  FILE* output = fopen_with_check("ET_Legacy_primitives_output.bin", "rb");

  key  = fread(rho_b_trusted, sizeof(double), arraylength, output);
  key += fread(press_trusted, sizeof(double), arraylength, output);
  key += fread(vx_trusted, sizeof(double), arraylength, output);
  key += fread(vy_trusted, sizeof(double), arraylength, output);
  key += fread(vz_trusted, sizeof(double), arraylength, output);

  if( key != 5*arraylength)
    ghl_error("An error has occured with reading in initial data. "
                 "Please check that comparison data "
                 "is up-to-date with current test version.\n");

  fclose(output);

  output = fopen_with_check("ET_Legacy_primitives_output_pert.bin", "rb");

  key  = fread(rho_b_pert, sizeof(double), arraylength, output);
  key += fread(press_pert, sizeof(double), arraylength, output);
  key += fread(vx_pert, sizeof(double), arraylength, output);
  key += fread(vy_pert, sizeof(double), arraylength, output);
  key += fread(vz_pert, sizeof(double), arraylength, output);

  if( key != 5*arraylength)
    ghl_error("An error has occured with reading in initial data. "
                 "Please check that comparison data "
                 "is up-to-date with current test version.\n");

  fclose(output);

  double prims_abs_error[5] = {0,0,0,0,0};
  double prims_trusted_abs_error[5] = {0,0,0,0,0};
  double prims_pert_abs_error[5] = {0,0,0,0,0};

  double prims_rel_error[5] = {0,0,0,0,0};
  double prims_trusted_rel_error[5] = {0,0,0,0,0};
  double prims_pert_rel_error[5] = {0,0,0,0,0};

  //Parallelizing this also needs parallel sum of abs/rel error arrays
  for(int index=0; index<arraylength; index++) {
    // Define the various GRHayL structs for the unit tests
    con2prim_diagnostics diagnostics;
    ghl_initialize_diagnostics(&diagnostics);
    metric_quantities ADM_metric;
    primitive_quantities prims;
    conservative_quantities cons, cons_undens;

    ghl_initialize_metric(lapse[index],
                      betax[index], betay[index], betaz[index],
                      gxx[index], gxy[index], gxz[index],
                      gyy[index], gyz[index], gzz[index],
                      &ADM_metric);

    ADM_aux_quantities metric_aux;
    ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

    ghl_initialize_primitives(rho_b[index], press[index], eps[index],
                          vx[index], vy[index], vz[index],
                          Bx[index], By[index], Bz[index],
                          poison, poison, poison, // entropy, Y_e, temp
                          &prims);

    ghl_initialize_conservatives(rho_star[index], tau[index],
                             S_x[index], S_y[index], S_z[index],
                             poison, poison, &cons);

    const primitive_quantities prims_orig = prims;
    int check = 0;
    if(cons.rho > 0.0) {

      //This applies the inequality (or "Faber") fixes on the conservatives
      if(eos.eos_type == 0) { //Hybrid-only
        if(index == arraylength-2 || index == arraylength-1) {
          params.psi6threshold = 1e-1; // Artificially triggering fix
          ghl_apply_conservative_limits(&params, &eos, &ADM_metric, &prims, &cons, &diagnostics);
          params.psi6threshold = Psi6threshold;
        } else {
          ghl_apply_conservative_limits(&params, &eos, &ADM_metric, &prims, &cons, &diagnostics);
        }
      }

      // The Con2Prim routines require the undensitized variables, but IGM evolves the densitized variables.
      ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);

      /************* Conservative-to-primitive recovery ************/
      check = ghl_con2prim_multi_method(&params, &eos, &ADM_metric, &metric_aux, &cons_undens, &prims, &diagnostics);

      /********** Artificial Font fix for code comparison **********
      This point corresponds to the second-to-last element of the edge
      cases for ghl_apply_conservative_limits function. Due to improvements
      in the Noble2D routine, the GRHayL code doesn't trigger Font
      fix while the old code does. This is also true for several of
      the more 'physically motivated' indices, but the edge case data
      isn't constructed to be physically reasonable, causing the Font
      fix to significantly affect the results. For the normal data, the
      Font fix data from IllinoisGRMHD and Noble2D data from GRHayL
      agree within tolerance, further validating the changes to Noble2D.
      Since this data doesn't match due to the code changing for the
      better, we just have to manually trigger Font fix to reproduce
      the behavior of IllinoisGRMHD.
      **************************************************************/
      if(index==arraylength-2)
        check = ghl_hybrid_Font_fix(&params, &eos, &ADM_metric, &metric_aux, &cons, &prims, &diagnostics);
      /*************************************************************/

      if(check)
        printf("Con2Prim failed!");
    } else {
      diagnostics.failure_checker+=1;
      ghl_set_prims_to_constant_atm(&eos, &prims);
      //TODO: Validate reset? (rhob press v)
      printf("Negative rho_* triggering atmospheric reset.\n");
    } // if rho_star > 0

    //Now we compute the difference between original & new primitives and conservatives, for diagnostic purposes:
    prims_abs_error[0] += fabs(prims.rho - prims_orig.rho);
    prims_abs_error[1] += fabs(prims.press - prims_orig.press);
    prims_abs_error[2] += fabs(prims.vU[0] - prims_orig.vU[0]);
    prims_abs_error[3] += fabs(prims.vU[1] - prims_orig.vU[1]);
    prims_abs_error[4] += fabs(prims.vU[2] - prims_orig.vU[2]);

    prims_rel_error[0] += relative_error(prims_orig.rho, prims.rho);
    prims_rel_error[1] += relative_error(prims_orig.press, prims.press);
    prims_rel_error[2] += relative_error(prims_orig.vU[0], prims.vU[0]);
    prims_rel_error[3] += relative_error(prims_orig.vU[1], prims.vU[1]);
    prims_rel_error[4] += relative_error(prims_orig.vU[2], prims.vU[2]);

    // Now, we load the trusted/perturbed data for this index and validate the computed results.
    primitive_quantities prims_trusted, prims_pert;

    ghl_initialize_primitives(rho_b_trusted[index], press_trusted[index], prims.eps, // Old code has no eps variable
                          vx_trusted[index], vy_trusted[index], vz_trusted[index],
                          poison, poison, poison, // B is C2P input, not output
                          poison, poison, poison, // entropy, Y_e, temp
                          &prims_trusted);

    ghl_initialize_primitives(rho_b_pert[index], press_pert[index], prims.eps, // Old code has no eps variable
                          vx_pert[index], vy_pert[index], vz_pert[index],
                          poison, poison, poison, // B is C2P input, not output
                          poison, poison, poison, // entropy, Y_e, temp
                          &prims_pert);

    prims_trusted_abs_error[0] += fabs(prims_trusted.rho - prims_orig.rho);
    prims_trusted_abs_error[1] += fabs(prims_trusted.press - prims_orig.press);
    prims_trusted_abs_error[2] += fabs(prims_trusted.vU[0] - prims_orig.vU[0]);
    prims_trusted_abs_error[3] += fabs(prims_trusted.vU[1] - prims_orig.vU[1]);
    prims_trusted_abs_error[4] += fabs(prims_trusted.vU[2] - prims_orig.vU[2]);

    prims_trusted_rel_error[0] += relative_error(prims_orig.rho, prims_trusted.rho);
    prims_trusted_rel_error[1] += relative_error(prims_orig.press, prims_trusted.press);
    prims_trusted_rel_error[2] += relative_error(prims_orig.vU[0], prims_trusted.vU[0]);
    prims_trusted_rel_error[3] += relative_error(prims_orig.vU[1], prims_trusted.vU[1]);
    prims_trusted_rel_error[4] += relative_error(prims_orig.vU[2], prims_trusted.vU[2]);

    prims_pert_abs_error[0] += fabs(prims_pert.rho - prims_orig.rho);
    prims_pert_abs_error[1] += fabs(prims_pert.press - prims_orig.press);
    prims_pert_abs_error[2] += fabs(prims_pert.vU[0] - prims_orig.vU[0]);
    prims_pert_abs_error[3] += fabs(prims_pert.vU[1] - prims_orig.vU[1]);
    prims_pert_abs_error[4] += fabs(prims_pert.vU[2] - prims_orig.vU[2]);

    prims_pert_rel_error[0] += relative_error(prims_orig.rho, prims_pert.rho);
    prims_pert_rel_error[1] += relative_error(prims_orig.press, prims_pert.press);
    prims_pert_rel_error[2] += relative_error(prims_orig.vU[0], prims_pert.vU[0]);
    prims_pert_rel_error[3] += relative_error(prims_orig.vU[1], prims_pert.vU[1]);
    prims_pert_rel_error[4] += relative_error(prims_orig.vU[2], prims_pert.vU[2]);

    if( validate(prims_trusted.rho, prims.rho, prims_pert.rho) )
      ghl_error("Test unit_test_ghl_hybrid_Noble2D has failed for variable rho.\n"
                   "  rho trusted %.14e computed %.14e perturbed %.14e\n"
                   "  rel.err. %.14e %.14e\n", prims_trusted.rho, prims.rho, prims_pert.rho,
                                               relative_error(prims_trusted.rho, prims.rho),
                                               relative_error(prims_trusted.rho, prims_pert.rho));

    const double min_rel = 8.0e-14; // This is the default relative tolerance cutoff used by validate()
    const double pressure_cutoff = 1.0e-18;
    // Pressure has an additional absolute difference check because the pressure can become very small depending on the
    // input values. The pressure coming out of HARM doesn't have the accuracy to preserve the stringent accuracy requirements
    // demanded elsewhere, so this relaxes the demands on the pressure for very small values.
    if( validate_with_tolerance(prims_trusted.press, prims.press, prims_pert.press, min_rel, pressure_cutoff))
      ghl_error("Test unit_test_ghl_hybrid_Noble2D has failed for variable press.\n"
                   "  press trusted %.14e computed %.14e perturbed %.14e\n"
                   "  rel.err. %.14e %.14e\n", prims_trusted.press, prims.press, prims_pert.press,
                                               relative_error(prims_trusted.press, prims.press),
                                               relative_error(prims_trusted.press, prims_pert.press));

    // Epsilon has a similar issue with pressure, so we compute a cutoff that is consistent with the above choice.
    //const double eps_cutoff = pressure_cutoff/(pow(pressure_cutoff/eos.K_ppoly[0], 1.0/eos.Gamma_ppoly[0]) * (eos.Gamma_ppoly[0] - 1.0));
    const double eps_cutoff = 1.0e-11; // Above computed 1e-9, which seemed too large to make sense as a cutoff
    if( validate_with_tolerance(prims_trusted.eps, prims.eps, prims_pert.eps, min_rel, eps_cutoff))
      ghl_error("Test unit_test_ghl_hybrid_Noble2D has failed for variable eps.\n"
                   "  eps trusted %.14e computed %.14e perturbed %.14e\n"
                   "  rel.err. %.14e %.14e\n", prims_trusted.eps, prims.eps, prims_pert.eps,
                                               relative_error(prims_trusted.eps, prims.eps),
                                               relative_error(prims_trusted.eps, prims_pert.eps));

    if( validate(prims_trusted.vU[0], prims.vU[0], prims_pert.vU[0]) )
      ghl_error("Test unit_test_ghl_hybrid_Noble2D has failed for variable vx.\n"
                   "  vx trusted %.14e computed %.14e perturbed %.14e\n"
                   "  rel.err. %.14e %.14e\n", prims_trusted.vU[0], prims.vU[0], prims_pert.vU[0],
                                               relative_error(prims_trusted.vU[0], prims.vU[0]),
                                               relative_error(prims_trusted.vU[0], prims_pert.vU[0]));

    if(validate(prims_trusted.vU[1], prims.vU[1], prims_pert.vU[1]))
      ghl_error("Test unit_test_ghl_hybrid_Noble2D has failed for variable vy.\n"
                   "  vy trusted %.14e computed %.14e perturbed %.14e\n"
                   "  rel.err. %.14e %.14e\n", prims_trusted.vU[1], prims.vU[1], prims_pert.vU[1],
                                               relative_error(prims_trusted.vU[1], prims.vU[1]),
                                               relative_error(prims_trusted.vU[1], prims_pert.vU[1]));

    if( validate(prims_trusted.vU[2], prims.vU[2], prims_pert.vU[2]) )
      ghl_error("Test unit_test_ghl_hybrid_Noble2D has failed for variable vz.\n"
                   "  vz trusted %.14e computed %.14e perturbed %.14e\n"
                   "  rel.err. %.14e %.14e\n", prims_trusted.vU[2], prims.vU[2], prims_pert.vU[2],
                                               relative_error(prims_trusted.vU[2], prims.vU[2]),
                                               relative_error(prims_trusted.vU[2], prims_pert.vU[2]));
  }

  output = fopen("prims_summary.asc","w");

  fprintf(output, "Absolute error\n"
          "rho_b %e %e %e\n"
          "press %e %e %e\n"
          "   vx %e %e %e\n"
          "   vy %e %e %e\n"
          "   vz %e %e %e\n\n"
          "Relative error\n"
          "rho_b %e %e %e\n"
          "press %e %e %e\n"
          "   vx %e %e %e\n"
          "   vy %e %e %e\n"
          "   vz %e %e %e\n",
          prims_abs_error[0], prims_trusted_abs_error[0], prims_pert_abs_error[0],
          prims_abs_error[1], prims_trusted_abs_error[1], prims_pert_abs_error[1],
          prims_abs_error[2], prims_trusted_abs_error[2], prims_pert_abs_error[2],
          prims_abs_error[3], prims_trusted_abs_error[3], prims_pert_abs_error[3],
          prims_abs_error[4], prims_trusted_abs_error[4], prims_pert_abs_error[4],
          prims_rel_error[0], prims_trusted_rel_error[0], prims_pert_rel_error[0],
          prims_rel_error[1], prims_trusted_rel_error[1], prims_pert_rel_error[1],
          prims_rel_error[2], prims_trusted_rel_error[2], prims_pert_rel_error[2],
          prims_rel_error[3], prims_trusted_rel_error[3], prims_pert_rel_error[3],
          prims_rel_error[4], prims_trusted_rel_error[4], prims_pert_rel_error[4]);

  fclose(output);
  ghl_info("ET_Legacy conservatives-to-primitives test has passed!\n");
  free(gxx); free(gxy); free(gxz);
  free(gyy); free(gyz); free(gzz);
  free(lapse);
  free(betax); free(betay); free(betaz);
  free(rho_b); free(press); free(eps);
  free(vx); free(vy); free(vz);
  free(Bx); free(By); free(Bz);
  free(rho_star); free(tau);
  free(S_x); free(S_y); free(S_z);
  free(rho_b_trusted); free(press_trusted);
  free(vx_trusted); free(vy_trusted); free(vz_trusted);
  free(rho_b_pert); free(press_pert);
  free(vx_pert); free(vy_pert); free(vz_pert);
}
