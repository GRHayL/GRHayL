// Thorn      : GRHayL
// File       : unit_test_data_con2prim.c
// Author(s)  : Leo Werneck & Samuel Cupp
// Description: In this file we provide an extensive unit test of
//              the Con2Prim gem.
#include "unit_tests.h"


int main(int argc, char **argv) {

  FILE* input = fopen("ET_Legacy_primitives_input.bin", "rb");
  check_file_was_successfully_open(input, "ET_Legacy_primitives_input.bin");

  int arraylength;
  int key = fread(&arraylength, sizeof(int), 1, input);
  if( key != 1 || arraylength < 1 )
    grhayl_error("An error has occured with reading the grid size. "
                 "Please check that Noble2D_initial_data.bin"
                 "is up-to-date with current test version.\n");

  const double poison = 1e300;
  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const int backup_routine[3] = {None,None,None};
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
  GRHayL_parameters params;
  grhayl_initialize(Noble2D, backup_routine, false /*evolve entropy*/, false /*evolve temperature*/, calc_prims_guess,
                    Psi6threshold, 0 /*Cupp Fix*/, 0 /*Lorenz damping factor*/, &params);

  eos_parameters eos;
  initialize_hybrid_eos_functions_and_params(W_max,
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
    grhayl_error("An error has occured with reading in initial data. "
                 "Please check that comparison data "
                 "is up-to-date with current test version.\n");

  fclose(input);

  // The output for this test is provided by IllinoisGRMHD via the ET
  // for validation with legacy code
  FILE* output = fopen("ET_Legacy_primitives_output.bin", "rb");
  check_file_was_successfully_open(output, "ET_Legacy_primitives_output.bin");

  key  = fread(rho_b_trusted, sizeof(double), arraylength, output);
  key += fread(press_trusted, sizeof(double), arraylength, output);
  key += fread(vx_trusted, sizeof(double), arraylength, output);
  key += fread(vy_trusted, sizeof(double), arraylength, output);
  key += fread(vz_trusted, sizeof(double), arraylength, output);

  if( key != 5*arraylength)
    grhayl_error("An error has occured with reading in initial data. "
                 "Please check that comparison data "
                 "is up-to-date with current test version.\n");

  fclose(output);

  output = fopen("ET_Legacy_primitives_output_pert.bin", "rb");
  check_file_was_successfully_open(output, "ET_Legacy_primitives_output_pert.bin");

  key  = fread(rho_b_pert, sizeof(double), arraylength, output);
  key += fread(press_pert, sizeof(double), arraylength, output);
  key += fread(vx_pert, sizeof(double), arraylength, output);
  key += fread(vy_pert, sizeof(double), arraylength, output);
  key += fread(vz_pert, sizeof(double), arraylength, output);

  if( key != 5*arraylength)
    grhayl_error("An error has occured with reading in initial data. "
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
    initialize_diagnostics(&diagnostics);
    metric_quantities metric;
    primitive_quantities prims;
    conservative_quantities cons, cons_undens;

    initialize_metric(lapse[index],
                      gxx[index], gxy[index], gxz[index],
                      gyy[index], gyz[index], gzz[index],
                      betax[index], betay[index], betaz[index],
                      &metric);

    initialize_primitives(rho_b[index], press[index], eps[index],
                          vx[index], vy[index], vz[index],
                          Bx[index], By[index], Bz[index],
                          poison, poison, poison, // entropy, Y_e, temp
                          &prims);

    initialize_conservatives(rho_star[index], tau[index],
                             S_x[index], S_y[index], S_z[index],
                             poison, poison, &cons);

    const primitive_quantities prims_orig = prims;
    int check = 0;
    if(cons.rho > 0.0) {

      //This applies the inequality (or "Faber") fixes on the conservatives
      if(eos.eos_type == 0) { //Hybrid-only
        if(index == arraylength-2 || index == arraylength-1) {
          params.psi6threshold = 1e-1; // Artificially triggering fix
          apply_inequality_fixes(&params, &eos, &metric, &prims, &cons, &diagnostics);
          params.psi6threshold = Psi6threshold;
        } else {
          apply_inequality_fixes(&params, &eos, &metric, &prims, &cons, &diagnostics);
        }
      }

      // The Con2Prim routines require the undensitized variables, but IGM evolves the densitized variables.
      undensitize_conservatives(&metric, &cons, &cons_undens);

      /************* Conservative-to-primitive recovery ************/
      check = grhayl_con2prim_multi_method(&params, &eos, &metric, &cons_undens, &prims, &diagnostics);

      if(check!=0)
        check = Hybrid_Font_Fix(&params, &eos, &metric, &cons, &prims, &diagnostics);
      /*************************************************************/

      /********** Artificial Font fix for code comparison **********
      This point corresponds to the second-to-last element of the edge
      cases for apply_inequality_fixes function. Due to improvements
      in the Noble2D routine, the GRHayL code doesn't trigger font
      fix while the old code does. This is also true for several of
      the more 'physically motivated' indices, but the edge case data
      isn't constructed to be physically reasonable, causing the font
      fix to significantly affect the results. For the normal data, the
      font fix data from IllinoisGRMHD and Noble2D data from GRHayL
      agree within tolerance, further validating the changes to Noble2D.
      Since this data doesn't match due to the code changing for the
      better, we just have to manually trigger font fix to reproduce
      the behavior of IllinoisGRMHD.
      **************************************************************/
      if(index==arraylength-2)
        check = Hybrid_Font_Fix(&params, &eos, &metric, &cons, &prims, &diagnostics);

      if(check)
        printf("Con2Prim and Font fix failed!");
    } else {
      diagnostics.failure_checker+=1;
      reset_prims_to_atmosphere(&eos, &prims);
      //TODO: Validate reset? (rhob press v)
      printf("Negative rho_* triggering atmospheric reset.\n");
    } // if rho_star > 0

    //Now we compute the difference between original & new primitives and conservatives, for diagnostic purposes:
    prims_abs_error[0] += fabs(prims.rho - prims_orig.rho);
    prims_abs_error[1] += fabs(prims.press - prims_orig.press);
    prims_abs_error[2] += fabs(prims.vx - prims_orig.vx);
    prims_abs_error[3] += fabs(prims.vy - prims_orig.vy);
    prims_abs_error[4] += fabs(prims.vz - prims_orig.vz);

    prims_rel_error[0] += relative_error(prims_orig.rho, prims.rho);
    prims_rel_error[1] += relative_error(prims_orig.press, prims.press);
    prims_rel_error[2] += relative_error(prims_orig.vx, prims.vx);
    prims_rel_error[3] += relative_error(prims_orig.vy, prims.vy);
    prims_rel_error[4] += relative_error(prims_orig.vz, prims.vz);

    // Now, we load the trusted/perturbed data for this index and validate the computed results.
    primitive_quantities prims_trusted, prims_pert;

    initialize_primitives(rho_b_trusted[index], press_trusted[index], prims.eps, // Old code has no eps variable
                          vx_trusted[index], vy_trusted[index], vz_trusted[index],
                          poison, poison, poison, // B is C2P input, not output
                          poison, poison, poison, // entropy, Y_e, temp
                          &prims_trusted);

    initialize_primitives(rho_b_pert[index], press_pert[index], prims.eps, // Old code has no eps variable
                          vx_pert[index], vy_pert[index], vz_pert[index],
                          poison, poison, poison, // B is C2P input, not output
                          poison, poison, poison, // entropy, Y_e, temp
                          &prims_pert);

    prims_trusted_abs_error[0] += fabs(prims_trusted.rho - prims_orig.rho);
    prims_trusted_abs_error[1] += fabs(prims_trusted.press - prims_orig.press);
    prims_trusted_abs_error[2] += fabs(prims_trusted.vx - prims_orig.vx);
    prims_trusted_abs_error[3] += fabs(prims_trusted.vy - prims_orig.vy);
    prims_trusted_abs_error[4] += fabs(prims_trusted.vz - prims_orig.vz);

    prims_trusted_rel_error[0] += relative_error(prims_orig.rho, prims_trusted.rho);
    prims_trusted_rel_error[1] += relative_error(prims_orig.press, prims_trusted.press);
    prims_trusted_rel_error[2] += relative_error(prims_orig.vx, prims_trusted.vx);
    prims_trusted_rel_error[3] += relative_error(prims_orig.vy, prims_trusted.vy);
    prims_trusted_rel_error[4] += relative_error(prims_orig.vz, prims_trusted.vz);

    prims_pert_abs_error[0] += fabs(prims_pert.rho - prims_orig.rho);
    prims_pert_abs_error[1] += fabs(prims_pert.press - prims_orig.press);
    prims_pert_abs_error[2] += fabs(prims_pert.vx - prims_orig.vx);
    prims_pert_abs_error[3] += fabs(prims_pert.vy - prims_orig.vy);
    prims_pert_abs_error[4] += fabs(prims_pert.vz - prims_orig.vz);

    prims_pert_rel_error[0] += relative_error(prims_orig.rho, prims_pert.rho);
    prims_pert_rel_error[1] += relative_error(prims_orig.press, prims_pert.press);
    prims_pert_rel_error[2] += relative_error(prims_orig.vx, prims_pert.vx);
    prims_pert_rel_error[3] += relative_error(prims_orig.vy, prims_pert.vy);
    prims_pert_rel_error[4] += relative_error(prims_orig.vz, prims_pert.vz);

    validate_primitives(params.evolve_entropy, &eos, &prims_trusted, &prims, &prims_pert);
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
  grhayl_info("ET_Legacy conservatives-to-primitives test has passed!\n");
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
