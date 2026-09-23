#include "ghl_con2prim.h"
#include "ghl_eos_functions_declaration.h"
#include "ghl_nrpyeos_hybrid.h"
#include "ghl_nrpyeos_tabulated.h"
#include <float.h>

#define init_common_eos_quantities \
  eos->rho_atm = rho_atm;          \
  eos->rho_min = rho_min;          \
  eos->rho_max = rho_max;

ghl_error_codes_t (*ghl_con2prim_multi_method)(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

/*
 * Function    : ghl_initialize_eos_functions()
 * Description : Initializes function pointers in EOS struct to NRPyEOS
 *
 * Input/Output: eos - ghl_eos_parameters struct with the function pointers
 *                     initialized
 */
void ghl_initialize_eos_functions(
      const ghl_eos_t eos_type) {

  // Step 1: Hybrid EOS functions (always available)
  NRPyEOS_initialize_hybrid_functions();

  // Step 2: Tabulated EOS functions require HDF5-backed NRPyEOS tables.
#ifndef GHL_DISABLE_HDF5
  NRPyEOS_initialize_tabulated_functions();
#endif

  // Step 3: General functions (same interface for all EOSs)
  if(eos_type == ghl_eos_hybrid || eos_type == ghl_eos_simple) {
    ghl_con2prim_multi_method = ghl_con2prim_hybrid_multi_method;
    ghl_compute_h_and_cs2 = NRPyEOS_hybrid_compute_enthalpy_and_cs2;
#ifndef GHL_DISABLE_HDF5
  } else if(eos_type == ghl_eos_tabulated) {
    ghl_con2prim_multi_method = ghl_con2prim_tabulated_multi_method;
    ghl_compute_h_and_cs2 = NRPyEOS_tabulated_compute_enthalpy_and_cs2;
#else
  } else if(eos_type == ghl_eos_tabulated) {
    GHL_HDF5_ERROR_IF_USED;
#endif
  }
}

/*
 * Function    : ghl_initialize_simple_eos()
 * Description : Initializes EOS struct elements for a simple EOS
*/
ghl_error_codes_t ghl_initialize_simple_eos(
      const double rho_atm,
      double rho_min,
      double rho_max,
      const double press_atm,
      double press_min,
      double press_max,
      const double Gamma,
      ghl_eos_parameters *eos) {

  if(eos == NULL) {
    return ghl_error_eos_struct_is_null;
  }
  if(!isfinite(rho_atm) || rho_atm <= 0.0) {
    return ghl_error_invalid_rho_atm;
  }
  if(!isfinite(press_atm) || press_atm < 0.0) {
    return ghl_error_invalid_press_atm;
  }
  if(!isfinite(rho_min) || !isfinite(rho_max) || !isfinite(press_min)
     || !isfinite(press_max) || !isfinite(Gamma) || Gamma == 0.0 || Gamma == 1.0) {
    return ghl_error_invalid_eos_parameters;
  }

  // Step 0: Enforce default values
  if(rho_min < 0) {
    ghl_warn("Minimum density not provided. Disabling density floor (rho_min = 0)\n");
    rho_min = 0.0;
  }
  if(rho_max < 0) {
    ghl_warn("Maximum density not provided. Disabling density ceiling (rho_max = 1e300)\n");
    rho_max = 1e300;
  }
  if(rho_max <= 0.0) {
    return ghl_error_invalid_eos_parameters;
  }
  if(rho_min > rho_max) return ghl_error_rho_min_gt_rho_max;

  if(press_min < 0) {
    ghl_warn("Minimum pressure not provided. Disabling pressure floor (press_min = 0)\n");
    press_min = 0.0;
  }
  if(press_max < 0) {
    ghl_warn("Maximum pressure not provided. Disabling pressure ceiling (press_max = 1e300)\n");
    press_max = 1e300;
  }
  if(press_min > press_max) return ghl_error_press_min_gt_press_max;

  ghl_eos_parameters candidate = { 0 };
  ghl_eos_parameters *const output = eos;
  eos = &candidate;

  // Step 1: Set EOS type to Ideal Fluid
  eos->eos_type = ghl_eos_simple;

  // Step 2: Initialize quantities which are common to all EOSs.
  init_common_eos_quantities;

  // Step 3: Set basic ideal fluid EOS parameters.
  eos->rho_atm = rho_atm;
  eos->rho_min = rho_min;
  eos->rho_max = rho_max;
  eos->press_atm = press_atm;
  eos->press_min = press_min;
  eos->press_max = press_max;
  eos->Gamma_th = eos->Gamma_ppoly[0] = Gamma;
  eos->neos = 1;
  eos->K_ppoly[0] = 1;
  eos->rho_ppoly[0] = 0.0;
  eos->p_ppoly[0] = 0.0;
  eos->eps_integ_const[0] = 0.0;
  // Unused-family atmosphere placeholders.
  eos->Y_e_atm = 0.0;
  eos->T_atm = 0.0;

  const double Gm1 = Gamma - 1.0;
  // For Gamma > 1, eps and entropy decrease with density at fixed pressure,
  // so the maxima pair press_max with rho_min and the minima press_min with
  // rho_max.
  // -------------- Ceilings --------------
  if(eos->rho_min == 0.0 || eos->press_max == 1e300) {
    eos->eps_max = eos->press_max > 0.0 ? DBL_MAX : 0.0;
    eos->entropy_max = eos->press_max > 0.0 ? DBL_MAX : 0.0;
  }
  else {
    eos->eps_max = eos->press_max / (eos->rho_min * Gm1);
    eos->entropy_max
          = ghl_hybrid_compute_entropy_function(eos, eos->rho_min, eos->press_max);
  }

  // --------------- Floors ---------------
  eos->eps_min = eos->press_min / (eos->rho_max * Gm1);
  eos->entropy_min
        = ghl_hybrid_compute_entropy_function(eos, eos->rho_max, eos->press_min);

  // --------- Atmospheric values ---------
  eos->eps_atm = eos->press_atm/(eos->rho_atm*Gm1);
  eos->entropy_atm = ghl_hybrid_compute_entropy_function(eos, eos->rho_atm, eos->press_atm);

  // Compute atmospheric tau
  eos->tau_atm = eos->rho_atm * eos->eps_atm;
  // --------------------------------------

  if(!isfinite(eos->rho_atm) || !isfinite(eos->rho_min) || !isfinite(eos->rho_max)
     || !isfinite(eos->press_atm) || !isfinite(eos->press_min)
     || !isfinite(eos->press_max) || !isfinite(eos->eps_atm) || !isfinite(eos->eps_min)
     || !isfinite(eos->eps_max) || !isfinite(eos->entropy_atm)
     || !isfinite(eos->entropy_min) || !isfinite(eos->entropy_max)
     || !isfinite(eos->tau_atm)) {
    return ghl_error_invalid_eos_parameters;
  }
  *output = candidate;
  return ghl_success;
}

/*
 * Function    : ghl_initialize_hybrid_eos()
 * Description : Initializes EOS struct elements for a hybrid EOS
*/
ghl_error_codes_t ghl_initialize_hybrid_eos(
      const double rho_atm,
      double rho_min,
      double rho_max,
      const int neos,
      const double *restrict rho_ppoly,
      const double *restrict Gamma_ppoly,
      const double K_ppoly0,
      const double Gamma_th,
      ghl_eos_parameters *eos) {

  if(eos == NULL) {
    return ghl_error_eos_struct_is_null;
  }
  if(neos < 1 || neos > MAX_EOS_PARAMS) {
    return ghl_error_invalid_neos;
  }
  if(!isfinite(rho_atm) || rho_atm <= 0.0) {
    return ghl_error_invalid_rho_atm;
  }
  if(!isfinite(rho_min) || !isfinite(rho_max) || !isfinite(K_ppoly0)
     || !isfinite(Gamma_th) || Gamma_th == 1.0 || Gamma_ppoly == NULL
     || (neos > 1 && rho_ppoly == NULL)) {
    return ghl_error_invalid_eos_parameters;
  }
  for(int j = 0; j < neos; j++) {
    if(!isfinite(Gamma_ppoly[j]) || Gamma_ppoly[j] == 0.0 || Gamma_ppoly[j] == 1.0) {
      return ghl_error_invalid_eos_parameters;
    }
  }
  for(int j = 0; j < neos - 1; j++) {
    if(!isfinite(rho_ppoly[j]) || rho_ppoly[j] <= 0.0
       || (j > 0 && rho_ppoly[j] <= rho_ppoly[j - 1])) {
      return ghl_error_invalid_eos_parameters;
    }
  }

  // Step 0: Enforce default values
  if(rho_min < 0) {
    ghl_warn("Minimum density not provided. Disabling density floor (rho_min = 0)\n");
    rho_min = 0.0;
  }
  if(rho_max < 0) {
    ghl_warn("Maximum density not provided. Disabling density ceiling (rho_max = 1e300)\n");
    rho_max = 1e300;
  }
  if(rho_max <= 0.0) {
    return ghl_error_invalid_eos_parameters;
  }
  if(rho_min > rho_max) return ghl_error_rho_min_gt_rho_max;

  ghl_eos_parameters candidate = { 0 };
  ghl_eos_parameters *const output = eos;
  eos = &candidate;

  // Step 1: Set EOS type to Hybrid
  eos->eos_type = ghl_eos_hybrid;

  // Step 2: Initialize quantities which are common to all EOSs.
  init_common_eos_quantities;

  // Step 3: Set basic Hybrid EOS parameters.
  eos->neos = neos;
  eos->Gamma_th = Gamma_th;
  eos->K_ppoly[0] = K_ppoly0;
  for(int j = 0; j < neos - 1; j++) {
    eos->rho_ppoly[j] = rho_ppoly[j];
  }
  for(int j = 0; j < neos; j++) {
    eos->Gamma_ppoly[j] = Gamma_ppoly[j];
  }

  // Step 4: Initialize {K_{j}}, j>=1, and {eps_integ_const_{j}}
  ghl_hybrid_set_K_ppoly_and_eps_integ_consts(eos);
  for(int j = 0; j < neos; j++) {
    if(!isfinite(eos->K_ppoly[j]) || (K_ppoly0 != 0.0 && eos->K_ppoly[j] == 0.0)
       || !isfinite(eos->eps_integ_const[j])) {
      return ghl_error_invalid_eos_parameters;
    }
  }

  // Initialize pressure breakpoints after the piece coefficients.
  for(int j = 0; j < eos->neos - 1; j++) {
    double P, eps;
    ghl_hybrid_compute_P_cold_and_eps_cold(eos, eos->rho_ppoly[j], &P, &eps);
    eos->p_ppoly[j] = P;
    if(!isfinite(eos->p_ppoly[j])) {
      return ghl_error_invalid_eos_parameters;
    }
  }

  // -------------- Ceilings --------------
  if(eos->rho_max == 1e300) {
    eos->press_max = DBL_MAX;
    eos->eps_max = DBL_MAX;
    eos->entropy_max = DBL_MAX;
  }
  else {
    ghl_hybrid_compute_P_cold_and_eps_cold(
          eos, eos->rho_max, &eos->press_max, &eos->eps_max);
    eos->entropy_max
          = ghl_hybrid_compute_entropy_function(eos, eos->rho_max, eos->press_max);
  }
  // --------------------------------------

  // --------------- Floors ---------------
  if(eos->rho_min == 0.0) {
    eos->press_min = -DBL_MAX;
    eos->eps_min = -DBL_MAX;
    eos->entropy_min = -DBL_MAX;
  }
  else {
    ghl_hybrid_compute_P_cold_and_eps_cold(
          eos, eos->rho_min, &eos->press_min, &eos->eps_min);
    eos->entropy_min
          = ghl_hybrid_compute_entropy_function(eos, eos->rho_min, eos->press_min);
  }
  // --------------------------------------

  // --------- Atmospheric values ---------
  // Compute atmospheric P and eps
  ghl_hybrid_compute_P_cold_and_eps_cold(eos, eos->rho_atm, &eos->press_atm, &eos->eps_atm);

  // Compute atmospheric entropy
  eos->entropy_atm = ghl_hybrid_compute_entropy_function(eos, eos->rho_atm, eos->press_atm);

  // Compute atmospheric tau
  eos->tau_atm = eos->rho_atm * eos->eps_atm;
  // --------------------------------------
  // Unused-family atmosphere placeholders.
  eos->Y_e_atm = 0.0;
  eos->T_atm = 0.0;

  if(!isfinite(eos->rho_atm) || !isfinite(eos->rho_min) || !isfinite(eos->rho_max)
     || !isfinite(eos->press_atm) || !isfinite(eos->press_min)
     || !isfinite(eos->press_max) || !isfinite(eos->eps_atm) || !isfinite(eos->eps_min)
     || !isfinite(eos->eps_max) || !isfinite(eos->entropy_atm)
     || !isfinite(eos->entropy_min) || !isfinite(eos->entropy_max)
     || !isfinite(eos->tau_atm)) {
    return ghl_error_invalid_eos_parameters;
  }
  *output = candidate;
  return ghl_success;
}

/*
 * Function    : ghl_initialize_tabulated_eos()
 * Description : Initializes EOS struct elements for tabulated EOS
*/
ghl_error_codes_t ghl_initialize_tabulated_eos(
      const char *table_filepath,
      const ghl_eos_table_t table_type,
      const bool clean_sound_speed,
      const bool enable_neural_net_c2p,
      const double rho_atm,
      double rho_min,
      double rho_max,
      const double Y_e_atm,
      double Y_e_min,
      double Y_e_max,
      const double T_atm,
      double T_min,
      double T_max,
      ghl_eos_parameters *eos) {

  if(eos == NULL) {
    return ghl_error_eos_struct_is_null;
  }
#ifdef GHL_DISABLE_HDF5
  *eos = (ghl_eos_parameters){ .eos_type = ghl_eos_tabulated };
  return ghl_error_used_disabled_hdf5;
#else
  ghl_eos_parameters candidate = { 0 };
  ghl_eos_parameters *const output = eos;
  eos = &candidate;

  // Step 1: Set EOS type and whether or not to clean sound speed
  eos->eos_type = ghl_eos_tabulated;
  eos->table_type = table_type;
  eos->clean_sound_speed = clean_sound_speed;
  eos->enable_neural_net_c2p = enable_neural_net_c2p;
  eos->c2p_nn = NULL;

  ghl_error_codes_t err;

  if(!isfinite(rho_atm) || !isfinite(rho_min) || !isfinite(rho_max) || !isfinite(Y_e_atm)
     || !isfinite(Y_e_min) || !isfinite(Y_e_max) || !isfinite(T_atm) || !isfinite(T_min)
     || !isfinite(T_max)) {
    err = ghl_error_invalid_eos_parameters;
    goto cleanup;
  }

  // Step 2: Read the EOS table
  err = ghl_tabulated_read_table_set_EOS_params(table_filepath, eos);
  if(err != ghl_success) {
    goto cleanup;
  }

  if(eos->enable_neural_net_c2p) {
    ghl_info("Loading neural-network parameters from '%s'\n", table_filepath);
    err = ghl_c2p_nn_load_from_eos_hdf5(table_filepath, eos);
    if(err != ghl_success) {
      goto cleanup;
    }
    ghl_info("Loaded neural-network Con2Prim model: %d hidden layer(s) of width %d\n",
             eos->c2p_nn->n_hidden, eos->c2p_nn->hidden_dim);
  }

  // Step 3: Enforce default values for (rho, Y_e, T) min, max, and atm
  // Step 3.a: Atmosphere values
  if(rho_atm < 0) {
    err = ghl_error_invalid_rho_atm;
    goto cleanup;
  }
  if(Y_e_atm < 0) {
    err = ghl_error_invalid_Y_e_atm;
    goto cleanup;
  }
  if(T_atm < 0) {
    err = ghl_error_invalid_T_atm;
    goto cleanup;
  }

  // Step 3.b: Minimum values
  if(rho_min < eos->table_rho_min) {
    ghl_warn("Minimum density is less than table bounds; using table bounds (%.15e)\n", eos->table_rho_min);
    rho_min = eos->table_rho_min;
  }
  if(Y_e_min < eos->table_Y_e_min) {
    ghl_warn("Minimum electron fraction is less than table bounds; using table bounds (%.15e)\n", eos->table_Y_e_min);
    Y_e_min = eos->table_Y_e_min;
  }
  if(T_min < eos->table_T_min) {
    ghl_warn("Minimum temperature is less than table bounds; using table bounds (%.15e)\n", eos->table_T_min);
    T_min = eos->table_T_min;
  }

  // Step 3.c: Maximum values
  if(rho_max < 0 || rho_max > eos->table_rho_max) {
    ghl_warn("Invalid maximum density; using table bounds (%.15e)\n", eos->table_rho_max);
    rho_max = eos->table_rho_max;
  }
  if(Y_e_max < 0 || Y_e_max > eos->table_Y_e_max) {
    ghl_warn("Invalid maximum electron fraction; using table bounds (%.15e)\n", eos->table_Y_e_max);
    Y_e_max = eos->table_Y_e_max;
  }
  if(T_max < 0 || T_max > eos->table_T_max) {
    ghl_warn("Invalid maximum temperature; using table bounds (%.15e)\n", eos->table_T_max);
    T_max = eos->table_T_max;
  }

  // Step 3.d: Sanity check mins and maxs
  if(rho_min > rho_max) {
    err = ghl_error_rho_min_gt_rho_max;
    goto cleanup;
  }
  if(Y_e_min > Y_e_max) {
    err = ghl_error_Y_e_min_gt_Y_e_max;
    goto cleanup;
  }
  if(T_min > T_max) {
    err = ghl_error_T_min_gt_T_max;
    goto cleanup;
  }

  // Step 4: Initialize quantities which are common to all EOSs.
  init_common_eos_quantities;

  // Step 5: Set parameters specific to Tabulated EOS.
  eos->Y_e_atm = Y_e_atm;
  eos->Y_e_min = Y_e_min;
  eos->Y_e_max = Y_e_max;
  eos->T_atm   = T_atm;
  eos->T_min   = T_min;
  eos->T_max   = T_max;
  err = ghl_tabulated_compute_P_eps_S_from_T(eos,
                                             eos->rho_atm,
                                             Y_e_atm, T_atm,
                                             &eos->press_atm,
                                             &eos->eps_atm,
                                             &eos->entropy_atm);
  if(err != ghl_success) {
    goto cleanup;
  }

  // Step 6: These parameters are manually set here, but
  //         can be overwritten later.
  eos->root_finding_precision = 1e-10;

  // Step 7: Set minimum values for eps, P, and S
  eos->press_min   = eos->table_P_min;
  eos->press_max   = eos->table_P_max;
  eos->eps_min     = eos->table_eps_min;
  eos->eps_max     = eos->table_eps_max;
  eos->entropy_min = eos->table_ent_min;
  eos->entropy_max = eos->table_ent_max;
  eos->tau_atm     = eos->rho_min * eos->eps_min;

  // Step 8: Initialize beta-equilibrium arrays to NULL
  eos->Ye_of_lr = NULL;
  eos->lp_of_lr = NULL;
  eos->le_of_lr = NULL;
  eos->lh_of_lr = NULL;
  if(!isfinite(eos->rho_atm) || !isfinite(eos->rho_min) || !isfinite(eos->rho_max)
     || !isfinite(eos->Y_e_atm) || !isfinite(eos->Y_e_min) || !isfinite(eos->Y_e_max)
     || !isfinite(eos->T_atm) || !isfinite(eos->T_min) || !isfinite(eos->T_max)
     || !isfinite(eos->press_atm) || !isfinite(eos->press_min)
     || !isfinite(eos->press_max) || !isfinite(eos->eps_atm) || !isfinite(eos->eps_min)
     || !isfinite(eos->eps_max) || !isfinite(eos->entropy_atm)
     || !isfinite(eos->entropy_min) || !isfinite(eos->entropy_max)
     || !isfinite(eos->tau_atm)) {
    err = ghl_error_invalid_eos_parameters;
    goto cleanup;
  }
  *output = candidate;
  return ghl_success;

cleanup:
  ghl_tabulated_free_memory(&candidate);
  *output = (ghl_eos_parameters){ .eos_type = ghl_eos_tabulated };
  return err;
#endif
}

/*
 * Function    : ghl_initialize_hybrid_eos_functions_and_params()
 * Description : Fully initializes EOS struct elements for a hybrid EOS
*/
ghl_error_codes_t ghl_initialize_simple_eos_functions_and_params(
      const double rho_atm,
      double rho_min,
      double rho_max,
      const double press_atm,
      double press_min,
      double press_max,
      const double Gamma,
      ghl_eos_parameters *restrict eos) {

  // Step 1: Initialize Hybrid EOS functions
  ghl_initialize_eos_functions(ghl_eos_simple);

  // Step 2: Initialize Hybrid EOS parameters
  return ghl_initialize_simple_eos(rho_atm, rho_min, rho_max,
                                   press_atm, press_min, press_max,
                                   Gamma, eos);
}

/*
 * Function    : ghl_initialize_hybrid_eos_functions_and_params()
 * Description : Fully initializes EOS struct elements for a hybrid EOS
*/
ghl_error_codes_t ghl_initialize_hybrid_eos_functions_and_params(
      const double rho_atm,
      double rho_min,
      double rho_max,
      const int neos,
      const double *restrict rho_ppoly,
      const double *restrict Gamma_ppoly,
      const double K_ppoly0,
      const double Gamma_th,
      ghl_eos_parameters *restrict eos) {

  // Step 1: Initialize Hybrid EOS functions
  ghl_initialize_eos_functions(ghl_eos_hybrid);

  // Step 2: Initialize Hybrid EOS parameters
  return ghl_initialize_hybrid_eos(rho_atm, rho_min, rho_max,
                                   neos, rho_ppoly, Gamma_ppoly,
                                   K_ppoly0, Gamma_th, eos);
}

/* Function    : ghl_initialize_tabulated_eos()
 * Description : Initializes EOS struct elements for tabulated EOS
*/
ghl_error_codes_t ghl_initialize_tabulated_eos_functions_and_params(
      const char *table_filepath,
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      const double Ye_atm,
      const double Ye_min,
      const double Ye_max,
      const double T_atm,
      const double T_min,
      const double T_max,
      ghl_eos_parameters *restrict eos) {

  if(eos == NULL) {
    return ghl_error_eos_struct_is_null;
  }
#ifdef GHL_DISABLE_HDF5
  *eos = (ghl_eos_parameters){ .eos_type = ghl_eos_tabulated };
  return ghl_error_used_disabled_hdf5;
#else
  // FIXME: these are hard-coded default values for now
  const ghl_eos_table_t default_table_type = ghl_eos_table_stellarcollapse;
  const bool default_clean_sound_speed = false;
  const bool default_enable_neural_net_c2p = false;

  // Step 1: Initialize Tabulated EOS functions
  ghl_initialize_eos_functions(ghl_eos_tabulated);

  // Step 2: Initialize Tabulated EOS parameters
  return ghl_initialize_tabulated_eos(table_filepath,
                                      default_table_type,
                                      default_clean_sound_speed,
                                      default_enable_neural_net_c2p,
                                      rho_atm, rho_min, rho_max,
                                      Ye_atm, Ye_min, Ye_max,
                                      T_atm, T_min, T_max,
                                      eos);
#endif
}
