#include "ghl_nrpyeos_tabulated.h"
#include "ghl_con2prim.h"

#include "stellarcollapse/NRPyEOS_stellarcollapse.h"
#include <limits.h>

#define GHL_CASE(name_) case name_: return #name_

static ghl_error_codes_t ghl_eos_read_stellarcollapse_table(
      const char *filepath,
      ghl_eos_parameters *restrict eos,
      bool *restrict cs2_is_relativistic) {

  NRPyEOS_stellarcollapse_t *sc = NULL;
  ghl_error_codes_t err = NRPyEOS_stellarcollapse_read_table(filepath, &sc);
  if(err != ghl_success) {
    return err;
  }
  *cs2_is_relativistic = sc->cs2_is_relativistic;
  err = NRPyEOS_stellarcollapse_to_ghl(sc, eos);
  NRPyEOS_stellarcollapse_free_table(sc);

  return err;
}

static inline double
get_EOS_table_max(const ghl_eos_parameters *restrict eos, const int var_key) {
  const int totalsize = eos->N_rho * eos->N_Ye * eos->N_T;
  double var_max_value = eos->table_all[var_key];

  for(int i = 1; i < totalsize; i++) {
    double var_value = eos->table_all[var_key + NRPyEOS_ntablekeys * i];
    if(var_value > var_max_value) {
      var_max_value = var_value;
    }
  }
  return var_max_value;
}

static inline double
get_EOS_table_min(const ghl_eos_parameters *restrict eos, const int var_key) {
  const int totalsize = eos->N_rho * eos->N_Ye * eos->N_T;
  double var_min_value = eos->table_all[var_key];

  for(int i = 1; i < totalsize; i++) {
    double var_value = eos->table_all[var_key + NRPyEOS_ntablekeys * i];
    if(var_value < var_min_value) {
      var_min_value = var_value;
    }
  }
  return var_min_value;
}

static ghl_error_codes_t validate_biased_axis(
      const double *restrict axis,
      const int n,
      const double inverse_spacing) {
  const double x0 = axis[0];
  const double values[2] = {
    (axis[0] - x0 - 1.0e-10) * inverse_spacing,
    (axis[n-1] - x0 - 1.0e-10) * inverse_spacing,
  };
  for(int i=0; i<2; i++) {
    if(!isfinite(values[i]) || values[i] < (double)INT_MIN
       || values[i] > (double)(INT_MAX - 1)) {
      return ghl_error_invalid_eos_table;
    }
  }
  return ghl_success;
}

ghl_error_codes_t NRPyEOS_read_table_set_EOS_params(
      const char *filepath,
      ghl_eos_parameters *restrict eos) {
#ifdef GHL_DISABLE_HDF5
  return ghl_error_used_disabled_hdf5;
#else

  if(eos == NULL) {
    return ghl_error_eos_struct_is_null;
  }

  if(eos->eos_type != ghl_eos_tabulated) {
    return ghl_error_invalid_eos_type;
  }

  bool cs2_is_relativistic = false;
  ghl_error_codes_t err = ghl_success;
  switch(eos->table_type) {
    case ghl_eos_table_stellarcollapse:
      err = ghl_eos_read_stellarcollapse_table(filepath, eos, &cs2_is_relativistic);
      break;
    default:
      return ghl_error_invalid_eos_table_type;
  }

  if(err != ghl_success) {
    NRPyEOS_free_memory(eos);
    return err;
  }

  const double dtemp = eos->table_logT[1] - eos->table_logT[0];
  const double drho = eos->table_logrho[1] - eos->table_logrho[0];
  const double dye = eos->table_Y_e[1] - eos->table_Y_e[0];

  eos->dtempi = 1.0 / dtemp;
  eos->drhoi = 1.0 / drho;
  eos->dyei = 1.0 / dye;
  eos->drhotempi = eos->drhoi * eos->dtempi;
  eos->drhoyei = eos->drhoi * eos->dyei;
  eos->dtempyei = eos->dtempi * eos->dyei;
  eos->drhotempyei = eos->drhoi * eos->dtempi * eos->dyei;

  eos->table_rho_max = exp(eos->table_logrho[eos->N_rho - 1]);
  eos->table_rho_min = exp(eos->table_logrho[0]);
  eos->table_T_max = exp(eos->table_logT[eos->N_T - 1]);
  eos->table_T_min = exp(eos->table_logT[0]);
  eos->table_Y_e_max = eos->table_Y_e[eos->N_Ye - 1];
  eos->table_Y_e_min = eos->table_Y_e[0];

  if(!isfinite(eos->table_rho_min) || eos->table_rho_min <= 0.0
     || !isfinite(eos->table_rho_max) || eos->table_rho_max <= 0.0
     || !isfinite(eos->table_T_min) || eos->table_T_min <= 0.0
     || !isfinite(eos->table_T_max) || eos->table_T_max <= 0.0
     || !isfinite(eos->dtempi) || eos->dtempi <= 0.0
     || !isfinite(eos->drhoi) || eos->drhoi <= 0.0
     || !isfinite(eos->dyei) || eos->dyei <= 0.0
     || !isfinite(eos->drhotempi) || eos->drhotempi <= 0.0
     || !isfinite(eos->drhoyei) || eos->drhoyei <= 0.0
     || !isfinite(eos->dtempyei) || eos->dtempyei <= 0.0
     || !isfinite(eos->drhotempyei) || eos->drhotempyei <= 0.0) {
    err = ghl_error_invalid_eos_table;
    goto cleanup;
  }
  err = validate_biased_axis(eos->table_logrho, eos->N_rho, eos->drhoi);
  if(err != ghl_success) goto cleanup;
  err = validate_biased_axis(eos->table_logT, eos->N_T, eos->dtempi);
  if(err != ghl_success) goto cleanup;
  err = validate_biased_axis(eos->table_Y_e, eos->N_Ye, eos->dyei);
  if(err != ghl_success) goto cleanup;

  const int npoints = eos->N_rho * eos->N_T * eos->N_Ye;
  for(int i = 0; i < npoints; i++) {
    int idx = NRPyEOS_press_key + NRPyEOS_ntablekeys * i;
    eos->table_all[idx] = eos->table_all[idx] * log(10.0) + log(CGS_TO_CODE_PRESSURE);

    idx = NRPyEOS_eps_key + NRPyEOS_ntablekeys * i;
    eos->table_all[idx] = eos->table_all[idx] * log(10.0) + log(CGS_TO_CODE_ENERGY);
    eos->table_eps[i] = exp(eos->table_all[idx]);

    idx = NRPyEOS_cs2_key + NRPyEOS_ntablekeys * i;
    eos->table_all[idx] *= CGS_TO_CODE_LENGTH * CGS_TO_CODE_LENGTH / CGS_TO_CODE_TIME
                           / CGS_TO_CODE_TIME;

    idx = NRPyEOS_depsdT_key + NRPyEOS_ntablekeys * i;
    eos->table_all[idx] *= CGS_TO_CODE_ENERGY;

    idx = NRPyEOS_dPdrho_key + NRPyEOS_ntablekeys * i;
    eos->table_all[idx] *= CGS_TO_CODE_PRESSURE / CGS_TO_CODE_DENSITY;

    idx = NRPyEOS_dPdeps_key + NRPyEOS_ntablekeys * i;
    eos->table_all[idx] *= CGS_TO_CODE_PRESSURE / CGS_TO_CODE_ENERGY;
  }

  err = NRPyEOS_tabulate_enthalpy(eos);
  if(err != ghl_success) goto cleanup;
  NRPyEOS_tabulated_adjust_sound_speed(eos, cs2_is_relativistic);

  eos->table_P_min = exp(get_EOS_table_min(eos, NRPyEOS_press_key));
  eos->table_P_max = exp(get_EOS_table_max(eos, NRPyEOS_press_key));
  eos->table_eps_min = exp(get_EOS_table_min(eos, NRPyEOS_eps_key)) - eos->energy_shift;
  eos->table_eps_max = exp(get_EOS_table_max(eos, NRPyEOS_eps_key)) - eos->energy_shift;
  eos->table_ent_min = get_EOS_table_min(eos, NRPyEOS_entropy_key);
  eos->table_ent_max = get_EOS_table_max(eos, NRPyEOS_entropy_key);

  return ghl_success;

cleanup:
  NRPyEOS_free_memory(eos);
  return err;
#endif
}
