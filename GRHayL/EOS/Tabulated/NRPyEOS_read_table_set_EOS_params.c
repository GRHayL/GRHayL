#include "ghl_nrpyeos_tabulated.h"
#include "stellarcollapse/NRPyEOS_stellarcollapse.h"

static const char *table_str(ghl_eos_table_t type) {

#define GHL_EOS_CASE(name_) \
  case name_:               \
    return #name_;

  switch(type) {
    GHL_EOS_CASE(ghl_eos_table_stellarcollapse);
    default:
      return "Unknown";
  }

#undef GHL_EOS_CASE
}

static void ghl_eos_read_stellarcollapse_table(
      const char *filepath,
      ghl_eos_parameters *restrict eos,
      bool *restrict cs2_is_relativistic) {
  NRPyEOS_stellarcollapse_t *sc = NRPyEOS_stellarcollapse_read_table(filepath);
  if(sc == NULL) {
    ghl_error("Could not read EOS table from '%s'\n", filepath);
  }
  NRPyEOS_stellarcollapse_to_ghl(sc, eos);
  *cs2_is_relativistic = sc->cs2_is_relativistic;
  NRPyEOS_stellarcollapse_free_table(sc);
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

void NRPyEOS_read_table_set_EOS_params(
      const char *filepath,
      ghl_eos_parameters *restrict eos) {
#ifndef GHL_USE_HDF5
  GHL_HDF5_ERROR_IF_USED;
#else

  if(eos == NULL) {
    ghl_error("Input EOS parameter struct is NULL\n");
  }

  if(eos->eos_type != ghl_eos_tabulated) {
    ghl_error("EOS type is not tabulated.\n");
  }

  bool cs2_is_relativistic = false;
  switch(eos->table_type) {
    case ghl_eos_table_stellarcollapse:
      ghl_eos_read_stellarcollapse_table(filepath, eos, &cs2_is_relativistic);
      break;
    default:
      ghl_info("Unsupported EOS table type. Supported types:\n");
      for(ghl_eos_table_t n = 0; n < ghl_eos_table_types; n++) {
        ghl_info("  - %s\n", table_str(n));
      }
      ghl_error("Please set the EOS type to one of the known values above\n");
  }

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

  NRPyEOS_tabulate_enthalpy(eos);
  if(eos->clean_sound_speed) {
    NRPyEOS_tabulated_clean_sound_speed(eos, cs2_is_relativistic);
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

  eos->table_P_min = exp(get_EOS_table_min(eos, NRPyEOS_press_key));
  eos->table_P_max = exp(get_EOS_table_max(eos, NRPyEOS_press_key));
  eos->table_eps_min = exp(get_EOS_table_min(eos, NRPyEOS_eps_key)) - eos->energy_shift;
  eos->table_eps_max = exp(get_EOS_table_max(eos, NRPyEOS_eps_key)) - eos->energy_shift;
  eos->table_ent_min = get_EOS_table_min(eos, NRPyEOS_entropy_key);
  eos->table_ent_max = get_EOS_table_max(eos, NRPyEOS_entropy_key);
#endif
}
