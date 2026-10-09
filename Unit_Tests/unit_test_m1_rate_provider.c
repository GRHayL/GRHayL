#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include "../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_rate_backend.h"
#include "../GRHayL/Radiation/Neutrinos/ghl_m1_nrpyleakage_kernel.h"
#include "ghl_neutrino_rate_provider.h"
#include "ghl_radiation.h"
#include "m1_neutrino_rate_provider_reference.h"
#ifndef GHL_DISABLE_HDF5
#include "ghl_nrpyeos_tabulated.h"
#include <hdf5.h>
#endif
#include "ghl_nrpyleakage_nucleon_blocking.h"
#include "ghl_nrpyleakage_rate_helpers.h"
#include "m1_test_prng.h"

#include <errno.h>
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define ghl_neutrino_rate_provider_context m1_test_reference_provider_context
#define ghl_neutrino_rate_provider_cache   m1_test_reference_provider_cache
#define ghl_neutrino_rate_provider_initialize_reference \
  m1_test_reference_provider_initialize
#define ghl_neutrino_rate_provider_cache_initialize \
  m1_test_reference_provider_cache_initialize
#define ghl_neutrino_rate_provider_compute_cell m1_test_reference_provider_compute_cell
#include "m1_neutrino_rate_provider_reference.inc"
#undef ghl_neutrino_rate_provider_compute_cell
#undef ghl_neutrino_rate_provider_cache_initialize
#undef ghl_neutrino_rate_provider_initialize_reference
#undef ghl_neutrino_rate_provider_cache
#undef ghl_neutrino_rate_provider_context

/*
 * Deterministic property coverage for the M1 frozen-rate provider boundary.
 * The synthetic reference model is deliberately exercised with generated primitive
 * keys, channel masks, cache reuse, context-generation changes, and each
 * supported recovery policy.  A table path is optional; when supplied in an
 * HDF5 build, the same test calls the installed production provider.
 */

enum { PROVIDER_RANDOM_CASES = 128, TABLE_RANDOM_CASES = 32 };

static char *owned_provider_fixture_directory;
static char *owned_provider_fixture_path;
#ifndef GHL_DISABLE_HDF5
static ghl_eos_parameters *active_provider_fixture_eos;
#endif

static void cleanup_provider_fixture(void) {
#ifndef GHL_DISABLE_HDF5
  if(active_provider_fixture_eos != NULL) {
    ghl_tabulated_free_memory(active_provider_fixture_eos);
    active_provider_fixture_eos = NULL;
  }
#endif
  if(owned_provider_fixture_path != NULL) {
    if(remove(owned_provider_fixture_path) == 0 || errno == ENOENT) {
      free(owned_provider_fixture_path);
      owned_provider_fixture_path = NULL;
    }
    else {
      fprintf(
            stderr, "Could not remove provider fixture file %s: %s\n",
            owned_provider_fixture_path, strerror(errno));
    }
  }
  if(owned_provider_fixture_directory != NULL) {
    if(rmdir(owned_provider_fixture_directory) == 0 || errno == ENOENT) {
      free(owned_provider_fixture_directory);
      owned_provider_fixture_directory = NULL;
    }
    else {
      fprintf(
            stderr, "Could not remove provider fixture directory %s: %s\n",
            owned_provider_fixture_directory, strerror(errno));
    }
  }
}

static void provider_test_error(const char *restrict message) {
  cleanup_provider_fixture();
  ghl_error("%s\n", message);
}

static void require_condition(
      const bool condition,
      const char *restrict message,
      const int case_index) {
  if(!condition) {
    cleanup_provider_fixture();
    ghl_error("M1 rate-provider case %d: %s\n", case_index, message);
  }
}

static void require_error(
      const ghl_error_codes_t actual,
      const ghl_error_codes_t expected,
      const char *restrict operation,
      const int case_index) {
  if(actual != expected) {
    cleanup_provider_fixture();
    ghl_error(
          "M1 rate-provider case %d: %s returned %d, expected %d\n", case_index,
          operation, (int)actual, (int)expected);
  }
}

static bool same_rates(
      const ghl_m1_neutrino_rates *restrict lhs,
      const ghl_m1_neutrino_rates *restrict rhs) {
  if(lhs->species != rhs->species) {
    return false;
  }
  const double lhs_scalars[]
        = { lhs->eta_N,       lhs->eta_E,         lhs->kappa_a_N, lhs->kappa_a_E,
            lhs->kappa_s,     lhs->kappa_tr,      lhs->n_eq,      lhs->J_eq,
            lhs->mean_energy, lhs->lepton_weight, lhs->eta_N_cc,  lhs->kappa_a_N_cc };
  const double rhs_scalars[]
        = { rhs->eta_N,       rhs->eta_E,         rhs->kappa_a_N, rhs->kappa_a_E,
            rhs->kappa_s,     rhs->kappa_tr,      rhs->n_eq,      rhs->J_eq,
            rhs->mean_energy, rhs->lepton_weight, rhs->eta_N_cc,  rhs->kappa_a_N_cc };
  for(size_t i = 0; i < sizeof(lhs_scalars) / sizeof(lhs_scalars[0]); ++i) {
    if(lhs_scalars[i] != rhs_scalars[i]) {
      return false;
    }
  }
  for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
    if(lhs->eta_N_pair[process] != rhs->eta_N_pair[process]
       || lhs->eta_E_pair[process] != rhs->eta_E_pair[process]) {
      return false;
    }
  }
  return true;
}

static bool same_rate_bundle(
      const ghl_m1_neutrino_rates lhs[ghl_m1_neutrino_species_count],
      const ghl_m1_neutrino_rates rhs[ghl_m1_neutrino_species_count]) {
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    if(!same_rates(&lhs[species], &rhs[species])) {
      return false;
    }
  }
  return true;
}

static bool provider_values_close(const double lhs, const double rhs) {
  return fabs(lhs - rhs) <= 4096.0 * DBL_EPSILON * fmax(fabs(lhs), fabs(rhs));
}

static bool same_nrpyleakage_raw_species(
      const ghl_m1_nrpyleakage_species_raw_rates *restrict lhs,
      const ghl_m1_nrpyleakage_species_raw_rates *restrict rhs) {
  const double lhs_values[] = { lhs->neutrino_degeneracy,
                                lhs->F2,
                                lhs->F3,
                                lhs->F4,
                                lhs->F5,
                                lhs->n_eq_cgs,
                                lhs->J_eq_mev_cgs,
                                lhs->mean_energy_mev,
                                lhs->eta_N_beta_cgs,
                                lhs->eta_N_pair_cgs,
                                lhs->eta_N_plasmon_cgs,
                                lhs->eta_N_brems_cgs,
                                lhs->eta_E_beta_mev_cgs,
                                lhs->eta_E_pair_mev_cgs,
                                lhs->eta_E_plasmon_mev_cgs,
                                lhs->eta_E_brems_mev_cgs,
                                lhs->kappa_a_N_cc_cgs,
                                lhs->kappa_a_E_cc_cgs,
                                lhs->kappa_s_N_neutron_cgs,
                                lhs->kappa_s_N_proton_cgs,
                                lhs->kappa_s_E_neutron_cgs,
                                lhs->kappa_s_E_proton_cgs };
  const double rhs_values[] = { rhs->neutrino_degeneracy,
                                rhs->F2,
                                rhs->F3,
                                rhs->F4,
                                rhs->F5,
                                rhs->n_eq_cgs,
                                rhs->J_eq_mev_cgs,
                                rhs->mean_energy_mev,
                                rhs->eta_N_beta_cgs,
                                rhs->eta_N_pair_cgs,
                                rhs->eta_N_plasmon_cgs,
                                rhs->eta_N_brems_cgs,
                                rhs->eta_E_beta_mev_cgs,
                                rhs->eta_E_pair_mev_cgs,
                                rhs->eta_E_plasmon_mev_cgs,
                                rhs->eta_E_brems_mev_cgs,
                                rhs->kappa_a_N_cc_cgs,
                                rhs->kappa_a_E_cc_cgs,
                                rhs->kappa_s_N_neutron_cgs,
                                rhs->kappa_s_N_proton_cgs,
                                rhs->kappa_s_E_neutron_cgs,
                                rhs->kappa_s_E_proton_cgs };
  for(size_t i = 0; i < sizeof(lhs_values) / sizeof(lhs_values[0]); ++i) {
    if(lhs_values[i] != rhs_values[i]) {
      return false;
    }
  }
  return true;
}

static bool same_nrpyleakage_raw_rates(
      const ghl_m1_nrpyleakage_raw_rates *restrict lhs,
      const ghl_m1_nrpyleakage_raw_rates *restrict rhs) {
  if(lhs->nux_single_species_multiplicity != rhs->nux_single_species_multiplicity) {
    return false;
  }
  for(int species = 0; species < ghl_m1_nrpyleakage_species_count; ++species) {
    if(!same_nrpyleakage_raw_species(&lhs->species[species], &rhs->species[species])) {
      return false;
    }
  }
  return true;
}

static double nrpyleakage_fraction_roundoff_envelope(void) {
  /* Keep this oracle identical to ghl_nrpyleakage_normalize_nucleon_fractions(). */
  const double gamma_64 = 64.0 * DBL_EPSILON / (1.0 - 64.0 * DBL_EPSILON);
  return 27.0 * gamma_64;
}

static void require_value_close(
      const double actual,
      const double expected,
      const char *restrict label,
      const int case_index) {
  require_condition(provider_values_close(actual, expected), label, case_index);
}

#ifndef GHL_DISABLE_HDF5
enum {
  PROVIDER_FIXTURE_NRHO = 3,
  PROVIDER_FIXTURE_NTEMP = 3,
  PROVIDER_FIXTURE_NYE = 3,
  PROVIDER_FIXTURE_CELL_COUNT = PROVIDER_FIXTURE_NRHO * PROVIDER_FIXTURE_NTEMP
        * PROVIDER_FIXTURE_NYE
};

/* Match the EOS interpolator's biased cell selection for the midpoint used
 * by the corruption witnesses below. The first cell is not the midpoint
 * cell of an arbitrary external table. No private EOS header is needed. */
static void provider_midpoint_table_cell(
      const ghl_eos_parameters *restrict eos, int lower[3]) {
  const double logrho = log(sqrt(eos->table_rho_min * eos->table_rho_max));
  const double logT = log(sqrt(eos->table_T_min * eos->table_T_max));
  const double Ye = 0.5 * (eos->table_Y_e_min + eos->table_Y_e_max);
  lower[0] = ghl_iclamp(
        1 + (int)((logrho - eos->table_logrho[0] - 1.e-10) * eos->drhoi),
        1, eos->N_rho - 1) - 1;
  lower[1] = ghl_iclamp(
        1 + (int)((logT - eos->table_logT[0] - 1.e-10) * eos->dtempi),
        1, eos->N_T - 1) - 1;
  lower[2] = ghl_iclamp(
        1 + (int)((Ye - eos->table_Y_e[0] - 1.e-10) * eos->dyei),
        1, eos->N_Ye - 1) - 1;
}

/* Keep malformed-table witnesses in memory only.  Each witness changes the
 * eight interpolation corners used by the authenticated interior point and
 * restores every byte before the next production-table assertion. */
static void provider_set_table_corners(
      ghl_eos_parameters *restrict eos,
      const int key,
      const double replacement,
      double saved[8]) {
  int lower[3];
  provider_midpoint_table_cell(eos, lower);
  int corner = 0;
  for(int ir = 0; ir < 2; ++ir) {
    for(int it = 0; it < 2; ++it) {
      for(int iy = 0; iy < 2; ++iy) {
        const size_t index = (size_t)NRPYEOS_IDX3D(
              eos, lower[0] + ir, lower[1] + it, lower[2] + iy, key);
        saved[corner++] = eos->table_all[index];
        eos->table_all[index] = replacement;
      }
    }
  }
}

static void provider_restore_table_corners(
      ghl_eos_parameters *restrict eos,
      const int key,
      const double saved[8]) {
  int lower[3];
  provider_midpoint_table_cell(eos, lower);
  int corner = 0;
  for(int ir = 0; ir < 2; ++ir) {
    for(int it = 0; it < 2; ++it) {
      for(int iy = 0; iy < 2; ++iy) {
        const size_t index = (size_t)NRPYEOS_IDX3D(
              eos, lower[0] + ir, lower[1] + it, lower[2] + iy, key);
        eos->table_all[index] = saved[corner++];
      }
    }
  }
}

static bool write_provider_fixture_scalar_int(
      const hid_t file,
      const char *restrict name,
      const int value) {
  const hsize_t dimensions[1] = { 1 };
  const hid_t space = H5Screate_simple(1, dimensions, NULL);
  if(space < 0) {
    return false;
  }
  const hid_t dataset = H5Dcreate2(
        file, name, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  const bool ok
        = dataset >= 0
          && H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value)
                   >= 0;
  if(dataset >= 0) {
    H5Dclose(dataset);
  }
  H5Sclose(space);
  return ok;
}

static bool write_provider_fixture_scalar_double(
      const hid_t file,
      const char *restrict name,
      const double value) {
  const hsize_t dimensions[1] = { 1 };
  const hid_t space = H5Screate_simple(1, dimensions, NULL);
  if(space < 0) {
    return false;
  }
  const hid_t dataset = H5Dcreate2(
        file, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  const bool ok
        = dataset >= 0
          && H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value)
                   >= 0;
  if(dataset >= 0) {
    H5Dclose(dataset);
  }
  H5Sclose(space);
  return ok;
}

static bool write_provider_fixture_vector(
      const hid_t file,
      const char *restrict name,
      const double *restrict values,
      const int count) {
  const hsize_t dimensions[1] = { (hsize_t)count };
  const hid_t space = H5Screate_simple(1, dimensions, NULL);
  if(space < 0) {
    return false;
  }
  const hid_t dataset = H5Dcreate2(
        file, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  const bool ok
        = dataset >= 0
          && H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values)
                   >= 0;
  if(dataset >= 0) {
    H5Dclose(dataset);
  }
  H5Sclose(space);
  return ok;
}

static bool write_provider_fixture_table(
      const hid_t file,
      const char *restrict name,
      const double values[PROVIDER_FIXTURE_CELL_COUNT]) {
  const hsize_t dimensions[3]
        = { PROVIDER_FIXTURE_NYE, PROVIDER_FIXTURE_NTEMP, PROVIDER_FIXTURE_NRHO };
  const hid_t space = H5Screate_simple(3, dimensions, NULL);
  if(space < 0) {
    return false;
  }
  const hid_t dataset = H5Dcreate2(
        file, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  const bool ok
        = dataset >= 0
          && H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values)
                   >= 0;
  if(dataset >= 0) {
    H5Dclose(dataset);
  }
  H5Sclose(space);
  return ok;
}

/* This is the small authenticated StellarCollapse fixture used by the
 * offline provider campaign, kept test-local so table coverage has no /work
 * or external-file dependency. */
static bool create_provider_fixture(
      const bool cold_degenerate_fixture,
      const bool high_temperature_fixture) {
  const char *const tmpdir = getenv("TMPDIR");
  const char *const base = (tmpdir == NULL || tmpdir[0] == '\0') ? "/tmp" : tmpdir;
  const size_t base_length = strlen(base);
  const bool needs_separator = base_length != 0 && base[base_length - 1] != '/';
  static const char directory_suffix[] = "ghl_m1_rate_provider_fixture.XXXXXX";
  const size_t separator_length = needs_separator ? 1 : 0;
  if(base_length > SIZE_MAX - separator_length
     || base_length + separator_length > SIZE_MAX - sizeof(directory_suffix)) {
    return false;
  }
  const size_t directory_size
        = base_length + separator_length + sizeof(directory_suffix);
  owned_provider_fixture_directory = malloc(directory_size);
  if(owned_provider_fixture_directory == NULL) {
    return false;
  }
  const int directory_characters = snprintf(
        owned_provider_fixture_directory, directory_size, "%s%s%s", base,
        needs_separator ? "/" : "", directory_suffix);
  if(directory_characters < 0 || (size_t)directory_characters >= directory_size
     || mkdtemp(owned_provider_fixture_directory) == NULL) {
    free(owned_provider_fixture_directory);
    owned_provider_fixture_directory = NULL;
    cleanup_provider_fixture();
    return false;
  }
  const size_t directory_length = strlen(owned_provider_fixture_directory);
  static const char file_suffix[] = "/table.h5";
  if(directory_length > SIZE_MAX - sizeof(file_suffix)) {
    cleanup_provider_fixture();
    return false;
  }
  const size_t path_size = directory_length + sizeof(file_suffix);
  owned_provider_fixture_path = malloc(path_size);
  if(owned_provider_fixture_path == NULL) {
    cleanup_provider_fixture();
    return false;
  }
  const int path_characters = snprintf(
        owned_provider_fixture_path, path_size, "%s%s", owned_provider_fixture_directory,
        file_suffix);
  if(path_characters < 0 || (size_t)path_characters >= path_size) {
    cleanup_provider_fixture();
    return false;
  }

  const double logrho[PROVIDER_FIXTURE_NRHO] = { 10.0, 11.0, 12.0 };
  const double logtemp[PROVIDER_FIXTURE_NTEMP] = {
    high_temperature_fixture ? 39.0 : (cold_degenerate_fixture ? log10(0.05) : 0.0),
    high_temperature_fixture ? 39.5 : (cold_degenerate_fixture ? log10(0.1) : 0.5),
    high_temperature_fixture ? 40.0 : (cold_degenerate_fixture ? log10(0.2) : 1.0)
  };
  const double ye[PROVIDER_FIXTURE_NYE] = { 0.1, 0.5, 0.9 };
  double abar[PROVIDER_FIXTURE_CELL_COUNT], xa[PROVIDER_FIXTURE_CELL_COUNT];
  double xh[PROVIDER_FIXTURE_CELL_COUNT], xn[PROVIDER_FIXTURE_CELL_COUNT];
  double xp[PROVIDER_FIXTURE_CELL_COUNT], zbar[PROVIDER_FIXTURE_CELL_COUNT];
  double cs2[PROVIDER_FIXTURE_CELL_COUNT], dedt[PROVIDER_FIXTURE_CELL_COUNT];
  double dpderho[PROVIDER_FIXTURE_CELL_COUNT];
  double dpdrhoe[PROVIDER_FIXTURE_CELL_COUNT];
  double entropy[PROVIDER_FIXTURE_CELL_COUNT], gamma[PROVIDER_FIXTURE_CELL_COUNT];
  double logenergy[PROVIDER_FIXTURE_CELL_COUNT];
  double logpress[PROVIDER_FIXTURE_CELL_COUNT];
  double mu_e[PROVIDER_FIXTURE_CELL_COUNT], mu_n[PROVIDER_FIXTURE_CELL_COUNT];
  double mu_p[PROVIDER_FIXTURE_CELL_COUNT], muhat[PROVIDER_FIXTURE_CELL_COUNT];
  double munu[PROVIDER_FIXTURE_CELL_COUNT];
  for(int iy = 0; iy < PROVIDER_FIXTURE_NYE; ++iy) {
    for(int it = 0; it < PROVIDER_FIXTURE_NTEMP; ++it) {
      for(int ir = 0; ir < PROVIDER_FIXTURE_NRHO; ++ir) {
        const int index
              = ir + PROVIDER_FIXTURE_NRHO * (it + PROVIDER_FIXTURE_NTEMP * iy);
        const double variation = 0.11 * ir + 0.23 * it + 0.37 * iy;
        abar[index] = 56.0 + variation;
        xa[index] = 0.1 + 0.01 * variation;
        xh[index] = 0.2 + 0.01 * variation;
        xn[index] = 1.0 - ye[iy];
        xp[index] = ye[iy];
        zbar[index] = 26.0 + 0.1 * variation;
        cs2[index] = 1.0e20 + 1.0e18 * variation;
        dedt[index] = 1.0e18 + 1.0e16 * variation;
        dpderho[index] = 1.0e15 + 1.0e13 * variation;
        dpdrhoe[index] = 1.0e10 + 1.0e8 * variation;
        entropy[index] = 1.0 + variation;
        gamma[index] = 1.5 + 0.01 * variation;
        logenergy[index] = 18.0 + 0.01 * variation;
        logpress[index] = 25.0 + 0.01 * variation;
        mu_e[index]
              = cold_degenerate_fixture ? 40.0 : 3.0 + 0.20 * ir + 0.30 * it + 0.40 * iy;
        mu_n[index]
              = cold_degenerate_fixture ? 40.0 : 1.0 + 0.10 * ir + 0.12 * it + 0.15 * iy;
        mu_p[index]
              = cold_degenerate_fixture ? 0.0 : 0.2 + 0.05 * ir + 0.06 * it + 0.08 * iy;
        muhat[index] = mu_n[index] - mu_p[index];
        munu[index] = 0.1 + 0.02 * ir + 0.03 * it + 0.04 * iy;
      }
    }
  }

  const hid_t file = H5Fcreate(
        owned_provider_fixture_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if(file < 0) {
    cleanup_provider_fixture();
    return false;
  }
  bool ok
        = write_provider_fixture_scalar_int(file, "have_rel_cs2", 1)
          && write_provider_fixture_scalar_int(file, "pointsrho", PROVIDER_FIXTURE_NRHO)
          && write_provider_fixture_scalar_int(
                file, "pointstemp", PROVIDER_FIXTURE_NTEMP)
          && write_provider_fixture_scalar_int(file, "pointsye", PROVIDER_FIXTURE_NYE)
          && write_provider_fixture_scalar_double(file, "energy_shift", 0.0)
          && write_provider_fixture_vector(file, "logrho", logrho, PROVIDER_FIXTURE_NRHO)
          && write_provider_fixture_vector(
                file, "logtemp", logtemp, PROVIDER_FIXTURE_NTEMP)
          && write_provider_fixture_vector(file, "ye", ye, PROVIDER_FIXTURE_NYE);
#define WRITE_PROVIDER_FIXTURE_TABLE(name, values) \
  ok = ok && write_provider_fixture_table(file, name, values)
  WRITE_PROVIDER_FIXTURE_TABLE("Abar", abar);
  WRITE_PROVIDER_FIXTURE_TABLE("Xa", xa);
  WRITE_PROVIDER_FIXTURE_TABLE("Xh", xh);
  WRITE_PROVIDER_FIXTURE_TABLE("Xn", xn);
  WRITE_PROVIDER_FIXTURE_TABLE("Xp", xp);
  WRITE_PROVIDER_FIXTURE_TABLE("Zbar", zbar);
  WRITE_PROVIDER_FIXTURE_TABLE("cs2", cs2);
  WRITE_PROVIDER_FIXTURE_TABLE("dedt", dedt);
  WRITE_PROVIDER_FIXTURE_TABLE("dpderho", dpderho);
  WRITE_PROVIDER_FIXTURE_TABLE("dpdrhoe", dpdrhoe);
  WRITE_PROVIDER_FIXTURE_TABLE("entropy", entropy);
  WRITE_PROVIDER_FIXTURE_TABLE("gamma", gamma);
  WRITE_PROVIDER_FIXTURE_TABLE("logenergy", logenergy);
  WRITE_PROVIDER_FIXTURE_TABLE("logpress", logpress);
  WRITE_PROVIDER_FIXTURE_TABLE("mu_e", mu_e);
  WRITE_PROVIDER_FIXTURE_TABLE("mu_n", mu_n);
  WRITE_PROVIDER_FIXTURE_TABLE("mu_p", mu_p);
  WRITE_PROVIDER_FIXTURE_TABLE("muhat", muhat);
  WRITE_PROVIDER_FIXTURE_TABLE("munu", munu);
#undef WRITE_PROVIDER_FIXTURE_TABLE
  if(H5Fclose(file) < 0) {
    ok = false;
  }
  if(!ok) {
    cleanup_provider_fixture();
    return false;
  }
  return true;
}

#endif

static void validate_rate_bundle(
      const ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count],
      const int case_index) {
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    const ghl_error_codes_t error
          = ghl_m1_validate_neutrino_rates(&rates[species], NULL);
    require_error(error, ghl_success, "published rate validation", case_index);
    require_condition(
          rates[species].species == (ghl_m1_neutrino_species_t)species,
          "published rate species is not indexed consistently", case_index);
    require_condition(
          isfinite(rates[species].mean_energy) && rates[species].mean_energy > 0.0,
          "published mean energy is not positive and finite", case_index);
  }
}

static void test_nrpyleakage_raw_kernel(void) {
  ghl_m1_nrpyleakage_thermo_state thermo;
  require_error(
        ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
              1.0e-4, 0.20, 8.0, 12.0, 18.0, 5.0, 20.0, 0.70, 0.20, &thermo),
        ghl_success, "direct NRPyLeakage thermo construction", 3070);
  const double eta[ghl_m1_nrpyleakage_species_count] = { 0.30, -0.20, 0.10 };
  ghl_m1_nrpyleakage_raw_rates raw;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(&thermo, eta, &raw),
        ghl_success, "direct NRPyLeakage strict raw rates", 3071);
  require_condition(
        raw.nux_single_species_multiplicity == 1,
        "direct NRPyLeakage raw multiplicity is not one", 3071);
  for(int species = 0; species < ghl_m1_nrpyleakage_species_count; ++species) {
    require_condition(
          isfinite(raw.species[species].mean_energy_mev)
                && raw.species[species].mean_energy_mev > 0.0,
          "direct NRPyLeakage mean energy is not positive", 3071);
  }

  ghl_m1_nrpyleakage_thermo_state roundoff_thermo = thermo;
  roundoff_thermo.X_n = -nrpyleakage_fraction_roundoff_envelope();
  ghl_m1_nrpyleakage_raw_rates roundoff_raw;
  ghl_m1_nrpyleakage_raw_rates normalized_raw;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(
              &roundoff_thermo, eta, &roundoff_raw),
        ghl_success, "direct NRPyLeakage roundoff composition", 3071);
  ghl_m1_nrpyleakage_thermo_state normalized_thermo = roundoff_thermo;
  normalized_thermo.X_n = 0.0;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(
              &normalized_thermo, eta, &normalized_raw),
        ghl_success, "direct NRPyLeakage exact normalized composition", 3071);
  require_condition(
        same_nrpyleakage_raw_rates(&roundoff_raw, &normalized_raw),
        "roundoff composition raw rates differ from exact normalized rates", 3071);

  /* Use the maintained NRPyLeakage helpers as the oracle for every changed
   * local microphysics component.  The M1 adapter is expected to retain only
   * its channel prefactors and unit bookkeeping around these helpers. */
  const double rho_cgs = thermo.rho * NRPyLeakage_units_geom_to_cgs_D;
  double B_n, B_p, Y_np, Y_pn, eta_n_minus_eta_p;
  require_error(
        NRPyLeakage_compute_nucleon_blocking(
              rho_cgs, thermo.T, thermo.X_n, thermo.X_p, &B_n, &B_p, &Y_np, &Y_pn,
              &eta_n_minus_eta_p),
        ghl_success, "canonical mixed-composition nucleon blocking", 3072);
  const double reaction_shift
        = nrpyl_compute_reaction_shift(thermo.T, thermo.muhat, eta_n_minus_eta_p);
  require_condition(
        isfinite(reaction_shift) && fabs(reaction_shift - thermo.muhat) > 1.0e-12,
        "canonical kinetic reaction shift was not exercised", 3072);

  nrpyl_beta_moments beta_nue_emission, beta_anue_emission;
  nrpyl_beta_moments beta_nue_absorption, beta_anue_absorption;
  require_error(
        nrpyl_compute_beta_emission_moments(
              thermo.T, thermo.mu_e, eta[ghl_m1_nrpyleakage_nue], 1, reaction_shift,
              &beta_nue_emission),
        ghl_success, "canonical nue emission moments", 3072);
  require_error(
        nrpyl_compute_beta_emission_moments(
              thermo.T, thermo.mu_e, eta[ghl_m1_nrpyleakage_anue], -1, reaction_shift,
              &beta_anue_emission),
        ghl_success, "canonical anue emission moments", 3072);
  require_error(
        nrpyl_compute_beta_absorption_moments(
              thermo.T, thermo.mu_e, eta[ghl_m1_nrpyleakage_nue], 1, reaction_shift,
              &beta_nue_absorption),
        ghl_success, "canonical nue absorption moments", 3072);
  require_error(
        nrpyl_compute_beta_absorption_moments(
              thermo.T, thermo.mu_e, eta[ghl_m1_nrpyleakage_anue], -1, reaction_shift,
              &beta_anue_absorption),
        ghl_success, "canonical anue absorption moments", 3072);

  const double beta_prefactor
        = 8.0 * NRPyLeakage_N_A * NRPyLeakage_beta * pow(thermo.T, 5) * rho_cgs
          * (M_PI / NRPyLeakage_hc3)
          * ((3.0 / 8.0) * NRPyLeakage_alpha * NRPyLeakage_alpha + 1.0 / 8.0);
  const double opacity_prefactor = NRPyLeakage_N_A * NRPyLeakage_sigma_0 * thermo.T
                                   * thermo.T * rho_cgs
                                   / (NRPyLeakage_m_e_c2 * NRPyLeakage_m_e_c2);
  const double cc_prefactor
        = opacity_prefactor
          * ((3.0 / 4.0) * NRPyLeakage_alpha * NRPyLeakage_alpha + 1.0 / 4.0);
  const double neutron_scattering_factor
        = B_n * ((5.0 / 24.0) * NRPyLeakage_alpha * NRPyLeakage_alpha + 1.0 / 24.0);
  const double proton_scattering_factor
        = B_p
          * ((1.0 / 6.0) * (NRPyLeakage_C_V - 1.0) * (NRPyLeakage_C_V - 1.0)
             + (5.0 / 24.0) * NRPyLeakage_alpha * NRPyLeakage_alpha);

  double nue_F2, nue_F3, nue_F4, nue_F5;
  double anue_F2, anue_F3, anue_F4, anue_F5;
  double nux_F2, nux_F3, nux_F4, nux_F5;
  require_error(
        NRPyLeakage_Fermi_Dirac_integrals(2, eta[0], &nue_F2), ghl_success,
        "canonical nue F2 oracle", 3072);
  require_error(
        NRPyLeakage_Fermi_Dirac_integrals(3, eta[0], &nue_F3), ghl_success,
        "canonical nue F3 oracle", 3072);
  require_error(
        NRPyLeakage_Fermi_Dirac_integrals(4, eta[0], &nue_F4), ghl_success,
        "canonical nue F4 oracle", 3072);
  require_error(
        NRPyLeakage_Fermi_Dirac_integrals(5, eta[0], &nue_F5), ghl_success,
        "canonical nue F5 oracle", 3072);
  require_error(
        NRPyLeakage_Fermi_Dirac_integrals(2, eta[1], &anue_F2), ghl_success,
        "canonical anue F2 oracle", 3072);
  require_error(
        NRPyLeakage_Fermi_Dirac_integrals(3, eta[1], &anue_F3), ghl_success,
        "canonical anue F3 oracle", 3072);
  require_error(
        NRPyLeakage_Fermi_Dirac_integrals(4, eta[1], &anue_F4), ghl_success,
        "canonical anue F4 oracle", 3072);
  require_error(
        NRPyLeakage_Fermi_Dirac_integrals(5, eta[1], &anue_F5), ghl_success,
        "canonical anue F5 oracle", 3072);
  require_error(
        NRPyLeakage_Fermi_Dirac_integrals(2, 0.0, &nux_F2), ghl_success,
        "canonical nux F2 oracle", 3072);
  require_error(
        NRPyLeakage_Fermi_Dirac_integrals(3, 0.0, &nux_F3), ghl_success,
        "canonical nux F3 oracle", 3072);
  require_error(
        NRPyLeakage_Fermi_Dirac_integrals(4, 0.0, &nux_F4), ghl_success,
        "canonical nux F4 oracle", 3072);
  require_error(
        NRPyLeakage_Fermi_Dirac_integrals(5, 0.0, &nux_F5), ghl_success,
        "canonical nux F5 oracle", 3072);

  const double expected_nue_scattering_n
        = opacity_prefactor * neutron_scattering_factor * nue_F4 / nue_F2;
  const double expected_nue_scattering_p
        = opacity_prefactor * proton_scattering_factor * nue_F4 / nue_F2;
  const double expected_anue_scattering_n
        = opacity_prefactor * neutron_scattering_factor * anue_F4 / anue_F2;
  const double expected_anue_scattering_p
        = opacity_prefactor * proton_scattering_factor * anue_F4 / anue_F2;
  const double expected_nux_scattering_n
        = opacity_prefactor * neutron_scattering_factor * nux_F4 / nux_F2;
  const double expected_nux_scattering_p
        = opacity_prefactor * proton_scattering_factor * nux_F4 / nux_F2;
  const double expected_nue_energy_n
        = opacity_prefactor * neutron_scattering_factor * nue_F5 / nue_F3;
  const double expected_nue_energy_p
        = opacity_prefactor * proton_scattering_factor * nue_F5 / nue_F3;
  const double expected_anue_energy_n
        = opacity_prefactor * neutron_scattering_factor * anue_F5 / anue_F3;
  const double expected_anue_energy_p
        = opacity_prefactor * proton_scattering_factor * anue_F5 / anue_F3;
  const double expected_nux_energy_n
        = opacity_prefactor * neutron_scattering_factor * nux_F5 / nux_F3;
  const double expected_nux_energy_p
        = opacity_prefactor * proton_scattering_factor * nux_F5 / nux_F3;
  require_value_close(
        raw.species[0].kappa_s_N_neutron_cgs, expected_nue_scattering_n,
        "canonical neutron blocking was not used for nue scattering", 3072);
  require_value_close(
        raw.species[0].kappa_s_N_proton_cgs, expected_nue_scattering_p,
        "canonical proton blocking was not used for nue scattering", 3072);
  require_value_close(
        raw.species[1].kappa_s_N_neutron_cgs, expected_anue_scattering_n,
        "canonical neutron blocking was not used for anue scattering", 3072);
  require_value_close(
        raw.species[1].kappa_s_N_proton_cgs, expected_anue_scattering_p,
        "canonical proton blocking was not used for anue scattering", 3072);
  require_value_close(
        raw.species[2].kappa_s_N_neutron_cgs, expected_nux_scattering_n,
        "canonical neutron blocking was not used for nux scattering", 3072);
  require_value_close(
        raw.species[2].kappa_s_N_proton_cgs, expected_nux_scattering_p,
        "canonical proton blocking was not used for nux scattering", 3072);
  require_value_close(
        raw.species[0].kappa_s_E_neutron_cgs, expected_nue_energy_n,
        "canonical neutron blocking was not used for nue energy scattering", 3072);
  require_value_close(
        raw.species[0].kappa_s_E_proton_cgs, expected_nue_energy_p,
        "canonical proton blocking was not used for nue energy scattering", 3072);
  require_value_close(
        raw.species[1].kappa_s_E_neutron_cgs, expected_anue_energy_n,
        "canonical neutron blocking was not used for anue energy scattering", 3072);
  require_value_close(
        raw.species[1].kappa_s_E_proton_cgs, expected_anue_energy_p,
        "canonical proton blocking was not used for anue energy scattering", 3072);
  require_value_close(
        raw.species[2].kappa_s_E_neutron_cgs, expected_nux_energy_n,
        "canonical neutron blocking was not used for nux energy scattering", 3072);
  require_value_close(
        raw.species[2].kappa_s_E_proton_cgs, expected_nux_energy_p,
        "canonical proton blocking was not used for nux energy scattering", 3072);

  require_value_close(
        raw.species[0].eta_N_beta_cgs, Y_pn * beta_prefactor * beta_nue_emission.number,
        "canonical nue emission number moment was not used", 3072);
  require_value_close(
        raw.species[1].eta_N_beta_cgs, Y_np * beta_prefactor * beta_anue_emission.number,
        "canonical anue emission number moment was not used", 3072);
  require_value_close(
        raw.species[0].eta_E_beta_mev_cgs,
        thermo.T * Y_pn * beta_prefactor * beta_nue_emission.energy,
        "canonical nue emission energy moment was not used", 3072);
  require_value_close(
        raw.species[1].eta_E_beta_mev_cgs,
        thermo.T * Y_np * beta_prefactor * beta_anue_emission.energy,
        "canonical anue emission energy moment was not used", 3072);
  require_value_close(
        raw.species[0].kappa_a_N_cc_cgs,
        cc_prefactor * Y_np * beta_nue_absorption.number,
        "canonical nue absorption number moment was not used", 3072);
  require_value_close(
        raw.species[1].kappa_a_N_cc_cgs,
        cc_prefactor * Y_pn * beta_anue_absorption.number,
        "canonical anue absorption number moment was not used", 3072);
  require_value_close(
        raw.species[0].kappa_a_E_cc_cgs,
        cc_prefactor * Y_np * beta_nue_absorption.energy,
        "canonical nue absorption energy moment was not used", 3072);
  require_value_close(
        raw.species[1].kappa_a_E_cc_cgs,
        cc_prefactor * Y_pn * beta_anue_absorption.energy,
        "canonical anue absorption energy moment was not used", 3072);

  const double expected_brems_number
        = nrpyl_bremsstrahlung_number_rate(thermo.T, rho_cgs, thermo.X_n, thermo.X_p);
  const double expected_brems_energy
        = nrpyl_bremsstrahlung_energy_rate(thermo.T, expected_brems_number);
  for(int species = 0; species < ghl_m1_nrpyleakage_species_count; ++species) {
    require_value_close(
          raw.species[species].eta_N_brems_cgs, expected_brems_number,
          "canonical bremsstrahlung number helper was not used", 3072);
    require_value_close(
          raw.species[species].eta_E_brems_mev_cgs, expected_brems_energy,
          "canonical bremsstrahlung energy helper was not used", 3072);
  }

  ghl_m1_nrpyleakage_thermo_state doubled_density = thermo;
  doubled_density.rho *= 2.0;
  ghl_m1_nrpyleakage_raw_rates doubled_raw;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(
              &doubled_density, eta, &doubled_raw),
        ghl_success, "doubled-density NRPyLeakage raw rates", 3073);
  require_value_close(
        doubled_raw.species[0].eta_N_brems_cgs / raw.species[0].eta_N_brems_cgs, 4.0,
        "bremsstrahlung number rate is not quadratic in density", 3073);
  require_value_close(
        doubled_raw.species[0].eta_E_brems_mev_cgs / raw.species[0].eta_E_brems_mev_cgs,
        4.0, "bremsstrahlung energy rate is not quadratic in density", 3073);

  /* Exact one-species endpoints must take the helper's analytic zero-product
   * limit.  DBL_MAX chemical shifts make any accidental reaction-shift
   * evaluation observable as a failure without changing the endpoint oracle. */
  const struct {
    const char *label;
    double X_n;
    double X_p;
    double Ye;
    double muhat;
  } endpoints[] = { { "pure-neutron beta endpoint", 1.0, 0.0, 0.0, DBL_MAX },
                    { "pure-proton beta endpoint", 0.0, 1.0, 1.0, -DBL_MAX } };
  for(size_t endpoint = 0; endpoint < sizeof(endpoints) / sizeof(endpoints[0]);
      ++endpoint) {
    ghl_m1_nrpyleakage_thermo_state endpoint_thermo;
    require_error(
          ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
                thermo.rho, endpoints[endpoint].Ye, thermo.T, endpoints[endpoint].muhat,
                thermo.mu_e, thermo.mu_p, thermo.mu_n, endpoints[endpoint].X_n,
                endpoints[endpoint].X_p, &endpoint_thermo),
          ghl_success, endpoints[endpoint].label, 3078 + (int)endpoint);
    double endpoint_B_n, endpoint_B_p, endpoint_Y_np, endpoint_Y_pn;
    double endpoint_eta_difference;
    require_error(
          NRPyLeakage_compute_nucleon_blocking(
                endpoint_thermo.rho * NRPyLeakage_units_geom_to_cgs_D, endpoint_thermo.T,
                endpoint_thermo.X_n, endpoint_thermo.X_p, &endpoint_B_n, &endpoint_B_p,
                &endpoint_Y_np, &endpoint_Y_pn, &endpoint_eta_difference),
          ghl_success, "canonical endpoint nucleon blocking", 3078 + (int)endpoint);
    require_condition(
          endpoint_eta_difference == 0.0 && endpoint_Y_np == endpoint_thermo.X_n
                && endpoint_Y_pn == endpoint_thermo.X_p,
          "canonical endpoint transition populations are wrong", 3078 + (int)endpoint);
    ghl_m1_nrpyleakage_raw_rates endpoint_raw;
    require_error(
          ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(
                &endpoint_thermo, eta, &endpoint_raw),
          ghl_success, "endpoint NRPyLeakage raw rates", 3078 + (int)endpoint);
    for(int species = ghl_m1_nrpyleakage_nue; species <= ghl_m1_nrpyleakage_anue;
        ++species) {
      require_condition(
            endpoint_raw.species[species].eta_N_beta_cgs == 0.0
                  && endpoint_raw.species[species].eta_E_beta_mev_cgs == 0.0
                  && endpoint_raw.species[species].kappa_a_N_cc_cgs == 0.0
                  && endpoint_raw.species[species].kappa_a_E_cc_cgs == 0.0,
            "endpoint beta channels are not exactly zero", 3078 + (int)endpoint);
    }
  }

  const ghl_m1_nrpyleakage_raw_rates untouched
        = { .nux_single_species_multiplicity = 17 };
  for(int invalid_species = 0; invalid_species < ghl_m1_nrpyleakage_species_count;
      ++invalid_species) {
    double invalid_eta[ghl_m1_nrpyleakage_species_count] = { eta[0], eta[1], eta[2] };
    invalid_eta[invalid_species] = NAN;
    ghl_m1_nrpyleakage_raw_rates staged = untouched;
    require_error(
          ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(
                &thermo, invalid_eta, &staged),
          ghl_error_m1_microphysics_failure, "direct NRPyLeakage invalid degeneracy",
          3074 + invalid_species);
    require_condition(
          memcmp(&staged, &untouched, sizeof(staged)) == 0,
          "invalid NRPyLeakage degeneracy changed output", 3074 + invalid_species);
  }
  ghl_m1_nrpyleakage_raw_rates staged = untouched;
  ghl_m1_nrpyleakage_thermo_state invalid_thermo = thermo;
  invalid_thermo.T = -1.0;
  staged = untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(&invalid_thermo, eta, &staged),
        ghl_error_m1_microphysics_failure, "direct NRPyLeakage invalid temperature",
        3075);
  require_condition(
        memcmp(&staged, &untouched, sizeof(staged)) == 0,
        "invalid NRPyLeakage temperature changed output", 3075);

  const double invalid_fractions[][2]
        = { { -0.01, 0.20 }, { 0.70, -0.01 }, { 1.01, 0.20 }, { 0.70, 1.01 } };
  for(size_t case_index = 0;
      case_index < sizeof(invalid_fractions) / sizeof(invalid_fractions[0]);
      ++case_index) {
    const ghl_m1_nrpyleakage_thermo_state thermo_before = thermo;
    require_error(
          ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
                1.0e-4, 0.20, 8.0, 12.0, 18.0, 5.0, 20.0,
                invalid_fractions[case_index][0], invalid_fractions[case_index][1],
                &thermo),
          ghl_error_m1_microphysics_failure,
          "NRPyLeakage out-of-range free-nucleon fraction", 3076 + (int)case_index);
    require_condition(
          memcmp(&thermo, &thermo_before, sizeof(thermo)) == 0,
          "invalid free-nucleon fraction changed thermodynamic output",
          3076 + (int)case_index);
  }
}

static void test_nrpyleakage_supported_thermo_boundaries(void) {
  /* The builder is the strict adapter boundary used after the provider owns
   * EOS lookup.  Exercise its input and EOS-quantity validation directly while
   * checking that rejected candidates never reach the caller-owned state. */
  const ghl_m1_nrpyleakage_thermo_state untouched = { .rho = 31.0,
                                                      .T = 33.0,
                                                      .Ye = 0.34,
                                                      .muhat = 35.0,
                                                      .mu_e = 36.0,
                                                      .mu_p = 37.0,
                                                      .mu_n = 38.0,
                                                      .X_n = 0.66,
                                                      .X_p = 0.34 };
  const struct {
    const char *label;
    double rho;
    double Ye;
    double T;
  } invalid_inputs[]
        = { { "strict adapter NaN density", NAN, 0.34, 1.0 },
            { "strict adapter nonpositive density", 0.0, 0.34, 1.0 },
            { "strict adapter NaN temperature", 1.0, 0.34, NAN },
            { "strict adapter nonpositive temperature", 1.0, 0.34, 0.0 },
            { "strict adapter NaN electron fraction", 1.0, NAN, 1.0 },
            { "strict adapter negative electron fraction", 1.0, -DBL_MIN, 1.0 },
            { "strict adapter electron fraction above one", 1.0, 1.0 + DBL_EPSILON,
              1.0 } };
  for(size_t i = 0; i < sizeof(invalid_inputs) / sizeof(invalid_inputs[0]); ++i) {
    ghl_m1_nrpyleakage_thermo_state staged = untouched;
    const ghl_error_codes_t error
          = ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
                invalid_inputs[i].rho, invalid_inputs[i].Ye, invalid_inputs[i].T, 12.0,
                8.0, 5.0, 20.0, 0.66, 0.34, &staged);
    require_error(
          error, ghl_error_m1_microphysics_failure, invalid_inputs[i].label,
          3130 + (int)i);
    require_condition(
          memcmp(&staged, &untouched, sizeof(staged)) == 0,
          "invalid strict adapter input changed thermodynamic output", 3130 + (int)i);
  }

  const struct {
    const char *label;
    double muhat;
    double mu_e;
    double mu_p;
    double mu_n;
    double X_n;
    double X_p;
  } invalid_eos_quantities[] = {
    { "strict adapter NaN muhat", NAN, 8.0, 5.0, 20.0, 0.66, 0.34 },
    { "strict adapter NaN electron chemical potential", 12.0, NAN, 5.0, 20.0, 0.66,
      0.34 },
    { "strict adapter NaN proton chemical potential", 12.0, 8.0, NAN, 20.0, 0.66, 0.34 },
    { "strict adapter NaN neutron chemical potential", 12.0, 8.0, 5.0, NAN, 0.66, 0.34 },
    { "strict adapter NaN neutron fraction", 12.0, 8.0, 5.0, 20.0, NAN, 0.34 },
    { "strict adapter positive-infinite neutron fraction", 12.0, 8.0, 5.0, 20.0,
      INFINITY, 0.34 },
    { "strict adapter negative-infinite neutron fraction", 12.0, 8.0, 5.0, 20.0,
      -INFINITY, 0.34 },
    { "strict adapter NaN proton fraction", 12.0, 8.0, 5.0, 20.0, 0.66, NAN },
    { "strict adapter positive-infinite proton fraction", 12.0, 8.0, 5.0, 20.0, 0.66,
      INFINITY },
    { "strict adapter negative-infinite proton fraction", 12.0, 8.0, 5.0, 20.0, 0.66,
      -INFINITY }
  };
  for(size_t i = 0;
      i < sizeof(invalid_eos_quantities) / sizeof(invalid_eos_quantities[0]); ++i) {
    ghl_m1_nrpyleakage_thermo_state staged = untouched;
    require_error(
          ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
                1.0e-4, 0.34, 8.0, invalid_eos_quantities[i].muhat,
                invalid_eos_quantities[i].mu_e, invalid_eos_quantities[i].mu_p,
                invalid_eos_quantities[i].mu_n, invalid_eos_quantities[i].X_n,
                invalid_eos_quantities[i].X_p, &staged),
          ghl_error_m1_microphysics_failure, invalid_eos_quantities[i].label,
          3140 + (int)i);
    require_condition(
          memcmp(&staged, &untouched, sizeof(staged)) == 0,
          "invalid strict adapter EOS quantity changed thermodynamic output",
          3140 + (int)i);
  }

  const double fraction_roundoff = nrpyleakage_fraction_roundoff_envelope();
  const struct {
    const char *label;
    double X_n;
    double X_p;
    double normalized_X_n;
    double normalized_X_p;
  } accepted_roundoff_compositions[]
        = { { "strict adapter lower composition roundoff", -fraction_roundoff, 0.34, 0.0,
              0.34 },
            { "strict adapter upper composition roundoff", 0.66, 1.0 + fraction_roundoff,
              0.66, 1.0 },
            { "strict adapter lower proton composition roundoff", 0.66,
              -fraction_roundoff, 0.66, 0.0 },
            { "strict adapter upper neutron composition roundoff",
              1.0 + fraction_roundoff, 0.34, 1.0, 0.34 } };
  for(size_t i = 0; i < sizeof(accepted_roundoff_compositions)
                              / sizeof(accepted_roundoff_compositions[0]);
      ++i) {
    ghl_m1_nrpyleakage_thermo_state staged = untouched;
    require_error(
          ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
                1.0e-4, 0.34, 8.0, 12.0, 8.0, 5.0, 20.0,
                accepted_roundoff_compositions[i].X_n,
                accepted_roundoff_compositions[i].X_p, &staged),
          ghl_success, accepted_roundoff_compositions[i].label, 3150 + (int)i);
    require_condition(
          staged.X_n == accepted_roundoff_compositions[i].normalized_X_n
                && staged.X_p == accepted_roundoff_compositions[i].normalized_X_p,
          "accepted composition roundoff was not normalized", 3150 + (int)i);
  }

  const struct {
    const char *label;
    double X_n;
    double X_p;
  } rejected_adjacent_compositions[]
        = { { "strict adapter below lower neutron envelope",
              nextafter(-fraction_roundoff, -INFINITY), 0.34 },
            { "strict adapter above upper neutron envelope",
              nextafter(1.0 + fraction_roundoff, INFINITY), 0.34 },
            { "strict adapter below lower proton envelope", 0.66,
              nextafter(-fraction_roundoff, -INFINITY) },
            { "strict adapter above upper proton envelope", 0.66,
              nextafter(1.0 + fraction_roundoff, INFINITY) } };
  for(size_t i = 0; i < sizeof(rejected_adjacent_compositions)
                              / sizeof(rejected_adjacent_compositions[0]);
      ++i) {
    ghl_m1_nrpyleakage_thermo_state staged = untouched;
    require_error(
          ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
                1.0e-4, 0.34, 8.0, 12.0, 8.0, 5.0, 20.0,
                rejected_adjacent_compositions[i].X_n,
                rejected_adjacent_compositions[i].X_p, &staged),
          ghl_error_m1_microphysics_failure, rejected_adjacent_compositions[i].label,
          3154 + (int)i);
    require_condition(
          memcmp(&staged, &untouched, sizeof(staged)) == 0,
          "adjacent out-of-envelope composition changed thermodynamic output",
          3154 + (int)i);
  }
}

static void test_nrpyleakage_boundary_paths(void) {
  ghl_m1_nrpyleakage_thermo_state thermo = { 0 };
  const double eta[ghl_m1_nrpyleakage_species_count] = { 0.0, 0.0, 0.0 };
  ghl_m1_nrpyleakage_raw_rates raw = { 0 };

  require_error(
        ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
              1.0e-4, 0.20, 8.0, 12.0, 18.0, 5.0, 20.0, 0.70, 0.20, &thermo),
        ghl_success, "strict NRPyLeakage thermo construction", 3090);
  require_condition(
        thermo.X_n == 0.70 && thermo.X_p == 0.20,
        "strict thermo construction changed composition inputs", 3090);

  const double composition_cases[][2]
        = { { 0.80, 12.0 }, { 0.50, 12.0 }, { 0.20, -1000.0 } };
  for(size_t case_index = 0;
      case_index < sizeof(composition_cases) / sizeof(composition_cases[0]);
      ++case_index) {
    const double Ye = composition_cases[case_index][0];
    const double muhat = composition_cases[case_index][1];
    require_error(
          ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
                1.0e-4, Ye, 8.0, muhat, 18.0, 5.0, 20.0, 1.0 - Ye, Ye, &thermo),
          ghl_success, "thermo composition branch", 3091 + (int)case_index);
    require_condition(
          thermo.Ye == Ye && thermo.X_n == 1.0 - Ye && thermo.X_p == Ye,
          "thermo composition branch changed primitive composition",
          3091 + (int)case_index);
  }

  require_error(
        ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
              1.0e-4, 0.20, 8.0, 12.0, 18.0, 5.0, 20.0, 0.70, 0.20, NULL),
        ghl_error_m1_null_pointer, "NULL thermo output", 3094);

  const ghl_m1_nrpyleakage_raw_rates untouched
        = { .nux_single_species_multiplicity = 23 };
  raw = untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(NULL, eta, &raw),
        ghl_error_m1_null_pointer, "NULL raw thermo input", 3095);
  require_condition(
        memcmp(&raw, &untouched, sizeof(raw)) == 0,
        "NULL raw thermo input changed output", 3095);
  raw = untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(&thermo, NULL, &raw),
        ghl_error_m1_null_pointer, "NULL raw degeneracy input", 3096);
  require_condition(
        memcmp(&raw, &untouched, sizeof(raw)) == 0,
        "NULL raw degeneracy input changed output", 3096);
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(&thermo, eta, NULL),
        ghl_error_m1_null_pointer, "NULL raw output", 3097);

  ghl_m1_nrpyleakage_thermo_state invalid_thermo = thermo;
  invalid_thermo.rho = NAN;
  raw = untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(&invalid_thermo, eta, &raw),
        ghl_error_m1_microphysics_failure, "nonfinite raw thermo state", 3104);
  require_condition(
        memcmp(&raw, &untouched, sizeof(raw)) == 0,
        "nonfinite raw thermo state changed output", 3104);

  /* The public raw interface accepts a finite thermodynamic record but still
   * owns the strict nucleon-fraction normalization at the kernel boundary. */
  ghl_m1_nrpyleakage_thermo_state invalid_composition = thermo;
  invalid_composition.X_n = 1.0 + 1.0e-6;
  raw = untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(
              &invalid_composition, eta, &raw),
        ghl_error_m1_microphysics_failure, "out-of-range raw nucleon composition", 3105);
  require_condition(
        memcmp(&raw, &untouched, sizeof(raw)) == 0,
        "out-of-range raw composition changed output", 3105);
}

static void test_masked_kernel_validation(void) {
  const double eta[] = { 0.0, 0.0, 0.0 };
  ghl_m1_nrpyleakage_thermo_state thermo
        = { .rho = 1.e-4, .T = 8.0, .Ye = 0.5, .X_n = 1.0, .X_p = 0.0 };
  ghl_m1_nrpyleakage_raw_rates raw = { .nux_single_species_multiplicity = 29 };
  const ghl_m1_nrpyleakage_raw_rates before = raw;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo_with_mask(NULL, eta, 0, &raw),
        ghl_error_m1_null_pointer, "masked NULL thermo", 3800);
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo_with_mask(
              &thermo, NULL, 0, &raw),
        ghl_error_m1_null_pointer, "masked NULL eta", 3801);
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo_with_mask(
              &thermo, eta, 0, NULL),
        ghl_error_m1_null_pointer, "masked NULL raw", 3802);
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo_with_mask(
              &thermo, eta, -1, &raw),
        ghl_error_m1_microphysics_failure, "unknown raw channels", 3803);
  const double invalid_states[][3] = {
    { 0.0, 8.0, 0.5 }, { 1.e-4, 0.0, 0.5 }, { 1.e-4, 8.0, -0.1 }, { 1.e-4, 8.0, 1.1 }
  };
  for(size_t i = 0; i < sizeof(invalid_states) / sizeof(invalid_states[0]); ++i) {
    ghl_m1_nrpyleakage_thermo_state bad = thermo;
    bad.rho = invalid_states[i][0];
    bad.T = invalid_states[i][1];
    bad.Ye = invalid_states[i][2];
    require_error(
          ghl_m1_nrpyleakage_compute_raw_rates_from_thermo_with_mask(&bad, eta, 0, &raw),
          ghl_error_m1_microphysics_failure, "invalid masked thermodynamics",
          3806 + (int)i);
  }
  /* With beta endpoint populations, density can be tiny while the positive
   * FD equilibrium moments overflow independently of every channel rate. */
  const double large_eta[] = { 1.e44, 1.e15 };
  for(size_t i = 0; i < sizeof(large_eta) / sizeof(large_eta[0]); ++i) {
    ghl_m1_nrpyleakage_thermo_state hot = thermo;
    hot.rho = 1.e-300;
    hot.T = 1.e60;
    const double hot_eta[] = { large_eta[i], 0.0, 0.0 };
    require_error(
          ghl_m1_nrpyleakage_compute_raw_rates_from_thermo_with_mask(
                &hot, hot_eta, 0, &raw),
          ghl_error_m1_microphysics_failure, "overflowing FD equilibrium moment",
          3810 + (int)i);
  }
  /* A finite electron potential can overflow its squared plasmon argument.
   * The pure-neutron endpoint avoids an unrelated beta shift failure. */
  thermo.mu_e = DBL_MAX;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo_with_mask(
              &thermo, eta, ghl_neutrino_rate_channel_plasmon, &raw),
        ghl_error_m1_microphysics_failure, "nonfinite plasmon factor", 3804);
  require_condition(
        memcmp(&raw, &before, sizeof(raw)) == 0, "masked input failure changed output",
        3804);
}

static void test_nrpyleakage_rate_overflow(void) {
  /* The strict raw interface must reject an intermediate rate overflow and
   * leave its caller-owned record untouched. */
  ghl_m1_nrpyleakage_thermo_state thermo;
  require_error(
        ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
              1.0e-4, 0.20, 1.0e100, 12.0, 8.0, 18.0, 5.0, 0.70, 0.20, &thermo),
        ghl_success, "large finite NRPyLeakage thermo construction", 3110);
  const double eta[ghl_m1_nrpyleakage_species_count]
        = { -(thermo.mu_e - thermo.muhat) / thermo.T, 0.0, 0.0 };
  const ghl_m1_nrpyleakage_raw_rates untouched
        = { .nux_single_species_multiplicity = 37 };
  ghl_m1_nrpyleakage_raw_rates strict_raw = untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(&thermo, eta, &strict_raw),
        ghl_error_m1_microphysics_failure, "strict NRPyLeakage rate-overflow rejection",
        3110);
  require_condition(
        memcmp(&strict_raw, &untouched, sizeof(strict_raw)) == 0,
        "strict NRPyLeakage rate overflow changed output", 3110);
}

static void test_nrpyleakage_kernel_representability_edges(void) {
  /* With mu_e/T = 800, the disabled pair-channel F3(-mu_e/T) tail is below
   * binary64 range. The full raw interface still rejects that requested
   * channel set, while a zero mask needs only the finite equilibrium moments
   * and unmasked beta/Kirchhoff diagnostics. */
  ghl_m1_nrpyleakage_thermo_state cold_degenerate_thermo;
  require_error(
        ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
              1.0e-4, 0.50, 0.05, 40.0, 40.0, 0.0, 40.0, 0.50, 0.50,
              &cold_degenerate_thermo),
        ghl_success, "cold degenerate NRPyLeakage thermodynamic state", 3128);
  const double cold_degenerate_eta[ghl_m1_nrpyleakage_species_count] = { 0.0, 0.0, 0.0 };
  const ghl_m1_nrpyleakage_raw_rates cold_degenerate_untouched
        = { .nux_single_species_multiplicity = 47 };
  ghl_m1_nrpyleakage_raw_rates cold_full = cold_degenerate_untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(
              &cold_degenerate_thermo, cold_degenerate_eta, &cold_full),
        ghl_error_m1_microphysics_failure, "enabled cold pair-channel tail", 3128);
  require_condition(
        memcmp(&cold_full, &cold_degenerate_untouched, sizeof(cold_full)) == 0,
        "enabled cold pair-channel failure changed raw output", 3128);

  ghl_m1_nrpyleakage_raw_rates cold_masked = cold_degenerate_untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo_with_mask(
              &cold_degenerate_thermo, cold_degenerate_eta, 0, &cold_masked),
        ghl_success, "zero-mask cold degenerate raw rates", 3128);
  require_condition(
        cold_masked.nux_single_species_multiplicity == 1,
        "zero-mask cold raw multiplicity changed", 3128);
  for(int species = 0; species < ghl_m1_nrpyleakage_species_count; ++species) {
    const ghl_m1_nrpyleakage_species_raw_rates *const rate
          = &cold_masked.species[species];
    require_condition(
          isfinite(rate->n_eq_cgs) && rate->n_eq_cgs > 0.0
                && isfinite(rate->J_eq_mev_cgs) && rate->J_eq_mev_cgs > 0.0
                && isfinite(rate->mean_energy_mev) && rate->mean_energy_mev > 0.0,
          "zero-mask cold equilibrium moments are invalid", 3128);
    require_condition(
          rate->eta_N_pair_cgs == 0.0 && rate->eta_E_pair_mev_cgs == 0.0
                && rate->eta_N_plasmon_cgs == 0.0 && rate->eta_E_plasmon_mev_cgs == 0.0
                && rate->eta_N_brems_cgs == 0.0 && rate->eta_E_brems_mev_cgs == 0.0
                && rate->kappa_s_N_neutron_cgs == 0.0
                && rate->kappa_s_N_proton_cgs == 0.0
                && rate->kappa_s_E_neutron_cgs == 0.0
                && rate->kappa_s_E_proton_cgs == 0.0,
          "zero-mask cold raw channel fields are nonzero", 3128);
  }

  /* The public state remains finite, but the cgs density conversion can
   * overflow before the canonical blocking helper is entered. */
  ghl_m1_nrpyleakage_thermo_state blocking_failure_thermo;
  require_error(
        ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
              DBL_MAX, 0.50, 1.0, 0.0, 0.0, 0.0, 0.0, 0.50, 0.50,
              &blocking_failure_thermo),
        ghl_success, "finite thermo with overflowing blocking density", 3127);
  const double blocking_failure_eta[ghl_m1_nrpyleakage_species_count]
        = { 0.0, 0.0, 0.0 };
  const ghl_m1_nrpyleakage_raw_rates blocking_untouched
        = { .nux_single_species_multiplicity = 43 };
  ghl_m1_nrpyleakage_raw_rates blocking_staged = blocking_untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(
              &blocking_failure_thermo, blocking_failure_eta, &blocking_staged),
        ghl_error_m1_microphysics_failure,
        "strict rejection of overflowing blocking density", 3127);
  require_condition(
        memcmp(&blocking_staged, &blocking_untouched, sizeof(blocking_staged)) == 0,
        "overflowing blocking density changed raw output", 3127);

  /* A positive FD moment can be representable while its product with T
   * underflows. The strict raw mean energy must stay positive, so reject this
   * finite-input endpoint without publishing the partially computed rates. */
  ghl_m1_nrpyleakage_thermo_state cold_thermo;
  require_error(
        ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
              1.0e-4, 0.50, 1.0e-30, 0.0, 0.0, 0.0, 0.0, 0.50, 0.50, &cold_thermo),
        ghl_success, "cold thermo with nonrepresentable mean energy", 3124);
  const double cold_eta[ghl_m1_nrpyleakage_species_count] = { -700.0, 0.0, 0.0 };
  const ghl_m1_nrpyleakage_raw_rates cold_untouched
        = { .nux_single_species_multiplicity = 41 };
  ghl_m1_nrpyleakage_raw_rates cold_staged = cold_untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(
              &cold_thermo, cold_eta, &cold_staged),
        ghl_error_m1_microphysics_failure,
        "strict rejection of underflowing mean energy", 3124);
  require_condition(
        memcmp(&cold_staged, &cold_untouched, sizeof(cold_staged)) == 0,
        "underflowing mean energy changed raw output", 3124);

  /* Keep the first FD group and all Fermi factors representable, then make a
   * later species moment underflow.  The adapter must collect that later
   * helper failure and preserve the caller-owned record. */
  ghl_m1_nrpyleakage_thermo_state late_fd_thermo;
  require_error(
        ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
              1.0e-4, 0.50, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, &late_fd_thermo),
        ghl_success, "thermo for later FD failure", 3126);
  const double late_fd_eta[ghl_m1_nrpyleakage_species_count] = { -1000.0, 0.0, 0.0 };
  ghl_m1_nrpyleakage_raw_rates late_fd_staged = cold_untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(
              &late_fd_thermo, late_fd_eta, &late_fd_staged),
        ghl_error_m1_microphysics_failure,
        "strict rejection of later nonpositive FD moment", 3126);
  require_condition(
        memcmp(&late_fd_staged, &cold_untouched, sizeof(late_fd_staged)) == 0,
        "later nonpositive FD moment changed raw output", 3126);

  /* The pair channel is disabled, leaving the plasmon Fermi factor as the
   * first unrepresentable helper result for this finite degeneracy. */
  const double plasmon_tail_eta[ghl_m1_nrpyleakage_species_count] = { 800.0, 0.0, 0.0 };
  ghl_m1_nrpyleakage_raw_rates plasmon_staged = cold_untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo_with_mask(
              &late_fd_thermo, plasmon_tail_eta, ghl_neutrino_rate_channel_plasmon,
              &plasmon_staged),
        ghl_error_m1_microphysics_failure, "plasmon Fermi tail underflow", 3129);
  require_condition(
        memcmp(&plasmon_staged, &cold_untouched, sizeof(plasmon_staged)) == 0,
        "plasmon helper failure changed raw output", 3129);

  /* These are finite caller-owned thermo states.  The first case makes the
   * generated FD polynomial overflow; the strict wrapper must report that
   * loss of a positive FD moment without publishing a partial record. */
  ghl_m1_nrpyleakage_thermo_state thermo;
  require_error(
        ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
              1.0e-4, 0.50, 1.0, 0.0, 0.0, 0.0, 0.0, 0.50, 0.50, &thermo),
        ghl_success, "finite thermo with overflowing FD moment", 3120);
  thermo.mu_e = 1.0e100;
  const double eta[ghl_m1_nrpyleakage_species_count] = { 0.0, 0.0, 0.0 };
  const ghl_m1_nrpyleakage_raw_rates untouched
        = { .nux_single_species_multiplicity = 41 };
  ghl_m1_nrpyleakage_raw_rates staged = untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(&thermo, eta, &staged),
        ghl_error_m1_microphysics_failure, "strict rejection of overflowing FD moment",
        3120);
  require_condition(
        memcmp(&staged, &untouched, sizeof(staged)) == 0,
        "overflowing FD moment changed raw output", 3120);

  /* The low-z approximation is a finite caller path too.  At z=-1000 its
   * positive FD approximation underflows to zero, so the strict evaluator
   * takes its nonpositive-integral guard before any raw field is published. */
  thermo.mu_e = -1000.0;
  staged = untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(&thermo, eta, &staged),
        ghl_error_m1_microphysics_failure, "strict rejection of nonpositive FD moment",
        3120);
  require_condition(
        memcmp(&staged, &untouched, sizeof(staged)) == 0,
        "nonpositive FD moment changed raw output", 3120);

  /* A finite density/temperature pair can still make mu_e/T nonrepresentable
   * at the first FD call.  This exercises the finite-input wrapper's
   * nonfinite FD-argument guard rather than the caller-side state validation. */
  require_error(
        ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
              1.0e-4, 0.50, 1.0e-4, 0.0, DBL_MAX, 0.0, 0.0, 0.50, 0.50, &thermo),
        ghl_success, "finite thermo with nonrepresentable FD argument", 3121);
  staged = untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(&thermo, eta, &staged),
        ghl_error_m1_microphysics_failure, "strict rejection of nonfinite FD argument",
        3121);
  require_condition(
        memcmp(&staged, &untouched, sizeof(staged)) == 0,
        "nonfinite FD argument changed raw output", 3121);

  /* eta_nue enters the first blocking factor directly.  At 1000 the
   * exp(-x) tail rounds to zero, so the strict kernel has to reject that
   * unsupported nonzero factor instead of silently changing the rate. */
  require_error(
        ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
              1.0e-4, 0.20, 8.0, 12.0, 18.0, 5.0, 20.0, 0.70, 0.20, &thermo),
        ghl_success, "finite thermo for FD-factor tail", 3122);
  const double extreme_eta[ghl_m1_nrpyleakage_species_count] = { 1000.0, 0.0, 0.0 };
  staged = untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(&thermo, extreme_eta, &staged),
        ghl_error_m1_microphysics_failure, "strict rejection of underflowing FD factor",
        3122);
  require_condition(
        memcmp(&staged, &untouched, sizeof(staged)) == 0,
        "underflowing FD factor changed raw output", 3122);

  /* Exercise the representable positive-x branch of the Fermi blocking
   * factor.  The extreme case above deliberately stops before assigning the
   * factor; this moderate degeneracy must reach the finite assignment and
   * publish a positive pair rate. */
  const double positive_factor_eta[ghl_m1_nrpyleakage_species_count]
        = { 10.0, 0.0, 0.0 };
  ghl_m1_nrpyleakage_raw_rates positive_factor_rates = untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(
              &thermo, positive_factor_eta, &positive_factor_rates),
        ghl_success, "strict acceptance of representable positive FD factor", 3122);
  require_condition(
        isfinite(positive_factor_rates.species[ghl_m1_nrpyleakage_nue].eta_N_pair_cgs)
              && positive_factor_rates.species[ghl_m1_nrpyleakage_nue].eta_N_pair_cgs
                       > 0.0,
        "representable positive FD factor did not produce a pair rate", 3122);

  /* The equilibrium moments are products of positive FD moments and powers
   * of T.  At this finite positive temperature, the generated pair-energy
   * ratio is 0/0 before the raw moments can be published; preserve the caller
   * record while that strict representability failure is reported. */
  require_error(
        ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
              1.0e-4, 0.50, 1.0e-100, 0.0, 0.0, 0.0, 0.0, 0.50, 0.50, &thermo),
        ghl_success, "finite thermo at pair-ratio underflow boundary", 3123);
  staged = untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(&thermo, eta, &staged),
        ghl_error_m1_microphysics_failure,
        "strict rejection of nonrepresentable pair ratio", 3123);
  require_condition(
        memcmp(&staged, &untouched, sizeof(staged)) == 0,
        "nonrepresentable pair ratio changed raw output", 3123);

  /* A finite EOS state can still make the canonical reaction-energy shift
   * nonfinite: use the smallest representable neutron fraction to produce a
   * large negative kinetic degeneracy difference, then add it to DBL_MAX.
   * The adapter must translate the beta-helper error and keep its output
   * transactionally untouched. */
  ghl_m1_nrpyleakage_thermo_state shift_failure_thermo;
  require_error(
        ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
              1.0e-4, 0.50, 9.0e304, DBL_MAX, 0.0, 0.0, 0.0, DBL_MIN, 1.0,
              &shift_failure_thermo),
        ghl_success, "finite thermo with overflowing reaction shift", 3125);
  const double shift_failure_eta[ghl_m1_nrpyleakage_species_count] = { 0.0, 0.0, 0.0 };
  ghl_m1_nrpyleakage_raw_rates shift_failure_staged = untouched;
  require_error(
        ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(
              &shift_failure_thermo, shift_failure_eta, &shift_failure_staged),
        ghl_error_m1_microphysics_failure,
        "strict rejection of nonfinite reaction shift", 3125);
  require_condition(
        memcmp(&shift_failure_staged, &untouched, sizeof(shift_failure_staged)) == 0,
        "nonfinite reaction shift changed raw output", 3125);
}

static void
make_primitives(m1_test_rng *restrict rng, ghl_primitive_quantities *restrict prims) {
  *prims = (ghl_primitive_quantities){ 0 };
  prims->rho = exp(m1_test_rng_between(rng, log(0.05), log(20.0)));
  prims->temperature = exp(m1_test_rng_between(rng, log(0.05), log(20.0)));
  prims->Y_e = m1_test_rng_between(rng, 0.02, 0.98);
  prims->eps = prims->temperature;
  prims->press = prims->rho * prims->eps * 0.1;
  prims->entropy = m1_test_rng_between(rng, 0.05, 2.0);
}

#define ghl_neutrino_rate_provider_context m1_test_reference_provider_context
#define ghl_neutrino_rate_provider_cache   m1_test_reference_provider_cache
#define ghl_neutrino_rate_provider_initialize_reference \
  m1_test_reference_provider_initialize
#define ghl_neutrino_rate_provider_cache_initialize \
  m1_test_reference_provider_cache_initialize
#define ghl_neutrino_rate_provider_compute_cell m1_test_reference_provider_compute_cell

static double reference_safe_exp(const double x) {
  const double bounded = fmin(40.0, fmax(-40.0, x));
  return exp(bounded);
}

/* Independent oracle for the table-free provider. Keep this formula-based
 * check separate from the implementation so channel masks cannot silently
 * produce merely valid, but physically empty, bundles. */
static void compute_expected_reference_rates(
      const ghl_neutrino_rate_provider_context *restrict provider,
      const ghl_primitive_quantities *restrict prims,
      ghl_m1_neutrino_rates expected[ghl_m1_neutrino_species_count]) {
  const double rho = prims->rho;
  const double T = prims->temperature;
  const double Ye = prims->Y_e;
  const double X_p = Ye;
  const double X_n = 1.0 - Ye;
  const double mu_e = T * log(Ye / (1.0 - Ye));
  const double mu_p = T * log(fmax(X_p, 1.0e-12));
  const double mu_n = T * log(fmax(X_n, 1.0e-12));
  const double mu_nue = mu_e + mu_p - mu_n;
  const double theta
        = fmax(T * provider->temperature_code_to_mev, provider->min_mean_energy);
  const double thermal_number = rho / fmax(provider->baryon_mass_code, DBL_MIN);
  const int mask = provider->channel_mask;

  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    const double lepton_weight = species == ghl_m1_neutrino_nue    ? 1.0
                                 : species == ghl_m1_neutrino_anue ? -1.0
                                                                   : 0.0;
    const double mu_nu = species == ghl_m1_neutrino_nue    ? mu_nue
                         : species == ghl_m1_neutrino_anue ? -mu_nue
                                                           : 0.0;
    const double degeneracy = mu_nu / theta;
    const double exp_degeneracy = reference_safe_exp(degeneracy);
    const double occupancy = exp_degeneracy / (1.0 + exp_degeneracy);
    const double multiplicity
          = species == ghl_m1_neutrino_nux ? provider->nu_x_multiplicity : 1.0;
    const double mean_energy
          = fmax(provider->min_mean_energy, theta * (3.1514 + 0.25 * fabs(degeneracy)));
    double n_eq = multiplicity * thermal_number * 1.0e-3 * occupancy;
    double kappa_a_N = 0.0;
    double kappa_a_E = 0.0;
    double kappa_a_N_cc = 0.0;
    double kappa_s = 0.0;

    expected[species] = (ghl_m1_neutrino_rates){ 0 };
    expected[species].species = (ghl_m1_neutrino_species_t)species;
    expected[species].mean_energy = mean_energy;
    expected[species].lepton_weight = lepton_weight;

    if((mask & ghl_neutrino_rate_channel_charged_current) != 0) {
      if(species == ghl_m1_neutrino_nue) {
        kappa_a_N_cc = provider->charged_current_scale * rho * X_n;
      }
      else if(species == ghl_m1_neutrino_anue) {
        kappa_a_N_cc = provider->charged_current_scale * rho * X_p;
      }
      kappa_a_N += kappa_a_N_cc;
      kappa_a_E += kappa_a_N_cc * mean_energy * mean_energy;
    }
    if((mask & ghl_neutrino_rate_channel_nucleon_scattering) != 0) {
      kappa_s = provider->scattering_scale * rho * fmax(X_n + X_p, 0.0) * mean_energy
                * mean_energy;
    }
    if((mask & ghl_neutrino_rate_channel_pair) != 0) {
      if(species == ghl_m1_neutrino_nux) {
        const double pair_kappa = provider->pair_scale * rho * theta * theta;
        kappa_a_N += pair_kappa;
        kappa_a_E += pair_kappa * mean_energy;
      }
      n_eq += multiplicity * thermal_number * 1.0e-4;
    }
    if((mask & ghl_neutrino_rate_channel_bremsstrahlung) != 0
       && species == ghl_m1_neutrino_nux) {
      const double brems_kappa = provider->bremsstrahlung_scale * rho * rho;
      kappa_a_N += brems_kappa;
      kappa_a_E += brems_kappa * mean_energy;
    }
    if((mask & ghl_neutrino_rate_channel_plasmon) != 0
       && species == ghl_m1_neutrino_nux) {
      const double plasmon_kappa = provider->plasmon_scale * theta * theta * theta;
      kappa_a_N += plasmon_kappa;
      kappa_a_E += plasmon_kappa * mean_energy;
    }

    expected[species].kappa_a_N = fmax(0.0, kappa_a_N);
    expected[species].kappa_a_E = fmax(0.0, kappa_a_E);
    expected[species].kappa_s = fmax(0.0, kappa_s);
    expected[species].kappa_tr = expected[species].kappa_a_E + expected[species].kappa_s;
    expected[species].n_eq = fmax(0.0, n_eq);
    expected[species].J_eq = mean_energy * expected[species].n_eq;
    expected[species].eta_N = expected[species].kappa_a_N * expected[species].n_eq;
    expected[species].eta_E = expected[species].kappa_a_E * expected[species].J_eq;
    expected[species].kappa_a_N_cc = lepton_weight == 0.0 ? 0.0 : kappa_a_N_cc;
    expected[species].eta_N_cc = expected[species].kappa_a_N_cc * expected[species].n_eq;
  }

  const int process_masks[]
        = { ghl_neutrino_rate_channel_pair, ghl_neutrino_rate_channel_plasmon,
            ghl_neutrino_rate_channel_bremsstrahlung };
  const double process_scales[] = { provider->pair_scale * rho * theta * theta,
                                    provider->plasmon_scale * theta * theta * theta,
                                    provider->bremsstrahlung_scale * rho * rho };
  const double electron_n_eq_geometric_mean
        = sqrt(expected[ghl_m1_neutrino_nue].n_eq)
          * sqrt(expected[ghl_m1_neutrino_anue].n_eq);
  for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
    if((mask & process_masks[process]) == 0) {
      continue;
    }
    const double kappa_N = process_scales[process];
    for(int species = ghl_m1_neutrino_nue; species <= ghl_m1_neutrino_anue; ++species) {
      expected[species].eta_N_pair[process] = kappa_N * electron_n_eq_geometric_mean;
      expected[species].eta_E_pair[process]
            = kappa_N * expected[species].mean_energy * expected[species].J_eq;
    }
  }
}

static void compare_reference_rates(
      const ghl_neutrino_rate_provider_context *restrict provider,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates actual[ghl_m1_neutrino_species_count],
      const int case_index) {
  ghl_m1_neutrino_rates expected[ghl_m1_neutrino_species_count];
  compute_expected_reference_rates(provider, prims, expected);
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    const double actual_scalars[]
          = { actual[species].eta_N,       actual[species].eta_E,
              actual[species].kappa_a_N,   actual[species].kappa_a_E,
              actual[species].kappa_s,     actual[species].kappa_tr,
              actual[species].n_eq,        actual[species].J_eq,
              actual[species].mean_energy, actual[species].lepton_weight,
              actual[species].eta_N_cc,    actual[species].kappa_a_N_cc };
    const double expected_scalars[]
          = { expected[species].eta_N,       expected[species].eta_E,
              expected[species].kappa_a_N,   expected[species].kappa_a_E,
              expected[species].kappa_s,     expected[species].kappa_tr,
              expected[species].n_eq,        expected[species].J_eq,
              expected[species].mean_energy, expected[species].lepton_weight,
              expected[species].eta_N_cc,    expected[species].kappa_a_N_cc };
    for(size_t i = 0; i < sizeof(actual_scalars) / sizeof(actual_scalars[0]); ++i) {
      require_condition(
            fabs(actual_scalars[i] - expected_scalars[i])
                  <= 2.0e-11 * fmax(1.0, fabs(expected_scalars[i])),
            "reference provider formula mismatch", case_index);
    }
    for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
      require_condition(
            fabs(actual[species].eta_N_pair[process]
                 - expected[species].eta_N_pair[process])
                        <= 2.0e-11
                                 * fmax(1.0, fabs(expected[species].eta_N_pair[process]))
                  && fabs(actual[species].eta_E_pair[process]
                          - expected[species].eta_E_pair[process])
                           <= 2.0e-11
                                    * fmax(
                                          1.0,
                                          fabs(expected[species].eta_E_pair[process])),
            "reference provider process mismatch", case_index);
    }
  }
}

static void
initialize_sentinel_rates(ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count]) {
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    rates[species] = (ghl_m1_neutrino_rates){ 0 };
    rates[species].species = (ghl_m1_neutrino_species_t)species;
    rates[species].eta_N = -11.0 - species;
    rates[species].eta_E = -21.0 - species;
    rates[species].kappa_a_N = -31.0 - species;
    rates[species].kappa_a_E = -41.0 - species;
    rates[species].kappa_s = -51.0 - species;
    rates[species].kappa_tr = -61.0 - species;
    rates[species].n_eq = -71.0 - species;
    rates[species].J_eq = -81.0 - species;
    rates[species].mean_energy = -91.0 - species;
    rates[species].lepton_weight = -101.0 - species;
    rates[species].eta_N_cc = -111.0 - species;
    rates[species].kappa_a_N_cc = -121.0 - species;
    for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
      rates[species].eta_N_pair[process] = -131.0 - process - species;
      rates[species].eta_E_pair[process] = -141.0 - process - species;
    }
  }
}

static void
test_invalid_provider_contexts(const ghl_primitive_quantities *restrict prims) {
  static const char *const labels[] = { "unknown channel bit",
                                        "failure policy below enum range",
                                        "failure policy above enum range",
                                        "table policy below enum range",
                                        "table policy above enum range",
                                        "nonpositive temperature conversion",
                                        "nonpositive baryon mass",
                                        "negative charged-current scale",
                                        "negative scattering scale",
                                        "negative pair scale",
                                        "negative bremsstrahlung scale",
                                        "negative plasmon scale",
                                        "nonpositive minimum mean energy",
                                        "nonpositive equilibrium recovery rate" };
  const int variant_count = (int)(sizeof(labels) / sizeof(labels[0]));

  for(int variant = 0; variant < variant_count; ++variant) {
    ghl_neutrino_rate_provider_context provider;
    require_error(
          ghl_neutrino_rate_provider_initialize_reference(&provider), ghl_success,
          "invalid-context baseline initialization", 4000 + variant);
    switch(variant) {
      case 0:
        provider.channel_mask = 1 << 8;
        break;
      case 1:
        provider.failure_policy = (ghl_neutrino_rate_failure_policy_t)-1;
        break;
      case 2:
        provider.failure_policy = (ghl_neutrino_rate_failure_policy_t)4;
        break;
      case 3:
        provider.table_bounds_policy = (ghl_neutrino_rate_table_bounds_policy_t)-1;
        break;
      case 4:
        provider.table_bounds_policy = (ghl_neutrino_rate_table_bounds_policy_t)2;
        break;
      case 5:
        provider.temperature_code_to_mev = 0.0;
        break;
      case 6:
        provider.baryon_mass_code = 0.0;
        break;
      case 7:
        provider.charged_current_scale = -1.0;
        break;
      case 8:
        provider.scattering_scale = -1.0;
        break;
      case 9:
        provider.pair_scale = -1.0;
        break;
      case 10:
        provider.bremsstrahlung_scale = -1.0;
        break;
      case 11:
        provider.plasmon_scale = -1.0;
        break;
      case 12:
        provider.min_mean_energy = 0.0;
        break;
      case 13:
        provider.equilibrium_recovery_rate = 0.0;
        break;
      default:
        provider_test_error(
              "M1 rate-provider invalid-context test has an unknown variant");
    }

    ghl_neutrino_rate_provider_cache cache;
    ghl_neutrino_rate_provider_cache_initialize(&cache);
    const ghl_neutrino_rate_provider_cache cache_before = cache;
    ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
    initialize_sentinel_rates(rates);
    const ghl_m1_neutrino_rates rates_before[ghl_m1_neutrino_species_count]
          = { rates[0], rates[1], rates[2] };
    ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
    const ghl_error_codes_t error = ghl_neutrino_rate_provider_compute_cell(
          &provider, &cache, &diagnostics, NULL, prims, rates);
    require_error(
          error, ghl_error_m1_microphysics_failure, labels[variant], 4000 + variant);
    require_condition(
          memcmp(&cache, &cache_before, sizeof(cache)) == 0
                && same_rate_bundle(rates, rates_before),
          "invalid provider context was not transactional", 4000 + variant);
    require_condition(
          diagnostics.failures == 1 && diagnostics.last_error == error,
          "invalid provider context diagnostics are incomplete", 4000 + variant);
  }
}

static void
test_nonfinite_provider_contexts(const ghl_primitive_quantities *restrict prims) {
  static const char *const labels[] = { "nonfinite multiplicity",
                                        "nonfinite charged-current scale",
                                        "nonfinite scattering scale",
                                        "nonfinite pair scale",
                                        "nonfinite bremsstrahlung scale",
                                        "nonfinite plasmon scale",
                                        "nonfinite temperature conversion",
                                        "nonfinite baryon mass",
                                        "nonfinite equilibrium recovery rate",
                                        "nonfinite minimum mean energy" };

  for(int variant = 0; variant < (int)(sizeof(labels) / sizeof(labels[0])); ++variant) {
    ghl_neutrino_rate_provider_context provider;
    require_error(
          ghl_neutrino_rate_provider_initialize_reference(&provider), ghl_success,
          "nonfinite-context baseline initialization", 4050 + variant);
    switch(variant) {
      case 0:
        provider.nu_x_multiplicity = NAN;
        break;
      case 1:
        provider.charged_current_scale = NAN;
        break;
      case 2:
        provider.scattering_scale = NAN;
        break;
      case 3:
        provider.pair_scale = NAN;
        break;
      case 4:
        provider.bremsstrahlung_scale = NAN;
        break;
      case 5:
        provider.plasmon_scale = NAN;
        break;
      case 6:
        provider.temperature_code_to_mev = NAN;
        break;
      case 7:
        provider.baryon_mass_code = NAN;
        break;
      case 8:
        provider.equilibrium_recovery_rate = NAN;
        break;
      case 9:
        provider.min_mean_energy = NAN;
        break;
      default:
        provider_test_error(
              "M1 rate-provider nonfinite-context test has an unknown variant");
    }

    ghl_neutrino_rate_provider_cache cache;
    ghl_neutrino_rate_provider_cache_initialize(&cache);
    const ghl_neutrino_rate_provider_cache cache_before = cache;
    ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
    initialize_sentinel_rates(rates);
    const ghl_m1_neutrino_rates rates_before[ghl_m1_neutrino_species_count]
          = { rates[0], rates[1], rates[2] };
    ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
    const ghl_error_codes_t error = ghl_neutrino_rate_provider_compute_cell(
          &provider, &cache, &diagnostics, NULL, prims, rates);
    require_error(
          error, ghl_error_m1_microphysics_failure, labels[variant], 4050 + variant);
    require_condition(
          memcmp(&cache, &cache_before, sizeof(cache)) == 0
                && same_rate_bundle(rates, rates_before),
          "nonfinite provider context was not transactional", 4050 + variant);
    require_condition(
          diagnostics.failures == 1 && diagnostics.last_error == error,
          "nonfinite provider context diagnostics are incomplete", 4050 + variant);
  }
}

static void test_invalid_primitive_keys(void) {
  ghl_neutrino_rate_provider_context provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_reference(&provider), ghl_success,
        "invalid-key provider initialization", 4100);
  ghl_neutrino_rate_provider_cache cache;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
  ghl_primitive_quantities prims = { 0 };
  prims.rho = 1.0;
  prims.temperature = 1.0;
  prims.Y_e = 0.5;
  prims.eps = 1.0;
  ghl_m1_neutrino_rates valid_rates[ghl_m1_neutrino_species_count];
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, NULL, &prims, valid_rates),
        ghl_success, "invalid-key cache seed", 4100);

  ghl_primitive_quantities invalid[] = { prims, prims, prims, prims, prims, prims };
  invalid[0].rho = 0.0;
  invalid[1].temperature = -1.0;
  invalid[1].eps = 0.0;
  invalid[2].Y_e = -0.01;
  invalid[3].Y_e = 1.01;
  invalid[4].rho = NAN;
  invalid[5].Y_e = NAN;
  for(int case_index = 0; case_index < 6; ++case_index) {
    ghl_neutrino_rate_provider_cache cache_before = cache;
    ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
    memcpy(rates, valid_rates, sizeof(rates));
    const ghl_m1_neutrino_rates rates_before[ghl_m1_neutrino_species_count]
          = { rates[0], rates[1], rates[2] };
    const ghl_error_codes_t error = ghl_neutrino_rate_provider_compute_cell(
          &provider, &cache, &diagnostics, NULL, &invalid[case_index], rates);
    require_error(
          error, ghl_error_m1_microphysics_failure, "invalid primitive key",
          4101 + case_index);
    require_condition(
          memcmp(&cache, &cache_before, sizeof(cache)) == 0
                && same_rate_bundle(rates, rates_before),
          "invalid primitive key was not transactional", 4101 + case_index);
  }
}

static void test_temperature_recovery(void) {
  ghl_neutrino_rate_provider_context provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_reference(&provider), ghl_success,
        "temperature-recovery provider initialization", 4200);
  ghl_primitive_quantities prims = { 0 };
  prims.rho = 0.8;
  prims.temperature = NAN;
  prims.eps = 0.65;
  prims.Y_e = 0.41;
  ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
  const ghl_error_codes_t error = ghl_neutrino_rate_provider_compute_cell(
        &provider, NULL, NULL, NULL, &prims, rates);
  require_error(error, ghl_success, "temperature recovered from epsilon", 4200);
  validate_rate_bundle(rates, 4200);
  ghl_primitive_quantities expected_prims = prims;
  expected_prims.temperature = prims.eps;
  compare_reference_rates(&provider, &expected_prims, rates, 4200);

  const double nonpositive_temperatures[] = { 0.0, -1.0 };
  for(size_t i = 0;
      i < sizeof(nonpositive_temperatures) / sizeof(nonpositive_temperatures[0]); ++i) {
    prims.temperature = nonpositive_temperatures[i];
    prims.eps = 0.65;
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, NULL, NULL, NULL, &prims, rates),
          ghl_success, "nonpositive temperature recovered from epsilon", 4201 + (int)i);
    validate_rate_bundle(rates, 4201 + (int)i);
    expected_prims = prims;
    expected_prims.temperature = prims.eps;
    compare_reference_rates(&provider, &expected_prims, rates, 4201 + (int)i);
  }

  const double invalid_eps[] = { 0.0, NAN, -1.0 };
  for(size_t i = 0; i < sizeof(invalid_eps) / sizeof(invalid_eps[0]); ++i) {
    prims.temperature = NAN;
    prims.eps = invalid_eps[i];
    initialize_sentinel_rates(rates);
    const ghl_m1_neutrino_rates rates_before[ghl_m1_neutrino_species_count]
          = { rates[0], rates[1], rates[2] };
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, NULL, NULL, NULL, &prims, rates),
          ghl_error_m1_microphysics_failure,
          "invalid epsilon rejected during temperature recovery", 4203 + (int)i);
    require_condition(
          same_rate_bundle(rates, rates_before),
          "invalid epsilon changed rates during temperature recovery", 4203 + (int)i);
  }
}

static void test_default_provider(m1_test_rng *restrict rng) {
  ghl_neutrino_rate_provider_context provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_reference(&provider), ghl_success,
        "default provider initialization", 0);
  require_error(
        ghl_neutrino_rate_provider_initialize_reference(NULL), ghl_error_m1_null_pointer,
        "NULL default provider initialization", 0);

  ghl_neutrino_rate_provider_cache cache;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  ghl_neutrino_rate_provider_cache_initialize(NULL);
  require_condition(
        !cache.thermo_valid && !cache.rates_valid,
        "cache initializer published a valid record", 0);

  static const int channel_masks[]
        = { 0,
            ghl_neutrino_rate_channel_charged_current
                  | ghl_neutrino_rate_channel_nucleon_scattering,
            ghl_neutrino_rate_channel_pair,
            ghl_neutrino_rate_channel_plasmon,
            ghl_neutrino_rate_channel_bremsstrahlung,
            ghl_neutrino_rate_channel_charged_current
                  | ghl_neutrino_rate_channel_nucleon_scattering
                  | ghl_neutrino_rate_channel_pair | ghl_neutrino_rate_channel_plasmon
                  | ghl_neutrino_rate_channel_bremsstrahlung };
  ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };

  for(int case_index = 0; case_index < PROVIDER_RANDOM_CASES; ++case_index) {
    ghl_primitive_quantities prims;
    make_primitives(rng, &prims);
    provider.channel_mask = channel_masks
          [case_index % (int)(sizeof(channel_masks) / sizeof(channel_masks[0]))];
    provider.eos_generation = (uint64_t)case_index;
    provider.charged_current_scale = m1_test_rng_between(rng, 1.0e-4, 2.0e-2);
    provider.scattering_scale = m1_test_rng_between(rng, 1.0e-4, 2.0e-2);
    provider.pair_scale = m1_test_rng_between(rng, 1.0e-6, 2.0e-4);
    provider.bremsstrahlung_scale = m1_test_rng_between(rng, 1.0e-6, 2.0e-4);
    provider.plasmon_scale = m1_test_rng_between(rng, 1.0e-7, 2.0e-5);

    ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
    const ghl_error_codes_t error = ghl_neutrino_rate_provider_compute_cell(
          &provider, &cache, &diagnostics, NULL, &prims, rates);
    require_error(error, ghl_success, "random reference provider call", case_index);
    validate_rate_bundle(rates, case_index);
    compare_reference_rates(&provider, &prims, rates, case_index);
    require_condition(
          diagnostics.active_channel_mask == provider.channel_mask,
          "diagnostics lost the active channel mask", case_index);
  }

  /* Repeated identical keys must hit the cache and publish exactly the same
   * frozen bundle; a context generation change must invalidate that hit. */
  require_error(
        ghl_neutrino_rate_provider_initialize_reference(&provider), ghl_success,
        "cache test provider initialization", 1000);
  provider.channel_mask
        = ghl_neutrino_rate_channel_charged_current
          | ghl_neutrino_rate_channel_nucleon_scattering | ghl_neutrino_rate_channel_pair
          | ghl_neutrino_rate_channel_plasmon | ghl_neutrino_rate_channel_bremsstrahlung;
  provider.eos_generation = 17;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  diagnostics = (ghl_neutrino_rate_provider_diagnostics){ 0 };
  ghl_primitive_quantities prims = { 0 };
  prims.rho = 1.25;
  prims.temperature = 0.8;
  prims.Y_e = 0.37;
  prims.eps = prims.temperature;
  ghl_m1_neutrino_rates first[ghl_m1_neutrino_species_count];
  ghl_error_codes_t error = ghl_neutrino_rate_provider_compute_cell(
        &provider, &cache, &diagnostics, NULL, &prims, first);
  require_error(error, ghl_success, "initial cache provider call", 1000);
  require_condition(
        diagnostics.cache_misses == 1 && diagnostics.cache_hits == 0,
        "initial provider call did not record a cache miss", 1000);
  validate_rate_bundle(first, 1000);

  ghl_m1_neutrino_rates second[ghl_m1_neutrino_species_count];
  error = ghl_neutrino_rate_provider_compute_cell(
        &provider, &cache, &diagnostics, NULL, &prims, second);
  require_error(error, ghl_success, "repeated cache provider call", 1001);
  require_condition(
        diagnostics.cache_hits == 1, "repeated provider call did not hit the cache",
        1001);
  require_condition(
        same_rate_bundle(first, second), "cache hit changed the frozen rate bundle",
        1001);

  prims.rho *= 1.01;
  error = ghl_neutrino_rate_provider_compute_cell(
        &provider, &cache, &diagnostics, NULL, &prims, second);
  require_error(error, ghl_success, "changed-key provider call", 1002);
  require_condition(
        diagnostics.cache_misses == 2, "changed primitive key did not miss the cache",
        1002);
  prims.rho = 1.25;
  prims.temperature *= 1.01;
  error = ghl_neutrino_rate_provider_compute_cell(
        &provider, &cache, &diagnostics, NULL, &prims, second);
  require_error(error, ghl_success, "changed-temperature provider call", 1005);
  require_condition(
        diagnostics.cache_misses == 3, "changed temperature key did not miss the cache",
        1005);
  prims.temperature = 0.8;
  prims.Y_e = 0.38;
  error = ghl_neutrino_rate_provider_compute_cell(
        &provider, &cache, &diagnostics, NULL, &prims, second);
  require_error(error, ghl_success, "changed-electron-fraction provider call", 1006);
  require_condition(
        diagnostics.cache_misses == 4,
        "changed electron-fraction key did not miss the cache", 1006);
  const uint64_t old_generation = provider.eos_generation;
  provider.eos_generation++;
  error = ghl_neutrino_rate_provider_compute_cell(
        &provider, &cache, &diagnostics, NULL, &prims, second);
  require_error(error, ghl_success, "changed-generation provider call", 1003);
  require_condition(
        provider.eos_generation != old_generation && diagnostics.cache_misses == 5,
        "EOS generation change did not invalidate the cache", 1003);
  validate_rate_bundle(second, 1003);

  /* The no-cache route is part of the public provider contract. */
  error = ghl_neutrino_rate_provider_compute_cell(
        &provider, NULL, NULL, NULL, &prims, second);
  require_error(error, ghl_success, "uncached reference provider call", 1004);
  validate_rate_bundle(second, 1004);
}

static void test_provider_cache_provenance(void) {
  static const char *const labels[]
        = { "channel-mask cache provenance",   "failure-policy cache provenance",
            "table-policy cache provenance",   "temperature conversion cache provenance",
            "baryon-mass cache provenance",    "charged-current cache provenance",
            "scattering cache provenance",     "pair cache provenance",
            "bremsstrahlung cache provenance", "plasmon cache provenance",
            "minimum-energy cache provenance", "equilibrium-rate cache provenance" };

  for(int variant = 0; variant < (int)(sizeof(labels) / sizeof(labels[0])); ++variant) {
    ghl_neutrino_rate_provider_context provider;
    require_error(
          ghl_neutrino_rate_provider_initialize_reference(&provider), ghl_success,
          "cache-provenance provider initialization", 1100 + variant);
    ghl_neutrino_rate_provider_cache cache;
    ghl_neutrino_rate_provider_cache_initialize(&cache);
    ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
    ghl_primitive_quantities prims = { 0 };
    prims.rho = 1.25;
    prims.temperature = 0.8;
    prims.Y_e = 0.37;
    prims.eps = prims.temperature;
    ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];

    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, &cache, &diagnostics, NULL, &prims, rates),
          ghl_success, "cache-provenance seed", 1100 + variant);
    switch(variant) {
      case 0:
        provider.channel_mask = 0;
        break;
      case 1:
        provider.failure_policy = ghl_neutrino_rate_failure_transparent;
        break;
      case 2:
        provider.table_bounds_policy = ghl_neutrino_rate_table_bounds_clamp;
        break;
      case 3:
        provider.temperature_code_to_mev = 2.0;
        break;
      case 4:
        provider.baryon_mass_code = 2.0;
        break;
      case 5:
        provider.charged_current_scale = 2.0e-2;
        break;
      case 6:
        provider.scattering_scale = 1.0e-2;
        break;
      case 7:
        provider.pair_scale = 2.0e-4;
        break;
      case 8:
        provider.bremsstrahlung_scale = 2.0e-4;
        break;
      case 9:
        provider.plasmon_scale = 2.0e-5;
        break;
      case 10:
        provider.min_mean_energy = 2.0e-12;
        break;
      case 11:
        provider.equilibrium_recovery_rate = 2.0;
        break;
      default:
        provider_test_error(
              "M1 rate-provider cache-provenance test has an unknown variant");
    }

    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, &cache, &diagnostics, NULL, &prims, rates),
          ghl_success, labels[variant], 1100 + variant);
    validate_rate_bundle(rates, 1100 + variant);
    require_condition(
          diagnostics.cache_hits == 0 && diagnostics.cache_misses == 2
                && memcmp(&cache.provider_snapshot, &provider, sizeof(provider)) == 0,
          "provider context change incorrectly reused cache", 1100 + variant);
  }
}

static void test_provider_cache_snapshot_mismatches(void) {
  /* The current provider context remains authenticated.  Only the cached
   * snapshot is stale, which lets same_provider_configuration evaluate each
   * field without being stopped by validate_provider_context. */
  static const char *const labels[]
        = { "cached temperature conversion mismatch", "cached baryon mass mismatch",
            "cached failure policy mismatch", "cached recovery rate mismatch",
            "cached EOS pointer mismatch" };
  for(int variant = 0; variant < (int)(sizeof(labels) / sizeof(labels[0])); ++variant) {
    ghl_neutrino_rate_provider_context provider;
    require_error(
          ghl_neutrino_rate_provider_initialize_reference(&provider), ghl_success,
          "snapshot-mismatch provider initialization", 1120 + variant);
    ghl_neutrino_rate_provider_cache cache;
    ghl_neutrino_rate_provider_cache_initialize(&cache);
    ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
    ghl_primitive_quantities prims = { 0 };
    prims.rho = 1.25;
    prims.temperature = 0.8;
    prims.Y_e = 0.37;
    prims.eps = prims.temperature;
    ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, &cache, &diagnostics, NULL, &prims, rates),
          ghl_success, "snapshot-mismatch cache seed", 1120 + variant);

    const ghl_eos_parameters *eos_argument = NULL;
    ghl_eos_parameters dummy_eos = { 0 };
    if(variant == 0) {
      cache.provider_snapshot.temperature_code_to_mev = 2.0;
    }
    else if(variant == 1) {
      cache.provider_snapshot.baryon_mass_code = 2.0;
    }
    else if(variant == 2) {
      cache.provider_snapshot.failure_policy = ghl_neutrino_rate_failure_transparent;
    }
    else if(variant == 3) {
      cache.provider_snapshot.equilibrium_recovery_rate = 2.0;
    }
    else {
      eos_argument = &dummy_eos;
      require_condition(
            cache.eos_snapshot == NULL,
            "snapshot-mismatch seed unexpectedly retained EOS", 1120 + variant);
    }

    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, &cache, &diagnostics, eos_argument, &prims, rates),
          ghl_success, labels[variant], 1120 + variant);
    validate_rate_bundle(rates, 1120 + variant);
    require_condition(
          diagnostics.cache_hits == 0 && diagnostics.cache_misses == 2
                && memcmp(&cache.provider_snapshot, &provider, sizeof(provider)) == 0
                && cache.eos_snapshot == eos_argument,
          "stale cache provenance was reused or not refreshed", 1120 + variant);
  }
}

static void test_provider_cache_same_rho_temperature_changed_ye(void) {
  /* Keep rho and temperature equal while changing only Ye.  This makes the
   * final equality in both cached-key predicates reachable; changing more
   * than one key would return through an earlier short-circuit operand. */
  ghl_neutrino_rate_provider_context provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_reference(&provider), ghl_success,
        "same-key provider initialization", 1130);
  ghl_neutrino_rate_provider_cache cache;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
  ghl_primitive_quantities prims = { 0 };
  prims.rho = 1.25;
  prims.temperature = 0.8;
  prims.Y_e = 0.37;
  prims.eps = prims.temperature;
  ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, NULL, &prims, rates),
        ghl_success, "same-key cache seed", 1130);
  prims.Y_e = 0.38;
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, NULL, &prims, rates),
        ghl_success, "same-rho-temperature changed-Ye call", 1131);
  validate_rate_bundle(rates, 1131);
  require_condition(
        diagnostics.cache_hits == 0 && diagnostics.cache_misses == 2
              && cache.thermo_Ye == prims.Y_e && cache.Ye == prims.Y_e,
        "changed Ye reused a cached key", 1131);
}

static void test_provider_cache_incomplete_rate_record(void) {
  /* A cache can carry valid thermodynamics before its rate bundle has been
   * published.  Keep the provenance and thermo key valid so same_rates_key()
   * must reject only the missing rates_valid flag. */
  ghl_neutrino_rate_provider_context provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_reference(&provider), ghl_success,
        "incomplete-cache provider initialization", 1140);
  ghl_neutrino_rate_provider_cache cache;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  cache.thermo_valid = true;
  cache.thermo_rho = 1.25;
  cache.thermo_T = 0.8;
  cache.thermo_Ye = 0.37;
  cache.muhat = 0.0;
  cache.mu_e = 0.0;
  cache.mu_p = 0.0;
  cache.mu_n = 0.0;
  cache.X_n = 0.63;
  cache.X_p = 0.37;
  cache.provider_snapshot = provider;
  cache.eos_snapshot = NULL;

  ghl_primitive_quantities prims = { 0 };
  prims.rho = cache.thermo_rho;
  prims.temperature = cache.thermo_T;
  prims.Y_e = cache.thermo_Ye;
  prims.eps = prims.temperature;
  ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
  ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, NULL, &prims, rates),
        ghl_success, "incomplete-cache rate publication", 1140);
  validate_rate_bundle(rates, 1140);
  require_condition(
        diagnostics.cache_hits == 0 && diagnostics.cache_misses == 1 && cache.rates_valid
              && cache.thermo_valid,
        "incomplete cache was treated as a complete rate record", 1140);
}

static void test_recovery_and_transactional_failures(void) {
  ghl_neutrino_rate_provider_context provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_reference(&provider), ghl_success,
        "recovery provider initialization", 2000);
  ghl_neutrino_rate_provider_cache cache;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
  ghl_primitive_quantities prims = { 0 };
  prims.rho = 1.0;
  prims.temperature = 1.0;
  prims.Y_e = 0.5;
  prims.eps = 1.0;
  ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(rates);
  const ghl_error_codes_t initial_error = ghl_neutrino_rate_provider_compute_cell(
        &provider, &cache, &diagnostics, NULL, &prims, rates);
  require_error(initial_error, ghl_success, "recovery cache seed", 2000);
  validate_rate_bundle(rates, 2000);

  const ghl_neutrino_rate_provider_cache cache_before = cache;
  const ghl_m1_neutrino_rates rates_before[ghl_m1_neutrino_species_count]
        = { rates[0], rates[1], rates[2] };
  ghl_primitive_quantities invalid_prims = prims;
  invalid_prims.rho = NAN;

  ghl_neutrino_rate_provider_diagnostics abort_diagnostics = { 0 };
  ghl_m1_neutrino_rates abort_rates[ghl_m1_neutrino_species_count];
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    abort_rates[species] = rates[species];
  }
  const ghl_error_codes_t abort_error = ghl_neutrino_rate_provider_compute_cell(
        &provider, &cache, &abort_diagnostics, NULL, &invalid_prims, abort_rates);
  require_error(
        abort_error, ghl_error_m1_microphysics_failure,
        "return-error-policy invalid primitive", 2001);
  require_condition(
        memcmp(&cache, &cache_before, sizeof(cache)) == 0,
        "return-error-policy failure changed the cache", 2001);
  require_condition(
        same_rate_bundle(abort_rates, rates_before),
        "return-error-policy failure changed output rates", 2001);
  require_condition(
        abort_diagnostics.failures == 1 && abort_diagnostics.last_error == abort_error,
        "return-error-policy failure diagnostics are incomplete", 2001);

  provider.failure_policy = ghl_neutrino_rate_failure_transparent;
  ghl_neutrino_rate_provider_cache transparent_cache;
  ghl_neutrino_rate_provider_cache_initialize(&transparent_cache);
  ghl_neutrino_rate_provider_diagnostics transparent_diagnostics = { 0 };
  ghl_m1_neutrino_rates transparent_rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(transparent_rates);
  const ghl_m1_neutrino_rates transparent_before[ghl_m1_neutrino_species_count]
        = { transparent_rates[0], transparent_rates[1], transparent_rates[2] };
  const ghl_neutrino_rate_provider_cache transparent_cache_before = transparent_cache;
  const ghl_error_codes_t transparent_error = ghl_neutrino_rate_provider_compute_cell(
        &provider, &transparent_cache, &transparent_diagnostics, NULL, &invalid_prims,
        transparent_rates);
  require_error(
        transparent_error, ghl_error_m1_microphysics_failure, "transparent recovery",
        2002);
  require_condition(
        same_rate_bundle(transparent_rates, transparent_before)
              && memcmp(
                       &transparent_cache, &transparent_cache_before,
                       sizeof(transparent_cache))
                       == 0
              && transparent_diagnostics.last_recovery == ghl_neutrino_rate_recovery_none
              && transparent_diagnostics.transparent_recoveries == 0
              && transparent_diagnostics.failures == 1
              && transparent_diagnostics.last_error == transparent_error,
        "transparent recovery fabricated targets without validated thermo", 2002);

  provider.failure_policy = ghl_neutrino_rate_failure_equilibrium;
  ghl_neutrino_rate_provider_diagnostics equilibrium_diagnostics = { 0 };
  ghl_m1_neutrino_rates equilibrium_rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(equilibrium_rates);
  const ghl_m1_neutrino_rates equilibrium_before[ghl_m1_neutrino_species_count]
        = { equilibrium_rates[0], equilibrium_rates[1], equilibrium_rates[2] };
  const ghl_error_codes_t equilibrium_error = ghl_neutrino_rate_provider_compute_cell(
        &provider, NULL, &equilibrium_diagnostics, NULL, &invalid_prims,
        equilibrium_rates);
  require_error(
        equilibrium_error, ghl_error_m1_microphysics_failure, "equilibrium recovery",
        2003);
  require_condition(
        same_rate_bundle(equilibrium_rates, equilibrium_before)
              && equilibrium_diagnostics.last_recovery == ghl_neutrino_rate_recovery_none
              && equilibrium_diagnostics.equilibrium_recoveries == 0
              && equilibrium_diagnostics.failures == 1
              && equilibrium_diagnostics.last_error == equilibrium_error,
        "equilibrium recovery fabricated targets without validated thermo", 2003);

  /* The recovery itself must not require diagnostics storage.  This also
   * exercises the successful equilibrium publication path's NULL optional
   * diagnostics arm. */
  ghl_m1_neutrino_rates equilibrium_no_diagnostics[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(equilibrium_no_diagnostics);
  const ghl_m1_neutrino_rates
        equilibrium_no_diagnostics_before[ghl_m1_neutrino_species_count]
        = { equilibrium_no_diagnostics[0], equilibrium_no_diagnostics[1],
            equilibrium_no_diagnostics[2] };
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, NULL, NULL, NULL, &invalid_prims, equilibrium_no_diagnostics),
        ghl_error_m1_microphysics_failure, "equilibrium recovery without diagnostics",
        2025);
  require_condition(
        same_rate_bundle(equilibrium_no_diagnostics, equilibrium_no_diagnostics_before),
        "recovery without physical thermo changed rates", 2025);

  /* Provider-context errors are non-recoverable and are checked before cache
   * lookup or any microphysics calculation. */
  require_error(
        ghl_neutrino_rate_provider_initialize_reference(&provider), ghl_success,
        "multiplicity provider initialization", 2006);
  provider.nu_x_multiplicity = 3.0;
  ghl_neutrino_rate_provider_cache bad_context_cache;
  ghl_neutrino_rate_provider_cache_initialize(&bad_context_cache);
  ghl_neutrino_rate_provider_diagnostics bad_context_diagnostics = { 0 };
  ghl_m1_neutrino_rates bad_context_rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(bad_context_rates);
  const ghl_neutrino_rate_provider_cache bad_context_cache_before = bad_context_cache;
  const ghl_m1_neutrino_rates bad_context_rates_before[ghl_m1_neutrino_species_count]
        = { bad_context_rates[0], bad_context_rates[1], bad_context_rates[2] };
  const ghl_error_codes_t context_error = ghl_neutrino_rate_provider_compute_cell(
        &provider, &bad_context_cache, &bad_context_diagnostics, NULL, &prims,
        bad_context_rates);
  require_error(
        context_error, ghl_error_m1_microphysics_failure, "invalid multiplicity context",
        2006);
  require_condition(
        memcmp(&bad_context_cache, &bad_context_cache_before, sizeof(bad_context_cache))
                    == 0
              && same_rate_bundle(bad_context_rates, bad_context_rates_before),
        "invalid context changed transactional outputs", 2006);

  require_error(
        ghl_neutrino_rate_provider_compute_cell(NULL, NULL, NULL, NULL, &prims, rates),
        ghl_error_m1_null_pointer, "NULL provider", 2007);
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, NULL, NULL, NULL, NULL, rates),
        ghl_error_m1_null_pointer, "NULL primitive input", 2008);
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, NULL, NULL, NULL, &prims, NULL),
        ghl_error_m1_null_pointer, "NULL rate output", 2009);

  test_invalid_provider_contexts(&prims);
  test_nonfinite_provider_contexts(&prims);
  test_invalid_primitive_keys();
}

static void test_recovery_publication_and_post_thermo_failures(void) {
  ghl_primitive_quantities prims = { 0 };
  prims.rho = 1.0;
  prims.temperature = 1.0;
  prims.Y_e = 0.5;
  prims.eps = 1.0;

  /* The synthetic channel scale overflows a rate product after finite
   * synthetic thermo has been validated. Recovery retains its same-cell
   * equilibrium moments and removes interactions without changing the error. */
  ghl_neutrino_rate_provider_context transparent_provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_reference(&transparent_provider),
        ghl_success, "same-cell transparent recovery initialization", 2020);
  transparent_provider.channel_mask = ghl_neutrino_rate_channel_charged_current;
  transparent_provider.charged_current_scale = DBL_MAX;
  transparent_provider.failure_policy = ghl_neutrino_rate_failure_transparent;
  m1_test_reference_provider_context target_provider = transparent_provider;
  target_provider.channel_mask = 0;
  ghl_m1_neutrino_rates expected_targets[ghl_m1_neutrino_species_count];
  compute_expected_reference_rates(&target_provider, &prims, expected_targets);

  m1_test_reference_provider_cache transparent_cache;
  m1_test_reference_provider_cache_initialize(&transparent_cache);
  const m1_test_reference_provider_cache transparent_cache_before = transparent_cache;
  ghl_neutrino_rate_provider_diagnostics transparent_diagnostics = { 0 };
  ghl_m1_neutrino_rates transparent_rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(transparent_rates);
  const ghl_error_codes_t transparent_error = ghl_neutrino_rate_provider_compute_cell(
        &transparent_provider, &transparent_cache, &transparent_diagnostics, NULL,
        &prims, transparent_rates);
  require_error(
        transparent_error, ghl_error_m1_microphysics_failure,
        "same-cell transparent recovery original error", 2020);
  validate_rate_bundle(transparent_rates, 2020);
  require_condition(
        memcmp(&transparent_cache, &transparent_cache_before, sizeof(transparent_cache))
                    == 0
              && transparent_diagnostics.failures == 1
              && transparent_diagnostics.last_error == transparent_error
              && transparent_diagnostics.transparent_recoveries == 1
              && transparent_diagnostics.last_recovery
                       == ghl_neutrino_rate_recovery_transparent,
        "same-cell transparent recovery did not preserve transaction metadata", 2020);
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    require_condition(
          transparent_rates[species].n_eq == expected_targets[species].n_eq
                && transparent_rates[species].J_eq == expected_targets[species].J_eq
                && transparent_rates[species].mean_energy
                         == expected_targets[species].mean_energy
                && transparent_rates[species].eta_N == 0.0
                && transparent_rates[species].eta_E == 0.0
                && transparent_rates[species].kappa_a_N == 0.0
                && transparent_rates[species].kappa_a_E == 0.0,
          "transparent recovery did not retain validated same-cell moments", 2020);
  }

  /* Equilibrium recovery uses the same physical cache targets and applies its
   * recovery opacity only after reconstructing those targets successfully. */
  ghl_neutrino_rate_provider_context equilibrium_provider = transparent_provider;
  equilibrium_provider.failure_policy = ghl_neutrino_rate_failure_equilibrium;
  equilibrium_provider.equilibrium_recovery_rate = 0.25;
  m1_test_reference_provider_context equilibrium_target_provider = equilibrium_provider;
  equilibrium_target_provider.channel_mask = 0;
  compute_expected_reference_rates(
        &equilibrium_target_provider, &prims, expected_targets);
  ghl_neutrino_rate_provider_diagnostics equilibrium_diagnostics = { 0 };
  ghl_m1_neutrino_rates equilibrium_rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(equilibrium_rates);
  const ghl_error_codes_t equilibrium_error = ghl_neutrino_rate_provider_compute_cell(
        &equilibrium_provider, NULL, &equilibrium_diagnostics, NULL, &prims,
        equilibrium_rates);
  require_error(
        equilibrium_error, ghl_error_m1_microphysics_failure,
        "same-cell equilibrium recovery original error", 2021);
  validate_rate_bundle(equilibrium_rates, 2021);
  require_condition(
        equilibrium_diagnostics.failures == 1
              && equilibrium_diagnostics.last_error == equilibrium_error
              && equilibrium_diagnostics.equilibrium_recoveries == 1
              && equilibrium_diagnostics.last_recovery
                       == ghl_neutrino_rate_recovery_equilibrium,
        "same-cell equilibrium recovery was not diagnosed", 2021);
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    require_condition(
          equilibrium_rates[species].n_eq == expected_targets[species].n_eq
                && equilibrium_rates[species].J_eq == expected_targets[species].J_eq
                && equilibrium_rates[species].mean_energy
                         == expected_targets[species].mean_energy
                && equilibrium_rates[species].kappa_a_N == 0.25
                && equilibrium_rates[species].kappa_a_E == 0.25
                && equilibrium_rates[species].kappa_tr == 0.25
                && equilibrium_rates[species].eta_N
                         == 0.25 * expected_targets[species].n_eq
                && equilibrium_rates[species].eta_E
                         == 0.25 * expected_targets[species].J_eq,
          "equilibrium recovery changed physical targets", 2021);
  }
}

#undef ghl_neutrino_rate_provider_compute_cell
#undef ghl_neutrino_rate_provider_cache_initialize
#undef ghl_neutrino_rate_provider_initialize_reference
#undef ghl_neutrino_rate_provider_cache
#undef ghl_neutrino_rate_provider_context

#ifndef GHL_DISABLE_HDF5
static void test_production_provider_context_validation(
      const ghl_neutrino_rate_provider_context *restrict baseline,
      const ghl_eos_parameters *restrict eos,
      const ghl_primitive_quantities *restrict prims) {
  static const char *const labels[] = { "production multiplicity",
                                        "production unknown channel mask",
                                        "production failure policy below enum",
                                        "production failure policy above enum",
                                        "production table bounds policy below enum",
                                        "production table bounds policy above enum",
                                        "production recovery rate zero",
                                        "production recovery rate NaN",
                                        "production NULL EOS",
                                        "production hybrid EOS",
                                        "production unknown table family" };
  for(int variant = 0; variant < (int)(sizeof(labels) / sizeof(labels[0])); ++variant) {
    ghl_neutrino_rate_provider_context provider = *baseline;
    ghl_eos_parameters bad_eos = *eos;
    const ghl_eos_parameters *eos_argument = &bad_eos;
    switch(variant) {
      case 0:
        provider.nu_x_multiplicity = 3.0;
        break;
      case 1:
        provider.channel_mask |= 1 << 20;
        break;
      case 2:
        provider.failure_policy = (ghl_neutrino_rate_failure_policy_t)-1;
        break;
      case 3:
        provider.failure_policy = (ghl_neutrino_rate_failure_policy_t)3;
        break;
      case 4:
        provider.table_bounds_policy = (ghl_neutrino_rate_table_bounds_policy_t)-1;
        break;
      case 5:
        provider.table_bounds_policy = (ghl_neutrino_rate_table_bounds_policy_t)2;
        break;
      case 6:
        provider.equilibrium_recovery_rate = 0.0;
        break;
      case 7:
        provider.equilibrium_recovery_rate = NAN;
        break;
      case 8:
        eos_argument = NULL;
        break;
      case 9:
        bad_eos.eos_type = ghl_eos_hybrid;
        break;
      case 10:
        bad_eos.table_type = ghl_eos_table_unknown;
        break;
      default:
        provider_test_error("Unknown production provider context variant");
    }
    ghl_neutrino_rate_provider_cache cache;
    ghl_neutrino_rate_provider_cache_initialize(&cache);
    const ghl_neutrino_rate_provider_cache cache_before = cache;
    ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
    initialize_sentinel_rates(rates);
    const ghl_m1_neutrino_rates rates_before[ghl_m1_neutrino_species_count]
          = { rates[0], rates[1], rates[2] };
    ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
    const ghl_error_codes_t error = ghl_neutrino_rate_provider_compute_cell(
          &provider, &cache, &diagnostics, eos_argument, prims, rates);
    require_error(
          error, ghl_error_m1_microphysics_failure, labels[variant], 3060 + variant);
    require_condition(
          memcmp(&cache, &cache_before, sizeof(cache)) == 0
                && same_rate_bundle(rates, rates_before) && diagnostics.failures == 1
                && diagnostics.last_error == error,
          "production context rejection changed transactional outputs", 3060 + variant);
  }
}

static void test_production_recovery_publication_overflow(void) {
  ghl_neutrino_rate_provider_context provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_default(&provider), ghl_success,
        "overflow recovery provider initialization", 3862);
  provider.channel_mask = ghl_neutrino_rate_channel_pair;
  provider.failure_policy = ghl_neutrino_rate_failure_equilibrium;
  ghl_eos_parameters eos = { 0 };
  eos.eos_type = ghl_eos_tabulated;
  eos.table_type = ghl_eos_table_stellarcollapse;
  eos.table_rho_min = 1.e-301;
  eos.table_rho_max = 1.0;
  eos.table_T_min = 1.e-120;
  eos.table_T_max = 1.e61;
  eos.table_Y_e_max = 1.0;
  const ghl_primitive_quantities prims
        = { .rho = 1.e-300, .temperature = 1.e60, .Y_e = 0.5 };
  const double recovery_rates[] = { 1.0, DBL_MAX };
  for(size_t i = 0; i < sizeof(recovery_rates) / sizeof(recovery_rates[0]); ++i) {
    provider.equilibrium_recovery_rate = recovery_rates[i];
    /* Existing public-cache fixture convention: thermodynamics are supplied
     * independently of table interpolation, while rate assembly is real. */
    ghl_neutrino_rate_provider_cache cache = { 0 };
    cache.thermo_valid = true;
    cache.thermo_rho = prims.rho;
    cache.thermo_T = prims.temperature;
    cache.thermo_Ye = prims.Y_e;
    cache.X_n = 1.0;
    cache.provider_snapshot = provider;
    cache.eos_snapshot = &eos;
    unsigned char cache_before[sizeof(cache)];
    memcpy(cache_before, &cache, sizeof(cache));
    ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
    initialize_sentinel_rates(rates);
    unsigned char rates_before[sizeof(rates)];
    memcpy(rates_before, rates, sizeof(rates));
    ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, &cache, &diagnostics, &eos, &prims, rates),
          ghl_error_m1_microphysics_failure, "pair failure with equilibrium recovery",
          3862);
    require_condition(
          memcmp(cache_before, &cache, sizeof(cache)) == 0 && diagnostics.failures == 1
                && diagnostics.last_error == ghl_error_m1_microphysics_failure,
          "failed primary assembly committed a cache record", 3862);
    if(i == 0) {
      validate_rate_bundle(rates, 3862);
      require_condition(
            diagnostics.last_recovery == ghl_neutrino_rate_recovery_equilibrium
                  && diagnostics.equilibrium_recoveries == 1,
            "representable equilibrium recovery was not published", 3862);
    }
    else {
      require_condition(
            memcmp(rates_before, rates, sizeof(rates)) == 0
                  && diagnostics.last_recovery == ghl_neutrino_rate_recovery_none
                  && diagnostics.equilibrium_recoveries == 0,
            "overflowed equilibrium recovery published rates or recovery status", 3862);
    }
  }
}

static void test_production_recovery_validation_failures(void) {
  ghl_neutrino_rate_provider_context provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_default(&provider), ghl_success,
        "validation recovery provider initialization", 3863);
  provider.channel_mask = 0;
  provider.failure_policy = ghl_neutrino_rate_failure_transparent;
  ghl_eos_parameters eos = { 0 };
  eos.eos_type = ghl_eos_tabulated;
  eos.table_type = ghl_eos_table_stellarcollapse;
  eos.table_rho_min = 1.e-301;
  eos.table_rho_max = 1.0;
  eos.table_T_min = 1.e-120;
  eos.table_T_max = 2.0;
  eos.table_Y_e_max = 1.0;
  const int saved_rounding = fegetround();
  /* Tiny T rejects the mask-zero recovery assembly itself. With ordinary T
   * and upward rounding, real assembly succeeds but the rate validator
   * rejects the unsupported rounding mode, including its recovery target. */
  const double temperatures[] = { 1.e-110, 1.0 };
  for(size_t i = 0; i < sizeof(temperatures) / sizeof(temperatures[0]); ++i) {
    const ghl_primitive_quantities prims
          = { .rho = 1.e-300, .temperature = temperatures[i], .Y_e = 0.5 };
    ghl_neutrino_rate_provider_cache cache = { 0 };
    cache.thermo_valid = true;
    cache.thermo_rho = prims.rho;
    cache.thermo_T = prims.temperature;
    cache.thermo_Ye = prims.Y_e;
    cache.X_n = 1.0;
    cache.provider_snapshot = provider;
    cache.eos_snapshot = &eos;
    unsigned char cache_before[sizeof(cache)];
    memcpy(cache_before, &cache, sizeof(cache));
    ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
    initialize_sentinel_rates(rates);
    unsigned char rates_before[sizeof(rates)];
    memcpy(rates_before, rates, sizeof(rates));
    ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
    require_condition(
          fesetround(i == 0 ? FE_TONEAREST : FE_UPWARD) == 0,
          "could not select provider validation rounding mode", 3863);
    const ghl_error_codes_t error = ghl_neutrino_rate_provider_compute_cell(
          &provider, &cache, &diagnostics, &eos, &prims, rates);
    require_condition(
          fesetround(saved_rounding) == 0, "could not restore provider rounding mode",
          3863);
    require_error(
          error, ghl_error_m1_microphysics_failure,
          "recovery assembly or candidate validation failure", 3863);
    require_condition(
          memcmp(cache_before, &cache, sizeof(cache)) == 0
                && memcmp(rates_before, rates, sizeof(rates)) == 0
                && diagnostics.failures == 1 && diagnostics.last_error == error
                && diagnostics.last_recovery == ghl_neutrino_rate_recovery_none
                && diagnostics.transparent_recoveries == 0,
          "failed recovery validation published rates, cache or status", 3863);
  }
}

static void test_production_conversion_boundaries(void) {
  ghl_neutrino_rate_provider_context provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_nrpyleakage(&provider), ghl_success,
        "conversion provider initialization", 3820);
  provider.channel_mask = 0;
  ghl_eos_parameters eos = { 0 };
  eos.eos_type = ghl_eos_tabulated;
  eos.table_type = ghl_eos_table_stellarcollapse;
  eos.table_rho_min = 1.e-301;
  eos.table_rho_max = 1.;
  eos.table_T_min = 1.e-120;
  eos.table_T_max = 1.e61;
  eos.table_Y_e_min = 0.;
  eos.table_Y_e_max = 1.;
  const double temperatures[] = { 1.e-110, 1.e-80, 8., 1.e60, 1.e-70 };
  for(size_t i = 0; i < sizeof(temperatures) / sizeof(temperatures[0]); ++i) {
    ghl_primitive_quantities prims = { 0 };
    prims.rho = 1.e-300;
    prims.temperature = temperatures[i];
    prims.Y_e = 0.5;
    ghl_neutrino_rate_provider_cache cache = { 0 };
    cache.thermo_valid = true;
    cache.thermo_rho = prims.rho;
    cache.thermo_T = prims.temperature;
    cache.thermo_Ye = prims.Y_e;
    cache.X_n = 1.;
    cache.X_p = 0.;
    cache.provider_snapshot = provider;
    cache.eos_snapshot = &eos;
    /* Public cached thermodynamics are validated again at the raw boundary. */
    if(i == 2) {
      cache.mu_e = NAN;
    }
    if(i == 3) {
      cache.mu_e = 700.0 * prims.temperature;
    }
    const ghl_neutrino_rate_provider_cache before = cache;
    ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
    initialize_sentinel_rates(rates);
    ghl_m1_neutrino_rates unchanged[ghl_m1_neutrino_species_count];
    memcpy(unchanged, rates, sizeof(rates));
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, &cache, NULL, &eos, &prims, rates),
          i < 3 ? ghl_error_m1_microphysics_failure : ghl_success,
          "converted equilibrium representability", 3821 + (int)i);
    if(i < 3) {
      require_condition(
            memcmp(&cache, &before, sizeof(cache)) == 0
                  && same_rate_bundle(rates, unchanged),
            "conversion failure changed outputs", 3821 + (int)i);
    }
    else {
      validate_rate_bundle(rates, 3821 + (int)i);
      for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
        require_condition(
              rates[species].n_eq > 0.0 && rates[species].J_eq > 0.0
                    && isfinite(rates[species].mean_energy)
                    && rates[species].mean_energy > 0.0,
              "extreme equilibrium moment is not representable", 3821 + (int)i);
      }
    }
  }
}

static void test_production_equilibrium_moments(
      const ghl_neutrino_rate_provider_context *restrict provider) {
  /* Keep this regression independent of the optional table's temperature
   * coverage.  The seeded cache is a valid thermodynamic fixture and makes
   * compute_cell exercise production assembly without an EOS interpolation. */
  static const double temperatures[] = { 0.1, 0.3, 0.4, 1.0 };
  const double rho = 1.0;
  const double Ye = 0.5;
  ghl_eos_parameters eos = { 0 };
  eos.eos_type = ghl_eos_tabulated;
  eos.table_type = ghl_eos_table_stellarcollapse;
  eos.table_rho_min = 0.5;
  eos.table_rho_max = 2.0;
  eos.table_T_min = 0.05;
  eos.table_T_max = 2.0;
  eos.table_Y_e_min = 0.25;
  eos.table_Y_e_max = 0.75;

  const double L0 = NRPyLeakage_units_geom_to_cgs_L;
  const double E0
        = NRPyLeakage_units_geom_to_cgs_M * NRPyLeakage_c_light * NRPyLeakage_c_light;
  const double volume = L0 * L0 * L0;
  const double time = NRPyLeakage_units_geom_to_cgs_T;
  const double energy_density_conversion = 1.602176634e-6 * volume / E0;
  const double number_emissivity_conversion = volume * time;
  const double energy_emissivity_conversion = 1.602176634e-6 * volume * time / E0;

  ghl_m1_parameters m1_params = { 0 };
  require_error(
        ghl_m1_initialize(
              1.0e-10, 1.0e-30, 1.0e-8, 1.0e-6, 1.0e-12, 20, 1.0e-10, &m1_params),
        ghl_success, "equilibrium-source M1 initialization", 3050);
  ghl_metric_quantities metric = { 0 };
  metric.lapse = 1.0;
  metric.lapseinv = 1.0;
  metric.lapseinv2 = 1.0;
  metric.detgamma = 1.0;
  metric.sqrt_detgamma = 1.0;
  for(int i = 0; i < 3; ++i) {
    metric.gammaDD[i][i] = 1.0;
    metric.gammaUU[i][i] = 1.0;
  }
  ghl_m1_neutrino_parameters nu_params = { 0 };

  for(size_t temperature_index = 0;
      temperature_index < sizeof(temperatures) / sizeof(temperatures[0]);
      ++temperature_index) {
    const double T = temperatures[temperature_index];
    const int case_index = 3051 + (int)temperature_index;
    ghl_primitive_quantities prims = { 0 };
    prims.rho = rho;
    prims.temperature = T;
    prims.Y_e = Ye;
    prims.eps = 1.0;

    ghl_m1_nrpyleakage_thermo_state thermo;
    require_error(
          ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
                rho, Ye, T, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, &thermo),
          ghl_success, "equilibrium fixture thermodynamic state", case_index);
    const double eta[ghl_m1_nrpyleakage_species_count] = { 0.0, 0.0, 0.0 };
    ghl_m1_nrpyleakage_raw_rates raw;
    require_error(
          ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(&thermo, eta, &raw),
          ghl_success, "independent equilibrium raw rates", case_index);

    ghl_neutrino_rate_provider_cache cache;
    ghl_neutrino_rate_provider_cache_initialize(&cache);
    cache.thermo_valid = true;
    cache.thermo_rho = rho;
    cache.thermo_T = T;
    cache.thermo_Ye = Ye;
    cache.muhat = 0.0;
    cache.mu_e = 0.0;
    cache.mu_p = 0.0;
    cache.mu_n = 0.0;
    cache.X_n = 0.5;
    cache.X_p = 0.5;
    cache.provider_snapshot = *provider;
    cache.eos_snapshot = &eos;
    ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
    ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                provider, &cache, &diagnostics, &eos, &prims, rates),
          ghl_success, "seeded-cache production provider call", case_index);
    validate_rate_bundle(rates, case_index);
    require_condition(
          diagnostics.cache_misses == 1 && diagnostics.cache_hits == 0
                && cache.thermo_valid && cache.rates_valid,
          "equilibrium fixture did not preserve its seeded thermodynamics", case_index);
    const int nux = ghl_m1_neutrino_nux;
    const double expected_n
          = provider->nu_x_multiplicity * raw.species[nux].n_eq_cgs * volume;
    const double expected_J = provider->nu_x_multiplicity * raw.species[nux].J_eq_mev_cgs
                              * energy_density_conversion;
    const double expected_mean = expected_J / expected_n;
    const double expected_eta_N
          = provider->nu_x_multiplicity
            * (raw.species[nux].eta_N_brems_cgs + raw.species[nux].eta_N_pair_cgs
               + raw.species[nux].eta_N_plasmon_cgs)
            * number_emissivity_conversion;
    const double expected_eta_E
          = provider->nu_x_multiplicity
            * (raw.species[nux].eta_E_brems_mev_cgs + raw.species[nux].eta_E_pair_mev_cgs
               + raw.species[nux].eta_E_plasmon_mev_cgs)
            * energy_emissivity_conversion;
    const double physical_mean_energy_mev = rates[nux].mean_energy * E0 / 1.602176634e-6;
    const bool below_physical_threshold = T <= 0.3;
    require_condition(
          (below_physical_threshold ? raw.species[nux].mean_energy_mev < 1.0
                                    : raw.species[nux].mean_energy_mev > 1.0)
                && provider_values_close(rates[nux].n_eq, expected_n)
                && provider_values_close(rates[nux].J_eq, expected_J)
                && provider_values_close(rates[nux].mean_energy, expected_mean),
          "production provider changed the independent raw FD moments", case_index);
    require_condition(
          isfinite(expected_eta_N) && expected_eta_N > 0.0 && isfinite(expected_eta_E)
                && expected_eta_E > 0.0,
          "independent production emissivity is not positive", case_index);
    require_condition(
          provider_values_close(
                physical_mean_energy_mev, raw.species[nux].mean_energy_mev),
          "production mean energy did not cross the physical 1-MeV threshold",
          case_index);
    if(T > 0.3) {
      require_condition(
            physical_mean_energy_mev > 1.0,
            "production mean energy was clipped at 1 MeV", case_index);
    }

    /* Construct the source state from the independent raw FD targets.  The
     * published bundle is used only as the source operator's frozen rates and
     * for checking the relative source residual. */
    ghl_m1_neutrino_state state
          = { .N = expected_n, .E = expected_J, .F = { 0.0, 0.0, 0.0 } };
    ghl_m1_sources sources = { 0 };
    double N_source = NAN;
    require_error(
          ghl_m1_compute_neutrino_interaction_sources(
                &m1_params, &nu_params, &metric, &prims, &state, &rates[nux], &sources,
                &N_source),
          ghl_success, "production equilibrium interaction source", case_index);
    const double source_tolerance = 1.0e-12;
    require_condition(
          fabs(sources.S_E) <= source_tolerance * fmax(expected_eta_E, DBL_MIN)
                && fabs(N_source) <= source_tolerance * fmax(expected_eta_N, DBL_MIN),
          "physical equilibrium produced a nonzero relative source", case_index);
    for(int i = 0; i < 3; ++i) {
      require_condition(
            fabs(sources.S[i]) <= source_tolerance * fmax(expected_eta_E, DBL_MIN),
            "physical equilibrium produced a momentum source", case_index);
    }
  }
}

static void test_production_zero_mask_cold_table(void) {
  require_condition(
        create_provider_fixture(true, false),
        "could not create cold degenerate provider fixture", 3129);

  ghl_eos_parameters eos = { 0 };
  eos.eos_type = ghl_eos_tabulated;
  eos.table_type = ghl_eos_table_stellarcollapse;
  eos.clean_sound_speed = true;
  active_provider_fixture_eos = &eos;
  require_error(
        ghl_initialize_tabulated_eos_functions_and_params(
              owned_provider_fixture_path, 1.0e-7, -1.0, -1.0, 0.5, -1.0, -1.0, 0.1,
              -1.0, -1.0, &eos),
        ghl_success, "cold degenerate EOS fixture initialization", 3129);
  require_condition(
        provider_values_close(eos.table_T_min, 0.05),
        "cold degenerate fixture does not begin at 0.05 MeV", 3129);

  ghl_neutrino_rate_provider_context provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_nrpyleakage(&provider), ghl_success,
        "cold table production provider initialization", 3129);
  provider.channel_mask = 0;
  ghl_primitive_quantities prims = { 0 };
  prims.rho = sqrt(eos.table_rho_min * eos.table_rho_max);
  prims.temperature = 0.05;
  prims.Y_e = 0.5;
  prims.eps = 1.0;

  ghl_neutrino_rate_provider_cache cache;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
  ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(rates);
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, &eos, &prims, rates),
        ghl_success, "cold degenerate zero-mask EOS-backed provider call", 3129);
  validate_rate_bundle(rates, 3129);
  require_condition(
        cache.thermo_valid && cache.rates_valid
              && provider_values_close(cache.mu_e, 40.0)
              && provider_values_close(cache.muhat, 40.0),
        "cold provider fixture did not reproduce mu_e=muhat=40 MeV", 3129);
  require_condition(
        diagnostics.active_channel_mask == 0
              && diagnostics.beta_kirchhoff_mismatch_valid[ghl_m1_neutrino_nue]
              && diagnostics.beta_kirchhoff_mismatch_valid[ghl_m1_neutrino_anue]
              && !diagnostics.beta_kirchhoff_mismatch_valid[ghl_m1_neutrino_nux]
              && isfinite(
                    diagnostics.beta_kirchhoff_relative_mismatch[ghl_m1_neutrino_nue])
              && isfinite(
                    diagnostics.beta_kirchhoff_relative_mismatch[ghl_m1_neutrino_anue]),
        "zero-mask production call lost beta/Kirchhoff diagnostics", 3129);
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    require_condition(
          isfinite(rates[species].n_eq) && rates[species].n_eq > 0.0
                && isfinite(rates[species].J_eq) && rates[species].J_eq > 0.0
                && isfinite(rates[species].mean_energy)
                && rates[species].mean_energy > 0.0,
          "zero-mask production equilibrium targets are invalid", 3129);
    require_condition(
          rates[species].eta_N == 0.0 && rates[species].eta_E == 0.0
                && rates[species].kappa_a_N == 0.0 && rates[species].kappa_a_E == 0.0
                && rates[species].kappa_a_N_cc == 0.0 && rates[species].kappa_s == 0.0
                && rates[species].kappa_tr == 0.0 && rates[species].eta_N_cc == 0.0,
          "zero-mask production call published a channel rate", 3129);
    for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
      require_condition(
            rates[species].eta_N_pair[process] == 0.0
                  && rates[species].eta_E_pair[process] == 0.0,
            "zero-mask production call published a pair-process rate", 3129);
    }
  }

  ghl_tabulated_free_memory(&eos);
  active_provider_fixture_eos = NULL;
  cleanup_provider_fixture();
}

static void test_production_recovery_from_generated_table(void) {
  require_condition(
        create_provider_fixture(false, true),
        "could not create high-temperature recovery fixture", 3130);
  ghl_eos_parameters eos = { 0 };
  eos.eos_type = ghl_eos_tabulated;
  eos.table_type = ghl_eos_table_stellarcollapse;
  eos.clean_sound_speed = true;
  active_provider_fixture_eos = &eos;
  require_error(
        ghl_initialize_tabulated_eos_functions_and_params(
              owned_provider_fixture_path, 1.0e-7, -1.0, -1.0, 0.5, -1.0, -1.0, 1.5e39,
              -1.0, -1.0, &eos),
        ghl_success, "high-temperature recovery EOS initialization", 3130);

  ghl_primitive_quantities prims = { 0 };
  prims.rho = sqrt(eos.table_rho_min * eos.table_rho_max);
  prims.temperature = sqrt(eos.table_T_min * eos.table_T_max);
  prims.Y_e = 0.5 * (eos.table_Y_e_min + eos.table_Y_e_max);
  prims.eps = 1.0;

  ghl_neutrino_rate_provider_context provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_nrpyleakage(&provider), ghl_success,
        "generated-table recovery provider initialization", 3130);
  provider.channel_mask = 0;
  ghl_neutrino_rate_provider_cache target_cache;
  ghl_neutrino_rate_provider_cache_initialize(&target_cache);
  ghl_m1_neutrino_rates physical_targets[ghl_m1_neutrino_species_count];
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &target_cache, NULL, &eos, &prims, physical_targets),
        ghl_success, "mask-zero generated-table physical targets", 3130);
  validate_rate_bundle(physical_targets, 3130);
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    require_condition(
          physical_targets[species].n_eq > 0.0 && physical_targets[species].J_eq > 0.0,
          "generated-table physical recovery targets are not positive", 3130);
  }

  /* A successful recovery publishes rates without a diagnostics output. */
  {
    ghl_m1_neutrino_rates silent_recovery[ghl_m1_neutrino_species_count];
    initialize_sentinel_rates(silent_recovery);
    ghl_neutrino_rate_provider_context silent_provider = provider;
    silent_provider.channel_mask = ghl_neutrino_rate_channel_pair;
    silent_provider.failure_policy = ghl_neutrino_rate_failure_transparent;
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &silent_provider, NULL, NULL, &eos, &prims, silent_recovery),
          ghl_error_m1_microphysics_failure,
          "generated-table recovery without diagnostics", 3134);
    validate_rate_bundle(silent_recovery, 3134);
  }

  static const ghl_neutrino_rate_failure_policy_t recovery_policies[]
        = { ghl_neutrino_rate_failure_transparent,
            ghl_neutrino_rate_failure_equilibrium };
  for(size_t policy_index = 0;
      policy_index < sizeof(recovery_policies) / sizeof(recovery_policies[0]);
      ++policy_index) {
    provider.channel_mask = ghl_neutrino_rate_channel_pair;
    provider.failure_policy = recovery_policies[policy_index];
    ghl_neutrino_rate_provider_cache cache;
    ghl_neutrino_rate_provider_cache_initialize(&cache);
    const ghl_neutrino_rate_provider_cache cache_before = cache;
    ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
    ghl_m1_neutrino_rates recovered[ghl_m1_neutrino_species_count];
    initialize_sentinel_rates(recovered);
    const ghl_error_codes_t error = ghl_neutrino_rate_provider_compute_cell(
          &provider, &cache, &diagnostics, &eos, &prims, recovered);
    const int case_index = 3131 + (int)policy_index;
    require_error(
          error, ghl_error_m1_microphysics_failure,
          "pair-kernel failure with generated-table recovery", case_index);
    validate_rate_bundle(recovered, case_index);
    require_condition(
          memcmp(&cache, &cache_before, sizeof(cache)) == 0 && diagnostics.failures == 1
                && diagnostics.last_error == error
                && diagnostics.last_recovery
                         == (policy_index == 0 ? ghl_neutrino_rate_recovery_transparent
                                               : ghl_neutrino_rate_recovery_equilibrium),
          "generated-table recovery changed cache or error metadata", case_index);
    if(policy_index == 0) {
      require_condition(
            diagnostics.transparent_recoveries == 1
                  && diagnostics.equilibrium_recoveries == 0,
            "transparent generated-table recovery was not counted", case_index);
    }
    else {
      require_condition(
            diagnostics.transparent_recoveries == 0
                  && diagnostics.equilibrium_recoveries == 1,
            "equilibrium generated-table recovery was not counted", case_index);
    }
    for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
      require_condition(
            recovered[species].n_eq == physical_targets[species].n_eq
                  && recovered[species].J_eq == physical_targets[species].J_eq
                  && recovered[species].mean_energy
                           == physical_targets[species].mean_energy,
            "generated-table recovery changed physical equilibrium moments", case_index);
      if(policy_index == 0) {
        require_condition(
              recovered[species].kappa_a_N == 0.0 && recovered[species].kappa_a_E == 0.0
                    && recovered[species].kappa_s == 0.0
                    && recovered[species].eta_N == 0.0
                    && recovered[species].eta_E == 0.0,
              "transparent generated-table recovery published interactions", case_index);
      }
      else {
        require_condition(
              recovered[species].kappa_a_N == provider.equilibrium_recovery_rate
                    && recovered[species].kappa_a_E == provider.equilibrium_recovery_rate
                    && recovered[species].eta_N
                             == provider.equilibrium_recovery_rate
                                      * physical_targets[species].n_eq
                    && recovered[species].eta_E
                             == provider.equilibrium_recovery_rate
                                      * physical_targets[species].J_eq,
              "equilibrium generated-table recovery lost its physical targets",
              case_index);
        if(recovered[species].lepton_weight != 0.0) {
          require_condition(
                recovered[species].kappa_a_N_cc == provider.equilibrium_recovery_rate
                      && recovered[species].eta_N_cc
                               == provider.equilibrium_recovery_rate
                                        * physical_targets[species].n_eq,
                "equilibrium recovery broke electron lepton exchange attribution",
                case_index);
        }
      }
    }
  }
  ghl_tabulated_free_memory(&eos);
  active_provider_fixture_eos = NULL;
  cleanup_provider_fixture();
}

static void test_production_representability_transaction(
      const ghl_neutrino_rate_provider_context *restrict provider) {
  /* Bypass EOS interpolation through a valid cached thermodynamic record so
   * the raw representability boundary is exercised without table access. */
  const double T = 1.0e-100;
  ghl_eos_parameters eos = { 0 };
  eos.eos_type = ghl_eos_tabulated;
  eos.table_type = ghl_eos_table_stellarcollapse;
  eos.table_rho_min = 0.5;
  eos.table_rho_max = 2.0;
  eos.table_T_min = 1.0e-101;
  eos.table_T_max = 1.0e-99;
  eos.table_Y_e_min = 0.25;
  eos.table_Y_e_max = 0.75;
  ghl_primitive_quantities prims = { 0 };
  prims.rho = 1.0;
  prims.temperature = T;
  prims.Y_e = 0.5;
  prims.eps = 1.0;

  ghl_neutrino_rate_provider_cache cache;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  cache.thermo_valid = true;
  cache.thermo_rho = prims.rho;
  cache.thermo_T = prims.temperature;
  cache.thermo_Ye = prims.Y_e;
  cache.muhat = 0.0;
  cache.mu_e = 0.0;
  cache.mu_p = 0.0;
  cache.mu_n = 0.0;
  cache.X_n = 0.5;
  cache.X_p = 0.5;
  cache.provider_snapshot = *provider;
  cache.eos_snapshot = &eos;
  const ghl_neutrino_rate_provider_cache cache_before = cache;

  ghl_m1_nrpyleakage_thermo_state thermo;
  require_error(
        ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
              prims.rho, prims.Y_e, T, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, &thermo),
        ghl_success, "representability thermodynamic state", 3060);
  const double eta[ghl_m1_nrpyleakage_species_count] = { 0.0, 0.0, 0.0 };
  /* The strict raw kernel rejects this unrepresentable positive FD target
   * before production assembly.  Keep that boundary explicit: this test
   * verifies public error propagation and transactionality, not an assembly
   * result that the kernel did not produce. */
  const ghl_error_codes_t raw_error = ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(
        &thermo, eta, &(ghl_m1_nrpyleakage_raw_rates){ 0 });
  require_condition(
        raw_error != ghl_success,
        "representability fixture did not fail at the raw boundary", 3060);

  ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(rates);
  const ghl_m1_neutrino_rates rates_before[ghl_m1_neutrino_species_count]
        = { rates[0], rates[1], rates[2] };
  ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
  const ghl_error_codes_t error = ghl_neutrino_rate_provider_compute_cell(
        provider, &cache, &diagnostics, &eos, &prims, rates);
  require_error(error, raw_error, "representability failure", 3060);
  require_condition(
        memcmp(&cache, &cache_before, sizeof(cache)) == 0
              && same_rate_bundle(rates, rates_before),
        "representability failure changed staged outputs", 3060);
  require_condition(
        diagnostics.failures == 1 && diagnostics.last_error == error,
        "representability failure diagnostics are incomplete", 3060);
}

static void test_production_provider_regressions(void) {
  ghl_neutrino_rate_provider_context provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_nrpyleakage(&provider), ghl_success,
        "unconditional production provider initialization", 3049);
  test_production_equilibrium_moments(&provider);
  test_production_zero_mask_cold_table();
  test_production_recovery_from_generated_table();
  test_production_representability_transaction(&provider);
}

/* This is a sufficient rejection proof, not a model of the rate backend.
 * For the negative FD branch, F2--F5 contain exp(eta) as their first
 * factor. A zero exponential therefore makes a required positive moment
 * unrepresentable, regardless of the enabled interaction channels.
 * Do not infer expected failures from the provider's return code. Leave the
 * rounding boundary and any other unexplained failures on the success path.
 */
static bool table_fd_tail_unrepresentable(const double eta) {
  const double zero_rounding_boundary = log(nextafter(0.0, 1.0)) - log(2.0);
  return isfinite(eta)
         && eta < nextafter(zero_rounding_boundary, -INFINITY) && exp(eta) == 0.0;
}

static const char *table_fd_tail_rejection_reason(
      const double T, const double muhat, const double mu_e, const int channel_mask) {
  if(table_fd_tail_unrepresentable(-fabs((mu_e - muhat) / T))) {
    return "neutrino equilibrium FD exponential underflow";
  }
  /* The pair kernel evaluates both signs of mu_e*(1/T), independently of
   * the neutrino degeneracy. Its negative electron/positron moment is also
   * required, but only when that channel is enabled. */
  if((channel_mask & ghl_neutrino_rate_channel_pair) != 0
     && table_fd_tail_unrepresentable(-fabs(mu_e * (1.0 / T)))) {
    return "pair electron/positron FD exponential underflow";
  }
  return NULL;
}

static void test_table_fd_tail_classifier(void) {
  const double min_subnormal = nextafter(0.0, 1.0);
  const double boundary = log(min_subnormal) - log(2.0);
  require_condition(
        table_fd_tail_unrepresentable(-818.5536117177104),
        "SFHo negative equilibrium tail was not classified", 3061);
  require_condition(
        !table_fd_tail_unrepresentable(0.0)
              && !table_fd_tail_unrepresentable(log(min_subnormal))
              && !table_fd_tail_unrepresentable(boundary)
              && !table_fd_tail_unrepresentable(nextafter(boundary, -INFINITY))
              && !table_fd_tail_unrepresentable(NAN),
        "equilibrium tail classifier excused an unproved rejection", 3061);
  require_condition(
        table_fd_tail_rejection_reason(1.0, 940.0, 940.0, ghl_neutrino_rate_channel_pair)
                    != NULL
              && table_fd_tail_rejection_reason(1.0, 940.0, 940.0, 0) == NULL
              && table_fd_tail_rejection_reason(1.0, 818.5536117177104, 0.0, 0) != NULL,
        "FD tail classifier ignored required-moment channel semantics", 3061);
}

static void test_table_provider(const char *restrict table_path) {
  test_table_fd_tail_classifier();
  ghl_eos_parameters eos = { 0 };
  active_provider_fixture_eos = &eos;
  eos.eos_type = ghl_eos_tabulated;
  eos.table_type = ghl_eos_table_stellarcollapse;
  eos.clean_sound_speed = true;
  const ghl_error_codes_t eos_error = ghl_initialize_tabulated_eos_functions_and_params(
        table_path, 1.0e-7, -1.0, -1.0, 0.5, -1.0, -1.0, 2.0, -1.0, -1.0, &eos);
  require_error(eos_error, ghl_success, "table initialization", 3000);
  require_condition(
        isfinite(eos.table_rho_min) && eos.table_rho_min > 0.0
              && eos.table_rho_min < eos.table_rho_max && isfinite(eos.table_T_min)
              && eos.table_T_min > 0.0 && eos.table_T_min < eos.table_T_max
              && eos.table_Y_e_min >= 0.0 && eos.table_Y_e_min < eos.table_Y_e_max,
        "table initialization published invalid bounds", 3000);

  ghl_neutrino_rate_provider_context provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_nrpyleakage(&provider), ghl_success,
        "table-backed provider initialization", 3000);
  ghl_neutrino_rate_provider_context default_provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_default(&default_provider), ghl_success,
        "default production provider initialization", 3000);
  require_condition(
        memcmp(&default_provider, &provider, sizeof(provider)) == 0,
        "default initializer did not select the installed production provider", 3000);
  ghl_primitive_quantities context_prims = { 0 };
  context_prims.rho = sqrt(eos.table_rho_min * eos.table_rho_max);
  context_prims.temperature = sqrt(eos.table_T_min * eos.table_T_max);
  context_prims.Y_e = 0.5 * (eos.table_Y_e_min + eos.table_Y_e_max);
  context_prims.eps = 1.0;
  test_production_provider_context_validation(&provider, &eos, &context_prims);

  /* Production-route input and cache-identity boundaries. */
  {
    ghl_m1_neutrino_rates boundary_rates[ghl_m1_neutrino_species_count];
    initialize_sentinel_rates(boundary_rates);
    const ghl_m1_neutrino_rates boundary_before[ghl_m1_neutrino_species_count]
          = { boundary_rates[0], boundary_rates[1], boundary_rates[2] };

    /* Optional-output guards accept NULL directly. */
    ghl_neutrino_rate_provider_cache_initialize(NULL);
    require_error(
          ghl_neutrino_rate_provider_initialize_nrpyleakage(NULL),
          ghl_error_m1_null_pointer, "production NULL provider initialization", 3002);

    /* NULL provider, primitives, and rates arms at the public boundary. */
    ghl_neutrino_rate_provider_cache boundary_cache;
    ghl_neutrino_rate_provider_cache_initialize(&boundary_cache);
    ghl_neutrino_rate_provider_diagnostics boundary_diagnostics = { 0 };
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                NULL, &boundary_cache, &boundary_diagnostics, &eos, &context_prims,
                boundary_rates),
          ghl_error_m1_null_pointer, "production NULL provider", 3002);
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, &boundary_cache, &boundary_diagnostics, &eos, NULL,
                boundary_rates),
          ghl_error_m1_null_pointer, "production NULL primitives", 3003);
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, &boundary_cache, &boundary_diagnostics, &eos, &context_prims,
                NULL),
          ghl_error_m1_null_pointer, "production NULL rates", 3004);
    require_condition(
          same_rate_bundle(boundary_rates, boundary_before),
          "production NULL boundary changed rates", 3004);

    /* A nonfinite nu_x multiplicity fails the context validation arm. */
    ghl_neutrino_rate_provider_context nan_multiplicity = provider;
    nan_multiplicity.nu_x_multiplicity = NAN;
    ghl_neutrino_rate_provider_diagnostics multiplicity_diagnostics = { 0 };
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &nan_multiplicity, &boundary_cache, &multiplicity_diagnostics, &eos,
                &context_prims, boundary_rates),
          ghl_error_m1_microphysics_failure, "production NaN multiplicity", 3005);

    /* Every invalid rho/Ye arm is a separate condition of the shared input
     * validation. */
    static const struct {
      double rho;
      double Ye;
      const char *label;
    } invalid_inputs[] = { { NAN, 0.5, "production NaN rho" },
                           { 0.0, 0.5, "production nonpositive rho" },
                           { -1.0, 0.5, "production negative rho" },
                           { 1.0, NAN, "production NaN electron fraction" },
                           { 1.0, -0.5, "production negative electron fraction" },
                           { 1.0, 1.5, "production super-unit electron fraction" } };
    for(size_t variant = 0; variant < sizeof(invalid_inputs) / sizeof(invalid_inputs[0]);
        ++variant) {
      ghl_primitive_quantities invalid_prims = context_prims;
      invalid_prims.rho = invalid_inputs[variant].rho;
      invalid_prims.Y_e = invalid_inputs[variant].Ye;
      ghl_neutrino_rate_provider_diagnostics invalid_diagnostics = { 0 };
      require_error(
            ghl_neutrino_rate_provider_compute_cell(
                  &provider, &boundary_cache, &invalid_diagnostics, &eos, &invalid_prims,
                  boundary_rates),
            ghl_error_m1_microphysics_failure, invalid_inputs[variant].label,
            3006 + (int)variant);
      require_condition(
            same_rate_bundle(boundary_rates, boundary_before)
                  && invalid_diagnostics.failures == 1,
            "production invalid inputs were not transactional", 3006 + (int)variant);
    }

    /* A productive call without diagnostics or cache exercises the optional
     * output arms. */
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, NULL, NULL, &eos, &context_prims, boundary_rates),
          ghl_success, "production diagnostics-free call", 3012);
    validate_rate_bundle(boundary_rates, 3012);

    /* Additional isolated coverage cases: one differing cached-snapshot field
     * per variant isolates every same_provider_configuration comparison, and
     * a coherent cloned cache with only thermo_valid cleared exercises the
     * thermodynamic-validity operand.  Mutating the passed provider instead
     * would fail context validation before the cache lookup; reseeding the
     * same cache would let the refreshed snapshot short-circuit each later
     * comparison at the first field. */
    {
      ghl_neutrino_rate_provider_context seed_provider = provider;
      ghl_m1_neutrino_rates seed_rates[ghl_m1_neutrino_species_count];
      static const enum {
        isolate_channel_mask,
        isolate_failure_policy,
        isolate_table_bounds_policy,
        isolate_nu_x_multiplicity,
        isolate_eos_generation,
        isolate_equilibrium_recovery_rate
      } isolate_kinds[]
            = { isolate_channel_mask,        isolate_failure_policy,
                isolate_table_bounds_policy, isolate_nu_x_multiplicity,
                isolate_eos_generation,      isolate_equilibrium_recovery_rate };
      for(int variant = 0;
          variant < (int)(sizeof(isolate_kinds) / sizeof(isolate_kinds[0])); ++variant) {
        ghl_neutrino_rate_provider_cache isolate_cache;
        ghl_neutrino_rate_provider_cache_initialize(&isolate_cache);
        ghl_neutrino_rate_provider_diagnostics isolate_seed_diagnostics = { 0 };
        require_error(
              ghl_neutrino_rate_provider_compute_cell(
                    &seed_provider, &isolate_cache, &isolate_seed_diagnostics, &eos,
                    &context_prims, seed_rates),
              ghl_success, "production isolate-field seed", 3030 + variant);
        switch(isolate_kinds[variant]) {
          case isolate_channel_mask:
            isolate_cache.provider_snapshot.channel_mask = 0;
            break;
          case isolate_failure_policy:
            isolate_cache.provider_snapshot.failure_policy
                  = ghl_neutrino_rate_failure_transparent;
            break;
          case isolate_table_bounds_policy:
            isolate_cache.provider_snapshot.table_bounds_policy
                  = ghl_neutrino_rate_table_bounds_clamp;
            break;
          case isolate_nu_x_multiplicity:
            isolate_cache.provider_snapshot.nu_x_multiplicity = 2.0;
            break;
          case isolate_eos_generation:
            isolate_cache.provider_snapshot.eos_generation += 1;
            break;
          default:
            isolate_cache.provider_snapshot.equilibrium_recovery_rate = 0.5;
            break;
        }
        ghl_m1_neutrino_rates isolate_rates[ghl_m1_neutrino_species_count];
        ghl_neutrino_rate_provider_diagnostics isolate_diagnostics = { 0 };
        require_error(
              ghl_neutrino_rate_provider_compute_cell(
                    &seed_provider, &isolate_cache, &isolate_diagnostics, &eos,
                    &context_prims, isolate_rates),
              ghl_success, "production isolate-field miss", 3030 + variant);
        require_condition(
              isolate_diagnostics.cache_hits == 0
                    && isolate_diagnostics.cache_misses == 1,
              "production snapshot-field change reused the cache", 3030 + variant);
        validate_rate_bundle(isolate_rates, 3030 + variant);
      }

      /* Only the thermodynamic validity flag is stale; the recomputed record
       * then matches the unchanged published rates key and is reused. */
      ghl_neutrino_rate_provider_cache validity_cache;
      ghl_neutrino_rate_provider_cache_initialize(&validity_cache);
      ghl_neutrino_rate_provider_diagnostics validity_seed_diagnostics = { 0 };
      require_error(
            ghl_neutrino_rate_provider_compute_cell(
                  &seed_provider, &validity_cache, &validity_seed_diagnostics, &eos,
                  &context_prims, seed_rates),
            ghl_success, "production validity-flag seed", 3040);
      ghl_neutrino_rate_provider_cache stale_thermo_cache = validity_cache;
      stale_thermo_cache.thermo_valid = false;
      ghl_m1_neutrino_rates stale_thermo_rates[ghl_m1_neutrino_species_count];
      ghl_neutrino_rate_provider_diagnostics stale_thermo_diagnostics = { 0 };
      require_error(
            ghl_neutrino_rate_provider_compute_cell(
                  &seed_provider, &stale_thermo_cache, &stale_thermo_diagnostics, &eos,
                  &context_prims, stale_thermo_rates),
            ghl_success, "production stale thermo record call", 3041);
      require_condition(
            stale_thermo_diagnostics.cache_hits == 1,
            "production stale thermo record was reused", 3041);
    }

    /* Trailing electron-fraction operands: seed on one key, then change only
     * Ye.  The first call evaluates rho and temperature before Ye rejects the
     * thermodynamic key; after the recomputed record moves the key, the
     * successful second call evaluates the same leading operands before Ye
     * rejects the published rates key. */
    {
      ghl_neutrino_rate_provider_cache ye_key_cache;
      ghl_neutrino_rate_provider_cache_initialize(&ye_key_cache);
      ghl_neutrino_rate_provider_diagnostics ye_seed_diagnostics = { 0 };
      ghl_neutrino_rate_provider_context ye_seed_provider = provider;
      ghl_m1_neutrino_rates ye_seed_rates[ghl_m1_neutrino_species_count];
      require_error(
            ghl_neutrino_rate_provider_compute_cell(
                  &ye_seed_provider, &ye_key_cache, &ye_seed_diagnostics, &eos,
                  &context_prims, ye_seed_rates),
            ghl_success, "production trailing-Ye seed", 3042);
      ghl_primitive_quantities changed_ye_prims = context_prims;
      changed_ye_prims.Y_e = 0.5 * (context_prims.Y_e + eos.table_Y_e_max);
      ghl_m1_neutrino_rates changed_ye_rates[ghl_m1_neutrino_species_count];
      ghl_neutrino_rate_provider_diagnostics changed_ye_diagnostics = { 0 };
      require_error(
            ghl_neutrino_rate_provider_compute_cell(
                  &ye_seed_provider, &ye_key_cache, &changed_ye_diagnostics, &eos,
                  &changed_ye_prims, changed_ye_rates),
            ghl_success, "production trailing-Ye thermo miss", 3043);
      require_condition(
            changed_ye_diagnostics.cache_hits == 0
                  && changed_ye_diagnostics.cache_misses == 1,
            "production trailing-Ye call reused the thermo record", 3043);
      require_error(
            ghl_neutrino_rate_provider_compute_cell(
                  &ye_seed_provider, &ye_key_cache, &changed_ye_diagnostics, &eos,
                  &changed_ye_prims, changed_ye_rates),
            ghl_success, "production trailing-Ye rates hit", 3044);
      validate_rate_bundle(changed_ye_rates, 3044);
    }

    /* Cache-identity false edges: same cache, changed thermodynamics key. */
    ghl_neutrino_rate_provider_cache identity_cache;
    ghl_neutrino_rate_provider_cache_initialize(&identity_cache);
    ghl_m1_neutrino_rates identity_first[ghl_m1_neutrino_species_count];
    ghl_neutrino_rate_provider_diagnostics identity_diagnostics = { 0 };
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, &identity_cache, &identity_diagnostics, &eos, &context_prims,
                identity_first),
          ghl_success, "production identity cache seed", 3013);
    ghl_primitive_quantities shifted_prims = context_prims;
    shifted_prims.rho = 0.8 * context_prims.rho;
    ghl_m1_neutrino_rates identity_second[ghl_m1_neutrino_species_count];
    ghl_neutrino_rate_provider_diagnostics shifted_diagnostics = { 0 };
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, &identity_cache, &shifted_diagnostics, &eos, &shifted_prims,
                identity_second),
          ghl_success, "production identity cache miss", 3013);
    require_condition(
          identity_diagnostics.cache_misses == 1 && identity_diagnostics.cache_hits == 0,
          "production changed-thermo call did not miss the cache", 3013);
    /* Same rho, changed temperature and electron fraction drive the
     * thermodynamic-key comparison arms independently. */
    ghl_primitive_quantities shifted_temperature = context_prims;
    shifted_temperature.temperature
          = 0.5 * (context_prims.temperature + eos.table_T_max);
    ghl_m1_neutrino_rates temperature_rates[ghl_m1_neutrino_species_count];
    ghl_neutrino_rate_provider_diagnostics temperature_diagnostics = { 0 };
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, &identity_cache, &temperature_diagnostics, &eos,
                &shifted_temperature, temperature_rates),
          ghl_success, "production changed-temperature cache miss", 3013);
    ghl_primitive_quantities shifted_ye = context_prims;
    shifted_ye.Y_e = 0.5 * (context_prims.Y_e + eos.table_Y_e_max);
    ghl_m1_neutrino_rates ye_rates[ghl_m1_neutrino_species_count];
    ghl_neutrino_rate_provider_diagnostics ye_diagnostics = { 0 };
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, &identity_cache, &ye_diagnostics, &eos, &shifted_ye,
                ye_rates),
          ghl_success, "production changed-Ye cache miss", 3013);

    /* Each differing provider-configuration field breaks cache provenance on
     * its own comparison term. */
    static const struct {
      const char *label;
      void (*mutate)(ghl_neutrino_rate_provider_context *restrict);
    } configuration_variants[] = { { 0 } };
    (void)configuration_variants;
    ghl_neutrino_rate_provider_context masked_provider = provider;
    masked_provider.channel_mask = 0;
    ghl_m1_neutrino_rates masked_rates[ghl_m1_neutrino_species_count];
    ghl_neutrino_rate_provider_diagnostics masked_diagnostics = { 0 };
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &masked_provider, &identity_cache, &masked_diagnostics, &eos,
                &context_prims, masked_rates),
          ghl_success, "production mask-changed provenance", 3014);
    require_condition(
          masked_diagnostics.cache_misses == 1,
          "production mask-changed call did not miss the cache", 3014);
    ghl_neutrino_rate_provider_context bounds_provider = provider;
    bounds_provider.table_bounds_policy = ghl_neutrino_rate_table_bounds_clamp;
    ghl_m1_neutrino_rates bounds_rates[ghl_m1_neutrino_species_count];
    ghl_neutrino_rate_provider_diagnostics bounds_diagnostics = { 0 };
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &bounds_provider, &identity_cache, &bounds_diagnostics, &eos,
                &context_prims, bounds_rates),
          ghl_success, "production bounds-changed provenance", 3014);
    ghl_neutrino_rate_provider_context multiplicity_provider = provider;
    multiplicity_provider.nu_x_multiplicity = 4.0;
    ghl_neutrino_rate_provider_context policy_provider = provider;
    policy_provider.failure_policy = ghl_neutrino_rate_failure_transparent;
    ghl_m1_neutrino_rates policy_rates[ghl_m1_neutrino_species_count];
    ghl_neutrino_rate_provider_diagnostics policy_diagnostics = { 0 };
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &policy_provider, &identity_cache, &policy_diagnostics, &eos,
                &context_prims, policy_rates),
          ghl_success, "production policy-changed provenance", 3014);
    ghl_neutrino_rate_provider_context recovery_provider = provider;
    recovery_provider.equilibrium_recovery_rate = 0.5;
    ghl_m1_neutrino_rates recovery_rates[ghl_m1_neutrino_species_count];
    ghl_neutrino_rate_provider_diagnostics recovery_diagnostics = { 0 };
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &recovery_provider, &identity_cache, &recovery_diagnostics, &eos,
                &context_prims, recovery_rates),
          ghl_success, "production recovery-changed provenance", 3014);

    /* A failed context validation with a NULL diagnostics output still
     * publishes nothing and fails. */
    ghl_neutrino_rate_provider_context invalid_context_provider = provider;
    invalid_context_provider.channel_mask |= 1 << 20;
    ghl_m1_neutrino_rates invalid_context_rates[ghl_m1_neutrino_species_count];
    initialize_sentinel_rates(invalid_context_rates);
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &invalid_context_provider, &identity_cache, NULL, &eos, &context_prims,
                invalid_context_rates),
          ghl_error_m1_microphysics_failure,
          "production invalid context without diagnostics", 3016);
  }

  ghl_neutrino_rate_provider_cache cache;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
  ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
  /* Keep a populated cache even when the first random cell is rejected, and
   * ensure this invocation performs a successful table-backed validation. */
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, &eos, &context_prims, rates),
        ghl_success, "random table provider control", 3000);
  validate_rate_bundle(rates, 3000);
  /* Exercise both rejection proofs even with the ordinary generated table.
   * Cached thermodynamics are caller-owned; existing production regressions
   * use this same seam without adding a fixture or a production hook. */
  for(int tail_case = 0; tail_case < 2; ++tail_case) {
    ghl_neutrino_rate_provider_cache tail_cache = cache;
    tail_cache.rates_valid = false;
    tail_cache.mu_e = tail_case == 0 ? 0.0 : 940.0 * context_prims.temperature;
    tail_cache.muhat = tail_case == 0 ? 818.5536117177104 * context_prims.temperature
                                    : tail_cache.mu_e;
    const ghl_neutrino_rate_provider_cache tail_cache_before = tail_cache;
    require_condition(
          table_fd_tail_rejection_reason(
                context_prims.temperature, tail_cache.muhat, tail_cache.mu_e,
                provider.channel_mask) != NULL,
          "cached FD tail fixture lacks an independent rejection proof", 3062 + tail_case);
    ghl_m1_neutrino_rates tail_rates[ghl_m1_neutrino_species_count];
    initialize_sentinel_rates(tail_rates);
    ghl_m1_neutrino_rates tail_rates_before[ghl_m1_neutrino_species_count];
    memcpy(tail_rates_before, tail_rates, sizeof(tail_rates));
    ghl_neutrino_rate_provider_diagnostics tail_diagnostics = { 0 };
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, &tail_cache, &tail_diagnostics, &eos, &context_prims, tail_rates),
          ghl_error_m1_microphysics_failure, "cached required FD tail rejection",
          3062 + tail_case);
    require_condition(
          memcmp(&tail_cache, &tail_cache_before, sizeof(tail_cache)) == 0
                && same_rate_bundle(tail_rates, tail_rates_before)
                && tail_diagnostics.failures == 1
                && tail_diagnostics.last_error == ghl_error_m1_microphysics_failure
                && tail_diagnostics.last_recovery == ghl_neutrino_rate_recovery_none,
          "cached FD tail rejection violated the failure transaction", 3062 + tail_case);
  }
  m1_test_rng table_rng = { .state = UINT64_C(0x4d315f5441424c45) };

  for(int case_index = 0; case_index < TABLE_RANDOM_CASES; ++case_index) {
    const double rho_fraction = m1_test_rng_between(&table_rng, 0.05, 0.95);
    const double temperature_fraction = m1_test_rng_between(&table_rng, 0.05, 0.95);
    const double ye_fraction = m1_test_rng_between(&table_rng, 0.05, 0.95);
    ghl_primitive_quantities prims = { 0 };
    prims.rho
          = exp(log(eos.table_rho_min)
                + rho_fraction * (log(eos.table_rho_max) - log(eos.table_rho_min)));
    prims.temperature
          = exp(log(eos.table_T_min)
                + temperature_fraction * (log(eos.table_T_max) - log(eos.table_T_min)));
    prims.Y_e
          = eos.table_Y_e_min + ye_fraction * (eos.table_Y_e_max - eos.table_Y_e_min);
    prims.eps = 1.0;
    double muhat, mu_e, mu_p, mu_n, X_n, X_p;
    require_error(
          ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T(
                &eos, prims.rho, prims.Y_e, prims.temperature, &muhat, &mu_e,
                &mu_p, &mu_n, &X_n, &X_p),
          ghl_success, "random table classifier thermodynamics", 3001 + case_index);
    const double negative_eta = -fabs((mu_e - muhat) / prims.temperature);
    const char *const rejection_reason = table_fd_tail_rejection_reason(
          prims.temperature, muhat, mu_e, provider.channel_mask);
    const bool expected_rejection = rejection_reason != NULL;
    initialize_sentinel_rates(rates);
    ghl_m1_neutrino_rates rates_before[ghl_m1_neutrino_species_count];
    memcpy(rates_before, rates, sizeof(rates));
    const ghl_neutrino_rate_provider_cache cache_before = cache;
    const ghl_neutrino_rate_provider_diagnostics diagnostics_before = diagnostics;
    const ghl_error_codes_t error = ghl_neutrino_rate_provider_compute_cell(
          &provider, &cache, &diagnostics, &eos, &prims, rates);
    if(error != (expected_rejection ? ghl_error_m1_microphysics_failure : ghl_success)) {
      fprintf(stderr,
              "Table sample %d: rho=%.17g T=%.17g Ye=%.17g eta_tail=%.17g (%s)\n",
              case_index, prims.rho, prims.temperature, prims.Y_e, negative_eta,
              expected_rejection ? rejection_reason : "no rejection proof");
    }
    require_error(
          error, expected_rejection ? ghl_error_m1_microphysics_failure : ghl_success,
          "random table provider call", 3001 + case_index);
    if(expected_rejection) {
      require_condition(
            memcmp(&cache, &cache_before, sizeof(cache)) == 0
                  && same_rate_bundle(rates, rates_before),
            "random table rejection changed staged outputs", 3001 + case_index);
      require_condition(
            diagnostics.failures == diagnostics_before.failures + 1
                  && diagnostics.last_error == error
                  && diagnostics.last_recovery == ghl_neutrino_rate_recovery_none,
            "random table rejection diagnostics are incomplete", 3001 + case_index);
      continue;
    }
    validate_rate_bundle(rates, 3001 + case_index);
    require_condition(
          diagnostics.beta_kirchhoff_mismatch_valid[ghl_m1_neutrino_nue]
                && diagnostics.beta_kirchhoff_mismatch_valid[ghl_m1_neutrino_anue],
          "table provider omitted beta Kirchhoff diagnostics", 3001 + case_index);
    for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
      require_condition(
            rates[ghl_m1_neutrino_nue].eta_N_pair[process] >= 0.0
                  && rates[ghl_m1_neutrino_anue].eta_N_pair[process] >= 0.0
                  && rates[ghl_m1_neutrino_nue].eta_E_pair[process] >= 0.0
                  && rates[ghl_m1_neutrino_anue].eta_E_pair[process] >= 0.0
                  && isfinite(rates[ghl_m1_neutrino_nue].eta_N_pair[process])
                  && isfinite(rates[ghl_m1_neutrino_anue].eta_N_pair[process])
                  && isfinite(rates[ghl_m1_neutrino_nue].eta_E_pair[process])
                  && isfinite(rates[ghl_m1_neutrino_anue].eta_E_pair[process]),
            "table provider published invalid electron pair process rates",
            3001 + case_index);
      require_condition(
            fabs(rates[ghl_m1_neutrino_nue].eta_N_pair[process]
                 - rates[ghl_m1_neutrino_anue].eta_N_pair[process])
                  <= 2.0e-12
                           * fmax(
                                 1.0, fmax(fabs(rates[ghl_m1_neutrino_nue]
                                                      .eta_N_pair[process]),
                                           fabs(rates[ghl_m1_neutrino_anue]
                                                      .eta_N_pair[process]))),
            "table provider broke shared electron number emissivity", 3001 + case_index);
      require_condition(
            rates[ghl_m1_neutrino_nux].eta_N_pair[process] == 0.0
                  && rates[ghl_m1_neutrino_nux].eta_E_pair[process] == 0.0,
            "table provider populated nu_x process arrays", 3001 + case_index);
    }
  }

  /* Each production channel is independently optional.  Exercise masks that
   * enter the raw-rate assembly with a single channel as well as the empty
   * mask, so every species/channel predicate is observed true and false. */
  const int production_channel_masks[]
        = { 0,
            ghl_neutrino_rate_channel_charged_current,
            ghl_neutrino_rate_channel_nucleon_scattering,
            ghl_neutrino_rate_channel_pair,
            ghl_neutrino_rate_channel_plasmon,
            ghl_neutrino_rate_channel_bremsstrahlung,
            ghl_neutrino_rate_channel_pair | ghl_neutrino_rate_channel_plasmon,
            ghl_neutrino_rate_channel_pair | ghl_neutrino_rate_channel_bremsstrahlung,
            ghl_neutrino_rate_channel_plasmon
                  | ghl_neutrino_rate_channel_bremsstrahlung };
  const int saved_channel_mask = provider.channel_mask;
  for(size_t mask_index = 0; mask_index < sizeof(production_channel_masks)
                                                / sizeof(production_channel_masks[0]);
      ++mask_index) {
    provider.channel_mask = production_channel_masks[mask_index];
    ghl_neutrino_rate_provider_cache mask_cache;
    ghl_neutrino_rate_provider_cache_initialize(&mask_cache);
    ghl_neutrino_rate_provider_diagnostics mask_diagnostics = { 0 };
    ghl_m1_neutrino_rates mask_rates[ghl_m1_neutrino_species_count];
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, &mask_cache, &mask_diagnostics, &eos, &context_prims,
                mask_rates),
          ghl_success, "production channel-mask assembly", 3030 + (int)mask_index);
    validate_rate_bundle(mask_rates, 3030 + (int)mask_index);
    require_condition(
          mask_diagnostics.cache_misses == 1 && mask_diagnostics.cache_hits == 0,
          "production channel-mask call did not publish a fresh bundle",
          3030 + (int)mask_index);
  }
  provider.channel_mask = saved_channel_mask;

  /* A repeated table key must reuse the cached thermodynamics and final rates. */
  ghl_primitive_quantities cached_prims = { 0 };
  cached_prims.rho = sqrt(eos.table_rho_min * eos.table_rho_max);
  cached_prims.temperature = sqrt(eos.table_T_min * eos.table_T_max);
  cached_prims.Y_e = 0.5 * (eos.table_Y_e_min + eos.table_Y_e_max);
  cached_prims.eps = 1.0;

  /* The public provider owns table lookup and then passes the resulting EOS
   * quantities through the strict adapter.  Corrupt only the interpolation
   * corners in memory, assert the rejected provider call is transactional, and
   * restore every authenticated fixture value before the next witness. */
  const int malformed_table_keys[]
        = { NRPyEOS_mu_e_key, NRPyEOS_X_n_key, NRPyEOS_X_p_key, NRPyEOS_mu_p_key,
            NRPyEOS_mu_n_key, NRPyEOS_X_n_key, NRPyEOS_X_p_key };
  const double malformed_table_values[]
        = { NAN, -1.0e-6, 1.0 + 1.0e-6, NAN, NAN, NAN, NAN };
  const char *const malformed_table_names[]
        = { "table NaN chemical potential",     "table negative neutron fraction",
            "table super-unit proton fraction", "table NaN proton potential",
            "table NaN neutron potential",      "table NaN neutron fraction",
            "table NaN proton fraction" };
  for(size_t malformed = 0;
      malformed < sizeof(malformed_table_keys) / sizeof(malformed_table_keys[0]);
      ++malformed) {
    double saved_table_corners[8];
    provider_set_table_corners(
          &eos, malformed_table_keys[malformed], malformed_table_values[malformed],
          saved_table_corners);
    ghl_neutrino_rate_provider_cache malformed_cache;
    ghl_neutrino_rate_provider_cache_initialize(&malformed_cache);
    const ghl_neutrino_rate_provider_cache malformed_cache_before = malformed_cache;
    ghl_m1_neutrino_rates malformed_rates[ghl_m1_neutrino_species_count];
    initialize_sentinel_rates(malformed_rates);
    const ghl_m1_neutrino_rates malformed_rates_before[ghl_m1_neutrino_species_count]
          = { malformed_rates[0], malformed_rates[1], malformed_rates[2] };
    ghl_neutrino_rate_provider_diagnostics malformed_diagnostics = { 0 };
    const ghl_error_codes_t malformed_table_error
          = ghl_neutrino_rate_provider_compute_cell(
                &provider, &malformed_cache, &malformed_diagnostics, &eos, &cached_prims,
                malformed_rates);
    require_error(
          malformed_table_error, ghl_error_m1_microphysics_failure,
          malformed_table_names[malformed], 3150 + malformed);
    require_condition(
          memcmp(&malformed_cache, &malformed_cache_before, sizeof(malformed_cache)) == 0
                && same_rate_bundle(malformed_rates, malformed_rates_before)
                && malformed_diagnostics.failures == 1
                && malformed_diagnostics.last_error == malformed_table_error,
          "malformed table changed public provider outputs", 3150 + malformed);
    provider_restore_table_corners(
          &eos, malformed_table_keys[malformed], saved_table_corners);
  }

  const double fraction_roundoff = nrpyleakage_fraction_roundoff_envelope();
  double saved_roundoff_table_corners[8];
  provider_set_table_corners(
        &eos, NRPyEOS_X_n_key, -fraction_roundoff, saved_roundoff_table_corners);
  ghl_neutrino_rate_provider_cache roundoff_cache;
  ghl_neutrino_rate_provider_cache_initialize(&roundoff_cache);
  ghl_neutrino_rate_provider_diagnostics roundoff_diagnostics = { 0 };
  ghl_m1_neutrino_rates roundoff_rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(roundoff_rates);
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &roundoff_cache, &roundoff_diagnostics, &eos, &cached_prims,
              roundoff_rates),
        ghl_success, "table roundoff composition", 3153);
  validate_rate_bundle(roundoff_rates, 3153);
  require_condition(
        roundoff_cache.thermo_valid && roundoff_cache.X_n == 0.0
              && roundoff_cache.X_p > 0.0 && roundoff_cache.X_p <= 1.0,
        "table roundoff composition was not normalized in the cache", 3153);
  provider_restore_table_corners(&eos, NRPyEOS_X_n_key, saved_roundoff_table_corners);

  /* Re-run the same table point with exact zero corners and an independent
   * cache.  The two published bundles must be identical after normalization. */
  double saved_exact_table_corners[8];
  provider_set_table_corners(&eos, NRPyEOS_X_n_key, 0.0, saved_exact_table_corners);
  ghl_neutrino_rate_provider_cache exact_cache;
  ghl_neutrino_rate_provider_cache_initialize(&exact_cache);
  ghl_neutrino_rate_provider_diagnostics exact_diagnostics = { 0 };
  ghl_m1_neutrino_rates exact_rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(exact_rates);
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &exact_cache, &exact_diagnostics, &eos, &cached_prims,
              exact_rates),
        ghl_success, "table exact normalized composition", 3154);
  validate_rate_bundle(exact_rates, 3154);
  require_condition(
        same_rate_bundle(roundoff_rates, exact_rates),
        "table roundoff rates differ from exact normalized rates", 3154);
  provider_restore_table_corners(&eos, NRPyEOS_X_n_key, saved_exact_table_corners);

  ghl_neutrino_rate_provider_cache_initialize(&cache);
  diagnostics = (ghl_neutrino_rate_provider_diagnostics){ 0 };
  ghl_m1_neutrino_rates first[ghl_m1_neutrino_species_count];
  ghl_m1_neutrino_rates second[ghl_m1_neutrino_species_count];
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, &eos, &cached_prims, first),
        ghl_success, "table cache seed", 3041);
  validate_rate_bundle(first, 3041);
  require_condition(
        diagnostics.cache_misses == 1 && diagnostics.cache_hits == 0 && cache.rates_valid
              && cache.thermo_valid,
        "table cache seed did not publish a complete production record", 3041);
  double seeded_beta_kirchhoff_relative_mismatch[ghl_m1_neutrino_species_count];
  bool seeded_beta_kirchhoff_mismatch_valid[ghl_m1_neutrino_species_count];
  memcpy(
        seeded_beta_kirchhoff_relative_mismatch,
        diagnostics.beta_kirchhoff_relative_mismatch,
        sizeof(seeded_beta_kirchhoff_relative_mismatch));
  memcpy(
        seeded_beta_kirchhoff_mismatch_valid, diagnostics.beta_kirchhoff_mismatch_valid,
        sizeof(seeded_beta_kirchhoff_mismatch_valid));
  diagnostics = (ghl_neutrino_rate_provider_diagnostics){ 0 };
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, &eos, &cached_prims, second),
        ghl_success, "table cache hit", 3041);
  require_condition(
        diagnostics.cache_hits == 1 && same_rate_bundle(first, second),
        "table cache hit was not exact", 3041);
  require_condition(
        memcmp(
              diagnostics.beta_kirchhoff_relative_mismatch,
              seeded_beta_kirchhoff_relative_mismatch,
              sizeof(seeded_beta_kirchhoff_relative_mismatch))
                    == 0
              && memcmp(
                       diagnostics.beta_kirchhoff_mismatch_valid,
                       seeded_beta_kirchhoff_mismatch_valid,
                       sizeof(seeded_beta_kirchhoff_mismatch_valid))
                       == 0,
        "production cache hit did not preserve its table diagnostics", 3041);

  /* A missing primitive temperature is recovered through the table's
   * inverse EOS using an independently obtained table-consistent epsilon. */
  double recovery_eps = NAN;
  require_error(
        NRPyEOS_eps_from_rho_Ye_T(
              &eos, cached_prims.rho, cached_prims.Y_e, cached_prims.temperature,
              &recovery_eps),
        ghl_success, "table-consistent recovery epsilon", 3044);
  ghl_primitive_quantities recovered_temperature = cached_prims;
  recovered_temperature.temperature = NAN;
  recovered_temperature.eps = recovery_eps;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  diagnostics = (ghl_neutrino_rate_provider_diagnostics){ 0 };
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, &eos, &recovered_temperature, rates),
        ghl_success, "table temperature recovery", 3044);
  validate_rate_bundle(rates, 3044);
  recovered_temperature.temperature = 0.0;
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, NULL, NULL, &eos, &recovered_temperature, rates),
        ghl_success, "zero table temperature recovery", 3812);
  require_condition(
        cache.thermo_valid && isfinite(cache.thermo_T) && cache.thermo_T > 0.0,
        "table temperature recovery did not publish thermodynamics", 3044);

  ghl_eos_parameters missing_ye_bounds = eos;
  missing_ye_bounds.table_Y_e_min = 0.0;
  missing_ye_bounds.table_Y_e_max = 0.0;
  missing_ye_bounds.Y_e_min = NAN;
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, NULL, NULL, &missing_ye_bounds, &cached_prims, rates),
        ghl_error_m1_microphysics_failure, "missing zero-based Ye table bounds", 3824);

  /* The legacy bound fields are a supported fallback when table-specific
   * bounds are absent.  Use the authenticated table limits as the fallback
   * values so the EOS interpolation still sees the same table. */
  ghl_eos_parameters fallback_bounds = eos;
  fallback_bounds.rho_min = eos.table_rho_min;
  fallback_bounds.rho_max = eos.table_rho_max;
  fallback_bounds.T_min = eos.table_T_min;
  fallback_bounds.T_max = eos.table_T_max;
  fallback_bounds.Y_e_min = eos.table_Y_e_min;
  fallback_bounds.Y_e_max = eos.table_Y_e_max;
  fallback_bounds.table_rho_min = 0.0;
  fallback_bounds.table_rho_max = 0.0;
  fallback_bounds.table_T_min = 0.0;
  fallback_bounds.table_T_max = 0.0;
  fallback_bounds.table_Y_e_min = -1.0;
  fallback_bounds.table_Y_e_max = 0.0;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  const ghl_neutrino_rate_provider_cache fallback_cache_before = cache;
  const ghl_m1_neutrino_rates fallback_rates_before[ghl_m1_neutrino_species_count]
        = { rates[0], rates[1], rates[2] };
  diagnostics = (ghl_neutrino_rate_provider_diagnostics){ 0 };
  const ghl_error_codes_t fallback_error = ghl_neutrino_rate_provider_compute_cell(
        &provider, &cache, &diagnostics, &fallback_bounds, &cached_prims, rates);
  require_error(
        fallback_error, ghl_error_table_max_rho, "legacy table-bound fallback", 3045);
  require_condition(
        memcmp(&cache, &fallback_cache_before, sizeof(cache)) == 0
              && same_rate_bundle(rates, fallback_rates_before),
        "invalid legacy table bounds changed transactional outputs", 3045);

  /* Malformed table metadata is rejected before the EOS callback and leaves
   * both cache and caller-owned rates unchanged. */
  ghl_eos_parameters malformed_bounds = eos;
  malformed_bounds.table_T_min = NAN;
  malformed_bounds.T_min = NAN;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  const ghl_neutrino_rate_provider_cache malformed_cache_before = cache;
  ghl_m1_neutrino_rates malformed_rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(malformed_rates);
  const ghl_m1_neutrino_rates malformed_rates_before[ghl_m1_neutrino_species_count]
        = { malformed_rates[0], malformed_rates[1], malformed_rates[2] };
  diagnostics = (ghl_neutrino_rate_provider_diagnostics){ 0 };
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, &malformed_bounds, &cached_prims,
              malformed_rates),
        ghl_error_m1_microphysics_failure, "malformed table bounds", 3046);
  require_condition(
        memcmp(&cache, &malformed_cache_before, sizeof(cache)) == 0
              && same_rate_bundle(malformed_rates, malformed_rates_before),
        "malformed table bounds changed transactional outputs", 3046);

  /* Temperature recovery validates table metadata before the inverse EOS. */
  ghl_eos_parameters malformed_recovery_bounds = eos;
  malformed_recovery_bounds.table_T_min = NAN;
  malformed_recovery_bounds.T_min = NAN;
  ghl_primitive_quantities malformed_recovery_prims = cached_prims;
  malformed_recovery_prims.temperature = NAN;
  malformed_recovery_prims.eps = 1.0;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  initialize_sentinel_rates(malformed_rates);
  diagnostics = (ghl_neutrino_rate_provider_diagnostics){ 0 };
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, &malformed_recovery_bounds,
              &malformed_recovery_prims, malformed_rates),
        ghl_error_m1_microphysics_failure, "malformed recovery temperature bounds",
        3047);

  static const char *const recovery_failure_names[]
        = { "nonfinite upper recovery temperature bound",
            "nonpositive lower recovery temperature bound",
            "reversed recovery temperature bounds", "nonfinite recovery epsilon" };
  for(size_t recovery_case = 0;
      recovery_case < sizeof(recovery_failure_names) / sizeof(recovery_failure_names[0]);
      ++recovery_case) {
    ghl_eos_parameters bad_recovery = eos;
    ghl_primitive_quantities bad_recovery_prims = cached_prims;
    bad_recovery_prims.temperature = NAN;
    bad_recovery_prims.eps = 1.0;
    switch(recovery_case) {
      case 0:
        bad_recovery.table_T_max = NAN;
        bad_recovery.T_max = NAN;
        break;
      case 1:
        bad_recovery.table_T_min = 0.0;
        bad_recovery.T_min = 0.0;
        break;
      case 2:
        bad_recovery.table_T_min = 2.0 * eos.table_T_max;
        break;
      case 3:
        bad_recovery_prims.eps = NAN;
        break;
      default:
        ghl_error("M1 rate-provider recovery test has an unknown variant\n");
    }
    ghl_neutrino_rate_provider_cache_initialize(&cache);
    initialize_sentinel_rates(malformed_rates);
    diagnostics = (ghl_neutrino_rate_provider_diagnostics){ 0 };
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                &provider, &cache, &diagnostics, &bad_recovery, &bad_recovery_prims,
                malformed_rates),
          ghl_error_m1_microphysics_failure, recovery_failure_names[recovery_case],
          3055 + (int)recovery_case);
  }

  /* A recovered temperature can also fail validation after the provisional
   * geometric temperature is formed, or fail inside the inverse EOS itself. */
  ghl_primitive_quantities out_of_bounds_recovery = cached_prims;
  out_of_bounds_recovery.temperature = NAN;
  out_of_bounds_recovery.eps = 1.0;
  out_of_bounds_recovery.rho = 0.5 * eos.table_rho_min;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  initialize_sentinel_rates(malformed_rates);
  diagnostics = (ghl_neutrino_rate_provider_diagnostics){ 0 };
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, &eos, &out_of_bounds_recovery,
              malformed_rates),
        ghl_error_m1_microphysics_failure, "out-of-bounds recovered temperature input",
        3048);

  ghl_primitive_quantities inverse_failure_prims = cached_prims;
  inverse_failure_prims.temperature = NAN;
  inverse_failure_prims.eps = DBL_MAX;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  initialize_sentinel_rates(malformed_rates);
  diagnostics = (ghl_neutrino_rate_provider_diagnostics){ 0 };
  const ghl_error_codes_t inverse_failure_error
        = ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, &eos, &inverse_failure_prims,
              malformed_rates);
  require_condition(
        inverse_failure_error != ghl_success,
        "inverse EOS failure was silently accepted", 3049);

  /* Exercise each table-bound metadata disjunct independently.  The table
   * lookup must reject malformed bounds before it can inspect any payload. */
  static const char *const malformed_bound_names[]
        = { "nonfinite lower density bound",
            "nonfinite upper density bound",
            "nonpositive lower density bound",
            "reversed density bounds",
            "nonfinite lower temperature bound",
            "nonfinite upper temperature bound",
            "nonpositive lower temperature bound",
            "reversed temperature bounds",
            "nonfinite lower electron-fraction bound",
            "nonfinite upper electron-fraction bound",
            "negative lower electron-fraction bound",
            "super-unit upper electron-fraction bound",
            "reversed electron-fraction bounds" };
  for(size_t bound_case = 0;
      bound_case < sizeof(malformed_bound_names) / sizeof(malformed_bound_names[0]);
      ++bound_case) {
    ghl_eos_parameters bad_bounds = eos;
    switch(bound_case) {
      case 0:
        bad_bounds.table_rho_min = NAN;
        bad_bounds.rho_min = NAN;
        break;
      case 1:
        bad_bounds.table_rho_max = NAN;
        bad_bounds.rho_max = NAN;
        break;
      case 2:
        bad_bounds.table_rho_min = 0.0;
        bad_bounds.rho_min = 0.0;
        break;
      case 3:
        bad_bounds.table_rho_min = 2.0 * eos.table_rho_max;
        break;
      case 4:
        bad_bounds.table_T_min = NAN;
        bad_bounds.T_min = NAN;
        break;
      case 5:
        bad_bounds.table_T_max = NAN;
        bad_bounds.T_max = NAN;
        break;
      case 6:
        bad_bounds.table_T_min = 0.0;
        bad_bounds.T_min = 0.0;
        break;
      case 7:
        bad_bounds.table_T_min = 2.0 * eos.table_T_max;
        break;
      case 8:
        bad_bounds.table_Y_e_min = NAN;
        bad_bounds.Y_e_min = NAN;
        break;
      case 9:
        bad_bounds.table_Y_e_max = NAN;
        bad_bounds.Y_e_max = NAN;
        break;
      case 10:
        bad_bounds.table_Y_e_min = -1.0;
        bad_bounds.Y_e_min = -1.0;
        break;
      case 11:
        bad_bounds.table_Y_e_max = 2.0;
        break;
      case 12:
        bad_bounds.table_Y_e_min = 0.8;
        bad_bounds.table_Y_e_max = 0.2;
        break;
      default:
        ghl_error("M1 rate-provider malformed-bound test has an unknown variant\n");
    }
    ghl_neutrino_rate_provider_cache bad_bounds_cache;
    ghl_neutrino_rate_provider_cache_initialize(&bad_bounds_cache);
    const ghl_neutrino_rate_provider_cache bad_bounds_cache_before = bad_bounds_cache;
    ghl_m1_neutrino_rates bad_bounds_rates[ghl_m1_neutrino_species_count];
    initialize_sentinel_rates(bad_bounds_rates);
    const ghl_m1_neutrino_rates bad_bounds_rates_before[ghl_m1_neutrino_species_count]
          = { bad_bounds_rates[0], bad_bounds_rates[1], bad_bounds_rates[2] };
    ghl_neutrino_rate_provider_diagnostics bad_bounds_diagnostics = { 0 };
    const ghl_error_codes_t bad_bounds_error = ghl_neutrino_rate_provider_compute_cell(
          &provider, &bad_bounds_cache, &bad_bounds_diagnostics, &bad_bounds,
          &cached_prims, bad_bounds_rates);
    require_error(
          bad_bounds_error, ghl_error_m1_microphysics_failure,
          malformed_bound_names[bound_case], 3048 + (int)bound_case);
    require_condition(
          memcmp(&bad_bounds_cache, &bad_bounds_cache_before, sizeof(bad_bounds_cache))
                      == 0
                && same_rate_bundle(bad_bounds_rates, bad_bounds_rates_before)
                && bad_bounds_diagnostics.failures == 1,
          "malformed table metadata changed transactional outputs",
          3048 + (int)bound_case);
  }

  const double out_of_bounds_values[3]
        = { 0.5 * eos.table_rho_min, 0.5 * eos.table_T_min,
            eos.table_Y_e_min > 0.0 ? 0.5 * eos.table_Y_e_min
                                    : 0.5 * (eos.table_Y_e_max + 1.0) };
  const char *const bound_names[3] = { "rho", "temperature", "Ye" };
  for(int bound = 0; bound < 3; ++bound) {
    ghl_primitive_quantities out_of_bounds = cached_prims;
    if(bound == 0) {
      out_of_bounds.rho = out_of_bounds_values[bound];
    }
    else if(bound == 1) {
      out_of_bounds.temperature = out_of_bounds_values[bound];
    }
    else {
      out_of_bounds.Y_e = out_of_bounds_values[bound];
    }
    ghl_m1_neutrino_rates unchanged[ghl_m1_neutrino_species_count];
    for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
      unchanged[species] = first[species];
    }
    const ghl_neutrino_rate_provider_cache cache_before = cache;
    diagnostics = (ghl_neutrino_rate_provider_diagnostics){ 0 };
    const ghl_error_codes_t bounds_error = ghl_neutrino_rate_provider_compute_cell(
          &provider, &cache, &diagnostics, &eos, &out_of_bounds, unchanged);
    require_error(
          bounds_error, ghl_error_m1_microphysics_failure, "abort-policy table bounds",
          3042 + bound);
    require_condition(
          memcmp(&cache, &cache_before, sizeof(cache)) == 0
                && same_rate_bundle(unchanged, first),
          "abort-policy table bounds changed outputs", 3042 + bound);
    require_condition(
          diagnostics.table_bound_hits == 1, bound_names[bound], 3042 + bound);
  }

  const double above_bounds[3] = { 2.0 * eos.table_rho_max, 2.0 * eos.table_T_max,
                                   0.5 * (eos.table_Y_e_max + 1.0) };
  require_condition(
        eos.table_Y_e_max < 1.0, "generated table cannot exercise an upper Ye bound",
        3047);
  for(int bound = 0; bound < 3; ++bound) {
    ghl_primitive_quantities out_of_bounds = cached_prims;
    if(bound == 0) {
      out_of_bounds.rho = above_bounds[bound];
    }
    else if(bound == 1) {
      out_of_bounds.temperature = above_bounds[bound];
    }
    else {
      out_of_bounds.Y_e = above_bounds[bound];
    }
    ghl_m1_neutrino_rates unchanged[ghl_m1_neutrino_species_count];
    for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
      unchanged[species] = first[species];
    }
    const ghl_neutrino_rate_provider_cache cache_before = cache;
    diagnostics = (ghl_neutrino_rate_provider_diagnostics){ 0 };
    const ghl_error_codes_t bounds_error = ghl_neutrino_rate_provider_compute_cell(
          &provider, &cache, &diagnostics, &eos, &out_of_bounds, unchanged);
    require_error(
          bounds_error, ghl_error_m1_microphysics_failure, "upper-bound table input",
          3047 + bound);
    require_condition(
          memcmp(&cache, &cache_before, sizeof(cache)) == 0
                && same_rate_bundle(unchanged, first),
          "upper-bound table input changed outputs", 3047 + bound);
    require_condition(
          diagnostics.table_bound_hits == 1, bound_names[bound], 3047 + bound);
  }

  provider.table_bounds_policy = ghl_neutrino_rate_table_bounds_clamp;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  diagnostics = (ghl_neutrino_rate_provider_diagnostics){ 0 };
  ghl_primitive_quantities out_of_bounds = cached_prims;
  out_of_bounds.rho = out_of_bounds_values[0];
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, &eos, &out_of_bounds, rates),
        ghl_success, "clamped table bounds", 3043);
  validate_rate_bundle(rates, 3043);
  require_condition(
        diagnostics.table_bound_hits > 0 && diagnostics.clamped_inputs > 0,
        "clamped table bounds were not diagnosed", 3043);

  ghl_neutrino_rate_provider_cache_initialize(&cache);
  ghl_primitive_quantities null_diagnostics_input = cached_prims;
  null_diagnostics_input.rho = out_of_bounds_values[0];
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, NULL, &eos, &null_diagnostics_input, rates),
        ghl_success, "clamped table bounds without diagnostics", 3049);
  validate_rate_bundle(rates, 3049);
  ghl_m1_neutrino_rates cached_rates[ghl_m1_neutrino_species_count];
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, NULL, &eos, &null_diagnostics_input, cached_rates),
        ghl_success, "cache hit without diagnostics", 3805);
  require_condition(
        same_rate_bundle(rates, cached_rates),
        "cache hit without diagnostics changed rates", 3805);

  ghl_tabulated_free_memory(&eos);
  active_provider_fixture_eos = NULL;
}
#endif

#ifdef GHL_DISABLE_HDF5
static void require_disabled_provider_failure(
      const ghl_neutrino_rate_provider_context *provider,
      const ghl_eos_parameters *eos,
      const ghl_primitive_quantities *prims,
      const ghl_error_codes_t expected_error,
      const int case_index) {
  /* Exercise every combination of optional cache and diagnostics pointers.
   * Neither a context error nor unavailable microphysics may publish rates. */
  for(int use_cache = 0; use_cache <= 1; ++use_cache) {
    for(int use_diagnostics = 0; use_diagnostics <= 1; ++use_diagnostics) {
      ghl_neutrino_rate_provider_cache cache;
      memset(&cache, 0xa5, sizeof(cache));
      unsigned char cache_before[sizeof(cache)];
      memcpy(cache_before, &cache, sizeof(cache));
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
      initialize_sentinel_rates(rates);
      unsigned char rates_before[sizeof(rates)];
      memcpy(rates_before, rates, sizeof(rates));
      ghl_neutrino_rate_provider_diagnostics diagnostics
            = { .failures = 7,
                .table_bound_hits = 3,
                .cache_hits = 5,
                .cache_misses = 9,
                .clamped_inputs = 4,
                .transparent_recoveries = 2,
                .equilibrium_recoveries = 6,
                .active_channel_mask = ghl_neutrino_rate_channel_pair,
                .last_error = ghl_error_m1_invalid_state,
                .last_recovery = ghl_neutrino_rate_recovery_equilibrium,
                .beta_kirchhoff_relative_mismatch = { 0.25, 0.5, 0.75 },
                .beta_kirchhoff_mismatch_valid = { true, false, true } };
      ghl_neutrino_rate_provider_diagnostics expected_diagnostics;
      memcpy(&expected_diagnostics, &diagnostics, sizeof(diagnostics));
      expected_diagnostics.failures++;
      expected_diagnostics.last_error = expected_error;
      expected_diagnostics.last_recovery = ghl_neutrino_rate_recovery_none;
      require_error(
            ghl_neutrino_rate_provider_compute_cell(
                  provider, use_cache ? &cache : NULL,
                  use_diagnostics ? &diagnostics : NULL, eos, prims, rates),
            expected_error, "disabled provider error precedence", case_index);
      require_condition(
            memcmp(cache_before, &cache, sizeof(cache)) == 0
                  && memcmp(rates_before, rates, sizeof(rates)) == 0,
            "disabled provider published cache or rates", case_index);
      if(use_diagnostics) {
        require_condition(
              memcmp(&expected_diagnostics, &diagnostics, sizeof(diagnostics)) == 0,
              "disabled provider changed unrelated diagnostics", case_index);
      }
    }
  }
}

static void test_disabled_provider_contract(void) {
  require_error(
        ghl_neutrino_rate_provider_initialize_default(NULL), ghl_error_m1_null_pointer,
        "disabled default initializer NULL output", 3840);
  require_error(
        ghl_neutrino_rate_provider_initialize_nrpyleakage(NULL),
        ghl_error_m1_null_pointer, "disabled NRPyLeakage initializer NULL output", 3841);
  ghl_neutrino_rate_provider_cache_initialize(NULL);
  ghl_neutrino_rate_provider_cache cache;
  memset(&cache, 0xa5, sizeof(cache));
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  const ghl_m1_neutrino_rates zero_rates[ghl_m1_neutrino_species_count] = { { 0 } };
  require_condition(
        !cache.thermo_valid && !cache.rates_valid && cache.eos_snapshot == NULL
              && cache.thermo_rho == 0.0 && cache.thermo_T == 0.0
              && cache.thermo_Ye == 0.0 && cache.rho == 0.0 && cache.T == 0.0
              && cache.Ye == 0.0 && cache.provider_snapshot.eos_generation == 0
              && same_rate_bundle(cache.rates, zero_rates),
        "disabled cache initializer did not clear its record", 3842);

  const ghl_neutrino_rate_provider_context valid
        = { .channel_mask = ghl_neutrino_rate_channel_charged_current,
            .failure_policy = ghl_neutrino_rate_failure_return_error,
            .table_bounds_policy = ghl_neutrino_rate_table_bounds_return_error,
            .nu_x_multiplicity = 4.0,
            .equilibrium_recovery_rate = 1.0 };
  const ghl_primitive_quantities prims = { .rho = 1.0, .temperature = 1.0, .Y_e = 0.5 };
  ghl_neutrino_rate_provider_context invalid[]
        = { valid, valid, valid, valid, valid, valid, valid, valid, valid };
  invalid[0].nu_x_multiplicity = NAN;
  invalid[1].nu_x_multiplicity = 3.0;
  invalid[2].channel_mask = ~0;
  invalid[3].failure_policy = (ghl_neutrino_rate_failure_policy_t)-1;
  invalid[4].failure_policy = ghl_neutrino_rate_failure_equilibrium + 1;
  invalid[5].table_bounds_policy = (ghl_neutrino_rate_table_bounds_policy_t)-1;
  invalid[6].table_bounds_policy = ghl_neutrino_rate_table_bounds_clamp + 1;
  invalid[7].equilibrium_recovery_rate = NAN;
  invalid[8].equilibrium_recovery_rate = 0.0;
  for(size_t i = 0; i < sizeof(invalid) / sizeof(invalid[0]); ++i) {
    require_disabled_provider_failure(
          &invalid[i], NULL, &prims, ghl_error_m1_microphysics_failure, 3843 + (int)i);
  }
  require_disabled_provider_failure(
        &valid, NULL, &prims, ghl_error_used_disabled_hdf5, 3852);
  const ghl_primitive_quantities invalid_prims
        = { .rho = NAN, .temperature = -1.0, .Y_e = NAN };
  const ghl_eos_parameters invalid_eos = { 0 };
  require_disabled_provider_failure(
        &valid, &invalid_eos, &invalid_prims, ghl_error_used_disabled_hdf5, 3853);
  require_disabled_provider_failure(
        &invalid[0], &invalid_eos, &invalid_prims, ghl_error_m1_microphysics_failure,
        3854);

  /* Each required-pointer guard precedes context checking and any diagnostic
   * mutation, including when other inputs would also be invalid. */
  for(int missing = 0; missing < 3; ++missing) {
    ghl_neutrino_rate_provider_diagnostics diagnostics
          = { .failures = 7,
              .last_error = ghl_error_m1_invalid_state,
              .last_recovery = ghl_neutrino_rate_recovery_equilibrium };
    unsigned char diagnostics_before[sizeof(diagnostics)];
    memcpy(diagnostics_before, &diagnostics, sizeof(diagnostics));
    ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
    initialize_sentinel_rates(rates);
    unsigned char rates_before[sizeof(rates)];
    memcpy(rates_before, rates, sizeof(rates));
    unsigned char cache_before[sizeof(cache)];
    memcpy(cache_before, &cache, sizeof(cache));
    require_error(
          ghl_neutrino_rate_provider_compute_cell(
                missing == 0 ? NULL : &invalid[0], &cache, &diagnostics, &invalid_eos,
                missing == 1 ? NULL : &invalid_prims, missing == 2 ? NULL : rates),
          ghl_error_m1_null_pointer, "disabled required pointer precedence", 3855);
    require_condition(
          memcmp(diagnostics_before, &diagnostics, sizeof(diagnostics)) == 0
                && memcmp(cache_before, &cache, sizeof(cache)) == 0
                && memcmp(rates_before, rates, sizeof(rates)) == 0,
          "disabled required-pointer rejection changed outputs", 3855);
  }
}

static void test_disabled_backend_entrypoints(void) {
  require_error(
        ghl_m1_neutrino_rate_backend_initialize(), ghl_error_used_disabled_hdf5,
        "disabled private backend initializer", 3856);
  double T = 12.0;
  require_error(
        ghl_m1_neutrino_rate_backend_temperature_from_eps(NULL, NAN, NAN, NAN, &T),
        ghl_error_used_disabled_hdf5, "disabled private temperature lookup", 3857);
  require_condition(T == 12.0, "disabled temperature lookup changed T", 3857);
  double thermo[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
  const double thermo_before[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
  require_error(
        ghl_m1_neutrino_rate_backend_thermo_from_T(
              NULL, NAN, NAN, NAN, &thermo[0], &thermo[1], &thermo[2], &thermo[3],
              &thermo[4], &thermo[5]),
        ghl_error_used_disabled_hdf5, "disabled private thermodynamics lookup", 3858);
  require_condition(
        memcmp(thermo, thermo_before, sizeof(thermo)) == 0,
        "disabled thermodynamics lookup changed outputs", 3858);
  ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(rates);
  unsigned char rates_before[sizeof(rates)];
  memcpy(rates_before, rates, sizeof(rates));
  double mismatch[] = { 0.25, 0.5, 0.75 };
  const double mismatch_before[] = { 0.25, 0.5, 0.75 };
  bool valid[] = { true, false, true };
  const bool valid_before[] = { true, false, true };
  require_error(
        ghl_m1_neutrino_rate_backend_assemble(
              0, 4.0, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, rates, mismatch,
              valid),
        ghl_error_used_disabled_hdf5, "disabled private rate assembly", 3859);
  require_condition(
        memcmp(rates_before, rates, sizeof(rates)) == 0
              && memcmp(mismatch_before, mismatch, sizeof(mismatch)) == 0
              && memcmp(valid_before, valid, sizeof(valid)) == 0,
        "disabled private rate assembly changed outputs", 3859);
}

static void test_disabled_provider_backend(void) {
  ghl_neutrino_rate_provider_context provider;
  memset(&provider, 0xa5, sizeof(provider));
  const ghl_neutrino_rate_provider_context provider_before = provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_default(&provider),
        ghl_error_used_disabled_hdf5, "disabled default-provider initializer", 3825);
  require_condition(
        memcmp(&provider, &provider_before, sizeof(provider)) == 0,
        "disabled default initializer changed its output", 3825);
  require_error(
        ghl_neutrino_rate_provider_initialize_nrpyleakage(&provider),
        ghl_error_used_disabled_hdf5, "disabled NRPyLeakage initializer", 3826);
  require_condition(
        memcmp(&provider, &provider_before, sizeof(provider)) == 0,
        "disabled NRPyLeakage initializer changed its output", 3826);

  provider = (ghl_neutrino_rate_provider_context){
    .channel_mask
    = ghl_neutrino_rate_channel_charged_current
      | ghl_neutrino_rate_channel_nucleon_scattering | ghl_neutrino_rate_channel_pair
      | ghl_neutrino_rate_channel_bremsstrahlung | ghl_neutrino_rate_channel_plasmon,
    .failure_policy = ghl_neutrino_rate_failure_return_error,
    .table_bounds_policy = ghl_neutrino_rate_table_bounds_return_error,
    .nu_x_multiplicity = 4.0,
    .equilibrium_recovery_rate = 1.0
  };
  ghl_neutrino_rate_provider_cache cache;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  const ghl_neutrino_rate_provider_cache cache_before = cache;
  ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(rates);
  const ghl_m1_neutrino_rates rates_before[ghl_m1_neutrino_species_count]
        = { rates[0], rates[1], rates[2] };
  const ghl_primitive_quantities prims = { .rho = 1.0, .temperature = 1.0, .Y_e = 0.5 };
  ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
  const ghl_error_codes_t error = ghl_neutrino_rate_provider_compute_cell(
        &provider, &cache, &diagnostics, NULL, &prims, rates);
  require_error(
        error, ghl_error_used_disabled_hdf5, "disabled production compute", 3827);
  require_condition(
        memcmp(&cache, &cache_before, sizeof(cache)) == 0
              && same_rate_bundle(rates, rates_before) && diagnostics.failures == 1
              && diagnostics.last_error == error,
        "disabled production call changed cache or rates", 3827);
}
#endif

#ifndef GHL_DISABLE_HDF5
/* Use the existing EOS dispatch boundary to inject successful lookups with
 * invalid results. These are defensive-boundary tests, separate from the
 * real generated-table production checks; disabled microphysics stays disabled. */
static ghl_error_codes_t provider_test_temperature_result(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Ye,
      const double eps,
      double *restrict T) {
  (void)eos;
  (void)rho;
  (void)Ye;
  *T = eps;
  return ghl_success;
}

static int provider_test_bad_chemical_potential;
static ghl_error_codes_t provider_test_thermo_result(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Ye,
      const double T,
      double *restrict muhat,
      double *restrict mu_e,
      double *restrict mu_p,
      double *restrict mu_n,
      double *restrict X_n,
      double *restrict X_p) {
  (void)eos;
  (void)rho;
  (void)Ye;
  (void)T;
  double *const chemical_potentials[] = { muhat, mu_e, mu_p, mu_n };
  for(size_t i = 0; i < sizeof(chemical_potentials) / sizeof(chemical_potentials[0]);
      ++i) {
    *chemical_potentials[i] = (int)i == provider_test_bad_chemical_potential ? NAN : 0.0;
  }
  *X_n = 0.5;
  *X_p = 0.5;
  return ghl_success;
}

static void test_backend_lookup_postconditions(void) {
  require_error(
        ghl_m1_neutrino_rate_backend_initialize(), ghl_success,
        "enabled private backend initializer", 3860);
  ghl_error_codes_t (*saved_temperature_lookup)(
        const ghl_eos_parameters *restrict, const double, const double, const double,
        double *restrict) = ghl_tabulated_compute_T_from_eps;
  ghl_tabulated_compute_T_from_eps = provider_test_temperature_result;
  const double invalid_temperatures[] = { NAN, 0.0, -1.0 };
  for(size_t i = 0; i < sizeof(invalid_temperatures) / sizeof(invalid_temperatures[0]);
      ++i) {
    double T = 12.0;
    require_error(
          ghl_m1_neutrino_rate_backend_temperature_from_eps(
                NULL, 1.0, 0.5, invalid_temperatures[i], &T),
          ghl_error_m1_microphysics_failure, "invalid successful EOS temperature", 3860);
    require_condition(T == 12.0, "invalid EOS temperature was published", 3860);
  }
  ghl_tabulated_compute_T_from_eps = saved_temperature_lookup;

  ghl_error_codes_t (*saved_thermo_lookup)(
        const ghl_eos_parameters *restrict, const double, const double, const double,
        double *restrict, double *restrict, double *restrict, double *restrict,
        double *restrict, double *restrict)
        = ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T;
  ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T = provider_test_thermo_result;
  for(provider_test_bad_chemical_potential = 0; provider_test_bad_chemical_potential < 4;
      ++provider_test_bad_chemical_potential) {
    double thermo[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
    const double before[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
    require_error(
          ghl_m1_neutrino_rate_backend_thermo_from_T(
                NULL, 1.0, 0.5, 1.0, &thermo[0], &thermo[1], &thermo[2], &thermo[3],
                &thermo[4], &thermo[5]),
          ghl_error_m1_microphysics_failure, "invalid successful EOS chemical potential",
          3861);
    require_condition(
          memcmp(thermo, before, sizeof(thermo)) == 0,
          "invalid EOS chemical potential was published", 3861);
  }
  ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T = saved_thermo_lookup;
}
#endif

#define ghl_neutrino_rate_provider_context m1_test_reference_provider_context
#define ghl_neutrino_rate_provider_cache   m1_test_reference_provider_cache
#define ghl_neutrino_rate_provider_initialize_reference \
  m1_test_reference_provider_initialize
#define ghl_neutrino_rate_provider_cache_initialize \
  m1_test_reference_provider_cache_initialize
#define ghl_neutrino_rate_provider_compute_cell m1_test_reference_provider_compute_cell

static void test_reference_cache_without_diagnostics(void) {
  ghl_neutrino_rate_provider_context provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_reference(&provider), ghl_success,
        "reference cache initialization", 3830);
  ghl_neutrino_rate_provider_cache cache = { 0 };
  ghl_primitive_quantities prims = { .rho = 1.e-4, .temperature = 1., .Y_e = 0.5 };
  ghl_m1_neutrino_rates first[ghl_m1_neutrino_species_count];
  ghl_m1_neutrino_rates second[ghl_m1_neutrino_species_count];
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, NULL, NULL, &prims, first),
        ghl_success, "reference cache seed without diagnostics", 3830);
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, NULL, NULL, &prims, second),
        ghl_success, "reference cache hit without diagnostics", 3830);
  require_condition(
        same_rate_bundle(first, second),
        "reference cache without diagnostics changed rates", 3830);
}

int main(int argc, char **argv) {
  if(argc > 2) {
    cleanup_provider_fixture();
    ghl_error(
          "Usage: %s [StellarCollapse EOS table path|--generated-fixture]\n", argv[0]);
  }

  m1_test_rng rng = { .state = UINT64_C(0x4d3150524f564944) };
  test_reference_cache_without_diagnostics();
  test_masked_kernel_validation();
  test_nrpyleakage_raw_kernel();
  test_nrpyleakage_boundary_paths();
  test_nrpyleakage_supported_thermo_boundaries();
  test_nrpyleakage_rate_overflow();
  test_nrpyleakage_kernel_representability_edges();
  test_default_provider(&rng);
  test_provider_cache_provenance();
  test_provider_cache_snapshot_mismatches();
  test_provider_cache_same_rho_temperature_changed_ye();
  test_provider_cache_incomplete_rate_record();
  test_recovery_and_transactional_failures();
  test_recovery_publication_and_post_thermo_failures();
  test_temperature_recovery();
  ghl_neutrino_rate_provider_context bad_context;
  (void)ghl_neutrino_rate_provider_initialize_reference(&bad_context);
  bad_context.nu_x_multiplicity = 0.0;
  ghl_primitive_quantities prim = { 0 };
  ghl_m1_neutrino_rates untouched_rates[ghl_m1_neutrino_species_count];
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &bad_context, NULL, NULL, NULL, &prim, untouched_rates),
        ghl_error_m1_microphysics_failure, "invalid context without diagnostics", 3813);

#ifndef GHL_DISABLE_HDF5
  test_backend_lookup_postconditions();
  test_production_conversion_boundaries();
  test_production_recovery_publication_overflow();
  test_production_recovery_validation_failures();
  test_production_provider_regressions();
  if(argc == 2) {
    const char *table_path = argv[1];
    if(strcmp(argv[1], "--generated-fixture") == 0) {
      if(!create_provider_fixture(false, false)) {
        provider_test_error("Could not create the generated provider EOS fixture");
      }
      table_path = owned_provider_fixture_path;
    }
    test_table_provider(table_path);
    if(owned_provider_fixture_path != NULL) {
      cleanup_provider_fixture();
    }
  }
#else
  test_disabled_provider_backend();
  test_disabled_provider_contract();
  test_disabled_backend_entrypoints();
#endif

  ghl_info(
        "M1 rate-provider randomized/property tests passed (%d reference cases%s)\n",
        PROVIDER_RANDOM_CASES,
#ifndef GHL_DISABLE_HDF5
        argc == 2 ? ", table-backed cases" : ""
#else
        ""
#endif
  );
  return 0;
}
