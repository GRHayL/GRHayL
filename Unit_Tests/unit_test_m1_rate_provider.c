#include "../GRHayL/Radiation/Neutrinos/ghl_m1_nrpyleakage_kernel.h"
#include "ghl_radiation.h"
#ifndef GHL_DISABLE_HDF5
#include "ghl_nrpyeos_tabulated.h"
#include <hdf5.h>
#include <unistd.h>
#endif
#include "../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h"
#include "../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_rate_helpers.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

/*
 * Deterministic property coverage for the M1 frozen-rate provider boundary.
 * The reference backend is deliberately exercised with generated primitive
 * keys, channel masks, cache reuse, context-generation changes, and each
 * recovery policy.  A table path is optional; when supplied in an HDF5 build,
 * the same test also enters the production NRPyLeakage backend.
 */

enum { PROVIDER_RANDOM_CASES = 128, TABLE_RANDOM_CASES = 32 };

typedef struct {
  uint64_t state;
} provider_rng;

static uint64_t provider_rng_next(provider_rng *restrict rng) {
  uint64_t z = (rng->state += UINT64_C(0x9e3779b97f4a7c15));
  z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
  return z ^ (z >> 31);
}

static double provider_rng_unit(provider_rng *restrict rng) {
  return (double)(provider_rng_next(rng) >> 11) * 0x1.0p-53;
}

static double provider_rng_between(
      provider_rng *restrict rng,
      const double lower,
      const double upper) {
  return lower + (upper - lower) * provider_rng_unit(rng);
}

static void require_condition(
      const bool condition,
      const char *restrict message,
      const int case_index) {
  if(!condition) {
    ghl_error("M1 rate-provider case %d: %s\n", case_index, message);
  }
}

static void require_error(
      const ghl_error_codes_t actual,
      const ghl_error_codes_t expected,
      const char *restrict operation,
      const int case_index) {
  if(actual != expected) {
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
  /* Keep this oracle identical to ghl_m1_nrpyleakage_normalize_nucleon_fractions(). */
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

/* Keep malformed-table witnesses in memory only.  Each witness changes the
 * eight interpolation corners used by the authenticated interior point and
 * restores every byte before the next production-table assertion. */
static void provider_set_table_corners(
      ghl_eos_parameters *restrict eos,
      const int key,
      const double replacement,
      double saved[8]) {
  int corner = 0;
  for(int ir = 0; ir < 2; ++ir) {
    for(int it = 0; it < 2; ++it) {
      for(int iy = 0; iy < 2; ++iy) {
        const size_t index = (size_t)NRPYEOS_IDX3D(eos, ir, it, iy, key);
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
  int corner = 0;
  for(int ir = 0; ir < 2; ++ir) {
    for(int it = 0; it < 2; ++it) {
      for(int iy = 0; iy < 2; ++iy) {
        const size_t index = (size_t)NRPYEOS_IDX3D(eos, ir, it, iy, key);
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
static bool create_provider_fixture(char *restrict path, const size_t path_size) {
  const int characters = snprintf(
        path, path_size, "/tmp/ghl_m1_rate_provider_fixture_%ld.h5", (long)getpid());
  if(characters < 0 || (size_t)characters >= path_size) {
    return false;
  }

  const double logrho[PROVIDER_FIXTURE_NRHO] = { 10.0, 11.0, 12.0 };
  const double logtemp[PROVIDER_FIXTURE_NTEMP] = { 0.0, 0.5, 1.0 };
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
        mu_e[index] = 3.0 + 0.20 * ir + 0.30 * it + 0.40 * iy;
        mu_n[index] = 1.0 + 0.10 * ir + 0.12 * it + 0.15 * iy;
        mu_p[index] = 0.2 + 0.05 * ir + 0.06 * it + 0.08 * iy;
        muhat[index] = mu_n[index] - mu_p[index];
        munu[index] = 0.1 + 0.02 * ir + 0.03 * it + 0.04 * iy;
      }
    }
  }

  const hid_t file = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if(file < 0) {
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
    remove(path);
  }
  return ok;
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
          rates[species].species == species,
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
make_primitives(provider_rng *restrict rng, ghl_primitive_quantities *restrict prims) {
  *prims = (ghl_primitive_quantities){ 0 };
  prims->rho = exp(provider_rng_between(rng, log(0.05), log(20.0)));
  prims->temperature = exp(provider_rng_between(rng, log(0.05), log(20.0)));
  prims->Y_e = provider_rng_between(rng, 0.02, 0.98);
  prims->eps = prims->temperature;
  prims->press = prims->rho * prims->eps * 0.1;
  prims->entropy = provider_rng_between(rng, 0.05, 2.0);
}

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
  static const char *const labels[] = { "backend below enum range",
                                        "backend above enum range",
                                        "unknown channel bit",
                                        "failure policy below enum range",
                                        "failure policy above enum range",
                                        "table policy below enum range",
                                        "table policy above enum range",
                                        "reference density conversion",
                                        "reference opacity conversion",
                                        "reference emissivity conversion",
                                        "nonpositive temperature conversion",
                                        "nonpositive baryon mass",
                                        "negative charged-current scale",
                                        "negative scattering scale",
                                        "negative pair scale",
                                        "negative bremsstrahlung scale",
                                        "negative plasmon scale",
                                        "nonpositive minimum mean energy",
                                        "nonpositive equilibrium recovery rate",
                                        "missing configured EOS" };
  const int variant_count = (int)(sizeof(labels) / sizeof(labels[0]));

  for(int variant = 0; variant < variant_count; ++variant) {
    ghl_neutrino_rate_provider_context provider;
    require_error(
          ghl_neutrino_rate_provider_initialize_default(&provider), ghl_success,
          "invalid-context baseline initialization", 4000 + variant);
    switch(variant) {
      case 0:
        provider.backend = (ghl_neutrino_rate_backend_t)-1;
        break;
      case 1:
        provider.backend = (ghl_neutrino_rate_backend_t)2;
        break;
      case 2:
        provider.channel_mask = 1 << 8;
        break;
      case 3:
        provider.failure_policy = (ghl_neutrino_rate_failure_policy_t)-1;
        break;
      case 4:
        provider.failure_policy = (ghl_neutrino_rate_failure_policy_t)4;
        break;
      case 5:
        provider.table_bounds_policy = (ghl_neutrino_rate_table_bounds_policy_t)-1;
        break;
      case 6:
        provider.table_bounds_policy = (ghl_neutrino_rate_table_bounds_policy_t)2;
        break;
      case 7:
        provider.rho_code_to_cgs = 2.0;
        break;
      case 8:
        provider.opacity_cgs_to_code = 2.0;
        break;
      case 9:
        provider.emissivity_cgs_to_code = 2.0;
        break;
      case 10:
        provider.temperature_code_to_mev = 0.0;
        break;
      case 11:
        provider.baryon_mass_code = 0.0;
        break;
      case 12:
        provider.charged_current_scale = -1.0;
        break;
      case 13:
        provider.scattering_scale = -1.0;
        break;
      case 14:
        provider.pair_scale = -1.0;
        break;
      case 15:
        provider.bremsstrahlung_scale = -1.0;
        break;
      case 16:
        provider.plasmon_scale = -1.0;
        break;
      case 17:
        provider.min_mean_energy = 0.0;
        break;
      case 18:
        provider.equilibrium_recovery_rate = 0.0;
        break;
      case 19:
        provider.use_tabulated_eos = true;
        break;
      default:
        ghl_error("M1 rate-provider invalid-context test has an unknown variant\n");
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
          ghl_neutrino_rate_provider_initialize_default(&provider), ghl_success,
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
        ghl_error("M1 rate-provider nonfinite-context test has an unknown variant\n");
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
        ghl_neutrino_rate_provider_initialize_default(&provider), ghl_success,
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
        ghl_neutrino_rate_provider_initialize_default(&provider), ghl_success,
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

static void test_default_provider(provider_rng *restrict rng) {
  ghl_neutrino_rate_provider_context provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_default(&provider), ghl_success,
        "default provider initialization", 0);
  require_error(
        ghl_neutrino_rate_provider_initialize_default(NULL), ghl_error_m1_null_pointer,
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
    provider.charged_current_scale = provider_rng_between(rng, 1.0e-4, 2.0e-2);
    provider.scattering_scale = provider_rng_between(rng, 1.0e-4, 2.0e-2);
    provider.pair_scale = provider_rng_between(rng, 1.0e-6, 2.0e-4);
    provider.bremsstrahlung_scale = provider_rng_between(rng, 1.0e-6, 2.0e-4);
    provider.plasmon_scale = provider_rng_between(rng, 1.0e-7, 2.0e-5);

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
        ghl_neutrino_rate_provider_initialize_default(&provider), ghl_success,
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
          ghl_neutrino_rate_provider_initialize_default(&provider), ghl_success,
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
        ghl_error("M1 rate-provider cache-provenance test has an unknown variant\n");
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
        = { "cached density conversion mismatch", "cached opacity conversion mismatch",
            "cached emissivity conversion mismatch", "cached backend mismatch",
            "cached EOS pointer mismatch" };
  for(int variant = 0; variant < (int)(sizeof(labels) / sizeof(labels[0])); ++variant) {
    ghl_neutrino_rate_provider_context provider;
    require_error(
          ghl_neutrino_rate_provider_initialize_default(&provider), ghl_success,
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
      cache.provider_snapshot.rho_code_to_cgs = 2.0;
    }
    else if(variant == 1) {
      cache.provider_snapshot.opacity_cgs_to_code = 2.0;
    }
    else if(variant == 2) {
      cache.provider_snapshot.emissivity_cgs_to_code = 2.0;
    }
    else if(variant == 3) {
      cache.provider_snapshot.backend = ghl_neutrino_rate_backend_nrpyleakage;
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
        ghl_neutrino_rate_provider_initialize_default(&provider), ghl_success,
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
        ghl_neutrino_rate_provider_initialize_default(&provider), ghl_success,
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
        ghl_neutrino_rate_provider_initialize_default(&provider), ghl_success,
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
        abort_error, ghl_error_m1_microphysics_failure, "abort-policy invalid primitive",
        2001);
  require_condition(
        memcmp(&cache, &cache_before, sizeof(cache)) == 0,
        "abort-policy failure changed the cache", 2001);
  require_condition(
        same_rate_bundle(abort_rates, rates_before),
        "abort-policy failure changed output rates", 2001);
  require_condition(
        abort_diagnostics.failures == 1 && abort_diagnostics.last_error == abort_error,
        "abort-policy failure diagnostics are incomplete", 2001);

  provider.failure_policy = ghl_neutrino_rate_failure_transparent;
  ghl_neutrino_rate_provider_cache transparent_cache;
  ghl_neutrino_rate_provider_cache_initialize(&transparent_cache);
  ghl_neutrino_rate_provider_diagnostics transparent_diagnostics = { 0 };
  ghl_m1_neutrino_rates transparent_rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(transparent_rates);
  const ghl_error_codes_t transparent_error = ghl_neutrino_rate_provider_compute_cell(
        &provider, &transparent_cache, &transparent_diagnostics, NULL, &invalid_prims,
        transparent_rates);
  require_error(transparent_error, ghl_success, "transparent recovery", 2002);
  validate_rate_bundle(transparent_rates, 2002);
  require_condition(
        transparent_diagnostics.last_recovery == ghl_neutrino_rate_recovery_transparent
              && transparent_diagnostics.transparent_recoveries == 1
              && !transparent_cache.thermo_valid && !transparent_cache.rates_valid,
        "transparent recovery did not preserve its transaction", 2002);

  provider.failure_policy = ghl_neutrino_rate_failure_equilibrium;
  ghl_neutrino_rate_provider_diagnostics equilibrium_diagnostics = { 0 };
  ghl_m1_neutrino_rates equilibrium_rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(equilibrium_rates);
  const ghl_error_codes_t equilibrium_error = ghl_neutrino_rate_provider_compute_cell(
        &provider, NULL, &equilibrium_diagnostics, NULL, &invalid_prims,
        equilibrium_rates);
  require_error(equilibrium_error, ghl_success, "equilibrium recovery", 2003);
  validate_rate_bundle(equilibrium_rates, 2003);
  require_condition(
        equilibrium_diagnostics.last_recovery == ghl_neutrino_rate_recovery_equilibrium
              && equilibrium_diagnostics.equilibrium_recoveries == 1,
        "equilibrium recovery diagnostics are incomplete", 2003);

  /* The recovery itself must not require diagnostics storage.  This also
   * exercises the successful equilibrium publication path's NULL optional
   * diagnostics arm. */
  ghl_m1_neutrino_rates equilibrium_no_diagnostics[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(equilibrium_no_diagnostics);
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, NULL, NULL, NULL, &invalid_prims, equilibrium_no_diagnostics),
        ghl_success, "equilibrium recovery without diagnostics", 2025);
  validate_rate_bundle(equilibrium_no_diagnostics, 2025);

  /* Hold-last requires an exact recovered key. An invalid primitive has no
   * recovered key, so the policy must fail closed and leave both records. */
  require_error(
        ghl_neutrino_rate_provider_initialize_default(&provider), ghl_success,
        "hold-last provider initialization", 2004);
  provider.failure_policy = ghl_neutrino_rate_failure_hold_last;
  ghl_neutrino_rate_provider_cache hold_last_cache;
  ghl_neutrino_rate_provider_cache_initialize(&hold_last_cache);
  ghl_neutrino_rate_provider_diagnostics hold_last_diagnostics = { 0 };
  ghl_m1_neutrino_rates hold_last_rates[ghl_m1_neutrino_species_count];
  const ghl_error_codes_t hold_seed_error = ghl_neutrino_rate_provider_compute_cell(
        &provider, &hold_last_cache, &hold_last_diagnostics, NULL, &prims,
        hold_last_rates);
  require_error(hold_seed_error, ghl_success, "hold-last cache seed", 2004);
  const ghl_neutrino_rate_provider_cache hold_cache_before = hold_last_cache;
  const ghl_m1_neutrino_rates hold_rates_before[ghl_m1_neutrino_species_count]
        = { hold_last_rates[0], hold_last_rates[1], hold_last_rates[2] };
  const ghl_error_codes_t hold_error = ghl_neutrino_rate_provider_compute_cell(
        &provider, &hold_last_cache, &hold_last_diagnostics, NULL, &invalid_prims,
        hold_last_rates);
  require_error(
        hold_error, ghl_error_m1_microphysics_failure, "hold-last invalid primitive",
        2005);
  require_condition(
        memcmp(&hold_last_cache, &hold_cache_before, sizeof(hold_last_cache)) == 0
              && same_rate_bundle(hold_last_rates, hold_rates_before),
        "hold-last rejection changed transactional outputs", 2005);

  /* Provider-context errors are non-recoverable and are checked before cache
   * lookup or any microphysics calculation. */
  require_error(
        ghl_neutrino_rate_provider_initialize_default(&provider), ghl_success,
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

  ghl_neutrino_rate_provider_context production = { 0 };
  require_error(
        ghl_neutrino_rate_provider_initialize_nrpyleakage(NULL),
        ghl_error_m1_null_pointer, "NULL production provider initialization", 2011);
  const ghl_error_codes_t production_init
        = ghl_neutrino_rate_provider_initialize_nrpyleakage(&production);
#ifdef GHL_DISABLE_HDF5
  require_error(
        production_init, ghl_error_used_disabled_hdf5,
        "disabled-HDF5 production provider initialization", 2010);
  ghl_neutrino_rate_provider_context disabled_backend;
  require_error(
        ghl_neutrino_rate_provider_initialize_default(&disabled_backend), ghl_success,
        "disabled-HDF5 backend-gate initialization", 2012);
  disabled_backend.backend = ghl_neutrino_rate_backend_nrpyleakage;
  ghl_neutrino_rate_provider_cache disabled_cache;
  ghl_neutrino_rate_provider_cache_initialize(&disabled_cache);
  const ghl_neutrino_rate_provider_cache disabled_cache_before = disabled_cache;
  ghl_m1_neutrino_rates disabled_rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(disabled_rates);
  const ghl_m1_neutrino_rates disabled_rates_before[ghl_m1_neutrino_species_count]
        = { disabled_rates[0], disabled_rates[1], disabled_rates[2] };
  ghl_neutrino_rate_provider_diagnostics disabled_diagnostics = { 0 };
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &disabled_backend, &disabled_cache, &disabled_diagnostics, NULL, &prims,
              disabled_rates),
        ghl_error_used_disabled_hdf5, "disabled-HDF5 NRPyLeakage backend gate", 2012);
  require_condition(
        memcmp(&disabled_cache, &disabled_cache_before, sizeof(disabled_cache)) == 0
              && same_rate_bundle(disabled_rates, disabled_rates_before)
              && disabled_diagnostics.failures == 1,
        "disabled-HDF5 backend gate changed transactional outputs", 2012);

  /* A reference context with the tabulated flag set is still a supported
   * representable context up to the HDF5 feature gate.  Supplying a non-NULL
   * EOS reaches that build-time gate instead of the earlier NULL-EOS guard. */
  disabled_backend.backend = ghl_neutrino_rate_backend_reference;
  disabled_backend.use_tabulated_eos = true;
  ghl_eos_parameters disabled_eos_argument = { 0 };
  ghl_neutrino_rate_provider_cache disabled_tabulated_cache;
  ghl_neutrino_rate_provider_cache_initialize(&disabled_tabulated_cache);
  ghl_neutrino_rate_provider_diagnostics disabled_tabulated_diagnostics = { 0 };
  ghl_m1_neutrino_rates disabled_tabulated_rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(disabled_tabulated_rates);
  const ghl_m1_neutrino_rates disabled_tabulated_before[ghl_m1_neutrino_species_count]
        = { disabled_tabulated_rates[0], disabled_tabulated_rates[1],
            disabled_tabulated_rates[2] };
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &disabled_backend, &disabled_tabulated_cache,
              &disabled_tabulated_diagnostics, &disabled_eos_argument, &prims,
              disabled_tabulated_rates),
        ghl_error_m1_microphysics_failure,
        "disabled-HDF5 non-NULL tabulated backend gate", 2013);
  require_condition(
        same_rate_bundle(disabled_tabulated_rates, disabled_tabulated_before)
              && disabled_tabulated_diagnostics.failures == 1,
        "disabled-HDF5 non-NULL backend gate changed rates", 2013);
#else
  require_error(
        production_init, ghl_success, "HDF5 production provider initialization", 2010);
  require_condition(
        production.backend == ghl_neutrino_rate_backend_nrpyleakage
              && production.use_tabulated_eos,
        "production initializer selected the wrong backend", 2010);
#endif
}

static void test_recovery_publication_and_post_thermo_failures(void) {
  ghl_primitive_quantities prims = { 0 };
  prims.rho = 1.0;
  prims.temperature = 1.0;
  prims.Y_e = 0.5;
  prims.eps = 1.0;

  /* A failed input with optional diagnostics must still publish a valid
   * transparent recovery, while no cache is available to update. */
  ghl_neutrino_rate_provider_context transparent_provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_default(&transparent_provider),
        ghl_success, "optional-diagnostics recovery initialization", 2020);
  transparent_provider.failure_policy = ghl_neutrino_rate_failure_transparent;
  ghl_primitive_quantities invalid_prims = prims;
  invalid_prims.rho = NAN;
  ghl_m1_neutrino_rates transparent_rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(transparent_rates);
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &transparent_provider, NULL, NULL, NULL, &invalid_prims,
              transparent_rates),
        ghl_success, "optional-diagnostics transparent recovery", 2020);
  validate_rate_bundle(transparent_rates, 2020);
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    require_condition(
          transparent_rates[species].eta_N == 0.0
                && transparent_rates[species].eta_E == 0.0
                && transparent_rates[species].kappa_a_N == 0.0
                && transparent_rates[species].kappa_a_E == 0.0,
          "transparent recovery published interaction rates", 2020);
  }

  /* Finite-valued reference inputs can still overflow an aggregate product
   * after thermodynamic reconstruction.  This is a supported caller path:
   * hold-last must fail closed when the failed key is not the cached key. */
  ghl_neutrino_rate_provider_context hold_provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_default(&hold_provider), ghl_success,
        "post-thermo hold-last initialization", 2021);
  hold_provider.failure_policy = ghl_neutrino_rate_failure_hold_last;
  ghl_neutrino_rate_provider_cache hold_cache;
  ghl_neutrino_rate_provider_cache_initialize(&hold_cache);
  ghl_neutrino_rate_provider_diagnostics hold_diagnostics = { 0 };
  ghl_m1_neutrino_rates hold_rates[ghl_m1_neutrino_species_count];
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &hold_provider, &hold_cache, &hold_diagnostics, NULL, &prims, hold_rates),
        ghl_success, "post-thermo hold-last seed", 2021);
  const ghl_neutrino_rate_provider_cache hold_cache_before = hold_cache;
  ghl_primitive_quantities overflow_prims = prims;
  overflow_prims.rho = DBL_MAX;
  initialize_sentinel_rates(hold_rates);
  const ghl_m1_neutrino_rates hold_rates_before_failure[ghl_m1_neutrino_species_count]
        = { hold_rates[0], hold_rates[1], hold_rates[2] };
  const ghl_error_codes_t hold_error = ghl_neutrino_rate_provider_compute_cell(
        &hold_provider, &hold_cache, &hold_diagnostics, NULL, &overflow_prims,
        hold_rates);
  require_error(
        hold_error, ghl_error_m1_microphysics_failure,
        "post-thermo hold-last key mismatch", 2022);
  require_condition(
        memcmp(&hold_cache, &hold_cache_before, sizeof(hold_cache)) == 0
              && same_rate_bundle(hold_rates, hold_rates_before_failure),
        "post-thermo failure changed transactional outputs", 2022);
  require_condition(
        hold_diagnostics.failures == 1 && hold_diagnostics.last_error == hold_error
              && hold_diagnostics.last_recovery == ghl_neutrino_rate_recovery_none
              && hold_diagnostics.hold_last_recoveries == 0,
        "post-thermo hold-last diagnostics are incomplete", 2022);

  /* The reference provider also has a supported pre-cache failure: a finite
   * DBL_MAX temperature makes its finite composition chemical potentials
   * overflow, after recovered_key_valid is set.  Keep a matching rates key
   * while invalidating only thermo_valid to exercise successful hold-last
   * recovery before the ordinary cache-hit check. */
  ghl_neutrino_rate_provider_context reference_hold_provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_default(&reference_hold_provider),
        ghl_success, "reference pre-cache hold-last initialization", 2026);
  reference_hold_provider.failure_policy = ghl_neutrino_rate_failure_hold_last;
  ghl_neutrino_rate_provider_cache reference_hold_cache;
  ghl_neutrino_rate_provider_cache_initialize(&reference_hold_cache);
  ghl_m1_neutrino_rates reference_hold_seed[ghl_m1_neutrino_species_count];
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &reference_hold_provider, &reference_hold_cache, NULL, NULL, &prims,
              reference_hold_seed),
        ghl_success, "reference pre-cache hold-last seed", 2026);
  reference_hold_cache.thermo_valid = false;
  reference_hold_cache.T = DBL_MAX;
  ghl_primitive_quantities reference_overflow_prims = prims;
  reference_overflow_prims.temperature = DBL_MAX;
  reference_overflow_prims.Y_e = 0.1;
  reference_hold_cache.Ye = reference_overflow_prims.Y_e;
  const ghl_neutrino_rate_provider_cache reference_hold_expected = reference_hold_cache;
  ghl_neutrino_rate_provider_diagnostics reference_hold_diagnostics = { 0 };
  ghl_m1_neutrino_rates reference_hold_rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(reference_hold_rates);
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &reference_hold_provider, &reference_hold_cache,
              &reference_hold_diagnostics, NULL, &reference_overflow_prims,
              reference_hold_rates),
        ghl_success, "reference pre-cache hold-last recovery", 2026);
  require_condition(
        same_rate_bundle(reference_hold_rates, reference_hold_cache.rates)
              && reference_hold_diagnostics.failures == 1
              && reference_hold_diagnostics.last_error
                       == ghl_error_m1_microphysics_failure
              && reference_hold_diagnostics.hold_last_recoveries == 1
              && reference_hold_diagnostics.last_recovery
                       == ghl_neutrino_rate_recovery_hold_last,
        "reference pre-cache hold-last recovery was not diagnosed", 2026);
  require_condition(
        memcmp(
              &reference_hold_cache, &reference_hold_expected,
              sizeof(reference_hold_cache))
              == 0,
        "reference pre-cache hold-last changed cache metadata", 2026);
  initialize_sentinel_rates(reference_hold_rates);
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &reference_hold_provider, &reference_hold_cache, NULL, NULL,
              &reference_overflow_prims, reference_hold_rates),
        ghl_success, "reference pre-cache hold-last without diagnostics", 2027);
  require_condition(
        same_rate_bundle(reference_hold_rates, reference_hold_cache.rates),
        "reference pre-cache hold-last without diagnostics changed rates", 2027);

  /* Equilibrium recovery builds finite intermediate fields here, but the
   * finite recovery rate times finite J_eq is not representable.  The
   * validation gate must reject publication transactionally and return the
   * original microphysics failure. */
  ghl_neutrino_rate_provider_context overflowing_recovery;
  require_error(
        ghl_neutrino_rate_provider_initialize_default(&overflowing_recovery),
        ghl_success, "overflowing recovery initialization", 2023);
  overflowing_recovery.failure_policy = ghl_neutrino_rate_failure_equilibrium;
  overflowing_recovery.min_mean_energy = DBL_MAX;
  overflowing_recovery.equilibrium_recovery_rate = DBL_MAX;
  ghl_neutrino_rate_provider_cache overflow_cache;
  ghl_neutrino_rate_provider_cache_initialize(&overflow_cache);
  const ghl_neutrino_rate_provider_cache overflow_cache_before = overflow_cache;
  ghl_m1_neutrino_rates overflow_rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(overflow_rates);
  const ghl_m1_neutrino_rates overflow_rates_before[ghl_m1_neutrino_species_count]
        = { overflow_rates[0], overflow_rates[1], overflow_rates[2] };
  ghl_neutrino_rate_provider_diagnostics overflow_diagnostics = { 0 };
  const ghl_error_codes_t overflow_error = ghl_neutrino_rate_provider_compute_cell(
        &overflowing_recovery, &overflow_cache, &overflow_diagnostics, NULL,
        &invalid_prims, overflow_rates);
  require_error(
        overflow_error, ghl_error_m1_microphysics_failure,
        "equilibrium recovery publication overflow", 2023);
  require_condition(
        memcmp(&overflow_cache, &overflow_cache_before, sizeof(overflow_cache)) == 0
              && same_rate_bundle(overflow_rates, overflow_rates_before),
        "overflowing equilibrium recovery changed outputs", 2023);
  require_condition(
        overflow_diagnostics.failures == 1
              && overflow_diagnostics.last_error == overflow_error
              && overflow_diagnostics.equilibrium_recoveries == 0
              && overflow_diagnostics.last_recovery == ghl_neutrino_rate_recovery_none,
        "overflowing equilibrium recovery diagnostics are incomplete", 2023);
}

#ifndef GHL_DISABLE_HDF5
static void test_reference_table_eos_validation(
      const ghl_eos_parameters *restrict eos,
      const ghl_primitive_quantities *restrict prims) {
  ghl_neutrino_rate_provider_context provider;
  require_error(
        ghl_neutrino_rate_provider_initialize_default(&provider), ghl_success,
        "reference table provider initialization", 3600);
  provider.use_tabulated_eos = true;

  static const char *const labels[] = { "reference table hybrid EOS",
                                        "reference table unknown type" };
  const int variant_count = (int)(sizeof(labels) / sizeof(labels[0]));
  for(int variant = 0; variant < variant_count; ++variant) {
    ghl_eos_parameters bad_eos = *eos;
    if(variant == 0) {
      bad_eos.eos_type = ghl_eos_hybrid;
    }
    else {
      bad_eos.table_type = ghl_eos_table_unknown;
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
          &provider, &cache, &diagnostics, &bad_eos, prims, rates);
    require_error(
          error, ghl_error_m1_microphysics_failure, labels[variant], 3601 + variant);
    require_condition(
          memcmp(&cache, &cache_before, sizeof(cache)) == 0
                && same_rate_bundle(rates, rates_before) && diagnostics.failures == 1
                && diagnostics.last_error == error,
          "reference table EOS rejection changed transactional outputs", 3601 + variant);
  }

  ghl_neutrino_rate_provider_cache cache;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
  ghl_m1_neutrino_rates first[ghl_m1_neutrino_species_count];
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, eos, prims, first),
        ghl_success, "valid reference table provider call", 3603);
  validate_rate_bundle(first, 3603);
  ghl_m1_neutrino_rates second[ghl_m1_neutrino_species_count];
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, eos, prims, second),
        ghl_success, "reference table provider cache hit", 3604);
  require_condition(
        diagnostics.cache_hits == 1 && same_rate_bundle(first, second),
        "reference table cache hit changed rates", 3604);
}

static void test_production_provider_context_validation(
      const ghl_neutrino_rate_provider_context *restrict baseline,
      const ghl_eos_parameters *restrict eos,
      const ghl_primitive_quantities *restrict prims) {
  /* Each mutation is a supported provider-boundary input, so a production
   * context must reject it before table lookup or raw-rate assembly.  Keeping
   * one mutation per call also identifies every disjunct in the production
   * unit's conversion and EOS contract. */
  static const char *const labels[] = { "production NaN minimum energy",
                                        "production nonzero minimum energy",
                                        "production reference backend without table EOS",
                                        "production NULL EOS",
                                        "production non-tabulated EOS type",
                                        "production non-StellarCollapse table",
                                        "production density conversion",
                                        "production temperature conversion",
                                        "production opacity conversion",
                                        "production emissivity conversion",
                                        "production baryon-mass conversion",
                                        "production charged-current scale",
                                        "production scattering scale",
                                        "production pair scale",
                                        "production bremsstrahlung scale",
                                        "production plasmon scale" };
  const int variant_count = (int)(sizeof(labels) / sizeof(labels[0]));
  for(int variant = 0; variant < variant_count; ++variant) {
    ghl_neutrino_rate_provider_context provider = *baseline;
    ghl_eos_parameters bad_eos = *eos;
    const ghl_eos_parameters *eos_argument = &bad_eos;
    switch(variant) {
      case 0:
        provider.min_mean_energy = NAN;
        break;
      case 1:
        provider.min_mean_energy = 1.0;
        break;
      case 2:
        provider.use_tabulated_eos = false;
        break;
      case 3:
        eos_argument = NULL;
        break;
      case 4:
        bad_eos.eos_type = ghl_eos_hybrid;
        break;
      case 5:
        bad_eos.table_type = ghl_eos_table_unknown;
        break;
      case 6:
        provider.rho_code_to_cgs *= 2.0;
        break;
      case 7:
        provider.temperature_code_to_mev = 2.0;
        break;
      case 8:
        provider.opacity_cgs_to_code *= 2.0;
        break;
      case 9:
        provider.emissivity_cgs_to_code *= 2.0;
        break;
      case 10:
        provider.baryon_mass_code *= 2.0;
        break;
      case 11:
        provider.charged_current_scale = 2.0;
        break;
      case 12:
        provider.scattering_scale = 2.0;
        break;
      case 13:
        provider.pair_scale = 2.0;
        break;
      case 14:
        provider.bremsstrahlung_scale = 2.0;
        break;
      case 15:
        provider.plasmon_scale = 2.0;
        break;
      default:
        ghl_error("M1 rate-provider production-context test has an unknown variant\n");
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
    require_condition(
          provider->min_mean_energy == 0.0,
          "production initializer retained a physical energy floor", case_index);

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
  test_production_representability_transaction(&provider);
}

static void test_table_provider(const char *restrict table_path) {
  ghl_eos_parameters eos = { 0 };
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
  require_condition(
        provider.min_mean_energy == 0.0,
        "production initializer selected a physical energy floor", 3000);
  ghl_primitive_quantities context_prims = { 0 };
  context_prims.rho = sqrt(eos.table_rho_min * eos.table_rho_max);
  context_prims.temperature = sqrt(eos.table_T_min * eos.table_T_max);
  context_prims.Y_e = 0.5 * (eos.table_Y_e_min + eos.table_Y_e_max);
  context_prims.eps = 1.0;
  test_reference_table_eos_validation(&eos, &context_prims);
  test_production_provider_context_validation(&provider, &eos, &context_prims);
  ghl_neutrino_rate_provider_cache cache;
  ghl_neutrino_rate_provider_cache_initialize(&cache);
  ghl_neutrino_rate_provider_diagnostics diagnostics = { 0 };
  ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
  provider_rng table_rng = { .state = UINT64_C(0x4d315f5441424c45) };

  for(int case_index = 0; case_index < TABLE_RANDOM_CASES; ++case_index) {
    const double rho_fraction = provider_rng_between(&table_rng, 0.05, 0.95);
    const double temperature_fraction = provider_rng_between(&table_rng, 0.05, 0.95);
    const double ye_fraction = provider_rng_between(&table_rng, 0.05, 0.95);
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
    const ghl_error_codes_t error = ghl_neutrino_rate_provider_compute_cell(
          &provider, &cache, &diagnostics, &eos, &prims, rates);
    require_error(error, ghl_success, "random table provider call", 3001 + case_index);
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
        = { NRPyEOS_mu_e_key, NRPyEOS_X_n_key, NRPyEOS_X_p_key };
  const double malformed_table_values[] = { NAN, -1.0e-6, 1.0 + 1.0e-6 };
  const char *const malformed_table_names[]
        = { "table NaN chemical potential", "table negative neutron fraction",
            "table super-unit proton fraction" };
  for(int malformed = 0; malformed < 3; ++malformed) {
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
  double seeded_beta_kirchhoff_relative_mismatch[ghl_m1_neutrino_species_count]
        = { 0.0 };
  bool seeded_beta_kirchhoff_mismatch_valid[ghl_m1_neutrino_species_count] = { false };
  const ghl_neutrino_rate_failure_policy_t saved_table_failure_policy
        = provider.failure_policy;
  /* Seed the rates key under the same policy used by the pre-cache witness;
   * changing provider configuration after seeding would intentionally make
   * same_rates_key reject the cache as stale. */
  provider.failure_policy = ghl_neutrino_rate_failure_hold_last;
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, &eos, &cached_prims, first),
        ghl_success, "table cache seed", 3040);
  memcpy(
        seeded_beta_kirchhoff_relative_mismatch,
        diagnostics.beta_kirchhoff_relative_mismatch,
        sizeof(seeded_beta_kirchhoff_relative_mismatch));
  memcpy(
        seeded_beta_kirchhoff_mismatch_valid, diagnostics.beta_kirchhoff_mismatch_valid,
        sizeof(seeded_beta_kirchhoff_mismatch_valid));
  require_condition(
        memcmp(
              cache.beta_kirchhoff_relative_mismatch,
              seeded_beta_kirchhoff_relative_mismatch,
              sizeof(seeded_beta_kirchhoff_relative_mismatch))
                    == 0
              && memcmp(
                       cache.beta_kirchhoff_mismatch_valid,
                       seeded_beta_kirchhoff_mismatch_valid,
                       sizeof(seeded_beta_kirchhoff_mismatch_valid))
                       == 0,
        "table cache seed did not retain production diagnostics", 3040);

  /* A failure while refreshing thermodynamics occurs before the ordinary
   * same-rates cache hit.  Keep a valid rates key but invalidate only the
   * thermo key, then make the public table interpolator return NaN.  This is
   * the supported pre-cache route to hold-last recovery; both diagnostics
   * conventions must preserve the cached bundle and cache metadata. */
  const ghl_neutrino_rate_provider_cache hold_cache_before = cache;
  cache.thermo_valid = false;
  double saved_hold_mu_e_corners[8];
  provider_set_table_corners(&eos, NRPyEOS_mu_e_key, NAN, saved_hold_mu_e_corners);
  ghl_neutrino_rate_provider_diagnostics hold_diagnostics = { 0 };
  ghl_m1_neutrino_rates held_rates[ghl_m1_neutrino_species_count];
  initialize_sentinel_rates(held_rates);
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &hold_diagnostics, &eos, &cached_prims, held_rates),
        ghl_success, "pre-cache hold-last recovery", 3160);
  require_condition(
        same_rate_bundle(held_rates, cache.rates) && hold_diagnostics.failures == 1
              && hold_diagnostics.last_error == ghl_error_m1_microphysics_failure
              && hold_diagnostics.hold_last_recoveries == 1
              && hold_diagnostics.last_recovery == ghl_neutrino_rate_recovery_hold_last
              && memcmp(
                       hold_diagnostics.beta_kirchhoff_relative_mismatch,
                       cache.beta_kirchhoff_relative_mismatch,
                       sizeof(hold_diagnostics.beta_kirchhoff_relative_mismatch))
                       == 0
              && memcmp(
                       hold_diagnostics.beta_kirchhoff_relative_mismatch,
                       seeded_beta_kirchhoff_relative_mismatch,
                       sizeof(seeded_beta_kirchhoff_relative_mismatch))
                       == 0
              && memcmp(
                       hold_diagnostics.beta_kirchhoff_mismatch_valid,
                       cache.beta_kirchhoff_mismatch_valid,
                       sizeof(hold_diagnostics.beta_kirchhoff_mismatch_valid))
                       == 0
              && memcmp(
                       hold_diagnostics.beta_kirchhoff_mismatch_valid,
                       seeded_beta_kirchhoff_mismatch_valid,
                       sizeof(seeded_beta_kirchhoff_mismatch_valid))
                       == 0,
        "pre-cache hold-last recovery was not diagnosed", 3160);
  require_condition(
        memcmp(&cache, &hold_cache_before, sizeof(cache)) != 0
              && cache.thermo_valid == false && cache.rates_valid
              && same_rate_bundle(cache.rates, first),
        "pre-cache hold-last did not retain the valid rates key", 3160);
  /* The only intended difference from the pre-call snapshot is the explicit
   * caller invalidation of thermo_valid; provider recovery itself is inert. */
  ghl_neutrino_rate_provider_cache hold_cache_expected = hold_cache_before;
  hold_cache_expected.thermo_valid = false;
  require_condition(
        memcmp(&cache, &hold_cache_expected, sizeof(cache)) == 0,
        "pre-cache hold-last changed cache metadata", 3160);

  initialize_sentinel_rates(held_rates);
  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, NULL, &eos, &cached_prims, held_rates),
        ghl_success, "pre-cache hold-last without diagnostics", 3161);
  require_condition(
        same_rate_bundle(held_rates, cache.rates) && same_rate_bundle(held_rates, first),
        "pre-cache hold-last without diagnostics changed rates", 3161);
  provider_restore_table_corners(&eos, NRPyEOS_mu_e_key, saved_hold_mu_e_corners);

  require_error(
        ghl_neutrino_rate_provider_compute_cell(
              &provider, &cache, &diagnostics, &eos, &cached_prims, second),
        ghl_success, "table cache hit", 3041);
  require_condition(
        diagnostics.cache_hits == 1 && same_rate_bundle(first, second),
        "table cache hit was not exact", 3041);
  provider.failure_policy = saved_table_failure_policy;

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
  require_condition(
        cache.thermo_valid && isfinite(cache.thermo_T) && cache.thermo_T > 0.0,
        "table temperature recovery did not publish thermodynamics", 3044);

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
  fallback_bounds.table_Y_e_min = 0.0;
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

  ghl_tabulated_free_memory(&eos);
}
#endif

int main(int argc, char **argv) {
  if(argc > 2) {
    ghl_error(
          "Usage: %s [StellarCollapse EOS table path|--generated-fixture]\n", argv[0]);
  }

  provider_rng rng = { .state = UINT64_C(0x4d3150524f564944) };
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

#ifndef GHL_DISABLE_HDF5
  test_production_provider_regressions();
  if(argc == 2) {
    char generated_fixture_path[128] = { 0 };
    const char *table_path = argv[1];
    if(strcmp(argv[1], "--generated-fixture") == 0) {
      if(!create_provider_fixture(
               generated_fixture_path, sizeof(generated_fixture_path))) {
        ghl_error("Could not create the generated provider EOS fixture\n");
      }
      table_path = generated_fixture_path;
    }
    test_table_provider(table_path);
    if(generated_fixture_path[0] != '\0') {
      remove(generated_fixture_path);
    }
  }
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
