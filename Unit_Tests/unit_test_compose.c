#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ghl.h"
#include "nrpyeos_tabulated.h"

#define MAX_FILENAME_SIZE            (1024)

#define MEV_OVER_C2_TO_GRAMS         (1.782661921627898e-27)
#define FM_TO_CM                     (1.000000000000000e-13)
#define MEV_OVER_FM3_TO_DYN_OVER_CM2 (1.602176633999999e+33)
#define MEV_OVER_C2_OVER_FM3_TO_GPCC (1.782661921627897e+12)
#define C2_TO_ERG_OVER_GRAMS         (8.987551787368175e+20)

#define MSUN_IN_MEV_OVER_C2          (1.115416135036007e+60)
#define MEV_OVER_C2_TO_CODE          (1.0 / MSUN_IN_MEV_OVER_C2)
#define FM3_TO_CODE                  (3.219664984827565e+54)
#define INV_FM3_TO_CODE              (1.0 / FM3_TO_CODE)
#define MEV_OVER_FM3_TO_CODE         (MEV_OVER_FM3_TO_DYN_OVER_CM2 * CGS_TO_CODE_PRESSURE)

/**
 * @brief Structure to hold CompOSE Equation of State (EOS) data.
 *
 * This structure stores information read from a CompOSE EOS file,
 * including metadata and thermodynamic quantities.
 * The thermodynamic quantities stored in the Q array are:
 *   - Q[0]: p/n_b               -> Pressure/n_b [MeV]
 *   - Q[1]: s/n_b               -> Entropy per baryon/n_b
 *   - Q[2]: mu_b/m_n - 1        -> Baryon chemical potential/m_n - 1
 *   - Q[3]: mu_q/m_n            -> Charge chemical potential/m_n
 *   - Q[4]: mu_l/m_n            -> Lepton chemical potential/m_n
 *   - Q[5]: f/(n_b * m_n) - 1   -> Helmholtz free energy/(n_b * m_n) - 1
 *   - Q[6]: eps                 -> Specific internal energy
 * where n_b is the baryon number density and m_n is the neutron mass. The only
 * dimensionful quantity is Q[0], which is given in MeV.
 */
typedef struct {
  char filename[MAX_FILENAME_SIZE]; ///< Aux: holds the name of the EOS file we are working on.
  int size[3];                      ///< Dimensions:  [0] n_b, [1] temperature, [2] charge fraction.
  double *grid[3];                  ///< Discretized: [0] n_b, [1] temperature, [2] charge fraction.
  double *Q[7];                     ///< Thermodynamic quantities described above.
  double m_n, m_p;                  ///< Neutron and proton masses.
  int contains_leptons;             ///< Whether or not the table contains leptons (0 if not).
} compose_eos_table;

static void *malloc_or_exit(const int size) {
  void *ptr = malloc(size);
  if(!ptr) {
    fprintf(stderr, "Error: could not allocate '%lu' bytes\n", size);
    exit(1);
  }
  return ptr;
}

static FILE *fopen_or_exit(const char *filename, const char *mode) {
  FILE *fp = fopen(filename, mode);
  if(!fp) {
    fprintf(stderr, "Error: could not open file '%s'\n", filename);
    exit(1);
  }
  return fp;
}

static void concatenate(const char *str1, const char *str2, char *out) {
  snprintf(out, MAX_FILENAME_SIZE, "%s%s", str1, str2);
}

static void compose_eos_table_free(compose_eos_table *eos) {
  if(eos == NULL) {
    return;
  }

  for(int n = 0; n < 3; n++) {
    if(eos->grid[n] != NULL) {
      free(eos->grid[n]);
    }
  }

  for(int n = 0; n < 7; n++) {
    if(eos->Q[n] != NULL) {
      free(eos->Q[n]);
    }
  }

  free(eos);
}

static void compose_eos_read_nb_t_ye(const char *basename, compose_eos_table *table) {
  const char *exts[3] = { ".nb", ".t", ".yq" };

  for(int n = 0; n < 3; n++) {
    concatenate(basename, exts[n], table->filename);
    FILE *fp = fopen_or_exit(table->filename, "r");

    if(fscanf(fp, "%*u\n%lu\n", &table->size[n]) != 1) {
      fprintf(stderr, "Error: could not read size from file '%s'\n", table->filename);
      fclose(fp);
      exit(1);
    }
    table->grid[n] = (double *)malloc_or_exit(sizeof(double) * table->size[n]);

    int count = 0;
    double value = 0;
    while(count < table->size[n]) {
      if(fscanf(fp, "%lf%*[^\n]", &value) != 1) {
        fprintf(stderr, "Error: Failed to parse line %lu in '%s'\n", count + 3, table->filename);
        compose_eos_table_free(table);
        fclose(fp);
        exit(1);
      }
      table->grid[n][count++] = value;
    }
    fclose(fp);

    printf("Finished reading file '%s'\n", table->filename);
  }
}

static int index_cp(const compose_eos_table *table, const int i_n, const int i_t, const int i_y) {
  return i_n + table->size[0] * (i_t + table->size[1] * i_y);
}

static int index_sc(const ghl_eos_parameters *table, const int i_r, const int i_t, const int i_y) {
  return i_r + table->N_rho * (i_t + table->N_T * i_y);
}

// Read file with n columns of doubles
static void compose_eos_read_thermo(const char *basename, compose_eos_table *table) {
  concatenate(basename, ".thermo", table->filename);
  FILE *fp = fopen_or_exit(table->filename, "r");

  if(fscanf(fp, "%lf %lf %d\n", &table->m_n, &table->m_p, &table->contains_leptons) != 3) {
    fprintf(stderr, "Error: Failed to read header from '%s'\n", table->filename);
    fclose(fp);
    compose_eos_table_free(table);
    exit(1);
  }

  const int size = table->size[0] * table->size[1] * table->size[2];
  for(int n = 0; n < 7; n++) {
    table->Q[n] = (double *)malloc_or_exit(sizeof(double) * size);
  }

  int n_lines = 0;
  while(n_lines < size) {

    int i_t, i_n, i_y;
    double Q[7];
    // clang-format off
    int read = fscanf(
          fp, "%d %d %d %lf %lf %lf %lf %lf %lf %lf%*[^\n]",
          &i_t, &i_n, &i_y,
          &Q[0], &Q[1], &Q[2], &Q[3], &Q[4], &Q[5], &Q[6]);
    // clang-format on
    if(read != 10) {
      fprintf(stderr, "Error: Failed to parse line %lu from '%s'\n", n_lines + 1, table->filename);
      fclose(fp);
      compose_eos_table_free(table);
      exit(1);
    }

    // CompOSE table indices adopt Fortran standards; fix it here.
    const int index = index_cp(table, --i_n, --i_t, --i_y);
    for(int n = 0; n < 7; n++) {
      table->Q[n][index] = Q[n];
    }

    n_lines++;
  }
  fclose(fp);

  printf("Finished reading file '%s'\n", table->filename);
}

static void compute_derivatives(ghl_eos_parameters *ghl_eos) {

  const double dlt = ghl_eos->table_logT[1] - ghl_eos->table_logT[0];

  // Shifts for finite differences points
  const int fd_shifts[2][2] = {
    { +0, +1}, // forward
    { -1, +0}  // backward
  };

  for(int iy = 0; iy < ghl_eos->N_Ye; iy++) {
    for(int it = 0; it < ghl_eos->N_T; it++) {
      const double temp = exp(ghl_eos->table_logT[it]);
      for(int ir = 0; ir < ghl_eos->N_rho; ir++) {

        const int fd_type = (it == ghl_eos->N_T - 1);

        const int it_l = it + fd_shifts[fd_type][0];
        const int it_r = it + fd_shifts[fd_type][1];

        const int index_tl = NRPyEOS_ntablekeys * index_sc(ghl_eos, ir, it_l, iy);
        const int index_tr = NRPyEOS_ntablekeys * index_sc(ghl_eos, ir, it_r, iy);

        const double le_tl = ghl_eos->table_all[index_tl + NRPyEOS_eps_key];
        const double le_tr = ghl_eos->table_all[index_tr + NRPyEOS_eps_key];

        // Compute dlog(eps + eps_0)
        const double dle = le_tr - le_tl;

        // Note the following identity:
        //
        // deps/dT = (eps + eps_0)/T * (dlog(eps + eps_0)/dlog(T))
        const int index = index_sc(ghl_eos, ir, it, iy);
        const double eps_shifted = ghl_eos->table_eps[index] + ghl_eos->energy_shift;
        const double dedt = (eps_shifted / temp) * (dle / dlt);

        // ghl_eos->table_all[index + NRPyEOS_dPdrho_key] = dpdrhoe;
        // ghl_eos->table_all[index + NRPyEOS_dPdeps_key] = dpderho;
        ghl_eos->table_all[NRPyEOS_ntablekeys * index + NRPyEOS_depsdT_key] = dedt;
      }
    }
  }
}

static void
compose_table_to_ghl_eos(const compose_eos_table *compose_table, ghl_eos_parameters *ghl_eos) {
  ghl_eos->N_rho = compose_table->size[0];
  ghl_eos->N_T = compose_table->size[1];
  ghl_eos->N_Ye = compose_table->size[2];

  // Default energy shift of 20 MeV = 1.91313e19 erg / g
  ghl_eos->energy_shift = 2.128643702854307e-02;

  // Log of baryonic density, converting from n_b to geometrized units
  const double n_b_to_rho = compose_table->m_n * MEV_OVER_C2_OVER_FM3_TO_GPCC * CGS_TO_CODE_DENSITY;
  ghl_eos->table_logrho = (double *restrict)malloc_or_exit(sizeof(double) * ghl_eos->N_rho);
  for(int ir = 0; ir < ghl_eos->N_rho; ir++) {
    ghl_eos->table_logrho[ir] = log(n_b_to_rho * compose_table->grid[0][ir]);
  }

  // Temperature, converting to log-scale
  ghl_eos->table_logT = (double *restrict)malloc_or_exit(sizeof(double) * ghl_eos->N_T);
  for(int it = 0; it < ghl_eos->N_T; it++) {
    ghl_eos->table_logT[it] = log(compose_table->grid[1][it]);
  }

  // Electron fraction, simple copy
  ghl_eos->table_Y_e = (double *restrict)malloc_or_exit(sizeof(double) * ghl_eos->N_Ye);
  for(int iy = 0; iy < ghl_eos->N_Ye; iy++) {
    ghl_eos->table_Y_e[iy] = compose_table->grid[2][iy];
  }

  // Neutron mass, in MeV
  const double m_n = compose_table->m_n;

  const int size = ghl_eos->N_rho * ghl_eos->N_T * ghl_eos->N_Ye;
  const int size_all = size * NRPyEOS_ntablekeys;
  ghl_eos->table_eps = (double *restrict)malloc_or_exit(sizeof(double) * size);
  ghl_eos->table_all = (double *restrict)malloc_or_exit(sizeof(double) * size_all);
  for(int iy = 0; iy < ghl_eos->N_Ye; iy++) {
    for(int it = 0; it < ghl_eos->N_T; it++) {
      for(int ir = 0; ir < ghl_eos->N_rho; ir++) {
        const int compose_index = index_cp(compose_table, ir, it, iy);
        const int ghl_eos_index = ir + ghl_eos->N_rho * (it + ghl_eos->N_T * iy);

        // Baryon number density, 1/fm^3
        const double n_b = compose_table->grid[0][ir];

        // Density, dimensionless
        const double rho = exp(ghl_eos->table_logrho[ir]);

        // Q[0] = P / n_b, in MeV
        const double press = compose_table->Q[0][compose_index] * n_b * MEV_OVER_FM3_TO_CODE;

        // Q[1] = s / n_b, dimensionless
        const double entropy = compose_table->Q[1][compose_index];

        // Q[2] = mu_b/m_n - 1, dimensionless. However, what we want is mu_b w.r.t. m_n,
        // so we must convert it using:
        // mu_b = m_n * (Q[2] + 1); mu_b -> mu_b - m_n = m_n * Q[2], which is now in MeV.
        const double mu_b = m_n * compose_table->Q[2][compose_index];

        // Q[3] = mu_q/m_n, dimensionless, so mu_q in MeV
        const double mu_q = m_n * compose_table->Q[3][compose_index];

        // Q[4] = mu_l/m_n, dimensionless, so mu_l in MeV
        const double mu_l = m_n * compose_table->Q[4][compose_index];

        // Q[5] = f/(n_b * m_n) - 1, dimensionless, so f dimensionless
        const double f = (compose_table->Q[5][compose_index] + 1) * rho;

        // Q[6] = eps, dimensionless
        const double eps = compose_table->Q[6][compose_index];
        const double eps_shifted = eps + ghl_eos->energy_shift;

        // Compute the chemical potentials
        const double mu_n = mu_b;
        const double mu_p = mu_b + mu_q;
        const double mu_e = mu_l - mu_q;
        const double muhat = mu_n - mu_p;
        const double munu = mu_e - mu_n + mu_p;

        // Write to the eos struct
        ghl_eos->table_eps[ghl_eos_index] = eps;
        ghl_eos->table_all[ghl_eos_index * NRPyEOS_ntablekeys + NRPyEOS_press_key] = log(press);
        ghl_eos->table_all[ghl_eos_index * NRPyEOS_ntablekeys + NRPyEOS_eps_key] = log(eps_shifted);
        ghl_eos->table_all[ghl_eos_index * NRPyEOS_ntablekeys + NRPyEOS_entropy_key] = entropy;
        ghl_eos->table_all[ghl_eos_index * NRPyEOS_ntablekeys + NRPyEOS_mu_n_key] = mu_n;
        ghl_eos->table_all[ghl_eos_index * NRPyEOS_ntablekeys + NRPyEOS_mu_p_key] = mu_p;
        ghl_eos->table_all[ghl_eos_index * NRPyEOS_ntablekeys + NRPyEOS_mu_e_key] = mu_e;
        ghl_eos->table_all[ghl_eos_index * NRPyEOS_ntablekeys + NRPyEOS_muhat_key] = muhat;
        ghl_eos->table_all[ghl_eos_index * NRPyEOS_ntablekeys + NRPyEOS_munu_key] = munu;
      }
    }
  }
  compute_derivatives(ghl_eos);
}

static ghl_eos_parameters *read_compose_eos_table_from_basename(const char *basename) {

  compose_eos_table *compose_table = (compose_eos_table *)malloc_or_exit(sizeof(compose_eos_table));
  memset(compose_table, 0, sizeof(ghl_eos_parameters));

  compose_eos_read_nb_t_ye(basename, compose_table);
  compose_eos_read_thermo(basename, compose_table);

  ghl_eos_parameters *ghl_eos = (ghl_eos_parameters *)malloc_or_exit(sizeof(ghl_eos_parameters));
  memset(ghl_eos, 0, sizeof(ghl_eos_parameters));
  compose_table_to_ghl_eos(compose_table, ghl_eos);

  compose_eos_table_free(compose_table);

  return ghl_eos;
}
static ghl_eos_parameters *read_stellarcollapse_eos_table_from_filepath(const char *filepath) {
  ghl_eos_parameters *eos = (ghl_eos_parameters *)malloc_or_exit(sizeof(ghl_eos_parameters));
  NRPyEOS_read_table_set_EOS_params(filepath, eos);
  return eos;
}

static double relative_error(const double x, const double y) {
  if(x != 0.0) {
    return fabs(1.0 - y / x);
  }
  if(y != 0.0) {
    return fabs(1.0 - x / y);
  }
  return 0.0;
}

#define COMPARE1D(n, qty)                                               \
  for(int i = 0; i < table1->n; i++) {                                  \
    const double qty1 = table1->table_##qty[i];                         \
    const double qty2 = table2->table_##qty[i];                         \
    if(relative_error(qty1, qty2) > 1e-7) {                             \
      fprintf(stderr, "%d, %s: %.15e != %.15e\n", i, #qty, qty1, qty2); \
    }                                                                   \
  }

#define COMPARE3D(index, qty)                                             \
  {                                                                       \
    const int idx = index + NRPyEOS_##qty##_key;                          \
    const double qty1 = table1->table_all[idx];                           \
    const double qty2 = table2->table_all[idx];                           \
    if(relative_error(qty1, qty2) > 1e-7) {                               \
      fprintf(stderr, "%d, %s: %.15e != %.15e\n", idx, #qty, qty1, qty2); \
    }                                                                     \
  }

static void compare_tables(const ghl_eos_parameters *table1, ghl_eos_parameters *table2) {
  assert(table1->N_rho == table2->N_rho);
  assert(table1->N_T == table2->N_T);
  assert(table1->N_Ye == table2->N_Ye);

  COMPARE1D(N_rho, logrho);
  COMPARE1D(N_T, logT);
  COMPARE1D(N_Ye, Y_e);

  for(int idx = 0; idx < table1->N_rho * table1->N_T * table1->N_Ye; idx++) {
    const int index = NRPyEOS_ntablekeys * idx;

    COMPARE3D(index, press);
    COMPARE3D(index, eps);
    COMPARE3D(index, entropy);
    COMPARE3D(index, mu_n);
    COMPARE3D(index, mu_p);
    COMPARE3D(index, mu_e);
    COMPARE3D(index, muhat);
    COMPARE3D(index, munu);
    COMPARE3D(index, depsdT);
  }
}

int main(int argc, char **argv) {
  if(argc != 3) {
    fprintf(stderr, "Usage: %s <CompOSE Table Basename> <StellarCollapse Table>\n", argv[0]);
    return 1;
  }
  ghl_eos_parameters *eos1 = read_compose_eos_table_from_basename(argv[1]);
  ghl_eos_parameters *eos2 = read_stellarcollapse_eos_table_from_filepath(argv[2]);

  compare_tables(eos1, eos2);

  NRPyEOS_free_memory(eos1);
  NRPyEOS_free_memory(eos2);

  return 0;
}
