#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_FILENAME_SIZE (1024)

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
 *   - Q[6]: eps/(n_b * m_n) - 1 -> Specific internal energy / (n_b * m_n) - 1
 * where n_b is the baryon number density and m_n is the neutron mass. The only
 * dimensionful quantity is Q[0], which is given in MeV.
 */
typedef struct {
  char filename[MAX_FILENAME_SIZE]; ///< Aux: holds the name of the EOS file we are working on.
  size_t size[3];                   ///< Dimensions:  [0] n_b, [1] temperature, [2] charge fraction.
  double *grid[3];                  ///< Discretized: [0] n_b, [1] temperature, [2] charge fraction.
  double *Q[7];                     ///< Thermodynamic quantities described above.
  double m_n, m_p;                  ///< Neutron and proton masses.
  int contains_leptons;             ///< Whether or not the table contains leptons (0 if not).
} compose_eos_table;

static void *malloc_or_exit(const size_t size) {
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

  for(size_t n = 0; n < 3; n++) {
    if(eos->grid[n] != NULL) {
      free(eos->grid[n]);
    }
  }

  for(size_t n = 0; n < 7; n++) {
    if(eos->Q[n] != NULL) {
      free(eos->Q[n]);
    }
  }

  free(eos);
}

static void compose_eos_read_nb_t_ye(const char *basename, compose_eos_table *table) {
  const char *exts[3] = { ".nb", ".t", ".yq" };

  for(size_t n = 0; n < 3; n++) {
    concatenate(basename, exts[n], table->filename);
    FILE *fp = fopen_or_exit(table->filename, "r");

    if(fscanf(fp, "%*u\n%lu\n", &table->size[n]) != 1) {
      fprintf(stderr, "Error: could not read size from file '%s'\n", table->filename);
      fclose(fp);
      exit(1);
    }
    table->grid[n] = (double *)malloc_or_exit(sizeof(double) * table->size[n]);

    size_t count = 0;
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

static size_t compose_eos_table_index(
      const compose_eos_table *table,
      const size_t i_n,
      const size_t i_t,
      const size_t i_y) {
  return i_n + table->size[0] * (i_t + table->size[1] * i_y);
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

  const size_t size = table->size[0] * table->size[1] * table->size[2];
  for(size_t n = 0; n < 7; n++) {
    table->Q[n] = (double *)malloc_or_exit(sizeof(double) * size);
  }

  size_t n_lines = 0;
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
    const size_t index = compose_eos_table_index(table, --i_n, --i_t, --i_y);
    for(size_t n = 0; n < 7; n++) {
      table->Q[n][index] = Q[n];
    }

    n_lines++;
  }
  fclose(fp);

  printf("Finished reading file '%s'\n", table->filename);
}

static compose_eos_table *compose_eos_table_from_basename(const char *basename) {
  compose_eos_table *table = (compose_eos_table *)malloc_or_exit(sizeof(compose_eos_table));
  memset(table, 0, sizeof(compose_eos_table));

  compose_eos_read_nb_t_ye(basename, table);
  compose_eos_read_thermo(basename, table);

  return table;
}

int main(int argc, char **argv) {
  if(argc != 2) {
    fprintf(stderr, "Usage: %s <eos basename with path, e.g., /my/path/eos>\n", argv[0]);
    return 1;
  }
  const char *basename = argv[1];
  compose_eos_table *table = compose_eos_table_from_basename(basename);

  compose_eos_table_free(table);

  return 0;
}
