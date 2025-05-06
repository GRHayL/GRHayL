#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_FILENAME_SIZE (1024)

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

static size_t read_size_from_file(const char *basename, const char *ext) {
  char filename[MAX_FILENAME_SIZE];
  concatenate(basename, ext, filename);

  FILE *fp = fopen_or_exit(filename, "r");

  size_t size = 0;
  if(fscanf(fp, "%*lu\n%lu\n", &size) != 1) {
    fprintf(stderr, "Error: could not read size from file '%s'\n", filename);
    fclose(fp);
    exit(1);
  }

  fclose(fp);
  return size;
}

// Read file with n columns of doubles
static void read_file_with_n_cols(
      const char *filename,
      const size_t n_lines_expected,
      const size_t n_lines_skip,
      const size_t n_cols,
      double **data_out) {
  FILE *fp = fopen_or_exit(filename, "r");

  char *line = NULL;
  size_t line_size = 0;
  for(size_t n = 0; n < n_lines_skip; n++) {
    if(getline(&line, &line_size, fp) == -1) {
      fprintf(stderr, "Error: Not enough lines to skip in file '%s'\n", filename);
      free(line);
      fclose(fp);
      exit(1);
    }
  }

  double *data = malloc_or_exit(sizeof(double) * n_cols * n_lines_expected);

  size_t n_lines = 0;
  while(getline(&line, &line_size, fp) > 0) {
    if(n_lines >= n_lines_expected) {
      fprintf(stderr, "Error: More lines in file than expected\n");
      free(data);
      free(line);
      fclose(fp);
      exit(1);
    }

    const size_t offset = n_lines * n_cols;

    char *token = strtok(line, " ");
    data[0 + offset] = strtod(token, NULL);

    for(size_t col = 1; col < n_cols; col++) {
      if(!token) {
        fprintf(stderr, "Error: Missing columns in line %lu\n", n_lines + 1);
        free(data);
        free(line);
        fclose(fp);
        exit(1);
      }
      char *token = strtok(NULL, " ");
      data[col + n_cols * n_lines] = strtod(token, NULL);
    }
  }

  *data_out = data;
  free(line);
  fclose(fp);
}

int main(int argc, char **argv) {
  if(argc != 2) {
    fprintf(stderr, "Usage: %s <eos basename with path, e.g., /my/path/eos>\n", argv[0]);
    return 1;
  }
  const char *basename = argv[1];

  const size_t n_rho = read_size_from_file(basename, ".nb");
  const size_t n_temp = read_size_from_file(basename, ".t");
  const size_t n_ye = read_size_from_file(basename, ".yq");

  printf("(%lu, %lu, %lu)\n", n_rho, n_temp, n_ye);

  return 0;
}
