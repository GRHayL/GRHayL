#include "m1_thcm1_fixture_utils.h"

#include "ghl_unit_tests.h"

static void m1_thcm1_fixture_set_errorf(
      char *restrict error,
      const size_t error_size,
      const char *restrict format,
      const char *restrict text,
      const size_t number) {
  if(error != NULL && error_size > 0) {
    snprintf(error, error_size, format, text != NULL ? text : "<unknown>", number);
  }
}

static void m1_thcm1_fixture_free_record(m1_thcm1_fixture_record *restrict record) {
  if(record == NULL) {
    return;
  }
  free(record->case_id);
  free(record->pair_id);
  free(record->origin);
  free(record->seed_id);
  free(record->family);
  free(record->perturbation);
  free(record->baseline_input);
  free(record->perturbed_input);
  free(record->baseline_output);
  free(record->perturbed_output);
  free(record->normalization);
  *record = (m1_thcm1_fixture_record){ 0 };
  record->sensitivity_start = SIZE_MAX;
}

void m1_thcm1_fixture_free(m1_thcm1_fixture_collection *restrict collection) {
  if(collection == NULL) {
    return;
  }
  if(collection->records != NULL) {
    for(size_t i = 0; i < collection->record_count; ++i) {
      m1_thcm1_fixture_free_record(&collection->records[i]);
    }
  }
  free(collection->records);
  free(collection->operation);
  free(collection->policy);
  *collection = (m1_thcm1_fixture_collection){ 0 };
}

static char *m1_thcm1_fixture_read_token(FILE *restrict file) {
  int character = fgetc(file);
  while(character != EOF && isspace((unsigned char)character)) {
    character = fgetc(file);
  }
  if(character == EOF) {
    return NULL;
  }
  size_t count = 0;
  size_t capacity = 16;
  char *token = (char *)malloc(capacity);
  if(token == NULL) {
    return NULL;
  }
  do {
    if(count + 1 >= capacity) {
      if(capacity > SIZE_MAX / 2) {
        free(token);
        return NULL;
      }
      const size_t next_capacity = capacity * 2;
      char *expanded = (char *)realloc(token, next_capacity);
      if(expanded == NULL) {
        free(token);
        return NULL;
      }
      token = expanded;
      capacity = next_capacity;
    }
    token[count++] = (char)character;
    character = fgetc(file);
  } while(character != EOF && !isspace((unsigned char)character));
  token[count] = '\0';
  return token;
}

static int m1_thcm1_fixture_expect(FILE *restrict file, const char *restrict expected) {
  char *token = m1_thcm1_fixture_read_token(file);
  const int matches = token != NULL && strcmp(token, expected) == 0;
  free(token);
  return matches;
}

static int m1_thcm1_fixture_string(FILE *restrict file, char **restrict destination) {
  char *token = m1_thcm1_fixture_read_token(file);
  if(token == NULL || token[0] == '\0') {
    free(token);
    return 0;
  }
  free(*destination);
  *destination = token;
  return 1;
}

static int m1_thcm1_fixture_vector(
      FILE *restrict file,
      const char *restrict key,
      double **restrict destination,
      size_t *restrict count) {
  size_t length = 0;
  if(!m1_thcm1_fixture_expect(file, key) || fscanf(file, "%zu", &length) != 1
     || length == 0) {
    return 0;
  }
  double *values = (double *)calloc(length, sizeof(*values));
  if(values == NULL) {
    return 0;
  }
  for(size_t i = 0; i < length; ++i) {
    if(fscanf(file, "%lf", &values[i]) != 1) {
      free(values);
      return 0;
    }
  }
  free(*destination);
  *destination = values;
  *count = length;
  return 1;
}

int m1_thcm1_fixture_record_valid(
      const m1_thcm1_fixture_record *restrict record,
      const size_t expected_input_count,
      const size_t expected_output_count) {
  if(record == NULL || record->case_id == NULL || record->pair_id == NULL
     || record->origin == NULL || record->seed_id == NULL || record->family == NULL
     || record->perturbation == NULL || record->case_id[0] == '\0'
     || record->pair_id[0] == '\0' || record->origin[0] == '\0'
     || record->seed_id[0] == '\0' || record->family[0] == '\0'
     || record->perturbation[0] == '\0' || record->input_count == 0
     || record->output_count == 0 || record->normalization_count != record->output_count
     || (expected_input_count != 0 && record->input_count != expected_input_count)
     || (expected_output_count != 0 && record->output_count != expected_output_count)
     || record->baseline_available != 1 || record->perturbed_available != 1
     || record->baseline_status != 0 || record->perturbed_status != 0
     || record->baseline_reference_status != 0 || record->perturbed_reference_status != 0
     || !m1_thcm1_fixture_finite_vector(record->baseline_input, record->input_count)
     || !m1_thcm1_fixture_finite_vector(record->perturbed_input, record->input_count)
     || !m1_thcm1_fixture_finite_vector(record->baseline_output, record->output_count)
     || !m1_thcm1_fixture_finite_vector(record->perturbed_output, record->output_count)
     || !m1_thcm1_fixture_finite_vector(
           record->normalization, record->normalization_count)) {
    return 0;
  }
  if(!m1_thcm1_fixture_pair_metadata_valid(record)) {
    return 0;
  }
  for(size_t i = 0; i < record->normalization_count; ++i) {
    if(!(record->normalization[i] > 0.0)) {
      return 0;
    }
  }
  if(record->sensitivity_count > 0) {
    if(record->sensitivity_start >= record->input_count
       || record->sensitivity_count > record->input_count - record->sensitivity_start) {
      return 0;
    }
    int changed = 0;
    for(size_t i = record->sensitivity_start;
        i < record->sensitivity_start + record->sensitivity_count; ++i) {
      changed |= record->baseline_input[i] != record->perturbed_input[i];
    }
    if(!changed) {
      return 0;
    }
    for(size_t i = 0; i < record->input_count; ++i) {
      if((i < record->sensitivity_start
          || i >= record->sensitivity_start + record->sensitivity_count)
         && record->baseline_input[i] != record->perturbed_input[i]) {
        return 0;
      }
    }
  }
  else {
    /* A control must actually have identical consumed inputs. */
    for(size_t i = 0; i < record->input_count; ++i) {
      if(record->baseline_input[i] != record->perturbed_input[i]) {
        return 0;
      }
    }
  }
  return 1;
}

static int m1_thcm1_fixture_compare_ids(const void *left, const void *right) {
  return strcmp(*(const char *const *)left, *(const char *const *)right);
}

static int m1_thcm1_fixture_unique_ids(
      const m1_thcm1_fixture_collection *restrict collection,
      char *restrict error,
      const size_t error_size) {
  if(collection->record_count < 2) {
    return 1;
  }
  const char **ids = NULL;
  if(collection->record_count <= SIZE_MAX / sizeof(*ids)) {
    ids = (const char **)malloc(collection->record_count * sizeof(*ids));
  }
  if(ids == NULL) {
    m1_thcm1_fixture_set_error(error, error_size, "could not allocate fixture ID index");
    return 0;
  }
  /* Sort borrowed IDs, not records: replay order stays unchanged. */
  for(int pair_ids = 0; pair_ids < 2; ++pair_ids) {
    for(size_t i = 0; i < collection->record_count; ++i) {
      ids[i] = pair_ids ? collection->records[i].pair_id : collection->records[i].case_id;
    }
    qsort(ids, collection->record_count, sizeof(*ids), m1_thcm1_fixture_compare_ids);
    for(size_t i = 1; i < collection->record_count; ++i) {
      if(strcmp(ids[i - 1], ids[i]) == 0) {
        free(ids);
        m1_thcm1_fixture_set_error(
              error, error_size, pair_ids ? "fixture has duplicate pair IDs"
                                         : "fixture has duplicate case IDs");
        return 0;
      }
    }
  }
  free(ids);
  return 1;
}

/* TestData/radiation/PROVENANCE.md defines this packed, little-endian layout.
 * Decode integers explicitly so fixture bytes do not depend on host endian or
 * compiler struct padding. */
static int m1_thcm1_fixture_read_u64(FILE *restrict file, uint64_t *restrict value) {
  unsigned char bytes[8];
  if(fread(bytes, 1, sizeof(bytes), file) != sizeof(bytes)) {
    return 0;
  }
  *value = 0;
  for(size_t i = 0; i < sizeof(bytes); ++i) {
    *value |= (uint64_t)bytes[i] << (8 * i);
  }
  return 1;
}

static int m1_thcm1_fixture_read_size(FILE *restrict file, size_t *restrict value) {
  uint64_t stored = 0;
  if(!m1_thcm1_fixture_read_u64(file, &stored) || stored > SIZE_MAX) {
    return 0;
  }
  *value = (size_t)stored;
  return 1;
}

static int m1_thcm1_fixture_read_string(FILE *restrict file, char **restrict value) {
  size_t length = 0;
  if(!m1_thcm1_fixture_read_size(file, &length) || length == 0 || length == SIZE_MAX) {
    return 0;
  }
  char *string = (char *)malloc(length + 1);
  if(string == NULL) {
    return 0;
  }
  if(fread(string, 1, length, file) != length) {
    free(string);
    return 0;
  }
  for(size_t i = 0; i < length; ++i) {
    if(string[i] == '\0' || isspace((unsigned char)string[i])) {
      free(string);
      return 0;
    }
  }
  string[length] = '\0';
  *value = string;
  return 1;
}

static int m1_thcm1_fixture_read_vector(
      FILE *restrict file,
      double **restrict values,
      size_t *restrict count) {
  size_t length = 0;
  if(!m1_thcm1_fixture_read_size(file, &length) || length == 0
     || length > SIZE_MAX / sizeof(double)) {
    return 0;
  }
  double *vector = (double *)malloc(length * sizeof(*vector));
  if(vector == NULL) {
    return 0;
  }
  for(size_t i = 0; i < length; ++i) {
    uint64_t bits = 0;
    if(!m1_thcm1_fixture_read_u64(file, &bits)) {
      free(vector);
      return 0;
    }
    memcpy(&vector[i], &bits, sizeof(bits));
  }
  *values = vector;
  *count = length;
  return 1;
}

static int m1_thcm1_fixture_load_binary(
      FILE *restrict file,
      const char *restrict operation,
      const size_t expected_input_count,
      const size_t expected_output_count,
      m1_thcm1_fixture_collection *restrict collection,
      char *restrict error,
      const size_t error_size) {
  if(sizeof(double) != sizeof(uint64_t) || DBL_MANT_DIG != 53 || DBL_MAX_EXP != 1024
     || !m1_thcm1_fixture_read_string(file, &collection->operation)
     || strcmp(collection->operation, operation) != 0
     || !m1_thcm1_fixture_read_string(file, &collection->policy)
     || !m1_thcm1_fixture_read_size(file, &collection->record_count)
     || collection->record_count == 0
     || collection->record_count > SIZE_MAX / sizeof(*collection->records)) {
    m1_thcm1_fixture_set_error(error, error_size, "invalid binary fixture header");
    m1_thcm1_fixture_free(collection);
    return 0;
  }
  collection->records = (m1_thcm1_fixture_record *)calloc(
        collection->record_count, sizeof(*collection->records));
  if(collection->records == NULL) {
    m1_thcm1_fixture_set_error(
          error, error_size, "binary fixture record allocation failed");
    m1_thcm1_fixture_free(collection);
    return 0;
  }
  for(size_t i = 0; i < collection->record_count; ++i) {
    m1_thcm1_fixture_record *record = &collection->records[i];
    record->sensitivity_start = SIZE_MAX;
    size_t baseline_input_count = 0;
    size_t perturbed_input_count = 0;
    size_t baseline_output_count = 0;
    size_t perturbed_output_count = 0;
    uint64_t status[6] = { 0 };
    int ok = m1_thcm1_fixture_read_string(file, &record->case_id)
             && m1_thcm1_fixture_read_string(file, &record->pair_id)
             && m1_thcm1_fixture_read_string(file, &record->origin)
             && m1_thcm1_fixture_read_string(file, &record->seed_id)
             && m1_thcm1_fixture_read_string(file, &record->family)
             && m1_thcm1_fixture_read_string(file, &record->perturbation)
             && m1_thcm1_fixture_read_size(file, &record->sensitivity_start)
             && m1_thcm1_fixture_read_size(file, &record->sensitivity_count)
             && m1_thcm1_fixture_read_vector(
                   file, &record->baseline_input, &baseline_input_count)
             && m1_thcm1_fixture_read_vector(
                   file, &record->perturbed_input, &perturbed_input_count)
             && m1_thcm1_fixture_read_vector(
                   file, &record->baseline_output, &baseline_output_count)
             && m1_thcm1_fixture_read_vector(
                   file, &record->perturbed_output, &perturbed_output_count)
             && m1_thcm1_fixture_read_vector(
                   file, &record->normalization, &record->normalization_count);
    for(size_t j = 0; ok && j < sizeof(status) / sizeof(status[0]); ++j) {
      ok = m1_thcm1_fixture_read_u64(file, &status[j]);
    }
    ok = ok && baseline_input_count == perturbed_input_count
         && baseline_output_count == perturbed_output_count && status[0] == 1
         && status[1] == 1 && status[2] == 0 && status[3] == 0 && status[4] == 0
         && status[5] == 0;
    if(ok) {
      record->input_count = baseline_input_count;
      record->output_count = baseline_output_count;
      record->baseline_available = (int)status[0];
      record->perturbed_available = (int)status[1];
      ok = m1_thcm1_fixture_record_valid(
            record, expected_input_count, expected_output_count);
    }
    if(!ok) {
      m1_thcm1_fixture_set_errorf(
            error, error_size, "invalid binary fixture record %s at index %zu",
            record->case_id, i);
      m1_thcm1_fixture_free(collection);
      return 0;
    }
  }
  const int valid_end = fgetc(file) == EOF && !ferror(file);
  if(!valid_end) {
    m1_thcm1_fixture_set_error(
          error, error_size, "binary fixture has trailing data or a read error");
    m1_thcm1_fixture_free(collection);
    return 0;
  }
  if(!m1_thcm1_fixture_unique_ids(collection, error, error_size)) {
    m1_thcm1_fixture_free(collection);
    return 0;
  }
  return 1;
}

int m1_thcm1_fixture_load(
      const char *restrict path,
      const char *restrict operation,
      const size_t expected_input_count,
      const size_t expected_output_count,
      m1_thcm1_fixture_collection *restrict collection,
      char *restrict error,
      const size_t error_size) {
  if(collection == NULL || path == NULL || operation == NULL) {
    m1_thcm1_fixture_set_error(error, error_size, "null fixture loader argument");
    return 0;
  }
  *collection = (m1_thcm1_fixture_collection){ 0 };
  FILE *file = fopen(path, "rb");
  if(file == NULL) {
    m1_thcm1_fixture_set_error(error, error_size, "fixture file could not be opened");
    return 0;
  }
  unsigned char binary_magic[8];
  if(fread(binary_magic, 1, sizeof(binary_magic), file) == sizeof(binary_magic)
     && memcmp(binary_magic, "M1THCMB1", sizeof(binary_magic)) == 0) {
    const int loaded = m1_thcm1_fixture_load_binary(
          file, operation, expected_input_count, expected_output_count, collection,
          error, error_size);
    fclose(file);
    return loaded;
  }
  if(fseek(file, 0, SEEK_SET) != 0) {
    fclose(file);
    m1_thcm1_fixture_set_error(error, error_size, "fixture file could not be rewound");
    return 0;
  }
  int version = 0;
  char *magic = m1_thcm1_fixture_read_token(file);
  const int valid_magic = magic != NULL && strcmp(magic, M1_THCM1_FIXTURE_MAGIC) == 0;
  free(magic);
  if(!valid_magic || fscanf(file, "%d", &version) != 1
     || version != M1_THCM1_FIXTURE_VERSION
     || !m1_thcm1_fixture_expect(file, "operation")
     || !m1_thcm1_fixture_string(file, &collection->operation)
     || strcmp(collection->operation, operation) != 0
     || !m1_thcm1_fixture_expect(file, "policy")
     || !m1_thcm1_fixture_string(file, &collection->policy)
     || !m1_thcm1_fixture_expect(file, "record_count")
     || fscanf(file, "%zu", &collection->record_count) != 1
     || collection->record_count == 0) {
    fclose(file);
    m1_thcm1_fixture_set_error(error, error_size, "invalid fixture header");
    m1_thcm1_fixture_free(collection);
    return 0;
  }
  collection->records = (m1_thcm1_fixture_record *)calloc(
        collection->record_count, sizeof(*collection->records));
  if(collection->records == NULL) {
    fclose(file);
    m1_thcm1_fixture_set_error(error, error_size, "fixture record allocation failed");
    m1_thcm1_fixture_free(collection);
    return 0;
  }
  for(size_t i = 0; i < collection->record_count; ++i) {
    m1_thcm1_fixture_record *record = &collection->records[i];
    record->sensitivity_start = SIZE_MAX;
    size_t baseline_input_count = 0;
    size_t perturbed_input_count = 0;
    size_t baseline_output_count = 0;
    size_t perturbed_output_count = 0;
    int ok
          = m1_thcm1_fixture_expect(file, "record_begin")
            && m1_thcm1_fixture_expect(file, "case_id")
            && m1_thcm1_fixture_string(file, &record->case_id)
            && m1_thcm1_fixture_expect(file, "pair_id")
            && m1_thcm1_fixture_string(file, &record->pair_id)
            && m1_thcm1_fixture_expect(file, "origin")
            && m1_thcm1_fixture_string(file, &record->origin)
            && m1_thcm1_fixture_expect(file, "seed_id")
            && m1_thcm1_fixture_string(file, &record->seed_id)
            && m1_thcm1_fixture_expect(file, "family")
            && m1_thcm1_fixture_string(file, &record->family)
            && m1_thcm1_fixture_expect(file, "perturbation")
            && m1_thcm1_fixture_string(file, &record->perturbation)
            && m1_thcm1_fixture_expect(file, "sensitivity_start")
            && fscanf(file, "%zu", &record->sensitivity_start) == 1
            && m1_thcm1_fixture_expect(file, "sensitivity_count")
            && fscanf(file, "%zu", &record->sensitivity_count) == 1
            && m1_thcm1_fixture_vector(
                  file, "baseline_input", &record->baseline_input, &baseline_input_count)
            && m1_thcm1_fixture_vector(
                  file, "perturbed_input", &record->perturbed_input,
                  &perturbed_input_count)
            && m1_thcm1_fixture_vector(
                  file, "baseline_output", &record->baseline_output,
                  &baseline_output_count)
            && m1_thcm1_fixture_vector(
                  file, "perturbed_output", &record->perturbed_output,
                  &perturbed_output_count)
            && m1_thcm1_fixture_vector(
                  file, "normalization", &record->normalization,
                  &record->normalization_count)
            && m1_thcm1_fixture_expect(file, "baseline_available")
            && fscanf(file, "%d", &record->baseline_available) == 1
            && m1_thcm1_fixture_expect(file, "perturbed_available")
            && fscanf(file, "%d", &record->perturbed_available) == 1
            && m1_thcm1_fixture_expect(file, "baseline_status")
            && fscanf(file, "%d", &record->baseline_status) == 1
            && m1_thcm1_fixture_expect(file, "perturbed_status")
            && fscanf(file, "%d", &record->perturbed_status) == 1
            && m1_thcm1_fixture_expect(file, "baseline_reference_status")
            && fscanf(file, "%d", &record->baseline_reference_status) == 1
            && m1_thcm1_fixture_expect(file, "perturbed_reference_status")
            && fscanf(file, "%d", &record->perturbed_reference_status) == 1
            && m1_thcm1_fixture_expect(file, "record_end") &&
            /* Do not assign record counts, or validate vectors through them,
             * until both sides of each pair have independently parsed counts. */
            baseline_input_count == perturbed_input_count
            && baseline_output_count == perturbed_output_count;
    if(ok) {
      record->input_count = baseline_input_count;
      record->output_count = baseline_output_count;
      ok = m1_thcm1_fixture_record_valid(
            record, expected_input_count, expected_output_count);
    }
    if(!ok) {
      fclose(file);
      m1_thcm1_fixture_set_errorf(
            error, error_size, "invalid fixture record %s at index %zu", record->case_id,
            i);
      m1_thcm1_fixture_free(collection);
      return 0;
    }
  }
  char *trailing = NULL;
  const int valid_end = m1_thcm1_fixture_expect(file, "end")
                        && (trailing = m1_thcm1_fixture_read_token(file)) == NULL
                        && !ferror(file);
  free(trailing);
  fclose(file);
  if(!valid_end) {
    m1_thcm1_fixture_set_error(
          error, error_size, "fixture has a missing end marker, trailing data or a read error");
    m1_thcm1_fixture_free(collection);
    return 0;
  }
  if(!m1_thcm1_fixture_unique_ids(collection, error, error_size)) {
    m1_thcm1_fixture_free(collection);
    return 0;
  }
  return 1;
}

/* The pointwise and strict-pair fixture policies use an inclusive bound on
 * error after input-defined normalization. The pointwise policy also adds an
 * absolute term to a symmetric relative term with a floor. The shared
 * ghl_pert_test_fail_with_tolerance() uses a strict absolute gate and a
 * trusted-value-relative response, so it cannot express these policies.
 * The baseline/response policy below uses that helper where its semantics
 * match, with a documented range fallback for unrepresentable scaling. */
static int m1_thcm1_fixture_compare_one(
      const double actual,
      const double expected,
      const double normalization,
      const double absolute_tolerance,
      const double relative_tolerance,
      const double floor_value,
      double *restrict scaled_error) {
  if(!isfinite(actual) || !isfinite(expected) || !isfinite(normalization)
     || !(normalization > 0.0)) {
    return 0;
  }
  const double error = fabs(actual - expected) / normalization;
  const double bound = absolute_tolerance
                       + relative_tolerance
                               * fmax(
                                     floor_value, fmax(fabs(actual) / normalization,
                                                       fabs(expected) / normalization));
  if(scaled_error != NULL) {
    *scaled_error = error;
  }
  return isfinite(error) && isfinite(bound) && error <= bound;
}

/* Return the relative response represented by the retained THC pair.  This
 * follows the trusted/perturbed convention used by ghl_pert_test_fail: the
 * perturbation response is a target-local numerical envelope, while a zero
 * trusted value uses the recorded normalization as its scale. */
static int m1_thcm1_fixture_response_relative(
      const double trusted,
      const double perturbed,
      const double normalization,
      double *restrict response_relative) {
  if(response_relative == NULL || !isfinite(trusted) || !isfinite(perturbed)
     || !isfinite(normalization) || !(normalization > 0.0)) {
    return 0;
  }
  const double response = fabs(perturbed - trusted);
  if(!isfinite(response)) {
    return 0;
  }
  if(trusted != 0.0) {
    *response_relative = response / fabs(trusted);
  }
  else {
    *response_relative = response / normalization;
  }
  return isfinite(*response_relative) && *response_relative >= 0.0;
}

/* This fallback is only for values that cannot be represented after exact
 * power-of-two normalization: normalization would overflow a value or erase
 * a nonzero subnormal.  Keep the original arithmetic for those range cases so
 * that the shared helper never sees collapsed inputs. */
static int m1_thcm1_fixture_compare_with_response_range_fallback(
      const double actual,
      const double trusted,
      const double perturbed,
      const double normalization,
      const double absolute_tolerance,
      const double relative_tolerance,
      double *restrict scaled_error,
      double *restrict response_relative) {
  if(!isfinite(actual) || !isfinite(trusted) || !isfinite(perturbed)
     || !isfinite(normalization) || !(normalization > 0.0)
     || !isfinite(absolute_tolerance) || absolute_tolerance < 0.0
     || !isfinite(relative_tolerance) || relative_tolerance < 0.0) {
    return 0;
  }
  const double absolute_error = fabs(actual - trusted) / normalization;
  if(!isfinite(absolute_error)) {
    return 0;
  }
  if(scaled_error != NULL) {
    *scaled_error = absolute_error;
  }
  double reference_response = 0.0;
  if(!m1_thcm1_fixture_response_relative(
           trusted, perturbed, normalization, &reference_response)) {
    return 0;
  }
  if(response_relative != NULL) {
    *response_relative = reference_response;
  }
  if(absolute_error <= absolute_tolerance) {
    return 1;
  }
  const double current_relative
        = trusted != 0.0 ? fabs(1.0 - actual / trusted) : absolute_error;
  const double bound = fmax(4.0 * reference_response, relative_tolerance);
  return isfinite(current_relative) && current_relative >= 0.0 && isfinite(bound)
         && current_relative <= bound;
}

static int m1_thcm1_fixture_scale_value(
      const double value,
      const int normalization_exponent,
      double *restrict scaled_value) {
  const double scaled = scalbn(value, -normalization_exponent);
  if(!isfinite(scaled) || (value != 0.0 && (scaled == 0.0 || fabs(scaled) < DBL_MIN))) {
    return 0;
  }
  *scaled_value = scaled;
  return 1;
}

/* Apply the common trusted-baseline response bar through the established
 * unit-test comparator.  Power-of-two scaling preserves the nonzero relative
 * comparisons while placing values near the normalization range.  The
 * shared helper has a strict `<` absolute gate, so its absolute tolerance is
 * disabled here and the retained inclusive `<=` gate is applied separately.
 * For a zero trusted value, scale the relative tolerance by the normalized
 * input scale so the shared helper's absolute zero-reference semantics match
 * the original normalized response rule. */
static int m1_thcm1_fixture_compare_with_response(
      const double actual,
      const double trusted,
      const double perturbed,
      const double normalization,
      const double absolute_tolerance,
      const double relative_tolerance,
      double *restrict scaled_error,
      double *restrict response_relative) {
  if(!isfinite(actual) || !isfinite(trusted) || !isfinite(perturbed)
     || !isfinite(normalization) || !(normalization > 0.0)
     || !isfinite(absolute_tolerance) || absolute_tolerance < 0.0
     || !isfinite(relative_tolerance) || relative_tolerance < 0.0) {
    return 0;
  }

  const double absolute_error = fabs(actual - trusted) / normalization;
  if(!isfinite(absolute_error)) {
    return 0;
  }
  if(scaled_error != NULL) {
    *scaled_error = absolute_error;
  }
  double reference_response = 0.0;
  if(!m1_thcm1_fixture_response_relative(
           trusted, perturbed, normalization, &reference_response)) {
    return 0;
  }
  if(response_relative != NULL) {
    *response_relative = reference_response;
  }
  if(absolute_error <= absolute_tolerance) {
    return 1;
  }

  const double relative_bound = fmax(4.0 * reference_response, relative_tolerance);
  const double current_relative
        = trusted != 0.0 ? fabs(1.0 - actual / trusted) : absolute_error;
  if(!isfinite(relative_bound) || !isfinite(current_relative)
     || current_relative < 0.0) {
    return 0;
  }

  double scaled_actual = 0.0;
  double scaled_trusted = 0.0;
  if(trusted == 0.0) {
    scaled_actual = actual / normalization;
    if(!isfinite(scaled_actual) || (actual != 0.0 && scaled_actual == 0.0)) {
      return m1_thcm1_fixture_compare_with_response_range_fallback(
            actual, trusted, perturbed, normalization, absolute_tolerance,
            relative_tolerance, scaled_error, response_relative);
    }
  }
  else {
    int normalization_exponent = 0;
    (void)frexp(normalization, &normalization_exponent);
    if(!m1_thcm1_fixture_scale_value(actual, normalization_exponent, &scaled_actual)
       || !m1_thcm1_fixture_scale_value(
             trusted, normalization_exponent, &scaled_trusted)) {
      return m1_thcm1_fixture_compare_with_response_range_fallback(
            actual, trusted, perturbed, normalization, absolute_tolerance,
            relative_tolerance, scaled_error, response_relative);
    }
  }

  /* The prior policy derives one exact response bound before comparing the
   * current value. Pass that bound as rel_tol and make the shared helper's
   * perturbation response zero, avoiding a second relative-response formula
   * with different rounding. For a zero trusted value, actual is normalized
   * directly above and the same bound is already in normalized units. */
  const double scaled_computed = scaled_actual;
  const double scaled_reference = trusted == 0.0 ? 0.0 : scaled_trusted;
  const int relative_failure = ghl_pert_test_fail_with_tolerance(
        scaled_reference, scaled_computed, scaled_reference, relative_bound, 0.0);
  return !relative_failure;
}

int m1_thcm1_fixture_compare_baseline_response(
      const m1_thcm1_fixture_record *restrict record,
      const char *restrict policy,
      const double *restrict computed_normalization,
      const double *restrict computed_baseline,
      m1_thcm1_fixture_comparison_report *restrict report,
      char *restrict error,
      const size_t error_size) {
  const int pointwise_policy
        = policy != NULL && strcmp(policy, "pointwise_a1_a2_v1") == 0;
  const int a1_a2_policy
        = policy != NULL && strcmp(policy, "source_a1_a2_rate_normalized_v1") == 0;
  const int strict_policy
        = policy != NULL
          && (strcmp(policy, "strict_relative_2e-12_propagated_response_v1") == 0
              || strcmp(policy, "strict_relative_2e-12_propagated_response_current_v2")
                       == 0);
  if(record == NULL || (!pointwise_policy && !a1_a2_policy && !strict_policy)
     || !m1_thcm1_fixture_record_valid(record, 0, 0) || computed_normalization == NULL
     || computed_baseline == NULL) {
    m1_thcm1_fixture_set_error(
          error, error_size, "invalid baseline-response comparison argument");
    return 0;
  }
  if(report != NULL) {
    *report = (m1_thcm1_fixture_comparison_report){ .max_scaled_error = 0.0,
                                                    .component = 0,
                                                    .gate = "baseline" };
  }

  /* Identical consumed inputs are controls.  They may certify fixed-input
   * agreement, but a nonzero retained output response is malformed metadata,
   * not sensitivity evidence. */
  if(record->sensitivity_count == 0) {
    for(size_t component = 0; component < record->output_count; ++component) {
      if(record->baseline_output[component] != record->perturbed_output[component]) {
        if(report != NULL) {
          report->gate = "response";
        }
        m1_thcm1_fixture_set_error(
              error, error_size,
              "input-invariant control has a nonzero reference response");
        return 0;
      }
    }
  }

  for(size_t component = 0; component < record->output_count; ++component) {
    const double normalization = computed_normalization[component];
    const double trusted = record->baseline_output[component];
    const double perturbed = record->perturbed_output[component];
    const double actual = computed_baseline[component];
    if(!isfinite(normalization) || !(normalization > 0.0)
       || normalization != record->normalization[component] || !isfinite(trusted)
       || !isfinite(perturbed) || !isfinite(actual)) {
      if(error != NULL && error_size > 0) {
        snprintf(
              error, error_size,
              "fixture %s component %zu has nonfinite or mismatched values",
              record->case_id, component);
      }
      return 0;
    }
    double scaled_error = 0.0;
    double response_relative = 0.0;
    const double absolute_tolerance = pointwise_policy || a1_a2_policy ? 2.0e-12 : 0.0;
    const double relative_tolerance
          = pointwise_policy || a1_a2_policy ? 2.0e-10 : 2.0e-12;
    const int baseline_ok = m1_thcm1_fixture_compare_with_response(
          actual, trusted, perturbed, normalization, absolute_tolerance,
          relative_tolerance, &scaled_error, &response_relative);
    if(!baseline_ok) {
      if(report != NULL) {
        report->max_scaled_error = fmax(report->max_scaled_error, scaled_error);
        report->component = component;
        report->gate = "baseline";
      }
      if(error != NULL && error_size > 0) {
        snprintf(
              error, error_size,
              "fixture policy %s case %s pair %s component %zu baseline comparison "
              "failed: actual=%.17g trusted=%.17g response_relative=%.17g",
              policy, record->case_id, record->pair_id, component, actual, trusted,
              response_relative);
      }
      return 0;
    }
    if(report != NULL) {
      report->max_scaled_error = fmax(report->max_scaled_error, scaled_error);
    }
  }
  return 1;
}

int m1_thcm1_fixture_compare_paired(
      const m1_thcm1_fixture_record *restrict record,
      const double *restrict computed_normalization,
      const double *restrict computed_baseline,
      const double *restrict computed_perturbed,
      m1_thcm1_fixture_comparison_report *restrict report,
      char *restrict error,
      const size_t error_size) {
  if(record == NULL || !m1_thcm1_fixture_pair_metadata_valid(record)
     || record->baseline_output == NULL || record->perturbed_output == NULL
     || record->normalization == NULL || computed_normalization == NULL
     || computed_baseline == NULL || computed_perturbed == NULL
     || record->output_count == 0
     || record->normalization_count != record->output_count) {
    m1_thcm1_fixture_set_error(error, error_size, "invalid paired comparison argument");
    return 0;
  }
  if(report != NULL) {
    *report = (m1_thcm1_fixture_comparison_report){ .max_scaled_error = 0.0,
                                                    .component = 0,
                                                    .gate = "baseline" };
  }
  if(record->sensitivity_count == 0) {
    for(size_t component = 0; component < record->output_count; ++component) {
      if(record->baseline_output[component] != record->perturbed_output[component]) {
        if(report != NULL) {
          report->gate = "response";
        }
        m1_thcm1_fixture_set_error(
              error, error_size,
              "input-invariant control has a nonzero reference response");
        return 0;
      }
    }
  }
  const double absolute_tolerance = 2.0e-12;
  const double relative_tolerance = 2.0e-10;
  const double floor_value = 1.0e-13;
  const double *expected[2] = { record->baseline_output, record->perturbed_output };
  const double *computed[2] = { computed_baseline, computed_perturbed };
  const char *gate_names[3] = { "baseline", "perturbed", "response" };
  for(size_t component = 0; component < record->output_count; ++component) {
    if(!isfinite(computed_normalization[component])
       || !(computed_normalization[component] > 0.0)
       || computed_normalization[component] != record->normalization[component]) {
      if(error != NULL && error_size > 0) {
        snprintf(
              error, error_size,
              "fixture %s component %zu has invalid input normalization",
              record->case_id, component);
      }
      return 0;
    }
    for(int gate = 0; gate < 3; ++gate) {
      double actual = 0.0, reference = 0.0;
      if(gate < 2) {
        actual = computed[gate][component];
        reference = expected[gate][component];
      }
      else {
        actual = computed[1][component] - computed[0][component];
        reference = expected[1][component] - expected[0][component];
      }
      double scaled_error = 0.0;
      if(!m1_thcm1_fixture_compare_one(
               actual, reference, computed_normalization[component], absolute_tolerance,
               relative_tolerance, floor_value, &scaled_error)) {
        if(report != NULL) {
          report->max_scaled_error = fmax(report->max_scaled_error, scaled_error);
          report->component = component;
          report->gate = gate_names[gate];
        }
        if(error != NULL && error_size > 0) {
          snprintf(
                error, error_size,
                "fixture %s component %zu %s comparison failed: actual=%.17g "
                "reference=%.17g normalization=%.17g",
                record->case_id, component, gate_names[gate], actual, reference,
                computed_normalization[component]);
        }
        return 0;
      }
      if(report != NULL) {
        report->max_scaled_error = fmax(report->max_scaled_error, scaled_error);
      }
    }
  }
  return 1;
}

int m1_thcm1_fixture_compare_paired_strict_relative(
      const m1_thcm1_fixture_record *restrict record,
      const double *restrict computed_normalization,
      const double *restrict computed_baseline,
      const double *restrict computed_perturbed,
      m1_thcm1_fixture_comparison_report *restrict report,
      char *restrict error,
      const size_t error_size) {
  if(record == NULL || computed_normalization == NULL || computed_baseline == NULL
     || computed_perturbed == NULL || record->normalization == NULL
     || record->baseline_output == NULL || record->perturbed_output == NULL
     || record->output_count == 0) {
    m1_thcm1_fixture_set_error(
          error, error_size, "invalid strict paired comparison argument");
    return 0;
  }
  if(report != NULL) {
    *report = (m1_thcm1_fixture_comparison_report){ .max_scaled_error = 0.0,
                                                    .component = 0,
                                                    .gate = "baseline" };
  }
  const double *expected[2] = { record->baseline_output, record->perturbed_output };
  const double *computed[2] = { computed_baseline, computed_perturbed };
  const char *gate_names[3] = { "baseline", "perturbed", "response" };
  for(size_t component = 0; component < record->output_count; ++component) {
    if(!isfinite(computed_normalization[component])
       || !(computed_normalization[component] > 0.0)
       || computed_normalization[component] != record->normalization[component]) {
      if(error != NULL && error_size > 0) {
        snprintf(
              error, error_size,
              "fixture %s component %zu has invalid input normalization",
              record->case_id, component);
      }
      return 0;
    }
    for(int gate = 0; gate < 3; ++gate) {
      const double actual = gate < 2 ? computed[gate][component]
                                     : computed[1][component] - computed[0][component];
      const double reference = gate < 2
                                     ? expected[gate][component]
                                     : expected[1][component] - expected[0][component];
      double denominator = fmax(1.0e-300, fmax(fabs(actual), fabs(reference)));
      double bound = 2.0e-12;
      if(gate == 2) {
        /* A difference of nearby outputs must not be tested relative to its
         * nearly zero value. Propagate the two retained endpoint bounds.
         * The three double subtractions contribute at most 4*DBL_EPSILON
         * times the summed endpoint scales (using epsilon conservatively
         * rather than half-epsilon). This is not independent sensitivity
         * evidence beyond the endpoint comparisons. */
        const double baseline_scale = fmax(
              1.0e-300,
              fmax(fabs(computed[0][component]), fabs(expected[0][component])));
        const double perturbed_scale = fmax(
              1.0e-300,
              fmax(fabs(computed[1][component]), fabs(expected[1][component])));
        denominator = baseline_scale + perturbed_scale;
        bound += 4.0 * DBL_EPSILON;
      }
      double scaled_error = 0.0;
      if(!m1_thcm1_fixture_compare_one(
               actual, reference, denominator, bound, 0.0, 0.0, &scaled_error)) {
        if(report != NULL) {
          report->max_scaled_error = fmax(report->max_scaled_error, scaled_error);
          report->component = component;
          report->gate = gate_names[gate];
        }
        if(error != NULL && error_size > 0) {
          snprintf(
                error, error_size, "fixture %s component %zu %s comparison failed",
                record->case_id, component, gate_names[gate]);
        }
        return 0;
      }
      if(report != NULL) {
        report->max_scaled_error = fmax(report->max_scaled_error, scaled_error);
      }
    }
  }
  return 1;
}

int m1_thcm1_fixture_compare_envelope(
      const m1_thcm1_fixture_record *restrict record,
      const char *restrict policy,
      const double *restrict computed_normalization,
      const double *restrict computed_baseline,
      m1_thcm1_fixture_comparison_report *restrict report,
      char *restrict error,
      const size_t error_size) {
  const char control_prefix[] = "input_invariant_control:";
  if(record == NULL || policy == NULL || strcmp(policy, "pointwise_a1_a2_v1") != 0) {
    m1_thcm1_fixture_set_error(
          error, error_size, "unsupported pointwise envelope policy");
    return 0;
  }
  if(!m1_thcm1_fixture_pair_id_matches_case_suffix(record, "-pair")) {
    m1_thcm1_fixture_set_error(
          error, error_size, "invalid pointwise envelope pair metadata");
    return 0;
  }
  const int is_control = record->sensitivity_count == 0;
  if((is_control
      && strncmp(record->perturbation, control_prefix, sizeof(control_prefix) - 1) != 0)
     || (!is_control
         && strncmp(record->perturbation, control_prefix, sizeof(control_prefix) - 1)
                  == 0)) {
    m1_thcm1_fixture_set_error(
          error, error_size, "pointwise perturbation role disagrees with input pair");
    return 0;
  }
  return m1_thcm1_fixture_compare_baseline_response(
        record, policy, computed_normalization, computed_baseline, report, error,
        error_size);
}
