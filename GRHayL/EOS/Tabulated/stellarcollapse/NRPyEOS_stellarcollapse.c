#include "NRPyEOS_stellarcollapse.h"
#include "ghl_nrpyeos_tabulated.h"
#include "ghl_io.h"

#ifndef GHL_DISABLE_HDF5
#include <limits.h>
#include <stdint.h>

ghl_error_codes_t NRPyEOS_stellarcollapse_check_dimensions(
      const int n_rho,
      const int n_temperature,
      const int n_ye,
      size_t *npoints) {
  if(n_rho < 2 || n_temperature < 2 || n_ye < 2) {
    return ghl_error_invalid_eos_table;
  }

  size_t limit = (size_t)INT_MAX / NRPyEOS_ntablekeys;
  const size_t allocation_limit = SIZE_MAX / sizeof(double) / NRPyEOS_ntablekeys;
  if(allocation_limit < limit) {
    limit = allocation_limit;
  }

  const int dimensions[3] = { n_rho, n_temperature, n_ye };
  size_t product = 1;
  for(int i = 0; i < 3; i++) {
    if(product > limit / (size_t)dimensions[i]) {
      return ghl_error_invalid_eos_table;
    }
    product *= (size_t)dimensions[i];
  }
  *npoints = product;
  return ghl_success;
}

#include <hdf5.h>
#include <math.h>
#include <stdlib.h>

#include "../NRPyEOS_hdf5_helpers.h"

static const char *dataset_names[NRPyEOS_sc_n_quantities] = {
  "Abar", "Xa",      "Xh",      "Xn",      "Xp",    "Zbar",      "cs2",
  "dedt", "dpderho", "dpdrhoe", "entropy", "gamma", "logenergy", "logpress",
  "mu_e", "mu_n",    "mu_p",    "muhat",   "munu",
};

// This function checks for the existence of the "have_rel_cs2" flag in
// the input HDF5 file. If it's not present, the sound speed squared is
// assumed to be non relativistic. Otherwise, the dataset value is used.
static ghl_error_codes_t table_has_rel_cs2(hid_t file, bool *has_rel_cs2) {
  *has_rel_cs2 = false;
  if(H5Lexists(file, "have_rel_cs2", H5P_DEFAULT) > 0) {
    int *flag = NULL;
    const ghl_error_codes_t error
          = NRPyEOS_hdf5_read_int_dataset(file, "have_rel_cs2", 1, (void **)&flag);
    if(error != ghl_success) {
      return error;
    }
    *has_rel_cs2 = (*flag != 0);
    free(flag);
  }
  return ghl_success;
}

static ghl_error_codes_t check_axis(
      const double *restrict axis,
      const int n,
      const double scale,
      const double offset,
      const char *name) {
  const bool transform_changes_values = scale != 1.0 || offset != 0.0;
  const double x0 = axis[0];
  const double dx = axis[1] - x0;
  const double dxi = 1.0 / dx;
  const double y0 = x0 * scale + offset;
  const double dy = axis[1] * scale + offset - y0;
  const double dyi = 1.0 / dy;
  if(!isfinite(x0) || !isfinite(dx) || dx <= 0.0 || !isfinite(dxi) || dxi <= 0.0
     || !isfinite(y0) || !isfinite(dy) || dy <= 0.0 || !isfinite(dyi) || dyi <= 0.0) {
    ghl_warn("Invalid %s axis origin or spacing in EOS table\n", name);
    return ghl_error_invalid_eos_table;
  }

  double previous_x = x0;
  double previous_y = y0;
  for(int i = 0; i < n; i++) {
    const double x = axis[i];
    const double y = x * scale + offset;
    if(!isfinite(x) || (i > 0 && x <= previous_x) || !isfinite((x - x0) * dxi)
       || fabs((x - x0) * dxi - (double)i) > 1.0e-8 || !isfinite(y)
       || (i > 0 && y <= previous_y)
       || (transform_changes_values
           && (!isfinite((y - y0) * dyi)
               || fabs((y - y0) * dyi - (double)i) > 1.0e-8))) {
      ghl_warn("Invalid %s axis at node %d in EOS table\n", name, i);
      return ghl_error_invalid_eos_table;
    }
    previous_x = x;
    previous_y = y;
  }
  return ghl_success;
}

#endif

#define GHL_GOTO_CLEANUP_IF_ERROR(call) \
  err = call;                           \
  if(err != ghl_success) {              \
    goto cleanup;                       \
  }                                     \

ghl_error_codes_t NRPyEOS_stellarcollapse_read_table(
      const char *filepath,
      NRPyEOS_stellarcollapse_t **sc) {
  *sc = NULL;
#ifdef GHL_DISABLE_HDF5
  return ghl_error_used_disabled_hdf5;
#else
  hid_t file_id = H5Fopen(filepath, H5F_ACC_RDONLY, H5P_DEFAULT);
  if(file_id < 0) {
    return ghl_error_could_not_open_file;
  }

  // We use calloc here because of the cleanup on failure. free() works fine
  // on NULL pointers, but table->data[n] would be uninitialized with malloc(),
  // which would cause free() to segfault if we encounter an error.
  NRPyEOS_stellarcollapse_t *table = (NRPyEOS_stellarcollapse_t *)calloc(1, sizeof(NRPyEOS_stellarcollapse_t));
  if(table == NULL) {
    H5Fclose(file_id);
    return ghl_error_out_of_memory;
  }

  ghl_error_codes_t err = table_has_rel_cs2(file_id, &table->cs2_is_relativistic);
  if(err != ghl_success) {
    free(table);
    H5Fclose(file_id);
    return err;
  }
  ghl_info("Table '%s' contains %srelativistic sound speed\n",
           filepath,
           table->cs2_is_relativistic ? "" : "non-");

  // Scalar quantities
  int *nr_ptr = NULL;
  int *ny_ptr = NULL;
  int *nt_ptr = NULL;
  double *es_ptr = NULL;

  GHL_GOTO_CLEANUP_IF_ERROR(NRPyEOS_hdf5_read_int_dataset(file_id, "pointsrho", 1, (void **)&nr_ptr));
  GHL_GOTO_CLEANUP_IF_ERROR(NRPyEOS_hdf5_read_int_dataset(file_id, "pointstemp", 1, (void **)&nt_ptr));
  GHL_GOTO_CLEANUP_IF_ERROR(NRPyEOS_hdf5_read_int_dataset(file_id, "pointsye", 1, (void **)&ny_ptr));
  GHL_GOTO_CLEANUP_IF_ERROR(NRPyEOS_hdf5_read_double_dataset(file_id, "energy_shift", 1, (void **)&es_ptr));

  table->n_rho = *nr_ptr;
  table->n_temperature = *nt_ptr;
  table->n_ye = *ny_ptr;
  table->energy_shift = *es_ptr;

  size_t size;
  GHL_GOTO_CLEANUP_IF_ERROR(NRPyEOS_stellarcollapse_check_dimensions(
        table->n_rho, table->n_temperature, table->n_ye, &size));

  // Basic tabulated quantities
  GHL_GOTO_CLEANUP_IF_ERROR(NRPyEOS_hdf5_read_double_dataset(file_id, "logrho", table->n_rho, (void **)&table->log10_rho));
  GHL_GOTO_CLEANUP_IF_ERROR(NRPyEOS_hdf5_read_double_dataset(file_id, "logtemp", table->n_temperature, (void **)&table->log10_temperature));
  GHL_GOTO_CLEANUP_IF_ERROR(NRPyEOS_hdf5_read_double_dataset(file_id, "ye", table->n_ye, (void **)&table->ye));

  const double log10_to_ln = log(10.0);
  GHL_GOTO_CLEANUP_IF_ERROR(check_axis(
        table->log10_rho, table->n_rho, log10_to_ln, log(CGS_TO_CODE_DENSITY),
        "density"));
  GHL_GOTO_CLEANUP_IF_ERROR(check_axis(
        table->log10_temperature, table->n_temperature, log10_to_ln, 0.0,
        "temperature"));
  GHL_GOTO_CLEANUP_IF_ERROR(
        check_axis(table->ye, table->n_ye, 1.0, 0.0, "electron-fraction"));

  // Tabulated data
  for(int n = 0; n < NRPyEOS_sc_n_quantities; n++) {
    GHL_GOTO_CLEANUP_IF_ERROR(NRPyEOS_hdf5_read_double_dataset(file_id, dataset_names[n], size, (void **)&table->data[n]));
  }

cleanup:

  free(nr_ptr);
  free(ny_ptr);
  free(nt_ptr);
  free(es_ptr);
  if(err != ghl_success) {
    free(table->log10_rho);
    free(table->log10_temperature);
    free(table->ye);
    for(int n = 0; n < NRPyEOS_sc_n_quantities; n++) {
      free(table->data[n]);
    }
    free(table);
    table = NULL;
  }

  H5Fclose(file_id);
  *sc = table;
  return err;
#endif
}

void NRPyEOS_stellarcollapse_free_table(NRPyEOS_stellarcollapse_t *table) {
#ifdef GHL_DISABLE_HDF5
  GHL_HDF5_ERROR_IF_USED;
#else
  if(!table) {
    return;
  }
  for(int n = 0; n < NRPyEOS_sc_n_quantities; n++) {
    free(table->data[n]);
  }
  free(table->log10_rho);
  free(table->log10_temperature);
  free(table->ye);
  free(table);
#endif
}
