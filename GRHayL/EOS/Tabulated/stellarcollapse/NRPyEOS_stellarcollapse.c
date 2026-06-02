#include <hdf5.h>
#include <math.h>
#include <stdlib.h>

#include "../NRPyEOS_hdf5_helpers.h"
#include "NRPyEOS_stellarcollapse.h"
#include "ghl_io.h"

static const char *dataset_names[NRPyEOS_sc_n_quantities] = {
  "Abar", "Xa",      "Xh",      "Xn",      "Xp",    "Zbar",      "cs2",
  "dedt", "dpderho", "dpdrhoe", "entropy", "gamma", "logenergy", "logpress",
  "mu_e", "mu_n",    "mu_p",    "muhat",   "munu",
};

// This function checks for the existence of the "have_rel_cs2" flag in
// the input HDF5 file. If it's not present, the sound speed squared is
// assumed to be non relativistic. Otherwise, the dataset value is used.
static bool table_has_rel_cs2(hid_t file) {
  int32_t flag = 0;

  if(H5Lexists(file, "have_rel_cs2", H5P_DEFAULT) > 0) {
    hid_t dset = H5Dopen2(file, "have_rel_cs2", H5P_DEFAULT);
    H5Dread(dset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, &flag);
    H5Dclose(dset);
  }

  return (flag != 0);
}

#define GHL_GOTO_CLEANUP_IF_ERROR(call) \
  err = call;                           \
  if(err != ghl_success) {              \
    goto cleanup;                       \
  }                                     \

ghl_error_codes_t NRPyEOS_stellarcollapse_read_table(
      const char *filepath,
      NRPyEOS_stellarcollapse_t **sc) {
#ifndef GHL_USE_HDF5
  GHL_HDF5_ERROR_IF_USED;
#else
  hid_t file_id = H5Fopen(filepath, H5F_ACC_RDONLY, H5P_DEFAULT);
  if(file_id < 0) {
    return ghl_error_could_not_open_file;
  }

  NRPyEOS_stellarcollapse_t *table = (NRPyEOS_stellarcollapse_t *)malloc(sizeof(NRPyEOS_stellarcollapse_t));
  if(table == NULL) {
    return ghl_error_out_of_memory;
  }

  table->cs2_is_relativistic = table_has_rel_cs2(file_id);
  ghl_info("Table '%s' contains %srelativistic sound speed\n",
           filepath,
           table->cs2_is_relativistic ? "" : "non-");

  ghl_error_codes_t err = ghl_success;

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

  // Basic tabulated quantities
  GHL_GOTO_CLEANUP_IF_ERROR(NRPyEOS_hdf5_read_double_dataset(file_id, "logrho", table->n_rho, (void **)&table->log10_rho));
  GHL_GOTO_CLEANUP_IF_ERROR(NRPyEOS_hdf5_read_double_dataset(file_id, "logtemp", table->n_temperature, (void **)&table->log10_temperature));
  GHL_GOTO_CLEANUP_IF_ERROR(NRPyEOS_hdf5_read_double_dataset(file_id, "ye", table->n_ye, (void **)&table->ye));

  // Tabulated data
  const size_t size = table->n_rho * table->n_temperature * table->n_ye;
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
#ifndef GHL_USE_HDF5
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
