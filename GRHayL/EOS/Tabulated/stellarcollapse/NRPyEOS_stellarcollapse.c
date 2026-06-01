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

static NRPyEOS_stellarcollapse_t *NRPyEOS_new_stellarcollapse_table() {
  NRPyEOS_stellarcollapse_t *t = NULL;
  t = (NRPyEOS_stellarcollapse_t *)malloc(sizeof(NRPyEOS_stellarcollapse_t));
  if(t == NULL) {
    ghl_error("Could not allocate memory for stellar collapse table\n");
  }
  return t;
}

NRPyEOS_stellarcollapse_t *NRPyEOS_stellarcollapse_read_table(const char *filepath) {
#ifndef GHL_USE_HDF5
  GHL_HDF5_ERROR_IF_USED;
#else
  hid_t file_id = H5Fopen(filepath, H5F_ACC_RDONLY, H5P_DEFAULT);
  if(file_id < 0) {
    ghl_error("Could not open file '%s'\n", filepath);
  }

  NRPyEOS_stellarcollapse_t *table = NRPyEOS_new_stellarcollapse_table();

  // Scalar quantities
  table->cs2_is_relativistic = table_has_rel_cs2(file_id);
  table->n_rho = *NRPyEOS_hdf5_read_int_dataset(file_id, "pointsrho");
  table->n_temperature = *NRPyEOS_hdf5_read_int_dataset(file_id, "pointstemp");
  table->n_ye = *NRPyEOS_hdf5_read_int_dataset(file_id, "pointsye");
  table->energy_shift = *NRPyEOS_hdf5_read_double_dataset(file_id, "energy_shift");

  ghl_info("Table '%s' contains %srelativistic sound speed\n",
           filepath,
           table->cs2_is_relativistic ? "" : "non-");

  // Basic tabulated quantities
  table->ye = NRPyEOS_hdf5_read_double_dataset(file_id, "ye");
  table->log10_temperature = NRPyEOS_hdf5_read_double_dataset(file_id, "logtemp");
  table->log10_rho = NRPyEOS_hdf5_read_double_dataset(file_id, "logrho");

  // Tabulated data
  for(int n = 0; n < NRPyEOS_sc_n_quantities; n++) {
    table->data[n] = NRPyEOS_hdf5_read_double_dataset(file_id, dataset_names[n]);
  }

  H5Fclose(file_id);

  return table;
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
    if(table->data[n]) {
      free(table->data[n]);
    }
  }
  free(table);
#endif
}
