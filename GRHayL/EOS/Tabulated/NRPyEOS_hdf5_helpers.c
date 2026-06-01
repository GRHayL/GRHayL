#include "ghl_nrpyeos_tabulated.h"
#include "NRPyEOS_hdf5_helpers.h"

void *
NRPyEOS_hdf5_read_dataset(hid_t file_id, ghl_hdf5_t dtype, const char *dataset_name) {
#ifndef GHL_USE_HDF5
  GHL_HDF5_ERROR_IF_USED;
#else
  const hid_t dataset_types[2] = { H5T_NATIVE_INT, H5T_NATIVE_DOUBLE };
  const size_t dtype_size[2] = { sizeof(int), sizeof(double) };

  const hid_t hdf5_dtype = dataset_types[dtype];

  // Open the dataset
  hid_t dataset_id = H5Dopen(file_id, dataset_name, H5P_DEFAULT);
  if(dataset_id < 0) {
    ghl_error("Dataset '%s' not found.\n", dataset_name);
  }

  // Get the dataspace and dimensions
  hid_t dataspace_id = H5Dget_space(dataset_id);

  int ndims = H5Sget_simple_extent_ndims(dataspace_id);
  if(ndims < 1) {
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    ghl_error("Invalid dimension for dataset '%s' (%d)\n", dataset_name, ndims);
  }

  hsize_t dims[ndims];
  H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

  // Allocate memory for the 3D array
  size_t total_size = 1;
  for(int i = 0; i < ndims; i++) {
    total_size *= dims[i];
  }

  // Read the 3D array
  void *array = malloc(total_size * dtype_size[dtype]);
  herr_t status = H5Dread(dataset_id, hdf5_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, array);
  if(status < 0) {
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    ghl_error("Problem reading dataset '%s'.\n", dataset_name);
  }

  // Close the dataspace and dataset
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);

  return array;
#endif
}

int *NRPyEOS_hdf5_read_int_dataset(hid_t file_id, const char *dataset_name) {
  return (int *)NRPyEOS_hdf5_read_dataset(file_id, GHL_HDF5_INT, dataset_name);
}

double *NRPyEOS_hdf5_read_double_dataset(hid_t file_id, const char *dataset_name) {
  return (double *)NRPyEOS_hdf5_read_dataset(file_id, GHL_HDF5_DOUBLE, dataset_name);
}
