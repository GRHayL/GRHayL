#include "NRPyEOS_hdf5_helpers.h"

ghl_error_codes_t NRPyEOS_hdf5_read_dataset(
      hid_t file_id,
      ghl_hdf5_t dtype,
      const char *dataset_name,
      const size_t expected_size,
      void **data) {
#ifndef GHL_USE_HDF5
  GHL_HDF5_ERROR_IF_USED;
#else
  ghl_error_codes_t error = ghl_success;

  const hid_t dataset_types[2] = { H5T_NATIVE_INT, H5T_NATIVE_DOUBLE };
  const size_t dtype_size[2] = { sizeof(int), sizeof(double) };

  hid_t dataset_id = -1;
  hid_t dataspace_id = -1;
  void *array = NULL;

  const hid_t hdf5_dtype = dataset_types[dtype];

  dataset_id = H5Dopen(file_id, dataset_name, H5P_DEFAULT);
  if(dataset_id < 0) {
    error = ghl_error_hdf5_dataset_could_not_open;
    goto cleanup_and_return;
  }

  dataspace_id = H5Dget_space(dataset_id);
  if(dataspace_id < 0) {
    error = ghl_error_hdf5_dataset_invalid_ndims;
    goto cleanup_and_return;
  }

  int ndims = H5Sget_simple_extent_ndims(dataspace_id);
  if(ndims < 0) {
    error = ghl_error_hdf5_dataset_invalid_ndims;
    goto cleanup_and_return;
  }

  size_t total_size = 1;
  if(ndims > 0) {
    hsize_t dims[ndims];

    if(H5Sget_simple_extent_dims(dataspace_id, dims, NULL) < 0) {
      error = ghl_error_hdf5_dataset_invalid_ndims;
      goto cleanup_and_return;
    }

    for(int i = 0; i < ndims; i++) {
      total_size *= dims[i];
    }
  }

  if(total_size != expected_size) {
    error = ghl_error_hdf5_dataset_size_mismatch;
    goto cleanup_and_return;
  }

  array = malloc(total_size * dtype_size[dtype]);
  if(array == NULL) {
    error = ghl_error_out_of_memory;
    goto cleanup_and_return;
  }

  herr_t status = H5Dread(dataset_id, hdf5_dtype,
                         H5S_ALL, H5S_ALL, H5P_DEFAULT, array);
  if(status < 0) {
    error = ghl_error_hdf5_dataset_could_not_open; /* ideally: could_not_read */
    goto cleanup_and_return;
  }

  *data = array;
  array = NULL;

cleanup_and_return:
  free(array);

  if(dataspace_id >= 0) {
    H5Sclose(dataspace_id);
  }

  if(dataset_id >= 0) {
    H5Dclose(dataset_id);
  }

 return error;
#endif
}

ghl_error_codes_t NRPyEOS_hdf5_read_int_dataset(
      hid_t file_id,
      const char *dataset_name,
      const size_t expected_size,
      void **data) {
  return NRPyEOS_hdf5_read_dataset(
    file_id, GHL_HDF5_INT, dataset_name, expected_size, data);
}

ghl_error_codes_t NRPyEOS_hdf5_read_double_dataset(
      hid_t file_id,
      const char *dataset_name,
      const size_t expected_size,
      void **data) {
  return NRPyEOS_hdf5_read_dataset(
    file_id, GHL_HDF5_DOUBLE, dataset_name, expected_size, data);
}
