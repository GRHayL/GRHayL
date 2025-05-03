#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"

int main() {
  hid_t file_id, dataset_id, dataspace_id;
  herr_t status;
  file_id = H5Fopen("test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    fprintf(stderr, "Failed to open HDF5 file!\\n");
    exit(1);
  }
  dataset_id = H5Dopen2(file_id, "test_dataset", H5P_DEFAULT);
  if (dataset_id < 0) {
    fprintf(stderr, "Failed to open HDF5 dataset!\\n");
    exit(1);
  }
  dataspace_id = H5Dget_space(dataset_id);
  if (H5Sget_simple_extent_ndims(dataspace_id) != 2) {
    fprintf(stderr, "Unexpected dataspace rank!\\n");
    exit(1);
  }
  hsize_t dims[2];
  status = H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
  if (status < 0 || dims[0] != 5 || dims[1] != 10) {
    fprintf(stderr, "Unexpected dataspace dimensions!\\n");
    exit(1);
  }
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);
  printf("HDF5 installation is OK!\\n");
  return 0;
}
