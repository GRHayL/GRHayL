/**
 * @file NRPyEOS_hdf5_helpers.h
 * @author Leo Werneck
 *
 * @brief Helper functions for reading and writing datasets in HDF5 files.
 */

#ifndef NRPYEOS_HDF5_HELPERS_H
#define NRPYEOS_HDF5_HELPERS_H

#include "ghl.h"
#include <hdf5.h>

/*
 * @brief Convenience enumeration to index HDF5 data types.
 */
typedef enum {
  GHL_HDF5_INT,
  GHL_HDF5_DOUBLE,
} ghl_hdf5_t;

/**
 * @brief Reads an HDF5 dataset from a file.
 *
 * @param file_id The HDF5 file identifier.
 * @param dtype The data type of the dataset (I32 or F64).
 * @param dataset_name The name of the dataset to read.
 * @param expected_size Expected dataset size.
 * @param data Pointer to the pointer that will hold the data.
 *
 * @return Return code, indicating success or failure.
 * @retval ghl_success Success.
 * @retval ghl_error_hdf5_dataset_could_not_open Could not open HDF5 dataset.
 * @retval ghl_error_hdf5_dataset_invalid_ndims HDF5 dataset has wrong number of
 *         dimensions.
 * @retval ghl_error_hdf5_dataset_size_mismatch HDF5 dataset size does not match
 *         expected size.
 * @retval ghl_error_out_of_memory Could not allocate memory to store the dataset
 *         data.
 */
ghl_error_codes_t NRPyEOS_hdf5_read_dataset(
      hid_t file_id,
      ghl_hdf5_t dtype,
      const char *dataset_name,
      const size_t expected_size,
      void **data);

ghl_error_codes_t NRPyEOS_hdf5_read_int_dataset(
      hid_t file_id,
      const char *dataset_name,
      const size_t expected_size,
      void **data);

ghl_error_codes_t NRPyEOS_hdf5_read_double_dataset(
      hid_t file_id,
      const char *dataset_name,
      const size_t expected_size,
      void **data);

#endif // NRPYEOS_HDF5_HELPERS_H
