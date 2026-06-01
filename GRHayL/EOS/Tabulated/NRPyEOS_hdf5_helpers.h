/**
 * @file NRPyEOS_hdf5_helpers.h
 * @author Leo Werneck
 *
 * @brief Helper functions for reading and writing datasets in HDF5 files.
 */

#ifndef NRPYEOS_HDF5_HELPERS_H
#define NRPYEOS_HDF5_HELPERS_H

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
 *
 * @return A pointer to the allocated memory containing the dataset data.
 *         The caller is responsible for freeing this memory. Returns NULL on error.
 */
void *NRPyEOS_hdf5_read_dataset(
      hid_t file_id,
      ghl_hdf5_t dtype,
      const char *dataset_name,
      const size_t expected_size);

// Convenience functions
int *NRPyEOS_hdf5_read_int_dataset(
      hid_t file_id,
      const char *dataset_name,
      const size_t expected_size);

double *NRPyEOS_hdf5_read_double_dataset(
      hid_t file_id,
      const char *dataset_name,
      const size_t expected_size);

#endif // NRPYEOS_HDF5_HELPERS_H
