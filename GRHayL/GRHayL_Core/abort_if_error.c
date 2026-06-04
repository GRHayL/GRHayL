#include "ghl.h"

#define GHL_CASE_ERROR(code, msg) case code: ghl_Error(code, msg); break;

void ghl_abort_if_error(const ghl_error_codes_t error) {
  switch(error) {
    case ghl_success: break;
    GHL_CASE_ERROR(ghl_error_u0_singular,
                   "u^0 evaluated to NaN while speed-limiting the velocity.\n")
    GHL_CASE_ERROR(ghl_error_unknown_eos_type,
                   "Unknown EOS found in struct element 'eos_type'.\n");
    GHL_CASE_ERROR(ghl_error_invalid_c2p_key,
                   "Con2Prim select_method function received an invalid con2prim method key.\n");
    GHL_CASE_ERROR(ghl_error_neg_rho,
                   "Negative rho returned from conservative-to-primitive process.\n");
    GHL_CASE_ERROR(ghl_error_neg_pressure,
                   "Negative pressure returned from conservative-to-primitive process.\n");
    GHL_CASE_ERROR(ghl_error_neg_vsq,
                   "Negative v^2 returned from conservative-to-primitive process.\n");
    GHL_CASE_ERROR(ghl_error_c2p_max_iter,
                   "Maximum iteration reached during the conservative-to-primitive process without finding a solution.\n");
    GHL_CASE_ERROR(ghl_error_c2p_singular,
                   "Infinite/singular quantities found after computing the primitives during the conservative-to-primitive process.\n");
    GHL_CASE_ERROR(ghl_error_root_not_bracketed,
                   "Interval in the root-finding method does not bracket a root.\n");
    GHL_CASE_ERROR(ghl_error_table_bisection,
                   "Table bisection failed while locating the interpolation interval.\n");
    GHL_CASE_ERROR(ghl_error_table_max_ye,
                   "Input Y_e is too large and outside table bounds. Table interpolation cannot be performed.\n");
    GHL_CASE_ERROR(ghl_error_table_min_ye,
                   "Input Y_e is too small and outside table bounds. Table interpolation cannot be performed.\n");
    GHL_CASE_ERROR(ghl_error_table_max_T,
                   "Input temperature is too large and outside table bounds. Table interpolation cannot be performed.\n");
    GHL_CASE_ERROR(ghl_error_table_min_T,
                   "Input temperature is too small and outside table bounds. Table interpolation cannot be performed.\n");
    GHL_CASE_ERROR(ghl_error_table_max_rho,
                   "Input rho is too large and outside table bounds. Table interpolation cannot be performed.\n");
    GHL_CASE_ERROR(ghl_error_table_min_rho,
                   "Input rho is too small and outside table bounds. Table interpolation cannot be performed.\n");
    GHL_CASE_ERROR(ghl_error_exceed_table_vars,
                   "The number of requested quantities for interpolation is greater than the number of table quantities.\n");
    GHL_CASE_ERROR(ghl_error_table_neg_energy,
                   "While interpolating temperature, found eps+energy_shift < 0.0. Interpolation cannot be performed.\n");
    GHL_CASE_ERROR(ghl_error_invalid_utsq,
                   "While computing u^0 squared in Con2Prim, found an invalid value (<0 or >>1)\n");
    GHL_CASE_ERROR(ghl_error_invalid_Z,
                   "While computing Z in Con2Prim, found an invalid value (<0 or >>1)\n");
    GHL_CASE_ERROR(ghl_error_newman_invalid_discriminant,
                   "Found negative discriminant inside Newman1D Con2Prim\n");
    GHL_CASE_ERROR(ghl_error_used_disabled_hdf5,
                   "Tried to use HDF5 function, but HDF5 is disabled\n");
    GHL_CASE_ERROR(ghl_error_out_of_memory,
                   "Failed to allocate memory\n");
    GHL_CASE_ERROR(ghl_error_eos_struct_is_null,
                   "Provided EOS struct pointer is NULL\n");
    GHL_CASE_ERROR(ghl_error_invalid_eos_type,
                   "Invalid EOS type\n");
    GHL_CASE_ERROR(ghl_error_invalid_eos_table_type,
                   "Invalid EOS table type\n");
    GHL_CASE_ERROR(ghl_error_could_not_open_file,
                   "Could not open file\n");
    GHL_CASE_ERROR(ghl_error_hdf5_dataset_could_not_open,
                   "Could not open HDF5 dataset\n");
    GHL_CASE_ERROR(ghl_error_hdf5_dataset_could_not_read,
                   "Could not read HDF5 dataset\n");
    GHL_CASE_ERROR(ghl_error_hdf5_dataset_invalid_ndims,
                   "HDF5 dataset had invalid number of dimensions\n");
    GHL_CASE_ERROR(ghl_error_hdf5_dataset_size_mismatch,
                   "HDF5 dataset size mismatch\n");
    GHL_CASE_ERROR(ghl_error_invalid_rho_atm,
                   "rho_atm must be specified.\n");
    GHL_CASE_ERROR(ghl_error_rho_min_gt_rho_max,
                   "rho_min cannot be greater than rho_max.\n");
    GHL_CASE_ERROR(ghl_error_invalid_press_atm,
                   "press_atm must be specified.\n");
    GHL_CASE_ERROR(ghl_error_press_min_gt_press_max,
                   "press_min cannot be greater than press_max.\n");
    GHL_CASE_ERROR(ghl_error_invalid_Y_e_atm,
                   "Y_e_atm must be specified.\n");
    GHL_CASE_ERROR(ghl_error_Y_e_min_gt_Y_e_max,
                   "Y_e_min cannot be greater than Y_e_max.\n");
    GHL_CASE_ERROR(ghl_error_invalid_T_atm,
                   "T_atm must be specified.\n");
    GHL_CASE_ERROR(ghl_error_T_min_gt_T_max,
                   "T_min cannot be greater than T_max.\n");
    GHL_CASE_ERROR(ghl_error_invalid_fermi_dirac_integral_key,
                   "Unsupported Fermi-Dirac integral key.\n");
  }
}
