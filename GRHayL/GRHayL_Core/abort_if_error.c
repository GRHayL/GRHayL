#include "ghl.h"

void ghl_abort_if_error(const ghl_error_codes_t error) {
  switch(error) {
    case ghl_success:
      break;
    case ghl_error_u0_singular:
      ghl_abort("u^0 evaluated to NaN while speed-limiting the velocity.");
      break;
    case ghl_error_unknown_eos_type:
      ghl_abort("Unknown EOS found in struct element 'eos_type'.");
      break;
    case ghl_error_invalid_c2p_key:
      ghl_abort("Con2Prim select_method function received an invalid con2prim method key.");
      break;
    case ghl_error_neg_rho:
      ghl_abort("Negative rho returned from conservative-to-primitive process.");
      break;
    case ghl_error_neg_pressure:
      ghl_abort("Negative pressure returned from conservative-to-primitive process.");
      break;
    case ghl_error_neg_vsq:
      ghl_abort("Negative v^2 returned from conservative-to-primitive process.");
      break;
    case ghl_error_c2p_max_iter:
      ghl_abort("Maximum iteration reached during the conservative-to-primitive process without finding a solution.");
      break;
    case ghl_error_c2p_singular:
      ghl_abort("Infinite/singular quantities found after computing the primitives during the conservative-to-primitive process.");
      break;
    case ghl_error_root_not_bracketed:
      ghl_abort("Interval in the root-finding method does not bracket a root.");
      break;
    case ghl_error_table_bisection:
      ghl_abort("Table bisection failed while locating the interpolation interval.");
      break;
    case ghl_error_table_max_ye:
      ghl_abort("Input Y_e is too large and outside table bounds. Table interpolation cannot be performed.");
      break;
    case ghl_error_table_min_ye:
      ghl_abort("Input Y_e is too small and outside table bounds. Table interpolation cannot be performed.");
      break;
    case ghl_error_table_max_T:
      ghl_abort("Input temperature is too large and outside table bounds. Table interpolation cannot be performed.");
      break;
    case ghl_error_table_min_T:
      ghl_abort("Input temperature is too small and outside table bounds. Table interpolation cannot be performed.");
      break;
    case ghl_error_table_max_rho:
      ghl_abort("Input rho is too large and outside table bounds. Table interpolation cannot be performed.");
      break;
    case ghl_error_table_min_rho:
      ghl_abort("Input rho is too small and outside table bounds. Table interpolation cannot be performed.");
      break;
    case ghl_error_exceed_table_vars:
      ghl_abort("The number of requested quantities for interpolation is greater than the number of table quantities.");
      break;
    case ghl_error_table_neg_energy:
      ghl_abort("While interpolating temperature, found eps+energy_shift < 0.0. Interpolation cannot be performed.\n");
      break;
    case ghl_error_invalid_utsq:
      ghl_abort("While computing u^0 squared in Con2Prim, found an invalid value (<0 or >>1)");
      break;
    case ghl_error_invalid_Z:
      ghl_abort("While computing Z in Con2Prim, found an invalid value (<0 or >>1)");
      break;
    case ghl_error_newman_invalid_discriminant:
      ghl_abort("Found negative discriminant inside Newman1D Con2Prim");
      break;
    case ghl_error_used_disabled_hdf5:
      ghl_abort("Tried to use HDF5 function, but HDF5 is disabled");
      break;
    case ghl_error_out_of_memory:
      ghl_abort("Failed to allocate memory");
      break;
    case ghl_error_eos_struct_is_null:
      ghl_abort("Provided EOS struct pointer is NULL");
      break;
    case ghl_error_invalid_eos_type:
      ghl_abort("Invalid EOS type");
      break;
    case ghl_error_invalid_eos_table_type:
      ghl_abort("Invalid EOS table type");
      break;
    case ghl_error_could_not_open_file:
      ghl_abort("Could not open file");
      break;
    case ghl_error_hdf5_dataset_could_not_open:
      ghl_abort("Could not open HDF5 dataset");
      break;
    case ghl_error_hdf5_dataset_could_not_read:
      ghl_abort("Could not read HDF5 dataset");
      break;
    case ghl_error_hdf5_dataset_invalid_ndims:
      ghl_abort("HDF5 dataset had invalid number of dimensions");
      break;
    case ghl_error_hdf5_dataset_size_mismatch:
      ghl_abort("HDF5 dataset size mismatch");
      break;
    case ghl_error_invalid_rho_atm:
      ghl_abort("rho_atm must be specified.");
      break;
    case ghl_error_rho_min_gt_rho_max:
      ghl_abort("rho_min cannot be greater than rho_max.");
      break;
    case ghl_error_invalid_press_atm:
      ghl_abort("press_atm must be specified.");
      break;
    case ghl_error_press_min_gt_press_max:
      ghl_abort("press_min cannot be greater than press_max.");
      break;
    case ghl_error_invalid_Y_e_atm:
      ghl_abort("Y_e_atm must be specified.");
      break;
    case ghl_error_Y_e_min_gt_Y_e_max:
      ghl_abort("Y_e_min cannot be greater than Y_e_max.");
      break;
    case ghl_error_invalid_T_atm:
      ghl_abort("T_atm must be specified.");
      break;
    case ghl_error_T_min_gt_T_max:
      ghl_abort("T_min cannot be greater than T_max.");
      break;
    case ghl_error_invalid_fermi_dirac_integral_key:
      ghl_abort("Unsupported Fermi-Dirac integral key.");
      break;
  }
}
