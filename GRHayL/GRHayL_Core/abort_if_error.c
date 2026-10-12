#include "ghl.h"

/* Error codes index immutable diagnostics; success has no fatal message. */
static const char *const error_messages[] = {
  [ghl_success] = NULL,
  [ghl_error_u0_singular] = "Velocity limiting could not produce a finite, subluminal, "
            "within-cap result, or u^0 was non-finite.\n",
  [ghl_error_unknown_eos_type] = "Unknown EOS found in struct element 'eos_type'.\n",
  [ghl_error_invalid_c2p_key] = "Con2Prim select_method function received an "
                                       "invalid con2prim method key.\n",
  [ghl_error_neg_rho] = "Negative rho returned from conservative-to-primitive process.\n",
  [ghl_error_neg_pressure] = "Negative pressure returned from conservative-to-primitive process.\n",
  [ghl_error_neg_vsq] = "Negative v^2 returned from conservative-to-primitive process.\n",
  [ghl_error_c2p_max_iter] = "Maximum iteration reached during the conservative-to-primitive process "
            "without finding a solution.\n",
  [ghl_error_c2p_singular] = "Infinite/singular quantities found after computing the primitives during "
            "the conservative-to-primitive process.\n",
  [ghl_error_root_not_bracketed] = "Interval in the root-finding method does not bracket a root.\n",
  [ghl_error_table_bisection] = "Table bisection failed while locating the interpolation interval.\n",
  [ghl_error_table_max_ye] = "Input Y_e is too large and outside table bounds. "
                                    "Table interpolation cannot be performed.\n",
  [ghl_error_table_min_ye] = "Input Y_e is too small and outside table bounds. "
                                    "Table interpolation cannot be performed.\n",
  [ghl_error_table_max_T] = "Input temperature is too large and outside table "
                                   "bounds. Table interpolation cannot be performed.\n",
  [ghl_error_table_min_T] = "Input temperature is too small and outside table "
                                   "bounds. Table interpolation cannot be performed.\n",
  [ghl_error_table_max_rho] = "Input rho is too large and outside table bounds. "
                                     "Table interpolation cannot be performed.\n",
  [ghl_error_table_min_rho] = "Input rho is too small and outside table bounds. "
                                     "Table interpolation cannot be performed.\n",
  [ghl_error_exceed_table_vars] = "The number of requested quantities for interpolation is greater than the "
            "number of table quantities.\n",
  [ghl_error_table_neg_energy] = "While interpolating temperature, found eps+energy_shift < 0.0. "
            "Interpolation cannot be performed.\n",
  [ghl_error_invalid_utsq] = "While computing u^0 squared in Con2Prim, found an "
                                    "invalid value (<0 or >>1)\n",
  [ghl_error_invalid_Z] = "While computing Z in Con2Prim, found an invalid value (<0 or >>1)\n",
  [ghl_error_newman_invalid_discriminant] = "Found negative discriminant inside Newman1D Con2Prim\n",
  [ghl_error_used_disabled_hdf5] = "Tried to use HDF5 function, but HDF5 is disabled\n",
  [ghl_error_out_of_memory] = "Failed to allocate memory\n",
  [ghl_error_eos_struct_is_null] = "Provided EOS struct pointer is NULL\n",
  [ghl_error_invalid_eos_type] = "Invalid EOS type\n",
  [ghl_error_invalid_eos_table_type] = "Invalid EOS table type\n",
  [ghl_error_could_not_open_file] = "Could not open file\n",
  [ghl_error_hdf5_dataset_could_not_open] = "Could not open HDF5 dataset\n",
  [ghl_error_hdf5_dataset_could_not_read] = "Could not read HDF5 dataset\n",
  [ghl_error_hdf5_dataset_invalid_ndims] = "HDF5 dataset had invalid number of dimensions\n",
  [ghl_error_hdf5_dataset_size_mismatch] = "HDF5 dataset size mismatch\n",
  [ghl_error_invalid_rho_atm] = "rho_atm must be specified.\n",
  [ghl_error_rho_min_gt_rho_max] = "rho_min cannot be greater than rho_max.\n",
  [ghl_error_invalid_press_atm] = "press_atm must be specified.\n",
  [ghl_error_press_min_gt_press_max] = "press_min cannot be greater than press_max.\n",
  [ghl_error_invalid_Y_e_atm] = "Y_e_atm must be specified.\n",
  [ghl_error_Y_e_min_gt_Y_e_max] = "Y_e_min cannot be greater than Y_e_max.\n",
  [ghl_error_invalid_T_atm] = "T_atm must be specified.\n",
  [ghl_error_T_min_gt_T_max] = "T_min cannot be greater than T_max.\n",
  [ghl_error_invalid_fermi_dirac_integral_key] = "Unsupported Fermi-Dirac integral key.\n",
  [ghl_error_nn_c2p_model_is_null] = "Neural-network Con2Prim model pointer is NULL.\n",
  [ghl_error_nn_c2p_invalid_dimensions] = "Neural-network Con2Prim model has invalid layer dimensions.\n",
  [ghl_error_nn_c2p_invalid_input_index] = "Neural-network Con2Prim model has an invalid input index.\n",
  [ghl_error_nn_c2p_missing_array] = "Neural-network Con2Prim model is missing one or more arrays.\n",
  [ghl_error_nn_c2p_invalid_kind] = "Neural-network Con2Prim model has an "
                                           "unsupported transform or output kind.\n",
  [ghl_error_nn_c2p_invalid_number] = "Neural-network Con2Prim model contains "
                                             "non-finite or invalid numeric data.\n",
  [ghl_error_nrpyleakage_blocking] = "NRPyLeakage nucleon-blocking evaluator received invalid data or failed "
            "numerically.\n",
  [ghl_error_nrpyleakage_nonfinite_output] = "NRPyLeakage replaced one or more non-finite outputs.\n",
  [ghl_error_invalid_neos] = "Hybrid EOS piece count is outside the supported range.\n",
  [ghl_error_invalid_eos_parameters] = "EOS parameters are invalid or produce non-finite state.\n",
  [ghl_error_invalid_eos_table] = "EOS table dimensions, grid, or derived data are invalid.\n",
  [ghl_error_invalid_hlle_wavespeeds] = "Invalid HLLE wave speeds.\n",
  [ghl_error_flux_source_invalid_input] = "Flux/source operation received invalid input.\n",
  [ghl_error_m1_null_pointer] = "M1 routine received a NULL pointer.\n",
  [ghl_error_m1_invalid_epsilon_c] = "M1 epsilon_c must be finite and strictly between zero and one.\n",
  [ghl_error_m1_invalid_E_floor] = "M1 E_floor must be finite and strictly positive.\n",
  [ghl_error_m1_invalid_zeta_min] = "M1 zeta_min must be finite and strictly positive.\n",
  [ghl_error_m1_invalid_fd_epsilon_rel] = "M1 relative finite-difference epsilon must be finite and strictly "
            "positive.\n",
  [ghl_error_m1_invalid_fd_epsilon_abs] = "M1 absolute finite-difference epsilon must be finite and strictly "
            "positive.\n",
  [ghl_error_m1_invalid_newton_max_iterations] = "M1 Newton maximum iteration count must be strictly positive.\n",
  [ghl_error_m1_invalid_newton_tolerance] = "M1 Newton relative tolerance must be finite and strictly positive.\n",
  [ghl_error_m1_invalid_metric] = "M1 operation received an invalid metric.\n",
  [ghl_error_m1_invalid_implicit_jacobian] = "M1 implicit solve produced an invalid Jacobian.\n",
  [ghl_error_m1_invalid_state] = "M1 operation received an invalid radiation state.\n",
  [ghl_error_m1_implicit_admissibility] = "M1 implicit trial state is outside the admissible domain.\n",
  [ghl_error_m1_implicit_solve_failure] = "M1 implicit Newton solve failed to converge.\n",
  [ghl_error_m1_implicit_terminal_fallback] = "M1 implicit solve exhausted its terminal fallback policy.\n",
  [ghl_error_m1_microphysics_failure] = "M1 neutrino microphysics rate evaluation failed.\n",
  [ghl_error_m1_invalid_newton_absolute_tolerance] = "M1 Newton absolute tolerance must be finite and strictly positive.\n",
  [ghl_error_m1_invalid_closure_tolerance] = "M1 closure tolerance must be finite and strictly positive.\n",
  [ghl_error_m1_invalid_closure_max_iterations] = "M1 closure maximum iteration count must be strictly positive.\n",
  [ghl_error_m1_closure_residual_too_large] = "M1 closure residual exceeded the configured acceptance tolerance.\n",
  [ghl_error_m1_source_double_application] = "M1 source update was applied more than once to the same state.\n",
};

void ghl_abort_if_error(const ghl_error_codes_t error) {
  const size_t index = (size_t)error;
  if(index >= sizeof(error_messages) / sizeof(error_messages[0])) {
    return;
  }
  const char *message = error_messages[index];
  if(message != NULL) {
    ghl_Error(error, "%s", message);
  }
}
