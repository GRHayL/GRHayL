#include "ghl.h"

GHL_HOST_DEVICE
void ghl_read_error_codes(
      const ghl_error_codes_t error) {
  switch(error) {
    case ghl_success:
      ghl_info("Success! No errors here.");
      break;
    case ghl_error_u0_singular:
      ghl_abort("u^0 evaluated to NaN while speed-limiting the velocity.");
      break;
    case ghl_error_unknown_eos_type:
      ghl_error("Unknown EOS found in struct element 'eos_type'.");
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
      ghl_abort("???");
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
  }
}
