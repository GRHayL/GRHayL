#include "ghl_unit_tests.h"

int main(int argc, char **argv) {

  // This is the SLy EOS with the values taken from the NRPyEOS python code

  const int neos = 4;
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300; //IGM default
  const double Gamma_th = 2.0; //Taken from magnetizedTOV.par
  const double rho_ppoly[3] = {2.44034e+07,3.78358e+11,2.62780e+12};
  const double Gamma_ppoly[4] = {1.58425,1.28733,0.62223,1.35692};
  const double k_ppoly0 = 6.80110e-9;

  ghl_eos_parameters eos = { 0 };
  ghl_error_codes_t error = ghl_initialize_hybrid_eos_functions_and_params(
        rho_b_min, rho_b_min, rho_b_max, neos, rho_ppoly, Gamma_ppoly, k_ppoly0,
        Gamma_th, &eos);
  if(error != ghl_success) {
    ghl_error("Four-piece EOS initialization failed with error %d.\n", error);
  }

  // Expected output values taken from NRPyEOS python code
  double k_comp[4] = {6.8010999999999996e-09, 1.0618444833278535e-06, 5.3275084583961501e+01, 3.9992069172431910e-08};
  double eps_comp[4] = {0.0, -2.490902788423614614e-04, 1.355533308962642014e-02, 7.652227660602063664e-03};
  double k_pert[4], eps_pert[4];
  for(int i=1; i<4; i++) {
    k_pert[i] = k_comp[i]*(1.0 + randf(-1,1)*1e-14);
    eps_pert[i] = eps_comp[i]*(1.0 + randf(-1,1)*1e-14);
  }

  for(int i=1; i<4; i++) {
    if(ghl_pert_test_fail(k_comp[i], eos.K_ppoly[i], k_pert[i]))
      ghl_error("unit_test_piecewise_polytrope has failed for K_ppoly.\n"
                   "For index %d, expected %e, computed %e, perturbed %e\n"
                   "%e\n", i, k_comp[i], eos.K_ppoly[i], k_pert[i],
                           ghl_pert_test_fail(k_comp[i], eos.K_ppoly[i], k_pert[i]));

    if(ghl_pert_test_fail(eps_comp[i], eos.eps_integ_const[i], eps_pert[i]))
      ghl_error("unit_test_piecewise_polytrope has failed for eps_integ_const.\n"
                   "For index %d, expected %e, computed %e, perturbed %e\n"
                   "relative error: %e\n", i, eps_comp[i], eos.eps_integ_const[i], eps_pert[i],
                           ghl_pert_test_fail(eps_comp[i], eos.eps_integ_const[i], eps_pert[i]));
  }

  for(int i = 0; i < neos - 1; i++) {
    const double expected_left = k_comp[i] * pow(rho_ppoly[i], Gamma_ppoly[i]);
    const double expected_right = k_comp[i + 1] * pow(rho_ppoly[i], Gamma_ppoly[i + 1]);
    if(!isfinite(eos.p_ppoly[i]) || relative_error(expected_left, eos.p_ppoly[i]) > 1e-14
       || relative_error(expected_left, expected_right) > 1e-14) {
      ghl_error(
            "Pressure transition %d is invalid: expected %.15e, got %.15e, adjacent "
            "%.15e.\n",
            i, expected_left, eos.p_ppoly[i], expected_right);
    }
  }
  if(eos.p_ppoly[neos - 1] != 0.0) {
    ghl_error("Unused four-piece pressure slot is not zero.\n");
  }

  const double one_rho_ppoly[1] = { 0.0 };
  const double one_Gamma_ppoly[1] = { 2.0 };
  ghl_eos_parameters one_piece_eos = { 0 };
  error = ghl_initialize_hybrid_eos_functions_and_params(
        1.0, 0.0, 100.0, 1, one_rho_ppoly, one_Gamma_ppoly, 1.0, 2.0, &one_piece_eos);
  if(error != ghl_success || !isfinite(one_piece_eos.p_ppoly[0])
     || one_piece_eos.p_ppoly[0] != 0.0) {
    ghl_error("One-piece EOS pressure sentinel was not initialized to zero.\n");
  }

  const double two_rho_ppoly[1] = { 2.0 };
  const double two_Gamma_ppoly[2] = { 2.0, 3.0 };
  ghl_eos_parameters two_piece_eos = { 0 };
  error = ghl_initialize_hybrid_eos_functions_and_params(
        1.0, 0.0, 100.0, 2, two_rho_ppoly, two_Gamma_ppoly, 1.0, 2.0, &two_piece_eos);
  const double expected_left = pow(two_rho_ppoly[0], two_Gamma_ppoly[0]);
  const double expected_K1
        = pow(two_rho_ppoly[0], two_Gamma_ppoly[0] - two_Gamma_ppoly[1]);
  const double expected_right = expected_K1 * pow(two_rho_ppoly[0], two_Gamma_ppoly[1]);
  if(error != ghl_success || !isfinite(two_piece_eos.p_ppoly[0])
     || relative_error(expected_left, two_piece_eos.p_ppoly[0]) > 1e-14
     || relative_error(expected_left, expected_right) > 1e-14
     || two_piece_eos.p_ppoly[1] != 0.0) {
    ghl_error("Two-piece EOS pressure transition is invalid.\n");
  }
}
