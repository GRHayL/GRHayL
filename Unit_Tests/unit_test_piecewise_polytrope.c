#include "unit_tests.h"

int main(int argc, char **argv) {

  // This is the SLy EOS with the values taken from the NRPyEOS python code

  const int neos = 4;
  const double W_max = 10.0; //IGM default
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300; //IGM default
  const double Gamma_th = 2.0; //Taken from magnetizedTOV.par
  const double rho_ppoly[3] = {2.44034e+07,3.78358e+11,2.62780e+12};
  const double Gamma_ppoly[4] = {1.58425,1.28733,0.62223,1.35692};
  const double k_ppoly0 = 6.80110e-9;

  eos_parameters eos;
  initialize_hybrid_eos_functions_and_params(W_max,
                                             rho_b_min, rho_b_min, rho_b_max,
                                             neos, rho_ppoly, Gamma_ppoly,
                                             k_ppoly0, Gamma_th, &eos);

  // Expected output values taken from NRPyEOS python code
  double k_comp[4] = {6.8010999999999996e-09, 1.0618444833278535e-06, 5.3275084583961501e+01, 3.9992069172431910e-08};
  double eps_comp[4] = {0.0, -2.490902788423614614e-04, 1.355533308962642014e-02, 7.652227660602063664e-03};
  double k_pert[4], eps_pert[4];
  for(int i=1; i<4; i++) {
    k_pert[i] = k_comp[i]*(1.0 + randf(-1,1)*1e-14);
    eps_pert[i] = eps_comp[i]*(1.0 + randf(-1,1)*1e-14);
  }

  for(int i=1; i<4; i++) {
    if(validate(k_comp[i], eos.K_ppoly[i], k_pert[i]))
      grhayl_error("unit_test_piecewise_polytrope has failed for K_ppoly.\n"
                   "For index %d, expected %e, computed %e, perturbed %e\n"
                   "%e\n", i, k_comp[i], eos.K_ppoly[i], k_pert[i],
                           validate(k_comp[i], eos.K_ppoly[i], k_pert[i]));

    if(validate(eps_comp[i], eos.eps_integ_const[i], eps_pert[i]))
      grhayl_error("unit_test_piecewise_polytrope has failed for eps_integ_const.\n"
                   "For index %d, expected %e, computed %e, perturbed %e\n"
                   "relative error: %e\n", i, eps_comp[i], eos.eps_integ_const[i], eps_pert[i],
                           validate(eps_comp[i], eos.eps_integ_const[i], eps_pert[i]));
  }
}
