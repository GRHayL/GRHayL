#include "ghl_radiation.h"
#include <math.h>

/*
 * (c) Leo Werneck
 * Compute Fermi-Dirac integrals according to the approximations
 * in Takahashi, Eid, and Hillebrandt (1978)
 * https://adsabs.harvard.edu/pdf/1978A%26A....67..185T
 *
 * For all functions below, see Eqs. (A2) in the reference above
 * for eta > 1e-3, and Eqs. (A3) otherwise. For k=0, see Eq. (A4).
 */
GHL_DEVICE
double NRPyLeakage_Fermi_Dirac_integrals_k0(const double eta) {
  // Note: Eq. (A4) reads eta + log(1 + exp(-eta)). But we have:
  //
  // eta + log(1 + exp(-eta)) = eta + log( (1 + exp(eta))/exp(eta) )
  //                          = eta + log(1 + exp(eta)) - log(exp(eta))
  //                          = eta + log(1 + exp(eta)) - eta
  //                          = log(1 + exp(eta))
  return log1p(exp(eta)); // log(1.0 + exp(eta))
}

GHL_DEVICE
double NRPyLeakage_Fermi_Dirac_integrals_k1(const double eta) {
  if(eta > 1e-3) {
    const double eta2 = eta * eta;
    const double numer = eta2 / 2.0 + 1.6449;
    const double denom = 1.0 + exp(-1.6855 * eta);
    return numer / denom;
  }
  return exp(eta) / (1.0 + 0.2159 * exp(0.8857 * eta));
}

GHL_DEVICE
double NRPyLeakage_Fermi_Dirac_integrals_k2(const double eta) {
  if(eta > 1e-3) {
    const double eta3 = eta * eta * eta;
    const double numer = eta3 / 3.0 + 3.2899 * eta;
    const double denom = - expm1(-1.8246 * eta); // 1.0 - exp(-1.8246 * eta)
    return numer / denom;
  }
  return 2.0 * exp(eta) / (1.0 + 0.1092 * exp(0.8908 * eta));
}

GHL_DEVICE
double NRPyLeakage_Fermi_Dirac_integrals_k3(const double eta) {
  if(eta > 1e-3) {
    const double eta2 = eta * eta;
    const double eta4 = eta2 * eta2;
    const double numer = eta4 / 4.0 + 4.9348 * eta2 + 11.3644;
    const double denom = 1.0 + exp(-1.9039 * eta);
    return numer / denom;
  }
  return 6.0 * exp(eta) / (1.0 + 0.0559 * exp(0.9069 * eta));
}

GHL_DEVICE
double NRPyLeakage_Fermi_Dirac_integrals_k4(const double eta) {
  if(eta > 1e-3) {
    const double eta2 = eta * eta;
    const double eta3 = eta * eta2;
    const double eta5 = eta2 * eta3;
    const double numer = eta5 / 5.0 + 6.5797 * eta3 + 45.4576 * eta;
    const double denom = - expm1(-1.9484 * eta); // 1.0 - exp(-1.9484 * eta)
    return numer / denom;
  }
  return 24.0 * exp(eta) / (1.0 + 0.0287 * exp(0.9257 * eta));
}

GHL_DEVICE
double NRPyLeakage_Fermi_Dirac_integrals_k5(const double eta) {
  if(eta > 1e-3) {
    const double eta2 = eta * eta;
    const double eta4 = eta2 * eta2;
    const double eta6 = eta2 * eta4;
    const double numer = eta6 / 6.0 + 8.2247 * eta4 + 113.6439 * eta2 + 236.5323;
    const double denom = 1.0 + exp(-1.9727 * eta);
    return numer / denom;
  }
  return 120.0 * exp(eta) / (1.0 + 0.0147 * exp(0.9431 * eta));
}
