#include "Neutrinos.h"

/*
 * (c) Leo Werneck
 * Compute Fermi-Dirac integrals according to the approximations
 * in Takahashi, Eid, and Hillebrandt (1978)
 * https://adsabs.harvard.edu/pdf/1978A%26A....67..185T
 */
double NRPyLeakage_Fermi_Dirac_integrals(const int k, const double z) {


  double Fermi_Dirac_integral = 0.0;
  if(z>1e-3) {
    switch(k) {
    case(0):
      Fermi_Dirac_integral = 1.0*log(exp(z) + 1);
      break;
    case(1):
      Fermi_Dirac_integral = ((1.0/2.0)*((z)*(z)) + 1.6449)/(1 + exp(-1.6855*z));
      break;
    case(2):
      Fermi_Dirac_integral = ((1.0/3.0)*((z)*(z)*(z)) + 3.2898999999999998*z)/(1 - exp(-1.8246*z));
      break;
    case(3):
      Fermi_Dirac_integral = ((1.0/4.0)*((z)*(z)*(z)*(z)) + 4.9348000000000001*((z)*(z)) + 11.3644)/(1 + exp(-1.9038999999999999*z));
      break;
    case(4):
      Fermi_Dirac_integral = ((1.0/5.0)*((z)*(z)*(z)*(z)*(z)) + 6.5796999999999999*((z)*(z)*(z)) + 45.457599999999999*z)/(1 - exp(-1.9483999999999999*z));
      break;
    case(5):
      Fermi_Dirac_integral = ((1.0/6.0)*pow(z, 6) + 8.2247000000000003*((z)*(z)*(z)*(z)) + 113.6439*((z)*(z)) + 236.53229999999999)/(1 + exp(-1.9726999999999999*z));
      break;
    default:
      fprintf(stderr, "Unsuported value of k: %d\n", k);
      exit(1);
    }
  }
  else {
    switch(k) {
    case(0):
      Fermi_Dirac_integral = 1.0*log(exp(z) + 1);
      break;
    case(1):
      Fermi_Dirac_integral = exp(z)/(0.21590000000000001*exp(0.88570000000000004*z) + 1);
      break;
    case(2):
      Fermi_Dirac_integral = 2*exp(z)/(0.10920000000000001*exp(0.89080000000000004*z) + 1);
      break;
    case(3):
      Fermi_Dirac_integral = 6*exp(z)/(0.055899999999999998*exp(0.90690000000000004*z) + 1);
      break;
    case(4):
      Fermi_Dirac_integral = 24*exp(z)/(0.0287*exp(0.92569999999999997*z) + 1);
      break;
    case(5):
      Fermi_Dirac_integral = 120*exp(z)/(0.0147*exp(0.94310000000000005*z) + 1);
      break;
    default:
      fprintf(stderr, "Unsuported value of k: %d\n", k);
      exit(1);
    }
  }

  return Fermi_Dirac_integral;
}
