#include "ghl.h"
#include "ghl_radiation.h"

void NRPyLeakage_compute_neutrino_opacities(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      const ghl_neutrino_optical_depths *restrict tau,
      ghl_neutrino_opacities *restrict kappa) {

  ghl_neutrino_opacities kappa_abs, kappa_scat;
  NRPyLeakage_compute_neutrino_absorption_and_scattering_opacities(
        eos, rho, Y_e, T, tau, &kappa_abs, &kappa_scat);

  for(int i = 0; i < 2; i++) {
    kappa->nue[i] = kappa_abs.nue[i] + kappa_scat.nue[i];
    if(!isfinite(kappa->nue[i])) {
      kappa->nue[i] = 1e-15;
    }

    kappa->anue[i] = kappa_abs.anue[i] + kappa_scat.anue[i];
    if(!isfinite(kappa->anue[i])) {
      kappa->anue[i] = 1e-15;
    }

    kappa->nux[i] = kappa_abs.nux[i] + kappa_scat.nux[i];
    if(!isfinite(kappa->nux[i])) {
      kappa->nux[i] = 1e-15;
    }
  }
}
