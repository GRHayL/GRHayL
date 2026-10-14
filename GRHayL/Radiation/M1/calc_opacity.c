#include "ghl.h"
#include "ghl_radiation.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

// This file translates from thc_M1_calc_opacity.cc
// The goal of this function is to calculate the flux term d/dx^i F_i and store it in rhs
// Currently only implementing this for grid-based codes
void ghl_M1_calc_opacity(
    const ghl_eos_parameters *restrict eos,
    const double rho,
    const double Y_e,
    const double T,
    double cdt
    //const ghl_neutrino_optical_depths *restrict tau,
    ){
    #define NGROUPS 1;
    #define NSPECIES 3;
    // Get the transport opacity (absorption + scattering)
    
    // setting tau to infinity?
    ghl_neutrino_optical_depths tau = {0};
    double tau_inf = 1e10;
    for (int i=0;i<2;i++) {
        tau.nue[i]  = tau_inf;
        tau.anue[i] = tau_inf;
        tau.nux[i]  = tau_inf;
    }
    // kappa      -> kappa_0_loc, kappa_1_loc in thc
    // kappa_abs  -> abs_0_loc, abs_1_loc in thc
    // kappa_scat -> scat_0_loc scat_1_loc in thc
    ghl_neutrino_opacities kappa = {0};
    ghl_neutrino_opacities kappa_abs = {0};
    ghl_neutrino_opacities kappa_scat = {0};
 
    NRPyLeakage_compute_neutrino_opacities(eos, rho, Y_e, T, &tau, &kappa);
    NRPyLeakage_compute_neutrino_absorption_and_scattering_opacities(
      eos, rho, Y_e, T, &tau, &kappa_abs, &kappa_scat);

    // assert(!ierr);
    // assert(isfinite(abs_0_loc[0]));
    // assert(isfinite(abs_0_loc[1]));
    // assert(isfinite(abs_0_loc[2]));
    // assert(isfinite(abs_1_loc[0]));
    // assert(isfinite(abs_1_loc[1]));
    // assert(isfinite(abs_1_loc[2]));

    // An effective optical depth used to decide whether to compute
    // the black body function for neutrinos assuming neutrino trapping
    // or at a fixed temperature and Ye
    tau = min(sqrt(abs_1_loc[0]*kappa_1_loc[0]),
            sqrt(abs_1_loc[1]*kappa_1_loc[1]))*cdt;
    // Compute the neutrino black body functions assuming trapped neutrinos
    double nudens_0_trap[3], nudens_1_trap[3];
    
    return GSL_SUCCESS;

} // ghl_M1_calc_opacity