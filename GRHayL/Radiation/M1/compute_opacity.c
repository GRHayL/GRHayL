#include "ghl.h"
#include "ghl_radiation.h"

//kappa
int neutrino_opacity_transport(
      double rho,
      double T,
      double Y_e,
      double *k0x,
      double *k0y,
      double *k0z,
      double *k1x,
      double *k1y,
      double *k1z) {return 0;}

//kappa_a
int neutrino_absorbtions(
      double rho,
      double T,
      double Y_e,
      double *a0x,
      double *a0y,
      double *a0z,
      double *a1x,
      double *a1y,
      double *a1z) {return 0;}

//eta
int neutrino_emissions(
      double rho,
      double T,
      double Y_e,
      double *n0x,
      double *n0y,
      double *n0z,
      double *n1x,
      double *n1y,
      double *n1z) {return 0;}

int neutrino_weak_equilibrium(){return 0;}

int calc_neutrino_density(
      double rho,
      double T,
      double Y_e,
      double *d0x,
      double *d0y,
      double *d0z,
      double *d1x,
      double *d1y,
      double *d1z){return 0;}

//heads up note: I think 0-index refers to electron neutrino and antineutrino species, and 1 refers to heavier neutrino and antineutrino (muon and tauon) species taken together?
void calc_opacity_and_blackbodies(double rho, double T, double T_trapping_limit, double Y_e, double Y_e_trapping_limit, double tau_trap, double dt) {

  //transport opacity
  double kappa_0_local[3];
  double kappa_1_local[3];
  int ierr = neutrino_opacity_transport(rho, T, Y_e, &kappa_0_local[0], &kappa_0_local[1], &kappa_0_local[2], &kappa_1_local[0], &kappa_1_local[1], &kappa_1_local[2]);
  //finite-ness and error checks on kappa_local
  
  
  //Absorbtion Opacity 
  double abs_0_local[3];
  double abs_1_local[3];
  int ierr = neutrino_absorbtions(rho, T, Y_e, &abs_0_local[0], &abs_0_local[1], &abs_0_local[2], &abs_1_local[0], &abs_1_local[1], &abs_1_local[2]);
  //finite-ness and error checks on abs_1_local

 
  double eta_0_local[3];
  double eta_1_local[3];
  int ierr = neutrino_emissions(rho, T, Y_e, &eta_0_local[0], &eta_0_local[1], &eta_0_local[2], &eta_1_local[0], &eta_1_local[1], &eta_1_local[2]);
  //finite-ness and error checks on eta_local
  
  
  //Now find the optical depth from these local values
  double tau = min(sqrt(abs_1_local[0]*kappa_1_local[0]),sqrt(abs_1_local[1]*kappa_1_local[1]))*dt;

  //now to begin the blackbody calcs, assuming neutrino trapping.
  double trapped_0_density[3];
  double trapped_1_density[3];


}
