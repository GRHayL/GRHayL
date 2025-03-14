#include "ghl_unit_tests.h"
#include "ghl.h"
#include "ghl_radiation.h"

int main(){

  double rN_e0 = 2./6.;
  double rN_ae0 = 1./6.;
  double rN_x0 = 2./6.;
  double rN_ax0 = 1./6.;
  double rE_e0 = 1.;
  double rE_ae0 = 1.;
  double rE_x0 = 2.;
  double rE_ax0 = 2.;

  ghl_radiation_flux_vector rF_e0;
  rF_e0.D[0] = 1.0;
  rF_e0.D[1] = 1.0;
  rF_e0.D[2] = 1.0;
  rF_e0.D[3] = 1.0;
  ghl_radiation_flux_vector rF_ae0; 
  rF_ae0.D[0] = 1.0;
  rF_ae0.D[1] = 1.0;
  rF_ae0.D[2] = 1.0;
  rF_ae0.D[3] = 1.0;
  ghl_radiation_flux_vector rF_x0;
  rF_x0.D[0] = 2.0;
  rF_x0.D[1] = 2.0;
  rF_x0.D[2] = 2.0;
  rF_x0.D[3] = 2.0;
  ghl_radiation_flux_vector rF_ax0;
  rF_ax0.D[0] = 2.0;
  rF_ax0.D[1] = 2.0;
  rF_ax0.D[2] = 2.0;
  rF_ax0.D[3] = 2.0;

  double Nfloor = 0.0;
  bool mix_type = 0;

  flavor_mix_adjustment(&rN_e0, &rN_ae0, &rN_x0, &rN_ax0, Nfloor,
                        &rE_e0, &rE_ae0, &rE_x0, &rE_ax0,
                        &rF_e0, &rF_ae0, &rF_x0, &rF_ax0, mix_type);
  
  printf("Mix Type: Equilibrium\n\n");
  printf("Number:\nElectron neu. = %e\nAntiElection neu. = %e\nHeavy neu. = %e\nAnitHeavy neu. = %e\n\n", rN_e0, rN_ae0, rN_x0, rN_ax0);
  printf("Energy:\nElectron neu. = %e\nAntiElection neu. = %e\nHeavy neu. = %e\nAnitHeavy neu. = %e\n\n", rE_e0, rE_ae0, rE_x0, rE_ax0);
  printf("Flux (0,1,2,3)):\nElectron neu. = %e, %e, %e, %e\nAntiElection neu. = %e, %e, %e, %e\nHeavy neu. = %e, %e, %e, %e\nAnitHeavy neu. = %e, %e, %e, %e\n\n", 
         rF_e0.D[0], rF_e0.D[1], rF_e0.D[2], rF_e0.D[3],
         rF_ae0.D[0], rF_ae0.D[1], rF_ae0.D[2], rF_ae0.D[3],
         rF_x0.D[0], rF_x0.D[1], rF_x0.D[2], rF_x0.D[3],
         rF_ax0.D[0], rF_ax0.D[1], rF_ax0.D[2], rF_ax0.D[3]);


  double rN_e1 = 2./6.;
  double rN_ae1 = 1./6.;
  double rN_x1 = 2./6.;
  double rN_ax1 = 1./6.;
  double rE_e1 = 1.0;
  double rE_ae1 = 1.0;
  double rE_x1 = 2.0;
  double rE_ax1 = 2.0;

  ghl_radiation_flux_vector rF_e1;
  rF_e1.D[0] = 1.0;
  rF_e1.D[1] = 1.0;
  rF_e1.D[2] = 1.0;
  rF_e1.D[3] = 1.0;
  ghl_radiation_flux_vector rF_ae1; 
  rF_ae1.D[0] = 1.0;
  rF_ae1.D[1] = 1.0;
  rF_ae1.D[2] = 1.0;
  rF_ae1.D[3] = 1.0;
  ghl_radiation_flux_vector rF_x1;
  rF_x1.D[0] = 2.0;
  rF_x1.D[1] = 2.0;
  rF_x1.D[2] = 2.0;
  rF_x1.D[3] = 2.0;
  ghl_radiation_flux_vector rF_ax1;
  rF_ax1.D[0] = 2.0;
  rF_ax1.D[1] = 2.0;
  rF_ax1.D[2] = 2.0;
  rF_ax1.D[3] = 2.0;

  mix_type = 1;

  flavor_mix_adjustment(&rN_e1, &rN_ae1, &rN_x1, &rN_ax1, Nfloor,
                        &rE_e1, &rE_ae1, &rE_x1, &rE_ax1,
                        &rF_e1, &rF_ae1, &rF_x1, &rF_ax1, mix_type);
  
  printf("Mix Type: Maximal\n\n");
  printf("Number:\nElectron neu. = %e\nAntiElection neu. = %e\nHeavy neu. = %e\nAnitHeavy neu. = %e\n\n", rN_e1, rN_ae1, rN_x1, rN_ax1);
  printf("Energy:\nElectron neu. = %e\nAntiElection neu. = %e\nHeavy neu. = %e\nAnitHeavy neu. = %e\n\n", rE_e1, rE_ae1, rE_x1, rE_ax1);
  printf("Flux (0,1,2,3)):\nElectron neu. = %e, %e, %e, %e\nAntiElection neu. = %e, %e, %e, %e\nHeavy neu. = %e, %e, %e, %e\nAnitHeavy neu. = %e, %e, %e, %e\n\n", 
         rF_e1.D[0], rF_e1.D[1], rF_e1.D[2], rF_e1.D[3],
         rF_ae1.D[0], rF_ae1.D[1], rF_ae1.D[2], rF_ae1.D[3],
         rF_x1.D[0], rF_x1.D[1], rF_x1.D[2], rF_x1.D[3],
         rF_ax1.D[0], rF_ax1.D[1], rF_ax1.D[2], rF_ax1.D[3]);

  return 1;

}
