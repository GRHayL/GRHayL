
#include "Neutrinos.h"

#define MINUS1 0
#define PLUS0  1
#define PLUS1  2

/*
 * (c) Leo Werneck
 * Compute GRMHD source terms following Ruffert et al. (1996)
 * https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
void NRPyLeakage_optical_depths_PathOfLeastResistance(
      const double *restrict dxx,
      const double *restrict stencil_gxx,
      const double *restrict stencil_gyy,
      const double *restrict stencil_gzz,
      const neutrino_opacities *restrict kappa_i_j_k,
      const neutrino_opacities *restrict kappa_im1_j_k,
      const neutrino_opacities *restrict kappa_ip1_j_k,
      const neutrino_opacities *restrict kappa_i_jm1_k,
      const neutrino_opacities *restrict kappa_i_jp1_k,
      const neutrino_opacities *restrict kappa_i_j_km1,
      const neutrino_opacities *restrict kappa_i_j_kp1,
      const neutrino_opacities *restrict tau_im1_j_k,
      const neutrino_opacities *restrict tau_ip1_j_k,
      const neutrino_opacities *restrict tau_i_jm1_k,
      const neutrino_opacities *restrict tau_i_jp1_k,
      const neutrino_opacities *restrict tau_i_j_km1,
      const neutrino_opacities *restrict tau_i_j_kp1,
      neutrino_optical_depths *restrict tau_i_j_k) {

  // Step 4: Compute metric at cell faces
  const double gxx_iphalf_j_k = 0.5*(stencil_gxx[PLUS0] + stencil_gxx[PLUS1] );
  const double gxx_imhalf_j_k = 0.5*(stencil_gxx[PLUS0] + stencil_gxx[MINUS1]);

  const double gyy_i_jphalf_k = 0.5*(stencil_gyy[PLUS0] + stencil_gyy[PLUS1] );
  const double gyy_i_jmhalf_k = 0.5*(stencil_gyy[PLUS0] + stencil_gyy[MINUS1]);

  const double gzz_i_j_kphalf = 0.5*(stencil_gzz[PLUS0] + stencil_gzz[PLUS1] );
  const double gzz_i_j_kmhalf = 0.5*(stencil_gzz[PLUS0] + stencil_gzz[MINUS1]);

  // Step 5: Compute ds^{i} = sqrt(gamma_{ii}dx^{i}dx^{i})
  const double ds_iphalf_j_k = sqrt(dxx[0]*dxx[0]*gxx_iphalf_j_k);
  const double ds_imhalf_j_k = sqrt(dxx[0]*dxx[0]*gxx_imhalf_j_k);
  const double ds_i_jphalf_k = sqrt(dxx[1]*dxx[1]*gyy_i_jphalf_k);
  const double ds_i_jmhalf_k = sqrt(dxx[1]*dxx[1]*gyy_i_jmhalf_k);
  const double ds_i_j_kphalf = sqrt(dxx[2]*dxx[2]*gzz_i_j_kphalf);
  const double ds_i_j_kmhalf = sqrt(dxx[2]*dxx[2]*gzz_i_j_kmhalf);

  // Step 6: Compute opacities at cell faces
  const double kappa_0_nue_iphalf_j_k = 0.5*(kappa_i_j_k->nue[0] + kappa_ip1_j_k->nue[0]);
  const double kappa_0_nue_imhalf_j_k = 0.5*(kappa_i_j_k->nue[0] + kappa_im1_j_k->nue[0]);
  const double kappa_0_nue_i_jphalf_k = 0.5*(kappa_i_j_k->nue[0] + kappa_i_jp1_k->nue[0]);
  const double kappa_0_nue_i_jmhalf_k = 0.5*(kappa_i_j_k->nue[0] + kappa_i_jm1_k->nue[0]);
  const double kappa_0_nue_i_j_kphalf = 0.5*(kappa_i_j_k->nue[0] + kappa_i_j_kp1->nue[0]);
  const double kappa_0_nue_i_j_kmhalf = 0.5*(kappa_i_j_k->nue[0] + kappa_i_j_km1->nue[0]);

  const double kappa_1_nue_iphalf_j_k = 0.5*(kappa_i_j_k->nue[1] + kappa_ip1_j_k->nue[1]);
  const double kappa_1_nue_imhalf_j_k = 0.5*(kappa_i_j_k->nue[1] + kappa_im1_j_k->nue[1]);
  const double kappa_1_nue_i_jphalf_k = 0.5*(kappa_i_j_k->nue[1] + kappa_i_jp1_k->nue[1]);
  const double kappa_1_nue_i_jmhalf_k = 0.5*(kappa_i_j_k->nue[1] + kappa_i_jm1_k->nue[1]);
  const double kappa_1_nue_i_j_kphalf = 0.5*(kappa_i_j_k->nue[1] + kappa_i_j_kp1->nue[1]);
  const double kappa_1_nue_i_j_kmhalf = 0.5*(kappa_i_j_k->nue[1] + kappa_i_j_km1->nue[1]);

  const double kappa_0_anue_iphalf_j_k = 0.5*(kappa_i_j_k->anue[0] + kappa_ip1_j_k->anue[0]);
  const double kappa_0_anue_imhalf_j_k = 0.5*(kappa_i_j_k->anue[0] + kappa_im1_j_k->anue[0]);
  const double kappa_0_anue_i_jphalf_k = 0.5*(kappa_i_j_k->anue[0] + kappa_i_jp1_k->anue[0]);
  const double kappa_0_anue_i_jmhalf_k = 0.5*(kappa_i_j_k->anue[0] + kappa_i_jm1_k->anue[0]);
  const double kappa_0_anue_i_j_kphalf = 0.5*(kappa_i_j_k->anue[0] + kappa_i_j_kp1->anue[0]);
  const double kappa_0_anue_i_j_kmhalf = 0.5*(kappa_i_j_k->anue[0] + kappa_i_j_km1->anue[0]);

  const double kappa_1_anue_iphalf_j_k = 0.5*(kappa_i_j_k->anue[1] + kappa_ip1_j_k->anue[1]);
  const double kappa_1_anue_imhalf_j_k = 0.5*(kappa_i_j_k->anue[1] + kappa_im1_j_k->anue[1]);
  const double kappa_1_anue_i_jphalf_k = 0.5*(kappa_i_j_k->anue[1] + kappa_i_jp1_k->anue[1]);
  const double kappa_1_anue_i_jmhalf_k = 0.5*(kappa_i_j_k->anue[1] + kappa_i_jm1_k->anue[1]);
  const double kappa_1_anue_i_j_kphalf = 0.5*(kappa_i_j_k->anue[1] + kappa_i_j_kp1->anue[1]);
  const double kappa_1_anue_i_j_kmhalf = 0.5*(kappa_i_j_k->anue[1] + kappa_i_j_km1->anue[1]);

  const double kappa_0_nux_iphalf_j_k = 0.5*(kappa_i_j_k->nux[0] + kappa_ip1_j_k->nux[0]);
  const double kappa_0_nux_imhalf_j_k = 0.5*(kappa_i_j_k->nux[0] + kappa_im1_j_k->nux[0]);
  const double kappa_0_nux_i_jphalf_k = 0.5*(kappa_i_j_k->nux[0] + kappa_i_jp1_k->nux[0]);
  const double kappa_0_nux_i_jmhalf_k = 0.5*(kappa_i_j_k->nux[0] + kappa_i_jm1_k->nux[0]);
  const double kappa_0_nux_i_j_kphalf = 0.5*(kappa_i_j_k->nux[0] + kappa_i_j_kp1->nux[0]);
  const double kappa_0_nux_i_j_kmhalf = 0.5*(kappa_i_j_k->nux[0] + kappa_i_j_km1->nux[0]);

  const double kappa_1_nux_iphalf_j_k = 0.5*(kappa_i_j_k->nux[1] + kappa_ip1_j_k->nux[1]);
  const double kappa_1_nux_imhalf_j_k = 0.5*(kappa_i_j_k->nux[1] + kappa_im1_j_k->nux[1]);
  const double kappa_1_nux_i_jphalf_k = 0.5*(kappa_i_j_k->nux[1] + kappa_i_jp1_k->nux[1]);
  const double kappa_1_nux_i_jmhalf_k = 0.5*(kappa_i_j_k->nux[1] + kappa_i_jm1_k->nux[1]);
  const double kappa_1_nux_i_j_kphalf = 0.5*(kappa_i_j_k->nux[1] + kappa_i_j_kp1->nux[1]);
  const double kappa_1_nux_i_j_kmhalf = 0.5*(kappa_i_j_k->nux[1] + kappa_i_j_km1->nux[1]);

  // Step 7: Compute optical depth at neighboring points
  const double tau_0_nue_ip1_j_k = tau_ip1_j_k->nue[0] + ds_iphalf_j_k*kappa_0_nue_iphalf_j_k;
  const double tau_0_nue_im1_j_k = tau_im1_j_k->nue[0] + ds_imhalf_j_k*kappa_0_nue_imhalf_j_k;
  const double tau_0_nue_i_jp1_k = tau_i_jp1_k->nue[0] + ds_i_jphalf_k*kappa_0_nue_i_jphalf_k;
  const double tau_0_nue_i_jm1_k = tau_i_jm1_k->nue[0] + ds_i_jmhalf_k*kappa_0_nue_i_jmhalf_k;
  const double tau_0_nue_i_j_kp1 = tau_i_j_kp1->nue[0] + ds_i_j_kphalf*kappa_0_nue_i_j_kphalf;
  const double tau_0_nue_i_j_km1 = tau_i_j_km1->nue[0] + ds_i_j_kmhalf*kappa_0_nue_i_j_kmhalf;

  const double tau_1_nue_ip1_j_k = tau_ip1_j_k->nue[1] + ds_iphalf_j_k*kappa_1_nue_iphalf_j_k;
  const double tau_1_nue_im1_j_k = tau_im1_j_k->nue[1] + ds_imhalf_j_k*kappa_1_nue_imhalf_j_k;
  const double tau_1_nue_i_jp1_k = tau_i_jp1_k->nue[1] + ds_i_jphalf_k*kappa_1_nue_i_jphalf_k;
  const double tau_1_nue_i_jm1_k = tau_i_jm1_k->nue[1] + ds_i_jmhalf_k*kappa_1_nue_i_jmhalf_k;
  const double tau_1_nue_i_j_kp1 = tau_i_j_kp1->nue[1] + ds_i_j_kphalf*kappa_1_nue_i_j_kphalf;
  const double tau_1_nue_i_j_km1 = tau_i_j_km1->nue[1] + ds_i_j_kmhalf*kappa_1_nue_i_j_kmhalf;

  const double tau_0_anue_ip1_j_k = tau_ip1_j_k->anue[0] + ds_iphalf_j_k*kappa_0_anue_iphalf_j_k;
  const double tau_0_anue_im1_j_k = tau_im1_j_k->anue[0] + ds_imhalf_j_k*kappa_0_anue_imhalf_j_k;
  const double tau_0_anue_i_jp1_k = tau_i_jp1_k->anue[0] + ds_i_jphalf_k*kappa_0_anue_i_jphalf_k;
  const double tau_0_anue_i_jm1_k = tau_i_jm1_k->anue[0] + ds_i_jmhalf_k*kappa_0_anue_i_jmhalf_k;
  const double tau_0_anue_i_j_kp1 = tau_i_j_kp1->anue[0] + ds_i_j_kphalf*kappa_0_anue_i_j_kphalf;
  const double tau_0_anue_i_j_km1 = tau_i_j_km1->anue[0] + ds_i_j_kmhalf*kappa_0_anue_i_j_kmhalf;

  const double tau_1_anue_ip1_j_k = tau_ip1_j_k->anue[1] + ds_iphalf_j_k*kappa_1_anue_iphalf_j_k;
  const double tau_1_anue_im1_j_k = tau_im1_j_k->anue[1] + ds_imhalf_j_k*kappa_1_anue_imhalf_j_k;
  const double tau_1_anue_i_jp1_k = tau_i_jp1_k->anue[1] + ds_i_jphalf_k*kappa_1_anue_i_jphalf_k;
  const double tau_1_anue_i_jm1_k = tau_i_jm1_k->anue[1] + ds_i_jmhalf_k*kappa_1_anue_i_jmhalf_k;
  const double tau_1_anue_i_j_kp1 = tau_i_j_kp1->anue[1] + ds_i_j_kphalf*kappa_1_anue_i_j_kphalf;
  const double tau_1_anue_i_j_km1 = tau_i_j_km1->anue[1] + ds_i_j_kmhalf*kappa_1_anue_i_j_kmhalf;

  const double tau_0_nux_ip1_j_k = tau_ip1_j_k->nux[0] + ds_iphalf_j_k*kappa_0_nux_iphalf_j_k;
  const double tau_0_nux_im1_j_k = tau_im1_j_k->nux[0] + ds_imhalf_j_k*kappa_0_nux_imhalf_j_k;
  const double tau_0_nux_i_jp1_k = tau_i_jp1_k->nux[0] + ds_i_jphalf_k*kappa_0_nux_i_jphalf_k;
  const double tau_0_nux_i_jm1_k = tau_i_jm1_k->nux[0] + ds_i_jmhalf_k*kappa_0_nux_i_jmhalf_k;
  const double tau_0_nux_i_j_kp1 = tau_i_j_kp1->nux[0] + ds_i_j_kphalf*kappa_0_nux_i_j_kphalf;
  const double tau_0_nux_i_j_km1 = tau_i_j_km1->nux[0] + ds_i_j_kmhalf*kappa_0_nux_i_j_kmhalf;

  const double tau_1_nux_ip1_j_k = tau_ip1_j_k->nux[1] + ds_iphalf_j_k*kappa_1_nux_iphalf_j_k;
  const double tau_1_nux_im1_j_k = tau_im1_j_k->nux[1] + ds_imhalf_j_k*kappa_1_nux_imhalf_j_k;
  const double tau_1_nux_i_jp1_k = tau_i_jp1_k->nux[1] + ds_i_jphalf_k*kappa_1_nux_i_jphalf_k;
  const double tau_1_nux_i_jm1_k = tau_i_jm1_k->nux[1] + ds_i_jmhalf_k*kappa_1_nux_i_jmhalf_k;
  const double tau_1_nux_i_j_kp1 = tau_i_j_kp1->nux[1] + ds_i_j_kphalf*kappa_1_nux_i_j_kphalf;
  const double tau_1_nux_i_j_km1 = tau_i_j_km1->nux[1] + ds_i_j_kmhalf*kappa_1_nux_i_j_kmhalf;

  // Step 8: Select path of least resistance
  const double new_tau_0_nue_i_j_k  = MIN(MIN(MIN(MIN(MIN(tau_0_nue_ip1_j_k,tau_0_nue_im1_j_k),tau_0_nue_i_jp1_k),tau_0_nue_i_jm1_k),tau_0_nue_i_j_kp1),tau_0_nue_i_j_km1);
  const double new_tau_1_nue_i_j_k  = MIN(MIN(MIN(MIN(MIN(tau_1_nue_ip1_j_k,tau_1_nue_im1_j_k),tau_1_nue_i_jp1_k),tau_1_nue_i_jm1_k),tau_1_nue_i_j_kp1),tau_1_nue_i_j_km1);
  const double new_tau_0_anue_i_j_k = MIN(MIN(MIN(MIN(MIN(tau_0_anue_ip1_j_k,tau_0_anue_im1_j_k),tau_0_anue_i_jp1_k),tau_0_anue_i_jm1_k),tau_0_anue_i_j_kp1),tau_0_anue_i_j_km1);
  const double new_tau_1_anue_i_j_k = MIN(MIN(MIN(MIN(MIN(tau_1_anue_ip1_j_k,tau_1_anue_im1_j_k),tau_1_anue_i_jp1_k),tau_1_anue_i_jm1_k),tau_1_anue_i_j_kp1),tau_1_anue_i_j_km1);
  const double new_tau_0_nux_i_j_k  = MIN(MIN(MIN(MIN(MIN(tau_0_nux_ip1_j_k,tau_0_nux_im1_j_k),tau_0_nux_i_jp1_k),tau_0_nux_i_jm1_k),tau_0_nux_i_j_kp1),tau_0_nux_i_j_km1);
  const double new_tau_1_nux_i_j_k  = MIN(MIN(MIN(MIN(MIN(tau_1_nux_ip1_j_k,tau_1_nux_im1_j_k),tau_1_nux_i_jp1_k),tau_1_nux_i_jm1_k),tau_1_nux_i_j_kp1),tau_1_nux_i_j_km1);

  // Step 9: Write results
  tau_i_j_k->nue[0]  = new_tau_0_nue_i_j_k;
  tau_i_j_k->nue[1]  = new_tau_1_nue_i_j_k;
  tau_i_j_k->anue[0] = new_tau_0_anue_i_j_k;
  tau_i_j_k->anue[1] = new_tau_1_anue_i_j_k;
  tau_i_j_k->nux[0]  = new_tau_0_nux_i_j_k;
  tau_i_j_k->nux[1]  = new_tau_1_nux_i_j_k;
}
