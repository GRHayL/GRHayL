#include "NRPyLeakageET.h"

static void set_opacity_struct_from_gfs(
      const int index,
      const CCTK_REAL *restrict kappa_0_nue,
      const CCTK_REAL *restrict kappa_1_nue,
      const CCTK_REAL *restrict kappa_0_anue,
      const CCTK_REAL *restrict kappa_1_anue,
      const CCTK_REAL *restrict kappa_0_nux,
      const CCTK_REAL *restrict kappa_1_nux,
      neutrino_opacities *restrict kappa) {
  kappa->nue [0] = kappa_0_nue [index];
  kappa->nue [1] = kappa_1_nue [index];
  kappa->anue[0] = kappa_0_anue[index];
  kappa->anue[1] = kappa_1_anue[index];
  kappa->nux [0] = kappa_0_nux [index];
  kappa->nux [1] = kappa_1_nux [index];
}

static void set_optical_depths_struct_from_gfs(
      const int index,
      const CCTK_REAL *restrict tau_0_nue,
      const CCTK_REAL *restrict tau_1_nue,
      const CCTK_REAL *restrict tau_0_anue,
      const CCTK_REAL *restrict tau_1_anue,
      const CCTK_REAL *restrict tau_0_nux,
      const CCTK_REAL *restrict tau_1_nux,
      neutrino_optical_depths *restrict tau) {
  tau->nue [0] = tau_0_nue [index];
  tau->nue [1] = tau_1_nue [index];
  tau->anue[0] = tau_0_anue[index];
  tau->anue[1] = tau_1_anue[index];
  tau->nux [0] = tau_0_nux [index];
  tau->nux [1] = tau_1_nux [index];
}

void NRPyLeakageET_optical_depths_PathOfLeastResistance(CCTK_ARGUMENTS) {
  
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(!NRPyLeakageET_ProcessOwnsData()) return;

  const CCTK_REAL dxx[3] = {CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1), CCTK_DELTA_SPACE(2)};

#pragma omp parallel
  for(int k=0;k<cctk_lsh[2];k++) {
    for(int j=0;j<cctk_lsh[1];j++) {
      for(int i=0;i<cctk_lsh[0];i++) {

        // Step 1: Set gridpoint indices
        const int i_j_k   = CCTK_GFINDEX3D(cctkGH, i  , j, k  );
        const int ip1_j_k = CCTK_GFINDEX3D(cctkGH, i+1, j, k  );
        const int im1_j_k = CCTK_GFINDEX3D(cctkGH, i-1, j, k  );
        const int i_jp1_k = CCTK_GFINDEX3D(cctkGH, i, j+1, k  );
        const int i_jm1_k = CCTK_GFINDEX3D(cctkGH, i, j-1, k  );
        const int i_j_kp1 = CCTK_GFINDEX3D(cctkGH, i, j  , k+1);
        const int i_j_km1 = CCTK_GFINDEX3D(cctkGH, i, j  , k-1);

        const CCTK_REAL rhoL = rho[i_j_k];
        // Only compute optical depths if density is above a threshold
        if( rhoL < rho_min_threshold || rhoL > rho_max_threshold ) {
          tau_0_nue [i_j_k] = 0.0;
          tau_1_nue [i_j_k] = 0.0;
          tau_0_anue[i_j_k] = 0.0;
          tau_1_anue[i_j_k] = 0.0;
          tau_0_nux [i_j_k] = 0.0;
          tau_1_nux [i_j_k] = 0.0;
        }
        else {
          const CCTK_REAL gxxL  = gxx[i_j_k];
          const CCTK_REAL gxyL  = gxy[i_j_k];
          const CCTK_REAL gxzL  = gxz[i_j_k];
          const CCTK_REAL gyyL  = gyy[i_j_k];
          const CCTK_REAL gyzL  = gyz[i_j_k];
          const CCTK_REAL gzzL  = gzz[i_j_k];
          const CCTK_REAL gdet  = fabs(gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL
                                     - gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL);
          const CCTK_REAL phiL  = (1.0/12.0) * log(gdet);
          const CCTK_REAL psiL  = exp(phiL);
          const CCTK_REAL psi2L = psiL *psiL;
          const CCTK_REAL psi4L = psi2L*psi2L;
          const CCTK_REAL psi6L = psi4L*psi2L;
          if( psi6L > psi6_threshold ) {
            tau_0_nue [i_j_k] = 0.0;
            tau_1_nue [i_j_k] = 0.0;
            tau_0_anue[i_j_k] = 0.0;
            tau_1_anue[i_j_k] = 0.0;
            tau_0_nux [i_j_k] = 0.0;
            tau_1_nux [i_j_k] = 0.0;
          }
          else {
            // Step 2: Read in metric gfs from main memory
            const CCTK_REAL stencil_gxx[3] = {gxx[im1_j_k], gxx[i_j_k], gxx[ip1_j_k]};
            const CCTK_REAL stencil_gyy[3] = {gyy[im1_j_k], gyy[i_j_k], gyy[ip1_j_k]};
            const CCTK_REAL stencil_gzz[3] = {gzz[im1_j_k], gzz[i_j_k], gzz[ip1_j_k]};

            // Step 3: Read in opacity gfs from main memory
            neutrino_opacities kappa_i_j_k;
            neutrino_opacities kappa_ip1_j_k, kappa_im1_j_k;
            neutrino_opacities kappa_i_jp1_k, kappa_i_jm1_k;
            neutrino_opacities kappa_i_j_kp1, kappa_i_j_km1;
            set_opacity_struct_from_gfs(i_j_k  , kappa_0_nue, kappa_1_nue, kappa_0_anue, kappa_1_anue, kappa_0_nux, kappa_1_nux, &kappa_i_j_k  );
            set_opacity_struct_from_gfs(ip1_j_k, kappa_0_nue, kappa_1_nue, kappa_0_anue, kappa_1_anue, kappa_0_nux, kappa_1_nux, &kappa_ip1_j_k);
            set_opacity_struct_from_gfs(im1_j_k, kappa_0_nue, kappa_1_nue, kappa_0_anue, kappa_1_anue, kappa_0_nux, kappa_1_nux, &kappa_im1_j_k);
            set_opacity_struct_from_gfs(i_jp1_k, kappa_0_nue, kappa_1_nue, kappa_0_anue, kappa_1_anue, kappa_0_nux, kappa_1_nux, &kappa_i_jp1_k);
            set_opacity_struct_from_gfs(i_jm1_k, kappa_0_nue, kappa_1_nue, kappa_0_anue, kappa_1_anue, kappa_0_nux, kappa_1_nux, &kappa_i_jm1_k);
            set_opacity_struct_from_gfs(i_j_kp1, kappa_0_nue, kappa_1_nue, kappa_0_anue, kappa_1_anue, kappa_0_nux, kappa_1_nux, &kappa_i_j_kp1);
            set_opacity_struct_from_gfs(i_j_km1, kappa_0_nue, kappa_1_nue, kappa_0_anue, kappa_1_anue, kappa_0_nux, kappa_1_nux, &kappa_i_j_km1);

            // Step 4: Read in optical depth gfs from main memory
            neutrino_optical_depths tau_ip1_j_k, tau_im1_j_k;
            neutrino_optical_depths tau_i_jp1_k, tau_i_jm1_k;
            neutrino_optical_depths tau_i_j_kp1, tau_i_j_km1;
            set_optical_depths_struct_from_gfs(ip1_j_k, tau_0_nue_p, tau_1_nue_p, tau_0_anue_p, tau_1_anue_p, tau_0_nux_p, tau_1_nux_p, &tau_ip1_j_k);
            set_optical_depths_struct_from_gfs(im1_j_k, tau_0_nue_p, tau_1_nue_p, tau_0_anue_p, tau_1_anue_p, tau_0_nux_p, tau_1_nux_p, &tau_im1_j_k);
            set_optical_depths_struct_from_gfs(i_jp1_k, tau_0_nue_p, tau_1_nue_p, tau_0_anue_p, tau_1_anue_p, tau_0_nux_p, tau_1_nux_p, &tau_i_jp1_k);
            set_optical_depths_struct_from_gfs(i_jm1_k, tau_0_nue_p, tau_1_nue_p, tau_0_anue_p, tau_1_anue_p, tau_0_nux_p, tau_1_nux_p, &tau_i_jm1_k);
            set_optical_depths_struct_from_gfs(i_j_kp1, tau_0_nue_p, tau_1_nue_p, tau_0_anue_p, tau_1_anue_p, tau_0_nux_p, tau_1_nux_p, &tau_i_j_kp1);
            set_optical_depths_struct_from_gfs(i_j_km1, tau_0_nue_p, tau_1_nue_p, tau_0_anue_p, tau_1_anue_p, tau_0_nux_p, tau_1_nux_p, &tau_i_j_km1);

            // Step 5: Compute the new optical depths
            neutrino_optical_depths tau_i_j_k;
            NRPyLeakage_optical_depths_PathOfLeastResistance(dxx, stencil_gxx, stencil_gyy, stencil_gzz,
                                                             &kappa_ip1_j_k, &kappa_im1_j_k,
                                                             &kappa_i_jp1_k, &kappa_i_jm1_k,
                                                             &kappa_i_j_kp1, &kappa_i_j_km1,
                                                             &tau_ip1_j_k  , &tau_im1_j_k,
                                                             &tau_i_jp1_k  , &tau_i_jm1_k,
                                                             &tau_i_j_kp1  , &tau_i_j_km1,
                                                             &kappa_i_j_k  , &tau_i_j_k);

            // Step 6: Write to main memory
            tau_0_nue [i_j_k] = tau_i_j_k.nue [0];
            tau_1_nue [i_j_k] = tau_i_j_k.nue [1];
            tau_0_anue[i_j_k] = tau_i_j_k.anue[0];
            tau_1_anue[i_j_k] = tau_i_j_k.anue[1];
            tau_0_nux [i_j_k] = tau_i_j_k.nux [0];
            tau_1_nux [i_j_k] = tau_i_j_k.nux [1];
          }
        }
      }
    }
  }
}
