#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "NRPyLeakageET.h"

/*
 * (c) Leo Werneck
 * Compute GRMHD source terms following Ruffert et al. (1996)
 * https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
void NRPyLeakageET_compute_neutrino_luminosities(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(!NRPyLeakageET_ProcessOwnsData()) return;

  if(verbosity_level>1) CCTK_VINFO("Computing neutrino luminosities at ref. lvl. %d...",GetRefinementLevel(cctkGH));

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) {
    for(int j=0;j<cctk_lsh[1];j++) {
      for(int i=0;i<cctk_lsh[0];i++) {
        // Step 1: Set gridpoint index
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        const CCTK_REAL rhoL = rho[index];
        if( rhoL < rho_min_threshold || rhoL > rho_max_threshold ) {
          lum_nue[index] = lum_anue[index] = lum_nux[index] = 0.0;
        }
        else {
          const CCTK_REAL gxxL = gxx[index];
          const CCTK_REAL gxyL = gxy[index];
          const CCTK_REAL gxzL = gxz[index];
          const CCTK_REAL gyyL = gyy[index];
          const CCTK_REAL gyzL = gyz[index];
          const CCTK_REAL gzzL = gzz[index];
          const CCTK_REAL gdet = fabs(gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL
                                      - gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL);
          const CCTK_REAL phiL  = (1.0/12.0) * log(gdet);
          const CCTK_REAL psiL  = exp(phiL);
          const CCTK_REAL psi2L = psiL *psiL;
          const CCTK_REAL psi4L = psi2L*psi2L;
          const CCTK_REAL psi6L = psi4L*psi2L;
          if( psi6L > psi6_threshold ) {
            lum_nue[index] = lum_anue[index] = lum_nux[index] = 0.0;
          }
          else {
            // Step 2: Read from main memory
            const CCTK_REAL alpL         = alp[index];
            const CCTK_REAL Y_eL         = Y_e[index];
            const CCTK_REAL temperatureL = temperature[index];
            const CCTK_REAL wL           = w_lorentz[index];
            neutrino_optical_depths tauL;
            tauL.nue [0] = tau_0_nue [index];
            tauL.nue [1] = tau_1_nue [index];
            tauL.anue[0] = tau_0_anue[index];
            tauL.anue[1] = tau_1_anue[index];
            tauL.nux [0] = tau_0_nux [index];
            tauL.nux [1] = tau_1_nux [index];

            // Step 3: Compute neutrino luminosities
            neutrino_luminosities lumL;
            NRPyLeakage_compute_neutrino_luminosities(ghl_eos,
                                                      alpL, gxxL, gxyL, gxzL, gyyL, gyzL, gzzL,
                                                      rhoL, Y_eL, temperatureL, wL,
                                                      &tauL, &lumL);

            // Step 4: Write to main memory
            lum_nue [index] = lumL.nue;
            lum_anue[index] = lumL.anue;
            lum_nux [index] = lumL.nux;
          }
        }
      }
    }
  }
  if(verbosity_level>1) CCTK_VINFO("Finished computing neutrino luminosities at ref. lvl. %d...",GetRefinementLevel(cctkGH));
}
