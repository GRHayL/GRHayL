#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "NRPyLeakageET.h"

void NRPyLeakageET_CopyOpticalDepthsToAux(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(!NRPyLeakageET_ProcessOwnsData()) return;

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) {
    for(int j=0;j<cctk_lsh[1];j++) {
      for(int i=0;i<cctk_lsh[0];i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        tau_0_nue_aux [index] = tau_0_nue [index];
        tau_1_nue_aux [index] = tau_1_nue [index];
        tau_0_anue_aux[index] = tau_0_anue[index];
        tau_1_anue_aux[index] = tau_1_anue[index];
        tau_0_nux_aux [index] = tau_0_nux [index];
        tau_1_nux_aux [index] = tau_1_nux [index];
      }
    }
  }
}

static inline CCTK_REAL relative_difference(const CCTK_REAL a, const CCTK_REAL b) {
  if     (a!=0.0) return 1.0-b/a;
  else if(b!=0.0) return 1.0-a/b;
  else            return 0.0;
}

void NRPyLeakageET_compute_optical_depth_change(CCTK_ARGUMENTS, const int it) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(!NRPyLeakageET_ProcessOwnsData()) return;

  if(verbosity_level>1) CCTK_VINFO("Beginning to compute changes in optical depths at Ref. Lev. %d",GetRefinementLevel(cctkGH));

#pragma omp parallel for
  for(int k=cctk_nghostzones[2];k<cctk_lsh[2]-cctk_nghostzones[2];k++) {
    for(int j=cctk_nghostzones[1];j<cctk_lsh[1]-cctk_nghostzones[1];j++) {
      for(int i=cctk_nghostzones[0];i<cctk_lsh[0]-cctk_nghostzones[0];i++) {
        // Step 1: Set local index
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        const CCTK_REAL rhoL = rho[index];
        if( rhoL < rho_min_threshold || rhoL > rho_max_threshold ) {
          tau_0_nue_aux [index] = 0.0;
          tau_1_nue_aux [index] = 0.0;
          tau_0_anue_aux[index] = 0.0;
          tau_1_anue_aux[index] = 0.0;
          tau_0_nux_aux [index] = 0.0;
          tau_1_nux_aux [index] = 0.0;
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
            tau_0_nue_aux [index] = 0.0;
            tau_1_nue_aux [index] = 0.0;
            tau_0_anue_aux[index] = 0.0;
            tau_1_anue_aux[index] = 0.0;
            tau_0_nux_aux [index] = 0.0;
            tau_1_nux_aux [index] = 0.0;
          }
          else {
            // Step 2: Compute relative differences in tau
            const CCTK_REAL rel_diff_tau_0_nue  = relative_difference(tau_0_nue [index],tau_0_nue_aux [index]);
            const CCTK_REAL rel_diff_tau_1_nue  = relative_difference(tau_1_nue [index],tau_1_nue_aux [index]);
            const CCTK_REAL rel_diff_tau_0_anue = relative_difference(tau_0_anue[index],tau_0_anue_aux[index]);
            const CCTK_REAL rel_diff_tau_1_anue = relative_difference(tau_1_anue[index],tau_1_anue_aux[index]);
            const CCTK_REAL rel_diff_tau_0_nux  = relative_difference(tau_0_nux [index],tau_0_nux_aux [index]);
            const CCTK_REAL rel_diff_tau_1_nux  = relative_difference(tau_1_nux [index],tau_1_nux_aux [index]);

            // Step 3: Write to main memory
            tau_0_nue_aux [index] = rel_diff_tau_0_nue *rel_diff_tau_0_nue;
            tau_1_nue_aux [index] = rel_diff_tau_1_nue *rel_diff_tau_1_nue;
            tau_0_anue_aux[index] = rel_diff_tau_0_anue*rel_diff_tau_0_anue;
            tau_1_anue_aux[index] = rel_diff_tau_1_anue*rel_diff_tau_1_anue;
            tau_0_nux_aux [index] = rel_diff_tau_0_nux *rel_diff_tau_0_nux;
            tau_1_nux_aux [index] = rel_diff_tau_1_nux *rel_diff_tau_1_nux;
          }
        }
      }
    }
  }

  if(verbosity_level>1) CCTK_VINFO("Finished computing changes in optical depths at Ref. Lev. %d",GetRefinementLevel(cctkGH));
}
