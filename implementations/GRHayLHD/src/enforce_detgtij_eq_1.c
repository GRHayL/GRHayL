#include "GRHayLHD.h"

void GRHayLHD_enforce_detgtij_eq_1(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLHD_enforce_detgtij_eq_1;
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel for
  for(int k=0; k<cctkGH->cctk_lsh[2]; k++) {
    for(int j=0; j<cctkGH->cctk_lsh[1]; j++) {
      for(int i=0; i<cctkGH->cctk_lsh[0]; i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
        double gxx_physL=gxx[index];
        double gxy_physL=gxy[index];
        double gxz_physL=gxz[index];
        double gyy_physL=gyy[index];
        double gyz_physL=gyz[index];
        double gzz_physL=gzz[index];

        /**********************************************************************
         * Compute \tilde{\gamma_{ij}}, phi, and psi (BSSN) from g_{ij} (ADM) *
         **********************************************************************/
        double gijdet = gxx_physL * gyy_physL * gzz_physL + gxy_physL * gyz_physL * gxz_physL + gxz_physL * gxy_physL * gyz_physL
          - gxz_physL * gyy_physL * gxz_physL - gxy_physL * gxy_physL * gzz_physL - gxx_physL * gyz_physL * gyz_physL;

        gijdet = fabs(gijdet);

        const double psiL = exp( (1.0/12.0) * log(gijdet) );

        const double Psim4 = 1.0/(psiL*psiL*psiL*psiL);
        double gtxxL = gxx_physL*Psim4;
        double gtxyL = gxy_physL*Psim4;
        double gtxzL = gxz_physL*Psim4;
        double gtyyL = gyy_physL*Psim4;
        double gtyzL = gyz_physL*Psim4;
        double gtzzL = gzz_physL*Psim4;

        /*********************************
         * Apply det gtij = 1 constraint *
         *********************************/
        double gtijdet = gtxxL * gtyyL * gtzzL + gtxyL * gtyzL * gtxzL + gtxzL * gtxyL * gtyzL -
          gtxzL * gtyyL * gtxzL - gtxyL * gtxyL * gtzzL - gtxxL * gtyzL * gtyzL;

        double gtijdet_Fm1o3 = fabs(1.0/cbrt(gtijdet));

        gtxxL = gtxxL * gtijdet_Fm1o3;
        gtxyL = gtxyL * gtijdet_Fm1o3;
        gtxzL = gtxzL * gtijdet_Fm1o3;
        gtyyL = gtyyL * gtijdet_Fm1o3;
        gtyzL = gtyzL * gtijdet_Fm1o3;
        gtzzL = gtzzL * gtijdet_Fm1o3;

        if(gtijdet<0.0) { CCTK_VWarn(CCTK_WARN_ALERT,__LINE__, __FILE__, CCTK_THORNSTRING,
                                     "WARNING: det[3-metric]<0.0 at point  %d %d %d | cctk_lsh: %d %d %d. Hopefully this is occurring in gz's! gtij_phys = %.2e %.2e %.2e %.2e %.2e %.2e gtij_new = %.2e %.2e %.2e %.2e %.2e %.2e | gijdet = %.2e | gtijdet = %.2e",
				     i,j,k,cctkGH->cctk_lsh[0],cctkGH->cctk_lsh[1],cctkGH->cctk_lsh[2],gxx_physL,gxy_physL,gxz_physL,gyy_physL,gyz_physL,gzz_physL,gtxxL,gtxyL,gtxzL,gtyyL,gtyzL,gtzzL,-gijdet,gtijdet); }

        /*******************************************
         * Set the ADM gridfunctions to new values *
         *******************************************/
        const double Psi4 = psiL*psiL*psiL*psiL;
        gxx[index] = gtxxL*Psi4;
        gxy[index] = gtxyL*Psi4;
        gxz[index] = gtxzL*Psi4;
        gyy[index] = gtyyL*Psi4;
        gyz[index] = gtyzL*Psi4;
        gzz[index] = gtzzL*Psi4;
      }
    }
  }
}
