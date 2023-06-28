#include "GRHayLHDX.h"

extern "C" void GRHayLHDX_compute_ccc_centered_spacetime_quantities(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_compute_ccc_centered_spacetime_quantities;
  DECLARE_CCTK_PARAMETERS;

  constexpr std::array<int, Loop::dim> indextype = {1, 1, 1};
  const Loop::GF3D2layout layout(cctkGH, indextype);

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);

    CCTK_REAL lapse_avg = 0.0;
    CCTK_REAL betax_avg = 0.0;
    CCTK_REAL betay_avg = 0.0;
    CCTK_REAL betaz_avg = 0.0;
    CCTK_REAL gxx_avg   = 0.0;
    CCTK_REAL gxy_avg   = 0.0;
    CCTK_REAL gxz_avg   = 0.0;
    CCTK_REAL gyy_avg   = 0.0;
    CCTK_REAL gyz_avg   = 0.0;
    CCTK_REAL gzz_avg   = 0.0;

    for (int k=0; k<2; k++) { 
      for (int j=0; j<2; j++) { 
        for (int i=0; i<2; i++) {
          lapse_avg += alp(p.I + i*p.DI[0] + j*p.DI[1] + k*p.DI[2]);
          betax_avg += betax(p.I + i*p.DI[0] + j*p.DI[1] + k*p.DI[2]);
          betay_avg += betay(p.I + i*p.DI[0] + j*p.DI[1] + k*p.DI[2]);
          betaz_avg += betaz(p.I + i*p.DI[0] + j*p.DI[1] + k*p.DI[2]);
          gxx_avg += gxx(p.I + i*p.DI[0] + j*p.DI[1] + k*p.DI[2]);
          gxy_avg += gxy(p.I + i*p.DI[0] + j*p.DI[1] + k*p.DI[2]);
          gxz_avg += gxz(p.I + i*p.DI[0] + j*p.DI[1] + k*p.DI[2]);
          gyy_avg += gyy(p.I + i*p.DI[0] + j*p.DI[1] + k*p.DI[2]);
          gyz_avg += gyz(p.I + i*p.DI[0] + j*p.DI[1] + k*p.DI[2]);
          gzz_avg += gzz(p.I + i*p.DI[0] + j*p.DI[1] + k*p.DI[2]);
        }
      }
    }

    ccc_lapse(index) = lapse_avg/8.0;
    ccc_betax(index) = betax_avg/8.0;
    ccc_betay(index) = betay_avg/8.0;
    ccc_betaz(index) = betaz_avg/8.0;

    // Now we enforce that the conformal metric determinant is one
    const CCTK_REAL gxx_physL = gxx_avg/8.0;
    const CCTK_REAL gxy_physL = gxy_avg/8.0;
    const CCTK_REAL gxz_physL = gxz_avg/8.0;
    const CCTK_REAL gyy_physL = gyy_avg/8.0;
    const CCTK_REAL gyz_physL = gyz_avg/8.0;
    const CCTK_REAL gzz_physL = gzz_avg/8.0;

    /**********************************************************************
     * Compute \tilde{\gamma_{ij}}, phi, and psi (BSSN) from g_{ij} (ADM) *
     **********************************************************************/
    const CCTK_REAL gijdet = fabs(gxx_physL * gyy_physL * gzz_physL
                                + gxy_physL * gyz_physL * gxz_physL
                                + gxz_physL * gxy_physL * gyz_physL
                                - gxz_physL * gyy_physL * gxz_physL
                                - gxy_physL * gxy_physL * gzz_physL
                                - gxx_physL * gyz_physL * gyz_physL);

    const CCTK_REAL psiL = exp( (1.0/12.0) * log(gijdet) );

    const CCTK_REAL Psim4 = 1.0/(psiL*psiL*psiL*psiL);
    CCTK_REAL gtxxL = gxx_physL*Psim4;
    CCTK_REAL gtxyL = gxy_physL*Psim4;
    CCTK_REAL gtxzL = gxz_physL*Psim4;
    CCTK_REAL gtyyL = gyy_physL*Psim4;
    CCTK_REAL gtyzL = gyz_physL*Psim4;
    CCTK_REAL gtzzL = gzz_physL*Psim4;

    /*********************************
     * Apply det gtij = 1 constraint *
     *********************************/
    const CCTK_REAL gtijdet = gtxxL * gtyyL * gtzzL
                            + gtxyL * gtyzL * gtxzL
                            + gtxzL * gtxyL * gtyzL
                            - gtxzL * gtyyL * gtxzL
                            - gtxyL * gtxyL * gtzzL
                            - gtxxL * gtyzL * gtyzL;

    const CCTK_REAL gtijdet_Fm1o3 = fabs(1.0/cbrt(gtijdet));

    gtxxL = gtxxL * gtijdet_Fm1o3;
    gtxyL = gtxyL * gtijdet_Fm1o3;
    gtxzL = gtxzL * gtijdet_Fm1o3;
    gtyyL = gtyyL * gtijdet_Fm1o3;
    gtyzL = gtyzL * gtijdet_Fm1o3;
    gtzzL = gtzzL * gtijdet_Fm1o3;

    if(gtijdet<0.0)
      CCTK_VWARN(CCTK_WARN_ALERT,
                 "WARNING: det[3-metric]<0.0 at coordinate %e %e %e | cctk_lsh: %d %d %d. "
                 "Hopefully this is occurring in gz's! gtij_phys = %.2e %.2e %.2e %.2e %.2e %.2e "
                 "gtij_new = %.2e %.2e %.2e %.2e %.2e %.2e | gijdet = %.2e | gtijdet = %.2e",
                 p.x, p.y, p.z, cctkGH->cctk_lsh[0], cctkGH->cctk_lsh[1], cctkGH->cctk_lsh[2],
                 gxx_physL, gxy_physL, gxz_physL, gyy_physL, gyz_physL, gzz_physL,
                 gtxxL, gtxyL, gtxzL, gtyyL, gtyzL, gtzzL, -gijdet, gtijdet);

    /*******************************************
     * Set the ADM gridfunctions to new values *
     *******************************************/
    const CCTK_REAL Psi4 = psiL*psiL*psiL*psiL;
    ccc_gxx(index) = gtxxL*Psi4;
    ccc_gxy(index) = gtxyL*Psi4;
    ccc_gxz(index) = gtxzL*Psi4;
    ccc_gyy(index) = gtyyL*Psi4;
    ccc_gyz(index) = gtyzL*Psi4;
    ccc_gzz(index) = gtzzL*Psi4;
  }); // ccc loop everywhere
}
