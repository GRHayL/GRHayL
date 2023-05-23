#include "GRHayLHDX.h"

extern "C" void GRHayLHDX_compute_ccc_centered_extrinsic_curvature(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_compute_ccc_centered_extrinsic_curvature;
  DECLARE_CCTK_PARAMETERS;

  constexpr std::array<int, Loop::dim> indextype = {1, 1, 1};
  const Loop::GF3D2layout layout(cctkGH, indextype);

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);

    CCTK_REAL kxx_avg   = 0.0;
    CCTK_REAL kxy_avg   = 0.0;
    CCTK_REAL kxz_avg   = 0.0;
    CCTK_REAL kyy_avg   = 0.0;
    CCTK_REAL kyz_avg   = 0.0;
    CCTK_REAL kzz_avg   = 0.0;

    for (int k=0; k<2; k++) { 
      for (int j=0; j<2; j++) { 
        for (int i=0; i<2; i++) {
          kxx_avg += kxx(p.I + i*p.DI[0] + j*p.DI[1] + k*p.DI[2]);
          kxy_avg += kxy(p.I + i*p.DI[0] + j*p.DI[1] + k*p.DI[2]);
          kxz_avg += kxz(p.I + i*p.DI[0] + j*p.DI[1] + k*p.DI[2]);
          kyy_avg += kyy(p.I + i*p.DI[0] + j*p.DI[1] + k*p.DI[2]);
          kyz_avg += kyz(p.I + i*p.DI[0] + j*p.DI[1] + k*p.DI[2]);
          kzz_avg += kzz(p.I + i*p.DI[0] + j*p.DI[1] + k*p.DI[2]);
        }
      }
    }
    ccc_kxx(index) = kxx_avg/8.0;
    ccc_kxy(index) = kxy_avg/8.0;
    ccc_kxz(index) = kxz_avg/8.0;
    ccc_kyy(index) = kyy_avg/8.0;
    ccc_kyz(index) = kyz_avg/8.0;
    ccc_kzz(index) = kzz_avg/8.0;
  }); // ccc loop everywhere
}
