#include "GRHayLHDX.h"

extern "C" void convert_GRHayLHDX_to_HydroBase(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_convert_GRHayLHDX_to_HydroBase;
  DECLARE_CCTK_PARAMETERS;

  // Generally, we only need the HydroBase variables for diagnostic purposes, so we run the below loop only at iterations in which diagnostics are run.
//  if(Convert_to_HydroBase_every==0 || cctk_iteration%Convert_to_HydroBase_every!=0) return;

  constexpr std::array<int, Loop::dim> indextype = {1, 1, 1};
  const Loop::GF3D2layout layout(cctkGH, indextype);

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);

    /* Note that we currently do not set Abar, Y_e, temperature, entropy, Avec[3], Aphi, Avec_stag[3], Aphi_stag */
    rho(index)   = rho_b(index);
    press(index) = pressure(index);

    // IllinoisGRMHD defines v^i = u^i/u^0.

    // Meanwhile, the ET/HydroBase formalism, called the Valencia
    // formalism, splits the 4 velocity into a purely spatial part
    // and a part that is normal to the spatial hypersurface:
    // u^a = G (n^a + U^a), (Eq. 14 of arXiv:1304.5544; G=W, U^a=v^a)
    // where n^a is the unit normal vector to the spatial hypersurface,
    // n_a = {-\alpha,0,0,0}, and U^a is the purely spatial part, which
    // is defined in HydroBase as the vel[] vector gridfunction.
    // Then u^a n_a = - \alpha u^0 = G n^a n_a = -G, and
    // of course \alpha u^0 = 1/sqrt(1+γ^ij u_i u_j) = \Gamma,
    // the standard Lorentz factor.

    // Note that n^i = - \beta^i / \alpha, so
    // u^a = \Gamma (n^a + U^a)
    // -> u^i = \Gamma ( U^i - \beta^i / \alpha )
    // which implies
    // v^i = u^i/u^0
    //     = \Gamma/u^0 ( U^i - \beta^i / \alpha ) <- \Gamma = \alpha u^0
    //     = \alpha ( U^i - \beta^i / \alpha )
    //     = \alpha U^i - \beta^i
    const CCTK_REAL lapseL = ccc_lapse(index);
    const CCTK_REAL lapseL_inv = 1.0/lapseL;

    velx(index) = lapseL_inv*(vx(index) + ccc_betax(index));
    vely(index) = lapseL_inv*(vy(index) + ccc_betay(index));
    velz(index) = lapseL_inv*(vz(index) + ccc_betaz(index));

    Bvecx(index) = 0.0;
    Bvecy(index) = 0.0;
    Bvecz(index) = 0.0;
  }); // ccc loop everywhere
}
