#include "GRHayLHDX.h"

extern "C" void GRHayLHDX_compute_Tmunu(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_compute_Tmunu;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL poison = 0.0/0.0;

  constexpr std::array<int, Loop::dim> vvvtype = {0, 0, 0};
  const Loop::GF3D2layout vvv_layout(cctkGH, vvvtype);

  constexpr std::array<int, Loop::dim> ccctype = {1, 1, 1};
  const Loop::GF3D2layout ccc_layout(cctkGH, ccctype);

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(ccc_layout, p.I);

    metric_quantities ADM_metric;
    ghl_initialize_metric(
          ccc_lapse(index),
          ccc_betax(index), ccc_betay(index), ccc_betaz(index),
          ccc_gxx(index), ccc_gxy(index), ccc_gxz(index),
          ccc_gyy(index), ccc_gyz(index), ccc_gzz(index),
          &ADM_metric);

    ADM_aux_quantities metric_aux;
    ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

    primitive_quantities prims;
    ghl_initialize_primitives(
          rho_b(index), pressure(index), eps(index),
          vx(index), vy(index), vz(index),
          0.0, 0.0, 0.0,
          poison, poison, poison, &prims);
    prims.u0 = u0(index);

    stress_energy Tmunu;
    ghl_compute_TDNmunu(
          &ADM_metric, &metric_aux, &prims, &Tmunu);

    ghl_return_stress_energy(
          &Tmunu, &ccc_Ttt(index), &ccc_Ttx(index),
          &ccc_Tty(index), &ccc_Ttz(index), &ccc_Txx(index),
          &ccc_Txy(index), &ccc_Txz(index), &ccc_Tyy(index),
          &ccc_Tyz(index), &ccc_Tzz(index));
  }); // ccc loop everywhere

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(vvv_layout, p.I);

    CCTK_REAL Ttt_avg   = 0.0;
    CCTK_REAL Ttx_avg   = 0.0;
    CCTK_REAL Tty_avg   = 0.0;
    CCTK_REAL Ttz_avg   = 0.0;
    CCTK_REAL Txx_avg   = 0.0;
    CCTK_REAL Txy_avg   = 0.0;
    CCTK_REAL Txz_avg   = 0.0;
    CCTK_REAL Tyy_avg   = 0.0;
    CCTK_REAL Tyz_avg   = 0.0;
    CCTK_REAL Tzz_avg   = 0.0;

    for (int k=0; k<2; k++) { 
      for (int j=0; j<2; j++) { 
        for (int i=0; i<2; i++) {
          Ttt_avg += ccc_Ttt(p.I - i*p.DI[0] - j*p.DI[1] - k*p.DI[2]);
          Ttx_avg += ccc_Ttx(p.I - i*p.DI[0] - j*p.DI[1] - k*p.DI[2]);
          Tty_avg += ccc_Tty(p.I - i*p.DI[0] - j*p.DI[1] - k*p.DI[2]);
          Ttz_avg += ccc_Ttz(p.I - i*p.DI[0] - j*p.DI[1] - k*p.DI[2]);
          Txx_avg += ccc_Txx(p.I - i*p.DI[0] - j*p.DI[1] - k*p.DI[2]);
          Txy_avg += ccc_Txy(p.I - i*p.DI[0] - j*p.DI[1] - k*p.DI[2]);
          Txz_avg += ccc_Txz(p.I - i*p.DI[0] - j*p.DI[1] - k*p.DI[2]);
          Tyy_avg += ccc_Tyy(p.I - i*p.DI[0] - j*p.DI[1] - k*p.DI[2]);
          Tyz_avg += ccc_Tyz(p.I - i*p.DI[0] - j*p.DI[1] - k*p.DI[2]);
          Tzz_avg += ccc_Tzz(p.I - i*p.DI[0] - j*p.DI[1] - k*p.DI[2]);
        }
      }
    }
    eTtt(index) = Ttt_avg/8.0;
    eTtx(index) = Ttx_avg/8.0;
    eTty(index) = Tty_avg/8.0;
    eTtz(index) = Ttz_avg/8.0;
    eTxx(index) = Txx_avg/8.0;
    eTxy(index) = Txy_avg/8.0;
    eTxz(index) = Txz_avg/8.0;
    eTyy(index) = Tyy_avg/8.0;
    eTyz(index) = Tyz_avg/8.0;
    eTzz(index) = Tzz_avg/8.0;
  }); // vvv loop interior
}
