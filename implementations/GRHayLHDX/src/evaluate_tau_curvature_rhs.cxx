#include "GRHayLHDX.h"

void GRHayLHDX_evaluate_tau_curvature_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_evaluate_tau_curvature_rhs;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL poison = 0.0/0.0;

  constexpr std::array<int, Loop::dim> indextype = {1, 1, 1};
  const Loop::GF3D2layout layout(cctkGH, indextype);

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);

    rho_star_rhs(index) = 0.0;
    tau_rhs(index)      = 0.0;
    Stildex_rhs(index)  = 0.0;
    Stildey_rhs(index)  = 0.0;
    Stildez_rhs(index)  = 0.0;

    metric_quantities ADM_metric;
    ghl_initialize_metric(ccc_lapse(index),
                      ccc_betax(index), ccc_betay(index), ccc_betaz(index),
                      ccc_gxx(index), ccc_gxy(index), ccc_gxz(index),
                      ccc_gyy(index), ccc_gyz(index), ccc_gzz(index),
                      &ADM_metric);

    extrinsic_curvature curv;
    ghl_initialize_extrinsic_curvature(
          ccc_kxx(index), ccc_kxy(index), ccc_kxz(index),
          ccc_kyy(index), ccc_kyz(index), ccc_kzz(index),
          &curv);

    primitive_quantities prims;
    ghl_initialize_primitives(
          rho_b(index), pressure(index), eps(index),
          vx(index), vy(index), vz(index),
          0.0, 0.0, 0.0,
          poison, poison, poison, // entropy, Y_e, temp
          &prims);

    int speed_limited = 0;
    ghl_limit_v_and_compute_u0(
          ghl_eos, &ADM_metric, &prims, &speed_limited);

    conservative_quantities cons_source;
    cons_source.tau = 0;
    ghl_calculate_tau_tilde_source_term_extrinsic_curv(&prims, ghl_eos, &ADM_metric, &curv, &cons_source);
    tau_rhs(index) += cons_source.tau;
  }); // ccc loop interior
}
