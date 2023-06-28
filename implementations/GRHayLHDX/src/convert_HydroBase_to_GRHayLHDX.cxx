/********************************
 * CONVERT ET ID TO IllinoisGRMHD
 *
 * Written in 2014 by Zachariah B. Etienne
 *
 * Sets metric & MHD variables needed
 * by IllinoisGRMHD, converting from
 * HydroBase and ADMBase.
 ********************************/

#include "GRHayLHDX.h"

extern "C" void convert_HydroBase_to_GRHayLHDX(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_convert_HydroBase_to_GRHayLHDX;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL poison = 0.0/0.0;

  constexpr std::array<int, Loop::dim> indextype = {1, 1, 1};
  const Loop::GF3D2layout layout(cctkGH, indextype);

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
   const Loop::GF3D2index index(layout, p.I);

   rho_b(index) = rho(index);
   pressure(index) = press(index);

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

   vx(index) = ccc_lapse(index)*velx(index) - ccc_betax(index);
   vy(index) = ccc_lapse(index)*vely(index) - ccc_betay(index);
   vz(index) = ccc_lapse(index)*velz(index) - ccc_betaz(index);
  }); // ccc loop everywhere

  // Neat feature for debugging: Add a roundoff-error perturbation
  //    to the initial data.
  // Set random_pert variable to ~1e-14 for a random 15th digit
  //    perturbation.
  if(random_pert > 1e-30) {
    srand(random_seed); // Use srand() as rand() is thread-safe.
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      const Loop::GF3D2index index(layout, p.I);
      const CCTK_REAL pert = (random_pert*(CCTK_REAL)rand() / RAND_MAX);
      const CCTK_REAL one_plus_pert=(1.0+pert);
      rho_b(index)*=one_plus_pert;
      vx(index)*=one_plus_pert;
      vy(index)*=one_plus_pert;
      vz(index)*=one_plus_pert;
    }); // ccc loop everywhere
  }

  // Finally, enforce limits on primitives & compute conservative variables.
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);

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
          poison, poison, poison,
          &prims);

    conservative_quantities cons;
    int speed_limited = 0;
    //This applies inequality fixes on the conservatives
    ghl_enforce_primitive_limits_and_compute_u0(
          ghl_params, ghl_eos, &ADM_metric,
          &prims, &speed_limited);
    //This computes the conservatives from the new primitives
    ghl_compute_conservs(
          &ADM_metric, &metric_aux, &prims, &cons);

    CCTK_REAL dummy1, dummy2, dummy3;
    CCTK_REAL dummy4, dummy5, dummy6;
    ghl_return_primitives(
          &prims,
          &rho_b(index), &pressure(index), &eps(index),
          &vx(index), &vy(index), &vz(index),
          &dummy1, &dummy2, &dummy3,
          &dummy4, &dummy5, &dummy6);

    ghl_return_conservatives(
          &cons,
          &rho_star(index), &tau(index),
          &Stildex(index), &Stildey(index), &Stildez(index),
          &dummy1, &dummy2);
  }); // ccc loop everywhere
}
