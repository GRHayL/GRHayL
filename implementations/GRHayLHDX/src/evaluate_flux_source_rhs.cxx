#include "GRHayLHDX.h"

template <int flux_dir>
         void GRHayLHDX_evaluate_flux_source_rhs_dir(
               cGH *restrict cctkGH,
               Loop::GF3D2<const CCTK_REAL> v_flux_dir,
               Loop::GF3D2<CCTK_REAL> rho_flux,
               Loop::GF3D2<CCTK_REAL> tau_flux,
               Loop::GF3D2<CCTK_REAL> Sx_flux,
               Loop::GF3D2<CCTK_REAL> Sy_flux,
               Loop::GF3D2<CCTK_REAL> Sz_flux) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_evaluate_flux_source_rhs;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL poison = 0.0/0.0;
  constexpr std::array<int, Loop::dim> ccctype = {1, 1, 1};
  const Loop::GF3D2layout ccc_layout(cctkGH, ccctype);

  constexpr std::array<int, Loop::dim> facetype = {flux_dir!=0, flux_dir!=1, flux_dir!=2};
  const Loop::GF3D2layout face_layout(cctkGH, facetype);

  // These nested condtional ternary operators let us tell the compiler that
  // these pointers can be set at compiler time while making the templated functions
  // instead of using a switch statement at runtime.
  constexpr void (*calculate_characteristic_speed)(const primitive_quantities *restrict prims_r,
                                         const primitive_quantities *restrict prims_l,
                                         struct eos_parameters const *restrict eos,
                                         const metric_quantities *restrict ADM_metric_face,
                                         CCTK_REAL *cmin, CCTK_REAL *cmax)
    = flux_dir==0 ? &ghl_calculate_characteristic_speed_dirn0 :
      flux_dir==1 ? &ghl_calculate_characteristic_speed_dirn1 :
                    &ghl_calculate_characteristic_speed_dirn2 ;

  constexpr void (*calculate_HLLE_fluxes)(const primitive_quantities *restrict prims_r,
                                const primitive_quantities *restrict prims_l,
                                const eos_parameters *restrict eos,
                                const metric_quantities *restrict ADM_metric_face,
                                const CCTK_REAL cmin,
                                const CCTK_REAL cmax,
                                conservative_quantities *restrict cons_fluxes)
    = flux_dir==0 ? &ghl_calculate_HLLE_fluxes_dirn0 :
      flux_dir==1 ? &ghl_calculate_HLLE_fluxes_dirn1 :
                    &ghl_calculate_HLLE_fluxes_dirn2 ;

  constexpr void (*calculate_source_terms)(const primitive_quantities *restrict prims,
                                 const eos_parameters *restrict eos,
                                 const metric_quantities *restrict ADM_metric,
                                 const metric_quantities *restrict metric_derivs,
                                 conservative_quantities *restrict cons_sources)
    = flux_dir==0 ? &ghl_calculate_source_terms_dirn0 :
      flux_dir==1 ? &ghl_calculate_source_terms_dirn1 :
                    &ghl_calculate_source_terms_dirn2 ;

  grid.loop_int_device<flux_dir!=0, flux_dir!=1, flux_dir!=2>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index indm2(ccc_layout, p.I + 2*p.DI[flux_dir]);
    const Loop::GF3D2index indm1(ccc_layout, p.I - p.DI[flux_dir]);
    const Loop::GF3D2index index(ccc_layout, p.I);
    const Loop::GF3D2index indp1(ccc_layout, p.I + p.DI[flux_dir]);

    const Loop::GF3D2index ind_flux(face_layout, p.I);
  
    CCTK_REAL rho_stencil[6], press_stencil[6], v_flux[6];
    CCTK_REAL rhor, rhol, pressr, pressl;
    CCTK_REAL vel_stencil[3][6], vel_r[3], vel_l[3];
  
    for(int ind=0; ind<6; ind++) {
      // Stencil from -3 to +2 reconstructs to e.g. i-1/2
      const Loop::GF3D2index stencil(ccc_layout, p.I + (ind-3)*p.DI[flux_dir]);
      v_flux[ind] = v_flux_dir(stencil); // Could be smaller; doesn't use full stencil
      rho_stencil[ind] = rho_b(stencil);
      press_stencil[ind] = pressure(stencil);
      vel_stencil[0][ind] = vx(stencil);
      vel_stencil[1][ind] = vy(stencil);
      vel_stencil[2][ind] = vz(stencil);
    }
  
    // Compute Gamma_eff
    CCTK_REAL K, Gamma;
    ghl_hybrid_get_K_and_Gamma(ghl_eos, rho_b(index), &K, &Gamma);
    const CCTK_REAL P_cold = K*pow(rho_b(index), Gamma);
    const CCTK_REAL Gamma_eff = ghl_eos->Gamma_th + (Gamma - ghl_eos->Gamma_th)*P_cold/pressure(index);
  
    ghl_simple_ppm(
          rho_stencil, press_stencil, vel_stencil,
          3, v_flux, Gamma_eff,
          &rhor, &rhol, &pressr, &pressl, vel_r, vel_l);
  
    metric_quantities ADM_metric_face;
    GRHayLHDX_interpolate_metric_to_face(
          indm2, indm1, index, indp1,
          ccc_lapse, ccc_betax, ccc_betay, ccc_betaz,
          ccc_gxx, ccc_gxy, ccc_gxz,
          ccc_gyy, ccc_gyz, ccc_gzz,
          &ADM_metric_face);
  
    primitive_quantities prims_r, prims_l;
    ghl_initialize_primitives(
          rhor, pressr, poison,
          vel_r[0], vel_r[1], vel_r[2],
          0.0, 0.0, 0.0,
          poison, poison, poison, // entropy, Y_e, temp
          &prims_r);
  
    ghl_initialize_primitives(
          rhol, pressl, poison,
          vel_l[0], vel_l[1], vel_l[2],
          0.0, 0.0, 0.0,
          poison, poison, poison, // entropy, Y_e, temp
          &prims_l);
  
    int speed_limited = 0;
    ghl_limit_v_and_compute_u0(
          ghl_eos, &ADM_metric_face, &prims_r, &speed_limited);
    ghl_limit_v_and_compute_u0(
          ghl_eos, &ADM_metric_face, &prims_l, &speed_limited);
  
    CCTK_REAL cmin, cmax;
    conservative_quantities cons_fluxes;
    calculate_characteristic_speed(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, &cmin, &cmax);
    calculate_HLLE_fluxes(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, cmin, cmax, &cons_fluxes);
  
    rho_flux(ind_flux) = cons_fluxes.rho;
    tau_flux(ind_flux) = cons_fluxes.tau;
    Sx_flux(ind_flux)  = cons_fluxes.SD[0];
    Sy_flux(ind_flux)  = cons_fluxes.SD[1];
    Sz_flux(ind_flux)  = cons_fluxes.SD[2];
  }); // staggered loop interior (e.g. flux_dir=0 gives vcc)

  const CCTK_REAL dxi = 1.0/CCTK_DELTA_SPACE(flux_dir);

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index indm2(ccc_layout, p.I - 2*p.DI[flux_dir]);
    const Loop::GF3D2index indm1(ccc_layout, p.I - p.DI[flux_dir]);
    const Loop::GF3D2index index(ccc_layout, p.I);
    const Loop::GF3D2index indp1(ccc_layout, p.I + p.DI[flux_dir]);
    const Loop::GF3D2index indp2(ccc_layout, p.I + 2*p.DI[flux_dir]);

    const Loop::GF3D2index ind_flux(face_layout, p.I);
    const Loop::GF3D2index ind_flp1(face_layout, p.I + p.DI[flux_dir]);

    rho_star_rhs(index) += dxi*(rho_flux(ind_flux) - rho_flux(ind_flp1));
    tau_rhs(index)      += dxi*(tau_flux(ind_flux) - tau_flux(ind_flp1));
    Stildex_rhs(index)  += dxi*(Sx_flux(ind_flux)  - Sx_flux(ind_flp1));
    Stildey_rhs(index)  += dxi*(Sy_flux(ind_flux)  - Sy_flux(ind_flp1));
    Stildez_rhs(index)  += dxi*(Sz_flux(ind_flux)  - Sz_flux(ind_flp1));

    metric_quantities ADM_metric;
    ghl_initialize_metric(ccc_lapse(index),
          ccc_betax(index), ccc_betay(index), ccc_betaz(index),
          ccc_gxx(index), ccc_gxy(index), ccc_gxz(index),
          ccc_gyy(index), ccc_gyz(index), ccc_gzz(index),
          &ADM_metric);

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

    metric_quantities ADM_metric_derivs;

    GRHayLHDX_compute_metric_derivs(
          dxi, indm2, indm1, indp1, indp2,
          ccc_lapse, ccc_betax, ccc_betay, ccc_betaz,
          ccc_gxx, ccc_gxy, ccc_gxz,
          ccc_gyy, ccc_gyz, ccc_gzz,
          &ADM_metric_derivs);

    conservative_quantities cons_source;
    cons_source.tau = 0.0;
    cons_source.SD[0] = 0.0;
    cons_source.SD[1] = 0.0;
    cons_source.SD[2] = 0.0;

    calculate_source_terms(&prims, ghl_eos, &ADM_metric, &ADM_metric_derivs, &cons_source);
    tau_rhs(index)     += cons_source.tau;
    Stildex_rhs(index) += cons_source.SD[0];
    Stildey_rhs(index) += cons_source.SD[1];
    Stildez_rhs(index) += cons_source.SD[2];
  }); // ccc loop interior
}

void GRHayLHDX_evaluate_flux_source_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_evaluate_flux_source_rhs;
  DECLARE_CCTK_PARAMETERS;

  /*
   *  Computation of \partial_i on RHS of \partial_t {rho_star,tau,Stilde{x,y,z}},
   *  via PPM reconstruction onto e.g. (i+1/2,j,k), so that
   *  \partial_x F = [ F(i+1/2,j,k) - F(i-1/2,j,k) ] / dx
   */
  GRHayLHDX_evaluate_flux_source_rhs_dir<0>(cctkGH, vx,
        rho_star_flux_x, tau_flux_x, Sx_flux_x, Sy_flux_x, Sz_flux_x);

  GRHayLHDX_evaluate_flux_source_rhs_dir<1>(cctkGH, vy,
        rho_star_flux_y, tau_flux_y, Sx_flux_y, Sy_flux_y, Sz_flux_y);

  GRHayLHDX_evaluate_flux_source_rhs_dir<2>(cctkGH, vz,
        rho_star_flux_z, tau_flux_z, Sx_flux_z, Sy_flux_z, Sz_flux_z);
}
