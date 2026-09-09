#ifndef GHL_FLUX_SOURCE_H_
#define GHL_FLUX_SOURCE_H_

#include "ghl.h"

static const double TINYDOUBLE = 1e-100;

static const double SQRT_4_PI = 1; //3.544907701811032054596334966682290365L;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Compute a component-wise symmetric Rusanov interface flux.
 * @ingroup Flux_Source
 *
 * The component arrays contain undensitized conserved quantities and their
 * corresponding physical fluxes on the left and right sides. The same
 * nonnegative speed is applied to every component. The output buffer must be
 * distinct from the input buffers. Inputs and candidates are validated before
 * publication, so an error leaves the output buffer unchanged.
 *
 * @param state_L Undensitized left state with @p component_count components.
 * @param state_R Undensitized right state with @p component_count components.
 * @param physical_flux_L Undensitized physical flux corresponding to
 *        @p state_L.
 * @param physical_flux_R Undensitized physical flux corresponding to
 *        @p state_R.
 * @param component_count Number of components in each input and output
 *        array; it must be positive.
 * @param speed Nonnegative interface speed applied componentwise.
 * @param flux Output numerical flux with @p component_count components. It is
 *        not densitized by this helper.
 * @return @c ghl_success on publication; otherwise @p flux is unchanged.
 */
ghl_error_codes_t ghl_calculate_Rusanov_flux(
      const double *restrict state_L,
      const double *restrict state_R,
      const double *restrict physical_flux_L,
      const double *restrict physical_flux_R,
      const int component_count,
      const double speed,
      double *restrict flux);

void ghl_calculate_source_terms(
      const ghl_eos_parameters *restrict eos,
      ghl_primitive_quantities *restrict prims,
      const ghl_metric_quantities *restrict metric,
      const ghl_metric_quantities *restrict metric_derivs_x,
      const ghl_metric_quantities *restrict metric_derivs_y,
      const ghl_metric_quantities *restrict metric_derivs_z,
      const ghl_extrinsic_curvature *restrict curv,
      ghl_conservative_quantities *restrict cons);

void ghl_calculate_characteristic_speed_dirn0(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      double *cmin_dirn0,
      double *cmax_dirn0);

void ghl_calculate_characteristic_speed_dirn1(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      double *cmin_dirn1,
      double *cmax_dirn1);

void ghl_calculate_characteristic_speed_dirn2(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      double *cmin_dirn2,
      double *cmax_dirn2);

void ghl_calculate_HLLE_fluxes_dirn0_hybrid(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn0,
      const double cmax_dirn0,
      ghl_conservative_quantities *restrict cons);

void ghl_calculate_HLLE_fluxes_dirn1_hybrid(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn1,
      const double cmax_dirn1,
      ghl_conservative_quantities *restrict cons);

void ghl_calculate_HLLE_fluxes_dirn2_hybrid(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn2,
      const double cmax_dirn2,
      ghl_conservative_quantities *restrict cons);

void ghl_calculate_HLLE_fluxes_dirn0_hybrid_entropy(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn0,
      const double cmax_dirn0,
      ghl_conservative_quantities *restrict cons);

void ghl_calculate_HLLE_fluxes_dirn1_hybrid_entropy(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn1,
      const double cmax_dirn1,
      ghl_conservative_quantities *restrict cons);

void ghl_calculate_HLLE_fluxes_dirn2_hybrid_entropy(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn2,
      const double cmax_dirn2,
      ghl_conservative_quantities *restrict cons);

void ghl_calculate_HLLE_fluxes_dirn0_tabulated(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn0,
      const double cmax_dirn0,
      ghl_conservative_quantities *restrict cons);

void ghl_calculate_HLLE_fluxes_dirn1_tabulated(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn1,
      const double cmax_dirn1,
      ghl_conservative_quantities *restrict cons);

void ghl_calculate_HLLE_fluxes_dirn2_tabulated(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn2,
      const double cmax_dirn2,
      ghl_conservative_quantities *restrict cons);

void ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn0,
      const double cmax_dirn0,
      ghl_conservative_quantities *restrict cons);

void ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn1,
      const double cmax_dirn1,
      ghl_conservative_quantities *restrict cons);

void ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn2,
      const double cmax_dirn2,
      ghl_conservative_quantities *restrict cons);

#ifdef __cplusplus
}
#endif

#endif // GHL_FLUX_SOURCE_H_
