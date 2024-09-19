#ifndef GHL_FLUX_SOURCE_H_
#define GHL_FLUX_SOURCE_H_

#include "ghl.h"


static const double TINYDOUBLE = 1e-100;

static const double SQRT_4_PI = 1; //3.544907701811032054596334966682290365L;

#ifdef __cplusplus
extern "C" {
#endif

void ghl_calculate_source_terms(
      const ghl_eos_parameters *restrict eos,
      ghl_primitive_quantities *restrict prims,
      const ghl_metric_quantities *restrict metric,
      const ghl_metric_quantities *restrict metric_derivs_x,
      const ghl_metric_quantities *restrict metric_derivs_y,
      const ghl_metric_quantities *restrict metric_derivs_z,
      const ghl_extrinsic_curvature *restrict curv,
      ghl_conservative_quantities *restrict cons);

void ghl_calculate_characteristic_speed(
      const int direction,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      double *cmin_dirn0,
      double *cmax_dirn0);

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
