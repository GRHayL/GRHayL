#ifndef FLUX_SOURCE_H_
#define FLUX_SOURCE_H_

#include "grhayl.h"


static const double TINYDOUBLE = 1e-100;

static const double SQRT_4_PI = 3.544907701811032054596334966682290365L;

#ifdef __cplusplus
extern "C" {
#endif

void ghl_calculate_characteristic_speed_dirn0(const primitive_quantities *restrict prims_r, 
      const primitive_quantities *restrict prims_l, 
      const eos_parameters *restrict eos, 
      const metric_quantities *restrict metric_face, 
      double *cmin_dirn0, 
      double *cmax_dirn0);

void ghl_calculate_characteristic_speed_dirn1(const primitive_quantities *restrict prims_r, 
      const primitive_quantities *restrict prims_l, 
      const eos_parameters *restrict eos, 
      const metric_quantities *restrict metric_face, 
      double *cmin_dirn1, 
      double *cmax_dirn1);

void ghl_calculate_characteristic_speed_dirn2(const primitive_quantities *restrict prims_r, 
      const primitive_quantities *restrict prims_l, 
      const eos_parameters *restrict eos, 
      const metric_quantities *restrict metric_face, 
      double *cmin_dirn2, 
      double *cmax_dirn2);

void ghl_calculate_HLLE_fluxes_dirn0(const primitive_quantities *restrict prims_r, 
      const primitive_quantities *restrict prims_l, 
      const eos_parameters *restrict eos, 
      const metric_quantities *restrict metric_face,
      const double cmin_dirn0,
      const double cmax_dirn0,
      conservative_quantities *restrict cons);

void ghl_calculate_HLLE_fluxes_dirn1(const primitive_quantities *restrict prims_r, 
      const primitive_quantities *restrict prims_l, 
      const eos_parameters *restrict eos, 
      const metric_quantities *restrict metric_face,
      const double cmin_dirn01,
      const double cmax_dirn01,
      conservative_quantities *restrict cons);

void ghl_calculate_HLLE_fluxes_dirn2(const primitive_quantities *restrict prims_r, 
      const primitive_quantities *restrict prims_l, 
      const eos_parameters *restrict eos, 
      const metric_quantities *restrict metric_face,
      const double cmin_dirn2,
      const double cmax_dirn2, 
      conservative_quantities *restrict cons);

void ghl_calculate_tau_tilde_source_term_extrinsic_curv(
      const primitive_quantities *restrict prims,
      struct eos_parameters const *restrict eos,
      const metric_quantities *restrict metric,
      const extrinsic_curvature *restrict curv,
      conservative_quantities *restrict cons);

void ghl_calculate_source_terms_dirn0(
      const primitive_quantities *restrict prims,
      struct eos_parameters const *restrict eos,
      const metric_quantities *restrict metric,
      const metric_quantities *restrict metric_derivs,
      conservative_quantities *restrict cons);

void ghl_calculate_source_terms_dirn1(
      const primitive_quantities *restrict prims,
      struct eos_parameters const *restrict eos,
      const metric_quantities *restrict metric,
      const metric_quantities *restrict metric_derivs,
      conservative_quantities *restrict cons);

void ghl_calculate_source_terms_dirn2(
      const primitive_quantities *restrict prims,
      struct eos_parameters const *restrict eos,
      const metric_quantities *restrict metric,
      const metric_quantities *restrict metric_derivs,
      conservative_quantities *restrict cons);

#ifdef __cplusplus
}
#endif

#endif // FLUX_SOURCE_H_
