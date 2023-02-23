#ifndef FLUX_SOURCE_H_
#define FLUX_SOURCE_H_

#include "GRHayL.h"


static const double TINYDOUBLE = 1e-100;

static const double SQRT_4_PI = 3.544907701811032054596334966682290365L;


typedef struct extrinsic_curvature {
  double Kxx, Kxy, Kxz;
  double Kyy, Kyz, Kzz;
} extrinsic_curvature;

typedef struct metric_derivatives {
  double lapse[3];
  double betax[3];
  double betay[3];
  double betaz[3];
  double adm_gxx[3];
  double adm_gxy[3];
  double adm_gxz[3];
  double adm_gyy[3];
  double adm_gyz[3];
  double adm_gzz[3];
} metric_derivatives;


void calculate_characteristic_speed_dirn0(const primitive_quantities *restrict prims_r, 
                                          const primitive_quantities *restrict prims_l, 
                                          const eos_parameters *restrict eos, 
                                          const metric_quantities *restrict metric_face, 
                                          double *cmin_dirn0, 
                                          double *cmax_dirn0);

void calculate_characteristic_speed_dirn1(const primitive_quantities *restrict prims_r, 
                                          const primitive_quantities *restrict prims_l, 
                                          const eos_parameters *restrict eos, 
                                          const metric_quantities *restrict metric_face, 
                                          double *cmin_dirn1, 
                                          double *cmax_dirn1);

void calculate_characteristic_speed_dirn2(const primitive_quantities *restrict prims_r, 
                                          const primitive_quantities *restrict prims_l, 
                                          const eos_parameters *restrict eos, 
                                          const metric_quantities *restrict metric_face, 
                                          double *cmin_dirn2, 
                                          double *cmax_dirn2);

void calculate_HLLE_fluxes_dirn0(const primitive_quantities *restrict prims_r, 
                                 const primitive_quantities *restrict prims_l, 
                                 const eos_parameters *restrict eos, 
                                 const metric_quantities *restrict metric_face, 
                                 conservative_quantities *restrict cons);

void calculate_HLLE_fluxes_dirn1(const primitive_quantities *restrict prims_r, 
                                 const primitive_quantities *restrict prims_l, 
                                 const eos_parameters *restrict eos, 
                                 const metric_quantities *restrict metric_face, 
                                 conservative_quantities *restrict cons);

void calculate_HLLE_fluxes_dirn2(const primitive_quantities *restrict prims_r, 
                                 const primitive_quantities *restrict prims_l, 
                                 const eos_parameters *restrict eos, 
                                 const metric_quantities *restrict metric_face, 
                                 conservative_quantities *restrict cons);

void calculate_all_source_terms(const primitive_quantities *restrict prims, 
                                const eos_parameters *restrict eos, 
                                const metric_quantities *restrict metric,
                                const extrinsic_curvature *restrict curv, 
                                const metric_derivatives *restrict metric_derivs, 
                                conservative_quantities *restrict cons);

#endif // FLUX_SOURCE_H_