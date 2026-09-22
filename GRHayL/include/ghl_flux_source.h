#ifndef GHL_FLUX_SOURCE_H_
#define GHL_FLUX_SOURCE_H_

#include "ghl.h"
#include <float.h>

/* GRHayL primitive magnetic fields already include the 1/sqrt(4 pi) rescaling. */
static const double SQRT_4_PI = 1;

#ifdef __cplusplus
extern "C" {
#endif

/** @addtogroup Flux_Source
 *
 * All routines use `ghl_compute_h_and_cs2`. Checked entry points return EOS
 * and wave-speed errors without changing outputs. Legacy entry points retain
 * their original `void` signatures and pass checked errors to
 * ghl_abort_if_error().
 *
 *  @{
 */

/** Compute GRMHD source terms.
 *
 * The EOS callback may update `prims`. On failure, `cons` is unchanged; the
 * checked counterpart returns the callback status and this legacy entry point
 * aborts. On success, the routine writes `tau` and all three `SD` components;
 * other conservative fields are not outputs.
 */
void ghl_calculate_source_terms(
      const ghl_eos_parameters *restrict eos,
      ghl_primitive_quantities *restrict prims,
      const ghl_metric_quantities *restrict metric,
      const ghl_metric_quantities *restrict metric_derivs_x,
      const ghl_metric_quantities *restrict metric_derivs_y,
      const ghl_metric_quantities *restrict metric_derivs_z,
      const ghl_extrinsic_curvature *restrict curv,
      ghl_conservative_quantities *restrict cons);

/** Compute direction-0 characteristic-speed magnitudes.
 *
 * The EOS callback may update either primitive state. On failure, both speed
 * outputs are unchanged; the checked counterpart returns the callback status
 * and this legacy entry point aborts.
 */
void ghl_calculate_characteristic_speed_dirn0(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      double *cmin_dirn0,
      double *cmax_dirn0);

/** Direction-1 counterpart of ghl_calculate_characteristic_speed_dirn0(). */
void ghl_calculate_characteristic_speed_dirn1(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      double *cmin_dirn1,
      double *cmax_dirn1);

/** Direction-2 counterpart of ghl_calculate_characteristic_speed_dirn0(). */
void ghl_calculate_characteristic_speed_dirn2(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      double *cmin_dirn2,
      double *cmax_dirn2);

/** Compute direction-0 hybrid-EOS HLLE fluxes.
 *
 * `cmin_dirn0` and `cmax_dirn0` are nonnegative wave-speed magnitudes.
 * Negative algebraic residue within `DBL_EPSILON` times the larger of one and
 * both magnitudes is floored at zero. Larger negative values are rejected.
 * The floored sum must be finite and at least `1/DBL_MAX`, and the floored
 * product must not overflow; the overflow test is exact only to within one
 * rounding. The EOS callback may update either primitive state. On
 * failure, `cons` is unchanged; the checked counterpart returns the error and
 * this legacy entry point aborts. On success, the
 * routine writes `rho`, `tau`, and all three `SD` components. Magnetic fields
 * in the primitive states are already rescaled by \f$1/\sqrt{4\pi}\f$.
 */
void ghl_calculate_HLLE_fluxes_dirn0_hybrid(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn0,
      const double cmax_dirn0,
      ghl_conservative_quantities *restrict cons);

/** Direction-1 counterpart of ghl_calculate_HLLE_fluxes_dirn0_hybrid(). */
void ghl_calculate_HLLE_fluxes_dirn1_hybrid(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn1,
      const double cmax_dirn1,
      ghl_conservative_quantities *restrict cons);

/** Direction-2 counterpart of ghl_calculate_HLLE_fluxes_dirn0_hybrid(). */
void ghl_calculate_HLLE_fluxes_dirn2_hybrid(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn2,
      const double cmax_dirn2,
      ghl_conservative_quantities *restrict cons);

/** Hybrid direction-0 HLLE flux including the `entropy` output. */
void ghl_calculate_HLLE_fluxes_dirn0_hybrid_entropy(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn0,
      const double cmax_dirn0,
      ghl_conservative_quantities *restrict cons);

/** Hybrid direction-1 HLLE flux including the `entropy` output. */
void ghl_calculate_HLLE_fluxes_dirn1_hybrid_entropy(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn1,
      const double cmax_dirn1,
      ghl_conservative_quantities *restrict cons);

/** Hybrid direction-2 HLLE flux including the `entropy` output. */
void ghl_calculate_HLLE_fluxes_dirn2_hybrid_entropy(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn2,
      const double cmax_dirn2,
      ghl_conservative_quantities *restrict cons);

/** Tabulated-EOS direction-0 HLLE flux; also writes `Y_e`.
 *
 * Wave-speed, mutation, error, and magnetic-rescaling contracts match
 * ghl_calculate_HLLE_fluxes_dirn0_hybrid().
 */
void ghl_calculate_HLLE_fluxes_dirn0_tabulated(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn0,
      const double cmax_dirn0,
      ghl_conservative_quantities *restrict cons);

/** Direction-1 counterpart of ghl_calculate_HLLE_fluxes_dirn0_tabulated(). */
void ghl_calculate_HLLE_fluxes_dirn1_tabulated(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn1,
      const double cmax_dirn1,
      ghl_conservative_quantities *restrict cons);

/** Direction-2 counterpart of ghl_calculate_HLLE_fluxes_dirn0_tabulated(). */
void ghl_calculate_HLLE_fluxes_dirn2_tabulated(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn2,
      const double cmax_dirn2,
      ghl_conservative_quantities *restrict cons);

/** Tabulated direction-0 HLLE flux including `Y_e` and `entropy` outputs. */
void ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn0,
      const double cmax_dirn0,
      ghl_conservative_quantities *restrict cons);

/** Tabulated direction-1 HLLE flux including `Y_e` and `entropy` outputs. */
void ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn1,
      const double cmax_dirn1,
      ghl_conservative_quantities *restrict cons);

/** Tabulated direction-2 HLLE flux including `Y_e` and `entropy` outputs. */
void ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn2,
      const double cmax_dirn2,
      ghl_conservative_quantities *restrict cons);

/** Checked source-term entry point. */
ghl_error_codes_t ghl_calculate_source_terms_checked(
      const ghl_eos_parameters *restrict eos,
      ghl_primitive_quantities *restrict prims,
      const ghl_metric_quantities *restrict metric,
      const ghl_metric_quantities *restrict metric_derivs_x,
      const ghl_metric_quantities *restrict metric_derivs_y,
      const ghl_metric_quantities *restrict metric_derivs_z,
      const ghl_extrinsic_curvature *restrict curv,
      ghl_conservative_quantities *restrict cons);

/** Checked characteristic-speed entry points. */
ghl_error_codes_t ghl_calculate_characteristic_speed_dirn0_checked(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      double *cmin_dirn0,
      double *cmax_dirn0);

ghl_error_codes_t ghl_calculate_characteristic_speed_dirn1_checked(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      double *cmin_dirn1,
      double *cmax_dirn1);

ghl_error_codes_t ghl_calculate_characteristic_speed_dirn2_checked(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      double *cmin_dirn2,
      double *cmax_dirn2);

/** Checked HLLE entry points. */
ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn0_hybrid_checked(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn0,
      const double cmax_dirn0,
      ghl_conservative_quantities *restrict cons);

ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn1_hybrid_checked(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn1,
      const double cmax_dirn1,
      ghl_conservative_quantities *restrict cons);

ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn2_hybrid_checked(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn2,
      const double cmax_dirn2,
      ghl_conservative_quantities *restrict cons);

ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn0_hybrid_entropy_checked(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn0,
      const double cmax_dirn0,
      ghl_conservative_quantities *restrict cons);

ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn1_hybrid_entropy_checked(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn1,
      const double cmax_dirn1,
      ghl_conservative_quantities *restrict cons);

ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn2_hybrid_entropy_checked(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn2,
      const double cmax_dirn2,
      ghl_conservative_quantities *restrict cons);

ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn0_tabulated_checked(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn0,
      const double cmax_dirn0,
      ghl_conservative_quantities *restrict cons);

ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn1_tabulated_checked(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn1,
      const double cmax_dirn1,
      ghl_conservative_quantities *restrict cons);

ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn2_tabulated_checked(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn2,
      const double cmax_dirn2,
      ghl_conservative_quantities *restrict cons);

ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy_checked(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn0,
      const double cmax_dirn0,
      ghl_conservative_quantities *restrict cons);

ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy_checked(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn1,
      const double cmax_dirn1,
      ghl_conservative_quantities *restrict cons);

ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy_checked(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn2,
      const double cmax_dirn2,
      ghl_conservative_quantities *restrict cons);

/** @} */

#ifdef __cplusplus
}
#endif

#endif // GHL_FLUX_SOURCE_H_
