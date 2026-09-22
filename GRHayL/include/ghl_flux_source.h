#ifndef GHL_FLUX_SOURCE_H_
#define GHL_FLUX_SOURCE_H_

#include <float.h>
#include "ghl.h"

/* GRHayL primitive magnetic fields already include the 1/sqrt(4 pi) rescaling. */
static const double SQRT_4_PI = 1;

#ifdef __cplusplus
extern "C" {
#endif

/** @addtogroup Flux_Source
 *
 * Flux and source routines use `ghl_compute_h`; characteristic-speed
 * routines use `ghl_compute_h_and_cs2`. A caller that replaces EOS dispatch
 * callbacks must install both pointers consistently.
 *
 *  @{
 */

/** Compute GRMHD source terms.
 *
 * The EOS callback may update `prims`. On failure, `cons` is unchanged and the
 * callback status is returned. On success, the routine writes `tau` and all
 * three `SD` components; other conservative fields are not outputs.
 */
ghl_error_codes_t ghl_calculate_source_terms(
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
 * outputs are unchanged and the callback status is returned.
 */
ghl_error_codes_t ghl_calculate_characteristic_speed_dirn0(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      double *cmin_dirn0,
      double *cmax_dirn0);

/** Direction-1 counterpart of ghl_calculate_characteristic_speed_dirn0(). */
ghl_error_codes_t ghl_calculate_characteristic_speed_dirn1(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      double *cmin_dirn1,
      double *cmax_dirn1);

/** Direction-2 counterpart of ghl_calculate_characteristic_speed_dirn0(). */
ghl_error_codes_t ghl_calculate_characteristic_speed_dirn2(
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
 * both magnitudes is clamped to zero. Larger negative values are rejected.
 * The clamped sum and reciprocal must be finite and positive, and the clamped
 * product must be representable. The EOS callback may update
 * either primitive state. On failure, `cons` is unchanged; on success, the
 * routine writes `rho`, `tau`, and all three `SD` components. Magnetic fields
 * in the primitive states are already rescaled by \f$1/\sqrt{4\pi}\f$.
 */
ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn0_hybrid(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn0,
      const double cmax_dirn0,
      ghl_conservative_quantities *restrict cons);

/** Direction-1 counterpart of ghl_calculate_HLLE_fluxes_dirn0_hybrid(). */
ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn1_hybrid(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn1,
      const double cmax_dirn1,
      ghl_conservative_quantities *restrict cons);

/** Direction-2 counterpart of ghl_calculate_HLLE_fluxes_dirn0_hybrid(). */
ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn2_hybrid(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn2,
      const double cmax_dirn2,
      ghl_conservative_quantities *restrict cons);

/** Hybrid direction-0 HLLE flux including the `entropy` output. */
ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn0_hybrid_entropy(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn0,
      const double cmax_dirn0,
      ghl_conservative_quantities *restrict cons);

/** Hybrid direction-1 HLLE flux including the `entropy` output. */
ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn1_hybrid_entropy(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn1,
      const double cmax_dirn1,
      ghl_conservative_quantities *restrict cons);

/** Hybrid direction-2 HLLE flux including the `entropy` output. */
ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn2_hybrid_entropy(
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
ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn0_tabulated(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn0,
      const double cmax_dirn0,
      ghl_conservative_quantities *restrict cons);

/** Direction-1 counterpart of ghl_calculate_HLLE_fluxes_dirn0_tabulated(). */
ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn1_tabulated(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn1,
      const double cmax_dirn1,
      ghl_conservative_quantities *restrict cons);

/** Direction-2 counterpart of ghl_calculate_HLLE_fluxes_dirn0_tabulated(). */
ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn2_tabulated(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn2,
      const double cmax_dirn2,
      ghl_conservative_quantities *restrict cons);

/** Tabulated direction-0 HLLE flux including `Y_e` and `entropy` outputs. */
ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn0,
      const double cmax_dirn0,
      ghl_conservative_quantities *restrict cons);

/** Tabulated direction-1 HLLE flux including `Y_e` and `entropy` outputs. */
ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy(
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_face,
      const double cmin_dirn1,
      const double cmax_dirn1,
      ghl_conservative_quantities *restrict cons);

/** Tabulated direction-2 HLLE flux including `Y_e` and `entropy` outputs. */
ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy(
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
