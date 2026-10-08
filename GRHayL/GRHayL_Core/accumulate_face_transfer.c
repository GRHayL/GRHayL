#include "ghl.h"

/**
 * @ingroup GRHayL_Core
 * @brief Applies one integrated face transfer to the two cells sharing the face
 *
 * @details
 * First-order flux correction replaces the flux on the faces of flagged cells by a
 * first-order flux (\cite Lemaster_2009; \cite Fields_2025, Sec. 3.3). The update
 * stays conservative because each face contributes one transfer, which leaves one
 * cell and enters the other: the integrated form of the conservation law
 * \f[
 *   \partial_t u + \partial_i f^i(u) = 0
 * \f]
 * (Eq. 10 of \cite Fields_2025) for every conservative field. The caller computes the
 * transfer, the donor-cell flux integrated over the face and the stage interval,
 * and selects it for flagged faces.
 *
 * For a face oriented from cell L toward cell R, the transfer \f$ I_f \f$ leaves L
 * and enters R in every conservative field:
 *
 * \f[
 * \Delta Q_L \leftarrow \Delta Q_L - I_f,
 * \qquad
 * \Delta Q_R \leftarrow \Delta Q_R + I_f
 * \f]
 *
 * All three arguments are integrated densitized quantities (\f$ Q \f$, not
 * \f$ Q/V \f$), and fields the caller does not evolve must be zero in all of them.
 * Both updated accumulators are formed first and checked, so neither changes unless
 * every updated field of both is finite: a transfer that is not finite, an accumulator
 * that already is not, and a finite transfer that overflows an accumulator are all
 * rejected. Call this once per face; the two cells must be different objects. A cell
 * that is its own neighbor, as across a periodic boundary one cell wide, has no net
 * change, so the caller skips such a face.
 *
 * @param[in] transfer pointer to the integrated transfer \f$ I_f \f$
 *
 * @param[in,out] delta_L pointer to the accumulator of the cell the face leaves
 *
 * @param[in,out] delta_R pointer to the accumulator of the cell the face enters
 *
 * @returns ghl_success, or ghl_error_invalid_face_transfer if an updated accumulator
 *          would not be finite, which leaves both accumulators unchanged
 */
ghl_error_codes_t ghl_accumulate_face_transfer(
      const ghl_conservative_quantities *restrict transfer,
      ghl_conservative_quantities *restrict delta_L,
      ghl_conservative_quantities *restrict delta_R) {

  // Eq. 10 of Fields et al. (2025) in integrated form: one transfer, two cells.
  ghl_conservative_quantities new_L = *delta_L;
  ghl_conservative_quantities new_R = *delta_R;
  new_L.rho -= transfer->rho;
  new_R.rho += transfer->rho;
  new_L.tau -= transfer->tau;
  new_R.tau += transfer->tau;
  new_L.Y_e -= transfer->Y_e;
  new_R.Y_e += transfer->Y_e;
  new_L.entropy -= transfer->entropy;
  new_R.entropy += transfer->entropy;
  for(int i = 0; i < 3; i++) {
    new_L.SD[i] -= transfer->SD[i];
    new_R.SD[i] += transfer->SD[i];
  }

  // Commit only if both updates are finite. A finite transfer can still overflow an
  // accumulator, and a transfer or accumulator that is not finite always gives a result
  // that is not.
  const ghl_conservative_quantities *const updated[2] = { &new_L, &new_R };
  for(int k = 0; k < 2; k++) {
    if(!(isfinite(updated[k]->rho) && isfinite(updated[k]->tau)
         && isfinite(updated[k]->Y_e) && isfinite(updated[k]->entropy)
         && isfinite(updated[k]->SD[0]) && isfinite(updated[k]->SD[1])
         && isfinite(updated[k]->SD[2]))) {
      return ghl_error_invalid_face_transfer;
    }
  }
  *delta_L = new_L;
  *delta_R = new_R;
  return ghl_success;
}
