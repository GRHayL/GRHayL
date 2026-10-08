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
 * The transfer is checked first, so nothing changes if it is not finite. Call this
 * once per face; the two cells must be different objects. A cell that is its own
 * neighbor, as across a periodic boundary one cell wide, has no net change, so the
 * caller skips such a face.
 *
 * @param[in] transfer pointer to the integrated transfer \f$ I_f \f$
 *
 * @param[in,out] delta_L pointer to the accumulator of the cell the face leaves
 *
 * @param[in,out] delta_R pointer to the accumulator of the cell the face enters
 *
 * @returns ghl_success, or ghl_error_invalid_face_transfer if transfer is not finite
 */
ghl_error_codes_t ghl_accumulate_face_transfer(
      const ghl_conservative_quantities *restrict transfer,
      ghl_conservative_quantities *restrict delta_L,
      ghl_conservative_quantities *restrict delta_R) {

  if(!(isfinite(transfer->rho) && isfinite(transfer->tau) && isfinite(transfer->Y_e)
       && isfinite(transfer->entropy) && isfinite(transfer->SD[0])
       && isfinite(transfer->SD[1]) && isfinite(transfer->SD[2]))) {
    return ghl_error_invalid_face_transfer;
  }

  // Eq. 10 of Fields et al. (2025) in integrated form: one transfer, two cells.
  delta_L->rho -= transfer->rho;
  delta_R->rho += transfer->rho;
  delta_L->tau -= transfer->tau;
  delta_R->tau += transfer->tau;
  delta_L->Y_e -= transfer->Y_e;
  delta_R->Y_e += transfer->Y_e;
  delta_L->entropy -= transfer->entropy;
  delta_R->entropy += transfer->entropy;
  for(int i = 0; i < 3; i++) {
    delta_L->SD[i] -= transfer->SD[i];
    delta_R->SD[i] += transfer->SD[i];
  }
  return ghl_success;
}
