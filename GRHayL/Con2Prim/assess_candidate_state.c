#include <float.h>

#include "ghl_con2prim.h"

/**
 * @ingroup Con2Prim
 * @brief Assesses whether recovering a candidate conservative state needs a repair
 *
 * @details
 * This is the assessment step of the first-order flux correction (FOFC) of
 * \cite Lemaster_2009, in the form used by AthenaK (\cite Stone_2024_AthenaK;
 * \cite Fields_2025, Sec. 3.3).
 * There, the next Runge-Kutta step is estimated and a cell is flagged if its new
 * state requires a floor or the conserved-to-primitive inversion fails. The
 * estimate is only a test: no source terms are applied and no matter field is
 * updated, so this function changes nothing the caller owns. The caller then
 * recomputes the fluxes of the flagged cells with first-order (donor-cell)
 * reconstruction and a more dissipative flux, LLF or HLLE in \cite Fields_2025.
 * The relaxed discrete maximum principle of \cite Fields_2025 (their Eq. 12) is not
 * applied here.
 *
 * The function runs the ordinary recovery sequence on copies of its inputs and
 * sets flagged when the candidate
 *
 * - has a \f$ \tilde{D} \f$, \f$ \tilde{\tau} \f$, or \f$ \tilde{S}_i \f$ that is not
 *   finite, or an entropy that is not finite when ghl_parameters::evolve_entropy is
 *   true, or is too large for the closure bound below to be formed,
 * - fails recovery numerically,
 * - made @ref ghl_apply_conservative_limits (the floors of \cite Faber_2007 and
 *   \cite Etienne_2012) or a velocity limit act
 *   (ghl_con2prim_diagnostics::tau_fix, ghl_con2prim_diagnostics::Stilde_fix,
 *   ghl_con2prim_diagnostics::speed_limited),
 * - was recovered by Font1D (\cite Font_2000), which replaces the energy equation
 *   by the cold EOS, or
 * - does not close: the conservatives rebuilt from the limited primitives differ
 *   from the original candidate by more than the closure bound,
 *   \f[
 *     \left| \tilde{u}_k^\mathrm{rebuilt} - \tilde{u}_k^\mathrm{candidate} \right|
 *     \le \max(\epsilon_\mathrm{C2P}, \epsilon_\mathrm{machine})
 *          \left( |\tilde{D}| + |\tilde{\tau}| \right),
 *   \f]
 *   for \f$ \tilde{u}_k = \tilde{D}, \tilde{\tau}, \tilde{S}_i \f$ and, for a
 *   tabulated EOS, \f$ \tilde{Y}_e \f$, where \f$ \epsilon_\mathrm{C2P} \f$ is
 *   ghl_parameters::con2prim_solver_tolerance. Primitive limits and solver
 *   clamps above the bound act through this comparison.
 *
 * With \f$ Q^\mathrm{test}_i \f$ the integrated candidate and \f$ V_i \f$ the
 * coordinate volume of the cell, cons_candidate is
 * \f$ \tilde{U}^\mathrm{test}_i = Q^\mathrm{test}_i / V_i \f$, and the undensitized
 * state
 * \f$ U^\mathrm{test}_i = \tilde{U}^\mathrm{test}_i / \sqrt{\gamma_i} \f$
 * is formed once, with metric_adm->sqrt_detgamma, by @ref ghl_undensitize_conservatives.
 * metric_adm and metric_aux describe the cell at the time of the candidate.
 * prims_guess supplies the magnetic field and, when
 * ghl_parameters::calc_prim_guess is false, the initial guess. A negative
 * \f$ \tilde{\tau} \f$ is not by itself a trigger for a tabulated EOS, which defines
 * its zero point through the table. For a tabulated EOS
 * @ref ghl_apply_conservative_limits is not applied, since its \f$ \tau_\mathrm{atm} \f$
 * can be negative, which makes the limiter's momentum rescaling NaN.
 *
 * Closure detects a floor, limit, or clamp only if it changes the state by more than the
 * closure bound. A smaller one is flagged only if one of the diagnostics listed above
 * also records it; a primitive floor with no such record, such as a pressure floor in a
 * nearly cold cell, goes unflagged. The evolved entropy is checked for finiteness but
 * not compared with the energy, because entropy generation at a shock is not a failure
 * of the energy recovery. With an entropy-based solver the entropy drives the recovery,
 * so a disagreement with the energy appears as a closure failure.
 *
 * Configuration errors are returned and leave flagged unchanged: disabled HDF5 for a
 * tabulated EOS, and an unknown or invalid EOS type or an invalid solver key that is the
 * final result of the recovery sequence. An error that a later backup recovers from is
 * not reported, a solver the sequence never calls is not checked, and a candidate that
 * is not finite is flagged before any solver runs, so validate the whole configured
 * solver list when initializing.
 *
 * @param[in] params pointer to ghl_parameters struct
 *
 * @param[in] eos pointer to ghl_eos_parameters struct
 *
 * @param[in] metric_adm pointer to ghl_metric_quantities struct with the ADM metric
 *                       of the cell
 *
 * @param[in] metric_aux pointer to ghl_ADM_aux_quantities struct computed from
 *                       metric_adm
 *
 * @param[in] cons_candidate pointer to ghl_conservative_quantities struct with
 *                           the **densitized** candidate per coordinate volume
 *
 * @param[in] prims_guess pointer to ghl_primitive_quantities struct with the
 *                        magnetic field and the initial guess
 *
 * @param[out] flagged true if the candidate needs a repair
 *
 * @returns ghl_success, or the configuration error
 */
ghl_error_codes_t ghl_assess_candidate_state(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_candidate,
      const ghl_primitive_quantities *restrict prims_guess,
      bool *restrict flagged) {

  const bool tabulated = eos->eos_type == ghl_eos_tabulated;
  // A tabulated EOS needs HDF5, whatever the candidate.
#ifdef GHL_DISABLE_HDF5
  if(tabulated) {
    return ghl_error_used_disabled_hdf5;
  }
#endif

  // An infinite candidate, or a sum |D| + |tau| that overflows, would make the closure
  // bound below infinite. An energy-based recovery would also discard a non-finite
  // evolved entropy, so it is checked here.
  const double scale = fabs(cons_candidate->rho) + fabs(cons_candidate->tau);
  if(!(isfinite(scale) && isfinite(cons_candidate->SD[0])
       && isfinite(cons_candidate->SD[1]) && isfinite(cons_candidate->SD[2])
       && (!params->evolve_entropy || isfinite(cons_candidate->entropy)))) {
    *flagged = true;
    return ghl_success;
  }

  // As in first-order flux correction (Lemaster & Stone 2009; Fields et al. 2025,
  // Sec. 3.3), the estimated state is only tested: everything below acts on copies,
  // applies no source terms, and updates no matter field.
  ghl_primitive_quantities prims = *prims_guess;
  ghl_conservative_quantities cons_limited = *cons_candidate;
  ghl_con2prim_diagnostics diagnostics;
  ghl_initialize_diagnostics(&diagnostics);

  // Conservative floors of Faber et al. (2007) and Etienne et al. (2012).
  if(!tabulated) {
    ghl_apply_conservative_limits(
          params, eos, metric_adm, &prims, &cons_limited, &diagnostics);
  }
  ghl_conservative_quantities cons_undens;
  ghl_undensitize_conservatives(metric_adm->sqrt_detgamma, &cons_limited, &cons_undens);

  ghl_error_codes_t error;
  if(tabulated) {
    error = ghl_con2prim_tabulated_multi_method(
          params, eos, metric_adm, metric_aux, &cons_undens, &prims, &diagnostics);
  }
  else {
    error = ghl_con2prim_hybrid_multi_method(
          params, eos, metric_adm, metric_aux, &cons_undens, &prims, &diagnostics);
  }
  if(error == ghl_success) {
    error = ghl_enforce_primitive_limits_and_compute_u0(
          params, eos, metric_adm, &prims, &diagnostics.speed_limited);
  }

  // Fields et al. (2025), Sec. 3.3: a cell is flagged if the conserved-to-primitive
  // inversion fails. A configuration error is returned instead; any other failure
  // flags the candidate.
  if(error == ghl_error_unknown_eos_type || error == ghl_error_invalid_c2p_key
     || error == ghl_error_invalid_eos_type) {
    return error;
  }
  if(error != ghl_success) {
    *flagged = true;
    return ghl_success;
  }

  // Fields et al. (2025), Sec. 3.3: a cell is flagged if its state requires a floor.
  // Closure sees a floor, limit, or clamp that changes the state by more than the bound.
  // A smaller one is flagged only if a diagnostic below records it (the cold-EOS energy
  // of Font1D, Font et al. 2000, by name); a primitive floor with no record goes
  // unflagged. Rebuild from the limited primitives and compare with the original
  // candidate. tau and S_i are rebuilt from terms of the size of D, so one scale serves
  // every component.
  ghl_conservative_quantities cons_rebuilt;
  ghl_compute_conservs(metric_adm, metric_aux, &prims, &cons_rebuilt);
  const double bound = fmax(params->con2prim_solver_tolerance, DBL_EPSILON) * scale;
  const double expected[6]
        = { cons_candidate->rho,   cons_candidate->tau,   cons_candidate->SD[0],
            cons_candidate->SD[1], cons_candidate->SD[2], cons_candidate->Y_e };
  const double actual[6] = { cons_rebuilt.rho,   cons_rebuilt.tau,   cons_rebuilt.SD[0],
                             cons_rebuilt.SD[1], cons_rebuilt.SD[2], cons_rebuilt.Y_e };
  bool closes = true;
  for(int i = 0; i < (tabulated ? 6 : 5); i++) {
    closes &= fabs(actual[i] - expected[i]) <= bound;
  }

  *flagged = diagnostics.tau_fix || diagnostics.Stilde_fix || diagnostics.speed_limited
             || diagnostics.which_routine == ghl_con2prim_id_Font1D || !closes;
  return ghl_success;
}
