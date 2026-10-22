#ifndef GHL_CON2PRIM_H_
#define GHL_CON2PRIM_H_

#include "ghl.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @ingroup Con2Prim
 * @struct ghl_con2prim_diagnostics
 * @brief Tracks @ref Con2Prim diagnostics
 *
 * @details
 * Initialize this struct with @ref ghl_initialize_diagnostics before each
 * logical recovery. tau_fix and Stilde_fix are sticky OR accumulators.
 * A successful direct solver writes its own speed_limited result; multi-method
 * drivers accumulate that result across attempts. backup and nn_guess_used
 * record attempted retry paths. which_routine and n_iter describe the
 * successful solver; n_iter is unspecified when which_routine is
 * ghl_con2prim_id_None.
 */
typedef struct ghl_con2prim_diagnostics {
  /** Whether any call limited \f$ \tilde{\tau} \f$ during this recovery */
  bool tau_fix;
  /** Whether any call limited \f$ \tilde{S}_i \f$ during this recovery */
  bool Stilde_fix;
  /** Direct-solver speed-limit result, accumulated by multi-method drivers */
  bool speed_limited;
  /** The Con2Prim routine which successfully found the primitive variables */
  ghl_con2prim_id_t which_routine;
  /** Whether a given backup routine was used (true) or not (false) */
  bool backup[3];
  /** Iterations used by which_routine. Noble totals all attempted outer,
   *  nested, retried, and follow-up Newton iterations; each invocation has its
   *  own iteration budget. Font1D totals outer iterations across its internal
   *  attempts and reports zero when its conservative momentum norm satisfies
   *  S_i gamma^ij S_j < 1e-300. */
  int n_iter;
  /** Whether the tabulated multi-method driver attempted a neural-network retry */
  bool nn_guess_used;
} ghl_con2prim_diagnostics;

/**
 * @ingroup c2p_internal
 * @brief Legacy Palenzuela workspace retained for source and ABI compatibility.
 *
 * The implementation uses private workspace types from Con2Prim/roots.h.
 * No supported solver API consumes this type; external users may retain its
 * existing layout without relying on it as the solver's active workspace.
 */
typedef struct ghl_palenzuela_quantities {
  double q, r, s, t, D, Y_e, S;
} ghl_palenzuela_quantities;

/**
 * @ingroup c2p_internal
 * @brief Stores shared auxiliary quantities for tabulated primitive guesses.
 */
typedef struct ghl_tabulated_primitive_guess_aux {
  double SU[3];
  double B_squared;
  double S_squared;
  double BdotS;
  double q, r, s, t;
} ghl_tabulated_primitive_guess_aux;

//--------- Initialization routines ----------------
void ghl_initialize_diagnostics(ghl_con2prim_diagnostics *restrict diagnostics);

//----------- Pre/Post-C2P routines ----------------
/**
 * Apply conservative limits. diagnostics must be initialized. tau_fix and
 * Stilde_fix are [in,out] OR accumulators: incoming true values are preserved.
 */
void ghl_apply_conservative_limits(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_primitive_quantities *restrict prims,
      ghl_conservative_quantities *restrict cons,
      ghl_con2prim_diagnostics *restrict diagnostics);

void ghl_undensitize_conservatives(
      const double psi6,
      const ghl_conservative_quantities *restrict cons,
      ghl_conservative_quantities *restrict cons_undens);

/**
 * @ingroup c2p_guess
 * @brief Build an initial guess from undensitized conservative variables.
 *
 * For tabulated EOS, atmosphere fallback preserves magnetic components and
 * uses bounded cons_undens->Y_e / cons_undens->rho when density is finite and
 * positive and the ratio is finite; otherwise it uses atmospheric composition.
 * A guess is a solver seed, not a successful recovery.
 */
void ghl_guess_primitives(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims);

/**
 * Enforce primitive limits and compute u0. speed_limited is an initialized
 * [in,out] OR accumulator: incoming true is preserved; start each logical
 * group false.
 */
ghl_error_codes_t ghl_enforce_primitive_limits_and_compute_u0(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      ghl_primitive_quantities *restrict prims,
      bool *restrict speed_limited);

/**
 * @ingroup Con2Prim
 * @brief Compute densitized conservatives and covariant stress-energy.
 *
 * Writes only Tmunu->T4[mu][nu] for mu <= nu; the lower triangle is untouched.
 * Use the symmetric access convention documented in ghl_stress_energy.
 */
void ghl_compute_conservs_and_Tmunu(
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_primitive_quantities *restrict prims,
      ghl_conservative_quantities *restrict cons,
      ghl_stress_energy *restrict Tmunu);

void ghl_compute_conservs(
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_primitive_quantities *restrict prims,
      ghl_conservative_quantities *restrict cons);

//-------------- Con2Prim routines -----------------
/**
 * @ingroup Con2Prim
 * @brief Recover hybrid primitives using the configured primary and backup methods.
 * See the common input, output, and initialization contract in @ref Con2Prim.
 * @return ghl_success on recovery; otherwise the last attempted method's
 *         error, or ghl_error_invalid_c2p_key for an unsupported selection.
 */
ghl_error_codes_t ghl_con2prim_hybrid_multi_method(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

/**
 * @ingroup Con2Prim
 * @brief Recover tabulated primitives using configured methods and optional NN guesses.
 * See the common input, output, and initialization contract in @ref Con2Prim.
 * @return ghl_success on recovery; otherwise the last attempted method's
 *         error, or ghl_error_invalid_c2p_key for an unsupported selection.
 *         Without HDF5, returns ghl_error_used_disabled_hdf5.
 */
ghl_error_codes_t ghl_con2prim_tabulated_multi_method(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

/**
 * @ingroup Con2Prim
 * @brief Dispatch a single hybrid recovery method on undensitized conservatives.
 * See the common input, output, and initialization contract in @ref Con2Prim.
 * @return ghl_success on recovery; otherwise the selected method's
 *         error, or ghl_error_invalid_c2p_key for an unsupported selection.
 */
ghl_error_codes_t ghl_con2prim_hybrid_select_method(
      const ghl_con2prim_id_t c2p_key,
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics );

/**
 * @ingroup Con2Prim
 * @brief Dispatch a single tabulated recovery method on undensitized conservatives.
 * See the common input, output, and initialization contract in @ref Con2Prim.
 * @return ghl_success on recovery; otherwise the selected method's
 *         error, or ghl_error_invalid_c2p_key for an unsupported selection.
 *         Without HDF5, returns ghl_error_used_disabled_hdf5.
 */
ghl_error_codes_t ghl_con2prim_tabulated_select_method(
      const ghl_con2prim_id_t c2p_key,
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics );

/**
 * @ingroup hyb_c2p
 * @brief Recover hybrid primitives with the Noble two-variable energy solve.
 */
ghl_error_codes_t ghl_hybrid_Noble2D(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

/**
 * @ingroup hyb_c2p
 * @brief Recover hybrid primitives with the Noble one-variable energy solve.
 */
ghl_error_codes_t ghl_hybrid_Noble1D(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

/**
 * @ingroup hyb_c2p
 * @brief Recover primitives with the Noble entropy solve for a one-piece, equal-Gamma EOS.
 */
ghl_error_codes_t ghl_hybrid_Noble1D_entropy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

/**
 * @ingroup hyb_c2p
 * @brief Recover primitives with the Noble density-based entropy solve.
 *
 * Density-based entropy recovery for a one-piece constant-Gamma EOS, requiring
 * neos == 1 and Gamma_th == Gamma_ppoly[0]. A successful direct solve returns
 * its positive mathematical root, which may lie outside configured EOS density
 * limits. Call ghl_enforce_primitive_limits_and_compute_u0 afterward to apply
 * those limits.
 */
ghl_error_codes_t ghl_hybrid_Noble1D_entropy2(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

/**
 * @ingroup hyb_c2p
 * @brief Recover hybrid primitives with the Palenzuela energy solve.
 * @param[in] params Recovery and velocity-limit parameters
 * @param[in] eos Initialized EOS parameters
 * @param[in] metric_adm ADM metric
 * @param[in] metric_aux Unused by this method; retained for the common interface
 * @param[in] cons_undens Undensitized conservatives
 * @param[in,out] prims Supplies magnetic fields and receives recovered primitives
 * @param[in,out] diagnostics Initialized recovery diagnostics
 * @return ghl_success on recovery; otherwise a root-finding or EOS error.
 */
ghl_error_codes_t ghl_hybrid_Palenzuela1D_energy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

/**
 * @ingroup hyb_c2p
 * @brief Recover hybrid primitives with the Palenzuela entropy solve.
 * @param[in] params Recovery and velocity-limit parameters
 * @param[in] eos Initialized EOS parameters
 * @param[in] metric_adm ADM metric
 * @param[in] metric_aux Unused by this method; retained for the common interface
 * @param[in] cons_undens Undensitized conservatives
 * @param[in,out] prims Supplies magnetic fields and receives recovered primitives
 * @param[in,out] diagnostics Initialized recovery diagnostics
 * @return ghl_success on recovery; otherwise a root-finding or EOS error.
 */
ghl_error_codes_t ghl_hybrid_Palenzuela1D_entropy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

/**
 * @ingroup hyb_c2p
 * @brief Recover hybrid primitives with the Font cold-EOS fallback.
 */
ghl_error_codes_t ghl_hybrid_Font1D(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

/**
 * @ingroup tab_c2p
 * @brief Recover tabulated primitives with the Noble two-variable energy solve.
 */
ghl_error_codes_t ghl_tabulated_Noble2D(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

/**
 * @ingroup tab_c2p
 * @brief Recover tabulated primitives with the Palenzuela energy solve.
 * @param[in] params Recovery and velocity-limit parameters
 * @param[in] eos Initialized EOS parameters
 * @param[in] metric_adm ADM metric
 * @param[in] metric_aux Unused by this method; retained for the common interface
 * @param[in] cons_undens Undensitized conservatives
 * @param[in,out] prims Supplies magnetic fields and receives recovered primitives
 * @param[in,out] diagnostics Initialized recovery diagnostics
 * @return ghl_success on recovery; otherwise a root-finding or EOS error.
 *
 * Direct callers supply prims->Y_e and a temperature seed. The incoming
 * composition is bounded by the EOS, not recomputed from conservatives. When
 * params->evolve_temp is false, supply the desired fixed temperature; an
 * internal retry may reset the seed to eos->T_min.
 */
ghl_error_codes_t ghl_tabulated_Palenzuela1D_energy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

/**
 * @ingroup tab_c2p
 * @brief Recover tabulated primitives with the Palenzuela entropy solve.
 * @param[in] params Recovery and velocity-limit parameters
 * @param[in] eos Initialized EOS parameters
 * @param[in] metric_adm ADM metric
 * @param[in] metric_aux Unused by this method; retained for the common interface
 * @param[in] cons_undens Undensitized conservatives
 * @param[in,out] prims Supplies magnetic fields and receives recovered primitives
 * @param[in,out] diagnostics Initialized recovery diagnostics
 * @return ghl_success on recovery; otherwise a root-finding or EOS error.
 *
 * Direct callers supply prims->Y_e and a temperature seed. The incoming
 * composition is bounded by the EOS, not recomputed from conservatives. When
 * params->evolve_temp is false, supply the desired fixed temperature; an
 * internal retry may reset the seed to eos->T_min.
 */
ghl_error_codes_t ghl_tabulated_Palenzuela1D_entropy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

/**
 * @ingroup tab_c2p
 * @brief Recover tabulated primitives with the Newman energy solve.
 */
ghl_error_codes_t ghl_tabulated_Newman1D_energy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

/**
 * @ingroup tab_c2p
 * @brief Recover tabulated primitives with the Newman entropy solve.
 */
ghl_error_codes_t ghl_tabulated_Newman1D_entropy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

//------------ Auxiliary Functions -----------------
bool ghl_limit_utilde_and_compute_v(
      const ghl_parameters *restrict params,
      const ghl_metric_quantities *restrict metric_adm,
      double utU[3],
      ghl_primitive_quantities *restrict prims);

/**
 * @ingroup c2p_internal
 * @brief Prepare momentum/magnetic contractions for tabulated primitive guesses.
 * @param[in] metric_adm ADM metric
 * @param[in] cons_undens Undensitized conservatives; finite positive density required
 * @param[in] prims Supplies the magnetic field
 * @param[out] aux Locally limited momentum, contractions, and normalized q,r,s,t
 *
 * Uses ghl_compute_SU_Bsq_Ssq_BdotS() without changing the conservative input.
 * This algebraic helper does not validate inputs or test EOS admissibility.
 */
void ghl_tabulated_compute_primitive_guess_auxiliaries(
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons_undens,
      const ghl_primitive_quantities *restrict prims,
      ghl_tabulated_primitive_guess_aux *restrict aux);

/**
 * @ingroup c2p_guess
 * @brief Complete a tabulated primitive initial guess from a Palenzuela variable.
 *
 * Complete a tabulated primitive initial guess from a finite Palenzuela
 * variable `x`. This helper accepts any algebraically usable finite `x`; its
 * sign is not an EOS admissibility test. In particular, tabulated EOS energy
 * reference conventions may admit negative internal energy and conservative
 * energy. Admissibility depends on EOS bounds local to density and
 * composition. This routine supplies a solver seed, not that decision.
 * Nonfinite `x`, invalid conservative input, or unusable intermediate/output
 * values return the initialized EOS atmosphere while preserving incoming
 * magnetic components. For finite positive conservative density and a finite
 * cons_undens->Y_e / cons_undens->rho, fallback retains that composition within
 * EOS bounds; otherwise it uses atmospheric composition. If the EOS inversion alone fails after the algebraic
 * state is valid, the helper preserves its bounded `rho`, `Y_e`, and `eps`,
 * uses `T_max` as the temperature seed, and supplies finite atmosphere
 * pressure and entropy before completing the velocity seed.
 *
 * @param[in] params Con2Prim parameters
 * @param[in] eos tabulated EOS parameters
 * @param[in] metric_adm ADM metric
 * @param[in] cons_undens undensitized conservatives
 * @param[in] aux precomputed algebraic contractions
 * @param[in] x finite Palenzuela variable; sign alone does not determine EOS
 *              admissibility
 * @param[in,out] prims primitive guess; incoming magnetic components are
 *                      preserved
 */
void ghl_tabulated_primitive_guess_from_x(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons_undens,
      const ghl_tabulated_primitive_guess_aux *restrict aux,
      double x,
      ghl_primitive_quantities *restrict prims);

/**
 * @ingroup Con2Prim
 * @brief EOS-selected multi-method driver; initialize EOS dispatch before use.
 * Accepts undensitized conservatives under the common @ref Con2Prim contract.
 */
extern ghl_error_codes_t (*ghl_con2prim_multi_method)(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

//------------ Neural Network Initial Guesses -----------------
/**
 * @ingroup c2p_guess
 * @brief Fixed NN input layout {q, r, s, t} from tabulated guess auxiliaries.
 * With D the undensitized conserved density, q=tau/D, r=S^2/D^2, s=B^2/D,
 * and t=(B.S)/D^(3/2), using the helper's locally limited momentum.
 */
typedef struct ghl_nn_c2p_input_t {
  float q;
  float r;
  float s;
  float t;
} ghl_nn_c2p_input_t;

/**
 * @ingroup c2p_guess
 * @brief Candidate Palenzuela variable x=hW; a seed rather than recovered primitives.
 */
typedef struct ghl_nn_c2p_guess_t {
  float x;
} ghl_nn_c2p_guess_t;

/**
 * @ingroup c2p_guess
 * @brief NN dimensions, preprocessing arrays, and layer weights for an x guess.
 *
 * Loaders allocate the model and its arrays, validate them, and transfer ownership
 * to eos->c2p_nn on success, freeing the previous model. EOS cleanup releases
 * that model. A manually constructed model must have correctly sized accessible
 * arrays and pass ghl_c2p_nn_validate_model() before inference. It may be passed
 * to ghl_c2p_nn_free() only if the struct and every non-NULL array can be freed.
 */
struct ghl_c2p_nn_model {
  int in_dim;
  int hidden_dim;
  int n_hidden;
  int out_dim;
  /** Fixed input layout is {q, r, s, t}; q_idx must be 0. */
  int q_idx;
  /** Fixed input layout is {q, r, s, t}; s_idx must be 2. */
  int s_idx;
  float x_eps;
  float y_eps;
  float dx_eps;

  int *x_kind;
  float *x_lo;
  float *x_hi;

  int *out_kind;
  float *x_invrng;
  float *out_lo;
  float *out_hi;
  float *out_invrng;

  float *W_in;
  float *b_in;
  float *W_hid;
  float *b_hid;
  float *W_out;
  float *b_out;
};

/** @ingroup c2p_guess
 * Informational version for the expected HDF5 schema/output-kind semantics.
 * The current loader validates fields but does not negotiate this value from
 * an on-disk version dataset. */
#define GHL_NN_C2P_API_VERSION 4u

/**
 * @ingroup c2p_guess
 * @brief Infer a candidate x from the fixed {q, r, s, t} input layout.
 * @param[in] model Validated model, or NULL; borrowed without ownership transfer
 * @param[in] input Normalized conservative/magnetic contractions
 * @return Candidate x, falling back to the bracket midpoint when usable or zero
 *         when a model/bracket is unavailable. No recovery status is implied.
 *
 * Inference does not validate model storage; callers must validate it beforehand.
 * This routine is available without HDF5.
 */
ghl_nn_c2p_guess_t ghl_c2p_nn_guess(
      const ghl_c2p_nn_model *restrict model,
      ghl_nn_c2p_input_t input);

/**
 * @ingroup c2p_guess
 * @brief Build a tabulated primitive initial guess using a validated NN model.
 *
 * Build a tabulated primitive initial guess. A present model must already be
 * validated and uses fixed {q, r, s, t} inference. A missing model or unusable
 * numerical state/candidate returns the atmosphere initial guess, preserving
 * magnetic components and bounded conservative composition when density is
 * finite and positive and the composition ratio is finite; otherwise it uses
 * atmospheric composition. Arguments must be non-NULL; parameters, the tabulated
 * EOS atmosphere, metric, and incoming magnetic fields must be valid and
 * finite. Inference additionally requires finite usable contractions. This
 * helper does not validate corrupted EOS/model objects, and its fallback is
 * not a recovery result.
 *
 * @param[in] params Recovery and velocity-limit parameters
 * @param[in] eos Tabulated EOS with atmosphere, bounds, and optional validated model
 * @param[in] metric_adm ADM metric
 * @param[in] cons_undens Undensitized conservatives
 * @param[in,out] prims Supplies magnetic fields and receives the primitive seed
 */
void ghl_c2p_nn_guess_primitives(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims);

/**
 * @ingroup c2p_guess
 * @brief Validate model metadata and numerical array contents without modifying it.
 * @param[in] model Borrowed model; NULL is rejected
 * @return ghl_success or a ghl_error_nn_c2p_* code describing invalid metadata,
 *         missing arrays, or invalid numerical entries.
 *
 * Array allocations must already match the declared dimensions: validation
 * cannot verify allocation extents. Available without HDF5.
 */
ghl_error_codes_t ghl_c2p_nn_validate_model(const ghl_c2p_nn_model *restrict model);

/**
 * @ingroup c2p_guess
 * @brief Free a heap-allocated model and every owned array; NULL is a no-op.
 * @param[in] model Model whose ownership is consumed
 *
 * Available without HDF5. This does not clear an owner's pointer. If freeing
 * eos->c2p_nn directly, set it to NULL to avoid a later double free by EOS cleanup.
 */
void ghl_c2p_nn_free(ghl_c2p_nn_model *restrict model);

/**
 * @ingroup c2p_guess
 * @brief Load and validate a model from the root of a standalone HDF5 file.
 * @param[in] nn_model_filepath Nonempty file path; the file is opened read-only
 * @param[in,out] eos EOS with c2p_nn initialized to NULL or an owned valid model
 * @return ghl_success, or an EOS-pointer, file, HDF5 dataset, allocation, or
 *         NN validation error. With GHL_DISABLE_HDF5, returns
 *         ghl_error_used_disabled_hdf5 without changing eos.
 *
 * On success, frees the old model and transfers the new model to eos->c2p_nn.
 * On failure, preserves the old model. Does not enable NN retries; configure
 * the associated runtime flags separately.
 */
ghl_error_codes_t ghl_c2p_nn_load_hdf5(
      const char *nn_model_filepath,
      ghl_eos_parameters *restrict eos);

/**
 * @ingroup c2p_guess
 * @brief Load and validate a model from the grhayl_nn_c2p group in an EOS HDF5 file.
 * @param[in] eos_table_filepath Nonempty file path; the file is opened read-only
 * @param[in,out] eos EOS with c2p_nn initialized to NULL or an owned valid model
 * @return ghl_success, or an EOS-pointer, file, HDF5 dataset, allocation, or
 *         NN validation error. With GHL_DISABLE_HDF5, returns
 *         ghl_error_used_disabled_hdf5 without changing eos.
 *
 * On success, frees the old model and transfers the new model to eos->c2p_nn.
 * On failure, preserves the old model. Does not enable NN retries; configure
 * the associated runtime flags separately.
 */
ghl_error_codes_t ghl_c2p_nn_load_from_eos_hdf5(
      const char *eos_table_filepath,
      ghl_eos_parameters *restrict eos);

#ifdef __cplusplus
}
#endif

#endif // GHL_CON2PRIM_H
