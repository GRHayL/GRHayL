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
 * This struct should be initialized with @ref ghl_initialize_diagnostics .
 */
typedef struct ghl_con2prim_diagnostics {
  /** Whether a limit was applied to \f$ \tilde{\tau} \f$ (true) or not (false) */
  bool tau_fix;
  /** Whether a limit was applied to \f$ \tilde{S}_i \f$ (true) or not (false) */
  bool Stilde_fix;
  /** Whether a speed limiter was triggered (true) or not (false) */
  bool speed_limited;
  /** The Con2Prim routine which successfully found the primitive variables */
  ghl_con2prim_id_t which_routine;
  /** Whether a given backup routine was used (true) or not (false) */
  bool backup[3];
  /** Number of iterations required to find the solution */
  int n_iter;
} ghl_con2prim_diagnostics;

/**
 * @ingroup c2p_internal
 * @brief Stores internal variables for Palenzuela1D routine
 */
typedef struct ghl_palenzuela_quantities {
  double q, r, s, t, D, Y_e, S;
} ghl_palenzuela_quantities;

//--------- Initialization routines ----------------
void ghl_initialize_diagnostics(ghl_con2prim_diagnostics *restrict diagnostics);

//----------- Pre/Post-C2P routines ----------------
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

void ghl_guess_primitives(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prims);

ghl_error_codes_t ghl_enforce_primitive_limits_and_compute_u0(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      ghl_primitive_quantities *restrict prims,
      bool *restrict speed_limited);

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
ghl_error_codes_t ghl_con2prim_hybrid_multi_method(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

ghl_error_codes_t ghl_con2prim_tabulated_multi_method(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

ghl_error_codes_t ghl_con2prim_hybrid_select_method(
      const ghl_con2prim_id_t c2p_key,
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics );

ghl_error_codes_t ghl_con2prim_tabulated_select_method(
      const ghl_con2prim_id_t c2p_key,
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics );

ghl_error_codes_t ghl_hybrid_Noble2D(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

ghl_error_codes_t ghl_hybrid_Noble1D(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

ghl_error_codes_t ghl_hybrid_Noble1D_entropy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

ghl_error_codes_t ghl_hybrid_Noble1D_entropy2(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

ghl_error_codes_t ghl_hybrid_Palenzuela1D_energy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

ghl_error_codes_t ghl_hybrid_Palenzuela1D_entropy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

ghl_error_codes_t ghl_hybrid_Font1D(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

ghl_error_codes_t ghl_tabulated_Noble2D(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

ghl_error_codes_t ghl_tabulated_Palenzuela1D_energy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

ghl_error_codes_t ghl_tabulated_Palenzuela1D_entropy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

ghl_error_codes_t ghl_tabulated_Newman1D_energy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

ghl_error_codes_t ghl_tabulated_Newman1D_entropy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

//------------ Auxiliary Functions -----------------
bool ghl_limit_utilde_and_compute_v(
      const ghl_parameters *restrict params,
      const ghl_metric_quantities *restrict metric,
      double utU[3],
      ghl_primitive_quantities *restrict prims);

extern ghl_error_codes_t (*ghl_con2prim_multi_method)(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

//------------ Neural Network Initial Guesses -----------------
typedef struct ghl_nn_c2p_input_t {
  float q;
  float r;
  float s;
  float t;
} ghl_nn_c2p_input_t;

typedef struct ghl_nn_c2p_guess_t {
  float x;
  float W;
} ghl_nn_c2p_guess_t;

typedef struct ghl_c2p_nn_model {
  int in_dim;
  int hidden_dim;
  int n_hidden;
  int out_dim;
  int q_idx;
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
} ghl_c2p_nn_model;

/* Public API version for the on-disk HDF5 schema/output-kind semantics. */
#define GHL_NN_C2P_API_VERSION 3u

ghl_nn_c2p_guess_t ghl_c2p_nn_guess(
      const ghl_c2p_nn_model *restrict model,
      ghl_nn_c2p_input_t input);

ghl_error_codes_t ghl_c2p_nn_validate_model(const ghl_c2p_nn_model *restrict model);

void ghl_c2p_nn_free(ghl_c2p_nn_model *restrict model);

ghl_error_codes_t ghl_c2p_nn_load_hdf5(
      const char *nn_model_filepath,
      ghl_eos_parameters *restrict eos);

ghl_error_codes_t ghl_c2p_nn_load_from_eos_hdf5(
      const char *eos_table_filepath,
      ghl_eos_parameters *restrict eos);

#ifdef __cplusplus
}
#endif

#endif // GHL_CON2PRIM_H
