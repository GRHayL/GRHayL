#ifndef CON2PRIM_H_
#define CON2PRIM_H_

#include "ghl.h"

//------------- Con2Prim struct --------------------
/*
   The struct con2prim_diagnostics contains variables for error-checking and
   diagnostic feedback. The struct elements are detailed below:

 --TODO

TODO: consider changing failure_checker to be bitwise; failure modes are currently
      1: atmosphere reset when rho_star < 0
      10: reseting P when P<P_min in enforce_...
      100: reseting P when P>P_max in enforce_...
      1k: Limiting velocity u~ after C2P/Font Fix or v in enforce_...
      10k: Font Fix was applied
      100k: Both C2P and Font Fix failed
      1M: tau~ was reset in ghl_apply_inequality_fixes
      10M: S~ was reset in ghl_apply_inequality_fixes via the first case
      100M: S~ was reset in ghl_apply_inequality_fixes via the second case
For bitwise, would become 1, 2, 4, 8, 16, 32. 64, 128, and 256
https://www.tutorialspoint.com/cprogramming/c_bitwise_operators.htm
*/

typedef enum {
  None=-1,
  Noble2D, Noble1D,
  Noble1D_entropy, Noble1D_entropy2,
  FontFix, CerdaDuran2D, CerdaDuran3D,
  Palenzuela1D, Palenzuela1D_entropy,
  Newman1D, Newman1D_entropy
} con2prim_method_t;

typedef struct con2prim_diagnostics {
  bool c2p_failed;
  int failures;
  int failure_checker;
  int speed_limited;
  int which_routine;
  int backup[3];
  int nan_found;
  int c2p_fail_flag;
  double error_int_numer;
  double error_int_denom;
  int n_iter;
} con2prim_diagnostics;

typedef struct palenzuela_quantities {
  double q, r, s, t, D, Y_e, S;
} palenzuela_quantities;

//--------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif

//--------- Initialization routines ----------------

void ghl_initialize_diagnostics(con2prim_diagnostics *restrict diagnostics);

//----------- Pre/Post-C2P routines ----------------

void ghl_apply_inequality_fixes(
      const ghl_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const primitive_quantities *restrict prims,
      conservative_quantities *restrict cons,
      con2prim_diagnostics *restrict diagnostics);

void ghl_undensitize_conservatives(
      const double psi6,
      const conservative_quantities *restrict cons,
      conservative_quantities *restrict cons_undens);

void ghl_guess_primitives(
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prims);

void ghl_enforce_primitive_limits_and_compute_u0(
      const ghl_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      primitive_quantities *restrict prims,
      int *restrict speed_limit);

void ghl_compute_conservs_and_Tmunu(
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const primitive_quantities *restrict prims,
      conservative_quantities *restrict cons,
      stress_energy *restrict Tmunu);

void ghl_compute_conservs(
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const primitive_quantities *restrict prims,
      conservative_quantities *restrict cons);

//--------------------------------------------------

//-------------- Con2Prim routines -----------------
int ghl_con2prim_multi_method(
      const ghl_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prim,
      con2prim_diagnostics *restrict diagnostics);

int ghl_con2prim_select_method(
      const con2prim_method_t c2p_key,
      const ghl_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prim,
      con2prim_diagnostics *restrict diagnostics);

int ghl_hybrid_Noble2D(
      const ghl_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prim,
      con2prim_diagnostics *restrict diagnostics);

int ghl_hybrid_Font_fix(
      const ghl_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const conservative_quantities *restrict cons_undens,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics);

int ghl_tabulated_Palenzuela1D_energy(
      const ghl_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prim,
      con2prim_diagnostics *restrict diagnostics);

int ghl_tabulated_Palenzuela1D_entropy(
      const ghl_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prim,
      con2prim_diagnostics *restrict diagnostics);

int ghl_tabulated_Newman1D_energy(
      const ghl_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prim,
      con2prim_diagnostics *restrict diagnostics);

int ghl_tabulated_Newman1D_entropy(
      const ghl_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prim,
      con2prim_diagnostics *restrict diagnostics);

//--------------------------------------------------

//------------ Auxiliary Functions -----------------

void ghl_limit_utilde_and_compute_v(
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      double utU[3],
      primitive_quantities *restrict prims,
      int *restrict speed_limit);

#ifdef __cplusplus
}
#endif

#endif // CON2PRIM_H
