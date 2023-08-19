#ifndef CON2PRIM_H_
#define CON2PRIM_H_

#include "ghl.h"

//------------- Con2Prim struct --------------------
/*
   The struct ghl_con2prim_diagnostics contains variables for error-checking and
   diagnostic feedback. The struct elements are detailed below:

 --TODO

TODO: consider changing failure_checker to be bitwise; failure modes are currently
      1: atmosphere reset when rho_star < 0
      10: reseting P when P<P_min in enforce_...
      100: reseting P when P>P_max in enforce_...
      1k: Limiting velocity u~ after C2P/Font Fix or v in enforce_...
      10k: Font Fix was applied
      100k: Both C2P and Font Fix failed
      1M: tau~ was reset in ghl_apply_conservative_limits
      10M: S~ was reset in ghl_apply_conservative_limits via the first case
      100M: S~ was reset in ghl_apply_conservative_limits via the second case
For bitwise, would become 1, 2, 4, 8, 16, 32. 64, 128, and 256
https://www.tutorialspoint.com/cprogramming/c_bitwise_operators.htm
*/

typedef struct ghl_con2prim_diagnostics {
  int tau_fix;
  int Stilde_fix;
  int speed_limited;
  int which_routine;
  bool backup[3];
  int n_iter;
} ghl_con2prim_diagnostics;

typedef struct ghl_palenzuela_quantities {
  double q, r, s, t, D, Y_e, S;
} ghl_palenzuela_quantities;

//--------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif

//--------- Initialization routines ----------------

void ghl_initialize_diagnostics(ghl_con2prim_diagnostics *restrict diagnostics);

//----------- Pre/Post-C2P routines ----------------

void ghl_apply_conservative_limits(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_primitive_quantities *restrict prims,
      ghl_conservative_quantities *restrict cons,
      ghl_con2prim_diagnostics *restrict diagnostics);

void ghl_undensitize_conservatives(
      const double psi6,
      const ghl_conservative_quantities *restrict cons,
      ghl_conservative_quantities *restrict cons_undens);

void ghl_guess_primitives(
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prims);

int ghl_enforce_primitive_limits_and_compute_u0(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      ghl_primitive_quantities *restrict prims);

void ghl_compute_conservs_and_Tmunu(
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_primitive_quantities *restrict prims,
      ghl_conservative_quantities *restrict cons,
      ghl_stress_energy *restrict Tmunu);

void ghl_compute_conservs(
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_primitive_quantities *restrict prims,
      ghl_conservative_quantities *restrict cons);

//--------------------------------------------------

//-------------- Con2Prim routines -----------------
int ghl_con2prim_multi_method(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

int ghl_con2prim_select_method(
      const ghl_con2prim_method_t c2p_key,
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

int ghl_hybrid_Noble2D(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

int ghl_hybrid_Font_fix(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

int ghl_tabulated_Palenzuela1D_energy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

int ghl_tabulated_Palenzuela1D_entropy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

int ghl_tabulated_Newman1D_energy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

int ghl_tabulated_Newman1D_entropy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

//--------------------------------------------------

//------------ Auxiliary Functions -----------------

int ghl_limit_utilde_and_compute_v(
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric,
      double utU[3],
      ghl_primitive_quantities *restrict prims);

#ifdef __cplusplus
}
#endif

#endif // CON2PRIM_H
