#include "ghl.h"

/**
 * @ingroup pack_struct
 * @brief Initialize the @grhayl parameter struct from user input
 *
 * @details
 * This function sets all the elements of the ghl_parameters
 * struct based on the inputs. These are used to control the
 * behavior of various aspects of the GRHayL library. In addition
 * to the input parameters, several others are set to their defaults.
 *
 * The Con2Prim parameters are set to the defaults of
 * ```
 * con2prim_max_iterations   = 30
 * con2prim_solver_tolerance = 1e-10
 * ```
 *
 * The PPM reconstruction parameters are set to the defaults of
 * ```
 * ppm_flattening_epsilon = 0.33
 * ppm_flattening_omega1  = 0.75
 * ppm_flattening_omega2  = 10.0
 * ppm_shock_k0           = 0.1
 * ppm_shock_eta1         = 20.0
 * ppm_shock_eta2         = 0.05
 * ppm_shock_epsilon      = 0.01
 * ```
 * which come from the original Colella and Woodward [paper](https://www.sciencedirect.com/science/article/abs/pii/0021999184901438?via%3Dihub).
 *
 * @param[in] main_routine:          selects the primary conservative-to-primitive routine for @ref ghl_con2prim_multi_method
 *                                   function; options are limited to @ref ghl_con2prim_method_t
 *
 * @param[in] backup_routine:        selects backup conservative-to-primitive routines for @ref ghl_con2prim_multi_method
 *                                   function; up to 3 backups can be selected, with no backup being -1; options are
 *                                   limited to @ref ghl_con2prim_method_t
 *
 * @param[in] evolve_entropy:        whether entropy should be evolved (True) or not (False)
 *
 * @param[in] evolve_temp:           whether temperature should be evolved (True) or not (False)
 *
 * @param[in] calc_prim_guess:       sets whether the provided ghl_primitive_quantities struct contains an initial
 *                                   Con2Prim guess for the ghl_con2prim_multi_method function
 *
 * @param[in] psi6threshold:         upper limit of \f$ \psi^6 = \sqrt{|\gamma|} \f$ above which the limits on conservatives and primitives are adjusted
 *
 * @param[in] max_Lorentz_factor:    maximum allowed Lorenz factor \f$ W \f$ in the simulation
 *
 * @param[in] Lorenz_damping_factor: sets the damping factor for the Lorenz gauge term in \f$ \tilde{\Phi}^\mathrm{RHS} \f$
 *
 * @param[out] params:               pointer to fully initialized ghl_parameters struct
 * 
 * @returns void
 */
void ghl_initialize_params(
      const ghl_con2prim_method_t main_routine,
      const ghl_con2prim_method_t backup_routine[3],
      const bool evolve_entropy,
      const bool evolve_temp,
      const bool calc_prim_guess,
      const double psi6threshold,
      const double max_Lorentz_factor,
      const double Lorenz_damping_factor,
      ghl_parameters *restrict params) {

  params->main_routine              = main_routine;
  params->backup_routine[0]         = backup_routine[0];
  params->backup_routine[1]         = backup_routine[1];
  params->backup_routine[2]         = backup_routine[2];
  params->evolve_entropy            = evolve_entropy;
  params->evolve_temp               = evolve_temp;
  params->calc_prim_guess           = calc_prim_guess;
  params->psi6threshold             = psi6threshold;
  params->max_Lorentz_factor         = max_Lorentz_factor;
  params->inv_sq_max_Lorentz_factor  = 1.0/SQR(max_Lorentz_factor);
  params->Lorenz_damping_factor     = Lorenz_damping_factor;

  // Initialize default Con2Prim values
  params->con2prim_max_iterations   = 30;
  params->con2prim_solver_tolerance = 1e-10;

  // Initialize default PPM values (from Colella & Woodward)
  params->ppm_flattening_epsilon = 0.33;
  params->ppm_flattening_omega1  = 0.75;
  params->ppm_flattening_omega2  = 10.0;
  params->ppm_shock_k0           = 0.1;
  params->ppm_shock_eta1         = 20.0;
  params->ppm_shock_eta2         = 0.05;
  params->ppm_shock_epsilon      = 0.01;
}
