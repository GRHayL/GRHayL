#ifndef GHL_M1_H
#define GHL_M1_H

#include "ghl.h"
#include <limits.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @addtogroup Radiation
 *  @{ */

typedef struct ghl_m1_rad_state {
  double E;
  double F[3];
} ghl_m1_rad_state;

typedef enum {
  ghl_m1_closure_solve_converged = 0,
  ghl_m1_closure_solve_endpoint_fallback,
  ghl_m1_closure_solve_iteration_exhausted,
  ghl_m1_closure_solve_invalid
} ghl_m1_closure_solve_status_t;

typedef struct ghl_m1_closure {
  double P[3][3];
  /** Eddington factor. */
  double chi;
  /** Comoving reduced-flux magnitude sqrt(H_mu H^mu)/J. */
  double xi;
  /** |J^2 xi^2 - H_mu H^mu| divided by the closure residual scale; zero for
   * the admissibility fallback. */
  double root_residual;
  /** Number of scalar-root iterations; zero for an endpoint or admissibility
   * fallback. */
  int root_iterations;
  ghl_m1_closure_solve_status_t solve_status;
  /**
   * True when the published tensor uses the primary full-four-dimensional
   * construction. False denotes the exceptional finite non-PSD admissibility
   * fallback. This field is diagnostic only and is not a runtime selector.
   */
  bool four_point_compatibility;
} ghl_m1_closure;

typedef struct ghl_m1_closure_counters {
  unsigned long long ordinary_convergence;
  unsigned long long endpoint_fallback;
  unsigned long long iteration_exhaustion;
  unsigned long long invalid_state;
  unsigned long long downstream_repair;
  /** Finite candidates rejected by the normalized-residual acceptance gate. */
  unsigned long long residual_rejection;
} ghl_m1_closure_counters;

typedef enum {
  ghl_m1_closure_failure_none = 0,
  ghl_m1_closure_failure_workspace = 1,
  ghl_m1_closure_failure_comoving_energy = 2,
  ghl_m1_closure_failure_comoving_flux_norm = 3,
  ghl_m1_closure_failure_residual = 4,
  ghl_m1_closure_failure_residual_gate = 5,
  ghl_m1_closure_failure_tensor_validation = 6
} ghl_m1_closure_failure_stage_t;

typedef enum {
  /**
   * Canonical realizability repair: when F^2 > E^2(1-epsilon_c), multiply
   * F_i by E^2(1-epsilon_c)/F^2. This is the only supported value.
   */
  ghl_m1_repair_linear_factor_compatibility = 1
} ghl_m1_repair_policy_t;

typedef enum {
  ghl_m1_repair_branch_none = 0,
  ghl_m1_repair_branch_energy_floor = 1,
  ghl_m1_repair_branch_flux_rescale = 2,
  ghl_m1_repair_branch_energy_and_flux = 3
} ghl_m1_repair_branch_t;

typedef struct ghl_m1_repair_diagnostics {
  /** Always ghl_m1_repair_linear_factor_compatibility. */
  ghl_m1_repair_policy_t policy;
  double original_flux_norm;
  double permitted_flux_norm;
  double scale;
  ghl_m1_rad_state repaired_state;
  ghl_m1_repair_branch_t branch;
} ghl_m1_repair_diagnostics;

typedef struct ghl_m1_closure_decomposition_diagnostic {
  double Pthin_dd[3][3];
  double Pthick_dd[3][3];
  double chi_minerbo_xi;
  double xi_HaHa_over_J2;
  double dthin_scalar;
  double dthick_scalar;
  double Pthick_minus_Pthin_dd[3];
  double Pth_dd_3_3_UU;
  double Pth_dd_0_0_DD;
} ghl_m1_closure_decomposition_diagnostic;

typedef struct ghl_m1_parameters ghl_m1_parameters;

typedef struct ghl_m1_comoving {
  double J;
  double HU[3];
  double HD[3];
  double Hn;
} ghl_m1_comoving;

typedef struct ghl_m1_sources {
  double S_E;
  double S[3];
} ghl_m1_sources;

typedef struct ghl_m1_diagnostics {
  double closure_xi;
  double closure_root_residual;
  int closure_root_iterations;
  ghl_m1_closure_solve_status_t closure_solve_status;
  /** Eulerian squared flux factor. */
  double r;
  /** Minerbo Eddington factor. */
  double chi_eddington;
  bool realizability_repaired;
  double Jthick;
  bool Jthick_is_valid;
  double diffusion_blend_factor[2][3];
} ghl_m1_diagnostics;

typedef struct ghl_m1_implicit_solve_diagnostics {
  int newton_iterations;
  int line_search_backtracks;
  int fallback_substeps;
  bool used_fallback_substepping;
  double residual_max_norm;
  /** Mixed absolute/relative weighted residual merit Phi. Success implies
   * residual_scaled_norm <= 1. */
  double residual_scaled_norm;
  /** Bitwise OR of ghl_m1_solution_path_flag_t values. */
  unsigned int solution_path_flags;
} ghl_m1_implicit_solve_diagnostics;

/**
 * Observable stages of the shared four-variable Newton driver.
 *
 * The driver reports only stages it owns. Residual events include the initial
 * iterate and every line-search trial evaluation, including a re-evaluation
 * after admissibility projection. Domain-specific callbacks remain responsible
 * for exposing trial-state construction, closure evaluation, rate/opacity
 * snapshots, and residual assembly. A solve that begins iteration emits a
 * completed-solve event on both success and failure. On iteration exhaustion,
 * that event carries the current iterate; its residual argument may be NULL
 * for that final event.
 */
typedef enum {
  ghl_m1_solver_stage_residual = 0,
  ghl_m1_solver_stage_jacobian,
  ghl_m1_solver_stage_newton_step,
  ghl_m1_solver_stage_completed_solve
} ghl_m1_solver_stage_t;

typedef struct ghl_m1_newton_diagnostics {
  int iterations;
  int backtracks;
  bool used_projection;
  double residual_max_norm;
  double residual_weighted_merit;
} ghl_m1_newton_diagnostics;

/**
 * Observer used by the public structured Newton diagnostic boundary.
 *
 * The callback is diagnostic-only and must not modify @p U, @p residual, or
 * @p diagnostics. @p residual is non-NULL for residual, Jacobian, and accepted
 * step events. It may be NULL only on the completed-solve event reported after
 * the iteration budget is exhausted. A completed-solve event is emitted for a
 * solve that entered iteration, including terminal callback or solve failure;
 * its @p U argument is the final accepted/current iterate.
 *
 * @param observer_context Opaque caller context supplied in
 *        ghl_m1_newton_callbacks.
 * @param stage Event kind reported by the driver.
 * @param status Callback or solve status associated with the event.
 * @param U Current or completed densitized vector in the ordering
 *        {tilde E, tilde F_x, tilde F_y, tilde F_z}.
 * @param residual Current four-component residual, or NULL only for the
 *        exhaustion completed-solve event.
 * @param diagnostics Snapshot of driver diagnostics for this event; it may be
 *        NULL when the caller did not request output diagnostics.
 */
typedef void (*ghl_m1_solver_stage_observer)(
      void *restrict observer_context,
      const ghl_m1_solver_stage_t stage,
      const ghl_error_codes_t status,
      const double U[4],
      const double residual[4],
      const ghl_m1_newton_diagnostics *restrict diagnostics);

/**
 * Evaluate the four-component residual at one densitized iterate.
 *
 * @param context Opaque caller context passed through the Newton driver.
 * @param U Densitized iterate in the ordering {tilde E, tilde F_x, tilde F_y,
 *        tilde F_z}; the callback must not modify it.
 * @param residual Output residual in the same four-component ordering.
 * @return Status consumed by the driver's admissibility and retry policy;
 *         @c ghl_success publishes a residual, while
 *         @c ghl_error_m1_implicit_admissibility permits the driver's
 *         admissibility projection during a line search.
 */
typedef ghl_error_codes_t (*ghl_m1_newton_residual_callback)(
      const void *restrict context,
      const double U[4],
      double residual[4]);

/**
 * Evaluate the 4x4 Jacobian at one densitized iterate.
 *
 * @param context Opaque caller context passed through the Newton driver.
 * @param U Densitized iterate in the ordering {tilde E, tilde F_x, tilde F_y,
 *        tilde F_z}; the callback must not modify it.
 * @param residual Residual at @p U, supplied by the driver and read-only.
 * @param jacobian Output Jacobian with row/column order matching @p U.
 * @return @c ghl_success when @p jacobian is valid; otherwise the driver
 *         reports the callback error and terminates the solve.
 */
typedef ghl_error_codes_t (*ghl_m1_newton_jacobian_callback)(
      const void *restrict context,
      const double U[4],
      const double residual[4],
      double jacobian[4][4]);

/**
 * Callback bundle for the public structured Newton diagnostic boundary.
 *
 * The @c residual and @c jacobian members are required. The @c observer
 * member is optional; when present, it receives the events described by
 * ghl_m1_solver_stage_observer. The driver does not interpret or own
 * @c observer_context.
 */
typedef struct ghl_m1_newton_callbacks {
  ghl_m1_newton_residual_callback residual;
  ghl_m1_newton_jacobian_callback jacobian;
  ghl_m1_solver_stage_observer observer;
  void *observer_context;
} ghl_m1_newton_callbacks;

typedef enum {
  ghl_m1_solution_path_primary_convergence = 1u << 0,
  ghl_m1_solution_path_closure_fallback = 1u << 1,
  ghl_m1_solution_path_projection = 1u << 2,
  ghl_m1_solution_path_line_search_backtracking = 1u << 3,
  ghl_m1_solution_path_substepping = 1u << 4,
  ghl_m1_solution_path_endpoint_acceptance = 1u << 5,
  ghl_m1_solution_path_terminal_failure = 1u << 6
} ghl_m1_solution_path_flag_t;

/** Explicit local neutrino-source route. The zero value preserves the
 * existing GRHayL implicit source behavior. */
typedef enum {
  ghl_m1_neutrino_source_grhayl_implicit = 0,
  ghl_m1_neutrino_source_branched_compatibility
} ghl_m1_neutrino_source_policy_t;

/** Matter-composition bookkeeping used by the source-policy dispatcher. */
typedef enum {
  /* Preserve GRHayL's established charged-current-only composition contract. */
  ghl_m1_neutrino_ye_from_charged_current = 0,
  /* Use the signed total-number composition source when requested. */
  ghl_m1_neutrino_ye_from_signed_total_number
} ghl_m1_neutrino_ye_policy_t;

/** Selected branch of the transactional local neutrino source update. */
typedef enum {
  ghl_m1_neutrino_source_path_none = 0,
  ghl_m1_neutrino_source_path_thin_explicit,
  ghl_m1_neutrino_source_path_thick_equilibrium,
  ghl_m1_neutrino_source_path_scattering_dominated,
  ghl_m1_neutrino_source_path_general_implicit,
  ghl_m1_neutrino_source_path_terminal_no_update,
  ghl_m1_neutrino_source_path_hard_failure
} ghl_m1_neutrino_source_path_t;

/** Host-supplied controls for the opt-in branched source policy.
 * A nonpositive thick or scattering threshold disables that shortcut.  A
 * negative thermalized-number threshold disables the optional equilibrium
 * mean-energy number update. A nonnegative threshold selects that projection
 * when dt_alpha*kappa_a_N is at least the threshold; zero therefore invokes
 * it even when the opacity or dt is zero, and may change N. The projection
 * uses the repaired endpoint E/F current and is distinct from both backward-
 * Euler number integration and the N_floor repair, which remains a separate
 * step. The zero-valued Y_e policy preserves GRHayL's charged-current-only contract;
 * the signed-total option is available when signed-total composition
 * bookkeeping is required. */
typedef struct {
  ghl_m1_neutrino_source_policy_t policy;
  double thick_equilibrium_threshold;
  double scattering_threshold;
  /**
   * Number-stiffness threshold for equilibrium mean-energy projection.
   * Negative disables the projection; zero selects it even for zero
   * opacity or zero timestep. The value is compared with the dimensionless
   * product dt_alpha*kappa_a_N.
   */
  double thermalized_number_threshold;
  bool allow_closure_fallback;
  bool interaction_sources_already_applied;
  ghl_m1_neutrino_ye_policy_t ye_policy;
} ghl_m1_neutrino_source_options;

typedef struct {
  ghl_m1_neutrino_source_path_t path;
  bool closure_fallback_used;
  bool terminal_no_update;
  ghl_m1_implicit_solve_diagnostics implicit;
} ghl_m1_neutrino_source_diagnostics;

typedef enum {
  ghl_m1_dirn0 = 0,
  ghl_m1_dirn1,
  ghl_m1_dirn2
} ghl_m1_direction_t;

struct ghl_m1_parameters {
  double epsilon_c;
  /* Historical name; stores the admissible squared reduced-flux limit
   * (F/E)^2 <= 1 - epsilon_c. */
  double one_minus_epsilon_c_sq;
  double E_floor;
  /**
   * The only supported value is
   * ghl_m1_repair_linear_factor_compatibility; direct mutation to any other
   * value is rejected by the M1 validation contract.
   */
  ghl_m1_repair_policy_t repair_policy;
  double closure_root_tolerance;
  int closure_root_max_iterations;
  double zeta_min;
  double fd_epsilon_rel;
  double fd_epsilon_abs;
  int newton_max_iterations;
  // Relative Newton tolerance used by the mixed solver.
  double newton_tolerance;
  // Absolute Newton tolerance in undensitized radiation-energy units.
  double newton_absolute_tolerance;
  /** Maximum normalized consistency residual accepted for publication. */
  double closure_root_residual_tolerance;
  /* Four-point transport controls. */
  /**
   * Four-point limiter parameter in [0, 2]. The initializer uses the
   * canonical value minmod_theta = 1.0.
   */
  double minmod_theta;
  /**
   * Lower bound for four-point opacity dissipation suppression in [0, 1].
   * The initializer uses the canonical value mindiss = 0.0.
   */
  double mindiss;
};

/**
 * Initialize the M1 runtime parameters used by the radiation routines.
 *
 * This routine initializes the canonical M1 method: metric light-cone
 * transport speeds, linear-factor realizability repair, and the primary full
 * four-dimensional closure with its flagged physical-PSD admissibility
 * fallback.
 */
ghl_error_codes_t ghl_m1_initialize(
      const double epsilon_c,
      const double E_floor,
      const double zeta_min,
      const double fd_epsilon_rel,
      const double fd_epsilon_abs,
      const int newton_max_iterations,
      const double newton_tolerance,
      ghl_m1_parameters *restrict m1_params);

/** Initialize M1 parameters with explicit relative and absolute Newton
 * tolerances. The absolute tolerance is expressed in undensitized radiation
 * energy-density units and is densitized at the solve point. */
ghl_error_codes_t ghl_m1_initialize_with_newton_tolerances(
      const double epsilon_c,
      const double E_floor,
      const double zeta_min,
      const double fd_epsilon_rel,
      const double fd_epsilon_abs,
      const int newton_max_iterations,
      const double newton_relative_tolerance,
      const double newton_absolute_tolerance,
      ghl_m1_parameters *restrict m1_params);

/** Update the mixed Newton tolerances on an initialized parameter bundle. */
ghl_error_codes_t ghl_m1_set_newton_tolerances(
      const double newton_relative_tolerance,
      const double newton_absolute_tolerance,
      ghl_m1_parameters *restrict m1_params);

/** Set the bracketed scalar Minerbo root-solver controls. */
ghl_error_codes_t ghl_m1_set_closure_solver_controls(
      const double root_interval_tolerance,
      const int root_max_iterations,
      ghl_m1_parameters *restrict m1_params);

/** Set the maximum normalized Minerbo consistency residual that may be
 * published. Finite candidates above this bound fail transactionally. */
ghl_error_codes_t ghl_m1_set_closure_residual_tolerance(
      const double max_normalized_residual,
      ghl_m1_parameters *restrict m1_params);

/** Reset and snapshot process-wide closure outcome counters. */
void ghl_m1_reset_closure_counters(void);
void ghl_m1_get_closure_counters(
      ghl_m1_closure_counters *restrict counters);

/** Return the stage at which the most recent closure call failed. */
void ghl_m1_get_last_closure_failure_stage(
      ghl_m1_closure_failure_stage_t *restrict stage);

void ghl_m1_get_last_closure_validation_reason(int *restrict reason);

/**
 * Runtime contract for public M1 compute kernels:
 *
 * - Initialization/registration routines validate configuration explicitly.
 * - Runtime kernels assume valid, non-NULL inputs in production builds.
 * - Defining GRHAYL_M1_DEBUG enables additional expensive input validation.
 */

/**
 * Apply only the scalar radiation-energy floor, without flux repair.
 * floor_applied may be NULL. Outputs are unchanged on error.
 */
ghl_error_codes_t ghl_m1_apply_energy_floor(
      const ghl_m1_parameters *restrict m1_params,
      const double E_in,
      double *restrict E_out,
      bool *restrict floor_applied);

/** Repair a radiation state using the canonical linear-factor bound. */
ghl_error_codes_t ghl_m1_realizability_repair(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      ghl_m1_rad_state *restrict rad_state);

/** Compute the relativistic Minerbo closure using fluid primitives. */
ghl_error_codes_t ghl_m1_compute_closure_minerbo(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      ghl_m1_closure *restrict closure);

/** Compute the production Minerbo closure. */
ghl_error_codes_t ghl_m1_compute_closure_with_primitives(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      ghl_m1_closure *restrict closure);

/**
 * Compute diagnostic closure-decomposition fields from the public closure path.
 *
 * This routine is diagnostic-only. It reports the production thin and
 * relativistic diffusion-limit tensors together with the solved weights.
 */
ghl_error_codes_t ghl_m1_compute_closure_decomposition_diagnostic(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      ghl_m1_closure_decomposition_diagnostic *restrict diagnostic);

/**
 * Compute comoving-frame radiation moments from the Eulerian state and closure.
 */
ghl_error_codes_t ghl_m1_compute_comoving_moments(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_comoving *restrict comoving);

/**
 * Compute geometric source terms for the M1 evolution equations.
 *
 * The directional metric derivatives are packed in ghl_metric_quantities as
 * lapse -> d_i(alpha), betaU[j] -> d_i(beta^j), and gammaDD[j][k] -> d_i(gamma_jk).
 */
ghl_error_codes_t ghl_m1_compute_geometry_sources(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_metric_quantities *restrict metric_derivs_x,
      const ghl_metric_quantities *restrict metric_derivs_y,
      const ghl_metric_quantities *restrict metric_derivs_z,
      const ghl_extrinsic_curvature *restrict curv,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_sources *restrict geometry_sources);

/**
 * Convert undensitized radiation interaction sources into equal-and-opposite
 * matter conservative RHS contributions.
 *
 * If the radiation equations receive
 *   + alpha * sqrt(gamma) * (S_E, S_i),
 * this routine returns the matching matter contributions
 *   - alpha * sqrt(gamma) * (S_E, S_i)
 * for the matter energy and covariant momentum conservative equations.
 */
ghl_error_codes_t ghl_m1_compute_matter_coupling_sources(
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_sources *restrict interaction_sources,
      double *restrict source_tildetau,
      double source_tildeS[3]);

/**
 * Compute the characteristic wavespeeds used by the M1 transport flux.
 *
 * This routine returns the canonical conservative light-cone speed estimates,
 * -beta^d +/- alpha sqrt(gamma^{dd}). It does not compute closure-dependent
 * M1 eigenvalue estimates or apply an optical-depth cap.
 *
 * metric_face must provide a coherent face metric with SPD and
 * inverse-consistent gammaDD and gammaUU tensors.
 */
ghl_error_codes_t ghl_m1_compute_raw_lightcone_speeds(
      const ghl_metric_quantities *restrict metric_face,
      const ghl_m1_direction_t direction,
      double *restrict s_minus_raw,
      double *restrict s_plus_raw);

ghl_error_codes_t ghl_m1_clip_hll_speeds(
      const double s_minus_raw,
      const double s_plus_raw,
      double *restrict s_minus,
      double *restrict s_plus);

ghl_error_codes_t ghl_m1_compute_wavespeeds(
      const ghl_metric_quantities *restrict metric_face,
      const ghl_m1_direction_t direction,
      double *restrict s_minus,
      double *restrict s_plus);

/** Undensitized pointwise physical E/F_i flux in one coordinate direction. */
ghl_error_codes_t ghl_m1_compute_physical_flux(
      const ghl_metric_quantities *restrict metric_face,
      const ghl_m1_direction_t direction,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      double *restrict flux_E,
      double flux_F[3]);

/** Symmetric Rusanov E/F_i flux from caller-supplied physical fluxes. */
ghl_error_codes_t ghl_m1_compute_rusanov_flux(
      const ghl_m1_rad_state *restrict state_L,
      const ghl_m1_rad_state *restrict state_R,
      const double physical_flux_E_L,
      const double physical_flux_F_L[3],
      const double physical_flux_E_R,
      const double physical_flux_F_R[3],
      const double speed,
      double *restrict flux_E,
      double flux_F[3]);

/** Scalar symmetric Rusanov flux, used by neutrino number transport. */
ghl_error_codes_t ghl_m1_compute_number_rusanov_flux(
      const double N_L,
      const double N_R,
      const double physical_flux_L,
      const double physical_flux_R,
      const double speed,
      double *restrict number_flux);

/** GRHayL ordering for the five grey-neutrino transport components. */
typedef enum {
  ghl_m1_neutrino_transport_N = 0,
  ghl_m1_neutrino_transport_E,
  ghl_m1_neutrino_transport_Fx,
  ghl_m1_neutrino_transport_Fy,
  ghl_m1_neutrino_transport_Fz,
  ghl_m1_neutrino_transport_component_count
} ghl_m1_neutrino_transport_component_t;

typedef struct {
  double phi[ghl_m1_neutrino_transport_component_count];
  bool sawtooth[ghl_m1_neutrino_transport_component_count];
  double opacity_suppression;
  double face_speed;
} ghl_m1_four_point_transport_diagnostics;

/** Evaluate one canonical four-point blended face flux.
 *
 * State and physical-flux inputs are undensitized and use
 * {N,E,Fx,Fy,Fz}. The stencil is {j-1,j,j+1,j+2}; physical_flux_L/R and
 * speed_L/R belong to j and j+1. delta_x is coordinate spacing, not proper
 * face-normal thickness. The returned flux is multiplied by the face
 * sqrt(det(gamma)) exactly once. The canonical caller supplies uncapped
 * light-cone speeds and disables the separate diffusion correction; requests
 * for either unsupported transport policy are rejected.
 *
 * @param m1_params Initialized M1 limiter parameters. The canonical
 *        four-point controls are read from this bundle.
 * @param metric_face Coherent face metric, including a positive
 *        @c sqrt_detgamma used for the single output densitization.
 * @param state_stencil Four states in stencil order {j-1,j,j+1,j+2}; each
 *        row uses the component order {N,E,Fx,Fy,Fz}.
 * @param physical_flux_L Undensitized physical flux at cell j, in the same
 *        five-component order.
 * @param physical_flux_R Undensitized physical flux at cell j+1, in the same
 *        five-component order.
 * @param speed_L Uncapped light-cone speed associated with cell j.
 * @param speed_R Uncapped light-cone speed associated with cell j+1.
 * @param kappa_face Face transport opacity in code inverse-length units,
 *        consistent with @p delta_x.
 * @param delta_x Coordinate spacing in the selected direction. This is not
 *        proper normal spacing.
 * @param diffusion_correction_enabled Must be false for the canonical
 *        operation; true is rejected.
 * @param flux_tilde Output densitized face flux in {N,E,Fx,Fy,Fz} order.
 * @param diagnostics Optional output diagnostics. Pass NULL when they are not
 *        needed.
 * @return @c ghl_success on publication; otherwise the output flux and
 *         diagnostics remain unchanged.
 */
ghl_error_codes_t ghl_m1_compute_neutrino_four_point_transport_flux(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric_face,
      const double state_stencil[4][ghl_m1_neutrino_transport_component_count],
      const double physical_flux_L[ghl_m1_neutrino_transport_component_count],
      const double physical_flux_R[ghl_m1_neutrino_transport_component_count],
      const double speed_L,
      const double speed_R,
      const double kappa_face,
      const double delta_x,
      const bool diffusion_correction_enabled,
      double flux_tilde[ghl_m1_neutrino_transport_component_count],
      ghl_m1_four_point_transport_diagnostics *restrict diagnostics);

/** Componentwise stages exposed for focused diagnostics and testing. */
ghl_error_codes_t ghl_m1_compute_four_point_flux_limiter(
      const ghl_m1_parameters *restrict m1_params,
      const double dum,
      const double duc,
      const double dup,
      double *restrict phi,
      bool *restrict sawtooth);

ghl_error_codes_t ghl_m1_compute_four_point_opacity_suppression(
      const ghl_m1_parameters *restrict m1_params,
      const double kappa_face,
      const double delta_x,
      double *restrict A);

ghl_error_codes_t ghl_m1_compute_four_point_blended_flux(
      const double flux_high,
      const double flux_low,
      const double phi,
      const bool sawtooth,
      const double A,
      double *restrict flux_num);

/**
 * Compute the proper face-normal thickness for a coordinate-direction face.
 *
 * For a face normal to x^d, this returns
 * delta_l = delta_x^d / sqrt(gamma^{dd}) using the supplied face metric.
 */
ghl_error_codes_t ghl_m1_compute_face_normal_delta_l(
      const ghl_metric_quantities *restrict metric_face,
      const ghl_m1_direction_t direction,
      const double delta_x_d,
      double *restrict delta_l);

/** Compute the thick-limit comoving radiation energy and validity flag. */
ghl_error_codes_t ghl_m1_compute_Jthick(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      double *restrict Jthick,
      bool *restrict Jthick_is_valid);

/** Compute the harmonic face diffusion coefficient from transport opacities. */
ghl_error_codes_t ghl_m1_compute_harmonic_diffusion_coefficient(
      const double chi_tr_L,
      const double chi_tr_R,
      double *restrict D_face);

/** Apply the optional Fick diffusion correction to an HLL energy flux.
 * The caller selects this route explicitly; the ordinary public HLL/Rusanov
 * fluxes remain unchanged. */
ghl_error_codes_t ghl_m1_compute_diffusion_flux(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric_face,
      const ghl_m1_direction_t direction,
      const double hll_flux_tildeE,
      const double E_star,
      const double Jthick_L,
      const bool Jthick_L_is_valid,
      const double Jthick_R,
      const bool Jthick_R_is_valid,
      const double gradJthick_face[3],
      const double W_face,
      const double V_face[3],
      const double chi_tr_face,
      const double D_face,
      const double delta_l,
      double *restrict corrected_flux_tildeE,
      double *restrict a_face);

/**
 * Run one host-independent four-variable Newton solve through the shared
 * callback boundary.
 *
 * U_base and U_out use the densitized ordering (tilde E, tilde F_i). The
 * optional observer receives every residual evaluation (initial and trial),
 * Jacobian, accepted-step, and completed-solve events. A residual pointer is
 * non-NULL for all normal events; it may be NULL only for terminal failure
 * after the final iteration. The observer is diagnostic-only: it cannot alter
 * a trial state, status, fallback policy, or production result.
 */
ghl_error_codes_t ghl_m1_newton_solve_4d(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_newton_callbacks *restrict callbacks,
      const void *restrict context,
      const double U_base[4],
      double U_out[4],
      ghl_m1_newton_diagnostics *restrict diagnostics);

/** Evaluate the shared mixed absolute/relative Newton merit. */
double ghl_m1_newton_weighted_merit(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const double U[4],
      const double U_base[4],
      const double residual[4]);

/** Project a densitized E/F_i vector into the configured admissible domain. */
ghl_error_codes_t ghl_m1_newton_project_admissible(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      double U[4]);

/**
 * Build the radiation stress-energy tensor from a state and its closure.
 *
 * This E/F-only kernel validates non-NULL inputs, the metric, finite E/F_i,
 * m1_params->E_floor, the shared M1 realizability cone, and the locally
 * observable supplied-closure tensor properties. It does not inspect a
 * neutrino number density or N_floor. On every error Rmunu is unchanged.
 */
ghl_error_codes_t ghl_m1_compute_stress_energy(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      ghl_stress_energy *restrict Rmunu);

/**
 * Compute diagnostic quantities associated with an M1 state and closure.
 *
 * The diagnostic r is reconstructed from the supplied state's scaled
 * covector norm, r = (sqrt(gamma^ij F_i F_j) / E)^2, then clamped to
 * [0, 1 - epsilon_c] before chi_eddington is reported from the supplied
 * closure.
 */
ghl_error_codes_t ghl_m1_compute_diagnostics(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_diagnostics *restrict diagnostics);

/* ============================================================================
 * Grey three-species neutrino M1 transport (one-group contract)
 *
 * Scope:
 *   - One transported number density N, one energy density E, one covariant
 *     flux F_i per species.
 *   - Three species: nu_e (lepton weight +1), nu_bar_e / anti-nu_e (-1),
 *     nu_x (lumped heavy lepton, 0).
 *   - Grey, one-group only; no group-indexed transport API is provided.
 *   - M1 closure and comoving moments for E/F_i are reused to construct the
 *     grey number current and its E/F-derived transport velocity.
 *   - Source solve uses frozen, provider-supplied neutrino rates
 *     (ghl_m1_neutrino_rates). The provider is responsible for EOS lookup, weak-equilibrium
 *     targets, unit conversion, and channel-specific microphysics.
 * ========================================================================== */

/**
 * Neutrino species enumeration for the grey three-species M1 framework.
 *
 * The naming follows the implementation whitepaper: nu_bar_e is the code-facing
 * name for the electron antineutrino (conceptual paper: anti-nu_e).
 */
typedef enum {
  ghl_m1_neutrino_nue  = 0,
  ghl_m1_neutrino_anue = 1,
  ghl_m1_neutrino_nux  = 2,
  ghl_m1_neutrino_species_count = 3
} ghl_m1_neutrino_species_t;

/** Grey pair-production channels carried separately from one-body rates. */
typedef enum {
  ghl_m1_neutrino_pair_process_pair = 0,
  ghl_m1_neutrino_pair_process_plasmon,
  ghl_m1_neutrino_pair_process_bremsstrahlung,
  ghl_m1_neutrino_pair_process_count
} ghl_m1_neutrino_pair_process_t;

/**
 * Per-species neutrino state carried through transport and source solves.
 *
 * N is the Eulerian neutrino number density, E is the Eulerian radiation
 * energy density, and F is the covariant Eulerian radiation flux. Densitized
 * host storage (sqrt(gamma) * N, sqrt(gamma) * E, sqrt(gamma) * F_i) must be
 * undensitized by the host with the appropriate cell or face sqrt_detgamma
 * before any Radiation pointwise call.
 */
typedef struct {
  double N;
  double E;
  double F[3];
} ghl_m1_neutrino_state;

typedef enum {
  ghl_m1_neutrino_terminal_fallback_no_update_all = 0
} ghl_m1_neutrino_terminal_fallback_policy_t;

/**
 * Per-species neutrino runtime parameters.
 *
 * Zero initialization permits N == 0, selects strict J > 0, uses
 * 64*DBL_EPSILON for the Gamma_N floor, disables optional comoving-mean-energy
 * bounds, and selects the transactional no-update terminal policy. When
 * enforce_mean_energy_bounds is nonzero, each positive bound is checked
 * independently against the final repaired endpoint mean energy
 * J*Gamma_N/N in the same code energy-per-number units as
 * ghl_m1_neutrino_rates::mean_energy. A nonpositive lower or upper bound
 * disables that bound. Out-of-bounds endpoints are rejected with no clamp;
 * N == 0 retains the existing skip of this ratio check. The check applies to
 * the ordinary implicit update, the thin/thick/scattering source paths, and
 * the public thin-update wrappers.
 */
typedef struct {
  /** Nonnegative neutrino-number floor used by state repair. */
  double N_floor;
  /** Positive lower mean-energy bound; nonpositive disables this bound. */
  double mean_energy_min;
  /** Positive upper mean-energy bound; nonpositive disables this bound. */
  double mean_energy_max;
  /** Nonzero enables validation of the positive mean-energy bounds. */
  int enforce_mean_energy_bounds;
  int terminal_fallback_policy;
  double J_floor;
  double Gamma_N_floor;
} ghl_m1_neutrino_parameters;

/**
 * Provider-supplied frozen neutrino rate bundle for a single cell and species.
 *
 * The provider computes this bundle from EOS state, composition, blocking,
 * degeneracy, weak equilibrium, and channel-specific microphysics. Radiation
 * consumes the bundle; it does not perform EOS lookup or evaluate production
 * weak-rate formulas. Field semantics:
 *   - eta_N, eta_E     : one-body grey number and energy emissivities (>= 0)
 *   - kappa_a_N,
 *     kappa_a_E        : one-body grey number and energy absorption opacities
 *                        (>= 0); for electron flavors these are the
 *                        independent charged-current/scattering channels.
 *   - kappa_s          : isoenergetic scattering opacity (>= 0)
 *   - kappa_tr         : kappa_a_E + kappa_s
 *   - n_eq, J_eq       : grey comoving-frame equilibrium number and energy
 *                        density targets (>= 0); the provider supplies
 *                        neutrino weak-equilibrium targets.
 *   - mean_energy      : positive grey mean neutrino energy used for the
 *                        reduced number-flux closure and as a diagnostic.
 *   - lepton_weight    : species lepton weight (+1, -1, or 0).
 *   - eta_N_cc, kappa_a_N_cc: charged-current subset used exclusively for
 *                        electron-lepton exchange; both are zero for nu_x.
 *   - eta_N_pair, eta_E_pair: pair, plasmon, and bremsstrahlung emissivities,
 *                        indexed by ghl_m1_neutrino_pair_process_t; number
 *                        emissivities are shared by the electron pair and
 *                        energy emissivities retain each raw spectrum.
 * For nu_x, the scalar coefficients retain the aggregate approximation and
 * the pair arrays are zero; its four-flavor multiplicity is already applied.
 * Initialize the complete struct, including absent process arrays. Electron
 * pair channels require ghl_m1_solve_neutrino_pair_source_update; single-species
 * source operations reject them and legacy non-CC electron number rates.
 * Electron transport opacity includes the partner-dependent inverse pair
 * energy opacity in addition to scalar kappa_tr (Radiation/PAIR_SOURCE_MODEL.md).
 * Validation enforces aggregate and charged-current Kirchhoff identities,
 * J_eq = n_eq*mean_energy, and exact species/lepton-weight consistency.
 */
typedef struct {
  ghl_m1_neutrino_species_t species;
  double eta_N;
  double eta_E;
  double kappa_a_N;
  double kappa_a_E;
  double kappa_s;
  double kappa_tr;
  double n_eq;
  double J_eq;
  double mean_energy;
  double lepton_weight;
  double eta_N_cc;
  double kappa_a_N_cc;
  /** Number emissivity for each shared electron-flavor pair process. */
  double eta_N_pair[ghl_m1_neutrino_pair_process_count];
  /** Raw energy emissivity for each electron-flavor pair process. */
  double eta_E_pair[ghl_m1_neutrino_pair_process_count];
} ghl_m1_neutrino_rates;

/**
 * Coupled radiation-matter exchange increments reported by the local
 * homogeneous implicit neutrino update.
 *
 * dN_rad_total, dE_rad, and dF_rad are the radiation increments. dL_rad_cc is
 * the charged-current radiation lepton-number increment retained for
 * diagnostics. The source-policy dispatcher normally uses it to compute
 * dYe_matter; its optional signed-total mode instead uses the signed total
 * number increment and species lepton weight. dTau_matter and dS_matter are
 * the equal-and-opposite densitized matter conservative
 * increments (radiation gets + alpha*sqrt(gamma)*(S_E,S_i); matter gets the
 * negation). dYe_matter is the host-applied Delta Y_e recommendation
 * computed via ghl_m1_compute_neutrino_lepton_increment. The host applies
 * one coupled limiter scalar theta to all six increments together to
 * preserve energy, momentum, and electron-lepton-number conservation.
 */
typedef struct {
  double dN_rad_total;
  double dL_rad_cc;
  double dE_rad;
  double dF_rad[3];
  double dTau_matter;
  double dS_matter[3];
  double dYe_matter;
} ghl_m1_neutrino_exchange;

/**
 * Local neutrino diagnostics counters.
 *
 * The comoving mean-energy fields report, per solve, whether
 * J_out * Gamma_N,out / N_out is consistent with the
 * provider-supplied rates->mean_energy and the equilibrium ratio
 * J_eq / n_eq, and whether the mean-energy diagnostic is invalid
 * because N_out is at or below the N floor. These diagnostics are
 * post-hoc and do not change the update; they are populated by
 * ghl_m1_solve_neutrino_implicit_homogeneous_update in both the
 * success and the terminal-fallback paths.
 *
 * The repair-budget fields accumulate, per species and over a
 * run, the absolute, component-wise magnitudes of state mutations introduced
 * by the N-floor and E/F_i realizability repairs. The historical repair_d*
 * field names are retained for source compatibility; these fields are not
 * signed conservation deltas and opposite-signed mutations do not cancel.
 * The host sets repair_stage and repair_lepton_weight as transient inputs to
 * each call of ghl_m1_repair_neutrino_state; repair_dL_e accumulates the
 * absolute magnitude of the corresponding weighted number mutation. The
 * Conservation accounting uses these budgets to ensure that repair leakage is
 * included in the total accounting.
 */
typedef struct {
  int provider_validation_failures;
  int source_converged;
  int source_terminal_fallbacks;
  int source_failures;
  int N_floor_repairs;
  int EF_repairs;
  int limiter_reductions;
  /* Comoving mean-energy diagnostics (post-update). */
  double mean_energy_diag;
  int    mean_energy_consistent;
  int    Jeq_over_neq_consistent;
  int    mean_energy_diag_invalid;
  /* Repair-budget accounting (accumulated over a run). */
  double repair_dN;
  double repair_dE;
  double repair_dF[3];
  double repair_dL_e;
  /* Transient inputs set by the host before each repair call.
   * repair_stage: 0=unknown, 1=post-transport, 2=post-source,
   *               3=face-state. repair_lepton_weight: the species
   *               lepton weight used to compute the repair_dL_e magnitude. */
  int    repair_stage;
  double repair_lepton_weight;
  /* Accepted positive rate-product identities that rounded to exact zero. */
  int    rate_product_underflows;
} ghl_m1_neutrino_diagnostics;

/**
 * Initialize neutrino diagnostics before their first use.
 *
 * Clears all counters, accumulated repair magnitudes, and transient repair
 * metadata. Passing NULL is a no-op. Value-initializing the struct to zero is
 * equivalent, but this function gives callers a stable public API for that
 * requirement.
 */
void ghl_m1_neutrino_diagnostics_initialize(
      ghl_m1_neutrino_diagnostics *restrict diagnostics);

/**
 * Apply only the scalar neutrino-number floor, without E/F_i repair.
 *
 * @param nu_params Neutrino parameters containing the nonnegative
 *        @c N_floor.
 * @param N_in Undensitized input number density.
 * @param N_out Output number density; unchanged on error.
 * @param floor_applied Optional flag set true when the floor changes the
 *        value; NULL is allowed.
 * @return @c ghl_success or an M1 validation error. A NULL required pointer
 *         returns @c ghl_error_m1_null_pointer.
 */
ghl_error_codes_t ghl_m1_apply_neutrino_number_floor(
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const double N_in,
      double *restrict N_out,
      bool *restrict floor_applied);

/**
 * Validate a provider-supplied ghl_m1_neutrino_rates bundle.
 *
 * Checks species range, lepton weight matching the species, finiteness and
 * nonnegativity of all emissivities, opacities, equilibrium targets, mean
 * energy > 0, and kappa_tr >= 0. Increments
 * diagnostics->provider_validation_failures when diagnostics is non-NULL and
 * a validation rule fails. Accepted positive products that round to exact
 * zero increment diagnostics->rate_product_underflows. The bundle is
 * intentionally not modified.
 *
 * @param rates Provider-supplied rate bundle for one species.
 * @param diagnostics Optional diagnostics record. A validation failure is
 *        counted when this pointer is non-NULL.
 * @return @c ghl_success when the bundle satisfies the aggregate and
 *         charged-current identities; otherwise
 *         @c ghl_error_m1_microphysics_failure or
 *         @c ghl_error_m1_null_pointer.
 */
ghl_error_codes_t ghl_m1_validate_neutrino_rates(
      const ghl_m1_neutrino_rates *restrict rates,
      ghl_m1_neutrino_diagnostics *restrict diagnostics);

/**
 * Repair a single neutrino state so it satisfies the per-species realizability
 * constraints.
 *
 * Applies the N floor (max(N, N_floor)) and delegates E/F_i repair to
 * ghl_m1_realizability_repair. Increments diagnostics counters for each
 * mutation when diagnostics is non-NULL. The repair_dN, repair_dE,
 * repair_dF[i], and repair_dL_e fields accumulate absolute component-wise
 * mutation magnitudes without signed cancellation. Nonfinite input state is
 * rejected.
 *
 * @param m1_params Shared M1 repair parameters.
 * @param nu_params Neutrino number-floor and mean-energy parameters.
 * @param metric Cell metric used by the E/F realizability repair.
 * @param state In/out undensitized state; it is modified only after the full
 *        repair succeeds.
 * @param diagnostics Optional diagnostics accumulator for repair accounting;
 *        NULL suppresses counter updates.
 * @return @c ghl_success on publication; otherwise the state remains
 *         unchanged and the returned error identifies invalid input or
 *         geometry.
 */
ghl_error_codes_t ghl_m1_repair_neutrino_state(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      ghl_m1_neutrino_state *restrict state,
      ghl_m1_neutrino_diagnostics *restrict diagnostics);

static inline ghl_m1_rad_state ghl_m1_neutrino_project_rad_state(
      const ghl_m1_neutrino_state *restrict state) {
  ghl_m1_rad_state rad_state;
  rad_state.E = state->E;
  rad_state.F[0] = state->F[0];
  rad_state.F[1] = state->F[1];
  rad_state.F[2] = state->F[2];
  return rad_state;
}

static inline ghl_error_codes_t ghl_m1_compute_neutrino_Jthick(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      double *restrict Jthick,
      bool *restrict Jthick_is_valid) {
  if(state == NULL)
    return ghl_error_m1_null_pointer;
  const ghl_m1_rad_state rad_state = ghl_m1_neutrino_project_rad_state(state);
  return ghl_m1_compute_Jthick(
      m1_params, metric, prims, &rad_state, Jthick, Jthick_is_valid);
}

/** Neutrino-state wrapper for the optional energy-flux diffusion correction. */
ghl_error_codes_t ghl_m1_compute_neutrino_diffusion_flux(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric_face,
      const ghl_m1_direction_t direction,
      const double hll_flux_tildeE,
      const double E_star,
      const double Jthick_L,
      const bool Jthick_L_is_valid,
      const double Jthick_R,
      const bool Jthick_R_is_valid,
      const double gradJthick_face[3],
      const double W_face,
      const double V_face[3],
      const ghl_m1_neutrino_rates *restrict rates_face,
      const double D_face,
      const double delta_l,
      double *restrict corrected_flux_tildeE,
      double *restrict a_face);

static inline ghl_error_codes_t ghl_m1_compute_neutrino_closure(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      ghl_m1_closure *restrict closure) {
  if(state == NULL)
    return ghl_error_m1_null_pointer;
  const ghl_m1_rad_state rad_state = ghl_m1_neutrino_project_rad_state(state);
  return ghl_m1_compute_closure_with_primitives(
      m1_params, metric, prims, &rad_state, closure);
}

static inline ghl_error_codes_t ghl_m1_compute_neutrino_comoving_moments(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_comoving *restrict comoving) {
  if(state == NULL)
    return ghl_error_m1_null_pointer;
  const ghl_m1_rad_state rad_state = ghl_m1_neutrino_project_rad_state(state);
  return ghl_m1_compute_comoving_moments(m1_params, metric, prims, &rad_state, closure, comoving);
}

static inline ghl_error_codes_t ghl_m1_compute_neutrino_wavespeeds(
      const ghl_metric_quantities *restrict metric_face,
      const ghl_m1_direction_t direction,
      double *restrict s_minus,
      double *restrict s_plus) {
  return ghl_m1_compute_wavespeeds(
      metric_face, direction, s_minus, s_plus);
}

static inline ghl_error_codes_t ghl_m1_compute_neutrino_stress_energy(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure,
      ghl_stress_energy *restrict Rmunu) {
  if(state == NULL)
    return ghl_error_m1_null_pointer;
  const ghl_m1_rad_state rad_state = ghl_m1_neutrino_project_rad_state(state);
  return ghl_m1_compute_stress_energy(m1_params, metric, &rad_state, closure, Rmunu);
}

static inline ghl_error_codes_t ghl_m1_compute_neutrino_geometry_sources(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_metric_quantities *restrict metric_derivs_x,
      const ghl_metric_quantities *restrict metric_derivs_y,
      const ghl_metric_quantities *restrict metric_derivs_z,
      const ghl_extrinsic_curvature *restrict curv,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_sources *restrict geometry_sources) {
  if(state == NULL)
    return ghl_error_m1_null_pointer;
  const ghl_m1_rad_state rad_state = ghl_m1_neutrino_project_rad_state(state);
  return ghl_m1_compute_geometry_sources(
      m1_params, metric, metric_derivs_x, metric_derivs_y, metric_derivs_z,
      curv, &rad_state, closure, geometry_sources);
}

static inline ghl_error_codes_t ghl_m1_compute_neutrino_matter_coupling_sources(
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_sources *restrict EF_sources,
      double *restrict dTau_matter,
      double dS_matter[3]) {
  return ghl_m1_compute_matter_coupling_sources(metric, EF_sources, dTau_matter, dS_matter);
}

/**
 * Compute the E/F_i and N interaction sources from frozen neutrino rates.
 *
 * EF_sources uses provider-supplied eta_E, kappa_a_E, kappa_s, and a fresh shared
 * comoving moment/current evaluation. The number source is
 * N_source = eta_N - kappa_a_N*N/Gamma_N. The closure is computed freshly.
 * E/F sources and N_source are published
 * transactionally.
 *
 * This is a manual-only explicit interaction operator for hosts that assemble
 * and limit their own source packet. The normal three-species IMEX stage must
 * instead use ghl_m1_solve_neutrino_implicit_homogeneous_update and must not
 * apply these explicit interaction sources in that same stage. GRHayL has no
 * global stage tracker because it cannot infer a host MoL schedule.
 *
 * @param m1_params Initialized shared M1 parameters.
 * @param nu_params Neutrino number-current parameters.
 * @param metric Frozen cell metric.
 * @param prims Frozen fluid primitives used for closure/current evaluation.
 * @param state Undensitized state in {N,E,Fx,Fy,Fz} order.
 * @param rates Validated frozen rates for one species.
 * @param EF_sources Output undensitized E/F interaction source components.
 * @param N_source Output undensitized number source.
 * @return @c ghl_success on publication; otherwise both output sources remain
 *         unchanged.
 */
ghl_error_codes_t ghl_m1_compute_neutrino_interaction_sources(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_neutrino_rates *restrict rates,
      ghl_m1_sources *restrict EF_sources,
      double *restrict N_source);

static inline ghl_error_codes_t ghl_m1_compute_neutrino_explicit_rhs_sources(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_metric_quantities *restrict metric_derivs_x,
      const ghl_metric_quantities *restrict metric_derivs_y,
      const ghl_metric_quantities *restrict metric_derivs_z,
      const ghl_extrinsic_curvature *restrict curv,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure,
      const bool include_interaction_sources,
      const ghl_m1_neutrino_rates *restrict rates,
      double *restrict source_tildeE,
      double source_tildeF[3],
      double *restrict source_tildeN) {
  if(state == NULL || source_tildeE == NULL || source_tildeF == NULL)
    return ghl_error_m1_null_pointer;
  if(include_interaction_sources && (prims == NULL || rates == NULL))
    return ghl_error_m1_null_pointer;

  ghl_m1_sources geometry_sources = {0};
  ghl_error_codes_t error = ghl_m1_compute_neutrino_geometry_sources(
      m1_params, metric, metric_derivs_x, metric_derivs_y, metric_derivs_z,
      curv, state, closure, &geometry_sources);
  if(error != ghl_success)
    return error;

  double candidate_E = geometry_sources.S_E;
  double candidate_F[3] = {
      geometry_sources.S[0], geometry_sources.S[1], geometry_sources.S[2]};
  double candidate_N = 0.0;

  if(include_interaction_sources) {
    ghl_m1_sources interaction_sources = {0};
    double N_source = 0.0;
    error = ghl_m1_compute_neutrino_interaction_sources(
        m1_params, nu_params, metric, prims, state, rates,
        &interaction_sources, &N_source);
    if(error != ghl_success)
      return error;
    const double alpha_sqrt_detgamma = metric->lapse * metric->sqrt_detgamma;
    if(!isfinite(alpha_sqrt_detgamma) || alpha_sqrt_detgamma <= 0.0)
      return ghl_error_m1_invalid_metric;
    candidate_E += alpha_sqrt_detgamma * interaction_sources.S_E;
    for(int i = 0; i < 3; i++)
      candidate_F[i] += alpha_sqrt_detgamma * interaction_sources.S[i];
    candidate_N = alpha_sqrt_detgamma * N_source;
  }

  if(!isfinite(candidate_E) || !isfinite(candidate_N))
    return ghl_error_m1_invalid_state;
  for(int i = 0; i < 3; i++) {
    if(!isfinite(candidate_F[i]))
      return ghl_error_m1_invalid_state;
  }

  *source_tildeE = candidate_E;
  for(int i = 0; i < 3; i++)
    source_tildeF[i] = candidate_F[i];
  if(source_tildeN != NULL)
    *source_tildeN = candidate_N;
  return ghl_success;
}

static inline ghl_error_codes_t ghl_m1_compute_neutrino_diagnostics(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_diagnostics *restrict diagnostics) {
  if(state == NULL)
    return ghl_error_m1_null_pointer;
  const ghl_m1_rad_state rad_state = ghl_m1_neutrino_project_rad_state(state);
  return ghl_m1_compute_diagnostics(m1_params, metric, &rad_state, closure, diagnostics);
}

/**
 * Compute the spatial contravariant number flux for one neutrino species.
 *
 * A fresh shared closure/current evaluation forms
 * Gamma_N = W - Hn/J, n_com = N/Gamma_N, and
 * number_flux^i = n_com*(W V^i + H_perp^i/J). The returned transport velocity
 * is number_flux^i/N, evaluated algebraically so it remains finite, causal,
 * and may be nonzero when N == 0. The number flux is exactly zero for N == 0.
 * If N == 0 and Gamma_N is singular, the number current is still returned as
 * zero with zero transport velocity; Gamma_N is undefined when no particles
 * are present. Nonzero-N states retain strict Gamma_N validation.
 * Outputs are meaningful only on success.
 *
 * Rejects invalid metric, invalid Lorentz factor, N < N_floor, E < E_floor,
 * non-realizable E/F_i, or nonfinite output.
 *
 * @param m1_params Initialized shared M1 parameters.
 * @param nu_params Neutrino number-current parameters, including N_floor and
 *        Gamma_N_floor.
 * @param metric Cell metric used by the current construction.
 * @param prims Fluid primitives used by the closure and comoving transform.
 * @param state Undensitized state in {N,E,Fx,Fy,Fz} order.
 * @param number_flux Output contravariant spatial number flux [3].
 * @param number_transport_velocity Output contravariant transport velocity
 *        [3]. Both output arrays are unchanged on error.
 * @return @c ghl_success on publication; otherwise an M1 validation error.
 */
ghl_error_codes_t ghl_m1_compute_neutrino_number_flux(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      double number_flux[3],
      double number_transport_velocity[3]);

/**
 * Compute the spatial number flux from a caller-supplied closure.
 *
 * This is the closure-reuse companion to
 * ghl_m1_compute_neutrino_number_flux. The closure must have been computed
 * from the supplied metric, primitives, and state. This routine never invokes
 * the closure solver. Outputs are meaningful only on success.
 *
 * @param m1_params Initialized shared M1 parameters.
 * @param nu_params Neutrino number-current parameters.
 * @param metric Cell metric used by the current construction.
 * @param prims Fluid primitives corresponding to @p closure and @p state.
 * @param state Undensitized state in {N,E,Fx,Fy,Fz} order.
 * @param closure Caller-supplied closure for the same metric, primitives, and
 *        E/F state; it is not recomputed.
 * @param number_flux Output contravariant spatial number flux [3].
 * @param number_transport_velocity Output contravariant transport velocity
 *        [3].
 * @return @c ghl_success on publication; otherwise outputs remain unchanged.
 */
ghl_error_codes_t ghl_m1_compute_neutrino_number_flux_from_closure(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure,
      double number_flux[3],
      double number_transport_velocity[3]);

/**
 * Compute the undensitized coordinate physical flux of N in one direction.
 *
 * This is the pointwise transport stage
 *   F_N^d = alpha n^d - beta^d N,
 * where n^i is obtained from the selected primitive-aware closure and number
 * current. It does not densitize the result or apply a numerical flux.
 */
ghl_error_codes_t ghl_m1_compute_neutrino_physical_number_flux(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_direction_t direction,
      double *restrict physical_number_flux);

/** Closure-reuse companion to
 * ghl_m1_compute_neutrino_physical_number_flux. */
ghl_error_codes_t ghl_m1_compute_neutrino_physical_number_flux_from_closure(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure,
      const ghl_m1_direction_t direction,
      double *restrict physical_number_flux);

/**
 * Compute an opt-in combined scalar N and E/F_i symmetric Rusanov flux.
 *
 * The caller supplies one nonnegative, undensitized speed that is applied to
 * all five components in the order {N, E, Fx, Fy, Fz}. The returned fluxes are
 * densitized with metric_face->sqrt_detgamma. No HLL envelope, speed cap, or
 * diffusion intermediate is constructed by this operation.
 *
 * @param m1_params Initialized shared M1 parameters.
 * @param nu_params Neutrino number-current parameters used to validate the
 *        number-current inputs.
 * @param metric_face Coherent face metric; its positive sqrt_detgamma is used
 *        for exactly one output densitization.
 * @param direction Coordinate direction of the E/F physical fluxes.
 * @param state_L Undensitized left state in {N,E,Fx,Fy,Fz} order.
 * @param state_R Undensitized right state in {N,E,Fx,Fy,Fz} order.
 * @param closure_L Closure corresponding to @p state_L.
 * @param closure_R Closure corresponding to @p state_R.
 * @param number_flux_L Undensitized contravariant number flux [3] at the left
 *        state.
 * @param number_flux_R Undensitized contravariant number flux [3] at the right
 *        state.
 * @param number_transport_velocity_L Contravariant number transport velocity
 *        [3] at the left state.
 * @param number_transport_velocity_R Contravariant number transport velocity
 *        [3] at the right state.
 * @param speed Nonnegative common Rusanov speed.
 * @param flux_tildeN Output densitized number flux.
 * @param flux_tildeE Output densitized energy flux.
 * @param flux_tildeF Output densitized covariant momentum flux [3].
 * @return @c ghl_success on publication; otherwise all output fluxes remain
 *         unchanged.
 */
ghl_error_codes_t ghl_m1_compute_neutrino_rusanov_flux(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric_face,
      const ghl_m1_direction_t direction,
      const ghl_m1_neutrino_state *restrict state_L,
      const ghl_m1_neutrino_state *restrict state_R,
      const ghl_m1_closure *restrict closure_L,
      const ghl_m1_closure *restrict closure_R,
      const double number_flux_L[3],
      const double number_flux_R[3],
      const double number_transport_velocity_L[3],
      const double number_transport_velocity_R[3],
      const double speed,
      double *restrict flux_tildeN,
      double *restrict flux_tildeE,
      double flux_tildeF[3]);

/**
 * Compute frozen-rate interaction sources from a caller-supplied closure.
 *
 * The closure must correspond to the supplied metric, primitives, and state.
 * This routine consumes the closure and never invokes the closure solver.
 * Outputs are published
 * transactionally.
 *
 * @param m1_params Initialized shared M1 parameters.
 * @param nu_params Neutrino number-current parameters.
 * @param metric Cell metric used for the comoving moments.
 * @param prims Frozen fluid primitives corresponding to @p state.
 * @param state Undensitized state in {N,E,Fx,Fy,Fz} order.
 * @param closure Caller-supplied closure for the same metric, primitives, and
 *        E/F state; it is not recomputed.
 * @param rates Validated frozen rates for one species.
 * @param EF_sources Output undensitized E/F interaction sources.
 * @param N_source Output undensitized number interaction source.
 * @return @c ghl_success on publication; otherwise both output sources remain
 *         unchanged.
 */
ghl_error_codes_t ghl_m1_compute_neutrino_interaction_sources_from_closure(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure,
      const ghl_m1_neutrino_rates *restrict rates,
      ghl_m1_sources *restrict EF_sources,
      double *restrict N_source);

/**
 * Closed-form backward-Euler update of the neutrino number density.
 *
 * With an endpoint current normalization Gamma_N, the update is
 *
 *   N_new = (N_old + dt_alpha*eta_N)
 *           / (1 + dt_alpha*kappa_a_N/Gamma_N).
 *
 * Gamma_N is not a fluid Lorentz factor and a valid value may lie below one.
 *
 * Scattering-only (kappa_a_N = 0) leaves N unchanged. Validation
 * rejects dt_alpha < 0, Gamma_N at/below its configured floor, nonfinite
 * inputs, and denominator <= 0. Floors
 * are NOT applied inside the formula; the host should run the repair helper
 * afterward if N_floor is configured.
 *
 * @param nu_params Neutrino number parameters, including N_floor and
 *        Gamma_N_floor.
 * @param rates Validated frozen rates for one species.
 * @param dt_alpha Lapse-weighted coordinate timestep in code-time units.
 * @param Gamma_N Endpoint number-current normalization; it is not a fluid
 *        Lorentz factor.
 * @param N_in Undensitized input number density.
 * @param N_out Output undensitized number density; unchanged on error.
 * @return @c ghl_success on publication; otherwise an M1 validation error.
 */
ghl_error_codes_t ghl_m1_update_neutrino_number_backward_euler(
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt_alpha,
      const double Gamma_N,
      const double N_in,
      double *restrict N_out);

/**
 * Try the opt-in thin explicit source update for one neutrino species.
 *
 * This validation-boundary helper is opt-in: it reports
 * `*thin_inequalities_hold = 0` and leaves a transactional no-op output when
 * `dt_alpha*kappa_a_E < 1` or `dt_alpha*kappa_s < 1` is false, where
 * `dt_alpha = metric->lapse * dt`. When both inequalities hold, it applies
 * the explicit E/F interaction source, then performs the backward-Euler N
 * update using the repaired E/F endpoint current, as in the branched source
 * policy's thin path, and assembles one equal-and-opposite exchange packet.
 * This compatibility wrapper intentionally does not report repair
 * diagnostics; use the `_with_diagnostics` variant when repair accounting is
 * required. Hosts remain responsible for their coupled source limiter and
 * stage policy; the ordinary implicit source solver is unchanged. When
 * enforce_mean_energy_bounds is enabled, the final repaired endpoint is
 * checked with the same rejection-only bounds as the other source paths; no
 * bound is clamped and N == 0 retains the existing skip.
 *
 * @param m1_params Initialized shared M1 parameters.
 * @param nu_params Neutrino number parameters, including N_floor and optional
 *        final mean-energy bounds.
 * @param metric Cell metric used for lapse weighting and exchange
 *        densitization.
 * @param prims Frozen fluid primitives for closure/current evaluation.
 * @param rates Validated frozen rates for one species.
 * @param dt Coordinate-time timestep in code-time units; the helper forms
 *        dt_alpha = metric->lapse * dt.
 * @param n_b_cons Positive conserved baryon-density normalization.
 * @param state_in Undensitized input state in {N,E,Fx,Fy,Fz} order.
 * @param thin_inequalities_hold Required output flag: zero when the thin
 *        inequalities do not select this branch, one when they do.
 * @param state_out Output state, initialized transactionally from @p state_in.
 * @param exchange Output exchange packet; zero when the branch is not selected
 *        or an error prevents publication.
 * @return @c ghl_success for a successful update or a valid non-selection;
 *         otherwise the output state and exchange remain transactional.
 */
ghl_error_codes_t ghl_m1_try_neutrino_explicit_thin_update(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double n_b_cons,
      const ghl_m1_neutrino_state *restrict state_in,
      int *restrict thin_inequalities_hold,
      ghl_m1_neutrino_state *restrict state_out,
      ghl_m1_neutrino_exchange *restrict exchange);

/**
 * Diagnostics-aware form of ghl_m1_try_neutrino_explicit_thin_update.
 *
 * The optional diagnostics record is updated transactionally with the
 * candidate repair. Pass NULL to obtain the same repair-reporting behavior as
 * the compatibility wrapper. Final repaired endpoint mean-energy bounds have
 * rejection semantics identical to the ordinary wrapper.
 *
 * @param m1_params Initialized shared M1 parameters.
 * @param nu_params Neutrino number parameters, including N_floor and optional
 *        final mean-energy bounds.
 * @param metric Cell metric used for lapse weighting and exchange
 *        densitization.
 * @param prims Frozen fluid primitives for closure/current evaluation.
 * @param rates Validated frozen rates for one species.
 * @param dt Coordinate-time timestep in code-time units; the helper forms
 *        dt_alpha = metric->lapse * dt.
 * @param n_b_cons Positive conserved baryon-density normalization.
 * @param state_in Undensitized input state in {N,E,Fx,Fy,Fz} order.
 * @param thin_inequalities_hold Required output flag: zero for a non-selection,
 *        one when the thin inequalities select this branch.
 * @param state_out Output state, initialized from @p state_in.
 * @param exchange Output exchange packet; zero when no update is published.
 * @param diagnostics Optional repair-budget accumulator; NULL is allowed.
 * @return @c ghl_success for a successful update or valid non-selection;
 *         otherwise the output state and exchange remain transactional.
 */
ghl_error_codes_t ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double n_b_cons,
      const ghl_m1_neutrino_state *restrict state_in,
      int *restrict thin_inequalities_hold,
      ghl_m1_neutrino_state *restrict state_out,
      ghl_m1_neutrino_exchange *restrict exchange,
      ghl_m1_neutrino_diagnostics *restrict diagnostics);

/** Transactional pointwise source-policy dispatcher.
 *
 * state_input is the pre-transport stage state and state_transport is the
 * source base. Returned radiation increments are always measured from
 * state_transport. Metric, primitives, and rates remain frozen. The caller
 * supplies the coordinate-time stage timestep (the implementation forms
 * `dt_alpha = metric->lapse * dt`) and conserved baryon-number normalization;
 * this operation owns neither stage counters nor matter publication. The
 * source options select whether dYe_matter follows the default charged-current
 * packet or signed total-number bookkeeping. For the branched policy,
 * thermalized_number_threshold controls the separate endpoint number
 * projection: negative disables it, while zero selects it even for zero
 * opacity or dt. The selected final repaired endpoint is then checked against
 * any enabled mean-energy bounds; violations are rejected rather than
 * clamped.
 *
 * @param options Optional source policy; NULL selects the established implicit
 *        path and its default options.
 * @param m1_params Initialized shared M1 parameters.
 * @param nu_params Neutrino number, floor, and mean-energy parameters.
 * @param metric Frozen cell metric.
 * @param prims_frozen Frozen fluid primitives.
 * @param rates Validated frozen rates for one species.
 * @param state_input Pre-transport input state used for public validation.
 * @param state_transport Transport-predicted state used as the source base.
 * @param dt Coordinate-time stage timestep in code-time units.
 * @param n_b_cons Positive conserved baryon-density normalization.
 * @param state_out Output state, initialized from @p state_transport.
 * @param exchange Output radiation/matter exchange packet.
 * @param diagnostics Required branch and implicit-solver diagnostics output.
 * @param neutrino_diagnostics Required caller-owned neutrino diagnostics
 *        accumulator.
 * @return @c ghl_success on publication; otherwise state and exchange remain
 *         at their transactional initialization and the error identifies the
 *         failed validation or source path.
 */
ghl_error_codes_t ghl_m1_solve_neutrino_source_update(
      const ghl_m1_neutrino_source_options *restrict options,
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_input,
      const ghl_m1_neutrino_state *restrict state_transport,
      const double dt,
      const double n_b_cons,
      ghl_m1_neutrino_state *restrict state_out,
      ghl_m1_neutrino_exchange *restrict exchange,
      ghl_m1_neutrino_source_diagnostics *restrict diagnostics,
      ghl_m1_neutrino_diagnostics *restrict neutrino_diagnostics);

/**
 * Transactional coupled electron-neutrino pair source update.
 *
 * The two array entries are ordered {ghl_m1_neutrino_nue,
 * ghl_m1_neutrino_anue}. Independent charged-current and scattering sources
 * are advanced first from the immutable state_transport source base (with
 * state_input retained for the ordinary stage validation), then the shared
 * pair reaction is advanced with one common number increment per substep.
 * Pair processes contribute no dYe_matter. The caller supplies n_b_cons as
 * the undensitized Eulerian baryon number density, W*rho/m_b; a densitized
 * host variable must be converted before this call.
 *
 * After required pointer validation, outputs are transactional: validation,
 * both source stages, endpoint
 * bounds, and exchange assembly must succeed before either species is
 * published. Number-floor repair is rejected rather than injected into the
 * conserving update.
 *
 * This is a first-order split grey collision model, with number-current
 * normalizations and partner occupancies frozen during their subsolves; see
 * Radiation/PAIR_SOURCE_MODEL.md. The call already includes independent CC
 * and scattering sources. The host must not apply these sources again.
 * Both diagnostics arrays are required; neutrino_diagnostics contains
 * caller-initialized accumulators. Exhausted Newton substep schedules return
 * ghl_error_m1_implicit_terminal_fallback with both source bases and zero
 * exchange packets.
 */
ghl_error_codes_t ghl_m1_solve_neutrino_pair_source_update(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters nu_params[2],
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates rates[2],
      const ghl_m1_neutrino_state state_input[2],
      const ghl_m1_neutrino_state state_transport[2],
      const double dt,
      const double n_b_cons,
      ghl_m1_neutrino_state state_out[2],
      ghl_m1_neutrino_exchange exchange[2],
      ghl_m1_neutrino_source_diagnostics diagnostics[2],
      ghl_m1_neutrino_diagnostics neutrino_diagnostics[2]);

/**
 * Compute the matter Delta Y_e recommendation from a radiation number
 * increment and the species lepton weight. The increment passed here must be
 * the charged-current electron-lepton exchange dL_rad_cc, not the total
 * radiation number update when pair/thermal channels are present.
 *
 * dL_rad_cc is already species-signed, so
 *   Delta Y_e_matter = - Delta L_rad_cc / baryon_density_conserved
 *
 * Validation rejects wrong lepton weight, nonfinite inputs, and
 * baryon_density_conserved <= 0. The argument is the undensitized Eulerian
 * baryon number density, W*rho/m_b. A host using a densitized baryon-mass
 * variable must convert its normalization before calling this helper.
 *
 * @param rates Validated rates whose lepton_weight supplies the species sign.
 * @param dL_rad_cc Species-signed charged-current radiation lepton-number
 *        increment.
 * @param baryon_density_conserved Positive conserved baryon-density
 *        normalization in the host's chosen compatible units.
 * @param dYe_matter Output matter composition increment recommendation;
 *        unchanged on error.
 * @return @c ghl_success on publication; otherwise an M1 validation error.
 */
ghl_error_codes_t ghl_m1_compute_neutrino_lepton_increment(
      const ghl_m1_neutrino_rates *restrict rates,
      const double dL_rad_cc,
      const double baryon_density_conserved,
      double *restrict dYe_matter);

/**
 * Solve the local homogeneous implicit update for one neutrino species.
 *
 * Uses frozen primitives and frozen rates; does not call Con2Prim and does
 * not update matter variables. Returns the updated state_out, populates the
 * exchange struct with radiation increments and equal-and-opposite matter
 * increments, and records solver and neutrino diagnostics. The host applies
 * the exchange increments after coupled limiting.
 *
 * The solve is transactional. After required pointers are validated, every
 * non-success return restores state_out to state_in and returns a zero exchange
 * packet. Diagnostics are observability state: success increments
 * source_converged, exhausted retry schedules increment
 * source_terminal_fallbacks, and hard failures increment source_failures.
 * The only accepted terminal_fallback_policy is no_update_all. n_b_cons is the
 * explicit conserved baryon-number density used for the signed Y_e packet.
 * The final repaired endpoint is checked against enabled mean-energy bounds;
 * an out-of-bounds endpoint is rejected without clamping, and N == 0 retains
 * the existing skip of that ratio check.
 *
 * @param m1_params Initialized shared M1 parameters.
 * @param nu_params Neutrino number, floor, and mean-energy parameters.
 * @param metric Frozen cell metric used for lapse weighting and densitization.
 * @param prims_frozen Frozen fluid primitives used by closure/current calls.
 * @param rates Validated frozen rates for one species.
 * @param dt Coordinate-time timestep in code-time units.
 * @param n_b_cons Positive conserved baryon-density normalization.
 * @param state_in Undensitized input state in {N,E,Fx,Fy,Fz} order.
 * @param state_out Output state; initialized from @p state_in and committed
 *        only on success.
 * @param exchange Output exchange packet; zero on a transactional failure or
 *        terminal no-update return.
 * @param solve_diagnostics Required output diagnostics for the implicit solve.
 * @param neutrino_diagnostics Required caller-owned neutrino diagnostics
 *        accumulator.
 * @return @c ghl_success on publication;
 *         @c ghl_error_m1_implicit_terminal_fallback when the retry schedule
 *         exhausts with no update; otherwise an error and transactional
 *         outputs.
 */
ghl_error_codes_t ghl_m1_solve_neutrino_implicit_homogeneous_update(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double n_b_cons,
      const ghl_m1_neutrino_state *restrict state_in,
      ghl_m1_neutrino_state *restrict state_out,
      ghl_m1_neutrino_exchange *restrict exchange,
      ghl_m1_implicit_solve_diagnostics *restrict solve_diagnostics,
      ghl_m1_neutrino_diagnostics *restrict neutrino_diagnostics);

/* Debug/public diagnostic neutrino implicit helpers. These declarations expose
 * the testing and source-policy diagnostic anchors without promoting the internal
 * _with_base fallback variants. */
ghl_error_codes_t ghl_m1_neutrino_compute_implicit_residual(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_in,
      const double dt,
      const double U[4],
      double residual[4]);

ghl_error_codes_t ghl_m1_neutrino_compute_implicit_jacobian(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_in,
      const double dt,
      const double U[4],
      const double residual_0[4],
      double jacobian[4][4]);

ghl_error_codes_t ghl_m1_neutrino_build_trial_state(
      const ghl_metric_quantities *restrict metric,
      const double U[4],
      ghl_m1_rad_state *restrict rad_state);

ghl_error_codes_t ghl_m1_neutrino_check_trial_admissibility(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_rad_state *restrict rad_state);

#ifdef __cplusplus
}
#endif

/** @} */

#endif // GHL_M1_H
