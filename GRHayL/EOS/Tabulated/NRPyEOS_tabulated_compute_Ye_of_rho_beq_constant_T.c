#include "nrpyeos_tabulated.h"
#include "../../Con2Prim/roots.h" // Assuming ghl_brent and roots_params definitions are here

#define munu_index(ir, it, iy) \
  (NRPyEOS_munu_key            \
   + NRPyEOS_ntablekeys * ((ir) + eos->N_rho * ((it) + eos->N_T * (iy))))

#define mu_p_index(ir, it, iy) \
  (NRPyEOS_mu_p_key            \
   + NRPyEOS_ntablekeys * ((ir) + eos->N_rho * ((it) + eos->N_T * (iy))))

#define mu_n_index(ir, it, iy) \
  (NRPyEOS_mu_n_key            \
  + NRPyEOS_ntablekeys * ((ir) + eos->N_rho * ((it) + eos->N_T * (iy))))

#define entropy_index(ir, it, iy) \
  (NRPyEOS_entropy_key         \
    + NRPyEOS_ntablekeys * ((ir) + eos->N_rho * ((it) + eos->N_T * (iy))))

// Modified struct to pass Ye result back from brent_objective
typedef struct {
  double target_entropy;
  int ir;
  double T;   // <--- new field for temperature // Original comment preserved
  double Ye;  // Field to store Ye found by inner solver
  ghl_eos_parameters *eos;
} brent_params;

static double find_Ye_st_munu_is_zero(
      const int n,
      const double *restrict Ye,
      const double *restrict munu) {

  int i0 = -1;
  for(int i = 0; i < n - 1; i++) {
    if(munu[i] * munu[i + 1] < 0) {
      i0 = i;
      break;
    }
  }

  // If munu does not cross zero, return the minimum Ye
  if(i0 == -1) {
    return Ye[0];
  }

  const int i1 = i0 + 1;
  const double x0 = Ye[i0];
  const double x1 = Ye[i1];
  const double y0 = munu[i0];
  const double y1 = munu[i1];

  return (y0 * x1 - y1 * x0) / (y0 - y1);
}

static int find_left_index_uniform_array(
      const int nx,
      const double *restrict x_arr,
      const double x) {
  // Ensure x is within bounds before calculating index directly
  // Although linterp checks bounds, find_left_index might be called elsewhere.
  if (x < x_arr[0]) return 0; // Or handle error, returning 0 can lead to extrapolation
  if (x > x_arr[nx-1]) return nx-1; // Or handle error

  // Calculate index based on uniform spacing
  int index = (int)((x - x_arr[0]) / (x_arr[1] - x_arr[0]));

  // Clamp index to valid range [0, nx-2] for interpolation interval [index, index+1]
  // If x == x_arr[nx-1], index will be nx-1, clamp to nx-2.
  if (index >= nx - 1) {
       index = nx - 2;
   }
   // If x == x_arr[0], index will be 0, which is correct.
   if (index < 0) { // Should not happen if x >= x_arr[0]
       index = 0;
   }

  // Original calculation (prone to slight offset issues, replaced above)
  // return (x - x_arr[0]) / (x_arr[1] - x_arr[0]) + 0.5;
  return index;
}

static int find_left_index_bisection(const int nx, const double *restrict x_arr, const double x) {

  int ia = 0;
  double a = x_arr[ia] - x;
  if(a == 0) { // Exact match at the lower bound
    return ia;
  }

  int ib = nx - 1;
  double b = x_arr[ib] - x;
  if(b == 0) {
    // If x matches the last element x_arr[nx-1], the interval should be
    // [nx-2, nx-1], so we return the index nx-2 (or ib-1).
    return ib - 1;
  }

  if(a * b >= 0) {
     // Check bounds before declaring failure
     if (x < x_arr[ia]) {
          ghl_error("Bisection error: Value %g is below the lower bound %g (index %d)\n", x, x_arr[ia], ia);
     } else if (x > x_arr[ib]) {
          ghl_error("Bisection error: Value %g is above the upper bound %g (index %d)\n", x, x_arr[ib], ib);
     } else {
          // Value is within bounds, but a*b >= 0. Could be duplicate points or zero function value at boundary.
          ghl_error("Interval [%g, %g] does not bracket the root %g\n", x_arr[ia], x_arr[ib], x);
     }
  }

  // Bisection loop
  do {
    int ic = ia + (ib - ia) / 2; // Avoid potential overflow
    double c = x_arr[ic] - x;

    if(c == 0) { // Exact match at midpoint
      return ic;
    }
    // Use the lower bound value 'a' for consistent interval halving logic
    else if(a * c < 0) { // Root is in [ia, ic]
      ib = ic;
      // b = c; // Not strictly needed for loop logic
    }
    else { // Root is in (ic, ib]
      ia = ic;
      a = c; // Update the lower bound value
    }
  } while(ib - ia > 1);

  // Final check (original logic) - verifies the loop terminated correctly
  // This check is primarily for debugging the bisection logic itself.
  a = x_arr[ia] - x;
  b = x_arr[ib] - x;
  if(a > 0 || b < 0) {
     // This error suggests the loop invariant or update logic failed.
     ghl_error("Bisection post-check failed: %g not in final interval [%g, %g]\n", x, x_arr[ia], x_arr[ib]);
  }

  // Return the lower index 'ia' of the final bracket [ia, ib] where ib = ia + 1
  return ia;
}


// This is a simple linear interpolation (linterp) function.
static double linterp(
      const int nx,
      const double *restrict x_arr,
      const double *restrict y_arr,
      const double x,
      int (*find_left_index)(const int, const double *restrict, const double)) {

  if(x < x_arr[0] || x > x_arr[nx - 1]) {
    ghl_error("Point (%e) out of array bounds [%e, %e]\n", x, x_arr[0], x_arr[nx - 1]);
  }

  // Set up basic quantities for the interpolation
  const int i0 = find_left_index(nx, x_arr, x);
   // Add check for valid index before proceeding
   if (i0 < 0 || i0 >= nx -1 ) { // i0 must be in [0, nx-2] for interval [i0, i0+1]
       // If i0==nx-1 it means x==x_arr[nx-1], handle this case specifically if find_left_index doesn't return nx-2.
       // With the fix in find_left_index_bisection, i0 should be nx-2 if x == x_arr[nx-1].
       if (i0 == nx - 1 && x == x_arr[nx-1]) {
           // This case should ideally not be reached if find_left_index returns nx-2 for x=x_arr[nx-1].
           // If it does, return the boundary value.
            ghl_warn("linterp: find_left_index returned boundary index %d for x=%e. Returning boundary value.\n", i0, x);
            return y_arr[i0];
       }
       ghl_error("Invalid index %d returned by find_left_index for x=%e (nx=%d)\n", i0, x, nx);
   }

  const int i1 = i0 + 1;
  const double x0 = x_arr[i0];
  const double x1 = x_arr[i1];
  const double y0 = y_arr[i0];
  const double y1 = y_arr[i1];

  // Perform the linear interpolation
  const double denom = x1 - x0;
  // Prevent division by zero if x points are identical
  if(denom == 0.0) {
      return y0; // Or y1, or average. Returning y0 is typical.
  }
  return (y0 * (x1 - x) + y1 * (x - x0)) / denom;
}


// Compute Ye(logrho) using linear interpolation
double NRPyEOS_tabulated_compute_Ye_from_rho(
      const ghl_eos_parameters *restrict eos,
      const double rho) {

  return linterp(
        eos->N_rho, eos->table_logrho, eos->Ye_of_lr, log(rho),
        find_left_index_uniform_array);
}

// Interpolate logP(logrho) using linear interpolation
double NRPyEOS_tabulated_compute_P_from_rho(
      const ghl_eos_parameters *restrict eos,
      const double rho) {

  return exp(linterp(
        eos->N_rho, eos->table_logrho, eos->lp_of_lr, log(rho),
        find_left_index_uniform_array));
}

// Interpolate logrho(logP) using linear interpolation
double NRPyEOS_tabulated_compute_rho_from_P(
      const ghl_eos_parameters *restrict eos,
      const double P) {

  return exp(linterp(
        eos->N_rho, eos->lp_of_lr, eos->table_logrho, log(P),
        find_left_index_bisection));
}

// Interpolate logeps(logrho) using linear interpolation
double NRPyEOS_tabulated_compute_eps_from_rho(
      const ghl_eos_parameters *restrict eos,
      const double rho) {

  return exp(linterp(
               eos->N_rho, eos->table_logrho, eos->le_of_lr, log(rho),
               find_left_index_bisection))
         - eos->energy_shift;
}

void NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T(
      const double T,
      ghl_eos_parameters *restrict eos) {

  const int it = ghl_tabulated_get_index_T(eos, T); // Assumes this function exists
  const int nr = eos->N_rho;
  const int ny = eos->N_Ye;

  if(eos->Ye_of_lr == NULL) {
    eos->Ye_of_lr = (double *)malloc(sizeof(double) * nr);
     if (!eos->Ye_of_lr) ghl_error("Failed to allocate memory for Ye_of_lr.\n");
  }
  double *munu_of_Ye = (double *)malloc(sizeof(double) * ny);
  if (!munu_of_Ye) ghl_error("Failed to allocate memory for munu_of_Ye.\n");

  for(int ir = 0; ir < nr; ir++) {
    for(int iy = 0; iy < ny; iy++) {
      munu_of_Ye[iy] = eos->table_all[munu_index(ir, it, iy)];
    }
    eos->Ye_of_lr[ir] = find_Ye_st_munu_is_zero(ny, eos->table_Y_e, munu_of_Ye);
  }
  free(munu_of_Ye);
}

void NRPyEOS_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(
      const double T,
      ghl_eos_parameters *restrict eos) {

  // Start by obtaining Ye(logrho)
  NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T(T, eos);

  // Now allocate memory for logP(logrho) and logeps(logrho)
  if(eos->lp_of_lr == NULL) {
    eos->lp_of_lr = (double *)malloc(sizeof(double) * eos->N_rho);
    if (!eos->lp_of_lr) ghl_error("Failed to allocate memory for lp_of_lr.\n");
  }
  if(eos->le_of_lr == NULL) {
    eos->le_of_lr = (double *)malloc(sizeof(double) * eos->N_rho);
    if (!eos->le_of_lr) ghl_error("Failed to allocate memory for le_of_lr.\n");
  }

  // Compute logP(logrho) and logeps(logrho)
  for(int ir = 0; ir < eos->N_rho; ir++) {
    const double rho = exp(eos->table_logrho[ir]);
    const double Y_e = eos->Ye_of_lr[ir];
    double P, eps;
    ghl_tabulated_compute_P_eps_from_T(eos, rho, Y_e, T, &P, &eps); // Assumes this function exists
    eos->lp_of_lr[ir] = log(P);
    eos->le_of_lr[ir] = log(eps + eos->energy_shift);
  }
}

// This function computes the entropy at a given Ye and T.
static double entropy_at_Ye_T(
    double Ye,
    double T,
    const int ir,
    ghl_eos_parameters *restrict eos) {

  Ye = MAX(Ye, eos->Y_e_min);
  Ye = MIN(Ye, eos->Y_e_max);

  T = MAX(T, eos->T_min);
  T = MIN(T, eos->T_max);

  const int nt = eos->N_T;
  const int ny = eos->N_Ye;

  const int it = find_left_index_bisection(nt, eos->table_logT, log(T));

  // Check if index is valid for interval [it, it+1]
  if (it < 0 || it >= nt - 1) {
     // Note: With the fix to find_left_index_bisection, it should now correctly be nt-2
     // if log(T) == eos->table_logT[nt-1]. So this error should only trigger
     // if T is truly out of bounds or bisection failed unexpectedly.
     ghl_error("Error: Temperature index %d out of bounds [0, %d) after bisection for T=%g (logT=%g).\n", it, nt-1, T, log(T));
     return NAN; // Indicate an error
   }

  double *entropy_of_Ye = (double *)malloc(sizeof(double) * ny);
  if (entropy_of_Ye == NULL) {
    ghl_error("Memory allocation failed.\n");
    return NAN; // Indicate an error
  }

  for (int iy = 0; iy < ny; iy++) {
    entropy_of_Ye[iy] = eos->table_all[entropy_index(ir, it, iy)];
  }

  double entropy_val = linterp(ny, eos->table_Y_e, entropy_of_Ye, Ye, find_left_index_uniform_array);
  free(entropy_of_Ye);

  return entropy_val;
}

// This function computes the difference between mu_p and mu_n at a given Ye and T.
static double mu_p__minus__mu_n__at_Ye_T(
    double Ye,
    double T,
    const int ir,
    ghl_eos_parameters *restrict eos) {

  Ye = MAX(Ye, eos->Y_e_min);
  Ye = MIN(Ye, eos->Y_e_max);

  T = MAX(T, eos->T_min);
  T = MIN(T, eos->T_max);

  const int nt = eos->N_T;
  const int ny = eos->N_Ye;

  const int it = find_left_index_bisection(nt, eos->table_logT, log(T));

  // Check if index is valid for interval [it, it+1]
  if (it < 0 || it >= nt - 1) {
      // See comment in entropy_at_Ye_T regarding this check after bisection fix.
      ghl_error("Error: Temperature index %d out of bounds [0, %d) after bisection for T=%g (logT=%g) in mu_p-mu_n.\n", it, nt-1, T, log(T));
      return NAN; // Indicate an error
   }

  double *mu_p__minus__mu_n__of_Ye = (double *)malloc(sizeof(double) * ny);
  if (mu_p__minus__mu_n__of_Ye == NULL) {
    ghl_error("Memory allocation failed.\n");
    return NAN; // Indicate an error
  }

  for (int iy = 0; iy < ny; iy++) {
    // munu_of_Ye[iy] = eos->table_all[munu_index(ir, it, iy)];
    // The above is commented out because this is for beta-equilibrium for a cold, neutrino-less star.
    // Given constant entropy, with a temperature gradient, we aim to model a newly formed, hot proto-neutron star.
    // Thus, we need beta-equilibrium, with mu_n + mu_nu = mu_p + mu_e.
    // --> mu_p - mu_n = 0.
    mu_p__minus__mu_n__of_Ye[iy] = eos->table_all[mu_p_index(ir, it, iy)] - eos->table_all[mu_n_index(ir, it, iy)];
  }

  double mu_p__minus__mu_n_val = linterp(ny, eos->table_Y_e, mu_p__minus__mu_n__of_Ye, Ye, find_left_index_uniform_array);
  free(mu_p__minus__mu_n__of_Ye);

  if (isnan(mu_p__minus__mu_n_val)) {
      ghl_warn("mu_p__minus__mu_n__at_Ye_T: linterp returned NaN for Ye=%g, T=%g, ir=%d, it=%d\n", Ye, T, ir, it);
   }

  return mu_p__minus__mu_n_val;
}

// Add a wrapper for mu_p__minus__mu_n__at_Ye_T so that it matches the prototype expected by ghl_brent.
static double mu_p__minus__mu_n_wrapper(const double Ye,
                           const ghl_parameters *restrict params,
                           const ghl_eos_parameters *restrict eos,
                           const ghl_conservative_quantities *restrict cons_undens,
                           fparams_struct *restrict fparams,
                           ghl_primitive_quantities *restrict prims) {
  // Cast back fparams to our brent_params struct.
  brent_params *bparams = (brent_params *)fparams;
  // Ensure eos is accessed correctly, preferably via bparams->eos if consistent
  return mu_p__minus__mu_n__at_Ye_T(Ye, bparams->T, bparams->ir, bparams->eos);
}

static double brent_objective(
  const double T, // Current Temperature guess from outer solver
  const ghl_parameters *restrict params,
  const ghl_eos_parameters *restrict eos, // Passed via bparams usually
  const ghl_conservative_quantities *restrict cons_undens,
  fparams_struct *restrict fparams,
  ghl_primitive_quantities *restrict prims)
{
  // Cast the fparams pointer to retrieve our brent_params structure
  brent_params *bparams = (brent_params *)fparams;
  bparams->T = T; // Store current T guess for inner solver

  // Find Ye_beq for the current T guess using the inner solver
  roots_params ye_root_params = {0};
  ye_root_params.tol = 1e-10; // Adjust tolerance as needed
  ye_root_params.max_iters = 100;

  // Use the wrapper function for mu_p__minus__mu_n_wrapper solve
  // Need the correct function pointer type for ghl_brent
  double (*ye_solver_func)(const double, const ghl_parameters *, const ghl_eos_parameters *,
                           const ghl_conservative_quantities *, fparams_struct *,
                           ghl_primitive_quantities *);
  ye_solver_func = mu_p__minus__mu_n_wrapper;

  // Call inner Brent solver for Ye
  int ye_status = ghl_brent(ye_solver_func, params, bparams->eos, cons_undens, fparams, prims,
                            // Use slightly inset bounds for robustness
                             bparams->eos->table_Y_e[1],
                             bparams->eos->table_Y_e[bparams->eos->N_Ye - 2], // Avoid exact boundary
                             &ye_root_params);

  // Check status of inner Ye solver (assuming 0 means success)
  if (ye_status != ghl_success) {
       ghl_warn("Inner Ye solver failed in brent_objective. T=%g, ir=%d, status=%d\n", T, bparams->ir, ye_status);
       return 1.0e100; // Return large residual to outer solver
   }
  double ye_beq = ye_root_params.root;
  bparams->Ye = ye_beq; // Store the found Ye in bparams


  double entropy_beq = entropy_at_Ye_T(ye_beq, T, bparams->ir, bparams->eos);
  // Check for NaN from entropy calculation
   if (isnan(entropy_beq)) {
       ghl_warn("entropy_at_Ye_T returned NaN in brent_objective. T=%g, Ye=%g, ir=%d\n", T, ye_beq, bparams->ir);
       return 1.0e100; // Return large residual
   }

  return entropy_beq - bparams->target_entropy;
}


// This function computes Ye, P, and eps at a given entropy and rho.
// It uses the outer Brent solver to find T and the inner solver to find Ye.
// The results are stored in the eos structure.
void NRPyEOS_tabulated_compute_Ye_P_eps_of_rho_beq_constant_entropy(
    const double entropy,
    ghl_eos_parameters *restrict eos) {

  const int nr = eos->N_rho;

  // Allocate memory with checks
  if (eos->Ye_of_lr == NULL) {
    eos->Ye_of_lr = (double *)malloc(sizeof(double) * nr);
     if (!eos->Ye_of_lr) ghl_error("Failed to allocate memory for Ye_of_lr.\n");
  }
  if (eos->lp_of_lr == NULL) {
    eos->lp_of_lr = (double *)malloc(sizeof(double) * nr);
    if (!eos->lp_of_lr) ghl_error("Failed to allocate memory for lp_of_lr.\n");
  }
  if (eos->le_of_lr == NULL) {
    eos->le_of_lr = (double *)malloc(sizeof(double) * nr);
    if (!eos->le_of_lr) ghl_error("Failed to allocate memory for le_of_lr.\n");
  }

  roots_params T_root_params = {0};
  T_root_params.tol = 1e-10; // Adjust tolerance as needed
  T_root_params.max_iters = 100;

  brent_params bparams = {0};
  bparams.target_entropy = entropy;
  bparams.eos = eos;

  // Define and assign function pointer outside loop
  double (*entropy_func)(const double, const ghl_parameters *, const ghl_eos_parameters *, const ghl_conservative_quantities *, fparams_struct *, ghl_primitive_quantities *);
  // Cast needed to match ghl_brent signature
  entropy_func = brent_objective;


  for (int ir = 0; ir < nr; ir++) {
    bparams.ir = ir; // Set current density index

    // Find Temperature T_beq using outer Brent solver
    int T_status = ghl_brent(entropy_func, NULL, eos, NULL, (fparams_struct*)&bparams, NULL,
                             // Use table temperature bounds for search
                              eos->T_min,
                              eos->T_max,
                              &T_root_params); // Pass T solver params

    // Check status of outer T solver (assuming 0 means success)
     if (T_status != ghl_success) {
         ghl_error("Outer T solver failed for ir=%d, target S=%g, status=%d\n", ir, entropy, T_status);
         // Consider alternative error handling (e.g., fill with NaN, break)
        //  return; // Exit on failure
     }

    double temp_beq = T_root_params.root;

    // Retrieve Ye_beq stored by the final call to brent_objective
    double ye_beq = bparams.Ye;

    eos->Ye_of_lr[ir] = ye_beq; // Store result

    double P, eps;
    // Assumes ghl_tabulated_compute_P_eps_from_T exists and works correctly
    ghl_tabulated_compute_P_eps_from_T(eos, exp(eos->table_logrho[ir]), ye_beq, temp_beq, &P, &eps);
    eos->lp_of_lr[ir] = log(P);
    eos->le_of_lr[ir] = log(eps + eos->energy_shift); // Remember energy shift
  }
}

void NRPyEOS_tabulated_free_beq_quantities(ghl_eos_parameters *restrict eos) {
  if (eos->Ye_of_lr) {
    free(eos->Ye_of_lr);
    eos->Ye_of_lr = NULL;
  }
  if (eos->lp_of_lr) {
    free(eos->lp_of_lr);
    eos->lp_of_lr = NULL;
  }
  if (eos->le_of_lr) {
    free(eos->le_of_lr);
    eos->le_of_lr = NULL;
  }
}