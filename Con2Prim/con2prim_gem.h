//#ifndef SMALLB_H_
//#define SMALLB_H_
//static const int SMALLBT=0,SMALLBX=1,SMALLBY=2,SMALLBZ=3,SMALLB2=4,NUMVARS_SMALLB=5;
//#endif //SMALLB_H_

#define smallb_defined = 1;
#define harm_params_defined = 1;

#ifndef CON2PRIM_GEM_H_
#define CON2PRIM_GEM_H_

#include "../GRHayL_Core/GRHayL.h"

//----------------- Ugly variables -----------------
// TODO: These things should be gotten rid of if possible.

static const int SMALLBT=0,SMALLBX=1,SMALLBY=2,SMALLBZ=3,SMALLB2=4,NUMVARS_SMALLB=5;
static const int NPR =8;
static const int NDIM=4;

/************* This version of vars are from Leo ****************/
//static const int MAX_NEWT_ITER       = 50;     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
//static const double NEWT_TOL      = 5e-9;    /* Min. of tolerance allowed for Newton-Raphson iterations */
//static const double MIN_NEWT_TOL  = 5e-9;    /* Max. of tolerance allowed for Newton-Raphson iterations */
/****************************************************************/

static const int MAX_NEWT_ITER    = 30;     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
//#define MAX_NEWT_ITER 300     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
static const double NEWT_TOL      = 1.0e-10;    /* Min. of tolerance allowed for Newton-Raphson iterations */
static const double MIN_NEWT_TOL  = 1.0e-10;    /* Max. of tolerance allowed for Newton-Raphson iterations */
static const int EXTRA_NEWT_ITER  = 0; /* ZACH SAYS: Original value = 2. But I don't think this parameter > 0 is warranted. Just slows the code for no reason, since our tolerances are fine. */

static const double NEWT_TOL2     = 1.0e-15;      /* TOL of new 1D^*_{v^2} gnr2 method */
static const double MIN_NEWT_TOL2 = 1.0e-10;  /* TOL of new 1D^*_{v^2} gnr2 method */

static const double W_TOO_BIG    = 1.e20;    /* \gamma^2 (\rho_0 + u + p) is assumed
                                                  to always be smaller than this.  This
                                                  is used to detect solver failures */
static const double UTSQ_TOO_BIG = 1.e20;    /* \tilde{u}^2 is assumed to be smaller
                                                  than this.  Used to detect solver
                                                  failures */

static const double FAIL_VAL     = 1.e30;    /* Generic value to which we set variables when a problem arises */

static const double NUMEPSILON   = 2.2204460492503131e-16;

//TODO: This had an error and I don't feel like debugging it rn.
//enum con2prim_routines{None=-1, Noble2D, Noble1D, Noble1D_entropy, Noble1D_entropy2,
//                       CerdaDuran2D, CerdaDuran3D, Palenzuela1D, Palenzuela1D_entropy,
//                       Newman1D};

static const int None = -1;
static const int Noble2D = 0;
static const int Noble1D = 1;
static const int Noble1D_entropy = 2;
static const int Noble1D_entropy2 = 3;
static const int CerdaDuran2D = 4;
static const int CerdaDuran3D = 5;
static const int Palenzuela1D = 6;
static const int Palenzuela1D_entropy = 7;
static const int Newman1D = 8;
static const int OldNoble2D = 20;

//--------------------------------------------------

//------------ Con2Prim structs --------------------
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
      1M: tau~ was reset in apply_inequality_fixes
      10M: S~ was reset in apply_inequality_fixes via the first case
      100M: S~ was reset in apply_inequality_fixes via the second case
For bitwise, would become 1, 2, 4, 8, 16, 32. 64, 128, and 256
https://www.tutorialspoint.com/cprogramming/c_bitwise_operators.htm
*/

typedef struct con2prim_diagnostics {
  int failures;
  int failure_checker;
  int font_fixes;
  int vel_limited_ptcount;
  int atm_resets;
  int which_routine;
  int rho_star_fix_applied;
  int pointcount;
  int failures_inhoriz;
  int pointcount_inhoriz;
  int backup[3];
  int nan_found;
  int c2p_fail_flag;
  double error_int_numer;
  double error_int_denom;
  int n_iter;
} con2prim_diagnostics;

//--------------------------------------------------

//------------- Misplaced functions ---------------
/* TODO: These functions come from the IllinoisGRMHD_header, EOS_header,
   or inlined_functions.h and
   are used by functions outside of con2prim. As such, they should be
   provided by something above c2p in the hierarchy. However, they use the
   cons/prim structs. The cons and prims structs should therefore be
   defined higher up in the code hierarchy. */

void reset_prims_to_atmosphere( 
             const GRHayL_parameters *restrict params,
             const eos_parameters *restrict eos,
             const metric_quantities *restrict metric,
             primitive_quantities *restrict prims,
             con2prim_diagnostics *restrict diagnostics);

//--------------------------------------------------

//--------- Initialization routines ----------------

void initialize_diagnostics(con2prim_diagnostics *restrict diagnostics);

//--------------------------------------------------

void con2prim_loop_kernel(
             const GRHayL_parameters *restrict params,
             const eos_parameters *restrict eos,
             metric_quantities *restrict metric,
             conservative_quantities *restrict cons,
             primitive_quantities *restrict prims,
             con2prim_diagnostics *restrict diagnostics,
             stress_energy *restrict Tmunu);

int C2P_Select_Hybrid_Method(
             const GRHayL_parameters *restrict params,
             const eos_parameters *restrict eos, const int c2p_key,
             const metric_quantities *restrict metric,
             const conservative_quantities *restrict cons,
             primitive_quantities *restrict prims,
             con2prim_diagnostics *restrict diagnostics);

int C2P_Hybrid_Noble2D(
             const GRHayL_parameters *restrict params,
             const eos_parameters *restrict eos,
             const metric_quantities *restrict metric,
             const conservative_quantities *restrict cons,
             primitive_quantities *restrict prim,
             con2prim_diagnostics *restrict diagnostics);

int C2P_Hybrid_OldNoble2D(
             const eos_parameters *restrict eos,
             const metric_quantities *restrict metric,
             const conservative_quantities *restrict cons,
             primitive_quantities *restrict prim,
             con2prim_diagnostics *restrict diagnostics);

int  apply_inequality_fixes(
             const GRHayL_parameters *restrict params,
             const eos_parameters *restrict eos,
             metric_quantities *restrict metric,
             const primitive_quantities *restrict prims,
             conservative_quantities *restrict cons,
             con2prim_diagnostics *restrict diagnostics);

void undensitize_conservatives(
             const metric_quantities *restrict metric,
             const conservative_quantities *restrict cons,
             conservative_quantities *restrict cons_undens);

void guess_primitives(
             const eos_parameters *restrict eos,
             const metric_quantities *restrict metric,
             const primitive_quantities *restrict prims,
             const conservative_quantities *restrict cons,
             primitive_quantities *restrict prims_guess);

void limit_utilde_and_compute_v(
             const eos_parameters *restrict eos,
             const metric_quantities *restrict metric,
             double *restrict utcon1_ptr, double *restrict u0_ptr,
             double *restrict utcon2_ptr, double *restrict utcon3_ptr,
             primitive_quantities *restrict prims,
             con2prim_diagnostics *restrict diagnostics);

void limit_v_and_output_u0(
             const eos_parameters *restrict eos,
             const metric_quantities *restrict metric,
             primitive_quantities *restrict prims,
             double *restrict u0,
             con2prim_diagnostics *restrict diagnostics);

void eigenvalues_3by3_real_sym_matrix(
             double *restrict  lam1, double *restrict  lam2, double *restrict  lam3,
             const double M11, const double M12, const double M13,
             const double M22, const double M23, const double M33);

void enforce_primitive_limits_and_output_u0(
             const GRHayL_parameters *restrict params,
             const eos_parameters *restrict eos,
             const metric_quantities *restrict metric,
             primitive_quantities *restrict prims,
             double *restrict u0,
             con2prim_diagnostics *restrict diagnostics);

void compute_conservs_and_Tmunu(
             const GRHayL_parameters *restrict params,
             const eos_parameters *restrict eos,
             const metric_quantities *restrict metric,
             primitive_quantities *restrict prims,
             const double u0,
             conservative_quantities *restrict cons,
             stress_energy *restrict Tmunu);

void initial_random_data(
             const double xrho,
             const double xpress,
             const bool random_metric,
             metric_quantities *restrict metric,
             primitive_quantities *restrict prims);

//--------------------------------------------------

//-------------- Font Fix routines -----------------
  
int font_fix(const eos_parameters *restrict eos,
             const metric_quantities *restrict metric,
             const conservative_quantities *restrict cons_undens,
             const primitive_quantities *restrict prims,
             primitive_quantities *restrict prims_guess,
             con2prim_diagnostics *restrict diagnostics);

int font_fix_hybrid_EOS(
             const eos_parameters *restrict eos,
             const metric_quantities *restrict metric,
             const conservative_quantities *restrict cons_undens,
             const primitive_quantities *restrict prims,
             double *restrict u_x_ptr, double *restrict u_y_ptr,
             double *restrict u_z_ptr );


//TODO: this was formerly inline. Include inside font_fix_hybrid_EOS directly?
int font_fix_rhob_loop(
             const eos_parameters *restrict eos,
             const int maxits, const double tol, const double W_in,
             const double Sf2_in, const double Psim6, const double sdots,
             const double BbardotS2, const double B2bar,
             const conservative_quantities *restrict cons,
             const double rhob_in, double *restrict rhob_out_ptr );

//--------------------------------------------------

/*TODO: use code below as basis to abstract con2prim calls
// A normal function with an int parameter
// and void return type
void fun(int a)
{
    printf("Value of a is %d\n", a);
}
  
int main()
{
    // fun_ptr is a pointer to function fun() 
    void (*fun_ptr)(int) = &fun;
  
    // The above line is equivalent of following two
    // void (*fun_ptr)(int);
    // fun_ptr = &fun; 
    
  
    // Invoking fun() using fun_ptr
    (*fun_ptr)(10);
  
    return 0;
}
*/

// TODO: The following functions are inside the functions that call them because they are only used by one function in c2p. They are used elsewhere inIGM,
// so we should consider a better solution.

// In enforce_limits_on_primitives_and_recompute_conservs.c: (split to enforce_primitive_limits... and compute_conservs...)
// impose_speed_limit_output_u0(const GRHayL_parameters *restrict params, const metric_quantities *restrict metric,
//                                   primitive_quantities *restrict prims, con2prim_diagnostics *restrict diagnostics,
//                                   double *restrict u0_out);
// compute_smallba_b2_and_u_i_over_u0_psi4(const metric_quantities *restrict metric, const primitive_quantities *restrict prims,
//                                                    const double u0L, const double ONE_OVER_LAPSE_SQRT_4PI, double *restrict u_x_over_u0_psi4,
//                                                    double *restrict u_y_over_u0_psi4, double *restrict u_z_over_u0_psi4, double *restrict smallb);

// In apply_tau_floor.c:
// eigenvalues_3by3_real_sym_matrix(double *restrict  lam1, double *restrict  lam2, double *restrict  lam3,
//                                      const double M11, const double M12, const double M13,
//                                      const double M22, const double M23, const double M33); 

#endif // CON2PRIM_GEM_H
