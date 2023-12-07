#include "../../../utils_Noble.h"

/*************************************************************************************

utoprim_1d_ee2.c:
---------------

  -- uses eq. (27) of Noble  et al. or the "momentum equation" and ignores
        the energy equation (29) in order to use the additional EOS, which
        is

             P = Sc rho^(GAMMA-1) / gamma

    Uses a method similiar to 1D_W method:
       -- solves for one independent variable (rho) via a 1D
          Newton-Raphson method
       -- by substituting
          W = Dc ( Dc + GAMMA Sc rho^(GAMMA-1) / (GAMMA-1) ) / rho
          into Qtsq equation, one can get one equation for
          one unknown (rho)

       -- can be used (in principle) with a general equation of state.

  -- Currently returns with an error state (>0) if a negative rest-mass
      density or internal energy density is calculated.  You may want
      to change this aspect of the code so that it still calculates the
      velocity and so that you can floor the densities.  If you want to
      change this aspect of the code please comment out the "return(retval)"
      statement after "retval = 5;" statement in Utoprim_new_body();

******************************************************************************/

/**********************************************************************************

  ghl_hybrid_Noble1D_entropy2():

     -- Attempt an inversion from U to prim using the initial guess prim.

     -- This is the main routine that calculates auxiliary quantities for the
        Newton-Raphson routine.

  -- assumes that
             /   rho gamma   \
         U = | alpha T^t_\mu |
             \   alpha B^i   /



             /     rho     \
      prim = |     uu      |
             | \tilde{u}^i |
             \  alpha B^i  /


return: i where
        i = 0 -> success
            1 -> initial v^2 < 0 with initial primitive guess;
            2 -> Newton-Raphson solver did not converge to a solution with the
                 given tolerances;
TODO: needs remaining error codes

**********************************************************************************/

int ghl_hybrid_Noble1D_entropy2(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics ) {

  double gnr_out[1];

  harm_aux_vars_struct harm_aux;

  double rho0, Z_last;
  if( ghl_initialize_Noble_entropy(params, eos, ADM_metric, metric_aux,
                                   cons_undens, prims, &harm_aux, &rho0, &Z_last) )
    return 1;


  gnr_out[0] = rho0;

  const int retval = ghl_general_newton_raphson(eos, &harm_aux, 1, Z_last, gnr_out, ghl_validate_1D_entropy, ghl_func_rho2);

  rho0 = gnr_out[0];

  /* Problem with solver, so return denoting error before doing anything further */
  if(retval != 0) {
    return retval;
  } else if(rho0 < 0 || rho0 > eos->rho_max) {
    return 4;
  }

  // Calculate v^2:

  const double rel_err = (harm_aux.D != 0.0) ? fabs((harm_aux.D-rho0)/harm_aux.D) :
                         (     (rho0 != 0.0) ? fabs((harm_aux.D-rho0)/rho0) : 0.0);
  const double utsq = ( rel_err > 1e-15 ) ? (harm_aux.D-rho0)*(harm_aux.D+rho0)/(rho0*rho0) : 0.0;

  if( utsq < 0. ) {
    return 5;
  }

  // Recover the primitive variables from the scalars and conserved variables:
  const double Wsq = 1.0+utsq;
  const double W   = sqrt(Wsq);

  const double Gamma_ppoly = eos->Gamma_ppoly[ghl_hybrid_find_polytropic_index(eos, rho0)];
  const double Gm1 = Gamma_ppoly - 1.0;
  const double rho_Gm1 = pow(rho0,Gm1);
  const double p_final = cons_undens->entropy * rho_Gm1/W;
  const double w = rho0 + p_final/Gm1 + p_final;
  const double Z = w * Wsq;

  prims->rho = rho0;

  if(prims->rho <= 0.0) {
    return 7;
  }

  ghl_finalize_Noble_entropy(params, eos, ADM_metric, metric_aux, cons_undens, &harm_aux, Z, W, prims);
  if(prims->press <= 0.0) {
    return 6;
  }

  /* Done! */
  diagnostics->n_iter = harm_aux.n_iter;
  diagnostics->which_routine = Noble1D_entropy2;
  return 0;
}
