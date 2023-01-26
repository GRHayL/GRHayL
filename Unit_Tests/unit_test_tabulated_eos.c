#include "NRPyEOS_Tabulated.h"
#include "unit_tests.h"

/*
 * (c) 2023 Leo Werneck
 */
int main(int argc, char **argv) {

  if( argc != 2 )
    grhayl_error("Correct usage is %s <table path>\n", argv[0]);

  grhayl_info("Beginning tabulated EOS unit test...\n");

  // Step 1: Initialize the EOS struct
  eos_parameters eos;
  initialize_tabulated_functions(&eos);
  eos.tabulated_read_table_set_EOS_params(argv[1], &eos);

  if( eos.N_rho != 7 || eos.N_T != 5 || eos.N_Ye != 3 )
    grhayl_error("Table dimension error: expected 7 x 5 x 3, but got %d x %d x %d\n",
                 eos.N_rho, eos.N_T, eos.N_Ye);
  grhayl_info("Table dimensions read in correctly\n");

  if( relative_error(eos.energy_shift, 123e-45) > 1e15 )
    grhayl_error("Error in energy shift exceeds round-off: %.15e vs. %.15e\n",
                 eos.energy_shift, 123e-45);
  grhayl_info("Energy shift read in correctly\n");

  // Step 2: Begin test
#pragma omp parallel for
  for(int i=0; i<eos.N_rho; i++) {
    const double logrho = eos.table_logrho[i];
    for(int j=0; j<eos.N_T; j++) {
      const double logT = eos.table_logT[j];
      for(int k=0; k<eos.N_Ye; k++) {
        const double Y_e = eos.table_Ye[k];
        const int index = i + eos.N_rho*( j + eos.N_T*k );
        for(int var_key=0; var_key<NRPyEOS_ntablekeys; var_key++) {
          const double var       = get_table_quantity(var_key, logrho, Y_e, logT);
          const double table_var = eos.table_all[var_key + NRPyEOS_ntablekeys*index];
          if( relative_error(var, table_var) > 2e-14 && fabs(var-table_var) > 1e-15 )
            grhayl_error("Error in variable %d exceeds round-off: %.15e vs. %.15e\n",
                         var_key, var, table_var);
        }
      }
    }
  }
  grhayl_info("All table quantities read in correctly!\n");

  grhayl_info("Beginning interpolation tests...\n");

  // Set number of points
  const int N_rho = 17;
  const int N_T   = 13;
  const int N_Ye  = 11;

  // Set useful quantities
  const double log_rho_min = log(eos.table_rho_min);
  const double log_rho_max = log(eos.table_rho_max);
  const double log_T_min   = log(eos.table_T_min);
  const double log_T_max   = log(eos.table_T_max);

  // Set step sizes
  const double dlogrho = (log_rho_max      - log_rho_min     )/N_rho;
  const double dlogT   = (log_T_max        - log_T_min       )/N_T;
  const double dYe     = (eos.table_Ye_max - eos.table_Ye_min)/N_Ye;

  // Begin test
  for(int k=0;k<N_Ye;k++) {
    const double Y_e = eos.table_Ye_min + k*dYe;
    for(int j=0;j<N_T;j++) {
      const double logT = log_T_min + j*dlogT;
      const double T    = exp(logT);
      for(int i=0;i<N_rho;i++) {
        const double logrho = log_rho_min + i*dlogrho;
        const double rho    = exp(logrho);

        // Compute all analytic quantities
        const double P      = exp(get_table_quantity(NRPyEOS_press_key  , logrho, Y_e, logT));
        const double eps    = exp(get_table_quantity(NRPyEOS_eps_key    , logrho, Y_e, logT));
        const double S      =     get_table_quantity(NRPyEOS_entropy_key, logrho, Y_e, logT);
        const double cs2    =     get_table_quantity(NRPyEOS_cs2_key    , logrho, Y_e, logT);
        const double depsdT =     get_table_quantity(NRPyEOS_depsdT_key , logrho, Y_e, logT);
        const double muhat  =     get_table_quantity(NRPyEOS_muhat_key  , logrho, Y_e, logT);
        const double mu_e   =     get_table_quantity(NRPyEOS_mu_e_key   , logrho, Y_e, logT);
        const double mu_p   =     get_table_quantity(NRPyEOS_mu_p_key   , logrho, Y_e, logT);
        const double mu_n   =     get_table_quantity(NRPyEOS_mu_n_key   , logrho, Y_e, logT);
        const double X_p    =     get_table_quantity(NRPyEOS_X_p_key    , logrho, Y_e, logT);
        const double X_n    =     get_table_quantity(NRPyEOS_X_n_key    , logrho, Y_e, logT);

        // Declare interpolated quantities
        double P_interp, eps_interp, T_interp, S_interp, cs2_interp;
        double muhat_interp, mu_e_interp, mu_p_interp, mu_n_interp;
        double depsdT_interp, X_p_interp, X_n_interp;

        // Now perform the interpolations, validating the results
        P_interp = 0.0/0.0;
        eos.tabulated_compute_P_from_T(&eos, rho, Y_e, T, &P_interp);
        if( relative_error(P, P_interp) > 2e-14 )
          grhayl_error("tabulated_compute_P_from_T validation failed:"
                       "Pressure : %22.15e %22.15e : %22.15e\n",
                       P, P_interp, relative_error(P, P_interp));

        eps_interp = 0.0/0.0;
        eos.tabulated_compute_eps_from_T(&eos, rho, Y_e, T, &eps_interp);
        if( relative_error(eps, eps_interp) > 2e-14 )
          grhayl_error("tabulated_compute_eps_from_T validation failed:"
                       "Energy   : %22.15e %22.15e : %22.15e\n",
                       eps, eps_interp, relative_error(eps, eps_interp));

        P_interp = eps_interp = 0.0/0.0;
        eos.tabulated_compute_P_eps_from_T(&eos, rho, Y_e, T, &P_interp, &eps_interp);
        if( relative_error(P  , P_interp  ) > 2e-14 ||
            relative_error(eps, eps_interp) > 2e-14 )
          grhayl_error("tabulated_compute_P_eps_from_T validation failed:\n"
                       "Pressure : %22.15e %22.15e : %22.15e\n"
                       "Energy   : %22.15e %22.15e : %22.15e\n",
                       P  , P_interp  , relative_error(P  , P_interp  ),
                       eps, eps_interp, relative_error(eps, eps_interp));

        P_interp = eps_interp = S_interp = 0.0/0.0;
        eos.tabulated_compute_P_eps_S_from_T(&eos, rho, Y_e, T, &P_interp, &eps_interp, &S_interp);
        if( relative_error(P  , P_interp  ) > 2e-14 ||
            relative_error(eps, eps_interp) > 2e-14 ||
            relative_error(S  , S_interp  ) > 2e-14 )
          grhayl_error("tabulated_compute_P_eps_S_from_T validation failed:\n"
                       "Pressure : %22.15e %22.15e : %22.15e\n"
                       "Energy   : %22.15e %22.15e : %22.15e\n"
                       "Entropy  : %22.15e %22.15e : %22.15e\n",
                       P  , P_interp  , relative_error(P  , P_interp  ),
                       eps, eps_interp, relative_error(eps, eps_interp),
                       S  , S_interp  , relative_error(S  , S_interp  ));

        P_interp = eps_interp = S_interp = cs2_interp = 0.0/0.0;
        eos.tabulated_compute_P_eps_S_cs2_from_T(&eos, rho, Y_e, T, &P_interp, &eps_interp, &S_interp, &cs2_interp);
        if( relative_error(P  , P_interp  ) > 2e-14 ||
            relative_error(eps, eps_interp) > 2e-14 ||
            relative_error(S  , S_interp  ) > 2e-14 ||
            relative_error(cs2, cs2_interp) > 2e-13 )
          grhayl_error("tabulated_compute_P_eps_S_cs2_from_T validation failed:\n"
                       "Pressure : %22.15e %22.15e : %22.15e\n"
                       "Energy   : %22.15e %22.15e : %22.15e\n"
                       "Entropy  : %22.15e %22.15e : %22.15e\n"
                       "cs2      : %22.15e %22.15e : %22.15e\n",
                       P  , P_interp  , relative_error(P  , P_interp  ),
                       eps, eps_interp, relative_error(eps, eps_interp),
                       S  , S_interp  , relative_error(S  , S_interp  ),
                       cs2, cs2_interp, relative_error(cs2, cs2_interp));

        P_interp = eps_interp = depsdT_interp = 0.0/0.0;
        eos.tabulated_compute_P_eps_depsdT_from_T(&eos, rho, Y_e, T, &P_interp, &eps_interp, &depsdT_interp);
        if( relative_error(P     , P_interp     ) > 2e-14 ||
            relative_error(eps   , eps_interp   ) > 2e-14 ||
            relative_error(depsdT, depsdT_interp) > 2e-14 )
          grhayl_error("tabulated_compute_P_eps_depsdT_from_T validation failed:\n"
                       "Pressure : %22.15e %22.15e : %22.15e\n"
                       "Energy   : %22.15e %22.15e : %22.15e\n"
                       "deps/dT  : %22.15e %22.15e : %22.15e\n",
                       P     , P_interp     , relative_error(P     , P_interp     ),
                       eps   , eps_interp   , relative_error(eps   , eps_interp   ),
                       depsdT, depsdT_interp, relative_error(depsdT, depsdT_interp));

        P_interp = eps_interp = muhat_interp = mu_e_interp = mu_p_interp = mu_n_interp = 0.0/0.0;
        eos.tabulated_compute_P_eps_muhat_mue_mup_mun_from_T(&eos, rho, Y_e, T,
                                                             &P_interp, &eps_interp,
                                                             &muhat_interp, &mu_e_interp,
                                                             &mu_p_interp, &mu_n_interp);
        if( relative_error(P    , P_interp    ) > 2e-14 ||
            relative_error(eps  , eps_interp  ) > 2e-14 ||
            relative_error(muhat, muhat_interp) > 2e-14 ||
            relative_error(mu_e , mu_e_interp ) > 6e-14 ||
            relative_error(mu_p , mu_p_interp ) > 2e-14 ||
            (relative_error(mu_n , mu_n_interp ) > 2e-14 && fabs(mu_n-mu_n_interp)>2e-14) )
          grhayl_error("tabulated_compute_P_eps_muhat_mue_mup_mun_from_T validation failed:\n"
                       "Pressure : %22.15e %22.15e : %22.15e\n"
                       "Energy   : %22.15e %22.15e : %22.15e\n"
                       "muhat    : %22.15e %22.15e : %22.15e\n"
                       "mu_e     : %22.15e %22.15e : %22.15e\n"
                       "mu_p     : %22.15e %22.15e : %22.15e\n"
                       "mu_n     : %22.15e %22.15e : %22.15e\n",
                       P    , P_interp    , relative_error(P    , P_interp    ),
                       eps  , eps_interp  , relative_error(eps  , eps_interp  ),
                       muhat, muhat_interp, relative_error(muhat, muhat_interp),
                       mu_e , mu_e_interp , relative_error(mu_e , mu_e_interp ),
                       mu_p , mu_p_interp , relative_error(mu_p , mu_p_interp ),
                       mu_n , mu_n_interp , relative_error(mu_n , mu_n_interp ));


        muhat_interp = mu_e_interp = mu_p_interp = mu_n_interp = X_n_interp = X_p_interp = 0.0/0.0;
        eos.tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T(&eos, rho, Y_e, T,
                                                             &muhat_interp, &mu_e_interp, &mu_p_interp,
                                                             &mu_n_interp, &X_n_interp, &X_p_interp);
        if( relative_error(muhat, muhat_interp) > 2e-14 ||
            relative_error(mu_e , mu_e_interp ) > 6e-14 ||
            relative_error(mu_p , mu_p_interp ) > 2e-14 ||
            (relative_error(mu_n, mu_n_interp ) > 2e-14 && fabs(mu_n-mu_n_interp)>2e-14) ||
            relative_error(X_n  , X_n_interp  ) > 6e-14 ||
            relative_error(X_p  , X_p_interp  ) > 2e-14 )
          grhayl_error("tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T validation failed:\n"
                       "muhat : %22.15e %22.15e : %22.15e\n"
                       "mu_e  : %22.15e %22.15e : %22.15e\n"
                       "mu_p  : %22.15e %22.15e : %22.15e\n"
                       "mu_n  : %22.15e %22.15e : %22.15e\n"
                       "X_n   : %22.15e %22.15e : %22.15e\n"
                       "X_p   : %22.15e %22.15e : %22.15e\n",
                       muhat, muhat_interp, relative_error(muhat, muhat_interp),
                       mu_e , mu_e_interp , relative_error(mu_e , mu_e_interp ),
                       mu_p , mu_p_interp , relative_error(mu_p , mu_p_interp ),
                       mu_n , mu_n_interp , relative_error(mu_n , mu_n_interp ),
                       X_n  , X_n_interp  , relative_error(X_n  , X_n_interp  ),
                       X_p  , X_p_interp  , relative_error(X_p  , X_p_interp  ));

        T_interp = eos.table_T_min; P_interp = 0.0/0.0;
        eos.tabulated_compute_P_T_from_eps(&eos, rho, Y_e, eps, &P_interp, &T_interp);
        if( relative_error(T, T_interp) > 3e-14 ||
            relative_error(P, P_interp) > 2e-14 )
          grhayl_error("tabulated_compute_P_T_from_eps validation failed:\n"
                       "Temperature : %22.15e %22.15e : %22.15e\n"
                       "Pressure    : %22.15e %22.15e : %22.15e\n",
                       T, T_interp, relative_error(T, T_interp),
                       P, P_interp, relative_error(P, P_interp));

        T_interp = eos.table_T_min; P_interp = S_interp = depsdT_interp = 0.0/0.0;
        eos.tabulated_compute_P_S_depsdT_T_from_eps(&eos, rho, Y_e, eps, &P_interp, &S_interp, &depsdT_interp, &T_interp);
        if( relative_error(T     , T_interp     ) > 3e-14 ||
            relative_error(P     , P_interp     ) > 2e-14 ||
            relative_error(S     , S_interp     ) > 7e-14 ||
            relative_error(depsdT, depsdT_interp) > 7e-14 )
          grhayl_error("tabulated_compute_P_S_depsdT_T_from_eps validation failed:\n"
                       "Temperature : %22.15e %22.15e : %22.15e\n"
                       "Pressure    : %22.15e %22.15e : %22.15e\n"
                       "Entropy     : %22.15e %22.15e : %22.15e\n"
                       "deps/dT     : %22.15e %22.15e : %22.15e\n",
                       T     , T_interp     , relative_error(T     , T_interp     ),
                       P     , P_interp     , relative_error(P     , P_interp     ),
                       S     , S_interp     , relative_error(S     , S_interp     ),
                       depsdT, depsdT_interp, relative_error(depsdT, depsdT_interp));

        T_interp = eos.table_T_min; eps_interp = S_interp = 0.0/0.0;
        eos.tabulated_compute_eps_S_T_from_P(&eos, rho, Y_e, P, &eps_interp, &S_interp, &T_interp);
        if( relative_error(T  , T_interp  ) > 3e-13 ||
            relative_error(eps, eps_interp) > 3e-13 ||
            relative_error(S  , S_interp  ) > 3e-13 )
          grhayl_error("tabulated_compute_eps_S_T_from_P validation failed:\n"
                       "Temperature : %22.15e %22.15e : %22.15e\n"
                       "Energy      : %22.15e %22.15e : %22.15e\n"
                       "Entropy     : %22.15e %22.15e : %22.15e\n",
                       T  , T_interp  , relative_error(T  , T_interp  ),
                       eps, eps_interp, relative_error(eps, eps_interp),
                       S  , S_interp  , relative_error(S  , S_interp  ));

        T_interp = eos.table_T_min; P_interp = eps_interp = 0.0/0.0;
        eos.tabulated_compute_P_eps_T_from_S(&eos, rho, Y_e, S, &P_interp, &eps_interp, &T_interp);
        if( relative_error(T  , T_interp  ) > 3e-14 ||
            relative_error(P  , P_interp  ) > 2e-14 ||
            relative_error(eps, eps_interp) > 2e-14 )
          grhayl_error("tabulated_compute_P_eps_T_from_S validation failed:\n"
                       "Temperature : %22.15e %22.15e : %22.15e\n"
                       "Pressure    : %22.15e %22.15e : %22.15e\n"
                       "Energy      : %22.15e %22.15e : %22.15e\n",
                       T  , T_interp  , relative_error(T  , T_interp  ),
                       P  , P_interp  , relative_error(P  , P_interp  ),
                       eps, eps_interp, relative_error(eps, eps_interp));
      }
    }
  }

  // Step 4: Free memory
  eos.tabulated_free_memory(&eos);

  grhayl_info("Test finished with no errors!\n");

  // All done!
  return 0;
}
