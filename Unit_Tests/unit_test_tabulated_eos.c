#include "nrpyeos_tabulated.h"
#include "unit_tests.h"

/*
 * (c) 2023 Leo Werneck
 */
int main(int argc, char **argv) {

  double rtol = 1.5e-12;
  double atol = 1e-14;

  if( argc < 2 || argc > 4 ) {
    ghl_error("Correct usage is %s <table path> [rel. err. tolerance] [abs. err. tolerance]\n", argv[0]);
  }
  else if( argc > 2 ) {
    rtol = strtod(argv[1], NULL);
    if( argc > 3 )
      atol = strtod(argv[2], NULL);
  }

  ghl_info("Beginning tabulated EOS unit test...\n");

  // Step 1: Initialize the EOS struct
  ghl_eos_parameters eos;
  ghl_initialize_tabulated_eos_functions_and_params(argv[1], exp(1), -1, -1, 1, -1, -1, exp(1), -1, -1, &eos);

  if( eos.N_rho != 7 || eos.N_T != 5 || eos.N_Ye != 3 )
    ghl_error("Table dimension error: expected 7 x 5 x 3, but got %d x %d x %d\n",
                 eos.N_rho, eos.N_T, eos.N_Ye);
  ghl_info("Table dimensions read in correctly\n");

  if( relative_error(eos.energy_shift, 123e-45) > rtol )
    ghl_error("Error in energy shift exceeds tolerance: %.15e vs. %.15e\n",
                 eos.energy_shift, 123e-45);
  ghl_info("Energy shift read in correctly\n");

  // Step 2: Begin test
  for(int i=0; i<eos.N_rho; i++) {
    const double logrho = eos.table_logrho[i];
    for(int j=0; j<eos.N_T; j++) {
      const double logT = eos.table_logT[j];
      for(int k=0; k<eos.N_Ye; k++) {
        const double Y_e = eos.table_Y_e[k];
        const int index = i + eos.N_rho*( j + eos.N_T*k );
        for(int var_key=0; var_key<NRPyEOS_ntablekeys; var_key++) {
          const double var       = get_table_quantity(var_key, logrho, Y_e, logT);
          const double table_var = eos.table_all[var_key + NRPyEOS_ntablekeys*index];
          if( relative_error(var, table_var) > rtol && fabs(var-table_var) > atol )
            ghl_error("Errors in variable %d exceed tolernaces: %.15e vs. %.15e\n",
                         var_key, var, table_var);
        }
      }
    }
  }
  ghl_info("All table quantities read in correctly!\n");

  ghl_info("Beginning interpolation tests...\n");

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
  const double dlogrho = (log_rho_max       - log_rho_min      )/N_rho;
  const double dlogT   = (log_T_max         - log_T_min        )/N_T;
  const double dYe     = (eos.table_Y_e_max - eos.table_Y_e_min)/N_Ye;

  // Begin test
  for(int k=0;k<N_Ye;k++) {
    const double Y_e = eos.table_Y_e_min + k*dYe;
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
        ghl_tabulated_compute_P_from_T(&eos, rho, Y_e, T, &P_interp);
        if( relative_error(P, P_interp) > rtol && fabs(P - P_interp) > atol )
          ghl_error("tabulated_compute_P_from_T validation failed:\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       " Varname :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       "Pressure : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n",
                       P, P_interp, relative_error(P, P_interp), fabs(P - P_interp));

        eps_interp = 0.0/0.0;
        ghl_tabulated_compute_eps_from_T(&eos, rho, Y_e, T, &eps_interp);
        if( relative_error(eps, eps_interp) > rtol && fabs(eps - eps_interp) > atol )
          ghl_error("tabulated_compute_eps_from_T validation failed:\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       " Varname :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       "Energy   : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n",
                       eps, eps_interp, relative_error(eps, eps_interp), fabs(eps - eps_interp));

        cs2_interp = 0.0/0.0;
        ghl_tabulated_compute_cs2_from_T(&eos, rho, Y_e, T, &cs2_interp);
        if( relative_error(cs2, cs2_interp) > rtol && fabs(cs2 - cs2_interp) > atol )
          ghl_error("tabulated_compute_cs2_from_T validation failed:\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       " Varname :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       "cs2      : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n",
                       cs2, cs2_interp, relative_error(cs2, cs2_interp), fabs(cs2 - cs2_interp));

        P_interp = eps_interp = 0.0/0.0;
        ghl_tabulated_compute_P_eps_from_T(&eos, rho, Y_e, T, &P_interp, &eps_interp);
        if( ( relative_error(P  , P_interp  ) > rtol && fabs(P   - P_interp  ) > atol ) ||
            ( relative_error(eps, eps_interp) > rtol && fabs(eps - eps_interp) > atol ))
          ghl_error("tabulated_compute_P_eps_from_T validation failed:\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       " Varname :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       "Pressure : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Energy   : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n",
                       P  , P_interp  , relative_error(P  , P_interp  ), fabs(P   - P_interp  ),
                       eps, eps_interp, relative_error(eps, eps_interp), fabs(eps - eps_interp));

        P_interp = eps_interp = S_interp = 0.0/0.0;
        ghl_tabulated_compute_P_eps_S_from_T(&eos, rho, Y_e, T, &P_interp, &eps_interp, &S_interp);
        if( ( relative_error(P  , P_interp  ) > rtol && fabs(P   - P_interp  ) > atol ) ||
            ( relative_error(eps, eps_interp) > rtol && fabs(eps - eps_interp) > atol ) ||
            ( relative_error(S  , S_interp  ) > rtol && fabs(S   - S_interp  ) > atol ) )
          ghl_error("tabulated_compute_P_eps_S_from_T validation failed:\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       " Varname :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       "Pressure : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Energy   : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Entropy  : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n",
                       P  , P_interp  , relative_error(P  , P_interp  ), fabs(P   - P_interp  ),
                       eps, eps_interp, relative_error(eps, eps_interp), fabs(eps - eps_interp),
                       S  , S_interp  , relative_error(S  , S_interp  ), fabs(S   - S_interp  ));

        P_interp = eps_interp = S_interp = cs2_interp = 0.0/0.0;
        ghl_tabulated_compute_P_eps_S_cs2_from_T(&eos, rho, Y_e, T, &P_interp, &eps_interp, &S_interp, &cs2_interp);
        if( ( relative_error(P  , P_interp  ) > rtol && fabs(P   - P_interp  ) > atol ) ||
            ( relative_error(eps, eps_interp) > rtol && fabs(eps - eps_interp) > atol ) ||
            ( relative_error(S  , S_interp  ) > rtol && fabs(S   - S_interp  ) > atol ) ||
            ( relative_error(cs2, cs2_interp) > rtol && fabs(cs2 - cs2_interp) > atol ) )
          ghl_error("tabulated_compute_P_eps_S_cs2_from_T validation failed:\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       " Varname :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       "Pressure : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Energy   : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Entropy  : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "cs2      : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n",
                       P  , P_interp  , relative_error(P  , P_interp  ), fabs(P   - P_interp  ),
                       eps, eps_interp, relative_error(eps, eps_interp), fabs(eps - eps_interp),
                       S  , S_interp  , relative_error(S  , S_interp  ), fabs(S   - S_interp  ),
                       cs2, cs2_interp, relative_error(cs2, cs2_interp), fabs(cs2 - cs2_interp));

        P_interp = eps_interp = depsdT_interp = 0.0/0.0;
        ghl_tabulated_compute_P_eps_depsdT_from_T(&eos, rho, Y_e, T, &P_interp, &eps_interp, &depsdT_interp);
        if( ( relative_error(P     , P_interp     ) > rtol && fabs(P      - P_interp     ) > atol ) ||
            ( relative_error(eps   , eps_interp   ) > rtol && fabs(eps    - eps_interp   ) > atol ) ||
            ( relative_error(depsdT, depsdT_interp) > rtol && fabs(depsdT - depsdT_interp) > atol ) )
          ghl_error("tabulated_compute_P_eps_depsdT_from_T validation failed:\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       " Varname :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       "Pressure : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Energy   : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "deps/dT  : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n",
                       P     , P_interp     , relative_error(P     , P_interp     ), fabs(P      - P_interp      ),
                       eps   , eps_interp   , relative_error(eps   , eps_interp   ), fabs(eps    - eps_interp    ),
                       depsdT, depsdT_interp, relative_error(depsdT, depsdT_interp), fabs(depsdT - depsdT_interp ));

        P_interp = eps_interp = muhat_interp = mu_e_interp = mu_p_interp = mu_n_interp = 0.0/0.0;
        ghl_tabulated_compute_P_eps_muhat_mue_mup_mun_from_T(&eos, rho, Y_e, T,
                                                             &P_interp, &eps_interp,
                                                             &muhat_interp, &mu_e_interp,
                                                             &mu_p_interp, &mu_n_interp);
        if( ( relative_error(P    , P_interp    ) > rtol && fabs(P     - P_interp    ) > atol ) ||
            ( relative_error(eps  , eps_interp  ) > rtol && fabs(eps   - eps_interp  ) > atol ) ||
            ( relative_error(muhat, muhat_interp) > rtol && fabs(muhat - muhat_interp) > atol ) ||
            ( relative_error(mu_e , mu_e_interp ) > rtol && fabs(mu_e  - mu_e_interp ) > atol ) ||
            ( relative_error(mu_p , mu_p_interp ) > rtol && fabs(mu_p  - mu_p_interp ) > atol ) ||
            ( relative_error(mu_n , mu_n_interp ) > rtol && fabs(mu_n  - mu_n_interp ) > atol ) )
          ghl_error("tabulated_compute_P_eps_muhat_mue_mup_mun_from_T validation failed:\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       " Varname :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       "Pressure : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Energy   : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "muhat    : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "mu_e     : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "mu_p     : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "mu_n     : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n",
                       P    , P_interp    , relative_error(P    , P_interp    ), fabs(P     - P_interp    ),
                       eps  , eps_interp  , relative_error(eps  , eps_interp  ), fabs(eps   - eps_interp  ),
                       muhat, muhat_interp, relative_error(muhat, muhat_interp), fabs(muhat - muhat_interp),
                       mu_e , mu_e_interp , relative_error(mu_e , mu_e_interp ), fabs(mu_e  - mu_e_interp ),
                       mu_p , mu_p_interp , relative_error(mu_p , mu_p_interp ), fabs(mu_p  - mu_p_interp ),
                       mu_n , mu_n_interp , relative_error(mu_n , mu_n_interp ), fabs(mu_n  - mu_n_interp ));


        muhat_interp = mu_e_interp = mu_p_interp = mu_n_interp = X_n_interp = X_p_interp = 0.0/0.0;
        ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T(&eos, rho, Y_e, T,
                                                             &muhat_interp, &mu_e_interp, &mu_p_interp,
                                                             &mu_n_interp, &X_n_interp, &X_p_interp);
        if( ( relative_error(muhat, muhat_interp) > rtol && fabs(muhat - muhat_interp) > atol ) ||
            ( relative_error(mu_e , mu_e_interp ) > rtol && fabs(mu_e  - mu_e_interp ) > atol ) ||
            ( relative_error(mu_p , mu_p_interp ) > rtol && fabs(mu_p  - mu_p_interp ) > atol ) ||
            ( relative_error(mu_n , mu_n_interp ) > rtol && fabs(mu_n  - mu_n_interp ) > atol ) ||
            ( relative_error(X_n  , X_n_interp  ) > rtol && fabs(X_n   - X_n_interp  ) > atol ) ||
            ( relative_error(X_p  , X_p_interp  ) > rtol && fabs(X_p   - X_p_interp  ) > atol ) )
          ghl_error("tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T validation failed:\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       " Varname :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n"
                       "muhat : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "mu_e  : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "mu_p  : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "mu_n  : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "X_n   : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "X_p   : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "---------:------------------------:------------------------:------------------------:-----------------------\n",
                       muhat, muhat_interp, relative_error(muhat, muhat_interp), fabs(muhat - muhat_interp),
                       mu_e , mu_e_interp , relative_error(mu_e , mu_e_interp ), fabs(mu_e  - mu_e_interp ),
                       mu_p , mu_p_interp , relative_error(mu_p , mu_p_interp ), fabs(mu_p  - mu_p_interp ),
                       mu_n , mu_n_interp , relative_error(mu_n , mu_n_interp ), fabs(mu_n  - mu_n_interp ),
                       X_n  , X_n_interp  , relative_error(X_n  , X_n_interp  ), fabs(X_n   - X_n_interp  ),
                       X_p  , X_p_interp  , relative_error(X_p  , X_p_interp  ), fabs(X_p   - X_p_interp  ));

        T_interp = eos.table_T_min; P_interp = 0.0/0.0;
        ghl_tabulated_compute_P_T_from_eps(&eos, rho, Y_e, eps, &P_interp, &T_interp);
        if( ( relative_error(T, T_interp) > rtol && fabs(T - T_interp) > atol ) ||
            ( relative_error(P, P_interp) > rtol && fabs(P - P_interp) > atol ) )
          ghl_error("tabulated_compute_P_T_from_eps validation failed:\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "  Varname   :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "Temperature : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Pressure    : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n",
                       T, T_interp, relative_error(T, T_interp), fabs(T - T_interp),
                       P, P_interp, relative_error(P, P_interp), fabs(P - P_interp));

        T_interp = eos.table_T_min; P_interp = S_interp = depsdT_interp = 0.0/0.0;
        ghl_tabulated_compute_P_S_depsdT_T_from_eps(&eos, rho, Y_e, eps, &P_interp, &S_interp, &depsdT_interp, &T_interp);
        if( ( relative_error(T     , T_interp     ) > rtol && fabs(T      - T_interp     ) > atol ) ||
            ( relative_error(P     , P_interp     ) > rtol && fabs(P      - P_interp     ) > atol ) ||
            ( relative_error(S     , S_interp     ) > rtol && fabs(S      - S_interp     ) > atol ) ||
            ( relative_error(depsdT, depsdT_interp) > rtol && fabs(depsdT - depsdT_interp) > atol ) )
          ghl_error("tabulated_compute_P_S_depsdT_T_from_eps validation failed:\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "  Varname   :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "Temperature : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Pressure    : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Entropy     : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "deps/dT     : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n",
                       T     , T_interp     , relative_error(T     , T_interp     ), fabs(T      - T_interp     ),
                       P     , P_interp     , relative_error(P     , P_interp     ), fabs(P      - P_interp     ),
                       S     , S_interp     , relative_error(S     , S_interp     ), fabs(S      - S_interp     ),
                       depsdT, depsdT_interp, relative_error(depsdT, depsdT_interp), fabs(depsdT - depsdT_interp));

        T_interp = eos.table_T_min; eps_interp = S_interp = 0.0/0.0;
        ghl_tabulated_compute_eps_S_T_from_P(&eos, rho, Y_e, P, &eps_interp, &S_interp, &T_interp);
        if( ( relative_error(T  , T_interp  ) > rtol && fabs(T   - T_interp  ) > atol ) ||
            ( relative_error(eps, eps_interp) > rtol && fabs(eps - eps_interp) > atol ) ||
            ( relative_error(S  , S_interp  ) > rtol && fabs(S   - S_interp  ) > atol ) )
          ghl_error("tabulated_compute_eps_S_T_from_P validation failed:\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "  Varname   :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "Temperature : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Energy      : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Entropy     : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n",
                       T  , T_interp  , relative_error(T  , T_interp  ), fabs(T   - T_interp  ),
                       eps, eps_interp, relative_error(eps, eps_interp), fabs(eps - eps_interp),
                       S  , S_interp  , relative_error(S  , S_interp  ), fabs(S   - S_interp  ));

        T_interp = eos.table_T_min; P_interp = eps_interp = 0.0/0.0;
        ghl_tabulated_compute_P_eps_T_from_S(&eos, rho, Y_e, S, &P_interp, &eps_interp, &T_interp);
        if( ( relative_error(T  , T_interp  ) > rtol && fabs(T   - T_interp  ) > atol ) ||
            ( relative_error(P  , P_interp  ) > rtol && fabs(P   - P_interp  ) > atol ) ||
            ( relative_error(eps, eps_interp) > rtol && fabs(eps - eps_interp) > atol ) )
          ghl_error("tabulated_compute_P_eps_T_from_S validation failed:\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "  Varname   :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "Temperature : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Pressure    : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Energy      : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n",
                       T  , T_interp  , relative_error(T  , T_interp  ), fabs(T   - T_interp  ),
                       P  , P_interp  , relative_error(P  , P_interp  ), fabs(P   - P_interp  ),
                       eps, eps_interp, relative_error(eps, eps_interp), fabs(eps - eps_interp));

        T_interp = eos.table_T_min; P_interp = 0.0/0.0;
        ghl_tabulated_compute_P_T_from_S(&eos, rho, Y_e, S, &P_interp, &T_interp);
        if( ( relative_error(T  , T_interp  ) > rtol && fabs(T   - T_interp  ) > atol ) ||
            ( relative_error(P  , P_interp  ) > rtol && fabs(P   - P_interp  ) > atol ) )
          ghl_error("tabulated_compute_P_eps_T_from_S validation failed:\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "  Varname   :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "Temperature : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Pressure    : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n",
                       T  , T_interp  , relative_error(T  , T_interp  ), fabs(T   - T_interp  ),
                       P  , P_interp  , relative_error(P  , P_interp  ), fabs(P   - P_interp  ));

        T_interp = eos.table_T_min; eps_interp = 0.0/0.0;
        ghl_tabulated_compute_eps_T_from_P(&eos, rho, Y_e, P, &eps_interp, &T_interp);
        if( ( relative_error(T  , T_interp  ) > rtol && fabs(T   - T_interp  ) > atol ) ||
            ( relative_error(eps, eps_interp) > rtol && fabs(eps - eps_interp) > atol ) )
          ghl_error("tabulated_compute_eps_T_from_P validation failed:\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "  Varname   :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "Temperature : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Energy      : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n",
                       T  , T_interp  , relative_error(T  , T_interp  ), fabs(T   - T_interp  ),
                       eps, eps_interp, relative_error(eps, eps_interp), fabs(eps - eps_interp));

        T_interp = eos.table_T_min; cs2_interp = P_interp = 0.0/0.0;
        ghl_tabulated_compute_P_cs2_T_from_eps(&eos, rho, Y_e, eps, &P_interp, &cs2_interp, &T_interp);
        if( ( relative_error(T  , T_interp  ) > rtol && fabs(T   - T_interp  ) > atol ) ||
            ( relative_error(cs2, cs2_interp) > rtol && fabs(cs2 - cs2_interp) > atol ) ||
            ( relative_error(P  , P_interp  ) > rtol && fabs(P   - P_interp  ) > atol ) )
          ghl_error("tabulated_compute_P_eps_T_from_S validation failed:\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "  Varname   :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "Temperature : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "cs2         : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Pressure    : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n",
                       T  , T_interp  , relative_error(T  , T_interp  ), fabs(T   - T_interp  ),
                       cs2, cs2_interp, relative_error(cs2, cs2_interp), fabs(cs2 - cs2_interp),
                       P  , P_interp  , relative_error(P  , P_interp  ), fabs(P   - P_interp  ));

        P_interp = cs2_interp = eps_interp = 0.0/0.0;
        ghl_tabulated_compute_P_eps_cs2_from_T(&eos, rho, Y_e, T, &P_interp, &eps_interp, &cs2_interp);
        if( ( relative_error(P  , P_interp  ) > rtol && fabs(P   - P_interp  ) > atol ) ||
            ( relative_error(cs2, cs2_interp) > rtol && fabs(cs2 - cs2_interp) > atol ) ||
            ( relative_error(eps, eps_interp) > rtol && fabs(eps - eps_interp) > atol ) )
          ghl_error("tabulated_compute_P_eps_T_from_S validation failed:\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "  Varname   :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "Pressure    : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "cs2         : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Energy      : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n",
                       P  , P_interp  , relative_error(P  , P_interp  ), fabs(P   - P_interp  ),
                       cs2, cs2_interp, relative_error(cs2, cs2_interp), fabs(cs2 - cs2_interp),
                       eps, eps_interp, relative_error(eps, eps_interp), fabs(eps - eps_interp));

        T_interp = eos.table_T_min; cs2_interp = eps_interp = 0.0/0.0;
        ghl_tabulated_compute_eps_cs2_T_from_P(&eos, rho, Y_e, P, &eps_interp, &cs2_interp, &T_interp);
        if( ( relative_error(T  , T_interp  ) > rtol && fabs(T   - T_interp  ) > atol ) ||
            ( relative_error(cs2, cs2_interp) > rtol && fabs(cs2 - cs2_interp) > atol ) ||
            ( relative_error(eps, eps_interp) > rtol && fabs(eps - eps_interp) > atol ) )
          ghl_error("tabulated_compute_P_eps_T_from_S validation failed:\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "  Varname   :     Analytic value     :  Interpolation Value   :     Relative Error     :     Absolute Error    \n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n"
                       "Temperature : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "cs2         : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "Energy      : %22.15e : %22.15e : %22.15e : %22.15e\n"
                       "------------:------------------------:------------------------:------------------------:-----------------------\n",
                       T  , T_interp  , relative_error(T  , T_interp  ), fabs(T   - T_interp  ),
                       cs2, cs2_interp, relative_error(cs2, cs2_interp), fabs(cs2 - cs2_interp),
                       eps, eps_interp, relative_error(eps, eps_interp), fabs(eps - eps_interp));
      }
    }
  }

  // Enforce limit tests
  {
    double rho = 0.9 * eos.rho_min;
    double Y_e = 0.9 * eos.Y_e_min;
    double P   = 0.9 * eos.press_min;
    ghl_tabulated_enforce_bounds_rho_Ye_P(&eos, &rho, &Y_e, &P);
    if(rho != eos.rho_min || Y_e != eos.Y_e_min || P != eos.press_min) {
      ghl_error("enforce bounds (rho, Y_e, P) failed for small values: %e != %e or %e != %e or %e != %e\n",
                rho, eos.rho_min, Y_e, eos.Y_e_min, P, eos.press_min);
    }
  }

  {
    double rho = 1.1 * eos.rho_max;
    double Y_e = 1.1 * eos.Y_e_max;
    double P   = 1.1 * eos.press_max;
    ghl_tabulated_enforce_bounds_rho_Ye_P(&eos, &rho, &Y_e, &P);
    if(rho != eos.rho_max || Y_e != eos.Y_e_max || P != eos.press_max) {
      ghl_error("enforce bounds (rho, Y_e, P) failed for large values");
    }
  }

  {
    double rho = 0.9 * eos.rho_max;
    double Y_e = 0.9 * eos.Y_e_max;
    double P   = 0.9 * eos.press_max;
    ghl_tabulated_enforce_bounds_rho_Ye_P(&eos, &rho, &Y_e, &P);
    if(rho != 0.9 * eos.rho_max || Y_e != 0.9 * eos.Y_e_max || P != 0.9 * eos.press_max) {
      ghl_error("enforce bounds (rho, Y_e, P) changed values that it shouldn't have");
    }
  }

  // Now test beta equilibrium stuff
  {
    ghl_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(exp(4), &eos);
    const double Y_e_expected[7] = {
      1.000000000000000e+00,
      1.000000000000000e+00,
      1.000000000000000e+00,
      1.000000000000000e+00,
      1.000000000000000e+00,
      1.000000000000000e+00,
      1.000000000000000e+00,
    };
    const double eps_expected[7] = {
      -2.000000000000000e+00,
      -9.999999999999933e-01,
      -4.440892098500627e-16,
      9.999999999999997e-01,
      2.000000000000000e+00,
      3.000000000000000e+00,
      4.000000000000006e+00,
    };
    const double P_expected[7] = {
      5.999999999999986e+00,
      7.000000000000014e+00,
      8.000000000000000e+00,
      9.000000000000000e+00,
      1.000000000000000e+01,
      1.100000000000000e+01,
      1.200000000000000e+01,
    };
    for(int i=0;i<eos.N_rho;i++) {
      if(eos.Ye_of_lr[i] == 0 || Y_e_expected[i] == 0) {
        if( fabs(eos.Ye_of_lr[i] - Y_e_expected[i]) > 1e-14 ) {
          ghl_error("Failed to impose beta equilibrium correctly: %d %.15e %.15e\n", i, eos.Ye_of_lr[i], Y_e_expected[i]);
        }
      }
      else {
        if(relative_error(eos.Ye_of_lr[i], Y_e_expected[i]) > 1e-14 ) {
          ghl_error("Failed to impose beta equilibrium correctly: %d %.15e %.15e\n", i, eos.Ye_of_lr[i], Y_e_expected[i]);
        }
      }
      if(eos.le_of_lr[i] == 0 || eps_expected[i] == 0) {
        if( fabs(eos.le_of_lr[i] - eps_expected[i]) > 1e-14 ) {
          ghl_error("Failed to impose beta equilibrium correctly: %d %.15e %.15e\n", i, eos.le_of_lr[i], eps_expected[i]);
        }
      }
      else {
        if(relative_error(eos.le_of_lr[i], eps_expected[i]) > 1e-14 ) {
          ghl_error("Failed to impose beta equilibrium correctly: %d %.15e %.15e\n", i, eos.le_of_lr[i], eps_expected[i]);
        }
      }
      if(eos.lp_of_lr[i] == 0 || P_expected[i] == 0) {
        if( fabs(eos.lp_of_lr[i] - P_expected[i]) > 1e-14 ) {
          ghl_error("Failed to impose beta equilibrium correctly: %d %.15e %.15e\n", i, eos.lp_of_lr[i], P_expected[i]);
        }
      }
      else {
        if(relative_error(eos.lp_of_lr[i], P_expected[i]) > 1e-14 ) {
          ghl_error("Failed to impose beta equilibrium correctly: %d %.15e %.15e\n", i, eos.lp_of_lr[i], P_expected[i]);
        }
      }
    }

    // Now test the interpolators
    const double rho_interp = exp(0.5*(eos.table_logrho[0] + eos.table_logrho[1]));

    const double Ye_interp = ghl_tabulated_compute_Ye_from_rho(&eos, rho_interp);
    const double Ye_expect = 0.5*(Y_e_expected[0] + Y_e_expected[1]);
    if(relative_error(Ye_interp, Ye_expect) > 1e-14) {
      ghl_error("Failed to interpolate Y_e(rho): %.15e %.15e\n", Ye_interp, Ye_expect);
    }

    const double P_interp   = log(ghl_tabulated_compute_P_from_rho(  &eos, rho_interp));
    const double P_expect   = 0.5*(P_expected[0] + P_expected[1]); 
    if(relative_error(P_interp, P_expect) > 1e-14) {
      ghl_error("Failed to interpolate P(rho): %.15e %.15e\n", P_interp, P_expect);
    }

    const double eps_interp = log(ghl_tabulated_compute_eps_from_rho(&eos, rho_interp));
    const double eps_expect = 0.5*(eps_expected[0] + eps_expected[1]); 
    if(relative_error(eps_interp, eps_expect) > 1e-14) {
      ghl_error("Failed to interpolate eps(rho): %.15e %.15e\n", eps_interp, eps_expect);
    }
  }

  // Step 4: Free memory
  ghl_tabulated_free_memory(&eos);

  ghl_info("Test finished with no errors!\n");

  // All done!
  return 0;
}
