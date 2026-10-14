#ifndef M1_PAIR_RANGE_TESTS_H
#define M1_PAIR_RANGE_TESTS_H

/* Exercise range rejection through the real paired solver, after the public
 * rate and current validation boundaries. No private-helper substitution. */
static void test_pair_quadratic_range_rejection(const ghl_m1_parameters *params) {
  const struct {
    double q, n_eq, N, dt;
  } cases[] = {
    { DBL_MIN, 1.0, 1.0, DBL_MIN }, /* H underflows. */
    { 1.0, DBL_MIN, 1.0, 1.0 },     /* D underflows. */
    { 1.0, 1.0, 1.0, 1.0e-310 },    /* D/H and B overflow. */
    { 1.0, 1.0e100, 1.0e200, 1.0 }  /* C overflows while B is finite. */
  };
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  ghl_primitive_quantities prims = { .rho = 1.0, .u0 = 1.0 };
  for(size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); ++c) {
    ghl_m1_neutrino_parameters nu[2] = { { .N_floor = 0.0 }, { .N_floor = 0.0 } };
    ghl_m1_neutrino_rates rates[2];
    ghl_m1_neutrino_state input[2], output[2];
    ghl_m1_neutrino_exchange exchange[2];
    ghl_m1_neutrino_source_diagnostics diagnostics[2];
    ghl_m1_neutrino_diagnostics nd[2] = { { 0 }, { 0 } };
    for(int s = 0; s < 2; ++s) {
      rates[s] = (ghl_m1_neutrino_rates){ .species = s == 0 ? ghl_m1_neutrino_nue
                                                            : ghl_m1_neutrino_anue,
                                          .lepton_weight = s == 0 ? 1.0 : -1.0,
                                          .mean_energy = 1.0,
                                          .n_eq = cases[c].n_eq,
                                          .J_eq = cases[c].n_eq };
      rates[s].eta_N_pair[0] = cases[c].q;
      require_error(
            ghl_m1_validate_neutrino_rates(&rates[s], NULL), ghl_success,
            "range-test pair rates", 790 + (int)c);
      input[s] = (ghl_m1_neutrino_state){ .E = 1.0, .N = cases[c].N };
    }
    require_error(
          ghl_m1_solve_neutrino_pair_source_update(
                params, nu, &metric, &prims, rates, input, input, cases[c].dt, 1.0,
                output, exchange, diagnostics, nd),
          ghl_error_m1_implicit_terminal_fallback, "quadratic range rejection",
          790 + (int)c);
    for(int s = 0; s < 2; ++s) {
      require_condition(
            memcmp(&output[s], &input[s], sizeof(input[s])) == 0,
            "range rejection changed state", 790 + (int)c);
      require_condition(
            diagnostics[s].terminal_no_update
                  && diagnostics[s].path
                           == ghl_m1_neutrino_source_path_terminal_no_update
                  && nd[s].source_terminal_fallbacks == 1,
            "range rejection lost terminal diagnostics", 790 + (int)c);
      require_condition(
            exchange[s].dN_rad_total == 0.0 && exchange[s].dL_rad_cc == 0.0
                  && exchange[s].dE_rad == 0.0 && exchange[s].dTau_matter == 0.0
                  && exchange[s].dYe_matter == 0.0,
            "range rejection published exchange", 790 + (int)c);
      for(int i = 0; i < 3; ++i) {
        require_condition(
              exchange[s].dF_rad[i] == 0.0 && exchange[s].dS_matter[i] == 0.0,
              "range rejection published momentum", 790 + (int)c);
      }
    }
  }
}

static void test_projected_number_range_rejection(void) {
  const ghl_m1_neutrino_parameters params = { .Gamma_N_floor = DBL_MIN * 0.5 };
  const ghl_m1_neutrino_rates rates
        = { .species = ghl_m1_neutrino_nux, .mean_energy = DBL_MIN };
  const ghl_m1_neutrino_state input = { .N = 1.0 };
  const ghl_m1_neutrino_current current = { .Gamma_N = 1.0, .J = DBL_MAX };
  double N = -7.0;
  bool projected = false;
  require_error(
        ghl_m1_validate_neutrino_rates(&rates, NULL), ghl_success,
        "projection overflow rates", 794);
  require_error(
        ghl_m1_neutrino_update_endpoint_number_with_policy(
              &params, &rates, 1.0, 1.0, 0.0, &input, &current, &N, &projected),
        ghl_error_m1_invalid_state, "projected quotient overflow", 794);
  require_condition(N == -7.0 && !projected, "failed projection modified outputs", 794);
  const struct {
    double gamma, J, mean, expected;
  } representable[] = { { DBL_MAX, 2.0, DBL_MAX, 2.0 },
                        { DBL_MIN, DBL_MIN, DBL_MIN, DBL_MIN },
                        { 1.0, 0.0, 1.0, 0.0 } };
  for(size_t c = 0; c < sizeof(representable) / sizeof(representable[0]); ++c) {
    const ghl_m1_neutrino_rates r
          = { .species = ghl_m1_neutrino_nux, .mean_energy = representable[c].mean };
    const ghl_m1_neutrino_current endpoint
          = { .Gamma_N = representable[c].gamma, .J = representable[c].J };
    N = -7.0;
    projected = false;
    require_error(
          ghl_m1_validate_neutrino_rates(&r, NULL), ghl_success,
          "representable projection rates", 795 + (int)c);
    require_error(
          ghl_m1_neutrino_update_endpoint_number_with_policy(
                &params, &r, 1.0, 1.0, 0.0, &input, &endpoint, &N, &projected),
          ghl_success, "representable scaled projection", 795 + (int)c);
    require_condition(
          N == representable[c].expected && projected,
          "scaled projection lost finite endpoint", 795 + (int)c);
  }
}

#endif
