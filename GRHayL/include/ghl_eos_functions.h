#ifndef GHL_EOS_FUNCTIONS_H_
#define GHL_EOS_FUNCTIONS_H_

extern ghl_error_codes_t (*ghl_compute_h_and_cs2)(
      const ghl_eos_parameters *restrict eos,
      ghl_primitive_quantities *restrict prims,
      double *restrict h,
      double *restrict cs2);

extern int (*ghl_hybrid_find_polytropic_index)(
      const ghl_eos_parameters *restrict eos,
      const double rho_in);

extern void (*ghl_hybrid_get_K_and_Gamma)(
      const ghl_eos_parameters *restrict eos,
      const double rho_in,
      double *restrict K,
      double *restrict Gamma);

extern void (*ghl_hybrid_set_K_ppoly_and_eps_integ_consts)(
      ghl_eos_parameters *restrict eos);

extern void (*ghl_hybrid_compute_P_cold)(
      const ghl_eos_parameters *restrict eos,
      const double rho_in,
      double *restrict P_cold_ptr);

extern void (*ghl_hybrid_compute_P_cold_and_eps_cold)(
      const ghl_eos_parameters *restrict eos,
      const double rho_in,
      double *restrict P_cold_ptr,
      double *restrict eps_cold_ptr);

extern double (*ghl_hybrid_compute_epsilon)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double press);

extern double (*ghl_hybrid_compute_entropy_function)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double press);

  // Function prototypes
extern void (*ghl_tabulated_read_table_set_EOS_params)(
      const char *nuceos_table_name,
      ghl_eos_parameters *restrict eos);

extern void (*ghl_tabulated_free_memory)(
      ghl_eos_parameters *restrict eos);

extern ghl_error_codes_t (*ghl_tabulated_compute_P_from_T)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P);

extern ghl_error_codes_t (*ghl_tabulated_compute_eps_from_T)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict eps);

extern ghl_error_codes_t (*ghl_tabulated_compute_cs2_from_T)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict cs2);

extern ghl_error_codes_t (*ghl_tabulated_compute_P_eps_from_T)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P,
      double *restrict eps);

extern ghl_error_codes_t (*ghl_tabulated_compute_P_eps_S_from_T)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P,
      double *restrict eps,
      double *restrict S);

extern ghl_error_codes_t (*ghl_tabulated_compute_P_eps_cs2_from_T)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P,
      double *restrict eps,
      double *restrict cs2);

extern ghl_error_codes_t (*ghl_tabulated_compute_P_eps_S_cs2_from_T)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P,
      double *restrict eps,
      double *restrict S,
      double *restrict cs2);

extern ghl_error_codes_t (*ghl_tabulated_compute_P_eps_depsdT_from_T)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P,
      double *restrict eps,
      double *restrict depsdT);

extern ghl_error_codes_t (*ghl_tabulated_compute_P_eps_muhat_mue_mup_mun_from_T)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P,
      double *restrict eps,
      double *restrict muhat,
      double *restrict mu_e,
      double *restrict mu_p,
      double *restrict mu_n);

extern ghl_error_codes_t (*ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict muhat,
      double *restrict mu_e,
      double *restrict mu_p,
      double *restrict mu_n,
      double *restrict X_n,
      double *restrict X_p);

extern ghl_error_codes_t (*ghl_tabulated_compute_T_from_eps)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double eps,
      double *restrict T);

extern ghl_error_codes_t (*ghl_tabulated_compute_P_T_from_eps)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double eps,
      double *restrict P,
      double *restrict T);

extern ghl_error_codes_t (*ghl_tabulated_compute_P_cs2_T_from_eps)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double eps,
      double *restrict P,
      double *restrict cs2,
      double *restrict T);

extern ghl_error_codes_t (*ghl_tabulated_compute_eps_T_from_P)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double P,
      double *restrict eps,
      double *restrict T);

extern ghl_error_codes_t (*ghl_tabulated_compute_eps_cs2_T_from_P)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double P,
      double *restrict eps,
      double *restrict cs2,
      double *restrict T);

extern ghl_error_codes_t (*ghl_tabulated_compute_P_T_from_S)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double S,
      double *restrict P,
      double *restrict T);

extern ghl_error_codes_t (*ghl_tabulated_compute_P_S_depsdT_T_from_eps)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double eps,
      double *restrict P,
      double *restrict S,
      double *restrict depsdT,
      double *restrict T);

extern ghl_error_codes_t (*ghl_tabulated_compute_eps_S_T_from_P)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double P,
      double *restrict eps,
      double *restrict S,
      double *restrict T);

extern ghl_error_codes_t (*ghl_tabulated_compute_P_eps_T_from_S)(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double S,
      double *restrict P,
      double *restrict eps,
      double *restrict T);

extern int (*ghl_tabulated_get_index_T)(
      const ghl_eos_parameters *restrict eos,
      const double T);

extern void (*ghl_tabulated_compute_Ye_of_rho_beq_constant_T)(
      const double T,
      ghl_eos_parameters *restrict eos );

extern void (*ghl_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T)(
      const double T,
      ghl_eos_parameters *restrict eos );

extern double (*ghl_tabulated_compute_Ye_from_rho)(
      const ghl_eos_parameters *restrict eos,
      const double rho );

extern double (*ghl_tabulated_compute_P_from_rho)(
      const ghl_eos_parameters *restrict eos,
      const double rho );

extern double (*ghl_tabulated_compute_rho_from_P)(
      const ghl_eos_parameters *restrict eos,
      const double P );

extern double (*ghl_tabulated_compute_eps_from_rho)(
      const ghl_eos_parameters *restrict eos,
      const double rho );

extern void (*ghl_tabulated_free_beq_quantities)(
      ghl_eos_parameters *restrict eos);

extern void (*ghl_calculate_HLLE_fluxes_dirn0)(
      const ghl_primitive_quantities *restrict prims_r,
      const ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric_face,
      const double cmin,
      const double cmax,
      ghl_conservative_quantities *restrict cons_fluxes);

extern void (*ghl_calculate_HLLE_fluxes_dirn1)(
      const ghl_primitive_quantities *restrict prims_r,
      const ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric_face,
      const double cmin,
      const double cmax,
      ghl_conservative_quantities *restrict cons_fluxes);

extern void (*ghl_calculate_HLLE_fluxes_dirn2)(
      const ghl_primitive_quantities *restrict prims_r,
      const ghl_primitive_quantities *restrict prims_l,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric_face,
      const double cmin,
      const double cmax,
      ghl_conservative_quantities *restrict cons_fluxes);

extern void (*ghl_tabulated_enforce_bounds_rho_Ye_T)(
      const ghl_eos_parameters *restrict eos,
      double *restrict rho,
      double *restrict Y_e,
      double *restrict T );

extern void (*ghl_tabulated_enforce_bounds_rho_Ye_eps)(
      const ghl_eos_parameters *restrict eos,
      double *restrict rho,
      double *restrict Y_e,
      double *restrict eps );

extern void (*ghl_tabulated_enforce_bounds_rho_Ye_S)(
      const ghl_eos_parameters *restrict eos,
      double *restrict rho,
      double *restrict Y_e,
      double *restrict S );

extern void (*ghl_tabulated_enforce_bounds_rho_Ye_P)(
      const ghl_eos_parameters *restrict eos,
      double *restrict rho,
      double *restrict Y_e,
      double *restrict P );

#endif // GHL_EOS_FUNCTIONS_H_
