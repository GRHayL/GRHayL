#ifndef NRPYLEAKAGE_H_
#define NRPYLEAKAGE_H_

#ifdef MIN
#undef MIN
#endif

#ifdef MAX
#undef MAX
#endif

#define MAX(a,b) ( (a) > (b) ? (a) : (b) )
#define MIN(a,b) ( (a) < (b) ? (a) : (b) )

// "Primary" parameters
#define NRPyLeakage_enable_beta_nue (1)
#define NRPyLeakage_enable_beta_anue (1)
#define NRPyLeakage_enable_pair_nue_anue (1)
#define NRPyLeakage_enable_pair_nux_anux (1)
#define NRPyLeakage_enable_plasmon_nue_anue (1)
#define NRPyLeakage_enable_plasmon_nux_anux (1)
#define NRPyLeakage_enable_brems_nui_anui (1)
#define NRPyLeakage_Q_npmass (1.2935)
#define NRPyLeakage_gamma_0 (5.565e-2)
#define NRPyLeakage_sigma_0 (1.76e-44)
#define NRPyLeakage_alpha (1.25)
#define NRPyLeakage_C_A (0.5)
#define NRPyLeakage_sinthw2 (0.23)
#define NRPyLeakage_Brems_C1 (2.9988e7)
#define NRPyLeakage_Brems_C2 (6.5428e7)
#define NRPyLeakage_Brems_zeta (0.5)
#define NRPyLeakage_eta_nue_0 (0.0)
#define NRPyLeakage_eta_anue_0 (0.0)
#define NRPyLeakage_c_light (2.997924580000000e+10)
#define NRPyLeakage_N_A (6.022140760000000e+23)
#define NRPyLeakage_alpha_fs (7.297352569300000e-03)
#define NRPyLeakage_amu (1.660539066600000e-24)
#define NRPyLeakage_hc3 (1.905895198207216e-30)
#define NRPyLeakage_m_e_c2 (5.109989499961642e-01)
#define NRPyLeakage_units_geom_to_cgs_D (6.175828479261933e+17)
#define NRPyLeakage_units_cgs_to_geom_D (1.619215953548485e-18)
#define NRPyLeakage_units_cgs_to_geom_R (7.975433521479384e-24)
#define NRPyLeakage_units_cgs_to_geom_Q (1.421750164721599e-50)
#define NRPyLeakage_units_geom_to_cgs_M (1.988409870698051e+33)
#define NRPyLeakage_units_geom_to_cgs_L (1.476625038050125e+05)
#define NRPyLeakage_units_geom_to_cgs_T (4.925490947641267e-06)
#define NRPyLeakage_units_cgs_to_geom_M (5.029144215870041e-34)
#define NRPyLeakage_units_cgs_to_geom_L (6.772199944005382e-06)
#define NRPyLeakage_units_cgs_to_geom_T (2.030254467280836e+05)
#define NRPyLeakage_ZL_Q_npmass (1.293333)
#define NRPyLeakage_ZL_alpha (1.23)
#define NRPyLeakage_ZL_N_A (6.0221367e+23)
#define NRPyLeakage_ZL_alpha_fs (7.297352520505561e-03)
#define NRPyLeakage_ZL_amu (1.674927211e-24)
#define NRPyLeakage_ZL_hc3 (1.905893979207552e-30)
#define NRPyLeakage_ZL_m_e_c2 (5.1099891e-01)
// "Derived" parameters
#define NRPyLeakage_C_V (NRPyLeakage_C_A + 2*NRPyLeakage_sinthw2)
#define NRPyLeakage_beta (NRPyLeakage_c_light*NRPyLeakage_sigma_0/((NRPyLeakage_m_e_c2)*(NRPyLeakage_m_e_c2)))
#define NRPyLeakage_C1pC2_nue_anue (((-NRPyLeakage_C_A + NRPyLeakage_C_V)*(-NRPyLeakage_C_A + NRPyLeakage_C_V)) + ((NRPyLeakage_C_A + NRPyLeakage_C_V)*(NRPyLeakage_C_A + NRPyLeakage_C_V)))
#define NRPyLeakage_C1pC2_nux_anux (((-NRPyLeakage_C_A + NRPyLeakage_C_V)*(-NRPyLeakage_C_A + NRPyLeakage_C_V)) + ((NRPyLeakage_C_A + NRPyLeakage_C_V - 2)*(NRPyLeakage_C_A + NRPyLeakage_C_V - 2)))

#ifdef __cplusplus
extern "C" {
#endif
// Function prototypes
double NRPyLeakage_Fermi_Dirac_integrals(const int k, const double z);

void NRPyLeakage_compute_ghl_neutrino_opacities(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      const ghl_neutrino_optical_depths *restrict tau,
      ghl_neutrino_opacities *restrict kappa );

void NRPyLeakage_compute_ghl_neutrino_luminosities(
      const ghl_eos_parameters *restrict eos,
      const double alpha,
      const double gammaxx,
      const double gammaxy,
      const double gammaxz,
      const double gammayy,
      const double gammayz,
      const double gammazz,
      const double rho,
      const double Y_e,
      const double T,
      const double W,
      const ghl_neutrino_optical_depths *restrict tau,
      ghl_neutrino_luminosities *restrict lum );

void NRPyLeakage_compute_ghl_neutrino_opacities_and_GRMHD_source_terms(
      const ghl_eos_parameters *restrict eos,
      const double rho_b,
      const double Y_e,
      const double T,
      const ghl_neutrino_optical_depths *restrict tau,
      ghl_neutrino_opacities *restrict kappa,
      double *restrict R_source,
      double *restrict Q_source );

void NRPyLeakage_optical_depths_PathOfLeastResistance(
      const double *restrict dxx,
      const double *restrict stencil_gxx,
      const double *restrict stencil_gyy,
      const double *restrict stencil_gzz,
      const ghl_neutrino_opacities *restrict kappa_im1_j_k,
      const ghl_neutrino_opacities *restrict kappa_ip1_j_k,
      const ghl_neutrino_opacities *restrict kappa_i_jm1_k,
      const ghl_neutrino_opacities *restrict kappa_i_jp1_k,
      const ghl_neutrino_opacities *restrict kappa_i_j_km1,
      const ghl_neutrino_opacities *restrict kappa_i_j_kp1,
      const ghl_neutrino_opacities *restrict tau_im1_j_k,
      const ghl_neutrino_opacities *restrict tau_ip1_j_k,
      const ghl_neutrino_opacities *restrict tau_i_jm1_k,
      const ghl_neutrino_opacities *restrict tau_i_jp1_k,
      const ghl_neutrino_opacities *restrict tau_i_j_km1,
      const ghl_neutrino_opacities *restrict tau_i_j_kp1,
      const ghl_neutrino_opacities *restrict kappa_i_j_k,
      ghl_neutrino_optical_depths *restrict tau_i_j_k );

static inline int robust_isnan(double x) {
  unsigned long *pbits = (unsigned long *)&x;
  return( (*pbits & 0x7ff0000000000000UL) == 0x7ff0000000000000UL &&
          (*pbits & 0x000fffffffffffffUL) );
}

static inline int robust_isfinite(double x) {
  unsigned long *pbits = (unsigned long *)&x;
  return( !((*pbits & 0x7ff0000000000000UL) == 0x7ff0000000000000UL &&
            ((*pbits & 0x7ff0000000000000UL) || (*pbits & 0xfff0000000000000UL))) );
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif // NRPYLEAKAGE_H_
