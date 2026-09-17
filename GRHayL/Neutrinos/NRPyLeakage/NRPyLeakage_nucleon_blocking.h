#ifndef NRPYLEAKAGE_NUCLEON_BLOCKING_H_
#define NRPYLEAKAGE_NUCLEON_BLOCKING_H_

#include <float.h>
#include <math.h>

#include "ghl_nrpyleakage.h"

/**
 * Free-nucleon blocking follows the density-derived kinetic-degeneracy
 * construction in ILEAS, Appendix B, Eqs. (69)--(71):
 * A. Ardevol-Pulpillo et al., MNRAS 485 (2019), 4754--4787,
 * doi:10.1093/mnras/stz613. The same-energy transition overlap below is
 * evaluated through the density-normalization identity. Charged-current grey
 * moments use the polynomial threshold expansion of ILEAS Appendix B,
 * Eqs. (78)--(88), with the reaction shift chosen so that the spectral parent
 * kernels obey the Kirchhoff pairing in Appendix C, Eqs. (100)--(109). This is
 * an algebraic approximation: it adds no quadrature, root solve, or table
 * lookup to the leakage hot path.
 * Exact subtraction of nearby positive floating-point populations uses
 * Sterbenz's lemma: P. H. Sterbenz, Floating-Point Computation,
 * Prentice-Hall, 1974, Sec. 4.3.
 *
 * The scalar fdm1h and ifd1h fits and coefficients below are ported
 * from Scott J. Maddox's FDINT implementation:
 * https://github.com/scott-maddox/fdint/blob/master/fdint/_fdint.pyx
 * They implement T. Fukushima's minimax approximations:
 * doi:10.1016/j.amc.2015.03.009 and doi:10.1016/j.amc.2015.03.015.
 * GRHayL's density normalization and transition-overlap evaluation are
 * adaptations in this file. FDINT does not provide or qualify them, and they
 * are not copied formulas from ILEAS.
 */

/*
 * Copyright (c) 2015, Scott J Maddox.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are
 * permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of
 * conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list
 * of conditions and the following disclaimer in the documentation and/or other materials
 * provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be
 * used to endorse or promote products derived from this software without specific prior
 * written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
 * WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 */

static inline double nrpyl_fdm1h_lt_m2(const double phi) {
  double exp_phi = exp(phi);
  double t = exp_phi * 7.38905609893065023;
  return exp_phi
         * (1.77245385090551603
            - exp_phi
                    * (40641.4537510284430
                       + t
                               * (9395.7080940846442
                                  + t
                                          * (649.96168315267301
                                             + t
                                                     * (12.7972295804758967
                                                        + t * 0.00153864350767585460))))
                    / (32427.1884765292940
                       + t
                               * (11079.9205661274782
                                  + t
                                          * (1322.96627001478859
                                             + t * (63.738361029333467 + t)))));
}

static inline double nrpyl_fdm1h_m2_to_0(const double phi) {
  double s = -0.5 * phi;
  double t = 1.0 - s;
  return (272.770092131932696
          + t
                  * (30.8845653844682850
                     + t
                             * (-6.43537632380366113
                                + t
                                        * (14.8747473098217879
                                           + t
                                                   * (4.86928862842142635
                                                      + t
                                                              * (-1.53265834550673654
                                                                 + t
                                                                         * (-1.02698898315597491
                                                                            + t
                                                                                    * (-0.177686820928605932
                                                                                       - t * 0.00377141325509246441))))))))
         / (293.075378187667857
            + s
                    * (305.818162686270816
                       + s
                               * (299.962395449297620
                                  + s
                                          * (207.640834087494249
                                             + s
                                                     * (92.0384803181851755
                                                        + s
                                                                * (37.0164914112791209
                                                                   + s
                                                                           * (7.88500950271420583
                                                                              + s)))))));
}

static inline double nrpyl_fdm1h_0_to_2(const double phi) {
  double t = 0.5 * phi;
  return (3531.50360568243046
          + t
                  * (6077.5339658420037
                     + t
                             * (6199.7700433981326
                                + t
                                        * (4412.78701919567594
                                           + t
                                                   * (2252.27343092810898
                                                      + t
                                                              * (811.84098649224085
                                                                 + t
                                                                         * (191.836401053637121
                                                                            + t * 23.2881838959183802)))))))
         / (3293.83702584796268
            + t
                    * (1528.97474029789098
                       + t
                               * (2568.48562814986046
                                  + t
                                          * (925.64264653555825
                                             + t
                                                     * (574.23248354035988
                                                        + t
                                                                * (132.803859320667262
                                                                   + t
                                                                           * (29.8447166552102115
                                                                              + t)))))));
}

static inline double nrpyl_fdm1h_2_to_5(const double phi) {
  double t = 0.3333333333333333333 * (phi - 2.0);
  return (4060.70753404118265
          + t
                  * (10812.7291333052766
                     + t
                             * (13897.5649482242583
                                + t
                                        * (10628.4749852740029
                                           + t
                                                   * (5107.70670190679021
                                                      + t
                                                              * (1540.84330126003381
                                                                 + t
                                                                         * (284.452720112970331
                                                                            + t * 29.5214417358484151)))))))
         / (1564.58195612633534
            + t
                    * (2825.75172277850406
                       + t
                               * (3189.16066169981562
                                  + t
                                          * (1955.03979069032571
                                             + t
                                                     * (828.000333691814748
                                                        + t
                                                                * (181.498111089518376
                                                                   + t
                                                                           * (32.0352857794803750
                                                                              + t)))))));
}

static inline double nrpyl_fdm1h_5_to_10(const double phi) {
  double t = 0.2 * phi - 1.0;
  return (1198.41719029557508
          + t
                  * (3263.51454554908654
                     + t
                             * (3874.97588471376487
                                + t
                                        * (2623.13060317199813
                                           + t
                                                   * (1100.41355637121217
                                                      + t
                                                              * (267.469532490503605
                                                                 + t
                                                                         * (25.4207671812718340
                                                                            + t * 0.389887754234555773)))))))
         / (273.407957792556998
            + t
                    * (595.918318952058643
                       + t
                               * (605.202452261660849
                                  + t
                                          * (343.183302735619981
                                             + t
                                                     * (122.187622015695729
                                                        + t
                                                                * (20.9016359079855933
                                                                   + t))))));
}

static inline double nrpyl_fdm1h_10_to_20(const double phi) {
  double t = 0.1 * phi - 1.0;
  return (9446.00169435237637
          + t
                  * (36843.4448474028632
                     + t
                             * (63710.1115419926191
                                + t
                                        * (62985.2197361074768
                                           + t
                                                   * (37634.5231395700921
                                                      + t
                                                              * (12810.9898627807754
                                                                 + t
                                                                         * (1981.56896138920963
                                                                            + t * 81.4930171897667580)))))))
         / (1500.04697810133666
            + t
                    * (5086.91381052794059
                       + t
                               * (7730.01593747621895
                                  + t
                                          * (6640.83376239360596
                                             + t
                                                     * (3338.99590300826393
                                                        + t
                                                                * (860.499043886802984
                                                                   + t
                                                                           * (78.8565824186926692
                                                                              + t)))))));
}

static inline double nrpyl_fdm1h_20_to_40(const double phi) {
  double t = 0.05 * phi - 1.0;
  return (22977.9657855367223
          + t
                  * (123416.616813887781
                     + t
                             * (261153.765172355107
                                + t
                                        * (274618.894514095795
                                           + t
                                                   * (149710.718389924860
                                                      + t
                                                              * (40129.3371700184546
                                                                 + t
                                                                         * (4470.46495881415076
                                                                            + t * 132.684346831002976)))))))
         / (2571.68842525335676
            + t
                    * (12521.4982290775358
                       + t
                               * (23268.1574325055341
                                  + t
                                          * (20477.2320119758141
                                             + t
                                                     * (8726.52577962268114
                                                        + t
                                                                * (1647.42896896769909
                                                                   + t
                                                                           * (106.475275142076623
                                                                              + t)))))));
}

static inline double nrpyl_fdm1h_gt_40(const double phi) {
  double factor = 2.0;
  double w = 1.0 / (phi * phi);
  double t = 1600.0 * w;
  return sqrt(phi) * factor
         * (1.0
            - w
                    * (0.411233516712009968
                       + t
                               * (0.00110980410034088951
                                  + t
                                          * (0.0000113689298990173683
                                             + t
                                                     * (2.56931790679436797e-7
                                                        + t
                                                                * (9.97897786755446178e-9
                                                                   + t * 8.67667698791108582e-10))))));
}

typedef struct {
  double eta;
} nrpyl_nucleon_population;

/**
 * @brief Dimensionless number and energy moments for one beta channel.
 *
 * The caller supplies the dimensional factors \f$T^5\f$ for emission and
 * \f$T^2\f$ for absorption. The energy member includes one additional factor
 * of \f$E/T\f$; callers multiply emission energy moments by \f$T\f$.
 */
typedef struct {
  double number;
  double energy;
} nrpyl_beta_moments;

static inline double nrpyl_fdm1h(const double phi) {
  if(phi < -2.0) {
    return nrpyl_fdm1h_lt_m2(phi);
  }
  if(phi < 0.0) {
    return nrpyl_fdm1h_m2_to_0(phi);
  }
  if(phi < 2.0) {
    return nrpyl_fdm1h_0_to_2(phi);
  }
  if(phi < 5.0) {
    return nrpyl_fdm1h_2_to_5(phi);
  }
  if(phi < 10.0) {
    return nrpyl_fdm1h_5_to_10(phi);
  }
  if(phi < 20.0) {
    return nrpyl_fdm1h_10_to_20(phi);
  }
  if(phi < 40.0) {
    return nrpyl_fdm1h_20_to_40(phi);
  }
  return nrpyl_fdm1h_gt_40(phi);
}

static inline double nrpyl_ifd1h(const double nu) {
  double v, s, t, w, y, z;
  if(nu < 1.17683303804380831) {
    t = nu * 0.849738210666018375;
    z = t
        * (156377.8333056294
           + t
                   * (48177.5705898287
                      + t
                              * (5847.07218383812
                                 + t * (335.3978079672194 + t * 7.84411868029912))))
        / (117762.02905535089
           + t
                   * (-19007.26938370368
                      + t * (1376.2936928453140 + t * (-54.11372698481717 + t))));
    y = log(z);
  }
  else if(nu < 3.82993088157949761) {
    t = 0.376917874490198033 * nu - 0.443569407329314587;
    y = (489.140447310410217
         + t
                 * (5335.07269317261966
                    + t
                            * (20169.0736140442509
                               + t
                                       * (35247.8115595510907
                                          + t
                                                  * (30462.3668614714761
                                                     + t
                                                             * (12567.9032426128967
                                                                + t
                                                                        * (2131.86789357398657
                                                                           + t * 93.6520172085419439)))))))
        / (656.826207643060606
           + t
                   * (4274.82831051941605
                      + t
                              * (10555.7581310151498
                                 + t
                                         * (12341.8742094611883
                                            + t
                                                    * (6949.18854413197094
                                                       + t
                                                               * (1692.19650634194002
                                                                  + t
                                                                          * (129.221772991589751
                                                                             + t)))))));
  }
  else if(nu < 13.3854493161866553) {
    t = 0.104651569335924949 * nu - 0.400808277205416960;
    y = (1019.84886406642351
         + t
                 * (9440.18255003922075
                    + t
                            * (33947.6616363762463
                               + t
                                       * (60256.7280980542786
                                          + t
                                                  * (55243.0045063055787
                                                     + t
                                                             * (24769.8354802210838
                                                                + t
                                                                        * (4511.77288617668292
                                                                           + t * 211.432806336150141)))))))
        / (350.502070353586442
           + t
                   * (2531.06296201234050
                      + t
                              * (6939.09850659439245
                                 + t
                                         * (9005.40197972396592
                                            + t
                                                    * (5606.73612994134056
                                                       + t
                                                               * (1488.76634564005075
                                                                  + t
                                                                          * (121.537028889412581
                                                                             + t)))))));
  }
  else if(nu < 53.2408277860982205) {
    t = 0.0250907164450825724 * nu - 0.335850513282463787;
    y = (11885.8779398399498
         + t
                 * (113220.250825178799
                    + t
                            * (408524.373881197840
                               + t
                                       * (695674.357483475952
                                          + t
                                                  * (569389.917088505552
                                                     + t
                                                             * (206433.082013681440
                                                                + t
                                                                        * (27307.2535671974100
                                                                           + t * 824.430826794730740)))))))
        / (1634.40491220861182
           + t
                   * (12218.1158551884025
                      + t
                              * (32911.7869957793233
                                 + t
                                         * (38934.6963039399331
                                            + t
                                                    * (20038.8358438225823
                                                       + t
                                                               * (3949.48380897796954
                                                                  + t
                                                                          * (215.607404890995706
                                                                             + t)))))));
  }
  else if(nu < 188.411871723022843) {
    t = 0.00739803415638806339 * nu - 0.393877462475929313;
    y = (11730.7011190435638
         + t
                 * (99421.7455796633651
                    + t
                            * (327706.968910706902
                               + t
                                       * (530425.668016563224
                                          + t
                                                  * (438631.900516555072
                                                     + t
                                                             * (175322.855662315845
                                                                + t
                                                                        * (28701.9605988813884
                                                                           + t * 1258.20914464286403)))))))
        / (634.080470383026173
           + t
                   * (4295.63159860265838
                      + t
                              * (10868.5260668911946
                                 + t
                                         * (12781.6871997977069
                                            + t
                                                    * (7093.80732100760563
                                                       + t
                                                               * (1675.06417056300026
                                                                  + t
                                                                          * (125.750901817759662
                                                                             + t)))))));
  }
  else {
    v = pow(nu, -4.0 / 3.0);
    s = 1080.13412050984017 * v;
    t = 1. - s;
    w = (1.12813495144821933e7
         + t * (420368.911157160874 + t * (1689.69475714536117 + t)))
        / (s
           * (6088.08350831295857
              + t * (221.445236759466761 + t * 0.718216708695397737)));
    y = sqrt(w);
  }
  return y;
}

static inline ghl_error_codes_t nrpyl_compute_population(
      const double log_y,
      nrpyl_nucleon_population *restrict population) {
  // Below DBL_MIN, the first fugacity term determines eta to much better than
  // double precision; keeping the normalization logarithmic avoids underflow.
  if(log_y < log(DBL_MIN)) {
    population->eta = log_y - log(0.886226925452758014);
    return robust_isfinite(population->eta) ? ghl_success : ghl_error_nrpyleakage_blocking;
  }

  if(log_y > log(DBL_MAX)) {
    return ghl_error_nrpyleakage_blocking;
  }

  const double y = exp(log_y);
  population->eta = nrpyl_ifd1h(y);
  if(!robust_isfinite(population->eta)) {
    return ghl_error_nrpyleakage_blocking;
  }
  return ghl_success;
}

static inline double nrpyl_log_fd1h_derivative(const double eta) {
  const double derivative = 0.5 * nrpyl_fdm1h(eta);
  if(derivative > 0.0 && robust_isfinite(derivative)) {
    return log(derivative);
  }

  // In the only underflowing branch, F'_{1/2} has the same leading
  // fugacity term as F_{1/2}.
  return eta + log(0.886226925452758014);
}

/**
 * @brief Evaluate the threshold-shifted Fermi moments used by beta reactions.
 *
 * For \f$a,b\geq0\f$ with \f$ab=0\f$, this evaluates
 * \f$\Phi_k=\int_0^\infty (x+a)^{2+k}(x+b)^2
 * [1+\exp(x-\eta)]^{-1}dx\f$ for \f$k=0,1\f$ by expanding the polynomial into
 * the existing complete Fermi integrals. This is the algebraic grey
 * approximation adapted from ILEAS Appendix B, Eqs. (78)--(88).
 *
 * @param[in] a Neutrino-energy threshold in units of temperature.
 * @param[in] b Charged-lepton-energy shift in units of temperature.
 * @param[in] eta Degeneracy parameter after the threshold change of variable.
 * @param[out] moments Number and energy phase-space moments.
 * @return `ghl_success` on success, otherwise a Fermi-integral error.
 */
static inline ghl_error_codes_t nrpyl_compute_shifted_fermi_moments(
      const double a,
      const double b,
      const double eta,
      nrpyl_beta_moments *restrict moments) {
  double F2, F3, F4, F5;
  ghl_error_codes_t error = NRPyLeakage_Fermi_Dirac_integrals(2, eta, &F2);
  if(error != ghl_success) {
    return error;
  }
  error = NRPyLeakage_Fermi_Dirac_integrals(3, eta, &F3);
  if(error != ghl_success) {
    return error;
  }
  error = NRPyLeakage_Fermi_Dirac_integrals(4, eta, &F4);
  if(error != ghl_success) {
    return error;
  }
  error = NRPyLeakage_Fermi_Dirac_integrals(5, eta, &F5);
  if(error != ghl_success) {
    return error;
  }

  /*
   * In a strongly blocked channel all complete moments can underflow
   * together. Their common Boltzmann factor then lies below the dynamic range
   * of double precision, so both raw phase-space moments have the
   * representable zero limit. Emission uses that limit directly; normalized
   * absorption cancels the common Boltzmann factor analytically below.
   */
  if(F2 == 0.0 && F3 == 0.0 && F4 == 0.0 && F5 == 0.0) {
    moments->number = 0.0;
    moments->energy = 0.0;
    return ghl_success;
  }

  moments->number = F4 + 2.0 * (a + b) * F3 + (a * a + b * b) * F2;
  moments->energy
        = F5 + (3.0 * a + 2.0 * b) * F4 + (3.0 * a * a + b * b) * F3 + a * a * a * F2;
  if(!robust_isfinite(moments->number) || !robust_isfinite(moments->energy) || moments->number < 0.0
     || moments->energy < 0.0) {
    return ghl_error_nrpyleakage_blocking;
  }
  return ghl_success;
}

/**
 * @brief Return the vacancy \f$1-f(E)\f$ without overflowing the exponential.
 *
 * @param[in] mu_minus_E_over_T The exponent \f$(\mu-E)/T\f$.
 * @return The Fermi vacancy in the closed interval \f$[0,1]\f$.
 */
static inline double nrpyl_fermi_vacancy(const double mu_minus_E_over_T) {
  if(mu_minus_E_over_T > 0.0) {
    const double inverse_exponential = exp(-mu_minus_E_over_T);
    return inverse_exponential / (1.0 + inverse_exponential);
  }
  return 1.0 / (1.0 + exp(mu_minus_E_over_T));
}

/**
 * @brief Reconcile kinetic nucleon occupations with the EOS energy convention.
 *
 * The shift \f$q=\hat\mu-T(\eta_n-\eta_p)\f$ makes the paired spectral beta
 * kernels obey detailed balance with the EOS chemical-potential difference.
 * No physical clamp is applied: the current EOS API supplies no authoritative
 * mean-field bound. Large finite \f$|q/T|\f$ can remain poorly conditioned;
 * downstream helpers reject nonfinite inputs and results.
 *
 * @param[in] T Temperature in MeV.
 * @param[in] muhat EOS neutron-minus-proton chemical potential in MeV.
 * @param[in] eta_n_minus_eta_p Kinetic nucleon degeneracy difference.
 * @return Reaction-energy shift in MeV.
 */
static inline double nrpyl_compute_reaction_shift(
      const double T,
      const double muhat,
      const double eta_n_minus_eta_p) {
  return muhat - T * eta_n_minus_eta_p;
}

/**
 * @brief Compute algebraic grey emission moments for one beta channel.
 *
 * The sign is \f$+1\f$ for electron-neutrino emission and \f$-1\f$ for
 * electron-antineutrino emission. The reaction shift is
 * \f$q=\hat\mu-T(\eta_n-\eta_p)\f$; this reconciles density-derived nucleon
 * occupations with the EOS chemical difference. The underlying spectral
 * kernels then satisfy detailed balance, while the returned grey neutrino
 * vacancy is evaluated at the mean emitted energy.
 *
 * @param[in] T Temperature in MeV.
 * @param[in] mu_e Electron chemical potential in MeV.
 * @param[in] eta_nu Neutrino degeneracy parameter.
 * @param[in] sign `+1` for nue or `-1` for anue.
 * @param[in] q Reaction-energy shift in MeV.
 * @param[out] moments Dimensionless blocked emission moments.
 * @return `ghl_success` on success, otherwise a Fermi-integral or input error.
 */
static inline ghl_error_codes_t nrpyl_compute_beta_emission_moments(
      const double T,
      const double mu_e,
      const double eta_nu,
      const int sign,
      const double q,
      nrpyl_beta_moments *restrict moments) {
  if(!(T > 0.0) || !robust_isfinite(T) || !robust_isfinite(mu_e)
     || !robust_isfinite(eta_nu) || !robust_isfinite(q)
     || (sign != 1 && sign != -1)) {
    return ghl_error_nrpyleakage_blocking;
  }

  const double shift = sign * q / T;
  const double a = fmax(0.0, -shift);
  const double b = fmax(0.0, shift);
  ghl_error_codes_t error
        = nrpyl_compute_shifted_fermi_moments(a, b, sign * mu_e / T - b, moments);
  if(error != ghl_success) {
    return error;
  }

  if(moments->number == 0.0) {
    return ghl_success;
  }

  const double mean_energy_over_T = moments->energy / moments->number;
  const double vacancy = nrpyl_fermi_vacancy(eta_nu - mean_energy_over_T);
  moments->number *= vacancy;
  moments->energy *= vacancy;
  return ghl_success;
}

/**
 * @brief Compute algebraic grey absorption moments for one beta channel.
 *
 * This uses the same reaction shift and threshold as the paired emission
 * kernel. Ordinary absorption includes final-state charged-lepton blocking;
 * stimulated neutrino absorption is not included in an opacity.
 *
 * @param[in] T Temperature in MeV.
 * @param[in] mu_e Electron chemical potential in MeV.
 * @param[in] eta_nu Neutrino degeneracy parameter.
 * @param[in] sign `+1` for nue or `-1` for anue.
 * @param[in] q Reaction-energy shift in MeV.
 * @param[out] moments Dimensionless blocked opacity moments.
 * @return `ghl_success` on success, otherwise a Fermi-integral or input error.
 */
static inline ghl_error_codes_t nrpyl_compute_beta_absorption_moments(
      const double T,
      const double mu_e,
      const double eta_nu,
      const int sign,
      const double q,
      nrpyl_beta_moments *restrict moments) {
  if(!(T > 0.0) || !robust_isfinite(T) || !robust_isfinite(mu_e)
     || !robust_isfinite(eta_nu) || !robust_isfinite(q)
     || (sign != 1 && sign != -1)) {
    return ghl_error_nrpyleakage_blocking;
  }

  const double shift = sign * q / T;
  const double a = fmax(0.0, -shift);
  const double b = fmax(0.0, shift);
  ghl_error_codes_t error
        = nrpyl_compute_shifted_fermi_moments(a, b, eta_nu - a, moments);
  if(error != ghl_success) {
    return error;
  }

  double F2, F3;
  error = NRPyLeakage_Fermi_Dirac_integrals(2, eta_nu, &F2);
  if(error != ghl_success) {
    return error;
  }
  error = NRPyLeakage_Fermi_Dirac_integrals(3, eta_nu, &F3);
  if(error != ghl_success) {
    return error;
  }

  if(moments->number == 0.0) {
    /*
     * Reconstruct a shifted numerator below double's range from its Boltzmann
     * limit. Retain a representable unshifted F2/F3 denominator; if it also
     * underflows, cancel the common exp(eta_nu) factor analytically. This
     * preserves every normalized opacity that remains representable.
     */
    const double number_polynomial = 12.0 + 6.0 * (a + b) + a * a + b * b;
    const double energy_polynomial
          = 20.0 + 12.0 * a + 8.0 * b + 3.0 * a * a + b * b + a * a * a / 3.0;
    if(!(number_polynomial > 0.0) || !(energy_polynomial > 0.0)
       || !robust_isfinite(number_polynomial) || !robust_isfinite(energy_polynomial)) {
      return ghl_error_nrpyleakage_blocking;
    }

    const double mean_energy_over_T = 3.0 * energy_polynomial / number_polynomial;
    const double lepton_energy_over_T = mean_energy_over_T + shift;
    const double vacancy = nrpyl_fermi_vacancy(sign * mu_e / T - lepton_energy_over_T);
    const double log_number_ratio
          = F2 > 0.0 ? log(2.0) + eta_nu - a + log(number_polynomial) - log(F2)
                     : log(number_polynomial) - a;
    const double log_energy_ratio
          = F3 > 0.0 ? log(6.0) + eta_nu - a + log(energy_polynomial) - log(F3)
                     : log(energy_polynomial) - a;
    moments->number = exp(log_number_ratio) * vacancy;
    moments->energy = exp(log_energy_ratio) * vacancy;
    if(!robust_isfinite(moments->number) || !robust_isfinite(moments->energy)
       || moments->number < 0.0
       || moments->energy < 0.0) {
      return ghl_error_nrpyleakage_blocking;
    }
    return ghl_success;
  }

  const double mean_energy_over_T = moments->energy / moments->number;
  const double lepton_energy_over_T = mean_energy_over_T + shift;
  const double vacancy = nrpyl_fermi_vacancy(sign * mu_e / T - lepton_energy_over_T);
  /* Divide the paired Boltzmann factors before applying the O(1) vacancy.
   * Forming vacancy/Fk first can overflow when both the shifted moment and Fk
   * are representable subnormals, even though their physical ratio is finite.
   */
  moments->number = (moments->number / F2) * vacancy;
  moments->energy = (moments->energy / F3) * vacancy;
  if(!robust_isfinite(moments->number) || !robust_isfinite(moments->energy)) {
    return ghl_error_nrpyleakage_blocking;
  }
  return ghl_success;
}

/**
 * @brief Compute scattering and charged-current free-nucleon populations.
 *
 * Density-derived kinetic occupations make the result independent of the
 * EOS chemical-energy zero. Stable normalization identities evaluate the
 * transition overlaps without a quotient singularity at equal populations.
 *
 * @param[in] rho_cgs Rest-mass density in g cm\f$^{-3}\f$.
 * @param[in] T Temperature in MeV.
 * @param[in] X_n Free-neutron mass fraction.
 * @param[in] X_p Free-proton mass fraction.
 * @param[out] B_n Effective neutron population for neutral-current scattering.
 * @param[out] B_p Effective proton population for neutral-current scattering.
 * @param[out] Y_np Effective neutron-to-proton transition population.
 * @param[out] Y_pn Effective proton-to-neutron transition population.
 * @param[out] eta_n_minus_eta_p Kinetic nucleon degeneracy difference. This
 * is zero and must not be used at an exact single-species endpoint, where
 * the charged-current products have their analytic zero limit.
 * @return `ghl_success` on success, otherwise a blocking-evaluator error.
 */
static inline ghl_error_codes_t NRPyLeakage_compute_nucleon_blocking(
      const double rho_cgs,
      const double T,
      const double X_n,
      const double X_p,
      double *restrict B_n,
      double *restrict B_p,
      double *restrict Y_np,
      double *restrict Y_pn,
      double *restrict eta_n_minus_eta_p) {
  /*
   * NRPyEOS evaluates an eight-corner trilinear polynomial in coefficient
   * form.  Expanding that expression gives 27 signed corner contributions,
   * each bounded by one for in-cell coordinates and fractions in [0,1].  Its
   * path has fewer than 64 rounded additions and multiplications, so
   * 27*gamma_64 is a conservative absolute forward-error bound.  Normalize
   * only that table/interpolation noise; larger excursions remain errors.
   * This includes the -8.24e-17 neutron fraction present in SLy4.
   */
  const double gamma_64 = 64.0 * DBL_EPSILON / (1.0 - 64.0 * DBL_EPSILON);
  const double fraction_roundoff = 27.0 * gamma_64;
  if(!robust_isfinite(rho_cgs) || !(rho_cgs > 0.0) || !robust_isfinite(T) || !(T > 0.0)
     || !robust_isfinite(X_n) || X_n < -fraction_roundoff || X_n > 1.0 + fraction_roundoff
     || !robust_isfinite(X_p) || X_p < -fraction_roundoff || X_p > 1.0 + fraction_roundoff) {
    return ghl_error_nrpyleakage_blocking;
  }

  const double physical_X_n = fmin(1.0, fmax(0.0, X_n));
  const double physical_X_p = fmin(1.0, fmax(0.0, X_p));

  *B_n = *B_p = *Y_np = *Y_pn = *eta_n_minus_eta_p = 0.0;
  if(physical_X_n == 0.0 && physical_X_p == 0.0) {
    return ghl_success;
  }

  // Common free-nucleon kinetic mass used by the ILEAS blocking model.
  const double nucleon_mass_MeV = 938.91872;
  const double log_C
        = log(4.0 * M_PI / NRPyLeakage_hc3) + 1.5 * log(2.0 * nucleon_mass_MeV * T);
  const double log_number_density = log(NRPyLeakage_N_A) + log(rho_cgs);
  if(!robust_isfinite(log_C) || !robust_isfinite(log_number_density)) {
    return ghl_error_nrpyleakage_blocking;
  }

  nrpyl_nucleon_population neutron = { 0.0 };
  nrpyl_nucleon_population proton = { 0.0 };
  ghl_error_codes_t error;
  if(physical_X_n > 0.0) {
    error = nrpyl_compute_population(
          log_number_density + log(physical_X_n) - log_C, &neutron);
    if(error != ghl_success) {
      return error;
    }
    *B_n = physical_X_n / (1.0 + (2.0 / 3.0) * fmax(neutron.eta, 0.0));
  }
  if(physical_X_p > 0.0) {
    error = nrpyl_compute_population(
          log_number_density + log(physical_X_p) - log_C, &proton);
    if(error != ghl_success) {
      return error;
    }
    *B_p = physical_X_p / (1.0 + (2.0 / 3.0) * fmax(proton.eta, 0.0));
  }

  /*
   * As one fraction tends to zero, its kinetic degeneracy tends to -infinity.
   * The occupied-to-empty overlap tends to the occupied fraction and the
   * reverse overlap tends to zero. The associated shifted beta moments and
   * reverse overlap make every charged-current product vanish exponentially,
   * so callers apply that analytic zero limit without forming an infinite
   * reaction shift. eta_n_minus_eta_p is consequently unused here.
   */
  if(physical_X_n == 0.0) {
    *Y_pn = physical_X_p;
    return ghl_success;
  }
  if(physical_X_p == 0.0) {
    *Y_np = physical_X_n;
    return ghl_success;
  }

  // Both populations have the same density normalization, so ordering their
  // mass fractions also orders F_{1/2}(eta) and eta.
  const bool neutron_is_high = physical_X_n >= physical_X_p;
  const nrpyl_nucleon_population high = neutron_is_high ? neutron : proton;
  const nrpyl_nucleon_population low = neutron_is_high ? proton : neutron;
  const double X_high = neutron_is_high ? physical_X_n : physical_X_p;
  const double X_low = neutron_is_high ? physical_X_p : physical_X_n;
  double a = high.eta - low.eta;
  if(!robust_isfinite(a)) {
    return ghl_error_nrpyleakage_blocking;
  }

  const double population_difference = X_high - X_low;
  if(population_difference <= sqrt(DBL_EPSILON) * X_high) {
    /*
     * Independent inverse-fit rounding dominates a direct subtraction as the
     * populations converge. A midpoint derivative gives the centered inverse
     * difference with second-order truncation error. The sqrt(DBL_EPSILON)
     * crossover follows by balancing that truncation against subtraction
     * roundoff; it is a machine-precision criterion, not a physical tuning
     * parameter.
     */
    nrpyl_nucleon_population midpoint;
    const double X_midpoint = 0.5 * (X_high + X_low);
    ghl_error_codes_t error = nrpyl_compute_population(
          log_number_density + log(X_midpoint) - log_C, &midpoint);
    if(error != ghl_success) {
      return error;
    }
    const double log_midpoint_overlap
          = log_C - log_number_density + nrpyl_log_fd1h_derivative(midpoint.eta);
    const double midpoint_overlap = exp(log_midpoint_overlap);
    if(!(midpoint_overlap > 0.0) || !robust_isfinite(midpoint_overlap)) {
      return ghl_error_nrpyleakage_blocking;
    }
    a = population_difference / midpoint_overlap;
  }
  else if(a < 0.0) {
    return ghl_error_nrpyleakage_blocking;
  }
  *eta_n_minus_eta_p = neutron_is_high ? a : -a;

  double Y_high_to_low, Y_low_to_high;
  if(a == 0.0) {
    // lim_{a->0} C*J/n_b = C*F'_{1/2}(eta)/n_b.
    const double log_overlap
          = log_C - log_number_density + nrpyl_log_fd1h_derivative(high.eta);
    if(robust_isnan(log_overlap)
       || (!robust_isfinite(log_overlap) && !signbit(log_overlap))) {
      return ghl_error_nrpyleakage_blocking;
    }
    Y_high_to_low = exp(log_overlap);
    Y_low_to_high = Y_high_to_low;
  }
  else {
    /*
     * The common density normalization gives
     *   C(F_high-F_low) = X_high-X_low.
     * Therefore the exact same-energy occupation overlap is
     *   Y_high->low = (X_high-X_low)/(1-exp(-a)),
     *   Y_low->high = exp(-a) Y_high->low.
     * Subtraction is exact for nearby positive doubles by Sterbenz's lemma,
     * while expm1 retains the close-degeneracy digits. No second Fermi
     * integral evaluation or user-chosen threshold is needed.
     */
    const double overlap = population_difference / (-expm1(-a));
    const double reverse_overlap = overlap * exp(-a);
    if(!(overlap >= 0.0) || !robust_isfinite(overlap) || !(reverse_overlap >= 0.0)
       || !robust_isfinite(reverse_overlap)) {
      return ghl_error_nrpyleakage_blocking;
    }
    Y_high_to_low = overlap;
    Y_low_to_high = reverse_overlap;
  }
  // Exact overlaps cannot exceed their initial populations. Project small fit
  // and rounding errors back onto that physical interval; independent
  // qualification, rather than this runtime projection, bounds fit error.
  Y_high_to_low = fmin(Y_high_to_low, X_high);
  Y_low_to_high = fmin(Y_low_to_high, X_low);
  if(neutron_is_high) {
    *Y_np = Y_high_to_low;
    *Y_pn = Y_low_to_high;
  }
  else {
    *Y_pn = Y_high_to_low;
    *Y_np = Y_low_to_high;
  }

  if(!robust_isfinite(*B_n) || !robust_isfinite(*B_p)
     || !robust_isfinite(*Y_np) || !robust_isfinite(*Y_pn)
     || *B_n < 0.0 || *B_n > physical_X_n || *B_p < 0.0 || *B_p > physical_X_p
     || *Y_np < 0.0 || *Y_np > physical_X_n || *Y_pn < 0.0 || *Y_pn > physical_X_p) {
    return ghl_error_nrpyleakage_blocking;
  }

  return ghl_success;
}

#endif // NRPYLEAKAGE_NUCLEON_BLOCKING_H_
