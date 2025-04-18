#ifndef GHL_METRIC_HELPERS_H_
#define GHL_METRIC_HELPERS_H_

/**
 * @ingroup metric_contract
 * @brief Raise or lower a 4-vector
 *
 * @details
 * This function can compute both raising or lowering:
 *
 * - \f$ v^\nu = g^{\mu \nu} v_\nu \f$
 * - \f$ v_\nu = g_{\mu \nu} v^\nu \f$
 *
 * The exact operation being done is simply determined by the input
 * quantities. This function works for any input arrays, though it
 * is defined with ghl_ADM_aux_quantities::g4DD and
 * ghl_ADM_aux_quantities::g4UU in mind.
 *
 * @param[in] g4:       2D array containing the 4-metric or its inverse
 *
 * @param[in] vec:      1D array containing a 4-vector to be raised/lowered
 *
 * @param[out] vec_inv: 1D array containing the dual of `vec`
 * 
 * @returns void
 */
static inline void ghl_raise_lower_vector_4D(
      const double g4[][4],
      const double vec[4],
      double vec_inv[4]) {

  vec_inv[0] = g4[0][0] * vec[0] + g4[0][1] * vec[1]
             + g4[0][2] * vec[2] + g4[0][3] * vec[3];
                                                     
  vec_inv[1] = g4[0][1] * vec[0] + g4[1][1] * vec[1]
             + g4[1][2] * vec[2] + g4[1][3] * vec[3];
                                                     
  vec_inv[2] = g4[0][2] * vec[0] + g4[1][2] * vec[1]
             + g4[2][2] * vec[2] + g4[2][3] * vec[3];
                                                     
  vec_inv[3] = g4[0][3] * vec[0] + g4[1][3] * vec[1]
             + g4[2][3] * vec[2] + g4[3][3] * vec[3];
}

/**
 * @ingroup metric_contract
 * @brief Raise or lower a 3-vector
 *
 * @details
 * This function can compute both raising or lowering:
 *
 * - \f$ v^i = g^{i j} v_j \f$
 * - \f$ v_i = g_{i j} v^j \f$
 *
 * The exact operation being done is simply determined by the input
 * quantities. This function works for any input arrays, though it
 * is defined with ghl_metric_quantities::gammaDD and
 * ghl_metric_quantities::gammaUU in mind.
 *
 * @param[in] gamma:    2D array containing the 3-metric or its inverse
 *
 * @param[in] vec:      1D array containing a 3-vector to be raised/lowered
 *
 * @param[out] vec_inv: 1D array containing the dual of `vec`
 * 
 * @returns void
 */
static inline void ghl_raise_lower_vector_3D(
      const double gamma[][3],
      const double vec[3],
      double vec_inv[3]) {

  vec_inv[0] = gamma[0][0] * vec[0]
             + gamma[0][1] * vec[1]
             + gamma[0][2] * vec[2];
                                    
  vec_inv[1] = gamma[0][1] * vec[0]
             + gamma[1][1] * vec[1]
             + gamma[1][2] * vec[2];
                                    
  vec_inv[2] = gamma[0][2] * vec[0]
             + gamma[1][2] * vec[1]
             + gamma[2][2] * vec[2];
}

/**
 * @ingroup metric_contract
 * @brief Compute square of 4-vector
 *
 * @details
 * This function computes the square of the 4-vector
 * 
 * \f[
 * v^2 = v^\mu v_\mu
 * \f]
 *
 * using one of
 *
 * - \f$ v^2 = g^{\mu \nu} v_\mu v_\nu \f$
 * - \f$ v^2 = g_{\mu \nu} v^\mu v^\nu \f$
 *
 * The exact operation being done is simply determined by the input
 * quantities. This function works for any input arrays, though it
 * is defined with ghl_ADM_aux_quantities::g4DD and
 * ghl_ADM_aux_quantities::g4UU in mind.
 *
 * @param[in] gamma:    2D array containing the 4-metric or its inverse
 *
 * @param[in] vec:      1D array containing a 4-vector to be squared
 * 
 * @returns the value \f$  v^2 \f$
 */
static inline double ghl_compute_vec2_from_vec4D(
      const double g4[][4],
      const double vec[4]) {

  return (g4[0][0] * vec[0] * vec[0] +
          g4[1][1] * vec[1] * vec[1] +
          g4[2][2] * vec[2] * vec[2] +
          g4[3][3] * vec[3] * vec[3] +
          2.0 * (g4[0][1] * vec[0] * vec[1] +
                 g4[0][2] * vec[0] * vec[2] +
                 g4[0][3] * vec[0] * vec[3] +
                 g4[1][2] * vec[1] * vec[2] +
                 g4[1][3] * vec[1] * vec[3] +
                 g4[2][3] * vec[2] * vec[3]));
}

/**
 * @ingroup metric_contract
 * @brief Compute square of 3-vector
 *
 * @details
 * This function computes the square of the 3-vector
 * 
 * \f[
 * v^2 = v^i v_i
 * \f]
 *
 * using one of
 *
 * - \f$ v^2 = g^{i j} v_i v_j \f$
 * - \f$ v^2 = g_{i j} v^i v^j \f$
 *
 * The exact operation being done is simply determined by the input
 * quantities. This function works for any input arrays, though it
 * is defined with ghl_metric_quantities::gammaDD and
 * ghl_metric_quantities::gammaUU in mind.
 *
 * @param[in] gamma:    2D array containing the 3-metric or its inverse
 *
 * @param[in] vec:      1D array containing a 3-vector to be squared
 * 
 * @returns the value \f$  v^2 \f$
 */
static inline double ghl_compute_vec2_from_vec3D(
      const double gamma[][3],
      const double vec[3]) {

  return (gamma[0][0] * vec[0] * vec[0] +
          gamma[1][1] * vec[1] * vec[1] +
          gamma[2][2] * vec[2] * vec[2] +
          2.0 * (gamma[0][1] * vec[0] * vec[1] +
                 gamma[0][2] * vec[0] * vec[2] +
                 gamma[1][2] * vec[1] * vec[2]));
}

#endif // GHL_METRIC_HELPERS_H
