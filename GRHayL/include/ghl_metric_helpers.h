#ifndef GHL_METRIC_HELPERS_H_
#define GHL_METRIC_HELPERS_H_

inline void ghl_raise_lower_vector_4D(
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

inline void ghl_raise_lower_vector_3D(
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

inline double ghl_compute_vec2_from_vec3D(
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
