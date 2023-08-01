#ifndef GHL_METRIC_HELPERS_H_
#define GHL_METRIC_HELPERS_H_

#define ghl_raise_lower_vector_4D(g4, vec, vec_inv)   \
  vec_inv[0] = g4[0][0] * vec[0] + g4[0][1] * vec[1]  \
             + g4[0][2] * vec[2] + g4[0][3] * vec[3]; \
                                                      \
  vec_inv[1] = g4[0][1] * vec[0] + g4[1][1] * vec[1]  \
             + g4[1][2] * vec[2] + g4[1][3] * vec[3]; \
                                                      \
  vec_inv[2] = g4[0][2] * vec[0] + g4[1][2] * vec[1]  \
             + g4[2][2] * vec[2] + g4[2][3] * vec[3]; \
                                                      \
  vec_inv[3] = g4[0][3] * vec[0] + g4[1][3] * vec[1]  \
             + g4[2][3] * vec[2] + g4[3][3] * vec[3];

#define ghl_raise_lower_vector_3D(gamma, vec, vec_inv) \
  vec_inv[0] = gamma[0][0] * vec[0]  \
             + gamma[0][1] * vec[1]  \
             + gamma[0][2] * vec[2]; \
                                     \
  vec_inv[1] = gamma[0][1] * vec[0]  \
             + gamma[1][1] * vec[1]  \
             + gamma[1][2] * vec[2]; \
                                     \
  vec_inv[2] = gamma[0][2] * vec[0]  \
             + gamma[1][2] * vec[1]  \
             + gamma[2][2] * vec[2];

#define ghl_compute_vec2_from_vec4(g4, vec) \
  (g4[0][0] * vec[0] * vec[0] +        \
   g4[1][1] * vec[1] * vec[1] +        \
   g4[2][2] * vec[2] * vec[2] +        \
   g4[3][3] * vec[3] * vec[3] +        \
   2.0 * (g4[0][1] * vec[0] * vec[1] + \
          g4[0][2] * vec[0] * vec[2] + \
          g4[0][3] * vec[0] * vec[3] + \
          g4[1][2] * vec[1] * vec[2] + \
          g4[1][3] * vec[1] * vec[3] + \
          g4[2][3] * vec[2] * vec[3]))

#define ghl_compute_vec2_from_vec(gamma, vec) \
  (gamma[0][0] * vec[0] * vec[0] +        \
   gamma[1][1] * vec[1] * vec[1] +        \
   gamma[2][2] * vec[2] * vec[2] +        \
   2.0 * (gamma[0][1] * vec[0] * vec[1] + \
          gamma[0][2] * vec[0] * vec[2] + \
          gamma[1][2] * vec[1] * vec[2]))

#endif // GHL_METRIC_HELPERS_H
