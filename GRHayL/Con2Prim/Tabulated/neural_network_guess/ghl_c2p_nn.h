#ifndef GHL_NN_C2P_H
#define GHL_NN_C2P_H

#include "ghl_con2prim.h"

#define GHL_NN_C2P_HANDLE_OTHER_KIND(kind, x_raw) (NAN)
#define GHL_NN_C2P_OUT_X_BOUNDED 0
#define GHL_NN_C2P_OUT_LINEAR 1
#define GHL_NN_C2P_OUT_LOG_LINEAR 2

static inline size_t nn__w_in_index(const ghl_c2p_nn_model *restrict model, size_t j, size_t i) {
  return j * (size_t)model->in_dim + i;
}

static inline size_t nn__w_hid_index(
      const ghl_c2p_nn_model *restrict model,
      size_t layer,
      size_t j,
      size_t k) {
  return ((layer * (size_t)model->hidden_dim) + j) * (size_t)model->hidden_dim + k;
}

static inline size_t nn__b_hid_index(const ghl_c2p_nn_model *restrict model, size_t layer, size_t j) {
  return layer * (size_t)model->hidden_dim + j;
}

static inline size_t nn__w_out_index(const ghl_c2p_nn_model *restrict model, size_t o, size_t k) {
  return o * (size_t)model->hidden_dim + k;
}

static inline float nn__hardtanhf(float x) {
  if(x < -1.0f) return -1.0f;
  if(x > 1.0f) return 1.0f;
  return x;
}

static inline float nn__sigmoidf(float z) {
  if(z >= 0.0f) {
    const float ez = expf(-z);
    return 1.0f / (1.0f + ez);
  }

  const float ez = expf(z);
  return ez / (1.0f + ez);
}

static inline float nn__clamp01_eps(float y, float eps) {
  const float hi = 1.0f - eps;
  if(y < eps) return eps;
  if(y > hi) return hi;
  return y;
}

static inline float nn__feature_transform(const ghl_c2p_nn_model *restrict model, int kind, float x_raw) {
  if(kind == 0) return x_raw;
  if(kind == 1) {
    float xc = x_raw;
    if(xc < model->x_eps) xc = model->x_eps;
    return log10f(xc);
  }
  return (float)GHL_NN_C2P_HANDLE_OTHER_KIND(kind, x_raw);
}

static inline float nn__clipf(float x, float lo, float hi) {
  if(x < lo) x = lo;
  if(x > hi) x = hi;
  return x;
}

static inline int nn__predict_y01_clamped_or_nan(
      const ghl_c2p_nn_model *restrict model,
      const float *restrict x_raw,
      float *restrict y01_out) {
  float x01[model->in_dim];
  float h0[model->hidden_dim];
  float h1[model->hidden_dim];
  float *h_prev = h0;
  float *h_curr = h1;

  for(size_t i = 0; i < (size_t)model->in_dim; ++i) {
    const float lo = model->x_lo[i];
    const float hi = model->x_hi[i];
    const float inv = model->x_invrng[i];

    if(!isfinite(lo) || !isfinite(hi) || !isfinite(inv) || !(inv > 0.0f)) {
      return 0;
    }

    const float xt = nn__feature_transform(model, model->x_kind[i], x_raw[i]);
    if(!isfinite(xt)) {
      return 0;
    }

    if(!(hi > lo)) {
      x01[i] = 0.0f;
      continue;
    }

    const float v = (nn__clipf(xt, lo, hi) - lo) * inv;
    if(!isfinite(v)) {
      return 0;
    }
    x01[i] = v;
  }

  for(size_t j = 0; j < (size_t)model->hidden_dim; ++j) {
    float acc = model->b_in[j];
    for(size_t i = 0; i < (size_t)model->in_dim; ++i) {
      acc += model->W_in[nn__w_in_index(model, j, i)] * x01[i];
    }
    h_prev[j] = nn__hardtanhf(acc);
  }

  for(size_t layer = 1; layer < (size_t)model->n_hidden; ++layer) {
    const size_t lh = layer - 1u;
    for(size_t j = 0; j < (size_t)model->hidden_dim; ++j) {
      float acc = model->b_hid[nn__b_hid_index(model, lh, j)];
      for(size_t k = 0; k < (size_t)model->hidden_dim; ++k) {
        acc += model->W_hid[nn__w_hid_index(model, lh, j, k)] * h_prev[k];
      }
      h_curr[j] = nn__hardtanhf(acc);
    }

    float *tmp = h_prev;
    h_prev = h_curr;
    h_curr = tmp;
  }

  for(size_t o = 0; o < (size_t)model->out_dim; ++o) {
    float z = model->b_out[o];
    for(size_t k = 0; k < (size_t)model->hidden_dim; ++k) {
      z += model->W_out[nn__w_out_index(model, o, k)] * h_prev[k];
    }

    const float y = nn__sigmoidf(z);
    if(!isfinite(y)) {
      return 0;
    }
    y01_out[o] = nn__clamp01_eps(y, model->y_eps);
  }
  return 1;
}

static inline float nn__decode_linear_output(
      const ghl_c2p_nn_model *restrict model,
      size_t output_idx,
      float y01) {
  const float lo = model->out_lo[output_idx];
  const float hi = model->out_hi[output_idx];
  const float inv = model->out_invrng[output_idx];
  if(!isfinite(lo) || !isfinite(hi) || !isfinite(inv) || !(inv > 0.0f) || !(hi > lo)) {
    return (float)NAN;
  }
  return lo + y01 / inv;
}

static inline float nn__decode_output(
      const ghl_c2p_nn_model *restrict model,
      size_t output_idx,
      float y01) {
  const int kind = model->out_kind[output_idx];
  if(kind == GHL_NN_C2P_OUT_LINEAR) {
    return nn__decode_linear_output(model, output_idx, y01);
  }
  if(kind == GHL_NN_C2P_OUT_LOG_LINEAR) {
    const float logv = nn__decode_linear_output(model, output_idx, y01);
    return isfinite(logv) ? expf(logv) : (float)NAN;
  }
  return (float)NAN;
}

#endif // GHL_NN_C2P_H
