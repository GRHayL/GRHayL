#include "ghl_c2p_nn.h"

ghl_nn_c2p_guess_t ghl_c2p_nn_guess(
      const ghl_c2p_nn_model *restrict model,
      ghl_nn_c2p_input_t input) {
  ghl_nn_c2p_guess_t guess = { 0.0f };
  const float q = input.q;
  const float r = input.r;
  const float s = input.s;
  const float t = input.t;

  if(model == NULL || !isfinite(q) || !isfinite(s)) {
    return guess;
  }

  const float x_lo = 1.0f + q - s;
  float width = 1.0f + q;

  if(!isfinite(x_lo) || !isfinite(width)) {
    return guess;
  }
  if(width < model->dx_eps) {
    width = model->dx_eps;
  }

  guess.x = x_lo + 0.5f * width;
  if(!isfinite(r) || !isfinite(t)) {
    return guess;
  }

  const float x_raw[4] = { q, r, s, t };
  float y01[model->out_dim];
  for(size_t o = 0; o < (size_t)model->out_dim; ++o) {
    y01[o] = (float)NAN;
  }
  if(!nn__predict_y01_clamped_or_nan(model, x_raw, y01)) {
    return guess;
  }

  const float x_pred = x_lo + y01[0] * width;
  if(isfinite(x_pred)) {
    guess.x = x_pred;
  }

  return guess;
}
