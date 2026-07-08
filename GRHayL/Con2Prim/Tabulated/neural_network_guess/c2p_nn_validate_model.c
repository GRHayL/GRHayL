#include "ghl_c2p_nn.h"

ghl_error_codes_t ghl_c2p_nn_validate_model(const ghl_c2p_nn_model *restrict model) {
  if(model == NULL) {
    return ghl_error_nn_c2p_model_is_null;
  }
  if(model->in_dim != 4) {
    return ghl_error_nn_c2p_invalid_dimensions;
  }
  if(model->hidden_dim <= 0 || model->n_hidden <= 0 || model->out_dim <= 0) {
    return ghl_error_nn_c2p_invalid_dimensions;
  }
  if(model->q_idx < 0 || model->q_idx >= model->in_dim) {
    return ghl_error_nn_c2p_invalid_input_index;
  }
  if(model->s_idx < 0 || model->s_idx >= model->in_dim) {
    return ghl_error_nn_c2p_invalid_input_index;
  }
  if(!isfinite(model->x_eps) || !(model->x_eps > 0.0f)) {
    return ghl_error_nn_c2p_invalid_number;
  }
  if(!isfinite(model->y_eps) || !(model->y_eps > 0.0f) || !(model->y_eps < 0.5f)) {
    return ghl_error_nn_c2p_invalid_number;
  }
  if(!isfinite(model->dx_eps) || !(model->dx_eps > 0.0f)) {
    return ghl_error_nn_c2p_invalid_number;
  }
  if(model->x_kind == NULL || model->x_lo == NULL || model->x_hi == NULL || model->x_invrng == NULL) {
    return ghl_error_nn_c2p_missing_array;
  }
  if(model->out_kind == NULL || model->out_lo == NULL || model->out_hi == NULL || model->out_invrng == NULL) {
    return ghl_error_nn_c2p_missing_array;
  }
  if(model->W_in == NULL || model->b_in == NULL || model->W_out == NULL || model->b_out == NULL) {
    return ghl_error_nn_c2p_missing_array;
  }
  if(model->n_hidden > 1 && (model->W_hid == NULL || model->b_hid == NULL)) {
    return ghl_error_nn_c2p_missing_array;
  }

  for(size_t i = 0; i < (size_t)model->in_dim; ++i) {
    const int kind = model->x_kind[i];
    if(kind != 0 && kind != 1) {
      return ghl_error_nn_c2p_invalid_kind;
    }
    if(!isfinite(model->x_lo[i]) || !isfinite(model->x_hi[i]) || !isfinite(model->x_invrng[i])) {
      return ghl_error_nn_c2p_invalid_number;
    }
    if(!(model->x_invrng[i] > 0.0f) || !(model->x_hi[i] > model->x_lo[i])) {
      return ghl_error_nn_c2p_invalid_number;
    }
  }

  for(size_t i = 0; i < (size_t)model->out_dim; ++i) {
    const int kind = model->out_kind[i];
    if(i == 0u) {
      if(kind != GHL_NN_C2P_OUT_X_BOUNDED) {
        return ghl_error_nn_c2p_invalid_kind;
      }
    }
    else {
      if(kind != GHL_NN_C2P_OUT_LINEAR && kind != GHL_NN_C2P_OUT_LOG_LINEAR) {
        return ghl_error_nn_c2p_invalid_kind;
      }
      if(!isfinite(model->out_lo[i]) || !isfinite(model->out_hi[i]) || !isfinite(model->out_invrng[i])) {
        return ghl_error_nn_c2p_invalid_number;
      }
      if(!(model->out_invrng[i] > 0.0f) || !(model->out_hi[i] > model->out_lo[i])) {
        return ghl_error_nn_c2p_invalid_number;
      }
    }
  }

  for(size_t i = 0; i < (size_t)model->hidden_dim; ++i) {
    if(!isfinite(model->b_in[i])) {
      return ghl_error_nn_c2p_invalid_number;
    }
    for(size_t j = 0; j < (size_t)model->in_dim; ++j) {
      if(!isfinite(model->W_in[nn__w_in_index(model, i, j)])) {
        return ghl_error_nn_c2p_invalid_number;
      }
    }
  }

  for(size_t i = 0; i < (size_t)(model->n_hidden - 1); ++i) {
    for(size_t j = 0; j < (size_t)model->hidden_dim; ++j) {
      if(!isfinite(model->b_hid[nn__b_hid_index(model, i, j)])) {
        return ghl_error_nn_c2p_invalid_number;
      }
      for(size_t k = 0; k < (size_t)model->hidden_dim; ++k) {
        if(!isfinite(model->W_hid[nn__w_hid_index(model, i, j, k)])) {
          return ghl_error_nn_c2p_invalid_number;
        }
      }
    }
  }

  for(size_t i = 0; i < (size_t)model->out_dim; ++i) {
    if(!isfinite(model->b_out[i])) {
      return ghl_error_nn_c2p_invalid_number;
    }
    for(size_t j = 0; j < (size_t)model->hidden_dim; ++j) {
      if(!isfinite(model->W_out[nn__w_out_index(model, i, j)])) {
        return ghl_error_nn_c2p_invalid_number;
      }
    }
  }

  return ghl_success;
}
