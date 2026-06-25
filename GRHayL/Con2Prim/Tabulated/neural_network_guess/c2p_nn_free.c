#include "ghl_con2prim.h"

void ghl_c2p_nn_free(ghl_c2p_nn_model *restrict model) {
  if(model == NULL) return;

  free(model->x_kind);
  free(model->x_lo);
  free(model->x_hi);
  free(model->x_invrng);
  free(model->out_kind);
  free(model->out_lo);
  free(model->out_hi);
  free(model->out_invrng);
  free(model->W_in);
  free(model->b_in);
  free(model->W_hid);
  free(model->b_hid);
  free(model->W_out);
  free(model->b_out);
  free(model);
}
