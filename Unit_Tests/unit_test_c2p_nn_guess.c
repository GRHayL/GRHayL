#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ghl_con2prim.h"

#ifndef GHL_DISABLE_HDF5
#include <hdf5.h>
#endif

#define CHECK(cond, ...)                                                       \
  do {                                                                        \
    if(!(cond)) {                                                             \
      fprintf(stderr, __VA_ARGS__);                                           \
      fprintf(stderr, "\n");                                                  \
      exit(1);                                                                \
    }                                                                         \
  } while(0)

#define CHECK_ERROR(got, expected)                                             \
  CHECK((got) == (expected), "expected error %d but got %d",                  \
        (int)(expected), (int)(got))

static float invrng(const float lo, const float hi) {
  return 1.0f / (hi - lo);
}

static void set_valid_model(
      ghl_c2p_nn_model *restrict model,
      int *restrict x_kind,
      float *restrict x_lo,
      float *restrict x_hi,
      float *restrict x_invrng,
      int *restrict out_kind,
      float *restrict out_lo,
      float *restrict out_hi,
      float *restrict out_invrng,
      float *restrict W_in,
      float *restrict b_in,
      float *restrict W_hid,
      float *restrict b_hid,
      float *restrict W_out,
      float *restrict b_out) {
  memset(model, 0, sizeof(*model));
  model->in_dim = 4;
  model->hidden_dim = 2;
  model->n_hidden = 2;
  model->out_dim = 2;
  model->q_idx = 0;
  model->s_idx = 2;
  model->x_eps = 1e-30f;
  model->y_eps = 1e-6f;
  model->dx_eps = 1e-5f;
  model->x_kind = x_kind;
  model->x_lo = x_lo;
  model->x_hi = x_hi;
  model->x_invrng = x_invrng;
  model->out_kind = out_kind;
  model->out_lo = out_lo;
  model->out_hi = out_hi;
  model->out_invrng = out_invrng;
  model->W_in = W_in;
  model->b_in = b_in;
  model->W_hid = W_hid;
  model->b_hid = b_hid;
  model->W_out = W_out;
  model->b_out = b_out;

  x_kind[0] = 0;
  x_kind[1] = 1;
  x_kind[2] = 0;
  x_kind[3] = 1;
  x_lo[0] = 0.0f;
  x_lo[1] = -2.0f;
  x_lo[2] = 0.0f;
  x_lo[3] = -1.0f;
  x_hi[0] = 10.0f;
  x_hi[1] = 2.0f;
  x_hi[2] = 10.0f;
  x_hi[3] = 1.0f;
  for(int i = 0; i < 4; ++i) {
    x_invrng[i] = invrng(x_lo[i], x_hi[i]);
  }

  out_kind[0] = 0;
  out_kind[1] = 1;
  out_lo[0] = 0.0f;
  out_hi[0] = 1.0f;
  out_invrng[0] = 1.0f;
  out_lo[1] = 2.0f;
  out_hi[1] = 6.0f;
  out_invrng[1] = invrng(out_lo[1], out_hi[1]);

  for(int i = 0; i < 8; ++i) {
    W_in[i] = 0.0f;
    W_out[i] = 0.0f;
  }
  b_in[0] = 0.25f;
  b_in[1] = -0.5f;
  for(int i = 0; i < 4; ++i) {
    W_hid[i] = 0.0f;
  }
  b_hid[0] = 0.125f;
  b_hid[1] = -0.25f;
  b_out[0] = 0.0f;
  b_out[1] = -1.0986122886681098f;
}

static ghl_c2p_nn_model valid_stack_model(void) {
  static int x_kind[4];
  static float x_lo[4], x_hi[4], x_invrng[4];
  static int out_kind[2];
  static float out_lo[2], out_hi[2], out_invrng[2];
  static float W_in[8], b_in[2], W_hid[4], b_hid[2], W_out[8], b_out[2];
  ghl_c2p_nn_model model;
  set_valid_model(&model, x_kind, x_lo, x_hi, x_invrng,
                  out_kind, out_lo, out_hi, out_invrng,
                  W_in, b_in, W_hid, b_hid, W_out, b_out);
  return model;
}

static void test_validate_model(void) {
  ghl_c2p_nn_model model = valid_stack_model();
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model), ghl_success);
  CHECK_ERROR(ghl_c2p_nn_validate_model(NULL), ghl_error_nn_c2p_model_is_null);

  model = valid_stack_model();
  model.in_dim = 3;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_invalid_dimensions);

  model = valid_stack_model();
  model.q_idx = model.in_dim;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_invalid_input_index);

  model = valid_stack_model();
  model.y_eps = 0.5f;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_invalid_number);

  model = valid_stack_model();
  model.x_lo = NULL;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_missing_array);

  model = valid_stack_model();
  model.x_kind[1] = 7;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_invalid_kind);

  model = valid_stack_model();
  model.x_invrng[0] = 0.0f;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_invalid_number);

  model = valid_stack_model();
  model.out_kind[0] = 1;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_invalid_kind);

  model = valid_stack_model();
  model.b_hid[1] = NAN;
  CHECK_ERROR(ghl_c2p_nn_validate_model(&model),
              ghl_error_nn_c2p_invalid_number);
}

static void test_guess_model(void) {
  ghl_c2p_nn_model model = valid_stack_model();
  ghl_nn_c2p_input_t input = { 2.0f, 0.25f, 0.5f, 0.1f };
  ghl_nn_c2p_guess_t guess = ghl_c2p_nn_guess(&model, input);
  CHECK(fabsf(guess.x - 4.0f) < 1e-6f,
        "x-bounded guess mismatch: %.9g", guess.x);

  model = valid_stack_model();
  model.out_kind[1] = 2;
  model.out_lo[1] = logf(2.0f);
  model.out_hi[1] = logf(8.0f);
  model.out_invrng[1] = invrng(model.out_lo[1], model.out_hi[1]);
  model.b_out[1] = 0.0f;
  guess = ghl_c2p_nn_guess(&model, input);
  CHECK(fabsf(guess.x - 4.0f) < 1e-6f,
        "x guess with extra log-linear output mismatch: %.9g", guess.x);

  model = valid_stack_model();
  input.r = NAN;
  guess = ghl_c2p_nn_guess(&model, input);
  CHECK(fabsf(guess.x - 4.0f) < 1e-6f,
        "finite fallback x mismatch: %.9g", guess.x);

  model = valid_stack_model();
  model.x_kind[1] = 9;
  input.r = 0.25f;
  guess = ghl_c2p_nn_guess(&model, input);
  CHECK(fabsf(guess.x - 4.0f) < 1e-6f,
        "invalid transform fallback x mismatch: %.9g", guess.x);

  guess = ghl_c2p_nn_guess(NULL, input);
  CHECK(guess.x == 0.0f, "NULL model fallback failed");
}

#ifndef GHL_DISABLE_HDF5
static void create_group_checked(hid_t file_id, const char *name) {
  hid_t group_id = H5Gcreate2(file_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  CHECK(group_id >= 0, "failed to create HDF5 group %s", name);
  CHECK(H5Gclose(group_id) >= 0, "failed to close HDF5 group %s", name);
}

static void write_scalar_i32(hid_t file_id, const char *name, int value) {
  hid_t space_id = H5Screate(H5S_SCALAR);
  CHECK(space_id >= 0, "failed to create scalar dataspace");
  hid_t dataset_id = H5Dcreate2(file_id, name, H5T_NATIVE_INT, space_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  CHECK(dataset_id >= 0, "failed to create HDF5 dataset %s", name);
  CHECK(H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, &value) >= 0,
        "failed to write HDF5 dataset %s", name);
  H5Dclose(dataset_id);
  H5Sclose(space_id);
}

static void write_scalar_f32(hid_t file_id, const char *name, float value) {
  hid_t space_id = H5Screate(H5S_SCALAR);
  CHECK(space_id >= 0, "failed to create scalar dataspace");
  hid_t dataset_id = H5Dcreate2(file_id, name, H5T_NATIVE_FLOAT, space_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  CHECK(dataset_id >= 0, "failed to create HDF5 dataset %s", name);
  CHECK(H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, &value) >= 0,
        "failed to write HDF5 dataset %s", name);
  H5Dclose(dataset_id);
  H5Sclose(space_id);
}

static void write_array_i32(
      hid_t file_id,
      const char *name,
      int rank,
      const hsize_t *restrict dims,
      const int *restrict values) {
  hid_t space_id = H5Screate_simple(rank, dims, NULL);
  CHECK(space_id >= 0, "failed to create array dataspace");
  hid_t dataset_id = H5Dcreate2(file_id, name, H5T_NATIVE_INT, space_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  CHECK(dataset_id >= 0, "failed to create HDF5 dataset %s", name);
  CHECK(H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, values) >= 0,
        "failed to write HDF5 dataset %s", name);
  H5Dclose(dataset_id);
  H5Sclose(space_id);
}

static void write_array_f32(
      hid_t file_id,
      const char *name,
      int rank,
      const hsize_t *restrict dims,
      const float *restrict values) {
  hid_t space_id = H5Screate_simple(rank, dims, NULL);
  CHECK(space_id >= 0, "failed to create array dataspace");
  hid_t dataset_id = H5Dcreate2(file_id, name, H5T_NATIVE_FLOAT, space_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  CHECK(dataset_id >= 0, "failed to create HDF5 dataset %s", name);
  CHECK(H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, values) >= 0,
        "failed to write HDF5 dataset %s", name);
  H5Dclose(dataset_id);
  H5Sclose(space_id);
}

static void prefixed_name(
      const char *restrict prefix,
      const char *restrict name,
      char *restrict out,
      size_t out_size) {
  const int written = (prefix == NULL || prefix[0] == '\0')
                    ? snprintf(out, out_size, "%s", name)
                    : snprintf(out, out_size, "%s/%s", prefix, name);
  CHECK(written > 0 && (size_t)written < out_size,
        "HDF5 path buffer too small");
}

static void create_nn_groups(hid_t file_id, const char *prefix) {
  const char *groups[] = {
    "dims", "meta", "scaling", "layers"
  };
  if(prefix != NULL && prefix[0] != '\0') {
    create_group_checked(file_id, prefix);
  }
  for(size_t i = 0; i < sizeof(groups) / sizeof(groups[0]); ++i) {
    char name[256];
    prefixed_name(prefix, groups[i], name, sizeof(name));
    create_group_checked(file_id, name);
  }
}

static void write_valid_hdf5_model(
      hid_t file_id,
      const char *prefix,
      const bool include_out_scaling,
      const int out_dim) {
  create_nn_groups(file_id, prefix);

  char name[256];
  prefixed_name(prefix, "dims/in_dim", name, sizeof(name));
  write_scalar_i32(file_id, name, 4);
  prefixed_name(prefix, "dims/hidden_dim", name, sizeof(name));
  write_scalar_i32(file_id, name, 2);
  prefixed_name(prefix, "dims/n_hidden", name, sizeof(name));
  write_scalar_i32(file_id, name, 2);
  prefixed_name(prefix, "dims/out_dim", name, sizeof(name));
  write_scalar_i32(file_id, name, out_dim);
  prefixed_name(prefix, "meta/q_idx", name, sizeof(name));
  write_scalar_i32(file_id, name, 0);
  prefixed_name(prefix, "meta/s_idx", name, sizeof(name));
  write_scalar_i32(file_id, name, 2);
  prefixed_name(prefix, "scaling/x_eps", name, sizeof(name));
  write_scalar_f32(file_id, name, 1e-30f);
  prefixed_name(prefix, "meta/y_eps", name, sizeof(name));
  write_scalar_f32(file_id, name, 1e-6f);
  prefixed_name(prefix, "meta/dx_eps", name, sizeof(name));
  write_scalar_f32(file_id, name, 1e-5f);

  const hsize_t dims_in[1] = { 4 };
  const hsize_t dims_hidden[1] = { 2 };
  const hsize_t dims_out[1] = { (hsize_t)out_dim };
  const hsize_t dims_w_in[2] = { 2, 4 };
  const hsize_t dims_b_hid[2] = { 1, 2 };
  const hsize_t dims_w_hid[3] = { 1, 2, 2 };
  const hsize_t dims_w_out[2] = { (hsize_t)out_dim, 2 };
  const int x_kind[4] = { 0, 1, 0, 1 };
  const float x_lo[4] = { 0.0f, -2.0f, 0.0f, -1.0f };
  const float x_hi[4] = { 10.0f, 2.0f, 10.0f, 1.0f };
  const float x_invrng[4] = { 0.1f, 0.25f, 0.1f, 0.5f };
  const int out_kind[2] = { 0, 1 };
  const float out_lo[2] = { 0.0f, 2.0f };
  const float out_hi[2] = { 1.0f, 6.0f };
  const float out_invrng[2] = { 1.0f, 0.25f };
  const float W_in[8] = { 0 };
  const float b_in[2] = { 0.25f, -0.5f };
  const float W_hid[4] = { 0 };
  const float b_hid[2] = { 0.125f, -0.25f };
  const float W_out[4] = { 0 };
  const float b_out[2] = { 0.0f, -1.0986122886681098f };

  prefixed_name(prefix, "scaling/x_kind", name, sizeof(name));
  write_array_i32(file_id, name, 1, dims_in, x_kind);
  prefixed_name(prefix, "scaling/x_lo", name, sizeof(name));
  write_array_f32(file_id, name, 1, dims_in, x_lo);
  prefixed_name(prefix, "scaling/x_hi", name, sizeof(name));
  write_array_f32(file_id, name, 1, dims_in, x_hi);
  prefixed_name(prefix, "scaling/x_invrng", name, sizeof(name));
  write_array_f32(file_id, name, 1, dims_in, x_invrng);
  if(include_out_scaling) {
    prefixed_name(prefix, "scaling/out_kind", name, sizeof(name));
    write_array_i32(file_id, name, 1, dims_out, out_kind);
    prefixed_name(prefix, "scaling/out_lo", name, sizeof(name));
    write_array_f32(file_id, name, 1, dims_out, out_lo);
    prefixed_name(prefix, "scaling/out_hi", name, sizeof(name));
    write_array_f32(file_id, name, 1, dims_out, out_hi);
    prefixed_name(prefix, "scaling/out_invrng", name, sizeof(name));
    write_array_f32(file_id, name, 1, dims_out, out_invrng);
  }
  prefixed_name(prefix, "layers/W_in", name, sizeof(name));
  write_array_f32(file_id, name, 2, dims_w_in, W_in);
  prefixed_name(prefix, "layers/b_in", name, sizeof(name));
  write_array_f32(file_id, name, 1, dims_hidden, b_in);
  prefixed_name(prefix, "layers/W_hid", name, sizeof(name));
  write_array_f32(file_id, name, 3, dims_w_hid, W_hid);
  prefixed_name(prefix, "layers/b_hid", name, sizeof(name));
  write_array_f32(file_id, name, 2, dims_b_hid, b_hid);
  prefixed_name(prefix, "layers/W_out", name, sizeof(name));
  write_array_f32(file_id, name, 2, dims_w_out, W_out);
  prefixed_name(prefix, "layers/b_out", name, sizeof(name));
  write_array_f32(file_id, name, 1, dims_out, b_out);
}

static void create_hdf5_file(
      const char *restrict path,
      const char *restrict prefix,
      const bool include_out_scaling,
      const int out_dim) {
  hid_t file_id = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  CHECK(file_id >= 0, "failed to create HDF5 file %s", path);
  write_valid_hdf5_model(file_id, prefix, include_out_scaling, out_dim);
  CHECK(H5Fclose(file_id) >= 0, "failed to close HDF5 file %s", path);
}

static void test_hdf5_loaders(void) {
  const ghl_nn_c2p_input_t input = { 2.0f, 0.25f, 0.5f, 0.1f };
  const int out_kind_linear = 1;

  const char root_path[] = "/tmp/unit_test_c2p_nn_root.h5";
  const char embedded_path[] = "/tmp/unit_test_c2p_nn_embedded.h5";
  const char legacy_path[] = "/tmp/unit_test_c2p_nn_legacy.h5";

  create_hdf5_file(root_path, "", true, 2);
  ghl_eos_parameters eos = { 0 };
  CHECK_ERROR(ghl_c2p_nn_load_hdf5(root_path, &eos), ghl_success);
  CHECK(eos.c2p_nn != NULL, "direct NN HDF5 load returned NULL model");
  CHECK(eos.c2p_nn->out_dim == 2, "direct HDF5 out_dim mismatch");
  CHECK(eos.c2p_nn->out_kind[1] == out_kind_linear,
        "direct HDF5 output-kind metadata mismatch");
  ghl_nn_c2p_guess_t guess = ghl_c2p_nn_guess(eos.c2p_nn, input);
  CHECK(fabsf(guess.x - 4.0f) < 1e-6f, "direct HDF5 x mismatch");
  ghl_c2p_nn_free(eos.c2p_nn);
  eos.c2p_nn = NULL;

  create_hdf5_file(embedded_path, "grhayl_nn_c2p", true, 2);
  CHECK_ERROR(ghl_c2p_nn_load_from_eos_hdf5(embedded_path, &eos),
              ghl_success);
  CHECK(eos.c2p_nn != NULL, "embedded NN HDF5 load returned NULL model");
  CHECK(eos.c2p_nn->out_dim == 2, "embedded HDF5 out_dim mismatch");
  CHECK(eos.c2p_nn->out_kind[1] == out_kind_linear,
        "embedded HDF5 output-kind metadata mismatch");
  guess = ghl_c2p_nn_guess(eos.c2p_nn, input);
  CHECK(fabsf(guess.x - 4.0f) < 1e-6f, "embedded HDF5 x mismatch");
  ghl_c2p_nn_free(eos.c2p_nn);
  eos.c2p_nn = NULL;

  create_hdf5_file(legacy_path, "", false, 1);
  CHECK_ERROR(ghl_c2p_nn_load_hdf5(legacy_path, &eos), ghl_success);
  CHECK(eos.c2p_nn != NULL, "legacy NN HDF5 load returned NULL model");
  CHECK(eos.c2p_nn->out_dim == 1, "legacy HDF5 out_dim mismatch");
  guess = ghl_c2p_nn_guess(eos.c2p_nn, input);
  CHECK(fabsf(guess.x - 4.0f) < 1e-6f, "legacy HDF5 x mismatch");
  ghl_c2p_nn_free(eos.c2p_nn);
  eos.c2p_nn = NULL;
}
#endif

int main(void) {
  test_validate_model();
  test_guess_model();
#ifndef GHL_DISABLE_HDF5
  test_hdf5_loaders();
#endif
  printf("All c2p neural-network guess tests succeeded\n");
  return 0;
}
