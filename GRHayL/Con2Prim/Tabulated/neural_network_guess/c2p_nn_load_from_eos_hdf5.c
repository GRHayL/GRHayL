#include "ghl_con2prim.h"

#ifdef GHL_DISABLE_HDF5
ghl_error_codes_t ghl_c2p_nn_load_hdf5(
      const char *nn_model_filepath,
      ghl_eos_parameters *restrict eos) {
  (void)nn_model_filepath;
  (void)eos;
  return ghl_error_used_disabled_hdf5;
}

ghl_error_codes_t ghl_c2p_nn_load_from_eos_hdf5(
      const char *tablepath,
      ghl_eos_parameters *restrict eos) {
  (void)tablepath;
  (void)eos;
  return ghl_error_used_disabled_hdf5;
}
#else

#include <hdf5.h>

#define GHL_NN_EOS_GROUP "grhayl_nn_c2p"
#define NN_C2P_OUT_X_BOUNDED 0
#define GHL_MAX_DIM_SIZE 8 // Leo says: out nets are tiny, so this is safe

static ghl_error_codes_t nn__build_dataset_path(
      const char *group_prefix,
      const char *dataset_name,
      char *buffer,
      const size_t buffer_size) {
  const int written = (group_prefix == NULL || group_prefix[0] == '\0')
                    ? snprintf(buffer, buffer_size, "%s", dataset_name)
                    : snprintf(buffer, buffer_size, "%s/%s", group_prefix, dataset_name);
  if(written < 0 || (size_t)written >= buffer_size) {
    return ghl_error_hdf5_dataset_could_not_open;
  }
  return ghl_success;
}

static ghl_error_codes_t nn__open_dataset(
      hid_t file_id,
      const char *dataset_name,
      hid_t *restrict dataset_id) {
  *dataset_id = H5Dopen(file_id, dataset_name, H5P_DEFAULT);
  if(*dataset_id < 0) {
    return ghl_error_hdf5_dataset_could_not_open;
  }
  return ghl_success;
}

static int nn__dataset_exists(hid_t file_id, const char *dataset_name) {
  return H5Lexists(file_id, dataset_name, H5P_DEFAULT) > 0;
}

static void nn__close_dataset(hid_t dataset_id, hid_t dataspace_id) {
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);
}

static ghl_error_codes_t nn__read_scalar(
      hid_t file_id,
      const char *dataset_name,
      hid_t mem_type,
      void *value) {
  hid_t dataset_id;
  ghl_error_codes_t err = nn__open_dataset(file_id, dataset_name, &dataset_id);
  if(err != ghl_success) {
    return err;
  }
  const hid_t dataspace_id = H5Dget_space(dataset_id);

  const int ndims = H5Sget_simple_extent_ndims(dataspace_id);
  if(ndims != 0) {
    nn__close_dataset(dataset_id, dataspace_id);
    return ghl_error_hdf5_dataset_invalid_ndims;
  }

  if(H5Dread(dataset_id, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value) < 0) {
    nn__close_dataset(dataset_id, dataspace_id);
    return ghl_error_hdf5_dataset_could_not_read;
  }

  nn__close_dataset(dataset_id, dataspace_id);
  return ghl_success;
}

static ghl_error_codes_t nn__read_array_exact(
      hid_t file_id,
      const char *dataset_name,
      hid_t mem_type,
      const int expected_rank,
      const hsize_t *restrict expected_dims,
      void *data) {
  hid_t dataset_id;
  ghl_error_codes_t err = nn__open_dataset(file_id, dataset_name, &dataset_id);
  if(err != ghl_success) {
    return err;
  }
  const hid_t dataspace_id = H5Dget_space(dataset_id);

  const int ndims = H5Sget_simple_extent_ndims(dataspace_id);
  if(ndims != expected_rank) {
    nn__close_dataset(dataset_id, dataspace_id);
    return ghl_error_hdf5_dataset_invalid_ndims;
  }

  hsize_t dims[expected_rank > 0 ? expected_rank : 1];
  H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
  for(int i = 0; i < expected_rank; ++i) {
    if(dims[i] != expected_dims[i]) {
      nn__close_dataset(dataset_id, dataspace_id);
      return ghl_error_hdf5_dataset_size_mismatch;
    }
  }

  if(H5Dread(dataset_id, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0) {
    nn__close_dataset(dataset_id, dataspace_id);
    return ghl_error_hdf5_dataset_could_not_read;
  }

  nn__close_dataset(dataset_id, dataspace_id);
  return ghl_success;
}

static int *nn__malloc_ints(size_t count) {
  return count == 0 ? NULL : malloc(count * sizeof(int));
}

static float *nn__malloc_floats(size_t count) {
  return count == 0 ? NULL : malloc(count * sizeof(float));
}

static ghl_error_codes_t nn__load_model_from_hdf5(
      const hid_t file_id,
      const char *group_prefix,
      const char *source_label,
      ghl_eos_parameters *restrict eos) {
  char dataset_path[256];
  ghl_error_codes_t err;

#define NN_BUILD_PATH(name) \
  do { \
    err = nn__build_dataset_path(group_prefix, (name), dataset_path, sizeof(dataset_path)); \
    if(err != ghl_success) goto cleanup_and_return; \
  } while(0)

  ghl_c2p_nn_model *model = malloc(sizeof(*model));
  if(model == NULL) {
    return ghl_error_out_of_memory;
  }
  memset(model, 0, sizeof(*model));

  NN_BUILD_PATH("dims/in_dim");
  err = nn__read_scalar(file_id, dataset_path, H5T_NATIVE_INT, &model->in_dim);
  if(err != ghl_success) goto cleanup_and_return;
  NN_BUILD_PATH("dims/hidden_dim");
  err = nn__read_scalar(file_id, dataset_path, H5T_NATIVE_INT, &model->hidden_dim);
  if(err != ghl_success) goto cleanup_and_return;
  NN_BUILD_PATH("dims/n_hidden");
  err = nn__read_scalar(file_id, dataset_path, H5T_NATIVE_INT, &model->n_hidden);
  if(err != ghl_success) goto cleanup_and_return;
  NN_BUILD_PATH("dims/out_dim");
  err = nn__read_scalar(file_id, dataset_path, H5T_NATIVE_INT, &model->out_dim);
  if(err != ghl_success) goto cleanup_and_return;
  NN_BUILD_PATH("meta/q_idx");
  err = nn__read_scalar(file_id, dataset_path, H5T_NATIVE_INT, &model->q_idx);
  if(err != ghl_success) goto cleanup_and_return;
  NN_BUILD_PATH("meta/s_idx");
  err = nn__read_scalar(file_id, dataset_path, H5T_NATIVE_INT, &model->s_idx);
  if(err != ghl_success) goto cleanup_and_return;
  NN_BUILD_PATH("scaling/x_eps");
  err = nn__read_scalar(file_id, dataset_path, H5T_NATIVE_FLOAT, &model->x_eps);
  if(err != ghl_success) goto cleanup_and_return;
  NN_BUILD_PATH("meta/y_eps");
  err = nn__read_scalar(file_id, dataset_path, H5T_NATIVE_FLOAT, &model->y_eps);
  if(err != ghl_success) goto cleanup_and_return;
  NN_BUILD_PATH("meta/dx_eps");
  err = nn__read_scalar(file_id, dataset_path, H5T_NATIVE_FLOAT, &model->dx_eps);
  if(err != ghl_success) goto cleanup_and_return;

  if(model->in_dim     <= 0 || model->in_dim     > GHL_MAX_DIM_SIZE ||
     model->hidden_dim <= 0 || model->hidden_dim > GHL_MAX_DIM_SIZE ||
     model->n_hidden   <= 0 || model->n_hidden   > GHL_MAX_DIM_SIZE ||
     model->out_dim    <= 0 || model->out_dim    > GHL_MAX_DIM_SIZE) {
    err = ghl_error_nn_c2p_invalid_dimensions;
    goto cleanup_and_return;
  }

  const hsize_t in_dim = (hsize_t)model->in_dim;
  const hsize_t hidden_dim = (hsize_t)model->hidden_dim;
  const hsize_t hidden_layers = (hsize_t)(model->n_hidden - 1);
  const hsize_t out_dim = (hsize_t)model->out_dim;

  model->x_kind = nn__malloc_ints((size_t)model->in_dim);
  model->x_lo = nn__malloc_floats((size_t)model->in_dim);
  model->x_hi = nn__malloc_floats((size_t)model->in_dim);
  model->x_invrng = nn__malloc_floats((size_t)model->in_dim);
  model->out_kind = nn__malloc_ints((size_t)model->out_dim);
  model->out_lo = nn__malloc_floats((size_t)model->out_dim);
  model->out_hi = nn__malloc_floats((size_t)model->out_dim);
  model->out_invrng = nn__malloc_floats((size_t)model->out_dim);
  model->W_in = nn__malloc_floats((size_t)(hidden_dim * in_dim));
  model->b_in = nn__malloc_floats((size_t)hidden_dim);
  model->W_hid = nn__malloc_floats((size_t)(hidden_layers * hidden_dim * hidden_dim));
  model->b_hid = nn__malloc_floats((size_t)(hidden_layers * hidden_dim));
  model->W_out = nn__malloc_floats((size_t)(out_dim * hidden_dim));
  model->b_out = nn__malloc_floats((size_t)out_dim);

  if(model->x_kind == NULL || model->x_lo == NULL || model->x_hi == NULL
        || model->x_invrng == NULL || model->out_kind == NULL
        || model->out_lo == NULL || model->out_hi == NULL
        || model->out_invrng == NULL || model->W_in == NULL
        || model->b_in == NULL || model->W_out == NULL || model->b_out == NULL
        || (hidden_layers > 0 && (model->W_hid == NULL || model->b_hid == NULL))) {
    err = ghl_error_out_of_memory;
    goto cleanup_and_return;
  }

  {
    const hsize_t dims1[1] = { in_dim };
    const hsize_t dims2_w_in[2] = { hidden_dim, in_dim };
    const hsize_t dims2_b_hid[2] = { hidden_layers, hidden_dim };
    const hsize_t dims3[3] = { hidden_layers, hidden_dim, hidden_dim };
    const hsize_t dims2_w_out[2] = { out_dim, hidden_dim };
    const hsize_t dims1_hidden[1] = { hidden_dim };
    const hsize_t dims1_out[1] = { out_dim };

    NN_BUILD_PATH("scaling/x_kind");
    err = nn__read_array_exact(file_id, dataset_path, H5T_NATIVE_INT, 1, dims1, model->x_kind);
    if(err != ghl_success) goto cleanup_and_return;
    NN_BUILD_PATH("scaling/x_lo");
    err = nn__read_array_exact(file_id, dataset_path, H5T_NATIVE_FLOAT, 1, dims1, model->x_lo);
    if(err != ghl_success) goto cleanup_and_return;
    NN_BUILD_PATH("scaling/x_hi");
    err = nn__read_array_exact(file_id, dataset_path, H5T_NATIVE_FLOAT, 1, dims1, model->x_hi);
    if(err != ghl_success) goto cleanup_and_return;
    NN_BUILD_PATH("scaling/x_invrng");
    err = nn__read_array_exact(file_id, dataset_path, H5T_NATIVE_FLOAT, 1, dims1, model->x_invrng);
    if(err != ghl_success) goto cleanup_and_return;
    NN_BUILD_PATH("scaling/out_kind");
    if(nn__dataset_exists(file_id, dataset_path)) {
      err = nn__read_array_exact(file_id, dataset_path, H5T_NATIVE_INT, 1, dims1_out, model->out_kind);
      if(err != ghl_success) goto cleanup_and_return;
      NN_BUILD_PATH("scaling/out_lo");
      err = nn__read_array_exact(file_id, dataset_path, H5T_NATIVE_FLOAT, 1, dims1_out, model->out_lo);
      if(err != ghl_success) goto cleanup_and_return;
      NN_BUILD_PATH("scaling/out_hi");
      err = nn__read_array_exact(file_id, dataset_path, H5T_NATIVE_FLOAT, 1, dims1_out, model->out_hi);
      if(err != ghl_success) goto cleanup_and_return;
      NN_BUILD_PATH("scaling/out_invrng");
      err = nn__read_array_exact(file_id, dataset_path, H5T_NATIVE_FLOAT, 1, dims1_out, model->out_invrng);
      if(err != ghl_success) goto cleanup_and_return;
    }
    else if(model->out_dim == 1) {
      model->out_kind[0] = NN_C2P_OUT_X_BOUNDED;
      model->out_lo[0] = 0.0f;
      model->out_hi[0] = 1.0f;
      model->out_invrng[0] = 1.0f;
    }
    else {
      (void)source_label;
      err = ghl_error_hdf5_dataset_could_not_open;
      goto cleanup_and_return;
    }
    NN_BUILD_PATH("layers/W_in");
    err = nn__read_array_exact(file_id, dataset_path, H5T_NATIVE_FLOAT, 2, dims2_w_in, model->W_in);
    if(err != ghl_success) goto cleanup_and_return;
    NN_BUILD_PATH("layers/b_in");
    err = nn__read_array_exact(file_id, dataset_path, H5T_NATIVE_FLOAT, 1, dims1_hidden, model->b_in);
    if(err != ghl_success) goto cleanup_and_return;
    if(hidden_layers > 0) {
      NN_BUILD_PATH("layers/W_hid");
      err = nn__read_array_exact(file_id, dataset_path, H5T_NATIVE_FLOAT, 3, dims3, model->W_hid);
      if(err != ghl_success) goto cleanup_and_return;
      NN_BUILD_PATH("layers/b_hid");
      err = nn__read_array_exact(file_id, dataset_path, H5T_NATIVE_FLOAT, 2, dims2_b_hid, model->b_hid);
      if(err != ghl_success) goto cleanup_and_return;
    }
    NN_BUILD_PATH("layers/W_out");
    err = nn__read_array_exact(file_id, dataset_path, H5T_NATIVE_FLOAT, 2, dims2_w_out, model->W_out);
    if(err != ghl_success) goto cleanup_and_return;
    NN_BUILD_PATH("layers/b_out");
    err = nn__read_array_exact(file_id, dataset_path, H5T_NATIVE_FLOAT, 1, dims1_out, model->b_out);
    if(err != ghl_success) goto cleanup_and_return;
  }

  err = ghl_c2p_nn_validate_model(model);
  if(err != ghl_success) {
    goto cleanup_and_return;
  }

  ghl_c2p_nn_free(eos->c2p_nn);
  eos->c2p_nn = model;
  return ghl_success;

cleanup_and_return:
  ghl_c2p_nn_free(model);
  return err;

#undef NN_BUILD_PATH
}

ghl_error_codes_t ghl_c2p_nn_load_hdf5(
      const char *nn_model_filepath,
      ghl_eos_parameters *restrict eos) {
  if(eos == NULL) {
    return ghl_error_eos_struct_is_null;
  }
  if(nn_model_filepath == NULL || nn_model_filepath[0] == '\0') {
    return ghl_error_could_not_open_file;
  }

  const hid_t file_id = H5Fopen(nn_model_filepath, H5F_ACC_RDONLY, H5P_DEFAULT);
  if(file_id < 0) {
    return ghl_error_could_not_open_file;
  }
  const ghl_error_codes_t err = nn__load_model_from_hdf5(file_id, NULL, nn_model_filepath, eos);
  H5Fclose(file_id);
  return err;
}

ghl_error_codes_t ghl_c2p_nn_load_from_eos_hdf5(
      const char *tablepath,
      ghl_eos_parameters *restrict eos) {
  if(eos == NULL) {
    return ghl_error_eos_struct_is_null;
  }
  if(tablepath == NULL || tablepath[0] == '\0') {
    return ghl_error_could_not_open_file;
  }

  const hid_t file_id = H5Fopen(tablepath, H5F_ACC_RDONLY, H5P_DEFAULT);
  if(file_id < 0) {
    return ghl_error_could_not_open_file;
  }
  if(H5Lexists(file_id, GHL_NN_EOS_GROUP, H5P_DEFAULT) <= 0) {
    H5Fclose(file_id);
    return ghl_error_hdf5_dataset_could_not_open;
  }
  const ghl_error_codes_t err = nn__load_model_from_hdf5(file_id, GHL_NN_EOS_GROUP, tablepath, eos);
  H5Fclose(file_id);
  return err;
}
#endif
