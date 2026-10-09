#!/bin/bash

set -Eeuxo pipefail

# Usage: .github/run_tests.sh [all|m1] [configure arguments...]
#   all (default)  build and run every unit test
#   m1             build and run only the Radiation M1 tests
# Any remaining arguments are passed to ./configure -r, for example
# --noomp, --disable-hdf5, or --cflags='-ftest-coverage -fprofile-arcs'.
# Runtime library paths follow the configured build directory.
suite=all
case "${1:-}" in
  all|m1) suite="$1"; shift ;;
esac

repo_root=$(pwd)
./configure -r "$@"
configured_build_dir=$(awk '
  /^BUILDDIR[[:space:]]*=/ {
    sub(/^BUILDDIR[[:space:]]*=[[:space:]]*/, "")
    print
    exit
  }
' Makefile)
if [[ -z "$configured_build_dir" ]]; then
  echo "Configured build directory is missing or empty in Makefile" >&2
  exit 1
fi
if [[ "$suite" == m1 ]]; then
  make tests
else
  make tests datagen
fi

if ! configured_lib_dir=$(cd "$repo_root" && cd -- "$configured_build_dir/lib" && pwd -P); then
  echo "Configured library directory is unavailable: $configured_build_dir/lib" >&2
  exit 1
fi
LD_LIBRARY_PATH="$configured_lib_dir${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH
DYLD_LIBRARY_PATH="$configured_lib_dir${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
export DYLD_LIBRARY_PATH

echo $LD_LIBRARY_PATH

created_paths=()
created_directories=()
cleanup_created_paths() {
  if [[ ${created_paths[*]-} ]]; then
    rm -f -- "${created_paths[@]}"
  fi
  if [[ ${created_directories[*]-} ]]; then
    rm -rf -- "${created_directories[@]}"
  fi
}
trap cleanup_created_paths EXIT

download_file() {
  url="$1"
  filename="${url##*/}"

  if [[ "$filename" == *.bz2 ]]; then
    uncompressed="${filename%.bz2}"

    if [ -f "$filename" ] || [ -f "$uncompressed" ]; then
      echo "File \"$filename\" or \"$uncompressed\" already exists; skipping it"
      return
    fi
  else
    if [ -f "$filename" ]; then
      echo "File \"$filename\" already exists; skipping it"
      return
    fi
  fi

  created_paths+=("$filename")
  curl -fL --retry 5 -O "$url"
}

download_test_data() {
  filepath="$1"
  url="${test_data_base_url}/${filepath}"
  download_file "$url"
}

run_m1_tests() {
  # Retained fixtures come from pinned TestData downloads unless the caller
  # supplies a raw or gzip directory. Both downloaded and decompressed data
  # are runner-owned temporary directories; caller inputs remain untouched.
  source .github/prepare_m1_fixtures.sh

  m1_fixture_source="${M1_FIXTURE_DIR:-}"
  m1_fixture_dir=""
  if [[ -z "$m1_fixture_source" ]]; then
    local radiation_testdata_ref
    radiation_testdata_ref=$(cat .github/radiation-testdata-ref)
    if [[ ! "$radiation_testdata_ref" =~ ^[0-9a-f]{40}$ ]]; then
      echo "Radiation TestData revision must be a published full commit SHA" >&2
      exit 1
    fi
    local test_data_base_url="https://raw.githubusercontent.com/GRHayL/TestData/${radiation_testdata_ref}"
    m1_fixture_source=$(mktemp -d "${TMPDIR:-/tmp}/grhayl_m1_downloads.XXXXXXXX")
    created_directories+=("$m1_fixture_source")
    for m1_fixture_member in "${m1_fixture_members[@]}"; do
      # The parent EXIT trap owns the whole download directory, including any
      # partial curl output; no subshell-created path needs separate cleanup.
      (cd "$m1_fixture_source" && download_test_data "radiation/$m1_fixture_member.gz")
    done
  fi
  if [[ -n "$m1_fixture_source" ]]; then
    if [[ ! -d "$m1_fixture_source" ]]; then
      echo "M1_FIXTURE_DIR is not a directory: $m1_fixture_source" >&2
      exit 1
    fi

    m1_fixture_is_raw=1
    for m1_fixture_member in "${m1_fixture_members[@]}"; do
      if [[ ! -f "$m1_fixture_source/$m1_fixture_member" ]]; then
        m1_fixture_is_raw=0
        break
      fi
    done

    if (( m1_fixture_is_raw )); then
      # Complete raw input is replayed directly from the supplied directory,
      # which stays caller-owned and untouched; any .gz counterparts are
      # redundant because raw members take precedence.
      m1_fixture_dir="$m1_fixture_source"
    else
      # Register the temporary destination with the existing EXIT cleanup
      # before preparation, and assign m1_fixture_dir only after complete
      # preparation so no partial directory can be replayed.
      m1_prepared_fixture_dir=$(mktemp -d "${TMPDIR:-/tmp}/grhayl_m1_fixtures.XXXXXXXX")
      created_directories+=("$m1_prepared_fixture_dir")
      prepare_m1_fixtures "$m1_fixture_source" "$m1_prepared_fixture_dir"
      m1_fixture_dir="$m1_prepared_fixture_dir"
    fi
  fi

  run_m1_fixture_test() {
    local test_executable="$1"
    local fixture_option="$2"
    local fixture_path="$3"
    "$test_executable" "$fixture_option" "$fixture_path"
  }

  ./test/unit_test_m1_closure_fallback
  run_m1_fixture_test ./test/unit_test_m1_diffusion_flux --fixture \
    "$m1_fixture_dir/jthick_thcm1.bin"
  ./test/unit_test_m1_error_handling
  ./test/unit_test_m1_fd_jacobian
  run_m1_fixture_test ./test/unit_test_m1_neutrino_rusanov_flux --fixture-dir \
    "$m1_fixture_dir"
  run_m1_fixture_test ./test/unit_test_m1_neutrino_seeded_invariants --fixture-dir \
    "$m1_fixture_dir"
  ./test/unit_test_m1_neutrino_source_update
  ./test/unit_test_m1_rate_provider --generated-fixture
  run_m1_fixture_test ./test/unit_test_m1_thcm1_blended_rusanov --fixture-dir \
    "$m1_fixture_dir"
  run_m1_fixture_test ./test/unit_test_rusanov_flux --fixture-dir \
    "$m1_fixture_dir"
}

if [[ "$suite" == m1 ]]; then
  run_m1_tests
  exit 0
fi

decompress_bz2() {
  archive="$1"
  uncompressed="${archive%.bz2}"
  if [ ! -f "$uncompressed" ]; then
    created_paths+=("$uncompressed")
    bunzip2 -k "$archive"
  fi
}

et_legacy_testdata_ref=$(cat .github/et-legacy-testdata-ref)
test_data_base_url="https://raw.githubusercontent.com/GRHayL/TestData/${et_legacy_testdata_ref}"
download_test_data ET_Legacy/ET_Legacy_conservs_input.bin
download_test_data ET_Legacy/ET_Legacy_conservs_output.bin
download_test_data ET_Legacy/ET_Legacy_conservs_output_pert.bin

download_test_data ET_Legacy/ET_Legacy_primitives_input.bin
download_test_data ET_Legacy/ET_Legacy_primitives_output.bin
download_test_data ET_Legacy/ET_Legacy_primitives_output_pert.bin

download_test_data ET_Legacy/ET_Legacy_induction_gauge_rhs_input.bin
download_test_data ET_Legacy/ET_Legacy_induction_gauge_rhs_output.bin
download_test_data ET_Legacy/ET_Legacy_induction_gauge_rhs_output_pert.bin

download_test_data ET_Legacy/ET_Legacy_HLL_flux_input.bin
download_test_data ET_Legacy/ET_Legacy_HLL_flux_output.bin
download_test_data ET_Legacy/ET_Legacy_HLL_flux_output_pert.bin

download_test_data ET_Legacy/ET_Legacy_reconstruction_input.bin
download_test_data ET_Legacy/ET_Legacy_reconstruction_output.bin
download_test_data ET_Legacy/ET_Legacy_reconstruction_output_pert.bin

download_test_data ET_Legacy/ET_Legacy_flux_source_input.bin
download_test_data ET_Legacy/ET_Legacy_flux_source_output.bin
download_test_data ET_Legacy/ET_Legacy_flux_source_output_pert.bin

./test/unit_test_ET_Legacy_conservs
./test/unit_test_ET_Legacy_primitives
./test/unit_test_ET_Legacy_induction_gauge_rhs
./test/unit_test_ET_Legacy_HLL_flux
./test/unit_test_ET_Legacy_reconstruction
./test/unit_test_ET_Legacy_flux_source

# These coupled fixtures must match the generator replay and CI consumers.
con2prim_testdata_ref=$(cat .github/con2prim-testdata-ref)
test_data_base_url="https://raw.githubusercontent.com/GRHayL/TestData/${con2prim_testdata_ref}"
download_test_data con2prim/metric_Bfield_initial_data.bin

download_test_data con2prim/apply_conservative_limits_input.bin
download_test_data con2prim/apply_conservative_limits_output.bin
download_test_data con2prim/apply_conservative_limits_output_pert.bin

download_test_data con2prim/con2prim_multi_method_hybrid_input.bin
download_test_data con2prim/con2prim_multi_method_hybrid_output.bin
download_test_data con2prim/con2prim_multi_method_hybrid_output_pert.bin

download_test_data con2prim/enforce_primitive_limits_and_compute_u0_input.bin
download_test_data con2prim/enforce_primitive_limits_and_compute_u0_output.bin
download_test_data con2prim/enforce_primitive_limits_and_compute_u0_output_pert.bin

download_test_data con2prim/compute_conservs_and_Tmunu_input.bin
download_test_data con2prim/compute_conservs_and_Tmunu_output.bin
download_test_data con2prim/compute_conservs_and_Tmunu_output_pert.bin
test_data_base_url="https://raw.githubusercontent.com/GRHayL/TestData/${et_legacy_testdata_ref}"
./test/unit_test_apply_conservative_limits
./test/unit_test_con2prim_multi_method_hybrid
./test/unit_test_enforce_primitive_limits_and_compute_u0
./test/unit_test_compute_conservs_and_Tmunu

./test/unit_test_hybrid_failure
./test/unit_test_c2p_nn_guess

download_test_data EOS/simple_table.h5
./test/unit_test_tabulated_eos simple_table.h5

./test/unit_test_piecewise_polytrope

download_test_data grhayl_core/grhayl_core_test_suite_input.bin
./test/unit_test_grhayl_core_test_suite

download_file https://stellarcollapse.org/EOS/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2
decompress_bz2 LS220_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2

download_test_data flux_source/hybrid_flux_input.bin
download_test_data flux_source/hybrid_flux_output.bin
download_test_data flux_source/hybrid_flux_output_pert.bin

download_test_data flux_source/tabulated_flux_input.bin
download_test_data flux_source/tabulated_flux_output.bin
download_test_data flux_source/tabulated_flux_output_pert.bin

./test/unit_test_hybrid_flux
./test/unit_test_tabulated_flux

download_test_data reconstruction/PLM_reconstruction_input.bin
download_test_data reconstruction/PLM_reconstruction_output.bin
download_test_data reconstruction/PLM_reconstruction_output_pert.bin

./test/unit_test_PLM_reconstruction

download_file https://stellarcollapse.org/EOS/SLy4_3335_rho391_temp163_ye66.h5.bz2
decompress_bz2 SLy4_3335_rho391_temp163_ye66.h5.bz2

download_test_data Neutrinos/nrpyleakage_optically_thin_gas_unperturbed.bin
download_test_data Neutrinos/nrpyleakage_optically_thin_gas_perturbed.bin

download_test_data Neutrinos/nrpyleakage_constant_density_sphere_unperturbed.bin
download_test_data Neutrinos/nrpyleakage_constant_density_sphere_perturbed.bin

download_test_data Neutrinos/nrpyleakage_luminosities_unperturbed.bin
download_test_data Neutrinos/nrpyleakage_luminosities_perturbed.bin

./test/unit_test_nrpyleakage_physics
./test/unit_test_nrpyleakage_classifier_fallback
./test/unit_test_nrpyleakage_optically_thin_gas SLy4_3335_rho391_temp163_ye66.h5 1
./test/unit_test_nrpyleakage_constant_density_sphere SLy4_3335_rho391_temp163_ye66.h5 1
./test/unit_test_nrpyleakage_luminosities SLy4_3335_rho391_temp163_ye66.h5 1

test_data_base_url="https://raw.githubusercontent.com/GRHayL/TestData/${con2prim_testdata_ref}"
download_test_data con2prim/con2prim_tabulated_Palenzuela1D_rho_vs_T_unperturbed.bin
download_test_data con2prim/con2prim_tabulated_Palenzuela1D_Pmag_vs_Wm1_unperturbed.bin
download_test_data con2prim/con2prim_tabulated_Palenzuela1D_rho_vs_T_perturbed.bin
download_test_data con2prim/con2prim_tabulated_Palenzuela1D_Pmag_vs_Wm1_perturbed.bin

download_test_data con2prim/con2prim_tabulated_Palenzuela1D_entropy_rho_vs_T_unperturbed.bin
download_test_data con2prim/con2prim_tabulated_Palenzuela1D_entropy_Pmag_vs_Wm1_unperturbed.bin
download_test_data con2prim/con2prim_tabulated_Palenzuela1D_entropy_rho_vs_T_perturbed.bin
download_test_data con2prim/con2prim_tabulated_Palenzuela1D_entropy_Pmag_vs_Wm1_perturbed.bin

download_test_data con2prim/con2prim_tabulated_Newman1D_rho_vs_T_unperturbed.bin
download_test_data con2prim/con2prim_tabulated_Newman1D_Pmag_vs_Wm1_unperturbed.bin
download_test_data con2prim/con2prim_tabulated_Newman1D_rho_vs_T_perturbed.bin
download_test_data con2prim/con2prim_tabulated_Newman1D_Pmag_vs_Wm1_perturbed.bin

download_test_data con2prim/con2prim_tabulated_Newman1D_entropy_rho_vs_T_unperturbed.bin
download_test_data con2prim/con2prim_tabulated_Newman1D_entropy_Pmag_vs_Wm1_unperturbed.bin
download_test_data con2prim/con2prim_tabulated_Newman1D_entropy_rho_vs_T_perturbed.bin
download_test_data con2prim/con2prim_tabulated_Newman1D_entropy_Pmag_vs_Wm1_perturbed.bin

download_test_data con2prim/con2prim_tabulated_Noble2D_rho_vs_T_unperturbed.bin
download_test_data con2prim/con2prim_tabulated_Noble2D_Pmag_vs_Wm1_unperturbed.bin
download_test_data con2prim/con2prim_tabulated_Noble2D_rho_vs_T_perturbed.bin
download_test_data con2prim/con2prim_tabulated_Noble2D_Pmag_vs_Wm1_perturbed.bin

code_error_workdir=$(mktemp -d)
created_directories+=("$code_error_workdir")
ln -s "$repo_root/SLy4_3335_rho391_temp163_ye66.h5" \
  "$code_error_workdir/SLy4_3335_rho391_temp163_ye66.h5"
for i in {0..90}; do
  if (cd "$code_error_workdir" && "$repo_root/test/unit_test_code_error" "$i"); then
    echo "Failed to fail!"
    exit 1
  else
    echo "Failed successfully!"
  fi
done

pyghl append SLy4_3335_rho391_temp163_ye66.h5
./test/unit_test_con2prim_tabulated SLy4_3335_rho391_temp163_ye66.h5 1

test_data_base_url="https://raw.githubusercontent.com/GRHayL/TestData/${et_legacy_testdata_ref}"
download_test_data induction/induction_interpolation_input.bin

download_test_data induction/induction_interpolation_ADM_input.bin
download_test_data induction/induction_interpolation_BSSN_input.bin

download_test_data induction/induction_interpolation_ccc_ADM_output.bin
download_test_data induction/induction_interpolation_ccc_ADM_output_pert.bin

download_test_data induction/induction_interpolation_vvv_ADM_output.bin
download_test_data induction/induction_interpolation_vvv_ADM_output_pert.bin

download_test_data induction/induction_interpolation_ccc_BSSN_output.bin
download_test_data induction/induction_interpolation_ccc_BSSN_output_pert.bin

./test/unit_test_induction_ccc_ADM
./test/unit_test_induction_vvv_ADM
./test/unit_test_induction_ccc_BSSN

download_test_data induction/HLL_flux_input.bin
download_test_data induction/HLL_flux_with_B_output.bin
download_test_data induction/HLL_flux_with_B_output_pert.bin
download_test_data induction/HLL_flux_with_Btilde_output.bin
download_test_data induction/HLL_flux_with_Btilde_output_pert.bin

./test/unit_test_HLL_flux

run_m1_tests
