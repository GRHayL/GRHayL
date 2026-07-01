#!/bin/bash

set -Eeuxo pipefail

./configure -r
make tests

LD_LIBRARY_PATH="$(pwd)/build/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH

echo $LD_LIBRARY_PATH

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

  curl -O "$url"
}

test_data_base_url="https://raw.githubusercontent.com/GRHayL/TestData/main"
download_test_data() {
  filepath="$1"
  url="${test_data_base_url}/${filepath}"
  download_file "$url"
}

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
if [ -f LS220_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2 ]; then
  bunzip2 LS220_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2
fi

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
if [ -f SLy4_3335_rho391_temp163_ye66.h5.bz2 ]; then
  bunzip2 SLy4_3335_rho391_temp163_ye66.h5.bz2
fi

download_test_data Neutrinos/nrpyleakage_optically_thin_gas_unperturbed.bin
download_test_data Neutrinos/nrpyleakage_optically_thin_gas_perturbed.bin

download_test_data Neutrinos/nrpyleakage_constant_density_sphere_unperturbed.bin
download_test_data Neutrinos/nrpyleakage_constant_density_sphere_perturbed.bin

download_test_data Neutrinos/nrpyleakage_luminosities_unperturbed.bin
download_test_data Neutrinos/nrpyleakage_luminosities_perturbed.bin

./test/unit_test_nrpyleakage_optically_thin_gas SLy4_3335_rho391_temp163_ye66.h5 1
./test/unit_test_nrpyleakage_constant_density_sphere SLy4_3335_rho391_temp163_ye66.h5 1
./test/unit_test_nrpyleakage_luminosities SLy4_3335_rho391_temp163_ye66.h5 1

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

for i in {0..85}; do
  if ./test/unit_test_code_error "$i"; then
    echo "Failed to fail!"
    exit 1
  else
    echo "Failed successfully!"
  fi
done

pyghl append SLy4_3335_rho391_temp163_ye66.h5
./test/unit_test_con2prim_tabulated SLy4_3335_rho391_temp163_ye66.h5 1

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

rm -f ./*.bin ./*.h5 ./*.bz2
