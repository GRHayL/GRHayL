export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`/lib

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_conservs_input.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_conservs_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_conservs_output_pert.bin

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_primitives_input.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_primitives_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_primitives_output_pert.bin

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_induction_gauge_rhs_input.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_induction_gauge_rhs_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_induction_gauge_rhs_output_pert.bin

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_HLL_flux_input.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_HLL_flux_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_HLL_flux_output_pert.bin

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_reconstruction_input.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_reconstruction_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_reconstruction_output_pert.bin

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_flux_source_input.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_flux_source_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/ET_Legacy/ET_Legacy_flux_source_output_pert.bin

./test/unit_test_ET_Legacy_conservs
./test/unit_test_ET_Legacy_primitives
./test/unit_test_ET_Legacy_induction_gauge_rhs
./test/unit_test_ET_Legacy_HLL_flux
./test/unit_test_ET_Legacy_reconstruction
./test/unit_test_ET_Legacy_flux_source

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/metric_Bfield_initial_data.bin

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/apply_conservative_limits_input.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/apply_conservative_limits_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/apply_conservative_limits_output_pert.bin

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_multi_method_hybrid_input.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_multi_method_hybrid_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_multi_method_hybrid_output_pert.bin

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/enforce_primitive_limits_and_compute_u0_input.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/enforce_primitive_limits_and_compute_u0_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/enforce_primitive_limits_and_compute_u0_output_pert.bin

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/compute_conservs_and_Tmunu_input.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/compute_conservs_and_Tmunu_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/compute_conservs_and_Tmunu_output_pert.bin

./test/unit_test_apply_conservative_limits
./test/unit_test_con2prim_multi_method_hybrid
./test/unit_test_enforce_primitive_limits_and_compute_u0
./test/unit_test_compute_conservs_and_Tmunu

./test/unit_test_hybrid_failure

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/EOS/simple_table.h5
./test/unit_test_tabulated_eos simple_table.h5

./test/unit_test_piecewise_polytrope

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/grhayl_core/grhayl_core_test_suite_input.bin
./test/unit_test_grhayl_core_test_suite

curl -O https://stellarcollapse.org/EOS/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2 && bunzip2 LS220_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/flux_source/hybrid_flux_input.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/flux_source/hybrid_flux_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/flux_source/hybrid_flux_output_pert.bin

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/flux_source/tabulated_flux_input.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/flux_source/tabulated_flux_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/flux_source/tabulated_flux_output_pert.bin

./test/unit_test_hybrid_flux
./test/unit_test_tabulated_flux

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/reconstruction/PLM_reconstruction_input.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/reconstruction/PLM_reconstruction_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/reconstruction/PLM_reconstruction_output_pert.bin

./test/unit_test_PLM_reconstruction

curl -O https://stellarcollapse.org/EOS/SLy4_3335_rho391_temp163_ye66.h5.bz2 && bunzip2 SLy4_3335_rho391_temp163_ye66.h5.bz2

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/Neutrinos/nrpyleakage_optically_thin_gas_unperturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/Neutrinos/nrpyleakage_optically_thin_gas_perturbed.bin

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/Neutrinos/nrpyleakage_constant_density_sphere_unperturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/Neutrinos/nrpyleakage_constant_density_sphere_perturbed.bin

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/Neutrinos/nrpyleakage_luminosities_unperturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/Neutrinos/nrpyleakage_luminosities_perturbed.bin

./test/unit_test_nrpyleakage_optically_thin_gas SLy4_3335_rho391_temp163_ye66.h5 1
./test/unit_test_nrpyleakage_constant_density_sphere SLy4_3335_rho391_temp163_ye66.h5 1
./test/unit_test_nrpyleakage_luminosities SLy4_3335_rho391_temp163_ye66.h5 1

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_tabulated_Palenzuela1D_rho_vs_T_unperturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_tabulated_Palenzuela1D_Pmag_vs_Wm1_unperturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_tabulated_Palenzuela1D_rho_vs_T_perturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_tabulated_Palenzuela1D_Pmag_vs_Wm1_perturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_tabulated_Palenzuela1D_entropy_rho_vs_T_unperturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_tabulated_Palenzuela1D_entropy_Pmag_vs_Wm1_unperturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_tabulated_Palenzuela1D_entropy_rho_vs_T_perturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_tabulated_Palenzuela1D_entropy_Pmag_vs_Wm1_perturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_tabulated_Newman1D_rho_vs_T_unperturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_tabulated_Newman1D_Pmag_vs_Wm1_unperturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_tabulated_Newman1D_rho_vs_T_perturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_tabulated_Newman1D_Pmag_vs_Wm1_perturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_tabulated_Newman1D_entropy_rho_vs_T_unperturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_tabulated_Newman1D_entropy_Pmag_vs_Wm1_unperturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_tabulated_Newman1D_entropy_rho_vs_T_perturbed.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/con2prim/con2prim_tabulated_Newman1D_entropy_Pmag_vs_Wm1_perturbed.bin
./test/unit_test_con2prim_tabulated SLy4_3335_rho391_temp163_ye66.h5 1

for i in {0..31}
do
  ./Unit_Tests/failure.sh "./test/unit_test_code_error "+$i
done

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/induction/induction_interpolation_input.bin

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/induction/induction_interpolation_ADM_input.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/induction/induction_interpolation_BSSN_input.bin

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/induction/induction_interpolation_ccc_ADM_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/induction/induction_interpolation_ccc_ADM_output_pert.bin

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/induction/induction_interpolation_vvv_ADM_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/induction/induction_interpolation_vvv_ADM_output_pert.bin

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/induction/induction_interpolation_ccc_BSSN_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/induction/induction_interpolation_ccc_BSSN_output_pert.bin

./test/unit_test_induction_ccc_ADM
./test/unit_test_induction_vvv_ADM
./test/unit_test_induction_ccc_BSSN

curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/induction/HLL_2D_flux_input.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/induction/HLL_2D_flux_with_B_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/induction/HLL_2D_flux_with_B_output_pert.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/induction/HLL_2D_flux_with_Btilde_output.bin
curl -O https://raw.githubusercontent.com/GRHayL/TestData/main/induction/HLL_2D_flux_with_Btilde_output_pert.bin

./test/unit_test_HLL_2D_flux

rm *.bin
rm *.h5
rm *.bz2
