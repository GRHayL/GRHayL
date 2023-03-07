common_sources = files('io_conservative_binary.c',
                       'io_flux_source_read_from_binary_file.c',
                       'io_metric_binary.c',
                       'io_primitive_binary.c',
                       'io_stress_energy_binary.c',
                       'randomize_metric.c',
                       'randomize_primitives.c',
                       'validate_conservatives.c',
                       'validate_primitives.c',
                       'validate_stress_energy.c')

dep_omp = dependency('openmp', required : false)
dep_without_hdf5 = [dep_grhayl, dep_math, dep_omp]
dep_with_hdf5 = [dep_grhayl, dep_math, dep_omp, dep_hdf5]

unit_test_apply_inequality_fixes = executable(
  'unit_test_apply_inequality_fixes',
  sources : ['unit_test_apply_inequality_fixes.c', common_sources],
  dependencies : dep_without_hdf5)
test('apply_inequality_fixes', unit_test_apply_inequality_fixes, args : ['myfile'])

unit_test_compute_conservs_and_Tmunu = executable(
  'unit_test_compute_conservs_and_Tmunu',
  sources : ['unit_test_compute_conservs_and_Tmunu.c', common_sources],
  dependencies : dep_without_hdf5)
test('compute_conservs_and_Tmunu', unit_test_compute_conservs_and_Tmunu)

unit_test_enforce_primitive_limits_and_compute_u0 = executable(
  'unit_test_enforce_primitive_limits_and_compute_u0',
  sources : ['unit_test_enforce_primitive_limits_and_compute_u0.c', common_sources],
  dependencies : dep_without_hdf5)
test('enforce_primitive_limits_and_compute_u0', unit_test_enforce_primitive_limits_and_compute_u0)

unit_test_flux_source = executable(
  'unit_test_flux_source',
  sources : ['unit_test_flux_source.c', common_sources],
  dependencies : dep_without_hdf5)
test('flux_source', unit_test_flux_source)

unit_test_font_fix = executable(
  'unit_test_font_fix',
  sources : ['unit_test_font_fix.c', common_sources],
  dependencies : dep_without_hdf5)
test('font_fix', unit_test_font_fix)

unit_test_ET_Legacy_HLL_flux = executable(
  'unit_test_ET_Legacy_HLL_flux',
  sources : ['unit_test_ET_Legacy_HLL_flux.c', common_sources],
  dependencies : dep_without_hdf5)
test('ET_Legacy_HLL_flux', unit_test_ET_Legacy_HLL_flux)

unit_test_Noble2D_Hybrid_Multi_Method = executable(
  'unit_test_Noble2D_Hybrid_Multi_Method',
  sources : ['unit_test_Noble2D_Hybrid_Multi_Method.c', common_sources],
  dependencies : dep_without_hdf5)
test('Noble2D_Hybrid_Multi_Method', unit_test_Noble2D_Hybrid_Multi_Method)

unit_test_ET_Legacy_induction_gauge_rhs = executable(
  'unit_test_ET_Legacy_induction_gauge_rhs',
  sources : ['unit_test_ET_Legacy_induction_gauge_rhs.c', common_sources],
  dependencies : dep_without_hdf5)
test('ET_Legacy_induction_gauge_rhs', unit_test_ET_Legacy_induction_gauge_rhs)

unit_test_ET_Legacy_reconstruction = executable(
  'unit_test_ET_Legacy_reconstruction',
  sources : ['unit_test_ET_Legacy_reconstruction.c', common_sources],
  dependencies : dep_without_hdf5)
test('ET_Legacy_reconstruction', unit_test_ET_Legacy_reconstruction)

unit_test_ET_Legacy_flux_source = executable(
  'unit_test_ET_Legacy_flux_source',
  sources : ['unit_test_ET_Legacy_flux_source.c', common_sources],
  dependencies : dep_without_hdf5)
test('ET_Legacy_flux_source', unit_test_ET_Legacy_flux_source)

unit_test_tabulated_eos = executable(
  'unit_test_tabulated_eos',
  sources : ['unit_test_tabulated_eos.c', 'tabulated_eos_unit_test_helpers.c', common_sources],
  dependencies : dep_with_hdf5)
test('tabulated_eos', unit_test_tabulated_eos, args : './simple_table.h5')
