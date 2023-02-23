# List each gem as a subdirectory
subdir('Con2Prim')
subdir('EOS')
subdir('Flux_Source')
subdir('GRHayL_Core')
subdir('Induction')
subdir('Neutrinos')
subdir('Reconstruction')

# List the include directory as a subdirectory
subdir('include')

# The variables below are defined inside each gem's subdirectory
sources_lib = [Con2Prim_sources,
               EOS_sources,
               Flux_Source_sources,
               GRHayL_Core_sources,
               Induction_sources,
               Neutrinos_sources,
               Reconstruction_sources]

# The include variable is defined in the include/ subdirectory
include_lib = [GRHayL_include_dir]

# The headers variable is defined in the include/ subdirectory
headers_lib = [headers]

# Define this as the base directory for the GRHayL library
base_lib_dir = meson.current_source_dir()

dep_extern  = [dep_hdf5, dep_math]

# Now define the GRHayL library file
lib_grhayl  = library('grhayl',
                      sources_lib,
                      include_directories : include_lib,
                      implicit_include_directories : true,
                      dependencies : dep_extern,
                      install : true)

# This dependency is used by the unit tests
dep_grhayl = declare_dependency(include_directories : include_lib,
                                link_with : lib_grhayl)
