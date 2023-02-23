project( 'GRHayL',
         'c',
         version : '0.0.0',
         license : 'BSD-2-Clause license',
         default_options : ['buildtype=plain',
                            'c_args=-O2 -march=native',
                            'c_std=c99'])
cc = meson.get_compiler('c')

dep_math = cc.find_library('m', required : true)
dep_hdf5 = dependency('hdf5', language: 'c')

subdir('GRHayL')
subdir('Unit_Tests')
