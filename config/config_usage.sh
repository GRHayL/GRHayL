#!/bin/sh

usage_ubuntu()
{
    cat <<EOF
--------------------------------------
Tested GRHayL configurations on Ubuntu
--------------------------------------

Ubuntu 22.04
------------
# sudo apt-get update
# sudo apt-get install -y gcc libhdf5-serial-dev
  ./configure --prefix=.
  make
  make install
EOF
}

usage_mac()
{
    cat <<EOF
-------------------------------------
Tested GRHayL configurations on macOS
-------------------------------------

macOS 13.2.1 (Intel Core i7)
----------------------------
# brew install gcc hdf5
  ./configure --prefix=. --hdf5dir=\$(brew --prefix)/opt/hdf5
  make
  make install
EOF
}

usage_falcon()
{
    cat <<EOF
--------------------------------------
Tested GRHayL configurations on Falcon
--------------------------------------

Falcon - GCC (tested on 03-29-2023)
-----------------------------------
# module list
# Currently Loaded Modules:
#   1) gcc/12.1.0   2) mpich/3.4.3   3) hdf5/1.12.2
  CC=mpicc ./configure --prefix=. --hdf5dir=\$HDF5_ROOT
  make
  make install

Falcon - Intel (tested on 03-29-2023)
-------------------------------------
# module list
# Currently Loaded Modules:
#   1) intel/2021.4.0   2) mpich/3.4.3   3) hdf5/1.12.2
  CC=mpicc ./configure --prefix=. --hdf5dir=\$HDF5_ROOT
  make
  make install
EOF
}

usage_sawtooth()
{
    cat <<EOF
----------------------------------------
Tested GRHayL configurations on Sawtooth
----------------------------------------

Sawtooth - GCC (tested on 03-29-2023)
-------------------------------------
# module list
# Currently Loaded Modules:
#   1) gcc/9.2.0-gcc-4.8.5-bxc7   2) openmpi/4.0.2-gcc-9.2.0-cuda-10.1-5xzs
#   3) cmake/3.16.2-gcc-9.2.0-r3q3   4) hdf5/1.12.0_gcc9.2.0
  CC=gcc ./configure --prefix=. --hdf5dir=\$(module show hdf5 2>&1 | awk '/"PATH"/ { match($0, "/"); s=substr($0, RSTART); sub(/\/*bin")$/, "", s); print s }')
  make
  make install

Sawtooth - Intel (tested on 03-29-2023)
---------------------------------------
# module list
# Currently Loaded Modules:
#  1) zlib/1.2.11-intel-19.0.5-p7mu   2) intel/19.0.5-gcc-9.2.0-kl4p
#  3) intel-mpi/2018.5.288-intel-19.0.5-alnu   4) hdf5/1.10.6-intel-19.0.5-m6j7
  CC=mpiicc ./configure --prefix=. --hdf5dir=\$(module show hdf5 2>&1 | awk '/"PATH"/ { match($0, "/"); s=substr($0, RSTART); sub(/\/*bin")$/, "", s); print s }')
  make
  make install
EOF
}

usage_lemhi()
{
    cat <<EOF
-------------------------------------
Tested GRHayL configurations on Lemhi
-------------------------------------

Lemhi - GCC (tested on 03-29-2023)
----------------------------------
# module list
# Currently Loaded Modules:
#   1) gcc/11.2.0-gcc-9.4.0-mm7z   2) openjdk/11.0.12_7-gcc-11.2.0-gzj7
#   3) openmpi/4.1.1-gcc-11.2.0-nt2l   4) zlib/1.2.11-gcc-11.2.0-lnxa
#   5) hdf5/1.12.2-gcc-11.2.0-ej5f
  CC=gcc ./configure --prefix=. --hdf5dir=\$(module show hdf5 2>&1 | awk '/"PATH"/ { match($0, "/"); s=substr($0, RSTART); sub(/\/*bin")$/, "", s); print s }')
  make
  make install
EOF
}

usage_message()
{
    cat <<EOF
-----------
Basic usage
-----------

To compile and install GRHayL on your local machine, run:

  ./configure
  make
  make install

To install GRHayL on the current directory instead of /usr/local, replace the
first line above by one of the following three, which should all be equivalent:

  ./configure --prefix=.
  ./configure --prefix=`pwd`
  ./configure --prefix=$(pwd)

To change the compiler, set the CC environment variable before configuring, i.e.

  CC=mpiicc ./configure

or

  export CC=mpiicc
  ./configure

The tabulated EOS gem in GRHayL depends on HDF5. By default, the configuration
script will attempt to find the HDF5 installation through pkg-config. If it
cannot be found, then an error will be produced. If HDF5 is installed in your
system but pkg-config cannot find it, you can specify the HDF5 installation
directory using

  ./configure --hdf5dir=<hdf5_directory>

Alternatively, if your installation is non-standard, you can manually specify
the include and lib directories of your HDF5 installation using

  ./configure --hdf5inc=<hdf5_include_directory> --hdf5lib=<hdf5_lib_directory>

If you are interested in hybrid EOS only, the HDF5 dependency might be
undesirable. To disable HDF5, simply run

  ./configure --nohdf5.

-----------------------
Specific configurations
-----------------------

To see specific configurations, please run one of the following commands:

  ./configure --usage ubuntu
  ./configure --usage mac
  ./configure --usage falcon
  ./configure --usage sawtooth
  ./configure --usage lemhi

In all cases you might want to run:

  export CFLAGS="\$CFLAGS -I$(pwd)/include"
  export LD_FLAGS="\$LD_FLAGS -L$(pwd)/lib"
  export LD_LIBRARY_PATH="\$LD_LIBRARY_PATH:\$(pwd)/lib"

so that compilers can find the GRHayL headers and so that executables that use
the GRHayL library are able to run without errors.
EOF
}
