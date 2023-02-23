# GRHayL

---

[![Ubuntu-gcc](https://github.com/SamuelCupp/GRHayL_beta/actions/workflows/github-actions-Ubuntu-gcc.yml/badge.svg)](https://github.com/SamuelCupp/GRHayL_beta/actions/workflows/github-actions-Ubuntu-gcc.yml)
[![Ubuntu-intel](https://github.com/SamuelCupp/GRHayL_beta/actions/workflows/github-actions-Ubuntu-intel.yml/badge.svg)](https://github.com/SamuelCupp/GRHayL_beta/actions/workflows/github-actions-Ubuntu-intel.yml)
[![Ubuntu-clang](https://github.com/SamuelCupp/GRHayL_beta/actions/workflows/github-actions-Ubuntu-clang.yml/badge.svg)](https://github.com/SamuelCupp/GRHayL_beta/actions/workflows/github-actions-Ubuntu-clang.yml)
[![MacOS-gcc](https://github.com/SamuelCupp/GRHayL_beta/actions/workflows/github-actions-MacOS-gcc.yml/badge.svg)](https://github.com/SamuelCupp/GRHayL_beta/actions/workflows/github-actions-MacOS-gcc.yml)
[![MacOS-clang](https://github.com/SamuelCupp/GRHayL_beta/actions/workflows/github-actions-MacOS-clang.yml/badge.svg)](https://github.com/SamuelCupp/GRHayL_beta/actions/workflows/github-actions-MacOS-clang.yml)

This is a beta version of the infrastructure-agnostic magnetohydrodynamics code
my collaborators and I are developing. This is a piece of the upcoming GRHayL
(General Relativistic Hydrodynamic Library) code. This repo currently contains
the core GRHayL "chalice" which provides the basic C structs and other fully
independent functions that act only on these structs.

In addtion, this is the current home of the Con2Prim "gem", which is designed to
simply take primitives and conservative structs, perform c2p routines, and
return the structs with the reconstructed primitive data. This code is being
developed with the goal of converting IllinoisGRMHD into a set of libraries
which are extensible and easy to plug-n-play with other codes.

In addition to the chalice and C2P gem, some elements of the EOS gem are also
present to facilitate the Con2Prim routines. This code does *not* depend on the
Einstein Toolkit or Cactus framework, though this repository is currently being
tested within that infrastructure. All GRHayL code except the included
con2prim_test_suite.c have no references to Cactus parameters or header files.

The test code provides a clean example of how to use the library's
functions. After the initial randomized data is set up and copies of the
original primitives are set, the remaining code is the C2P routine in its
entirety.

## Building and Installing `GRHayL`

`GRHayL` uses the [meson build system](https://mesonbuild.com). For the
convenience of users we package `GRHayL` with a `configure` script wrapper,
which is a `Python3` script that handles the `meson` configuration. The builds
happen out of the source tree, with all object files placed in the `build/`
directory (by default).

### Dependencies

Currently, `GRHayL` requires `HDF5` C bindings in order to compile. On Ubuntu,
you can install `HDF5` using

```shell
$ sudo apt-get install libhdf5-serial-dev
```

`GRHayL` is also actively tested on MacOS using the [Homebrew](https://brew.sh/)
package manager. To install `HDF5` on MacOS we thus recommend using

```shell
$ brew install hdf5
```

We have not tested `GRHayL` with [macports](https://www.macports.org/), although
we do not expect users to have issues when using it.

The `HDF5` dependency arises from the Tabulated Equation of State (EOS) gem, as
it currently only supports EOS tables in that format. We are working on adding
an option to compile `GRHayL` with Hybrid EOS support only, thus removing the
`HDF5` dependency for those who do not wish to use Tabulated EOSs. We are also
working supporting tables in the [CompOSE](https://compose.obspm.fr/table)
format, which would also remove the `HDF5` dependency for those who wish to use
these larger ASCII tables.

### System-wide Installation (default)

As an example, suppose you wish to install `GRHayL` on your system using
`gcc`. To do so, run the following commands:

```shell
$ CC=gcc ./configure
$ make
$ make install
```

### Local Installation

As an example, a local installation that uses the next-generation intel
compilers is achieved as follows:

```shell
$ CC=icx ./configure --prefix=<local_path>
$ make
$ make install
```

### Legacy Development Build System

`GRHayL` currently contains an "old" `Makefile` that is still used by its
developers, but is likely to be deprecated and/or removed in the
future. Although this `Makefile` does take *some* measures that allow `GRHayL`
to be compiled on a few different systems, it is *not* the recommended way to
compile the library. One cannot specify the installation (`./lib`) and build
(`./build`) paths when using this `Makefile`, and it also does not produce an
`include` directory alongside the `lib` directory. To use it, simply run

```shell
$ make
```

without running `./configure` first. When `./configure` is run, the `Makefile`
is replaced by a different one that makes use of the `meson` build system. To
restore the old `Makefile` and clean up all files produced by the `meson` build,
use

```shell
$ make realclean
```
