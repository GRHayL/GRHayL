# GRHayL

---

[![Ubuntu-gcc](https://github.com/SamuelCupp/GRHayL_beta/actions/workflows/github-actions-Ubuntu.yml/badge.svg)](https://github.com/SamuelCupp/GRHayL_beta/actions/workflows/github-actions-Ubuntu.yml)
[![MacOS-gcc](https://github.com/SamuelCupp/GRHayL_beta/actions/workflows/github-actions-MacOS-gcc.yml/badge.svg)](https://github.com/SamuelCupp/GRHayL_beta/actions/workflows/github-actions-MacOS-gcc.yml)
[![MacOS-clang](https://github.com/SamuelCupp/GRHayL_beta/actions/workflows/github-actions-MacOS-clang.yml/badge.svg)](https://github.com/SamuelCupp/GRHayL_beta/actions/workflows/github-actions-MacOS-clang.yml)

This is a beta version of the infrastructure-agnostic magnetohydrodynamics code
my collaborators and I are developing. This is a piece of the upcoming GRHayL
(General Relativistic Hydrodynamic Library) code. This repo currently contains the core
GRHayL "chalice" which provides the basic C structs and other fully independent functions
that act only on these structs.

In addtion, this is the current home of the Con2Prim "gem", which is designed to simply
take primitives and conservative structs, perform c2p routines, and return the structs
with the reconstructed primitive data. This code is being developed with the  goal of
converting IllinoisGRMHD into a set of libraries which are extensible and easy to
plug-n-play with other codes.

In addition to the chalice and C2P gem, some elements of the EOS gem are also present
to facilitate the Con2Prim routines. This code does *not* depend on the Einstein
Toolkit or Cactus framework, though this repository is currently being tested within
that infrastructure. All GRHayL code except the included con2prim_test_suite.c have no
references to Cactus parameters or header files.

The test code provides a clean example of how to use the library's functions. After
the initial randomized data is set up and copies of the original primitives are set,
the remaining code is the C2P routine in its entirety.
