// Thorn      : GRHayL
// File       : con2prim_unit_test.c
// Author(s)  : Leo Werneck & Samuel Cupp
// Description: In this file we provide an extensive unit test of
//              the Con2Prim gem.
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "unit_tests.h"

void unit_tests( CCTK_ARGUMENTS ) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  con2prim_unit_test();

  exit(1);
}
