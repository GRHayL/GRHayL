#include "ghl.h"
#include "ghl_radiation.h"

// This is the global will generally point to a closure of choice.
double (*ghl_m1_closure)(double) = NULL;

// Initialize m1_closure
ghl_error_codes_t ghl_initialize_m1_closure(ghl_m1_closure_t ghl_m1_closure_type) {
  switch(ghl_m1_closure_type) {
    case Eddington:
      ghl_m1_closure = eddington;
      break;
    case Kershaw:
      ghl_m1_closure = kershaw;
      break;
    case Minerbo:
      ghl_m1_closure = minerbo;
      break;
    case Thin:
      ghl_m1_closure = thin;
      break;
    default:
      return ghl_error_invalid_m1_closure;
  }
  return ghl_success;
}
