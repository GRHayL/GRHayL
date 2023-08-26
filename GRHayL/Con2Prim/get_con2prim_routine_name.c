#include "ghl.h"

char *ghl_get_con2prim_routine_name(const ghl_con2prim_method_t key) {
  switch(key) {
    case Noble2D:
      return "Noble2D";
    case Noble1D:
      return "Noble1D";
    case Noble1D_entropy:
      return "Noble1D_entropy";
    case Noble1D_entropy2:
      return "Noble1D_entropy2";
    case Font1D:
      return "Font1D";
    case CerdaDuran2D:
      return "CerdaDuran2D";
    case CerdaDuran3D:
      return "CerdaDuran3D";
    case Palenzuela1D:
      return "Palenzuela1D";
    case Palenzuela1D_entropy:
      return "Palenzuela1D_entropy";
    case Newman1D:
      return "Newman1D";
    case Newman1D_entropy:
      return "Newman1D_entropy";
    default:
      ghl_error("Unknown routine key %d\n", key);
  }
  return NULL;
}
