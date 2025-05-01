#include "ghl.h"

/**
 * @ingroup Con2Prim
 * @brief Gets @ref Con2Prim routine name from method key
 *
 * @details
 * This function returns the name of a conservative-to-primitive solver method
 * given a @ref ghl_con2prim_method_t key. If key matches nothing in the enum,
 * then the code will fail with a message giving the value of the failed key.
 *
 * @param[in] key: key from @ref ghl_con2prim_method_t selecting a method
 *
 * @returns the name of the conservative-to-primitive solver method
 */
char *ghl_get_con2prim_routine_name(const ghl_con2prim_method_t key) {
  switch(key) {
    case None:
      return "None";
    case Noble2D:
      return "Noble2D";
    case Noble1D:
      return "Noble1D";
    case Noble1D_entropy:
      return "Noble1D_entropy";
    case Font1D:
      return "Font1D";
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
