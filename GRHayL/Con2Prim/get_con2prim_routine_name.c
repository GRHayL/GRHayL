#include "ghl.h"

/**
 * @ingroup Con2Prim
 * @brief Gets @ref Con2Prim routine name from method key
 *
 * @details
 * This function returns the name of a conservative-to-primitive solver method
 * given a @ref ghl_con2prim_id_t key. If key matches nothing in the enum,
 * then the function returns NULL.
 *
 * @param[in] key key from @ref ghl_con2prim_id_t selecting a method
 *
 * @returns the name of the conservative-to-primitive solver method
 */
const char *ghl_get_con2prim_routine_name(const ghl_con2prim_id_t key) {
  switch(key) {
    case ghl_con2prim_id_None:
      return "None";
    case ghl_con2prim_id_Noble2D:
      return "Noble2D";
    case ghl_con2prim_id_Noble1D:
      return "Noble1D";
    case ghl_con2prim_id_Noble1D_entropy:
      return "Noble1D_entropy";
    case ghl_con2prim_id_Noble1D_entropy2:
      return "Noble1D_entropy2";
    case ghl_con2prim_id_Font1D:
      return "Font1D";
    case ghl_con2prim_id_Palenzuela1D:
      return "Palenzuela1D";
    case ghl_con2prim_id_Palenzuela1D_entropy:
      return "Palenzuela1D_entropy";
    case ghl_con2prim_id_Newman1D:
      return "Newman1D";
    case ghl_con2prim_id_Newman1D_entropy:
      return "Newman1D_entropy";
    default:
      return NULL;
  }
}
