#include "ghl.h"

/* Function    : ghl_initialize_conservatives()
 * Description : Initialize the conservative struct from user input
 *
 * Inputs      : rho            - value of rho_star (densitized density)
 *             : tau            - value of \tilde{tau} (densitized energy
 *                                variable)
 *             : S_x            - value of \tilde{S}_1 (densitized momentum variable)
 *             : S_y            - value of \tilde{S}_2 (densitized momentum variable)
 *             : S_z            - value of \tilde{S}_3 (densitized momentum variable)
 *             : entropy        - value of densitized entropy
 *
 * Outputs     : cons           - returns ghl_conservative_quantities struct containing the input
 */

void ghl_initialize_conservatives(
             const double rho, const double tau,
             const double S_x, const double S_y, const double S_z,
             const double Y_e, const double entropy,
             ghl_conservative_quantities *restrict cons) {

  cons->rho = rho;
  cons->SD[0] = S_x;
  cons->SD[1] = S_y;
  cons->SD[2] = S_z;
  cons->tau = tau;
  cons->Y_e = Y_e;
  cons->entropy = entropy;
}
