#include "nrpyeos_hybrid.h"

/*
 * Function     : NRPyEOS_find_polytropic_index_from_P()
 * Description  : For a given value of pressure, find the appropriate polytropic index;
 *                usually aliased as ghl_hybrid_find_polytropic_index
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_hybrid_find_polytropic_index_from_P
*/

int NRPyEOS_find_polytropic_index_from_P(
    const ghl_eos_parameters *restrict eos,
    const double P_in) {

/* We want to find the appropriate polytropic EOS for the
 * input value P_in. Remember that:
 *
 * if P < P_{0}:                P_{0} , index: 0
 * if P >= P_{0} but < P_{1}: P_{1} , index: 1
 * if P >= P_{1} but < P_{2}: P_{2} , index: 2
 *                      ...
 * if P >= P_{j-1} but < P_{j}: P_{j} , index: j
 *
 * Then, a simple way of determining the index is through
 * the formula:
 *  ---------------------------------------------------------------------------
 * | index = (P >= P_{0}) + (P >= P_{1}) + ... + (P >= P_{neos-2}) |
 *  ---------------------------------------------------------------------------
 */
if(eos->neos == 1) return 0;

int polytropic_index = 0;
for(int j=0; j <= eos->neos-2; j++) polytropic_index += (P_in >= eos->p_ppoly[j]);

return polytropic_index;
}