#include "ghl_nrpyeos_tabulated.h"

int NRPyEOS_tabulated_get_index_T(
    const ghl_eos_parameters *restrict eos,
    const double T) {
  const double lt = log(T);
  const double ltmin = eos->table_logT[0];
  const double ltmax = eos->table_logT[eos->N_T - 1];
  if(lt < ltmin || lt > ltmax) {
    return -1;
  }
  const double dlt = eos->table_logT[1] - ltmin;
  return (lt - ltmin) / dlt;
}
