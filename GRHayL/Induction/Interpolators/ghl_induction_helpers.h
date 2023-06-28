#ifndef INDUCTION_HELPERS_H_
#define INDUCTION_HELPERS_H_

// This variable size stuff only works in C,
// not C++. D:

void ghl_A_i_avg(
      const int size,
      const double Ax_stencil[size][size][size],
      const double Ay_stencil[size][size][size],
      const double Az_stencil[size][size][size],
      double A_to_phitilde[3],
      double A_to_Ax[3],
      double A_to_Ay[3],
      double A_to_Az[3]);

void ghl_BSSN_cell_interp(
      const int size,
      const metric_quantities metric_stencil[size][size][size],
      const double psi_stencil[size][size][size],
      metric_quantities *restrict metric_interp,
      double lapse_psi2_interp[3],
      double *restrict lapse_over_psi6_interp);

void ghl_ADM_cell_interp(
      const int size,
      const metric_quantities metric_stencil[size][size][size],
      metric_quantities *restrict metric_interp,
      double *restrict lapse_over_psi6_interp);

void ghl_ADM_vertex_interp(
      const int size,
      const metric_quantities metric_stencil[size][size][size],
      double gammaUU_interp[3][3]);

#endif // INDUCTION_HELPERS_H_
