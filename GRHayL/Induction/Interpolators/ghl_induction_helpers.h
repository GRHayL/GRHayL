#ifndef INDUCTION_HELPERS_H_
#define INDUCTION_HELPERS_H_

void ghl_A_i_avg(
      const double Ax_stencil[3][3][3],
      const double Ay_stencil[3][3][3],
      const double Az_stencil[3][3][3],
      double A_to_phitilde[3],
      double A_to_Ax[3],
      double A_to_Ay[3],
      double A_to_Az[3]);

void ghl_BSSN_cell_interp(
      const ghl_metric_quantities metric_stencil[2][2][2],
      const double psi_stencil[2][2][2],
      ghl_metric_quantities *restrict metric_interp,
      double lapse_psi2_interp[3],
      double *restrict lapse_over_psi6_interp);

void ghl_ADM_cell_interp(
      const ghl_metric_quantities metric_stencil[2][2][2],
      ghl_metric_quantities *restrict metric_interp,
      double *restrict lapse_over_psi6_interp);

void ghl_ADM_vertex_interp(
      const ghl_metric_quantities metric_stencil[2][2][2],
      double gammaUU_interp[3][3]);

#endif // INDUCTION_HELPERS_H_
