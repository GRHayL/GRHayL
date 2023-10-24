#include "ghl.h"

void ghl_enforce_detgtij_and_initialize_ADM_metric(
      const double lapse,
      const double betax,
      const double betay,
      const double betaz,
      const double gxx,
      const double gxy,
      const double gxz,
      const double gyy,
      const double gyz,
      const double gzz,
      ghl_metric_quantities *restrict ADM_metric) {

  ghl_initialize_metric(
        lapse, betax, betay, betaz,
        gxx, gxy, gxz,
        gyy, gyz, gzz,
        ADM_metric);

  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(ADM_metric, &metric_aux);
  const double phi = (1.0/12.0) * log(ADM_metric->detgamma);

  const double psi2 = exp(2.0*phi);
  const double psi4 = SQR(psi2);
  const double psi4inv = 1.0/psi4;


        /**********************************************************************
         * Compute \tilde{\gamma_{ij}}, phi, and psi (BSSN) from g_{ij} (ADM) *
         **********************************************************************/

  const double gtxx = gxx*psi4inv;
  const double gtxy = gxy*psi4inv;
  const double gtxz = gxz*psi4inv;
  const double gtyy = gyy*psi4inv;
  const double gtyz = gyz*psi4inv;
  const double gtzz = gzz*psi4inv;

  /*********************************
   * Apply det gtij = 1 constraint *
   *********************************/
  const double gtijdet = gtxx * gtyy * gtzz
                       + gtxy * gtyz * gtxz
                       + gtxz * gtxy * gtyz
                       - gtxz * gtyy * gtxz
                       - gtxy * gtxy * gtzz
                       - gtxx * gtyz * gtyz;

  const double gtijdet_Fm1o3 = fabs(1.0/cbrt(gtijdet));

  const double gxx_new = gtxx * gtijdet_Fm1o3 * psi4;
  const double gxy_new = gtxy * gtijdet_Fm1o3 * psi4;
  const double gxz_new = gtxz * gtijdet_Fm1o3 * psi4;
  const double gyy_new = gtyy * gtijdet_Fm1o3 * psi4;
  const double gyz_new = gtyz * gtijdet_Fm1o3 * psi4;
  const double gzz_new = gtzz * gtijdet_Fm1o3 * psi4;

  ghl_initialize_metric(lapse,
                    betax, betay, betaz,
                    gxx_new, gxy_new, gxz_new,
                    gyy_new, gyz_new, gzz_new,
                    ADM_metric);

  if(gtijdet<0.0) ghl_warn(
                      "WARNING: det[3-metric]<0.0. Hopefully this is occurring in gz's! "
                      "gtij_phys = %.2e %.2e %.2e %.2e %.2e %.2e gtij_new = %.2e %.2e %.2e %.2e %.2e %.2e | gijdet = %.2e | gtijdet = %.2e\n",
  			     gxx, gxy, gxz, gyy, gyz, gzz, gtxx, gtxy, gtxz, gtyy, gtyz, gtzz, ADM_metric->detgamma, gtijdet);
}
