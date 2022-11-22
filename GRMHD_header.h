#ifndef GRMHD_HEADER_H_
#define GRMHD_HEADER_H_

#include "math.h"
#include "stdbool.h"

#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define SQR(x) ((x) * (x))
#define ONE_OVER_SQRT_4PI 0.282094791773878143474039725780

//------------------ temporaries -------------------
/* This section should be housed somewhere higher up in the hierarchy and not defined
   in con2prim. The parameters, particularly, should be defined once and kept static. */

/*
   The struct GRMHD_parameters contains parameters for controlling
   the behavior of the con2prim routines. The struct elements are detailed below:

 --main_routine: selects which con2prim routine to use. The
   available con2prim routines are given by the enum con2prim_routines.

 --backup_routine[3]: stores up to three backup routines using the values
   in con2prim_routines. The "None" option is provided if backups are not
   desired or if there are less than 3 backup routines.

 --evolve_entropy: tells the code whether or not entropy is being evolved.
   If true, then the code will set the conservative entropy when needed.

 --calc_prim_guess: currently does nothing. However, some codes might be able to
   directly provide guesses for the initial primitive values for con2prim. If they
   can, this should be supported.
*/

typedef struct GRMHD_parameters {
  int main_routine, backup_routine[3];
  bool evolve_entropy;
  bool evolve_temp;
  bool calc_prim_guess;
  double psi6threshold;
  bool update_Tmunu;
} GRMHD_parameters;

/* The struct metric_quantities contains variables for storing the (point-wise)
   metric data. The struct elements are detailed below:

 --bssn_phi: TODO

 --bssn_psi: TODO

 --bssn_gij: the BSSN conformal metric g_{i j}. These variables are set via the
   initialize_metric function.

 --betai: the shift vector (beta^x, beta^y, beta^z) from the 3+1
   decomposition. These variables are set via the initialize_metric function.

 --lapse: the lapse from the 3+1 decomposition. These variables are set via the
   initialize_metric function.

   All of the following varibles are calculated from the previous quantities in the initialize_metric
   function.

 --bssn_gupij: the BSSN conformal metric inverse g^{i j}. They are calculated from
   the bssn_gij variables using the relations
     gupxx =   ( gyy * gzz - gyz * gyz );
     gupxy = - ( gxy * gzz - gyz * gxz );
     gupxz =   ( gxy * gyz - gyy * gxz );
     gupyy =   ( gxx * gzz - gxz * gxz );
     gupyz = - ( gxx * gyz - gxy * gxz );
     gupzz =   ( gxx * gyy - gxy * gxy );

 --lapseinv: this variable stores the inverse of the lapse.

 --lapm1: this variable stores the quantity (lapse-1).

 --psi2: this variable stores the quantity exp(2.0*metric.bssn_phi).

 --psi4: this variable stores the quantity (psi2^2).

 --psi6: this variable stores the quantity (psi4*psi2).

 --psi4inv: this variable stores the inverse of psi4.

 --lapseinv2: this variables stores the quantity (lapse^(-2))

 --adm_gij: the ADM metric g_{i j} They are computed by the initialize_metric
   function using the relation adm_gij = psi4*bssn_gij.

 --adm_gupij: the ADM metric inverse g^{i j}. They are computed by the
   initialize_metric function using the relation adm_gupij = psi4inv*bssn_gupij.

 --g4dn: the 4-metric g_{\mu \nu}. This quantity is needed for computing T_{\mu \nu} and T^{\mu \nu}
   and the HARM con2prim lowlevel functions.

 --g4up: the 4-metric inverse g^{\mu \nu}. This quantity is needed for computing T_{\mu \nu} and T^{\mu \nu}
   and the HARM con2prim lowlevel functions.
*/

typedef struct metric_quantities {
  double adm_gxx, adm_gxy, adm_gxz;
  double adm_gyy, adm_gyz, adm_gzz;
  double adm_gupxx, adm_gupxy, adm_gupxz;
  double adm_gupyy, adm_gupyz, adm_gupzz;
  double betax, betay, betaz;
  double lapse, lapseinv;
  double psi2, psi4, psi6;
  double psi4inv, lapseinv2;
  double g4dn[4][4],g4up[4][4];
} metric_quantities;


//TODO: comment/add to this
typedef struct stress_energy {
  double Ttt, Ttx, Tty, Ttz;
  double Txx, Txy, Txz;
  double Tyy, Tyz, Tzz;
} stress_energy;

void initialize_parameters(GRMHD_parameters *restrict params,
             const int main, const int backup[3], const int evolve_entropy,
             const int evolve_temp, const int calc_prim_guess,
             const double psi6threshold, const int update_Tmunu);

void initialize_metric(metric_quantities *restrict metric, const double lapse,
             const double gxx, const double gxy, const double gxz,
             const double gyy, const double gyz, const double gzz,
             const double betax, const double betay, const double betaz);

//--------------------------------------------------

#endif // GRMHD_HEADER_H
