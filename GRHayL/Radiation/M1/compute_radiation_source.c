#include "ghl.h"
#include "ghl_m1.h"

/* return dthin */
inline double radM1_set_dthin(double chi)
{
  return 1.5*chi-0.5;
}


/* return dthick */
inline double radM1_set_dthick(double chi)
{
  return 1.5*(1-chi);
}

// Low level kernel computing the Jacobian matrix
void __source_jacobian_low_level(
        double *qpre, double Fup[4], double F2,
        double chi,
        double kapa, double kaps,
        double vup[4], double vdown[4], double v2,
        double W,
        double alpha,
        double cdt,
        double *qstar,
        ghl_radiation_con_source_tensor * J)
{
  const double kapas = kapa+kaps;
  const double alpW  = alpha * W;

  const double dthin = radM1_set_dthin(chi);
  const double dthick = radM1_set_dthick(chi);

  const double vx = vdown[1];
  const double vy = vdown[2];
  const double vz = vdown[3];
  const double W2 = SQ(W);
  const double W3 = W2*W;

  const double vdotF = Fup[1]*vdown[1] + Fup[2]*vdown[2] + Fup[3]*vdown[3];
  const double normF = sqrt(F2);
  const double inormF = (normF > 0 ? 1/normF : 0);
  const double vdothatf = vdotF*inormF;
  const double vdothatf2 = SQ(vdothatf);
  const double hatfx = qpre[1]*inormF; // hatf_i
  const double hatfy = qpre[2]*inormF;
  const double hatfz = qpre[3]*inormF;
  const double hatfupx = Fup[1]*inormF; // hatf^i
  const double hatfupy = Fup[2]*inormF;
  const double hatfupz = Fup[3]*inormF;
  const double e = qpre[0];
  const double eonormF = min(e*inormF, 1.0); // with factor dthin ...

  // drvts of J 
  // Eq (65)
  double JdE = W2 + dthin*vdothatf2*W2 + (dthick*(3 - 2*W2)*(-1 + W2))/(1 + 2*W2);
  // Eq (69)
  double JdFv = 2*W2*(-1 + (dthin*eonormF*vdothatf) + (2*dthick*(-1 + W2))/(1 + 2*W2));
  // Eq (70)
  double JdFf = (-2*dthin*eonormF*vdothatf2*W2);
  // Eq (66)
  double JdFx = JdFv * vup[1] + JdFf * hatfupx;
  double JdFy = JdFv * vup[2] + JdFf * hatfupy;
  double JdFz = JdFv * vup[3] + JdFf * hatfupz;

  // drvts of Hi
  // Eq (71)
  double HdEv = W3*(-1 - dthin*vdothatf2 + (dthick*(-3 + 2*W2))/(1 + 2*W2));
  // Eq (72)
  double HdEf = -(dthin*vdothatf*W);
  // Eq (67)
  double HxdE = HdEv * vx + HdEf * hatfx;
  double HydE = HdEv * vy + HdEf * hatfy;
  double HzdE = HdEv * vz + HdEf * hatfz;
  // Eq (73)
  double HdFdelta = (1 - dthick*v2 - (dthin*eonormF*vdothatf))*W;
  // Eq (74)
  double HdFvv = (2*(1 - dthin*eonormF*vdothatf)*W3) + dthick*W*(2 - 2*W2 + 1/(-1 - 2*W2));
  // Eq (75)
  double HdFff = (2*dthin*eonormF*vdothatf*W);
  // Eq (76)
  double HdFvf = (2*dthin*eonormF*vdothatf2*W3);
  // Eq (77)
  double HdFfv = -(dthin*eonormF*W);
  // Eq (68)
  double HxdFx = HdFdelta + HdFvv * vx * vup[1] + HdFff * hatfx * hatfupx + HdFvf * vx * hatfupx + HdFfv * hatfx * vup[1];
  double HydFx = HdFvv * vy * vup[1] + HdFff * hatfy * hatfupx + HdFvf * vy * hatfupx + HdFfv * hatfy * vup[1];
  double HzdFx = HdFvv * vz * vup[1] + HdFff * hatfz * hatfupx + HdFvf * vz * hatfupx + HdFfv * hatfz * vup[1];

  double HxdFy = HdFvv * vx * vup[2] + HdFff * hatfx * hatfupy + HdFvf * vx * hatfupy + HdFfv * hatfx * vup[2];
  double HydFy = HdFdelta + HdFvv * vy * vup[2] + HdFff * hatfy * hatfupy + HdFvf * vy * hatfupy + HdFfv * hatfy * vup[2];
  double HzdFy = HdFvv * vz * vup[2] + HdFff * hatfz * hatfupy + HdFvf * vz * hatfupy + HdFfv * hatfz * vup[2];

  double HxdFz = HdFvv * vx * vup[3] + HdFff * hatfx * hatfupz + HdFvf * vx * hatfupz + HdFfv * hatfx * vup[3];
  double HydFz = HdFvv * vy * vup[3] + HdFff * hatfy * hatfupz + HdFvf * vy * hatfupz + HdFfv * hatfy * vup[3];
  double HzdFz = HdFdelta + HdFvv * vz * vup[3] + HdFff * hatfz * hatfupz + HdFvf * vz * hatfupz + HdFfv * hatfz * vup[3];

  // Build the Jacobian
  // Eq (61)
  double J00 = - alpW * ( kapas - kaps * JdE);
  // Eq (62)
  double J0x = + alpW * kaps * JdFx + alpW * kapas * vup[1];
  double J0y = + alpW * kaps * JdFy + alpW * kapas * vup[2];
  double J0z = + alpW * kaps * JdFz + alpW * kapas * vup[3];
  // Eq (63)
  double Jx0 = - alpha * ( kapas * HxdE + W * kapa * vx * JdE );
  double Jy0 = - alpha * ( kapas * HydE + W * kapa * vy * JdE );
  double Jz0 = - alpha * ( kapas * HzdE + W * kapa * vz * JdE );
  // Eq (64)
  double Jxx = - alpha * ( kapas * HxdFx + W * kapa * vx * JdFx );
  double Jxy = - alpha * ( kapas * HxdFy + W * kapa * vx * JdFy );
  double Jxz = - alpha * ( kapas * HxdFz + W * kapa * vx * JdFz );

  double Jyy = - alpha * ( kapas * HydFx + W * kapa * vy * JdFx );
  double Jyx = - alpha * ( kapas * HydFy + W * kapa * vy * JdFy );
  double Jyz = - alpha * ( kapas * HydFz + W * kapa * vy * JdFz );

  double Jzx = - alpha * ( kapas * HzdFx + W * kapa * vz * JdFx );
  double Jzy = - alpha * ( kapas * HzdFy + W * kapa * vz * JdFy );
  double Jzz = - alpha * ( kapas * HzdFz + W * kapa * vz * JdFz );

  // Store Jacobian into J
  double A_data[4][4] = { 1 - cdt*J00, - cdt*J0x, - cdt*J0y, - cdt*J0z,
		      - cdt*Jx0, 1 - cdt*Jxx, - cdt*Jxy, - cdt*Jxz,
		      - cdt*Jy0, - cdt*Jyx, 1 - cdt*Jyy, - cdt*Jyz,
		      - cdt*Jz0, - cdt*Jzx, - cdt*Jzy, 1 - cdt*Jzz, };
  for (int a = 0; a < 4; ++a)
  for (int b = 0; b < 4; ++b) {
    // gsl_matrix_set(J, a, b, A_data[a][b]);
    J->DD[a][b] = A_data[a][b];
  }
}