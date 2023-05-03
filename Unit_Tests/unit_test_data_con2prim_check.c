#include "unit_tests.h"

int nchecks=0;

FILE *
open_file(const char *filename, const char *type) {
  FILE *fp = fopen(filename, type);
  if( !fp ) grhayl_error("Failed to open file %s\n", filename);
  return fp;
}

#define check(s1, s2, q, tol)                          \
  {                                                    \
    if( (relative_error(s1.q, s2.q) > tol) &&          \
        (fabs(s1.q-s2.q) > tol) )                      \
      grhayl_info("%s mismatch (rel. err.: %e)\n",     \
                   #q, relative_error(s1.q, s2.q));    \
  }

int
initial_metric_data_check(
    const char *dirname,
    const char *routine,
    const double tol ) {

  char filename1[512], filename2[512];
  sprintf(filename1, "%s/%s_initial_data.bin", dirname, routine);
  sprintf(filename2, "%s_initial_data.bin"   , routine);

  FILE *fp1 = open_file(filename1, "rb");
  FILE *fp2 = open_file(filename2, "rb");

  int  __attribute__((unused)) err, npoints1, npoints2;
  err = fread(&npoints1, sizeof(int), 1, fp1);
  err = fread(&npoints2, sizeof(int), 1, fp2);
  if( npoints1 != npoints2 )
    grhayl_error("npoints mismatch: %d %d\n", npoints1, npoints2);

  grhayl_info("Number of points in files: %d\n", npoints1);

  for(int i=0;i<npoints1;i++) {
    metric_quantities m1, m2;
    read_metric_struct_binary(&m1, fp1);
    read_metric_struct_binary(&m2, fp2);
    check(m1, m2, lapse, tol);
    check(m1, m2, betax, tol);
    check(m1, m2, betay, tol);
    check(m1, m2, betaz, tol);
    check(m1, m2, adm_gxx, tol);
    check(m1, m2, adm_gxy, tol);
    check(m1, m2, adm_gxz, tol);
    check(m1, m2, adm_gyy, tol);
    check(m1, m2, adm_gyz, tol);
    check(m1, m2, adm_gzz, tol);
  }

  fclose(fp1);
  fclose(fp2);

  grhayl_info("Files %s and %s agree\n", filename1, filename2);

  nchecks++;
  return npoints1;
}

void
inequality_fixes_data_check(
    const int npoints,
    const char *dirname,
    const char *filename,
    const double tol ) {

  char filename1[512], filename2[512];
  sprintf(filename1, "%s/%s", dirname, filename);
  sprintf(filename2, "%s"   , filename);

  FILE *fp1 = open_file(filename1, "rb");
  FILE *fp2 = open_file(filename2, "rb");

  for(int i=0;i<npoints;i++) {
    primitive_quantities p1, p2;
    read_primitive_struct_binary(0, 0, &p1, fp1);
    read_primitive_struct_binary(0, 0, &p2, fp2);
    check(p1, p2, rho  , tol);
    check(p1, p2, press, tol);
    check(p1, p2, vx   , tol);
    check(p1, p2, vy   , tol);
    check(p1, p2, vz   , tol);
    check(p1, p2, eps  , tol);
    check(p1, p2, Bx   , tol);
    check(p1, p2, By   , tol);
    check(p1, p2, Bz   , tol);
  }

  fclose(fp1);
  fclose(fp2);
  grhayl_info("Files %s and %s agree\n", filename1, filename2);
  nchecks++;
}

void
con2prim_data_check(
    const int npoints,
    const bool perturbed,
    const char *dirname,
    const char *filename,
    const double tol ) {

  char filename1[512], filename2[512];
  sprintf(filename1, "%s/%s", dirname, filename);
  sprintf(filename2, "%s"   , filename);

  FILE *fp1 = open_file(filename1, "rb");
  FILE *fp2 = open_file(filename2, "rb");

  for(int i=0;i<npoints;i++) {

    if( !perturbed ) {
      primitive_quantities p1, p2;
      read_primitive_struct_binary(0, 0, &p1, fp1);
      read_primitive_struct_binary(0, 0, &p2, fp2);
      check(p1, p2, rho  , tol);
      check(p1, p2, press, 1e-12);
      check(p1, p2, vx   , tol);
      check(p1, p2, vy   , tol);
      check(p1, p2, vz   , tol);
      check(p1, p2, eps  , 1e-11);
      check(p1, p2, Bx   , tol);
      check(p1, p2, By   , tol);
      check(p1, p2, Bz   , tol);

      conservative_quantities c1, c2;
      read_conservative_struct_binary(0, &c1, fp1);
      read_conservative_struct_binary(0, &c2, fp2);
      check(c1, c2, rho, tol);
      check(c1, c2, tau, tol);
      check(c1, c2, S_x, tol);
      check(c1, c2, S_y, tol);
      check(c1, c2, S_z, tol);
    }

    primitive_quantities p1, p2;
    read_primitive_struct_binary(0, 0, &p1, fp1);
    read_primitive_struct_binary(0, 0, &p2, fp2);
    check(p1, p2, rho  , tol);
    check(p1, p2, press, 1e-12);
    check(p1, p2, vx   , tol);
    check(p1, p2, vy   , tol);
    check(p1, p2, vz   , tol);
    check(p1, p2, eps  , 1e-11);
    check(p1, p2, Bx   , tol);
    check(p1, p2, By   , tol);
    check(p1, p2, Bz   , tol);
  }

  fclose(fp1);
  fclose(fp2);
  grhayl_info("Files %s and %s agree\n", filename1, filename2);
  nchecks++;
}

void
enforce_limits_data_check(
    const int npoints,
    const bool perturbed,
    const char *dirname,
    const char *filename,
    const double tol ) {

  char filename1[512], filename2[512];
  sprintf(filename1, "%s/%s", dirname, filename);
  sprintf(filename2, "%s"   , filename);

  FILE *fp1 = open_file(filename1, "rb");
  FILE *fp2 = open_file(filename2, "rb");

  for(int i=0;i<npoints;i++) {

    if( !perturbed ) {
      primitive_quantities p1, p2;
      read_primitive_struct_binary(0, 0, &p1, fp1);
      read_primitive_struct_binary(0, 0, &p2, fp2);
      check(p1, p2, rho  , tol);
      check(p1, p2, press, 1e-12);
      check(p1, p2, vx   , tol);
      check(p1, p2, vy   , tol);
      check(p1, p2, vz   , tol);
      check(p1, p2, eps  , 1e-11);
      check(p1, p2, Bx   , tol);
      check(p1, p2, By   , tol);
      check(p1, p2, Bz   , tol);
    }

    primitive_quantities p1, p2;
    read_primitive_struct_binary(0, 0, &p1, fp1);
    read_primitive_struct_binary(0, 0, &p2, fp2);
    check(p1, p2, rho  , tol);
    check(p1, p2, press, 1e-12);
    check(p1, p2, vx   , tol);
    check(p1, p2, vy   , tol);
    check(p1, p2, vz   , tol);
    check(p1, p2, eps  , 1e-11);
    check(p1, p2, Bx   , tol);
    check(p1, p2, By   , tol);
    check(p1, p2, Bz   , tol);

    int err;
    double u01, u02;
    err  = fread(&u01, sizeof(double), 1, fp1);
    err += fread(&u02, sizeof(double), 1, fp2);
    if( err != 2 )
      grhayl_error("Failed to read from enforce limits file\n");
    if( relative_error(u01, u02) > tol )
      grhayl_error("u^0 mismatch (rel. err.: %e)\n", relative_error(u01, u02));
  }

  fclose(fp1);
  fclose(fp2);
  grhayl_info("Files %s and %s agree\n", filename1, filename2);
  nchecks++;
}

void
conservs_and_Tmunu_data_check(
    const int npoints,
    const bool perturbed,
    const char *dirname,
    const char *filename,
    const double tol ) {

  char filename1[512], filename2[512];
  sprintf(filename1, "%s/%s", dirname, filename);
  sprintf(filename2, "%s"   , filename);

  FILE *fp1 = open_file(filename1, "rb");
  FILE *fp2 = open_file(filename2, "rb");

  for(int i=0;i<npoints;i++) {
    int __attribute__((unused)) err;
    primitive_quantities p1, p2;
    conservative_quantities c1, c2;
    stress_energy T1, T2;
    double u01, u02;

    if( !perturbed ) {
      read_primitive_struct_binary(0, 0, &p1, fp1);
      read_primitive_struct_binary(0, 0, &p2, fp2);
      check(p1, p2, rho  , tol);
      check(p1, p2, press, 1e-12);
      check(p1, p2, vx   , tol);
      check(p1, p2, vy   , tol);
      check(p1, p2, vz   , tol);
      check(p1, p2, eps  , 1e-11);
      check(p1, p2, Bx   , tol);
      check(p1, p2, By   , tol);
      check(p1, p2, Bz   , tol);

      err = fread(&u01, sizeof(double), 1, fp1);
      err = fread(&u02, sizeof(double), 1, fp2);
      if( relative_error(u01, u02) > tol )
        grhayl_error("u^0 mismatch (rel. err.: %e)\n", relative_error(u01, u02));
    }

    read_conservative_struct_binary(0, &c1, fp1);
    read_conservative_struct_binary(0, &c2, fp2);
    check(c1, c2, rho, tol);
    check(c1, c2, tau, tol);
    check(c1, c2, S_x, tol);
    check(c1, c2, S_y, tol);
    check(c1, c2, S_z, tol);

    read_stress_energy_struct_binary(&T1, fp1);
    read_stress_energy_struct_binary(&T2, fp2);
    check(T1, T2, Ttt, tol);
    check(T1, T2, Ttx, tol);
    check(T1, T2, Tty, tol);
    check(T1, T2, Ttz, tol);
    check(T1, T2, Txx, tol);
    check(T1, T2, Txy, tol);
    check(T1, T2, Txz, tol);
    check(T1, T2, Tyy, tol);
    check(T1, T2, Tyz, tol);
    check(T1, T2, Tzz, tol);
  }

  fclose(fp1);
  fclose(fp2);
  grhayl_info("Files %s and %s agree\n", filename1, filename2);
  nchecks++;
}

int main(int argc, char **argv) {

  const double tol  = 1e-14;
  const int npoints = initial_metric_data_check("../GRHayL_TestData/con2prim",
                                                "Noble2D", tol);
  inequality_fixes_data_check(npoints,
                              "../GRHayL_TestData/con2prim",
                              "apply_inequality_fixes.bin",
                              tol);
  inequality_fixes_data_check(npoints,
                              "../GRHayL_TestData/con2prim",
                              "apply_inequality_fixes_pert.bin",
                              tol);

  con2prim_data_check(npoints, false,
                      "../GRHayL_TestData/con2prim",
                      "Noble2D_Hybrid_Multi_Method.bin",
                      tol);

  con2prim_data_check(npoints, true,
                      "../GRHayL_TestData/con2prim",
                      "Noble2D_Hybrid_Multi_Method_pert.bin",
                      tol);

  con2prim_data_check(npoints, false,
                      "../GRHayL_TestData/con2prim",
                      "font_fix.bin",
                      tol);

  con2prim_data_check(npoints, true,
                      "../GRHayL_TestData/con2prim",
                      "font_fix_pert.bin",
                      tol);

  enforce_limits_data_check(npoints, false,
                            "../GRHayL_TestData/con2prim",
                            "enforce_primitive_limits_and_compute_u0.bin",
                            tol);

  enforce_limits_data_check(npoints, true,
                            "../GRHayL_TestData/con2prim",
                            "enforce_primitive_limits_and_compute_u0_pert.bin",
                            tol);

  conservs_and_Tmunu_data_check(npoints, true,
                                "../GRHayL_TestData/con2prim",
                                "compute_conservs_and_Tmunu_pert.bin",
                                tol);

  conservs_and_Tmunu_data_check(npoints, false,
                                "../GRHayL_TestData/con2prim",
                                "compute_conservs_and_Tmunu.bin",
                                tol);

  grhayl_info("Number of successful comparisons: %d\n", nchecks);

  return 0;
}
