#include "unit_tests.h"

void read_c2p_initial_data_binary(
      const int eos_type,
      const bool evolve_entropy,
      metric_quantities *restrict metric,
      primitive_quantities *restrict prims,
      FILE *restrict infile) {

  int key;
  key += fread(&metric->lapse,   sizeof(double), 1, infile);
  key += fread(&metric->adm_gxx, sizeof(double), 1, infile);
  key += fread(&metric->adm_gxy, sizeof(double), 1, infile);
  key += fread(&metric->adm_gxz, sizeof(double), 1, infile);
  key += fread(&metric->adm_gyy, sizeof(double), 1, infile);
  key += fread(&metric->adm_gyz, sizeof(double), 1, infile);
  key += fread(&metric->adm_gzz, sizeof(double), 1, infile);
  key += fread(&metric->betax,   sizeof(double), 1, infile);
  key += fread(&metric->betay,   sizeof(double), 1, infile);
  key += fread(&metric->betaz,   sizeof(double), 1, infile);

  initialize_metric(metric->lapse,
                    metric->adm_gxx, metric->adm_gxy, metric->adm_gxz,
                    metric->adm_gyy, metric->adm_gyz, metric->adm_gzz,
                    metric->betax, metric->betay, metric->betaz,
                    metric);

  key += fread(&prims->rho,   sizeof(double), 1, infile);
  key += fread(&prims->press, sizeof(double), 1, infile);
  key += fread(&prims->vx,    sizeof(double), 1, infile);
  key += fread(&prims->vy,    sizeof(double), 1, infile);
  key += fread(&prims->vz,    sizeof(double), 1, infile);
  key += fread(&prims->Bx,    sizeof(double), 1, infile);
  key += fread(&prims->By,    sizeof(double), 1, infile);
  key += fread(&prims->Bz,    sizeof(double), 1, infile);

  if(evolve_entropy)
    key += fread(&prims->entropy, sizeof(double), 1, infile);

  if(eos_type == 2) { //Tabulated
    key += fread(&prims->Y_e, sizeof(double), 1, infile);
    key += fread(&prims->temperature, sizeof(double), 1, infile);
  }

  // Since each read only reads a single double, the key should just be a sum of every read
  // that happens. Hence, 18 for the metric and primitives, +1 for entropy, and
  // +2 for tabulated (Y_e and temperature).
  const int correct_key = 18 + evolve_entropy + (eos_type == 2)*2;
  if( key != correct_key) {
    printf("An error has occured with reading in initial data."
           "Please check that comparison data"
           "is up-to-date with current test version.\n");
    exit(1);
  }
}

void write_c2p_initial_data_binary(
      const int eos_type,
      const bool evolve_entropy,
      const metric_quantities *restrict metric,
      const primitive_quantities *restrict prims,
      FILE *restrict outfile) {

  fwrite(&metric->lapse,   sizeof(double), 1, outfile);
  fwrite(&metric->adm_gxx, sizeof(double), 1, outfile);
  fwrite(&metric->adm_gxy, sizeof(double), 1, outfile);
  fwrite(&metric->adm_gxz, sizeof(double), 1, outfile);
  fwrite(&metric->adm_gyy, sizeof(double), 1, outfile);
  fwrite(&metric->adm_gyz, sizeof(double), 1, outfile);
  fwrite(&metric->adm_gzz, sizeof(double), 1, outfile);
  fwrite(&metric->betax,   sizeof(double), 1, outfile);
  fwrite(&metric->betay,   sizeof(double), 1, outfile);
  fwrite(&metric->betaz,   sizeof(double), 1, outfile);

  fwrite(&prims->rho,   sizeof(double), 1, outfile);
  fwrite(&prims->press, sizeof(double), 1, outfile);
  fwrite(&prims->vx,    sizeof(double), 1, outfile);
  fwrite(&prims->vy,    sizeof(double), 1, outfile);
  fwrite(&prims->vz,    sizeof(double), 1, outfile);
  fwrite(&prims->Bx,    sizeof(double), 1, outfile);
  fwrite(&prims->By,    sizeof(double), 1, outfile);
  fwrite(&prims->Bz,    sizeof(double), 1, outfile);

  if(evolve_entropy)
    fwrite(&prims->entropy, sizeof(double), 1, outfile);

  if(eos_type == 2) { //Tabulated
    fwrite(&prims->Y_e, sizeof(double), 1, outfile);
    fwrite(&prims->temperature, sizeof(double), 1, outfile);
  }
}
