#include "unit_tests.h"

void read_metric_binary(
      metric_quantities *restrict metric,
      FILE *restrict infile) {

  int key;
  key  = fread(&metric->lapse,   sizeof(double), 1, infile);
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

  // Since each read only reads a single double, the key should just be a sum of every read
  // that happens.
  if( key != 10)
    grhayl_error("An error has occured with reading in initial data."
                 "Please check that comparison data"
                 "is up-to-date with current test version.\n");
}

void write_metric_binary(
      const metric_quantities *restrict metric,
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
}
