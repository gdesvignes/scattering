
typedef struct param {
  // Sampler
  int sampler; // 0: MultiNest, 1: Polychord

  // Multinest
  int IS;
  int nlive;
  int ceff;
  double efr;
  char *basename;
  
  // config
  int have_efac;  
  double threshold;
  int margin_phi0;  
  
  // Params
  double *sigma;
  double *tau;
  double *DM;
  double *t0_inf;
  double *efac;
  double *amp;
  double *b;
  double win;
  int scatter_index_fixed;
  // Data
  int nchan;
  int bscrunch;

} param;

int readParameters(param *p, char *paramFile);
