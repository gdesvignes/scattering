
typedef struct param {
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
    // Data
    int nchan;
    int bscrunch;

} param;

int readParameters(param *p, char *paramFile);
