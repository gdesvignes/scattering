#include<vector>

using namespace std;

typedef struct {

  int nchan;    /* Number of channels */
  int nbin;
  
  double period;
  double DM;
  double tau1;
  double max_phase;
  int nfiles;  
  int *chan_idx;
  double **phase;
  double **I;
  
  double *rmsI;
  double *freq;
  double *scale;
  double *efac;
  double *cfreq;
  double *sigma;
  double *t0_inf;
  double *A;
  double *b;
  
  double *r_sigma;
  double *r_tau;
  double *r_DM;
  double *r_t0_inf;
  double *r_gamma;
  double *r_amp;
  double *r_b;
  double *lo, *hi;


} MNStruct;

MNStruct* init_struct(vector<int> chan_idx, double period, double DM, int nbin, vector< vector<double> > phase, vector<double> freq, vector< vector<double> >  I, vector<double> rmsI, vector<double> scale, vector<double> cfreq );
void ScatterLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void ScatterLike2(double *Cube, int &ndim, int &npars, double &lnew, void *context);
double get_scatter_chi(MNStruct *par);
