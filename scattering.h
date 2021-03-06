
#ifndef SCATTERING_H_INCLUDED
#define SCATTERING_H_INCLUDED

#include<vector>

using namespace std;

typedef struct {
  int sampler;
  int ndims;

  int nchan;    /* Number of channels */
  int* nbin;
  int do_plot;
  
  double period;
  double DM;
  double tau1;
  double gamma;
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
  double win;
  int scatter_index_fixed;
  double *r_sigma;
  double *r_tau;
  double *r_DM;
  double *r_t0_inf;
  double *r_gamma;
  double *r_amp;
  double *r_b;
  double *lo, *hi;


} MNStruct;

MNStruct* init_struct(vector<int> chan_idx, double period, double DM, vector<int> nbin, vector< vector<double> > phase, vector<double> freq, vector< vector<double> >  I, vector<double> rmsI, vector<double> scale, vector<double> cfreq );
void ScatterLike_MN(double *Cube, int &ndim, int &npars, double &lnew, void *context);
//void ScatterLike2(double *Cube, int &ndim, int &npars, double &lnew, void *context);
double get_scatter_chi(MNStruct *par);
void read_stats(char *filename, int npar, MNStruct *p);

#endif
