#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "scattering.h"
#include <complex>

using namespace std;

double get_scatter_chi(MNStruct *par) {
  
  double chi = 0.0;
  int ii, ichan=0,nchans=0;
  double A, b, sigma, t0_inf;
  double gamma = par->gamma;
  double tau1 = par->tau1;
  double DM = par->DM;
  tau1 =  tau1 * 360./par->period;
  double period = par->period;
  int *chan_idx = par->chan_idx;
  double *freq = par->freq;
  double *rmsI = par->rmsI;
  int *nbin = par->nbin;
  int npts=0;
  double *cfreq = par->cfreq;
  double tau;
  double t1, t3, t5, t9, t10, t19, delay, T0, result;
  double win_lo, win_hi;
  ofstream output_file;

  if (par->do_plot) {
    output_file.open("profiles.ascii");
  }
  
  for (unsigned int jj = 0; jj < par->nfiles; jj++) { // Loop over files
    sigma = par->sigma[jj];
    sigma *= 360./period;
    t0_inf = par->t0_inf[jj];
    nchans += chan_idx[jj];
    for (unsigned int i = 0; i < chan_idx[jj]; i++) {
      double *I = par->I[ichan];
      delay = DM * 4.14879e3 * (1./(par->freq[ichan]*par->freq[ichan]) - 1./(cfreq[jj]*cfreq[jj]) ); // delay in s
      T0 = t0_inf + delay/period * 360; // T0 in deg

      tau = tau1 * pow((par->freq[ichan]*1e-3), gamma);

      //cout << jj << " " << i << " "<< delay/period*360. << endl;

      while (T0 < 100) T0+=360.;
      while (T0 > 620.) T0-=360;
      
      A = pow(10, par->A[ichan]);
      b = par->b[ichan];

      win_lo = T0 - par->win / 3.; win_hi = T0 + par->win / 2.;
      if (par->freq[ichan] < 2000.) {
	win_lo = T0 - 160/3; win_hi = T0 + 80;
      }
      while (win_lo > 540) win_lo -= 360;  while (win_hi > 540) win_hi -= 360;
      int lobin = (int)(win_lo * nbin[i]/360.);
      int hibin = (int)(win_hi * nbin[i]/360.);
      //printf("lobin = %d hibin=%d\n", lobin,hibin);
      
      if (lobin>hibin && lobin - nbin[i] > 0) lobin-= nbin[i]; 
      
      for (unsigned int j = lobin; j < hibin; j++) {
	ii = j;
	while(ii>=nbin[i]) ii -= nbin[i];
	double dt = (j*360./(double)nbin[i] - T0) ;
	while (dt < 100) dt+=360.;
	while (dt > 180.) dt-=360.;
	
	t1 = tau * dt;
	t3 = sigma * sigma;
	t5 = tau * tau;
	t9 = exp(-(0.4e1 * t1 - t3) / t5 / 0.4e1);
	t10 = 1.77245385090552;
	t19 = erf((0.2e1 * t1 - t3) / sigma / tau / 0.2e1);
	result = (A*t9 * t10 * sigma * (t19 + 0.1e1) / 0.2e1+b);
	chi += pow((I[ii] - result)* par->scale[ichan] / par->rmsI[ichan], 2); 
	npts++;
	//cout << jj <<  " " << i << " "<< ii << " " << I[ii] << " " << chi <<endl; 

	if (par->do_plot) output_file << ichan << " "<< j/(float)nbin[i] << " " << I[ii] << "  " << result << " " << chi << endl;
	
	//cout << ichan << " " << ii << " " << I[ii] << " " << result << endl;
      }
      ichan++;
    }
  }
  if (par->do_plot) {
    cout << "Chi**2 = "<< chi << " Npts = " << npts << endl;
    output_file.close();
  }
  return chi;
}
