#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "scattering.h"
#include <complex>

using namespace std;



/******************************************** loglikelihood routine ****************************************************/

// Now an example, sample an egg box likelihood

// Input arguments
// ndim                                                 = dimensionality (total number of free parameters) of the problem
// npars                                                = total number of free plus derived parameters
// context                                              void pointer, any additional information
//
// Input/Output arguments
// Cube[npars]                                          = on entry has the ndim parameters in unit-hypercube
//                                                      on exit, the physical parameters plus copy any derived parameters you want to store with the free parameters
//       
// Output arguments
// lnew                                                 = loglikelihood

void ScatterLike(double *Cube, int &ndim, int &npars, double &lnew, void *context) {

	MNStruct *par = ((MNStruct *)context);

        double chi = 0.0;
// 	cerr << "Here"<< endl;
// 	cerr << par->r_sigma[0] << " "<< par->r_sigma[1]<< endl;
// 	cerr << par->r_tau[0] << " "<< par->r_tau[1]<< endl;
// 	cerr << par->r_DM[0] << " "<< par->r_DM[1]<< endl;
// 	cerr << par->r_t0_inf[0] << " "<< par->r_t0_inf[1]<< endl;
	double sigma, t0_inf;
	Cube[0] = Cube[0] * (par->r_tau[1] - par->r_tau[0]) + par->r_tau[0];
	Cube[1] = Cube[1] * (par->r_DM[1] - par->r_DM[0]) + par->r_DM[0];
	double A,b ;
	double gamma = -4.0;
	double tau1 = Cube[0];
	double DM = Cube[1];
	tau1 =  tau1 * 360./par->period;
	int ii, ichan=0,nchans=0;
	double period = par->period;
	int *chan_idx = par->chan_idx;
	double *freq = par->freq;
	double *rmsI = par->rmsI;
	int nbin = par->nbin;
	double cfreq = par->cfreq[0];
	double tau;
	double t1, t3, t5, t9, t10, t19, delay, T0, result;
	double win_lo, win_hi;

	for (unsigned int jj = 0; jj < par->nfiles; jj++) { // Loop over files
	  Cube[2 + nchans*2 + jj*2] = Cube[2 + nchans*2 + jj*2] * (par->r_sigma[1] - par->r_sigma[0]) + par->r_sigma[0];
	  sigma = Cube[2 + nchans*2 + jj*2];
	  sigma *= 360/period;
          Cube[3 + nchans*2 + jj*2] = Cube[3 + nchans*2 + jj*2] * (par->r_t0_inf[1] - par->r_t0_inf[0]) + par->r_t0_inf[0];
	  t0_inf = Cube[3 + nchans*2 + jj*2];
	  nchans += chan_idx[jj];

	  for (unsigned int i = 0; i < chan_idx[jj]; i++) {
	    double *I = par->I[ichan];
	    delay = DM * 4.14879e3 * (1./(par->freq[ichan]*par->freq[ichan])); // delay in s
	    tau = tau1 * pow((par->freq[ichan]*1e-3), gamma);
	    
	    T0 = t0_inf + delay/period * 360.; // T0 in deg
	    T0 -= DM * 4.14879e3 * (1./(cfreq*cfreq)) / period * 360.; // TODO

	    while (T0 < 100) T0+=360.;
	    while (T0 > 620.) T0-=360;

	    Cube[4 + ichan*2 + jj*2] = Cube[4 + ichan*2 + jj*2] * (par->r_amp[1] - par->r_amp[0]) + par->r_amp[0];
	    A = Cube[4 + ichan*2 + jj*2];
       	    Cube[5  + ichan*2 + jj*2] = Cube[5 + ichan*2 + jj*2] * (par->r_b[1] - par->r_b[0]) + par->r_b[0];
	    b = Cube[5 + ichan*2 + jj*2];

	    double win = 80; // Set a restricted window for chi**2 calculation 
	    win_lo = T0 - win / 3.; win_hi = T0 + win / 2.;
	    while (win_lo > 540) win_lo -= 360;  while (win_hi > 540) win_hi -= 360;
	    int lobin = (int)(win_lo * nbin/360.);
	    int hibin = (int)(win_hi * nbin/360.);

	    for (unsigned int j = lobin; j < hibin; j++) {
	      ii = j;
	      while(ii>=nbin) ii -= nbin;
	      double dt = (j*360./(double)nbin - T0) ;
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
	      //cout << jj <<  " " << i << " "<< ii << " " << I[ii] << " " << chi <<endl; 
	      //printf(" %lf %lf %lf %lf\n", t1,t3,t5,t9);
	      //cout << I[ii] << "  " << result << " " << chi << endl;
	      //cout << ichan << " " << ii << " " << I[ii] << " " << result << endl;
	    }
	    ichan++;
	  }
	}

	if (isnan(chi) || isinf(chi))
	  lnew=-pow(10.0,200);
	else
	  lnew = -chi/2;

}
