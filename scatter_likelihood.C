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

	double sigma = Cube[0] * (par->r_sigma[1] - par->r_sigma[0]) + par->r_sigma[0];
	double tau1 = Cube[1] * (par->r_tau[1] - par->r_tau[0]) + par->r_tau[0];
	double DM = Cube[2] * (par->r_DM[1] - par->r_DM[0]) + par->r_DM[0];
	double t0_inf = Cube[3] * (par->r_t0_inf[1] - par->r_t0_inf[0]) + par->r_t0_inf[0];
	double gamma = -4.0;

	tau1 *= 360./par->period;
	sigma *= 360./par->period;

	int ii;
	int nchan = par->nchan;
	double period = par->period;
	double *freq = par->freq;
	double *rmsI = par->rmsI;
	int nbin = par->nbin;


	//cout << sigma << " " << tau1 << " " << DM << " " << t0_inf << " "<< nchan << " "<<period<<endl;

	double tau;
	double t1, t3, t5, t9, t10, t19, delay, T0, result;
	double win_lo, win_hi;

	for (unsigned int i = 0; i < nchan; i++) {
	    double *I = par->I[i];

	    delay = DM * 4.14879e3 * (1./(par->freq[i]*par->freq[i])); // delay in s
	    tau = tau1 * pow((par->freq[i]*1e-3), gamma);

	    // cerr << sigma << " " << tau << " " << DM << " " << t0_inf << " "<< nchan << " "<<period<<endl;
	
	    T0 = t0_inf + delay/period * 360.; // T0 in deg
	    //cerr << T0 << endl;
	    T0 -= DM * 4.14879e3 * (1./(par->freq[nchan-1]*par->freq[nchan-1])) / period * 360.;
	    //cerr << T0 << endl;

	    while (T0 < 100) T0+=360.;
	    while (T0 > 620.) T0-=360;
	    //cerr << delay << " "<< tau << " " << T0 << endl;

	    double A = Cube[4+i*2] * (par->r_amp[1] - par->r_amp[0]) + par->r_amp[0];
	    double b = Cube[5+i*2] * (par->r_b[1] - par->r_b[0]) + par->r_b[0];


	    //delay = DM * 4.14879e3 * (1./(par->freq[i]*par->freq[i])) / period * 360.; // Delay in deg
	    double win = 80;
	    win_lo = T0 - win / 3.; win_hi = T0 + win / 2.;
	    //cerr << "Win: "<< win_lo << " "<< win_hi << endl;
	    while (win_lo > 540) win_lo -= 360;  while (win_hi > 540) win_hi -= 360;
	    
	    //cerr << "Nbins: "<< (int)win_lo * nbin/360. << " "<< (int)win_hi* nbin/360 << endl;

	    //cerr << (int)(win_lo * nbin/360.) << " "<< (int)(win_hi * nbin/360.) << endl;
	    int lobin = (int)(win_lo * nbin/360.);
	    int hibin = (int)(win_hi * nbin/360.);

	    

	    for (unsigned int j = lobin; j < hibin; j++) {
	    //for (unsigned int j = 0; j < 2* nbin; j++) {
      	      ii = j;
	      while(ii>=nbin) ii -= nbin;
	      //cerr << "j: "<< j << endl;
	      //cerr << par->rmsI[i] << endl;
	      //cout << x[0]  << " " << T0 << endl;
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

		chi += pow((I[ii] - result)* par->scale[i] / par->rmsI[i], 2); 
		//cerr << t1 << " " << t3 << " " << t5 << " " << t9 << endl;
		//cerr <<i << " "<< j << " " << result << " "<< I[ii]<<  " Freq: "<< par->freq[i] << " Chi: "<< chi << endl;
		//exit(0);
	    }
	}

	
	Cube[0] = Cube[0] * (par->r_sigma[1] - par->r_sigma[0]) + par->r_sigma[0];
	Cube[1] = Cube[1] * (par->r_tau[1] - par->r_tau[0]) + par->r_tau[0];
	Cube[2] = Cube[2] * (par->r_DM[1] - par->r_DM[0]) + par->r_DM[0];
	Cube[3] = Cube[3] * (par->r_t0_inf[1] - par->r_t0_inf[0]) + par->r_t0_inf[0];
	for (unsigned int i = 0; i < par->nchan; i++) {
	    Cube[4+i*2] = Cube[4+i*2] * (par->r_amp[1] - par->r_amp[0]) + par->r_amp[0];
	    Cube[5+i*2] = Cube[5+i*2] * (par->r_b[1] - par->r_b[0]) + par->r_b[0];
        }

	if (isnan(chi) || isinf(chi))
	  lnew=-pow(10.0,200);
	else
	  lnew = -chi/2;

	//cerr << lnew << endl;
	//exit(0);
}
