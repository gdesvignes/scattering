#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "scattering.h"
#include <complex>

using namespace std;

void ScatterLike_MN(double *Cube, int &ndim, int &npars, double &lnew, void *context) {

	MNStruct *par = ((MNStruct *)context);

        double chi = 0.0;
	double sigma, t0_inf;
	int ichan=0,nchans=0;
	int *chan_idx = par->chan_idx;

	Cube[0] = Cube[0] * (par->r_tau[1] - par->r_tau[0]) + par->r_tau[0];
	par->tau1 = Cube[0] ;
	Cube[1] = Cube[1] * (par->r_DM[1] - par->r_DM[0]) + par->r_DM[0];
	par->DM = Cube[1];

	for (unsigned int jj = 0; jj < par->nfiles; jj++) { // Loop over files
	  Cube[2 + nchans*2 + jj*2] = Cube[2 + nchans*2 + jj*2] * (par->r_sigma[1] - par->r_sigma[0]) + par->r_sigma[0];
	  par->sigma[jj] = Cube[2 + nchans*2 + jj*2];

          Cube[3 + nchans*2 + jj*2] = Cube[3 + nchans*2 + jj*2] * (par->r_t0_inf[1] - par->r_t0_inf[0]) + par->r_t0_inf[0];
	  par->t0_inf[jj] = Cube[3 + nchans*2 + jj*2];
	  nchans += chan_idx[jj];

	  for (unsigned int i = 0; i < chan_idx[jj]; i++) {
	    Cube[4 + ichan*2 + jj*2] = Cube[4 + ichan*2 + jj*2] * (par->r_amp[1] - par->r_amp[0]) + par->r_amp[0];
	    par->A[ichan] = Cube[4 + ichan*2 + jj*2];
	    par->A[ichan] = pow(10, par->A[ichan]);
       	    Cube[5  + ichan*2 + jj*2] = Cube[5 + ichan*2 + jj*2] * (par->r_b[1] - par->r_b[0]) + par->r_b[0];
	    par->b[ichan] = Cube[5 + ichan*2 + jj*2];

	    ichan++;
	  }
	}

	// Call function
	chi = get_scatter_chi(par);

	if (isnan(chi) || isinf(chi))
	  lnew=-pow(10.0,200);
	else
	  lnew = -chi/2;

}
