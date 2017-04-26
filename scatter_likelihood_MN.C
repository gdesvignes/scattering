#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "scattering.h"
#include <complex>

using namespace std;

void ScatterLike_MN(double *Cube, int &ndim, int &npars, double &lnew, void *context) {

	MNStruct *par = ((MNStruct *)context);
	int ipar = 0;
        double chi = 0.0;
	double sigma, t0_inf;
	int ichan=0,nchans=0;
	int *chan_idx = par->chan_idx;

	Cube[ipar] = Cube[ipar] * (par->r_tau[1] - par->r_tau[0]) + par->r_tau[0];
	par->tau1 = Cube[ipar];
	ipar++;

	Cube[ipar] = Cube[ipar] * (par->r_DM[1] - par->r_DM[0]) + par->r_DM[0];
	par->DM = Cube[ipar];
	ipar++;

	if(!par->scatter_index_fixed) {
	  Cube[ipar] = Cube[ipar] * (1.6) - 4.8;
	  par->gamma = Cube[ipar];
	  ipar++;
	} else {
	  par->gamma = -4.;
	}

	for (unsigned int jj = 0; jj < par->nfiles; jj++) { // Loop over files
	  Cube[ipar] = Cube[ipar] * (par->r_sigma[1] - par->r_sigma[0]) + par->r_sigma[0];
	  par->sigma[jj] = Cube[ipar + nchans*2 + jj*2];
	  ipar++;

          Cube[ipar] = Cube[ipar] * (par->r_t0_inf[1] - par->r_t0_inf[0]) + par->r_t0_inf[0];
	  par->t0_inf[jj] = Cube[ipar];
	  ipar++;
	  nchans += chan_idx[jj];

	  for (unsigned int i = 0; i < chan_idx[jj]; i++) {
	    Cube[ipar] = Cube[ipar] * (par->r_amp[1] - par->r_amp[0]) + par->r_amp[0];
	    par->A[ichan] = Cube[ipar];
	    ipar++;

       	    Cube[ipar] = Cube[ipar] * (par->r_b[1] - par->r_b[0]) + par->r_b[0];
	    par->b[ichan] = Cube[ipar];
	    ipar++;

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
