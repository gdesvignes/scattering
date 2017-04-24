#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "scattering.h"
#include <complex>

#include "Parameters.h"
#include "scatter_likelihood_PC.h"
using namespace std;

// Prior function                                                                                                                      
//                                                                                                                                     
// Either write your prior code directly into this function, or call an                                                                
// external library from it. This should transform a coordinate in the unit hypercube                                                  
// stored in cube (of size nDims) to a coordinate in the physical system stored in theta                                               
//                                                                                                                                     
// This function is called from likelihoods/fortran_cpp_wrapper.f90                                                                    
// If you would like to adjust the signature of this call, then you should adjust it there,                                            
// as well as in likelihoods/my_cpp_likelihood.hpp                                                                                     
//                                                                                                                                     
void prior (double cube[], double theta[], int nDims,  void *context)
{
  //MNStruct *par = ((MNStruct *)context);
  double chi = 0.0;
  double sigma, t0_inf;
  int ichan=0,nchans=0;
  //int *chan_idx = sp->chan_idx;
  //printf("Nfiles = %d\n", sp->nfiles); fflush(stdout);

  int *chan_idx = sp->chan_idx; 
                                                                 
  theta[0] = cube[0] * (sp->r_tau[1] - sp->r_tau[0]) + sp->r_tau[0];
  theta[1] = cube[1] * (sp->r_DM[1] - sp->r_DM[0]) + sp->r_DM[0];

  for (unsigned int jj = 0; jj < sp->nfiles; jj++) { // Loop over files                   
    theta[2 + nchans*2 + jj*2] = cube[2 + nchans*2 + jj*2] * (sp->r_sigma[1] - sp->r_sigma[0]) + sp->r_sigma[0];
    theta[3 + nchans*2 + jj*2] = cube[3 + nchans*2 + jj*2] * (sp->r_t0_inf[1] - sp->r_t0_inf[0]) + sp->r_t0_inf[0];
    nchans += chan_idx[jj];

    for (unsigned int i = 0; i < chan_idx[jj]; i++) {
      theta[4 + ichan*2 + jj*2] = cube[4 + ichan*2 + jj*2] * (sp->r_amp[1] - sp->r_amp[0]) + sp->r_amp[0];
      theta[5  + ichan*2 + jj*2] = cube[5 + ichan*2 + jj*2] * (sp->r_b[1] - sp->r_b[0]) + sp->r_b[0];

      ichan++;
    }
  }
}

double ScatterLike_PC(double theta[], int nDims, double phi[], int nDerived, void *context) {


        double chi = 0.0, lnew;
	double sigma, t0_inf;
	int ichan=0,nchans=0;
	int *chan_idx = sp->chan_idx;

	sp->tau1 = theta[0] ;
	sp->DM = theta[1];

	for (unsigned int jj = 0; jj < sp->nfiles; jj++) { // Loop over files
	  sp->sigma[jj] = theta[2 + nchans*2 + jj*2];
	  sp->t0_inf[jj] = theta[3 + nchans*2 + jj*2];
	  nchans += chan_idx[jj];

	  for (unsigned int i = 0; i < chan_idx[jj]; i++) {
	    sp->A[ichan] = theta[4 + ichan*2 + jj*2];
	    sp->b[ichan] = theta[5 + ichan*2 + jj*2];

	    ichan++;
	  }
	}

	// Call function
	chi = get_scatter_chi(sp);

	if (isnan(chi) || isinf(chi))
	  lnew=-pow(10.0,200);
	else
	  lnew = -chi/2;
	return lnew;
}
