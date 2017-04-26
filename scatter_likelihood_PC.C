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

  double chi = 0.0;
  double sigma, t0_inf;
  int ichan=0,nchans=0,ipar=0;
  int *chan_idx = sp->chan_idx; 
                                                                 
  theta[ipar] = cube[ipar] * (sp->r_tau[1] - sp->r_tau[0]) + sp->r_tau[0];
  ipar++;
  theta[ipar] = cube[ipar] * (sp->r_DM[1] - sp->r_DM[0]) + sp->r_DM[0];
  ipar++;
  
  if(!sp->scatter_index_fixed) {
    theta[ipar] = cube[ipar] * (1.6) - 4.8;
    ipar++;
  }


  for (unsigned int jj = 0; jj < sp->nfiles; jj++) { // Loop over files                   
    theta[ipar] = cube[ipar] * (sp->r_sigma[1] - sp->r_sigma[0]) + sp->r_sigma[0];
    ipar++;
    theta[ipar] = cube[ipar] * (sp->r_t0_inf[1] - sp->r_t0_inf[0]) + sp->r_t0_inf[0];
    ipar++;
    nchans += chan_idx[jj];

    for (unsigned int i = 0; i < chan_idx[jj]; i++) {
      theta[ipar] = cube[ipar] * (sp->r_amp[1] - sp->r_amp[0]) + sp->r_amp[0];
      ipar++;
      theta[ipar] = cube[ipar] * (sp->r_b[1] - sp->r_b[0]) + sp->r_b[0];
      ipar++;

      ichan++;
    }
  }
}

double ScatterLike_PC(double theta[], int nDims, double phi[], int nDerived, void *context) {


        double chi = 0.0, lnew;
	double sigma, t0_inf;
	int ichan=0,nchans=0,ipar=0;
	int *chan_idx = sp->chan_idx;

	sp->tau1 = theta[ipar];
	ipar++;
	sp->DM = theta[ipar];
	ipar++;
	if(!sp->scatter_index_fixed) {
	  sp->gamma = theta[ipar];
	  ipar++;
	} else {sp->gamma = -4;}


	for (unsigned int jj = 0; jj < sp->nfiles; jj++) { // Loop over files
	  sp->sigma[jj] = theta[ipar]; ipar++;
	  sp->t0_inf[jj] = theta[ipar]; ipar++;
	  nchans += chan_idx[jj];

	  for (unsigned int i = 0; i < chan_idx[jj]; i++) {
	    sp->A[ichan] = theta[ipar]; ipar++;
	    sp->b[ichan] = theta[ipar]; ipar++;

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
