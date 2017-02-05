#include "scattering.h"
#include <stdlib.h>
#include <math.h>

#include <iostream>

using namespace std;

MNStruct* init_struct(int nchan, double period, double DM, int nbin, vector< vector<double> > phase, vector<double> freq, vector< vector<double> > I, vector<double> rmsI, vector<double> scale, vector<double> cfreq) {

	// Init struct
	MNStruct* MNS = (MNStruct*)malloc(sizeof(MNStruct));
	
	// Malloc
	//MNS->npts = (int *)malloc(nfiles * sizeof(int));
	//MNS->efac = (double *)malloc(nfiles * sizeof(double));
	//MNS->epoch = (double *)malloc(nfiles * sizeof(double));
	MNS->period = period;
	MNS->DM = DM;
	MNS->nchan = nchan;

	MNS->freq = (double *)malloc(nchan * sizeof(double));
	MNS->rmsI = (double *)malloc(nchan * sizeof(double));
	MNS->scale = (double *)malloc(nchan * sizeof(double));
	MNS->lo = (double *)malloc(nchan * sizeof(double));
	MNS->hi = (double *)malloc(nchan * sizeof(double));
	MNS->cfreq = (double *)malloc(cfreq.size() * sizeof(double));

	MNS->phase = (double **)malloc(nchan * sizeof(double *));
	MNS->I = (double **)malloc(nchan * sizeof(double *));

	for (int i=0; i<nchan; i++) {
		//cout << i << " "<<  nfiles << " "<<  MNS->npts[i] <<endl;
		MNS->phase[i] = (double *)malloc(nbin * sizeof(double));
		MNS->I[i] = (double *)malloc(nbin * sizeof(double));
		MNS->rmsI[i] = rmsI[i];
		MNS->scale[i] = scale[i];

	}
		
	cout << "Malloc finished" << endl;

	for (int i=0; i<cfreq.size() ; i++) MNS->cfreq[i] = cfreq[i];
	
	// Assign data
	for (int i=0; i<nchan; i++) {
	  MNS->freq[i]=freq[i];
	    for (int j=0; j<nbin; j++) {
		MNS->phase[i][j] = phase[i][j]; // TBD
		MNS->I[i][j] = I[i][j];
		//cout << "Data: "<< MNS->I[i][j] << endl;
	    }
	}

	cout << "Params finished" << endl;

	return MNS;
}
