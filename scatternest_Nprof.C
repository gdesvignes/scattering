#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <mpi.h>
#include "multinest.h"
#include "scattering.h"
#include "Parameters.h"

// psrchive stuff
#include "Pulsar/Archive.h"
#include "Pulsar/Integration.h"
#include "Pulsar/Profile.h"
//#include "Pulsar/FaradayRotation.h"
//#include "Pulsar/PolnProfileStats.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;
using namespace Pulsar;



/************************************************* dumper routine ******************************************************/

// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants
//
//
// Arguments:
//
// nSamples 						= total number of samples in posterior distribution
// nlive 						= total number of live points
// nPar 						= total number of parameters (free + derived)
// physLive[1][nlive * (nPar + 1)] 			= 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
// posterior[1][nSamples * (nPar + 2)] 			= posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
// paramConstr[1][4*nPar]:
// paramConstr[0][0] to paramConstr[0][nPar - 1] 	= mean values of the parameters
// paramConstr[0][nPar] to paramConstr[0][2*nPar - 1] 	= standard deviation of the parameters
// paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters
// paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters
// maxLogLike						= maximum loglikelihood value
// logZ							= log evidence value from the default (non-INS) mode
// INSlogZ						= log evidence value from the INS mode
// logZerr						= error on log evidence value
// context						void pointer, any additional information

void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context)
{
	// convert the 2D Fortran arrays to C++ arrays
	
	
	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
	
	int i, j;
	
	double postdist[nSamples][nPar + 2];
	for( i = 0; i < nPar + 2; i++ )
		for( j = 0; j < nSamples; j++ )
			postdist[j][i] = posterior[0][i * nSamples + j];
	
	
	
	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
	
	double pLivePts[nlive][nPar + 1];
	for( i = 0; i < nPar + 1; i++ )
		for( j = 0; j < nlive; j++ )
			pLivePts[j][i] = physLive[0][i * nlive + j];
}

/***********************************************************************************************************************/




/************************************************** Main program *******************************************************/



int main(int argc, char *argv[])
{
	// set the MultiNest sampling parameters
	
	int IS = 0;					// do Nested Importance Sampling?
	int mmodal = 0;					// do mode separation?
	int ceff = 0;					// run in constant efficiency mode?
	int nlive = 2000;				// number of live points
	double efr = 0.1;				// set the required efficiency
	double tol = 0.5;				// tol, defines the stopping criteria
	int ndims = 4;					// dimensionality (no. of free parameters)
	int nPar = 4;					// total no. of parameters including free & derived parameters
	int nClsPar = 2;				// no. of parameters to do mode separation on
	int updInt = 1000;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(int i = 0; i < ndims; i++) pWrap[i] = 0;
	//pWrap[0] = 1; pWrap[1] = 1; pWrap[2] = 1; pWrap[3] = 1;
	
	char confname[128];
	char filename[128];
	char root[100] = "chains/eggboxCC-";		// root for output files
	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	int fb = 1;					// need feedback on standard output?
	int resume = 1;					// resume from a previous job?
	int outfile = 1;				// write output files?
	int initMPI = 0;				// initialize MPI routines?, relevant only if compiling with MPI
							// set it to F if you want your main program to handle MPI initialization
	double logZero = -1E90;				// points with loglike < logZero will be ignored by MultiNest
	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
							// has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
	void *context = 0;				// not required by MultiNest, any additional information user wants to pass
	int nbin;
	int ichan = 0;
	int bscrunch = 1;
	int nchan = 8;
	double DM, period;
	vector< vector<double> > phase, I;

	int rank, size;
        MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

	//if (rank==0) {

	int nfiles = argc - 1;
	cerr << "Nfiles = "<< nfiles<< endl;	
	param p;
        cerr << "Reading parameters" << endl;
        int rv = readParameters(&p, "config.txt");

        if (rv == EXIT_SUCCESS) {
          IS = p.IS;
          nlive = p.nlive;
          ceff = p.ceff;
          efr = p.efr;
          nchan = p.nchan;
          bscrunch = p.bscrunch;
	}

	vector <int> chan_idx;
	vector <double> RMS_I;
	vector <double> scale;
	vector <double> cfreq, freq;

	chan_idx.resize(nfiles);
	phase.resize(nchan * nfiles);
	I.resize(nchan * nfiles);

	for (int jj=0; jj<nfiles; jj++) {

	  nPar = 2 + (2 + nchan * 2) * nfiles; // DM and scattering are assumed to be constant here; width, phase are profile dependent 
	  ndims = nPar;
	  
	  strcpy(filename, argv[jj+1]);
	  cerr << "Filename: "<< filename<< " " << nchan << endl;


	  Reference::To< Archive > archive = Archive::load( filename );
	  if( !archive ) return -1;
	  cerr << "Archive loaded "<< endl;
	  
	  archive->tscrunch();
	  archive->fscrunch_to_nchan(nchan);
	  if (bscrunch > 1) archive->bscrunch(bscrunch);
	  //int nchan = archive->get_nchan();
	  nbin = archive->get_nbin();
	  DM = archive->get_dispersion_measure();
	  cerr << nbin<< " "<< DM << endl;
	  if( archive->get_state() != Signal::Stokes) archive->convert_state(Signal::Stokes);
	  archive->remove_baseline();
	  cfreq.push_back(archive->get_centre_frequency());
	  cerr << "Finished reading file"<< endl;
	
	  // Get Data
	  Pulsar::Integration* integration = archive->get_Integration(0);
	  period = integration->get_folding_period();

	  vector< vector< double > > std_variance;
	  integration->baseline_stats (0, &std_variance);

	  if (jj == 0 && rv == EXIT_SUCCESS) { // Create directory on first profile only
	    strcpy(root, p.basename);
	    sprintf(root, "%d_%d", (int)integration->get_epoch().intday(), (int)archive->get_centre_frequency());
	    mkdir(root, 0700);
	    sprintf(root, "%d_%d/chains-", (int)integration->get_epoch().intday(), (int)archive->get_centre_frequency());
	    // Copy config file                                                                                           
	    sprintf(confname,"%s.config", root);
	    ifstream  src("config.txt", ios::binary);
	    ofstream  dst(confname,   ios::binary);
	    dst << src.rdbuf();
	  }

	  chan_idx[jj] = nchan;
	  for (int ii=0; ii<archive->get_nchan(); ii++) {
	    //Skip profile with zero-weight
	    if (integration->get_Profile(0,ii)->get_weight() == 0.0) {
	      chan_idx[jj]--;
	      nPar -=2;
	      ndims -=2;
	      continue;
	    }
	    freq.push_back(integration->get_Profile(0, ii)->get_centre_frequency());
	    for (int ibin=0; ibin< archive->get_nbin(); ibin++) {
	      phase[ichan].push_back((ibin+.5)*(period/(double) archive->get_nbin()));
	      I[ichan].push_back(integration->get_Profile(0,ii)->get_amps()[ibin%nbin] / integration->get_Profile(0,ii)->max());
	      //cout << ibin << " " << phase[ichan][ibin] << " " << I[ichan][ibin] << endl;
	      //cout << ibin << " " <<I.back() << " " << 0.5 * atan(U.back(),  Q.back()) <<endl;
	      //if (L.back() > max_L) max_L = L.back();
	    }
	    RMS_I.push_back(sqrt(std_variance[0][ii]));
	    scale.push_back(integration->get_Profile(0,ii)->max());
	    ichan ++;
	  }
	}


	// Init struct
	context = init_struct(chan_idx, period, DM, nbin, phase, freq, I, RMS_I, scale, cfreq);

	MNStruct *par = ((MNStruct *)context);
	
	par->r_sigma = p.sigma;
	par->r_tau = p.tau;
	par->r_DM = p.DM;
	par->r_t0_inf = p.t0_inf;
	par->r_amp = p.amp;
	par->nbin = nbin;
	par->r_b = p.b;

	
	// calling MultiNest
	nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, ScatterLike, dumper, context);

	MPI_Finalize();

}

/***********************************************************************************************************************/
