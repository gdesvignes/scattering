#include "string.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include "scattering.h"

using namespace std;

double get_val(ifstream *s) {
  int ipar;
  double val, val_err;
  string line;
  getline(*s, line);
  istringstream iss(line);
  iss >> ipar; iss >> val;  
  return val;
}

void read_stats(char *root, int npar, MNStruct *p)
{

  int i, ipar,ichan, nmodes;
  string line, a,b,c;
  char filename[256];

  double val, val_err;

  sprintf(filename, "%sstats.dat", root);

  ifstream stats(filename);

  ichan = 0;
  // Ask to output the result
  p->do_plot = 1;

  if (stats.is_open())
    {
      //getline(stats, line);

      for(i=0; i<3; i++) getline(stats, line);
      // Read the number of modes
      //istringstream iss(line);
      //iss >> a; iss>>b; iss >> c; iss>> nmodes;
      //cout << "Found Nmodes: " << nmodes << endl;
      
      //for(i=0; i<7; i++) getline(stats, line);

      for(i=0; i<npar; i++) getline(stats, line);
      
      // Skip 3 lines to go to ML results
      for(i=0; i<3; i++) getline(stats, line);

      p->tau1 = get_val(&stats);
      p->DM = get_val(&stats);
            
      for (unsigned int j = 0; j < p->nfiles; j++) {
	p->sigma[j] = get_val(&stats);
	p->t0_inf[j] = get_val(&stats);
	
	for (unsigned int i=0; i<p->chan_idx[j]; i++) {
	  p->A[ichan] = get_val(&stats);
	  p->b[ichan] = get_val(&stats);
	  ichan ++;
	}
      }
      
      printf("Tau1 = %lf  DM = %lf\n", p->tau1, p->DM);

    }
  else
    {
      cout << "Unable to open file "<< filename << endl;
    }
  cout << "Finished reading stats file "<< filename << endl;
}
