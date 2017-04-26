#include "string.h"

#include <iostream>
#include <algorithm> 
#include <fstream>
#include <string>
#include <iterator>
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

  int i, ipar=0,ichan, nmodes;
  string line;
  char filename[256];
  double maxlike = -1.0*pow(10.0,10);
  double *tmp_cube = new double [p->ndims];

  if (p->sampler==0) // MultiNest
    sprintf(filename, "%s/chainsMN-phys_live.points", root);
  if (p->sampler==1) // PolyChord
    sprintf(filename, "%s/chainsPC_phys_live.txt", root);

  //TODO
  ifstream physlive_file(filename);

  if (!physlive_file.is_open()) {cerr << "Error opening: " << filename << endl; return;}

  while (getline(physlive_file,line)) {
    std::istringstream myStream( line );
    std::istream_iterator< double > begin(myStream),eof;
    std::vector<double> paramlist(begin,eof);
    
    double like = paramlist[p->ndims+1];
    if(like > maxlike) {
      maxlike = like;
      for (int j=0; j< p->ndims; j++) 
	tmp_cube[j] = paramlist[j];
    }
  } 

  ichan = 0;
  // Ask to output the result
  p->do_plot = 1;

  p->tau1 = tmp_cube[ipar]; ipar++;
  p->DM = tmp_cube[ipar]; ipar++;
            
  for (unsigned int j = 0; j < p->nfiles; j++) {
    p->sigma[j] = tmp_cube[ipar];ipar++;
    p->t0_inf[j] = tmp_cube[ipar];ipar++;
    
    for (unsigned int i=0; i<p->chan_idx[j]; i++) {
      p->A[ichan] = tmp_cube[ipar];ipar++;
      p->b[ichan] = tmp_cube[ipar];ipar++;
      ichan ++;
    }
  }
  
  printf(" ML results: Tau1 = %lf  DM = %lf\n", p->tau1, p->DM);

  delete[] tmp_cube;
  cout << "Finished reading stats file "<< filename << endl;
}
