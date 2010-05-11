/*******************************************************************
Exact Diagonalization program (Lanczos and Complete)
for the fully frustrated honeycomb XXZ model
Roger Melko, October 2008
********************************************************************/

#define DO_LANCZOS //Lanczos or LAPACK full ED

#include <fstream> 
#include <vector> 
#include <math.h>
using namespace std;

#include <iostream>
#include <limits.h>
#include <stdio.h>
#include <time.h>
#include <iomanip>
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

#include "Lanczos_07.h"
#include "GenHam.h"
#include "lapack.h"
#include "simparam.h"

int main(){

//  cout<<sizeof(long double)<<endl;
//  cout<<sizeof(float)<<endl;

  PARAMS prm;
  double J;  //bonds a
  int Sz;
  
  J=prm.JJ_;
  Sz=prm.Sz_;

  GENHAM HV(18,J,Sz); 
  HV.Bonds_18p(); 

#ifdef DO_LANCZOS    //EVERYTHING IN THIS BLOCK OF CODE FOR LANCZOS

  LANCZOS lancz(HV.Vdim);  //dimension of reduced Hilbert space (Sz sector)

  HV.SparseHam();  //generates sparse matrix Hamiltonian for Lanczos
  //HV.printg();
  lancz.Diag(HV,prm.Neigen_,prm.valvec_); // second parameter: # of eigenvalues to converge
                      // third parameter: 1 for -values only, 2 for vals AND vectors

#else             //EVERTHING IN THIS BLOCK OF CODE FOR FULL LAPACK ED

  HV.FullHamJQ();  //generates full Sz sector Hamiltonian

  cout<<"Done making the Hamiltonian \n";
  cout.flush();
  vector<double> dd;
  diagWithLapack_R(HV.Ham,dd);  //*** LAPACK DIAG
  cout<<"Done diagonalizting the Hamiltonian \n";
  cout.flush();

  double min=999.0;
  for (int ii=0; ii<dd.size(); ii++) {
    cout<<setprecision(12)<<dd.at(ii)<<"\n";     
    if (dd.at(ii)<min) min=dd.at(ii);
  }
  cout<<endl;
  cout<<"HH energy: "<<min<<endl;

#endif

  return 0;

}
