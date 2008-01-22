/*******************************************************************
Exact Diagonalization program (Lanczos and Complete)
for AWS' J-Q model
Roger Melko, November 2007
********************************************************************/


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
  double J;
  double Q;
  int Sz;
  
  J=prm.JJ_;
  Q=-prm.QQ_;// SIGN HAS TO BE FLIPPED: NOW SIGN AT INPUT IS SAME AS PAPER 
  Sz=prm.Sz_;

  GENHAM HV(24,J,Q,Sz); 
  HV.Bonds_24R(); 

  LANCZOS lancz(HV.Vdim);  //dimension of reduced Hilbert space (Sz sector)

  HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos
//  HV.printg();
  lancz.Diag(HV,prm.Neigen_,prm.valvec_); // second parameter: # of eigenvalues to converge
                      // third parameter: 1 for -values only, 2 for vals AND vectors

//  HV.FullHamJQ();  //generates full Sz sector Hamiltonian

//  Array<double,1> dd(HV.Vdim);         //************
//  Array<double,1> ee(HV.Vdim);         //  RGM HOUSEHOLDER
//  H.tred3(HV.Ham, dd, ee, HV.Vdim);    //************

//  vector<double> dd;
//  diagWithLapack_R(HV.Ham,dd);  //*** LAPACK DIAG
//
//  double min=999.0;
//  for (int ii=0; ii<dd.size(); ii++) {
//    cout<<setprecision(12)<<dd.at(ii)<<"\n";     
//    if (dd.at(ii)<min) min=dd.at(ii);
//  }
//  cout<<endl;
//  cout<<"HH energy: "<<min<<endl;

  return 0;

}