//c++ class for creating a general Hamiltonian in c++ vectors
//Roger Melko, November 2007
#ifndef GenHam_H
#define GenHam_H

#include <iostream>
#include <vector>
using namespace std;

#include <blitz/array.h>
BZ_USING_NAMESPACE(blitz)

typedef long double h_float;  //precision for Hamiltonian storage

class GENHAM{

  public:
    int Fdim; //"full" Hilbert space
    int Vdim; //dimenson of reduced Hilbert space

    vector<vector<long> > PosHam;
    vector<vector<h_float> > ValHam;
    //vector<double> DiagHam;

    vector<long> Basis;
    vector<long> BasPos;

    Array<int,2> Bond;
    
    Array<double,2> Ham;  //full hamiltonian

    GENHAM(const int,const h_float J_, const h_float J2_,const int Sz); 
    
    Array<double,1> apply(const Array<double,1>&);

    void Bonds_8();
    void SparseHamJQ();
    void FullHamJQ();

  private:
    int Nsite;

    h_float JJ; //heisenberg exchange value
    h_float J2; //exchange value for bond b

    double HdiagPart(const long);
    double HOFFdBond0(const int, const long);
    double HOFFdBond1(const int, const long);
    double HOFFdBond2(const int, const long);

};

#endif
