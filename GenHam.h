//c++ class for creating a general Hamiltonian in c++ vectors
//Roger Melko, November 2007
#ifndef GenHam_H
#define GenHam_H

#include <iostream>
#include <vector>
using namespace std;

#include <blitz/array.h>
BZ_USING_NAMESPACE(blitz)

typedef float h_float;  //precision for Hamiltonian storage

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
    Array<int,2> OtherTwoX; //sites on the plaquette that are not in bond
    Array<int,2> OtherTwoY;
    Array<int,2> PlaqX;
    Array<int,2> PlaqY;
    
    Array<double,2> Ham;  //full hamiltonian

    //Added TF
	 vector< vector<bool> > daughters; //vector of a parent and its daughter states
	 vector<long> stateVec; //workspace for functions, since I can't figure out how to pass vectors

    bool equal(); //poor workaround, see GenHam.cpp for details

	 void expandVec(int);
	 void collapseVec(int&);
	 double calc_Sz();
	 void translate(bool, bool, int);
	 //End of Added TF

    GENHAM(const int,const h_float J_, const h_float Q_,const int Sz); 
    void printg();
    //double at(const int , const int );
    Array<double,1> apply(const Array<double,1>&);

    void Bonds_16A();
    void Bonds_16B();
    void Bonds_18A();
    void Bonds_20A();
    void Bonds_24A();
    void Bonds_24R();
    void Bonds_26A();
    void SparseHamJQ();
    void FullHamJQ();

  private:
    int Nsite;

    h_float JJ; //heisenberg exchange value
    h_float QQ; //ring-exchange value

    double HdiagPart(const long);
    double HOFFdBondX(const int, const long);
    double HOFFdBondY(const int, const long);
    double HOFFdPlaq(const int, const long);

};

#endif
