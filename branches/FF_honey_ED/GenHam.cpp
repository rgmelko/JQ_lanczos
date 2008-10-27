#include "GenHam.h"

//----------------------------------------------------------
GENHAM::GENHAM(const int Ns, const h_float J_, const h_float J2_, const int Sz)  
               : JJ(J_), J2(J2_)
//create bases and determine dim of full Hilbert space
{
  int Dim;
  Nsite = Ns;

  Dim = 2;  //  S=1/2 models : two states
  for (int ch=1; ch<Nsite; ch++) Dim *=2;
  Fdim = Dim;

  BasPos.resize(Dim,-1); //initialization 

  Vdim=0;
  unsigned long temp;    //create basis (16 site cluster)
  for (int i1=0; i1<Dim; i1++) 
  {
    temp = 0;
    for (int sp =0; sp<Nsite; sp++)
       temp += (i1>>sp)&1;  //unpack bra
    if (temp==(Nsite/2+Sz) ){ 
       Basis.push_back(i1);
       BasPos.at(i1)=Basis.size()-1;
       Vdim++;
     }
  }//Dim

  cout<<"Vdim "<<Vdim<<" "<<Dim<<endl;

}//constructor



//----------------------------------------------------------
void GENHAM::FullHamJQ(){

  int ii, tempI;
  vector<long> revBas(Fdim,-1);

  for (ii=0; ii<Basis.size(); ii++) { //reverse lookup
    tempI = Basis.at(ii);
    revBas.at(tempI) = ii;
  }
    
  Ham.resize(Vdim,Vdim);
  Ham = 0;

  unsigned long tempi, tempj, tempod;
  double tempD;
  int si,sj,sk,sl, revPos;
  for (ii=0; ii<Basis.size(); ii++){
    tempi = Basis.at(ii);

    //Hamiltonian for diagonal
    tempD = (*this).HdiagPart(tempi);
    Ham(ii,ii) = tempD;

//    for (int T0=0; T0<Nsite; T0++){ // Now generate off-diagonal part
//
//      //si = PlaqX(T0,0); //si = Bond(T0,0);
//      //if (si != T0) cout<<"Square error \n";
//      // X Bond
//      tempod = tempi;
//      //sj = PlaqX(T0,1); //sj = Bond(T0,1);
//      tempod ^= (1<<si);   //toggle bit 
//      tempod ^= (1<<sj);   //toggle bit 
//      revPos = revBas.at(tempod);
//      if (revPos != -1){
//        tempD = (*this).HOFFdBond0(T0,tempi);
//        Ham(ii,revPos) = tempD;
//      }
//
//      // Y Bond
//      tempod = tempi;
//      //sj = PlaqX(T0,3); //sj = Bond(T0,2);
//      tempod ^= (1<<si);   //toggle bit 
//      tempod ^= (1<<sj);   //toggle bit 
//      revPos = revBas.at(tempod);
//      if (revPos != -1){
//        tempD = (*this).HOFFdBond1(T0,tempi);
//        Ham(ii,revPos) = tempD;
//      }
//
//
//    }//si

  }//ii

}//FullHamHeis


//----------------------------------------------------------
double GENHAM::HdiagPart(const long bra){

  int S0b,S1b ;  //spins (bra 
  int T0,T1;  //site
  int P0, P1, P2, P3; //sites for plaquette (Q)
  int s0p, s1p, s2p, s3p;
  double valH = 0;

  for (int Ti=0; Ti<Nsite; Ti++){
    //***HEISENBERG PART

    //T0 = PlaqX(Ti,0); //T0 = Bond(Ti,0); //lower left spin
    S0b = (bra>>T0)&1;  
    //if (T0 != Ti) cout<<"Square error 3\n";
    //T1 = PlaqX(Ti,1); //T1 = Bond(Ti,1); //first bond
    S1b = (bra>>T1)&1;  //unpack bra
    valH += JJ*(S0b-0.5)*(S1b-0.5);
    //T1 = PlaqX(Ti,3); //T1 = Bond(Ti,2); //second bond
    S1b = (bra>>T1)&1;  //unpack bra
    valH += JJ*(S0b-0.5)*(S1b-0.5);

    //Next-Nearest Neighbor part
     //bond 0,2 
    //T0 = PlaqX(Ti,0); 
    S0b = (bra>>T0)&1;
    //T1 = PlaqX(Ti,2);
    S1b = (bra>>T1)&1; 
    valH += J2*(S0b-0.5)*(S1b-0.5);
     //bond 1,3
    //T0 = PlaqX(Ti,1); 
    S0b = (bra>>T0)&1;
    //T1 = PlaqX(Ti,3);
    S1b = (bra>>T1)&1;
    valH += J2*(S0b-0.5)*(S1b-0.5);


  }//T0

  //cout<<bra<<" "<<valH<<endl;

  return valH;

}//HdiagPart

//----------------------------------------------------------
double GENHAM::HOFFdBond0(const int si, const long bra){

  double valH;
  int S0, S1;
  int T0, T1;

  valH = JJ*0.5; //contribution from the J part of the Hamiltonian

//  T0 = OtherTwoX(si,0); //first other set (diagonal - Q) spin
//  T1 = OtherTwoX(si,1); 
//  S0 = (bra>>T0)&1;    //spin values base 0
//  S1 = (bra>>T1)&1;
//  valH += QQ*0.5*( (S0-0.5)*(S1-0.5) - 0.25);

//  T0 = OtherTwoX(si,2); //second other set (diagonal - Q) spin
//  T1 = OtherTwoX(si,3); 
//  S0 = (bra>>T0)&1;    //spin values base 0
//  S1 = (bra>>T1)&1;
//  valH += QQ*0.5*( (S0-0.5)*(S1-0.5) - 0.25);

  return valH;

}//HOFFdPart

//----------------------------------------------------------
double GENHAM::HOFFdBond1(const int si, const long bra){

  double valH;
  int S0, S1;
  int T0, T1;

  valH = JJ*0.5; //contribution from the J part of the Hamiltonian


  return valH;

}//HOFFdPart


