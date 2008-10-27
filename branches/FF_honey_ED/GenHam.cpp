#include "GenHam.h"

//----------------------------------------------------------
GENHAM::GENHAM(const int Ns, const h_float J_, const int Sz)  
               : JJ(J_)
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

    for (int T0=0; T0<Nsite; T0++){ // Now generate off-diagonal part

      si = Bond(T0,0);
      if (si != T0) cout<<"Square error \n";
      // X Bond
      tempod = tempi;
      sj = Bond(T0,1);
      tempod ^= (1<<si);   //toggle bit 
      tempod ^= (1<<sj);   //toggle bit 
      revPos = revBas.at(tempod);
      if (revPos != -1){
        //tempD = (*this).HOFFdBond0(T0,tempi);
        tempD = -0.5; //FM in plane exchange
        Ham(ii,revPos) = tempD;
      }

      // Y Bond
      tempod = tempi;
      sj = Bond(T0,2);
      if (sj != -99) {  //if the second bond exists
        tempod ^= (1<<si);   //toggle bit 
        tempod ^= (1<<sj);   //toggle bit 
        revPos = revBas.at(tempod);
        if (revPos != -1){
          //tempD = (*this).HOFFdBond1(T0,tempi);
          tempD = -0.5; //FM in plane exchange
          Ham(ii,revPos) = tempD;
        }
      }


    }//si

  }//ii

}//FullHamHeis


//----------------------------------------------------------
double GENHAM::HdiagPart(const long bra){

  int S0b,S1b ;  //spins (bra 
  int T0,T1,F1;  //site
//  int s0p, s1p, s2p, s3p;
  double valH = 0;

  for (int Ti=0; Ti<Nsite; Ti++){
    //***HEISENBERG PART

    T0 = Bond(Ti,0); //lower left spin
    S0b = (bra>>T0)&1;  
    if (T0 != Ti) cout<<"Square error 3\n";
    T1 = Bond(Ti,1); //first bond
    S1b = (bra>>T1)&1;  //unpack bra
    valH += JJ*(S0b-0.5)*(S1b-0.5);
    T1 = Bond(Ti,2); //second bond
    if (T1 != -99) { //if there exists a 2nd bond
      S1b = (bra>>T1)&1;  //unpack bra
      F1 = Bond(Ti,3);
      if (F1 == 0)
        valH += JJ*(S0b-0.5)*(S1b-0.5);
      else
        valH -= JJ*(S0b-0.5)*(S1b-0.5);
    }

  }//T0

  //cout<<bra<<" "<<valH<<endl;

  return valH;

}//HdiagPart

//----------------------------------------------------------
double GENHAM::HOFFdBond(const int si, const long bra){

  double valH;

  valH = JJ*0.5; //contribution from the J part of the Hamiltonian


  return valH;

}//HOFFdPart



