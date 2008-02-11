#include "GenHam.h"

//----------------------------------------------------------
//Added TF
void GENHAM::expandVec(int stateNum)
//expands integer into binary string
{
	bool digit; //selected digit

	stateVec.clear(); //to make sure no memory leaking

	do
	{
		digit = stateNum % 2;
		stateVec.push_back(digit);
		stateNum/=2;
	}
	while(stateNum>0);
}

void GENHAM::collapseVec(int& stateNum)
//collapses binary string into an integer
{
	stateNum=0;

	for (int i=(stateVec.size()-1); i>=0; i--)
	{
		stateNum*=2;
		stateNum+=stateVec.back();
		stateVec.pop_back();//eliminates element
	}
}

double GENHAM::calc_Sz()
//returns sum of up-bits minus sum of down bits
{
	double total; //running total of Sz

	for (int i=0; i<stateVec.size(); i++)
	{
		total+=stateVec[i];
	}
	total-=0.5*stateVec.size();//effectively subtract half off each value

	return total;
}

void GENHAM::translate(bool x, bool y, int n)
//translate statevector of nxn lattice right x and down y
//TODO generalize this to x,y>1
// will require as many temp's as there are shifts in that direction
{
	bool temp; //for bit shifting

	if (x==true)
	{  //x-translations
		for (int row=0; row<n; row++)
		{
			temp=stateVec[n*row+n-1]; //save the last bit in the row
			for (int column=n-1; column>=(0+1); column--)//shift everything in the row right one bit
			{
				stateVec[n*row+column]=stateVec[n*row+column-1];
			}
			stateVec[n*row]=temp;
		}
	}

	if (y==true)
	{
		//y-translations
		for (int column=0; column<n; column++)
		{
			temp=stateVec[n*(n-1)+column]; //save the last bit in each column
			for (int row=n-1; row>=(0+1); row--) //shift everything down one bit
			{
				stateVec[n*row+column]=stateVec[n*(row-1)+column];
			}
			stateVec[column]=temp;
		}
	}
}
//TODO Vector operators = and == are not working
//this is a workaround for the latter
bool GENHAM::equal()
{
	bool same=true;

	for (int i=0; i<(daughters.back())[i]; i++)
	{
		if (((daughters.back())[i])!=stateVec[i])
		{
			same=false;
		}
	}
	return same;
}

//End Added TF
//----------------------------------------------------------
GENHAM::GENHAM(const int Ns, const h_float J_, const h_float Q_, const int Sz)  
					: JJ(J_), QQ(Q_) 
//constructor, TF's
{
	Nsite=Ns; //dimension of lattice

	vector<int> Basis; //define vector to hold parent states in integer form
	vector<bool> test; //expanded version of an integer as an array of booleans
	vector<bool> temp; //temporary version of a test for comparisons

	bool found; //true if some translation of a state has been found in the basis set on a previous translation (so stop checking)

	for (int state=0; state<=pow(2.,Nsite); state++)//cycle through all states
	{
		//clear daughters and reset found to false
		daughters.clear();
		found=false;

		expandVec(state); //convert state from integer to bit string
		test.assign(stateVec.begin(), stateVec.end());

		if (calc_Sz()==Sz) //only pay attention to those in the correct spin sector
		{
			for (int i=0; i<Nsite; i++) //for the application of the 1-right operator to get to all possibilities
			{
				for (int j=0; j<Nsite; j++) //for the application of the 1-down operator to get to all possibilities
				{
					//all states will be in daughters, except parent state
					daughters.push_back(test);

					if (found==false)
					{
						for (int k=0; k<Basis.size(); k++) //compare translated state to all basis vectors
						{
							expandVec(Basis[k]);//puts expansion into stateVec
							if (equal())
							{
								found=true;
							}
						}
					}
					translate(0, 1, Nsite); //translate one down
				}
				translate(1, 0, Nsite); //translate one right
			}

			//if this point is reached, and found!=true, then no translation of state is in the basis set so far
			if (found==false)
			{
				Basis.push_back(state);
			}

			//Also if found=false, generate hamiltonian elements with all other states
		}
	}
}
/*//create bases and determine dim of full Hilbert space
{
  int Dim;
  Nsite = Ns;

  Dim = 2;  //  S=1/2 models : two states
  for (int ch=1; ch<Nsite; ch++) Nsite *=2;
  FNsite = Nsite;

  BasPos.resize(Nsite,-1); //initialization 

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

//  cout<<"Vdim "<<Vdim<<" "<<Dim<<endl;

}//constructor*/


//----------------------------------------------------------
void GENHAM::printg()
{
  int i,j;
  vector<int> tempP;
  vector<h_float> tempV;

  for (i=0; i<PosHam.size(); i++){
    //cout<<PosHam[i][0]<<" * ";
    cout<<i+1<<" * ";
    for (j=0; j<=PosHam[i][0]; j++){
     cout<<"("<<PosHam[i][j]+1<<","<<ValHam[i][j]<<") ";
   }
   cout<<endl;
  }

}//print


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

      si = PlaqX(T0,0); //si = Bond(T0,0);
      //if (si != T0) cout<<"Square error \n";
      // X Bond
      tempod = tempi;
      sj = PlaqX(T0,1); //sj = Bond(T0,1);
      tempod ^= (1<<si);   //toggle bit 
      tempod ^= (1<<sj);   //toggle bit 
      revPos = revBas.at(tempod);
      if (revPos != -1){
        tempD = (*this).HOFFdBondX(T0,tempi);
        Ham(ii,revPos) = tempD;
      }

      // Y Bond
      tempod = tempi;
      sj = PlaqX(T0,3); //sj = Bond(T0,2);
      tempod ^= (1<<si);   //toggle bit 
      tempod ^= (1<<sj);   //toggle bit 
      revPos = revBas.at(tempod);
      if (revPos != -1){
        tempD = (*this).HOFFdBondY(T0,tempi);
        Ham(ii,revPos) = tempD;
      }

      // Square Plaquette
      tempod = tempi;
      si = PlaqX(T0,0);
      sj = PlaqX(T0,1);
      sk = PlaqX(T0,2);
      sl = PlaqX(T0,3);
      tempod ^= (1<<si);   //toggle bits 
      tempod ^= (1<<sj);   
      tempod ^= (1<<sk);   
      tempod ^= (1<<sl);   
      revPos = revBas.at(tempod);
      if (revPos != -1){ 
        tempD = (*this).HOFFdPlaq(T0,tempi);
        Ham(ii,revPos) = tempD;
      }

    }//si

  }//ii

}//FullHamHeis

//----------------------------------------------------------
void GENHAM::SparseHamJQ()
{
  int ii, jj;

  int Rsize;
  vector<long> tempBas;
  //vector<long> tempBas2;
  vector<h_float> tempH;
  unsigned long tempi, tempj, tempod;
  int si, sj,sk,sl;
  double tempD;

  for (ii=0; ii<Basis.size(); ii++){
    tempH.clear(); 
    tempBas.clear();

    tempi = Basis.at(ii);
    tempBas.push_back(0); //first element (Row size)
    tempH.push_back(0); //make base 0

    //-----1:   diagonal 
    tempBas.push_back(BasPos.at(tempi));  
    tempD = (*this).HdiagPart(tempi);
    tempH.push_back(tempD); 

    for (int T0=0; T0<Nsite; T0++){ //T0 is your square index

      si = PlaqX(T0,0); //si = Bond(T0,0); //the lower left bond spin is not always T0
      //if (si != T0) cout<<"Square error 2\n";
      //-----2:   first bond (Horizontal)
      tempod = tempi;
      sj = PlaqX(T0,1); //sj = Bond(T0,1);
      tempod ^= (1<<si);   //toggle bit 
      tempod ^= (1<<sj);   //toggle bit 
      if (BasPos.at(tempod) != -1 && BasPos.at(tempod) > ii){ //build only upper half of matrix
        tempBas.push_back(BasPos.at(tempod));
        tempD = (*this).HOFFdBondX(T0,tempi);
        tempH.push_back(tempD); 
      }
 
       //-----3:   second bond (Vertical)
      tempod = tempi;
      sj = PlaqX(T0,3); //sj = Bond(T0,2);
      tempod ^= (1<<si);   //toggle bit 
      tempod ^= (1<<sj);   //toggle bit 
      if (BasPos.at(tempod) != -1 && BasPos.at(tempod) > ii){ 
        tempBas.push_back(BasPos.at(tempod));
        tempD = (*this).HOFFdBondY(T0,tempi);
        tempH.push_back(tempD); 
      }

       //-----4:   plaquette 
      tempod = tempi;
      si = PlaqX(T0,0);   //not always redundant here
      sj = PlaqX(T0,1);  
      sk = PlaqX(T0,2);
      sl = PlaqX(T0,3);
      tempod ^= (1<<si);   //toggle bits 
      tempod ^= (1<<sj);   
      tempod ^= (1<<sk);   
      tempod ^= (1<<sl);   
      if (BasPos.at(tempod) != -1 && BasPos.at(tempod) > ii){ 
        tempBas.push_back(BasPos.at(tempod));
        tempD = (*this).HOFFdPlaq(T0,tempi);
        tempH.push_back(tempD); 
      }

    }//si

    tempBas.at(0) = tempBas.size()-1;
    //cout<<tempBas.at(0)<<" "<<tempBas.size()<<" "<<tempH.size()<<endl;

    //bubble sort (slow) 
    long stemp;
    bool noswap = false;
    while (noswap == false){
      noswap = true; 
      for (int i2=1; i2<tempBas.size()-1; i2++){ //ignore 0 element
        if (tempBas.at(i2) > tempBas.at(i2+1) ) {
          stemp = tempBas.at(i2);
          tempBas.at(i2) = tempBas.at(i2+1);
          tempBas.at(i2+1) = stemp;
          tempD = tempH.at(i2);
          tempH.at(i2) = tempH.at(i2+1);
          tempH.at(i2+1) = tempD;
          noswap = false;
        }
      }//i2
    }//while

    PosHam.push_back(tempBas);
    ValHam.push_back(tempH);

  }//ii       

}//Heisenberg

//----------------------------------------------------------
double GENHAM::HdiagPart(const long bra){

  int S0b,S1b ;  //spins (bra 
  int T0,T1;  //site
  int P0, P1, P2, P3; //sites for plaquette (Q)
  int s0p, s1p, s2p, s3p;
  double valH = 0;

  for (int Ti=0; Ti<Nsite; Ti++){
    //***HEISENBERG PART

    T0 = PlaqX(Ti,0); //T0 = Bond(Ti,0); //lower left spin
    S0b = (bra>>T0)&1;  
    //if (T0 != Ti) cout<<"Square error 3\n";
    T1 = PlaqX(Ti,1); //T1 = Bond(Ti,1); //first bond
    S1b = (bra>>T1)&1;  //unpack bra
    valH += JJ*(S0b-0.5)*(S1b-0.5);
    T1 = PlaqX(Ti,3); //T1 = Bond(Ti,2); //second bond
    S1b = (bra>>T1)&1;  //unpack bra
    valH += JJ*(S0b-0.5)*(S1b-0.5);

    //X Plaquettes
    P0 = PlaqX(Ti,0);  //if (P0 != Ti) cout<<"ERROR \n";
    s0p = (bra>>P0)&1;
    P1 = PlaqX(Ti,1); 
    s1p = (bra>>P1)&1;
    P2 = PlaqX(Ti,2); 
    s2p = (bra>>P2)&1;
    P3 = PlaqX(Ti,3); 
    s3p = (bra>>P3)&1;
    valH += QQ* ((s0p-0.5)*(s1p-0.5) - 0.25) * ((s2p-0.5)*(s3p-0.5) - 0.25);

    //Y Plaquettes  -rotate PlaqX 
    P0 = PlaqX(Ti,1);
    s0p = (bra>>P0)&1;
    P1 = PlaqX(Ti,2); 
    s1p = (bra>>P1)&1;
    P2 = PlaqX(Ti,3); 
    s2p = (bra>>P2)&1;
    P3 = PlaqX(Ti,0); 
    s3p = (bra>>P3)&1;
    valH += QQ* ((s0p-0.5)*(s1p-0.5) - 0.25) * ((s2p-0.5)*(s3p-0.5) - 0.25);

  }//T0

  //cout<<bra<<" "<<valH<<endl;

  return valH;

}//HdiagPart

//----------------------------------------------------------
double GENHAM::HOFFdBondX(const int si, const long bra){

  double valH;
  int S0, S1;
  int T0, T1;

  valH = JJ*0.5; //contribution from the J part of the Hamiltonian

  T0 = OtherTwoX(si,0); //first other set (diagonal - Q) spin
  T1 = OtherTwoX(si,1); 
  S0 = (bra>>T0)&1;    //spin values base 0
  S1 = (bra>>T1)&1;
  valH += QQ*0.5*( (S0-0.5)*(S1-0.5) - 0.25);

  T0 = OtherTwoX(si,2); //second other set (diagonal - Q) spin
  T1 = OtherTwoX(si,3); 
  S0 = (bra>>T0)&1;    //spin values base 0
  S1 = (bra>>T1)&1;
  valH += QQ*0.5*( (S0-0.5)*(S1-0.5) - 0.25);

  return valH;

}//HOFFdPart

//----------------------------------------------------------
double GENHAM::HOFFdBondY(const int si, const long bra){

  double valH;
  int S0, S1;
  int T0, T1;

  valH = JJ*0.5; //contribution from the J part of the Hamiltonian

  T0 = OtherTwoY(si,0); //first other set (diagonal - Q) spin
  T1 = OtherTwoY(si,1); 
  S0 = (bra>>T0)&1;    //spin values base 0
  S1 = (bra>>T1)&1;
  valH += QQ*0.5*( (S0-0.5)*(S1-0.5) - 0.25);

  T0 = OtherTwoY(si,2); //second other set (diagonal - Q) spin
  T1 = OtherTwoY(si,3); 
  S0 = (bra>>T0)&1;    //spin values base 0
  S1 = (bra>>T1)&1;
  valH += QQ*0.5*( (S0-0.5)*(S1-0.5) - 0.25);

  return valH;

}//HOFFdPart

//----------------------------------------------------------
double GENHAM::HOFFdPlaq(const int si, const long bra){

  int S0, S1, S2, S3;
  int T0, T1, T2, T3;

  double valH=0.0;

  T0 = PlaqX(si,0);   //if (T0 != si) cout<<"ERROR \n";
  S0 = (bra>>T0)&1;
  T1 = PlaqX(si,1);
  S1 = (bra>>T1)&1; // 
  T2 = PlaqX(si,2); //   3  2 
  S2 = (bra>>T2)&1; //   0  1
  T3 = PlaqX(si,3); // 
  S3 = (bra>>T3)&1;

  //X-contribution to energy
  //if ( (S0 != S1) && (S2 != S3) ) valH += QQ*0.25;
  if  (S0 != S1) {
    if  (S2 == S3) cout<<"P1 error \n";
    valH += QQ*0.25;
   }
  //Y-contribution to energy
  if (S0 != S3) {
    if (S2 == S1) cout<<"P2 error \n";
    valH += QQ*0.25;
  }

  return valH;

}//HOFFdPart

int main()
//Added TF, so this would compile on its own without worrying about all the Lanczoos stuff yet
{
	//Testing routing for new GenHam class
}
