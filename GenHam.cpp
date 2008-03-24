#include "GenHam.h"
#include <iostream>
#include <fstream>
#include "vectorMath.h"

//----------------------------------------------------------
//Added TF

void printLongVector(const vector <long>& aV)
{
	for (int i=0; i<aV.size(); i++)
	{cout << aV[i] << " " ;}
	cout << endl;
}
void GENHAM::printBoolVector(const vector <bool>& aV)
{
	for (int i=0; i<Nsite; i++)
	{
		for (int j=0; j<Nsite; j++)
		{cout << aV[i*Nsite+j] << " " ;}
		cout << endl;
	}
	cout << endl << endl;
}

void GENHAM::printBoolArray(const vector <vector<bool> >& aV, const vector<double>& aCoeffs)
//Prints a set of daughter states
{
	for (int i=0; i<aV.size(); i++)
	{
		cout << aCoeffs[i] << " | ";
		for (int j=0; j<aV[i].size(); j++)
		{cout << aV[i][j] << " ";}
		cout << " >" << endl;
	}
	//cout << endl;
}

void GENHAM::expandVec(int stateNum, vector <bool>& stateVec)
//expands integer into binary string
{
	bool digit; //selected digit

	stateVec.clear(); //to make sure no memory leaking

	for (int place=0; place<(Nsite*Nsite); place++) //dimension of expanded vector must be that of the hilbert space
	{
		digit = stateNum % 2;
		stateVec.push_back(digit);
		stateNum=(stateNum-digit)/2;
	}
}

void GENHAM::collapseVec(int& stateNum, vector <bool>& stateVec)
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

double GENHAM::calc_Sz(const vector <bool>& stateVec)
//returns sum of up-bits minus sum of down bits
{
	double total; //running total of Sz

	for (int i=0; i<stateVec.size(); i++)
	{total+=stateVec[i];}
	total-=0.5*stateVec.size();//effectively subtract half off each value

	return total;
}

void GENHAM::translate(int x, int y, int n, vector <bool>& stateVec, double& coeff, const vector<double>& K_sect)
//translate statevector of nxn lattice right x and down y
{
	bool temp; //for bit shifting

	do
	{
		if (x<0)
		{	x=Nsite+x;}
		if (y<0)
		{	y=Nsite+y;}
	}while ((x<0) || (y<0));

	//x-translations
		for (int i=0; i<x; i++)//perform translation x times
		{
			for (int row=0; row<n; row++)
			{
				temp=stateVec[n*row+n-1]; //save the last bit in the row
				for (int column=n-1; column>=(0+1); column--)//shift everything in the row right one bit
				{stateVec[n*row+column]=stateVec[n*row+column-1];}
				stateVec[n*row]=temp;
			}
			if (K_sect[0]==1)
			{	coeff*=(-1);}
		}

	//y-translations
		for (int i=0; i<y; i++)//perform translation x times
		{
			for (int column=0; column<n; column++)
			{
				temp=stateVec[n*(n-1)+column]; //save the last bit in each column
				for (int row=n-1; row>=(0+1); row--) //shift everything down one bit
				{stateVec[n*row+column]=stateVec[n*(row-1)+column];}
				stateVec[column]=temp;
			}
			if (K_sect[1]==1)
			{	coeff*=(-1);}
		}
}

bool GENHAM::equal(const vector <bool>& v1, const vector <bool>& v2)
//tests equality of v1 and v2
//TODO integrate comparison of coefficients too
{
	bool same=false;

	if (v1.size()==v2.size())
	{
		same=true;
		for (int i=0; i<v1.size(); i++)
		{
			if (v1[i]!=v2[i])
			{same=false;}
		}
	}
	return same;
}

void GENHAM::gen_daughters(long aBasis, vector <vector <bool> >& daughters, vector<double>& coeffs, const vector<double>& K_sect)
//accepts integer and finds all translated versions thereof
{
	daughters.clear();
	coeffs.clear();

	vector <bool> test;
	double testCoeff;
	expandVec(aBasis,test);
	testCoeff=1;
	bool foundInDaughters;

	for (int i=0; i<Nsite; i++) //for the application of the 1-right operator to get to all possibilities
	{
		for (int j=0; j<Nsite; j++) //for the application of the 1-down operator to get to all possibilities
		{
			//all states will be in daughters, including parent state
			//compare translated state to all previous daughter vectors
			foundInDaughters=false;

			for (int k=0; k<daughters.size(); k++)
			{
				if (equal(test, daughters[k]))
				{
					foundInDaughters=true;
				}

			}
			if (foundInDaughters==false)
			{
				daughters.push_back(test);
				coeffs.push_back(testCoeff);
			}

			if ((j+1)<Nsite) //to avoid unnessessary last rotation
			{	translate(0,1,Nsite,test,testCoeff,K_sect);} //translate one down
		}
		if ((i+1)<Nsite) //to avoid unnessessary last rotation
		{	translate(1,0,Nsite,test,testCoeff,K_sect);} //translate one right
	}
}

void GENHAM::scanTranslation(	bool& foundInBasis, vector<long> Basis, vector<vector <bool> >& daughters, vector<double>& dCoeffs,
										vector<bool>& test, double& testCoeff)
//Tests translated vector to see if its already in the basis and puts it in daughters
{
	vector<bool> temp;//temporary version of a test for comparisons

	bool foundInDaughters=false;//true if daughter state already in daughters array
	for (int k=0; k<daughters.size(); k++)
	{
		if (equal(test, daughters[k]))
		{foundInDaughters=true;}
	}
	if (foundInDaughters==false)
	{
		daughters.push_back(test);
		dCoeffs.push_back(testCoeff);
	}

	if (foundInBasis==false)
	{
		//compare translated state to all basis vectors
		for (int k=0; k<Basis.size(); k++)
		{
			expandVec(Basis[k],temp);//puts expansion into stateVec
			if (equal(test, temp))
			{foundInBasis=true;}
		}
	}
}
//End Added TF
//----------------------------------------------------------
GENHAM::GENHAM(const int Ns, const h_float J_, const h_float Q_, const double Sz, const vector<double>& K_sect)  
					: JJ(J_), QQ(Q_) 
//constructor, TF's
{
	Fdim=(int)pow(2.,(Ns*Ns));
	Nsite=Ns; //dimension of lattice, this is a slightly different implementation
				 //before, Nsite was the dimension of a state
	
	cout << "Hamiltonian Generator Starting." << endl;
	cout << "Important parameters: Ns: " << Ns << " Sz: " << Sz << " K_sect: (" << K_sect[0] << "," << K_sect[1] << ")" << endl;
	getchar();

	vector<long> Basis; //define vector to hold parent states in integer form
	vector<bool> test; //expanded version of an integer as an array of booleans
	double testCoeff; //coefficient of daughter function to be tested

	bool foundInBasis; //true if some translation of a state has been found in the basis set on a previous translation (so stop checking)
	
	vector<bool> temp;//temporary version of a test for comparisons

	//Ham; //initialize hamiltonian to zero
	double element; //running total of hamiltonian matrix element

	for (int state=0; state<Fdim; state++)//cycle through all states
	{
		//clear daughters and reset foundInBasis to false
		daughters.clear();
		dCoeffs.clear();
		foundInBasis=false;

		expandVec(state,test); //convert state from integer to bit string
		testCoeff=1;//assumed coefficient, will be inverted by some translations

		if (calc_Sz(test)==Sz) //only pay attention to those in the correct spin sector
		{
			//cout << "Original Vector" << endl;
			//printBoolVector(test);
			//getchar();
			for (int i=0; i<Nsite; i++) //for the application of the 1-right operator
			{
				for (int j=0; j<Nsite; j++) //for the application of the 1-down operator
				{
					//if unique among daughters, put it in, if found in the basis, rule out this one's parent
					scanTranslation(foundInBasis, Basis, daughters, dCoeffs, test, testCoeff);

					//cout << "Translating by y=" << j << endl;
					if ((j+1)<Nsite) //to avoid unnessessary last rotation
					{	translate(0,1,Nsite,test,testCoeff,K_sect);} //translate one down
					//printBoolVector(test);
				}
				//cout << "Translating by x=" << i << endl;
				if ((i+1)<Nsite) //to avoid unnessessary last rotation
				{	translate(1,0,Nsite,test,testCoeff,K_sect);} //translate one right
				//printBoolVector(test);
			}

			//if this point is reached, and found!=true, then no translation of state is in the basis set so far
			if (foundInBasis==false)
			{
				Basis.push_back(state);
				int BasisSize=Basis.size();

				//generate column of hamiltonian elements with all other states
	cout << "Before resizing: " << endl << Ham << endl;
				Ham.resizeAndPreserve(Basis.size(),Basis.size());//grow the array by one in each direction
	cout << "After resizing: " << endl << Ham << endl;
				//off diagonal elements
				for (int i=0; i<(Basis.size()-1); i++)// for all basis not including last
				{
					Ham(i,BasisSize)=0; //set element to 0
					gen_daughters(Basis[i],daughters1,d1Coeffs,K_sect);//generate daughters of second state, since first state was generated in discovery process, and not yet cleared
					for (int d1=0; d1<daughters.size(); d1++) //for all combinations of daughter and daughter1 members
					{
						for (int d2=d1; d2<daughters1.size(); d2++)
						{
							cout << "States : " << endl;
							printBoolVector(daughters1[d2]);
							printBoolVector(daughters[d1]);
							Ham(i,BasisSize)+=hamElementGenerator(daughters1[d2],d1Coeffs[d2],daughters[d1],dCoeffs[d1]);// add the contribution of this daughter-pair
							cout << "Should give hamElement: " << hamElementGenerator(daughters1[d2],d1Coeffs[d2],daughters[d1],dCoeffs[d1]) << " at " << i << "x" << Basis.size() << endl;
						}
					}
					Ham(BasisSize,i)=Ham(i,BasisSize); //duplicate in lower half of matrix
	cout << "After off-diagonal generation: " << endl << Ham << endl;
				}
				//diagonal element - set aside since daughters already generated
						int d2=0;
						printBoolVector(daughters[d2]);
						printBoolVector(daughters[0]);
				cout << "Fine until here." << endl;
				Ham(BasisSize,BasisSize)=0; //set element to 0
				cout << "Made it to here?" << endl;
				for (int d1=0; d1<daughters.size(); d1++) //for all combinations of daughter and daughter1 members
				{
					for (int d2=d1; d2<daughters.size(); d2++)
					{
						cout << "States : " << endl;
						cout << "Attempting to access daughters element d2=" << d2 << " where the size of daughters is " << daughters.size() << endl;
						printBoolVector(daughters[d2]);
						printBoolVector(daughters[d1]);
						cout << "Should give hamElement: " << hamElementGenerator(daughters[d2],d1Coeffs[d2],daughters[d1],dCoeffs[d1])<< " at " << state << "x" << state << endl;
						Ham(BasisSize,BasisSize)+=hamElementGenerator(daughters[d2],d1Coeffs[d2],daughters[d1],dCoeffs[d1]); //add the contribution of this daughter-pair
					}
				}
	cout << "After diagonal generation: " << endl << Ham << endl;

			}

		}
	}
//	cout << "At end: " << Ham << endl;
	
	/*printLongVector(Basis);
	for (int i=0; i<Basis.size(); i++)
	{
		cout << "Basis state " << Basis[i] << endl;
		expandVec(Basis[i],temp);
		printBoolVector(temp);
		cout << "Daughters1" << endl;
		gen_daughters(Basis[i],daughters1,d1Coeffs,K_sect);
		printBoolArray(daughters1,d1Coeffs);
	}*/

	//testing section for hamiltonian element generator

	/*cout << "Here's a quick attempt to test hamiltonian generator.  Test on " << Basis[0] << "&" << Basis[2] << endl;

	expandVec(Basis[0],temp);
	printBoolVector(temp);
	vector<bool> temp1;
	expandVec(Basis[2],temp1);
	printBoolVector(temp1);

	element=hamElementGenerator(temp, temp1);

	cout << "The result is " << element << endl;
	getchar();
	*/

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

/*void GENHAM::printg()
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
*/
#include "simparam.h"

int main()
//Added TF, so this would compile on its own without worrying about all the Lanczoos stuff yet
{
	//copied from ed_Lan_1107.cpp
	PARAMS prm;
	double J;
	double Q;
	vector <double> K_Sect(2); //Two element array of K-sector, in multiples of pi TODO Generalize to more than just (0,0),(0,pi),(pi,0),(pi,pi)
	int Sz;
 
	J=prm.JJ_;
	Q=-prm.QQ_;// SIGN HAS TO BE FLIPPED: NOW SIGN AT INPUT IS SAME AS PAPER 
	Sz=prm.Sz_;
	K_Sect[0]=0;//prm.K_SectX;
	K_Sect[1]=1;//prm.K_SectY;

	//slightly different implementation of dimensionality
	//first parameter is now the length of each dimension, so, the square root of the total enumber of elements
	GENHAM HV(2,J,Q,Sz,K_Sect);

}
