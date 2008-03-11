//Vector Operations Library
//Tim Foster
//Library for inner product definition

#ifndef VECTORMATH
#define VECTORMATH


bool innerProductStates(vector <bool> aV1, vector <bool> aV2)
//returns a 1 if the states are identical
{
	bool result=true;

	for (int i=0; i<aV1.size; i++)
	{
		if (aV1[i]!=aV2[i])
		{result=false;}
	}
	return result;
}

double GENHAM::hamElementGenerator(two states)
//wrapper for inner product across a matrix
	double element=0; //running total of the matrix element
	int site2; //index of the next element

	vector <bool> temp1;
	vector <bool> temp2;

//calculate spin-z component - always there for states identical, just a sign issuee
//	loop over all bonds to right and down
	for (int i=0; i<Nsite; i++)
	{
		for (int j=0; j<Nsite; j++)
		{
			for (int dir=0; dir<2; dir++)//adding dir*Nsite will move to right or dowwn site
			{
				site1=i*Nsite+j;
				if (dir==0)
				{
					site2=i*Nsite+( (j+1)%Nsite );//mod wraps it around
				}
				else
				{
					site2=( (i+1)%Nsite )*Nsite+j;
				}
				//xor them, such that true output gives 1 and false gives -1
				element+=2*state1[site1]^state1[site2]-1;

				//bond flipping term
				//test to see if up can be made down on (i,j) and partner
				if ( (state1[site1]==1) && (state1[site2]==0))
				{
					//implement switching, get prefactor 1/2
					//create temp1 the switched state
					temp1.assign(state1.begin(),state1.end());
					temp1[site1]=0;
					temp1[site2]=1;
					element+=0.5;
				}
				//testing opposite
				else if ( (state1[site1]==0) && (state1[site2]==1))
				{
					//implement switching, get prefactor 1/2
					//create temp2 the switched state
					temp2.assign(state1.begin(),state1.end());
					temp2[site1]=1;
					temp2[site2]=0;
					element+=0.5;
				}
				//otherwise done
//TODO I don't think this is a proper implementation.  Investigate further
//TODO still need to corrent translation operator for only one-site hop in
//Hilbert space reduction
		}
	}
//multiply the whole thing by -J parameter


#endif

