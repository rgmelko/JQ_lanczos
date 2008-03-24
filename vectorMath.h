//Vector Operations Library
//Tim Foster
//Library for inner product definition

#ifndef VECTORMATH
#define VECTORMATH


bool innerProductStates(vector <bool> aV1, vector <bool> aV2)
//returns a 1 if the states are identical
{
	bool result=true;

	for (int i=0; i<aV1.size(); i++)
	{
		if (aV1[i]!=aV2[i])
		{result=false;}
	}
	return result;
}

double GENHAM::hamElementGenerator(vector <bool> state1, double coeff1, vector <bool> state2, double coeff2)
//wrapper for inner product across a matrix
{
	double element=0; //running total of the matrix element
	double tempElement=0; //increment added by each bond found
	int site1, site2; //index of two elements under consideration

	vector <bool> temp1;
	vector <bool> temp2;

//	loop over all bonds to right and down
	for (int i=0; i<Nsite; i++)
	{
		for (int j=0; j<Nsite; j++)
		{
			//calculate spin-z component - always there for states identical, just a sign issuee
			//cout << "Considering interactions of site " << i << j << endl;
			for (int dir=0; dir<2; dir++)//adding dir*Nsite will move to right or dowwn site
			{
				site1=i*Nsite+j;
				//cout << "Site 1 = " << site1 << " with bit equal to " << state1[site1] << endl;
				if (dir==0)
				{
					//cout << "Horizonal bond" << endl;
					site2=i*Nsite+( (j+1)%Nsite );//mod wraps it around
				}
				else
				{
					//cout << "Vertical bond" << endl;
					site2=( (i+1)%Nsite )*Nsite+j;
				}
				//cout << "Site 2 = " << site2 << " with bit equal to " << state1[site2] << endl;
				
				//xor them, such that true output gives 1 (if states are parallel) and false gives -1, term only survives if states are equal
				//element+= -(2*(state1[site1]^state1[site2])-1)*equal(state1,state2);
				tempElement= -(2*(state1[site1]^state1[site2])-1)*equal(state1,state2);
				element+=tempElement;

				//cout << "From Sz overlap term, " << tempElement << endl;
				tempElement=0;

				//bond flipping term
				//test to see if up can be made down on (i,j) and partner
				if ( (state1[site1]==1) && (state1[site2]==0))
				{
					//implement switching, get prefactor 1/2
					//create temp1 the switched state
					//cout << "They are switchable." << endl;
					temp1.assign(state1.begin(),state1.end());
					temp1[site1]=0;
					temp1[site2]=1;
					tempElement=0.5*equal(state2,temp1)*coeff1*coeff2;
					element+=tempElement;
				}
				//testing opposite
				else if ( (state1[site1]==0) && (state1[site2]==1))
				{
					//implement switching, get prefactor 1/2
					//create temp2 the switched state
					//cout << "They are switchable." << endl;
					temp1.assign(state1.begin(),state1.end());
					temp1[site1]=1;
					temp1[site2]=0;
					tempElement=0.5*equal(state2,temp1)*coeff1*coeff2;
					element+=tempElement;
				}
				//cout << "Kinetic term added: " << tempElement << endl << endl;
				//cout << "Total element is " << element << endl << endl;
			}
				//otherwise done
//TODO still need to corrent translation operator for only one-site hop in
		}
	}
	//element*=JJ;//multiply the whole thing by -J parameter
	return element;
}

#endif

