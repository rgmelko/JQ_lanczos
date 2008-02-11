void expand (int stateNum, vector<bool> stateVec)
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

void collapse (int &stateNum, vector<bool> stateVec)
//collapses binary string into an integer
{
	stateNum=0;

	for (int i=(stateVec.size()-1); i>=0; i--)
	{
		stateNum*=2;
		stateNum+=stateVec.back();
		stateVec.pop_back();//eliminates element
	}

	//TODO: Check this at some point, to make sure it's appropriate to start from size-1 rather than size
}

double calc_Sz(vector<bool> stateVec)
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

void translate(vector<bool> stateVec, bool x,y; int n)
//translate statevector of nxn lattice right x and down y
//TODO generalize this to x,y>1
//	will require as many temp's as there are shifts in that direction
{
	bool temp; //for bit shifting

	if (x==true)
	{	//x-translations
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

int main();
{
	int dim; //dimension of lattice
	int Sz; //total Z-angular momentum to examine

	cout << "Enter dimension of square lattice" << endl;
	cin >> dim;
	//Defined states represented by integers 0..2^n-1

	cout << "Enter Spin-Z value." << endl;
	cin >> Sz;
	//This is total Sz, for an even number of sites, this is an integer

	vector<int> basis; //define vector to hold parent states in integer form
	vector<bool> test; //expanded version of an integer as an array of booleans
	vector<bool> temp; //temporary version of a test for comparisons
	vector<vector <bool>> daughters; //Vector of parent and daughter states

	bool found; //true if some translation of state i has been found in the basis set (so stop checking)

	for (int state=0; state<=pow(2,dim); state++)//cycle through all states
	{
		//clear daughters and reset found to false
		daughters.clear();
		found=false;

		expand(state, test); //convert state from integer to bit string

		if (calc(test)==Sz) //only pay attention to those in the correct spin sector
		{
			for (int i=0; i<dim; i++) //for the application of the 1-right operator to get to all possibilities
			{
				for (int j=0; j<dim; j++) //for the application of the 1-down operator to get to all possibilities
				{
					//all states will be in daughters, except parent state
					daughters.push_back(test);

					if (found==false)
					{
						for (int k=0; k<basis.size(); k++) //compare translated state to all basis vectors 
						{
							expand(basis[k], temp);
							if (daughters.back() == temp)
							{
								found=true;
							}
						}
					}
					translate(test, 0, 1, dim); //translate one down
				}
				translate(test, 1, 0, dim); //translate one right
			}

			//if this point is reached, and found!=true, then no translation of state is in the basis set so far
			if (found==false)
			{
				basis.push_back(state);
			}

			//Also if found=false, generate hamiltonian elements with all other states
		}
	}
}





