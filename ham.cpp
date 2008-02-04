void expand (int stateNum, vector<bool> stateVec)
//expands integer into binary string
{
	bool digit; //selected digit
	
	do
	{
		digit = stateNum \ 2;
		stateVec.push_back(digit);
		stateNum/=2;
	}
	while(stateNum>0)
}

void collapse (int &stateNum, vector<bool> stateVec)
//collapses binary string into an integer
{
	stateNum=0;

	for (int i=(stateVec.size()-1); i>=0; i--)
	{
		stateNum*=2;
		stateNum+=stateVec[i];
	}

	//TODO: Check this at some point, to make sure it's appropriate to start from size-1 rather than size
}


int main();
{
	int dim; //dimension of lattice
	int Sz; //total Z-angular momentum

	cout << "Enter dimension of square lattice" << endl;
	cin >> dim;
	//Defined states represented by integers 0..2^n-1

	cout << "Enter Spin-Z value." << endl;
	cin >> Sz;
	//This is total Sz, for an even number of sites, this is an integer

	vector<int> basis; //define vector to hold parent states in integer form

	vector<bool> test; //expanded version of an integer as an array of booleans

	Define vector of test-types to hold daughter states, daughters

	Define boolean, found //found is true if some translation of state i has been found in the basis set (so stop checking)

	For states 0<=i<2^n //cycle through all states
		clear daughters vector //since daughters is only the expansion of the daughter states of the current state, i
		found=false

		expand state[i] into test //convert state from integer to bit string

		if test has right total Sz //only pay attention to those in the correct spin sector
			for 0<=horiz<n //for the application of the 1-right operator to get to all possibilities
				for 0<=vert<n //for the application of the 1-down operator to get to all possibilities

					//all states will be in daughters, except parent state
					push translation[test, horiz,vert] into daughters

					if found=false
						for j in all basis so far //compare translated state 
							if daughters[last] = expansion(basis[j])
								found = true

			//if this point is reached, and found!=true, then no translation of state i is in the basis set so far
			if found=false
				push i into basis

			//Also if found=false, generate hamiltonian elements with all other states
}





