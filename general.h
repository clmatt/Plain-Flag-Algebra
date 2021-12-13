//---------------------------
//---------------------------
//-----General Functions-----
//---------------------------
//---------------------------

//For permutations use stl::next_permutation

//---------------------
//-----Next Subset-----
//---------------------

//Requires to be a subset from 0,1,...,n-1 of size k
//Initialize with 0,...,k-1
bool nextSubset(vector<int> &subset, int n, int k) {
	int index,i;
	bool val;

	val = true;
	index = k-1;

	while (val == true) {

		if (index < 0) {
			return false;
			val = false;
		}

		else if (subset[index] < n - k + index) {
			val = false;

			for (i = index; i < k; ++i) {
				if (i == index) {
					++subset[i];
				}
				else {
					subset[i] = subset[i - 1] + 1;
				}
			}
		}

		else if (subset[index] > n - k + index) {
			cout << "Error in nextSubset" << endl << endl;
			throw exception();
		}

		else {
			index = index - 1;
		}
	}
	
	return true;
}


//-------------------
//-----Next List-----
//-------------------

//maxVals is the maximum each element can take in the list
//unlike subset we allow repeats
//Intialize list with all zeros
bool nextList(vector<int> &list, const vector<int> &maxVals) {
	if(list.size() != maxVals.size()) {
		cout << "In next_list the list size must be the same as maxVals size." << endl;
		cout << "List size = " << list.size() << endl;
		cout << "maxVals size = " << maxVals.size() << endl << endl;
		throw exception();
	}
	
	for(int i = 0; i < (int)list.size(); ++i) {
		if(list[i] > maxVals[i]) {
			cout << "In next_list the list entries must be <= the maxVal entries." << endl;
			cout << "Index = " << i << endl;
			cout << "List entry = " << list[i] << endl;
			cout << "maxVals entry = " << maxVals[i] << endl << endl;
			throw exception();
		}
		
		if(maxVals[i] < 0) {
			cout << "In next_list, maxVals can't have negative entries." << endl;
			cout << "Index = " << i << endl;
			cout << "maxVals entry = " << maxVals[i] << endl << endl;
			throw exception();
		}
		
		if(list[i] < 0) {
			cout << "In next_list, list can't have negative entries." << endl;
			cout << "Index = " << i << endl;
			cout << "maxVals entry = " << list[i] << endl << endl;
			throw exception();
		}
	}
	
	bool output = false;
	
	for(int i = 0; i < (int)list.size(); ++i) {
		if(list[i] == maxVals[i]) {
			list[i] = 0;
		}
		
		else {
			list[i] = list[i] + 1;
			i = list.size();
			output = true;
		}
	}
	
	return output;
}

//----------------------------
//-----Euclid's Algorithm-----
//----------------------------

//Returns positive number
long long int gcd(long long int a, long long int b) {
	if(a < 0) {
		a = -a;
	}

	if(b < 0) {
		b = -b;
	}

	while (b != 0) {
		int s = b;
		b = a % b;
		a = s;
	}
	
	return a;
}


//-------------------------------
//-----Least Common Multiple-----
//-------------------------------

long long int lcm(long long int a, long long int b) {
	if(a < 0) {
		a = -a;
	}

	if(b < 0) {
		b = -b;
	}
	
	if(b > a) {
		return a*(b/gcd(a,b));
	}
	
	return b*(a/gcd(a,b));
}


//------------------
//-----My Power-----
//------------------

//C++ has a stupid power function, this is better if everything is positive integers
long long int myPow(long long int a, int b) {
	long long int val = a;

	if (b < 0) {
		cout << "Only using myPow if power is >= 0." << endl << endl;
		throw exception();
	}
	
	else if(b == 0) {
		return 1;
	}
	
	else { 
		for(int i = 0; i < b-1; ++i) {
			val = val*a;
		}
	}
	
	return val;
}

//----------------
//-----Choose-----
//----------------

long long int choose(long long int n, long long int k) {
	if(n < 0 || k < 0) {
		cout << "Choose can't take negative numbers." << endl;
		cout << "n = " << n << endl;
		cout << "k = " << k << endl;
		throw exception();
	}
	
	if(n < k) {
		return 0;
	}
	
	if(k == 0) {
		return 1;
	}
	
	if(n == k) {
		return 1;
	}
	
	return choose(n-1,k-1) + choose (n-1,k);
	
}

//----------------------
//-----Binary Digit-----
//----------------------

//Find the kth binary digit in i 
int binaryDigit(int n, int k) {
	return ((n & (1 << (k - 1))) >> (k - 1));
}


//-------------------
//-----Factorial-----
//-------------------

int factorial(int n) {
	if(n < 0) {
		cout << "Factorial doesn't take negative numbers." << endl;
		cout << "n = " << n << endl << endl;
		throw exception();
	}
	
	if((n == 1) || (n == 0)) {
		return 1;
	}
	
	return n*factorial(n-1);
}


