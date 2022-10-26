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
bool nextSubset(std::vector<int> &subset, int n, int k) {
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
			std::cout << "Error in nextSubset" << std::endl << std::endl;
			throw std::exception();
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
bool nextList(std::vector<int> &list, const std::vector<int> &maxVals) {
	if(list.size() != maxVals.size()) {
		std::cout << "In nextList the list size must be the same as maxVals size." << std::endl;
		std::cout << "List size = " << list.size() << std::endl;
		std::cout << "maxVals size = " << maxVals.size() << std::endl << std::endl;
		throw std::exception();
	}
	
	for(int i = 0; i < (int)list.size(); ++i) {
		if(list[i] > maxVals[i]) {
			std::cout << "In nextList the list entries must be <= the maxVal entries." << std::endl;
			std::cout << "Index = " << i << std::endl;
			std::cout << "List entry = " << list[i] << std::endl;
			std::cout << "maxVals entry = " << maxVals[i] << std::endl << std::endl;
			throw std::exception();
		}
		
		if(maxVals[i] < 0) {
			std::cout << "In nextList, maxVals can't have negative entries." << std::endl;
			std::cout << "Index = " << i << std::endl;
			std::cout << "maxVals entry = " << maxVals[i] << std::endl << std::endl;
			throw std::exception();
		}
		
		if(list[i] < 0) {
			std::cout << "In nextList, list can't have negative entries." << std::endl;
			std::cout << "Index = " << i << std::endl;
			std::cout << "maxVals entry = " << list[i] << std::endl << std::endl;
			throw std::exception();
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


//--------------------------
//-----Next Ordered Sum-----
//--------------------------

//Outputs all ordered sets of size k such that they sum to equal n
//Intialize with n,0,0,...,0 goes until 0,0,...,0,n
bool nextOrderedSum(std::vector<int> &sum) {
	int k = sum.size();
	
	int n = 0;
	for(int i = 0; i < k; ++i) {
		n += sum[i];
	}
	
	bool output = false;
	for(int i = 0; i < k-1; ++i) {
		if(sum[i] != 0) {
			output = true;
			i = k-1;
		}
	}
	
	if(!output) {
		return output;
	}
	
	if(k == 2) {
		--sum[0];
		++sum[1];
	}
	
	else {
		std::vector<int> partialSum(k-1);
		for(int i = 0; i < k-1; ++i) {
			partialSum[i] = sum[i+1];
		}
		
		bool test = nextOrderedSum(partialSum);
		
		if(test) {
			for(int i = 0; i < k-1; ++i) {
				sum[i+1] = partialSum[i];
			}
		} 
		
		else {
			--sum[0];
			sum[1] = n - sum[0];
			for(int i = 2; i < k; ++i) {
				sum[i] = 0;
			}	
		}
	}
	
	return true;
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
	
	long long int ans = b*(a/gcd(a,b));
	
	if(ans < 0) { 
		std::cout << "Overflow in lcm." << std::endl << std::endl;
		throw std::exception();
	}
	
	return ans;
}


//------------------
//-----My Power-----
//------------------

//C++ has a stupid power function, this is better if everything is positive integers
long long int myPow(long long int a, int b) {
	long long int val = a;

	if (b < 0) {
		std::cout << "Only using myPow if power is >= 0." << std::endl << std::endl;
		throw std::exception();
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
		std::cout << "Choose can't take negative numbers." << std::endl;
		std::cout << "n = " << n << std::endl;
		std::cout << "k = " << k << std::endl;
		throw std::exception();
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

long long int factorial(int n) {
	if(n < 0) {
		std::cout << "Factorial doesn't take negative numbers." << std::endl;
		std::cout << "n = " << n << std::endl << std::endl;
		throw std::exception();
	}
	
	if((n == 1) || (n == 0)) {
		return 1;
	}
	
	return n*factorial(n-1);
}


//---------------------
//-----Multinomial-----
//---------------------

long long int multinomial(const int n, const std::vector<int> &input) {
	if(n < 0) {
		std::cout << "In multinomial n >= 0. " << std::endl << std::endl;
		throw std::exception();
	}
	
	if(input.size() == 0) {
		std::cout << "In multinomial, input.size() must be > 0." << std::endl << std::endl;
		throw std::exception();
	}
	
	int tot = 0;
	for(int i = 0; i < input.size(); ++i) {
		if(input[i] < 0) {
			std::cout << "In multinomial, everything in vector must be >= 0." << std::endl << std::endl;
			throw std::exception();
		}
	
		tot = tot + input[i];
	}
	
	if(tot != n) {
		std::cout << "In multinomial, input must sum to n." << std::endl << std::endl;
		throw std::exception(); 
	}
	
	long long int output = 1;
	long long int sum = input[0];
	
	for(int i = 1; i < input.size(); ++i) {
		sum += input[i];
		output = output*choose(sum,input[i]);
		
		if (output < 0) {
			std::cout << "Overflow in multinomial." << std::endl << std::endl;
			throw std::exception();
		}
	}
	
	
	
	return output;
}


//--------------------------
//-----Hash for Vectors-----
//--------------------------

//Who knows why C++ doesn't have built in hash for vectors
//Stolen for stackoverflow
#include <limits>
#include <cstdint>

template<typename T>
T xorshift(const T& n,int i){
  return n^(n>>i);
}

// a hash function with another name as to not confuse with std::hash
uint32_t distribute(const uint32_t& n){
  uint32_t p = 0x55555555ul; // pattern of alternating 0 and 1
  uint32_t c = 3423571495ul; // random uneven integer constant; 
  return c*xorshift(p*xorshift(n,16),16);
}

// a hash function with another name as to not confuse with std::hash
uint64_t distribute(const uint64_t& n){
  uint64_t p = 0x5555555555555555ull; // pattern of alternating 0 and 1
  uint64_t c = 17316035218449499591ull;// random uneven integer constant; 
  return c*xorshift(p*xorshift(n,32),32);
}

// if c++20 rotl is not available:
template <typename T,typename S>
typename std::enable_if<std::is_unsigned<T>::value,T>::type
constexpr rotl(const T n, const S i){
  const T m = (std::numeric_limits<T>::digits-1);
  const T c = i&m;
  return (n<<c)|(n>>((T(0)-c)&m)); // this is usually recognized by the compiler to mean rotation, also c++20 now gives us rotl directly
}

// call this function with the old seed and the new key to be hashed and combined into the new seed value, respectively the final hash
template <class T>
inline size_t hash_combine(std::size_t& seed, const T& v)
{
    return rotl(seed,std::numeric_limits<size_t>::digits/3) ^ distribute(std::hash<T>{}(v));
}

//Could also do this for containers which aren't vectors
template <typename T> 
struct vectorHash {
    std::size_t operator()(const std::vector<T> &vec) const {
        std::size_t seed = vec.size();
        for(auto& i : vec) {
			 seed = hash_combine(seed, i);
		  }
		  return seed;
    }
};


