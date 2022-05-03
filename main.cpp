//using namespace std;
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <limits> 
#include <queue>
#include <tuple>
#include <fusion.h>
#include <omp.h>
#include <threads.h>
#include <boost/functional/hash.hpp>



//using namespace std::chrono;
//using namespace mosek::fusion;
//using namespace monty;


//My headers
#include "general.h"
#include "fractions.h"
#include "flags.h"


//Nauty
extern "C" {
	#include "nauty.h"   
	#include "naututil.h"
	#include "gtools.h"
}

//This is going to be extremely confusing, but graph is from Nauty and Graph is from me


int main() {
	omp_set_num_threads(1); //Use for debugging
	auto start=std::chrono::high_resolution_clock::now();
	
	std::vector<Graph> f;
	std::vector<Graph> zeros;
	std::vector<Equation> known;
	
	f.push_back(Graph({{0,1,1},{1,2,1},{0,2,1}},3,2));
	f.push_back(Graph({{}},3,2));
	
	//plainFlagAlgebraApprox(f,7,1,zeros,known,false);
	plainFlagAlgebra(f,5,zeros,known,false);

	auto end=std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	std::cout << "Running time in seconds: " << duration.count() << std::endl << std::endl;
	
	return 0;
}

