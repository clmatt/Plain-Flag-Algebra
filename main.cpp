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
	std::cout.precision(17);

	std::vector<Graph> f;
	std::vector<Graph> zeros;
	std::vector<Equation> known;
	std::vector<Edge> edges;
	const int numColors = 2;
	
	//vector<Edge> C5edges = {{0,1,1},{1,2,1},{2,3,1},{3,4,1},{4,0,1}};
	
	//Graph C5(C5edges,5,numColors);
	
	//f.push_back(C5);
	
	Graph K3({{0,1,1},{0,2,1},{1,2,1}},3,numColors);
	Graph K3c({{}},3,numColors);
	
	f.push_back(K3);
	f.push_back(K3c);
	
	Graph K2({{0,1,1}},2,numColors);
	
	K2.setCoefficient(Frac(-1,1));
	Equation K2Dens({K2},zeros,Frac(-9,10),0);
	known.push_back(K2Dens);
	
	Equation K3cDens({K3c},zeros,Frac(1,100),0);
	known.push_back(K3cDens);
	
	//Minimize for false
	plainFlagAlgebra(f,6,zeros,known,false);
	
	auto end=std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	std::cout << "Running time in seconds: " << duration.count() << std::endl << std::endl;
	
	return 0;
}

