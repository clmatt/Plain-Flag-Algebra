using namespace std;
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


using namespace std::chrono;
using namespace mosek::fusion;
using namespace monty;


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
	auto start=high_resolution_clock::now();
	cout.precision(17);

	vector<Graph> f;
	vector<Graph> zeros;
	vector<Equation> known;
	vector<Edge> edges;
	const int numColors = 2;
	
	//vector<Edge> C5edges = {{0,1,1},{1,2,1},{2,3,1},{3,4,1},{4,0,1}};
	
	//Graph C5(C5edges,5,numColors);
	
	//f.push_back(C5);
	
	//Graph K2({{0,1,1}},2,numColors);
	
	//K2.setCoefficient(Frac(-1,1));
	//Equation K2Dens({K2},zeros,Frac(-9,10),1);
	//known.push_back(K2Dens);
	
	Graph K3({{0,1,1},{0,2,1},{1,2,1}},3,numColors);
	Graph K3c({{}},3,numColors);
	
	f.push_back(K3);
	f.push_back(K3c);
	
	//Minimize for false
	plainFlagAlgebra(f,7,zeros,known,false);
	
	auto end=high_resolution_clock::now();
	auto duration = duration_cast<seconds>(end - start);
	cout << "Running time in seconds: " << duration.count() << endl << endl;
	
	return 0;
}

