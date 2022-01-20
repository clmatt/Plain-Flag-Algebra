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
	auto time1=high_resolution_clock::now();		
	
	vector<Graph> f;
	vector<Graph> zeros;
	vector<Equation> known;
	
	f.push_back(Graph({{0,1,1},{0,2,1},{1,2,1}},3,2));
	f.push_back(Graph({{}},3,2));

	
	plainFlagAlgebra(f,6,zeros,known,false);
	
	
	auto time2=high_resolution_clock::now();
	auto duration1 = duration_cast<milliseconds>(time2 - time1);
	cout << "The running time is: " << duration1.count() << endl << endl;

	return 0;
}

