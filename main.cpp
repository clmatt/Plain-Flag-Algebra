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

using namespace std::chrono;


//My headers
#include "general.h"
#include "fractions.h"
#include "flags.h"

//#define MAXN 100    /* Define this before including nauty.h */

//Nauty
extern "C" {
	#include "nauty.h"   
	#include "naututil.h"
	#include "gtools.h"
}

//This is going to be extremely confusing, but graph is from Nauty and Graph is from me


int main() {		
	vector<Graph> f;
	vector<Graph> zeros;
	vector<Equation> known;
	
	//Graph G({{}},2,2);
	//known.push_back(Equation({G},zeros,Frac(1,1),1));
	
	f.push_back(Graph({{0,1,1},{0,2,1},{1,2,1}},3,2));
	f.push_back(Graph({{}},3,2));
	
	auto time1=high_resolution_clock::now();
	NEWplainFlagAlgebra(f,5,zeros,known);
	auto time2=high_resolution_clock::now();
	
	//auto time3=high_resolution_clock::now();
	//plainFlagAlgebra(f,7,zeros,known);
	//auto time4=high_resolution_clock::now();
	
	auto duration1 = duration_cast<milliseconds>(time2 - time1);
	//auto duration2 = duration_cast<milliseconds>(time4 - time3);
	
	cout << "The running time of NEWplainFlagAlgebra is: " << duration1.count() << endl << endl;
	//cout << "The running time of plainFlagAlgebra is: " << duration2.count() << endl << endl;
	
	return 0;
}

