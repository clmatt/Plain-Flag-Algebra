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
	omp_set_num_threads(1);
	auto start=std::chrono::high_resolution_clock::now();
	
	std::vector<Graph> f;
	std::vector<Graph> zeros;
	std::vector<Equation> known;
	
	/*Graph K4({{0,1,1},{0,2,1},{0,3,1},{1,2,1},{1,3,1},{2,3,1}},4,2);
	Graph G3({{0,1,1},{1,2,1},{0,2,1},{2,3,1},{3,0,1}},4,2);
	G3.setCoefficient(Frac(1,2));
	Graph G2({{0,1,1},{1,2,1},{2,3,1},{3,0,1}},4,2);
	G2.setCoefficient(Frac(1,3));
	Graph G1({{0,1,1},{1,2,1},{0,2,1},{2,3,1}},4,2);
	G1.setCoefficient(Frac(1,6));
	Graph P4({{0,1,1},{1,2,1},{2,3,1}},4,2);
	P4.setCoefficient(Frac(1,12));
	Graph E({{0,1,1}},2,2);
	
	Equation mainEQ({K4,G3,G2,G1,P4});
	Equation Eeq({E});
	
	Equation ans = mainEQ + (-1)*Eeq*Eeq*Eeq;
	
	for(int i = 0; i < ans.getNumVariables(); ++i) {
		f.push_back(ans.getVariable(i));
	}*/
	
	Graph K3({{0,1,1},{0,2,1},{1,2,1}},3,2);
	Graph K3c({},3,2);
	
	f = {K3,K3c};
	
	//Max = true
	plainFlagAlgebra(f,4,zeros,known,false);
	
	auto end=std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	std::cout << "Running time in seconds: " << duration.count() << std::endl << std::endl;
	
	return 0;
}

