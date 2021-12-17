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
	auto start=high_resolution_clock::now();
	
	
	
	//R(4,6) <= 39
	/*vector<Graph> f;
	vector<Graph> zeros;
	vector<Equation> known;
	vector<Graph> totalZeros;
	vector<Edge> edges;
	const int numColors = 3;
	
	//zeros
	edges = {{0,1,1}};
	zeros.push_back(Graph(edges,3,numColors));
	
	edges = {{0,1,2}};
	zeros.push_back(Graph(edges,3,numColors));
	
	edges = {{0,1,1},{1,2,2}};
	zeros.push_back(Graph(edges,3,numColors));
	
	//K4
	edges.clear();
	for(int i = 0; i < 3; ++i) {
		for(int j = i+1; j < 4; ++j) {
			edges.push_back({i,j,1});
		}
	}
	zeros.push_back(Graph(edges,4,numColors));
	
	//K6
	edges.clear();
	for(int i = 0; i < 5; ++i) {
		for(int j = i+1; j < 6; ++j) {
			edges.push_back({i,j,2});
		}
	}
	zeros.push_back(Graph(edges,6,numColors));
	
	//Total Zeros
	edges = {{0,1,1}}; 
	totalZeros.push_back(Graph(edges,3,2));
	
	//f-Maximize number of edges
	f.push_back(Graph({{0,1,2}},2,numColors));
	
	//Make known
	//CHANGE if running on more vertices
	for(int i = 2; i <= 6; ++i) {
		Graph Kempty({{}},i,numColors);
		Equation knownEmpty({Kempty},zeros,Frac(1,myPow(39,i-1)),0);
		known.push_back(knownEmpty);
	}
	
	//Edge density bounds based on R(3,6)
	Graph K21({{0,1,1}},2,numColors);
	Equation known1({K21},zeros,Frac(17,39),1);
	known.push_back(known1);
	
	K21.setCoefficient({-1});
	Equation known2({K21},zeros,Frac(-14,39),1);
	known.push_back(known2);
	
	Graph K21f({{0,1,1}},2,numColors);
	K21f.setFlag({0});
	Graph K11f({{}},1,numColors);
	K11f.setFlag({0});
	K11f.setCoefficient(Frac(-14,39));
	Equation eq1({K21f,K11f},zeros,Frac(0,1),0);
	K11f.setCoefficient(Frac(-15,39));
	Equation eq2({K21f,K11f},zeros,Frac(0,1),0);
	K11f.setCoefficient(Frac(-16,39));
	Equation eq3({K21f,K11f},zeros,Frac(0,1),0);
	K11f.setCoefficient(Frac(-17,39));
	Equation eq4({K21f,K11f},zeros,Frac(0,1),0);
	
	Equation known3 = eq1*eq1*eq2*eq3*eq4;
	known3.averageAll();
	known.push_back(known3);*/
	
	
	/*Equation known4 = eq1*eq2*eq2*eq3*eq4;
	known4.averageAll();
	known.push_back(known4);
	Equation known5 = eq1*eq2*eq3*eq3*eq4;
	known5.averageAll();
	known.push_back(known5);
	Equation known6 = eq1*eq2*eq3*eq4*eq4;
	known6.averageAll();
	known.push_back(known6);*/
	
	//plainFlagAlgebra(f,7,zeros,known,totalZeros);
	
	
	auto end=high_resolution_clock::now();
	auto duration = duration_cast<seconds>(end - start);
	cout << "Running time in seconds: " << duration.count() << endl << endl;
	
	return 0;
}

