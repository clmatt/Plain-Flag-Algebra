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
	
	//R(4,4) <= 17
	vector<Graph> f;
	vector<Graph> zeros;
	vector<Equation> known;
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
	
	//K4
	edges.clear();
	for(int i = 0; i < 3; ++i) {
		for(int j = i+1; j < 4; ++j) {
			edges.push_back({i,j,2});
		}
	}
	zeros.push_back(Graph(edges,4,numColors));
	
	//f-Maximize number of edges
	f.push_back(Graph({{0,1,2}},2,numColors));
	
	//Make known
	//CHANGE if running on more vertices
	for(int i = 2; i <= 3; ++i) {
		Graph Kempty({{}},i,numColors);
		Equation knownEmpty({Kempty},zeros,Frac(1,myPow(17,i-1)),0);
		known.push_back(knownEmpty);
	}
	
	//Edge density bounds based on R(3,4)
	Graph K21({{0,1,1}},2,numColors);
	Equation known1({K21},zeros,Frac(8,17),1);
	known.push_back(known1);
	
	K21.setCoefficient(Frac(-1,1));
	Equation known2({K21},zeros,Frac(-8,17),1);
	known.push_back(known2);
	
	//Equation doesn't play nicely with constants or single vertex flags
	vector<Graph> vOne; //Identically equal to one when we have a single flag
	vOne.push_back(Graph({{}},2,numColors));
	vOne[0].setFlag({0});
	vOne.push_back(Graph({{0,1,1}},2,numColors));
	vOne[1].setFlag({0});
	vOne.push_back(Graph({{0,1,2}},2,numColors));
	vOne[2].setFlag({0});
	
	Graph K21f({{0,1,1}},2,numColors);
	K21f.setFlag({0});
	
	vector<Graph> deg8Vec = {K21f};
	for(int i = 0; i < vOne.size(); ++i) {
		vOne[i].setCoefficient(Frac(-8,17));
		deg8Vec.push_back(vOne[i]);
	}
	
	Equation deg8(deg8Vec,zeros,Frac(0,1),0);
	
	Equation known3 = deg8*deg8;
	known3.averageAll();
	known.push_back(known3);
	
	
	plainFlagAlgebra(f,6,zeros,known);
	
	
	auto end=high_resolution_clock::now();
	auto duration = duration_cast<seconds>(end - start);
	cout << "Running time in seconds: " << duration.count() << endl << endl;
	
	return 0;
}
