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

//Doesn't work yet!
int main() {
	auto start=high_resolution_clock::now();
		
	//R(4,6) <= 39
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
	
	//K6
	edges.clear();
	for(int i = 0; i < 5; ++i) {
		for(int j = i+1; j < 6; ++j) {
			edges.push_back({i,j,2});
		}
	}
	zeros.push_back(Graph(edges,6,numColors));
	
	//f-Maximize number of edges
	f.push_back(Graph({{0,1,2}},2,numColors));
	
	//Make known
	//CHANGE if running on more vertices
	for(int i = 2; i <= 6; ++i) {
		Graph Kempty({{}},i,numColors);
		Kempty.setCoefficient(Frac(myPow(39,i-1),1));
		Equation knownEmpty({Kempty},zeros,Frac(1,1),0);
		known.push_back(knownEmpty);
	}
	
	//Edge density bounds based on R(3,6)
	Graph K21({{0,1,1}},2,numColors);
	Equation known1({K21},zeros,Frac(17,39),1);
	known.push_back(known1);
	
	K21.setCoefficient(Frac(-1,1));
	Equation known2({K21},zeros,Frac(-14,39),1);
	known.push_back(known2);
	
	Graph K21f({{0,1,1}},2,numColors);
	K21f.setFlag({0});

	vector<Graph> vOne; //Identically equal to one when we have a single flag
	vOne.push_back(Graph({{}},2,numColors));
	vOne[0].setFlag({0});
	vOne.push_back(Graph({{0,1,1}},2,numColors));
	vOne[1].setFlag({0});
	vOne.push_back(Graph({{0,1,2}},2,numColors));
	vOne[2].setFlag({0});
	
	K21f.setCoefficient(Frac(39,1));
	
	//Degree 14
	vector<Graph> deg14Vec = {K21f};
	for(int i = 0; i < vOne.size(); ++i) {
		vOne[i].setCoefficient(Frac(-14,1));
		deg14Vec.push_back(vOne[i]);
	}
	Equation deg14(deg15Vec,zeros,Frac(0,1),0);

	//Degree 15
	vector<Graph> deg15Vec = {K21f};
	for(int i = 0; i < vOne.size(); ++i) {
		vOne[i].setCoefficient(Frac(-15,1));
		deg15Vec.push_back(vOne[i]);
	}
	Equation deg15(deg15Vec,zeros,Frac(0,1),0);
	
	//Degree 16
	vector<Graph> deg16Vec = {K21f};
	for(int i = 0; i < vOne.size(); ++i) {
		vOne[i].setCoefficient(Frac(-16,1));
		deg16Vec.push_back(vOne[i]);
	}
	Equation deg16(deg16Vec,zeros,Frac(0,1),0);
	
	//Degree 17
	vector<Graph> deg17Vec = {K21f};
	for(int i = 0; i < vOne.size(); ++i) {
		vOne[i].setCoefficient(Frac(-17,1));
		deg17Vec.push_back(vOne[i]);
	}
	Equation deg17(deg17Vec,zeros,Frac(0,1),0);
	
	Equation known3 = deg14*deg14*deg15*deg16*deg17;
	known3.averageAll();
	known.push_back(known3);
	
	Equation known4 = deg14*deg15*deg15*deg16*deg17;
	known4.averageAll();
	known.push_back(known4);
	
	Equation known5 = deg14*deg15*deg16*deg16*deg17;
	known5.averageAll();
	known.push_back(known5);
	
	Equation known6 = deg14*deg15*deg16*deg17*deg17;
	known6.averageAll();
	known.push_back(known6);
	
	
	plainFlagAlgebra(f,6,zeros,known);
	
	
	auto end=high_resolution_clock::now();
	auto duration = duration_cast<seconds>(end - start);
	cout << "Running time in seconds: " << duration.count() << endl << endl;
	
	return 0;
}

