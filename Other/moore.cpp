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
		
	//Moore Graph	- color 2 are original edges
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
	
	//Girth 5
	edges = {{0,1,2},{1,2,2},{2,0,2}};
	zeros.push_back(Graph(edges,3,numColors));
	
	edges = {{0,1,2},{1,2,2},{2,3,2},{3,0,2},{0,2,1},{1,3,1}};
	zeros.push_back(Graph(edges,4,numColors));
	
	//f-Maximize number of non-triangles
	f.push_back(Graph({{0,1,1},{0,2,1},{1,2,1}},3,numColors));
	
	//Make known
	//CHANGE if running on more vertices
	for(int i = 2; i <= 4; ++i) {
		Graph Kempty({{}},i,numColors);
		Kempty.setCoefficient(Frac(myPow(3250,i-1),1));
		Equation knownEmpty({Kempty},zeros,Frac(1,1),0);
		known.push_back(knownEmpty);
	}
	
	//Know edge density
	Graph K21({{0,1,2}},2,numColors);
	Equation known1({K21},zeros,Frac(57,3250),1);
	known.push_back(known1);
	
	K21.setCoefficient(Frac(-1,1));
	Equation known2({K21},zeros,Frac(-57,3250),1);
	known.push_back(known2);
	
	Graph K21f({{0,1,2}},2,numColors,{0});

	vector<Graph> vOne; //Identically equal to one when we have a single flag
	vOne.push_back(Graph({{}},2,numColors,{0}));
	vOne.push_back(Graph({{0,1,1}},2,numColors,{0}));
	vOne.push_back(Graph({{0,1,2}},2,numColors,{0}));
	
	K21f.setCoefficient(Frac(3250,1));

	//Degree
	vector<Graph> degVec = {K21f};
	for(int i = 0; i < vOne.size(); ++i) {
		vOne[i].setCoefficient(Frac(-57,1));
		degVec.push_back(vOne[i]);
	}
	Equation deg(degVec,zeros,Frac(0,1),0);

	
	Equation known3 = deg*deg;
	known3.averageAll();

	known.push_back(known3);

	//False = minimize
	plainFlagAlgebra(f,6,zeros,known,true);
	
	
	auto end=high_resolution_clock::now();
	auto duration = duration_cast<seconds>(end - start);
	cout << "Running time in seconds: " << duration.count() << endl << endl;
	
	return 0;
}

