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
	auto start=std::chrono::high_resolution_clock::now();
	//omp_set_num_threads(1);	
	//R(4,6) <= 40
	std::vector<Graph> f;
	std::vector<Graph> zeros;
	std::vector<Equation> known;
	std::vector<Edge> edges;
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
			edges.push_back({i,j,2});
		}
	}
	zeros.push_back(Graph(edges,4,numColors));
	
	//K6
	edges.clear();
	for(int i = 0; i < 5; ++i) {
		for(int j = i+1; j < 6; ++j) {
			edges.push_back({i,j,1});
		}
	}
	zeros.push_back(Graph(edges,6,numColors));
	
	//f-Maximize number of non-edges
	f.push_back(Graph({{0,1,1}},2,numColors));
	
	//Make known
	//CHANGE if running on more vertices
	for(int i = 2; i <= 5; ++i) {
		Graph Kempty({{}},i,numColors);
		Kempty.setCoefficient(Frac(myPow(40,i-1),1));
		Equation knownEmpty({Kempty},zeros,Frac(1,1),0);
		known.push_back(knownEmpty);
	}
	
	//Edge density bounds based on R(4,5),R(3,6)
	Graph K21({{0,1,1}},2,numColors);
	Equation known1({K21},zeros,Frac(24,40),1);
	known.push_back(known1);
	
	K21.setCoefficient(Frac(-1,1));
	Equation known2({K21},zeros,Frac(-22,40),1);
	known.push_back(known2);
	
	Graph K21f({{0,1,1}},2,numColors,{0});

	std::vector<Graph> vOne; //Identically equal to one when we have a single flag
	vOne.push_back(Graph({{}},2,numColors,{0}));
	vOne.push_back(Graph({{0,1,1}},2,numColors,{0}));
	vOne.push_back(Graph({{0,1,2}},2,numColors,{0}));
	
	K21f.setCoefficient(Frac(40,1));

	//Degree 22
	std::vector<Graph> deg22Vec = {K21f};
	for(int i = 0; i < vOne.size(); ++i) {
		vOne[i].setCoefficient(Frac(-22,1));
		deg22Vec.push_back(vOne[i]);
	}
	Equation deg22(deg22Vec,zeros,Frac(0,1),0);
	
	//Degree 23
	std::vector<Graph> deg23Vec = {K21f};
	for(int i = 0; i < vOne.size(); ++i) {
		vOne[i].setCoefficient(Frac(-23,1));
		deg23Vec.push_back(vOne[i]);
	}
	Equation deg23(deg23Vec,zeros,Frac(0,1),0);
	
	//Degree 24
	std::vector<Graph> deg24Vec = {K21f};
	for(int i = 0; i < vOne.size(); ++i) {
		vOne[i].setCoefficient(Frac(-24,1));
		deg24Vec.push_back(vOne[i]);
	}
	Equation deg24(deg23Vec,zeros,Frac(0,1),0);
	
	Equation known3 = deg22*deg23*deg24;
	known3.averageAll();
	known.push_back(known3);
	
	Equation known4 = deg22*deg23*deg24*deg22;
	known4.averageAll();
	known.push_back(known4);
	
	Equation known5 = deg22*deg23*deg24*deg23;
	known5.averageAll();
	known.push_back(known5);
	
	Equation known6 = deg22*deg23*deg24*deg24;
	known6.averageAll();
	known.push_back(known6);
	
	Equation known7 = deg22*deg23*deg24*deg22*deg23;
	known7.averageAll();
	known.push_back(known7);
	
	Equation known8 = deg22*deg23*deg24*deg22*deg24;
	known8.averageAll();
	known.push_back(known8);
	
	Equation known9 = deg22*deg23*deg24*deg23*deg23;
	known9.averageAll();
	known.push_back(known9);

	
	fastPlainFlagAlgebra(f,8,zeros,known,false);
	
	
	auto end=std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	std::cout << "Running time in seconds: " << duration.count() << std::endl << std::endl;
	
	return 0;
}

