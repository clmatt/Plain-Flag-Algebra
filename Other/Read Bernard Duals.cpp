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
		
	ifstream myFile;
	myFile.open("dual.txt", ios::in);
	
	string fileLine;
	vector<Graph> allGraphs;
	vector<double> density;
	const int numColors = 3;
	vector<Edge> edges;
	
	do {
		double densTemp;
		myFile >> densTemp;
		
		density.push_back(densTemp);
		
		int dummy;
		myFile >> dummy;
		myFile >> dummy;
		
		int edgeColor;
		for(int i = 0; i < 7; ++i) {
			for(int j = i+1; j < 8; ++j) {
				myFile >> edgeColor;
				edges.push_back({i,j,edgeColor-1});
			}
		}
		
		Graph G(edges,8,numColors);
		allGraphs.push_back(G);
		
		edges.clear();
		
		getline(myFile,fileLine);
	} while(!myFile.eof());
	
	vector<Graph> zeros;
	
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
	
	Graph K21f({{0,1,1}},2,numColors,{0});

	vector<Graph> vOne; //Identically equal to one when we have a single flag
	vOne.push_back(Graph({{}},2,numColors,{0}));
	vOne.push_back(Graph({{0,1,1}},2,numColors,{0}));
	vOne.push_back(Graph({{0,1,2}},2,numColors,{0}));
	
	K21f.setCoefficient(Frac(40,1));
	
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
	
	Equation degree = deg15*deg16*deg15*deg16*deg17;
	degree.averageAll();
	
	double total = 0.; 
	
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		cout << i << " out of " << allGraphs.size() << endl;
		cout << "Total = " << total << endl << endl;
		for(int j = 0; j < degree.getNumVariables(); ++j) {
			Graph H = degree.getVariable(j);
		
			int mult = numSubgraphs(H,allGraphs[i]);
			//cout << mult << endl;
			total += (H.getCoefficient().getNum()*density[i]*(double)mult)/((double)choose(8,H.getN())*H.getCoefficient().getDen());  
		}
	}
	
	cout << total << endl << endl;
	
	myFile.close();
	
	auto end=high_resolution_clock::now();
	auto duration = duration_cast<seconds>(end - start);
	cout << "Running time in seconds: " << duration.count() << endl << endl;
	
	return 0;
}

