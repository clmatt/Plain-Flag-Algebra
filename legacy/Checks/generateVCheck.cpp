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
#include <string>

using namespace std::chrono;


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
	vector<Graph> zeros;
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

	auto start1=high_resolution_clock::now();
	
	vector < vector <Graph> > test1 = NEWgenerateV(7,3,zeros);
	
	auto start2=high_resolution_clock::now();
	auto duration1 = duration_cast<milliseconds>(start2 - start1);
	
	
	auto start3=high_resolution_clock::now();
	
	vector < vector <Graph> > test2 = generateV(7,3,zeros);

	
	
	
	auto start4=high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(start4 - start3);
	cout << "Running time of NEWgenerateV: " << duration1.count() << endl << endl;
	cout << "Running time in generateV: " << duration.count() << endl << endl;
	
	for(int i = 0; i < test1.size(); ++i) {
		cout << test1[i].size() << " ";
	}
	cout << endl << endl << endl;
	for(int i = 0; i < test2.size(); ++i) {
		cout << test2[i].size() << " ";
	}
	cout << endl;
	
	return 0;
}

