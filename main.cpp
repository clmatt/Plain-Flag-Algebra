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
	auto start1=high_resolution_clock::now();
	
	Graph G({{0,2,1},{1,3,1}},4,2);
	vector<Graph> output = expandGraphs({G},{});
	G.printOrbits();
	
	for(int i = 0; i < output.size(); ++i) {
		output[i].printEdges();
	}
	
	return -1;
	
	vector<Graph> test1 = NEWgenerate(5,2,{});
	
	auto start2=high_resolution_clock::now();
	auto duration1 = duration_cast<milliseconds>(start2 - start1);
	
	
	auto start3=high_resolution_clock::now();
	
	vector<Graph> test2 = generate(5,2,{});

	
	
	
	auto start4=high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(start4 - start3);
	cout << "Running time of NEWgenerate: " << duration1.count() << endl << endl;
	cout << "Running time in generate: " << duration.count() << endl << endl;
	
	bool test;
	for(int i = 0; i < test2.size(); ++i) {
		test = false;
		for(int j = 0; j < test1.size(); ++j) {
			if(isomorphic(test2[i],test1[j])) {
				test = true;
			} 
		}
		
		if(!test) {
			test2[i].printEdges();
		}
	}
	
	cout << test1.size() << endl;
	cout << test2.size() << endl;
	
	return 0;
}

