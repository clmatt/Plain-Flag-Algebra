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
	auto time1=high_resolution_clock::now();
	
	for(int i = 1; i <= 4; ++i) {
		for(int j = i; j <= 4; ++j) {
			cout << i << " " << j << endl;
			for(int k = 0; k < 100; ++k) {
				for(int f = 0; f <= 0; ++f) {
					Graph H = uniformRandomGraph(i,2,f);
					Graph G = uniformRandomGraph(j,2,f);
					
					NEWmultiply(H,G,{});
				}
			}
		}
	}
	
	auto time2=high_resolution_clock::now();
	auto duration1 = duration_cast<milliseconds>(time2-time1);
	cout << "Running time of new multiply: " << duration1.count() << endl << endl;
	
	auto time3=high_resolution_clock::now();
	
	for(int i = 1; i <= 4; ++i) {
		for(int j = i; j <= 4; ++j) {
			cout << i << " " << j << endl;
			for(int k = 0; k < 100; ++k) {
				Graph H = uniformRandomGraph(i,2,0);
				Graph G = uniformRandomGraph(j,2,0);
				multiply(H,G,{});
			}
		}
	}
	
	
	auto time4=high_resolution_clock::now();
	auto duration2 = duration_cast<milliseconds>(time4 - time3);
	cout << "Running time of old multiply: " << duration2.count() << endl << endl;
	
	return 0;
}

