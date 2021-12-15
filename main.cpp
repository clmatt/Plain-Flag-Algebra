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
	auto start=high_resolution_clock::now();
	
	
	for(int j = 2; j < 5; ++j) {
		for(int k = j; k < 10; ++k) {
			for(int i = 0; i < 1000; ++i) {
				cout << i << " " << j << " " << k << endl;
				Graph H = randomGraph(j, 2, {0.5,0.5}, 0);
				Graph G = randomGraph(k, 2, {0.5,0.5}, 0);
				
				vector< vector<int> > temp;
				returnSubgraphsNoFlags(H,G,temp);
				
				if(numSubgraphs(H,G) != (int)temp.size()) {
					cout << numSubgraphs(H,G) << " " << (int)temp.size() << endl;
					for(int l = 0; l < (int)temp.size(); ++l) {
						for(int n = 0; n < (int)temp[l].size(); ++n) {
							cout << temp[l][n] << " ";
						}
						cout << endl;
					}
					cout << endl;
					H.printEdges();
					G.printEdges();
					return 0;
				}
			}
		}
	}
	
	/*ifstream myFile;
	myFile.open("r44_14.g6");
	string line;
	char* canonLabel;
	const int subgraphSize = 4;
	
	vector< vector< Graph> > allGraphs;
	
	allGraphs.push_back({});
	allGraphs.push_back({});
	for(int i = 2; i <= subgraphSize; ++i) {
		allGraphs.push_back(generate(i, 2, {}, {}));
	}
	
	vector<Graph> fromFile;
	
	cout << endl << "Reading in file." << endl << endl;
	int counter = 0;
	
	if (myFile.is_open()) {
		while (getline(myFile,line)) {
			if(counter %10000 == 0) {
				cout << "Iteration number " << counter << endl;
			}
			++counter;
			canonLabel = &line[0];
			sparsegraph sg;
    		SG_INIT(sg);
    		int num_loops;
			
			stringtosparsegraph(&canonLabel[0], &sg, &num_loops);
    		Graph G = convertFromNauty(sg);
    		fromFile.push_back(G);
		}
   	myFile.close();
  	}
  	cout << endl;
  	
  	ofstream outputFile;
  	outputFile.open("r44_14_subgraphBounds.txt");
  	
  	for(int i = 2; i <= subgraphSize; ++i) {
  		for(int j = 0; j < (int)allGraphs[i].size(); ++j) {
  			cout << "(" << i << ", " << j << ") out of (" << subgraphSize << ", " << (int)allGraphs[i].size() << ")" << endl;
  			int temp = numSubgraphs(allGraphs[i][j],fromFile[0]);
  			int minNum = temp;
  			int maxNum = temp;
  			int counter = 1;
  			
  			for(int k = 1; k < (int)fromFile.size(); ++k) {
  				if(counter %10000 == 0) {
					cout << "Iteration number " << counter << endl;
				}
				++counter;
				temp = numSubgraphs(allGraphs[i][j],fromFile[k]);
  				minNum = min(minNum,temp);
  				maxNum = max(maxNum,temp);
  			}
  			cout << endl;
  			outputFile << "Adjacency Matrix:" << endl;
  			allGraphs[i][j].printAdjMatToFile(outputFile);
  			outputFile << "Minimum number: " << minNum << endl;
  			outputFile << "Maximum number: " << maxNum << endl << endl;
  		}
  	}*/
	
	
	auto end=high_resolution_clock::now();
	auto duration = duration_cast<seconds>(end - start);
	cout << "Running time in seconds: " << duration.count() << endl << endl;
	
	return 0;
}

