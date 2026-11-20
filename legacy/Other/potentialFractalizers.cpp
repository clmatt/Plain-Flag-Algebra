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
	//omp_set_num_threads(1); //Use for debugging
	auto start=std::chrono::high_resolution_clock::now();
	
	//Change this	
	const int n = 8;
	std::vector<std::vector<int> > num(18, std::vector<int>(18,0));
	
	/*num[6][6] = 2; //n = 5
	num[7][6] = 2;
	num[7][7] = 4;
	num[8][6] = 1;
	num[8][7] = 3;
	num[8][8] = 8;
	num[9][6] = 1;
	num[9][7] = 3;
	num[9][8] = 8;
	num[9][9] = 16;*/
	
	/*num[7][7] = 2; //n = 6
	num[8][7] = 1;
	num[8][8] = 4;
	num[9][7] = 1;
	num[9][8] = 3;
	num[9][9] = 8;
	num[10][7] = 1;
	num[10][8] = 3;
	num[10][9] = 7;
	num[10][10] = 16;*/
	
	/*num[8][8] = 2; //n = 7
	
	num[9][8] = 1;
	num[9][9] = 4;
	
	num[10][8] = 1;
	num[10][9] = 3;
	num[10][10] = 8;
	
	num[11][8] = 1;
	num[11][9] = 2;
	num[11][10] = 6;
	num[11][11] = 16;
	
	num[12][8] = 1;
	num[12][9] = 2;
	num[12][10] = 6;
	num[12][11] = 14;
	num[12][12] = 32;
	
	num[13][8] = 1;
	num[13][9] = 2;
	num[13][10] = 5;
	num[13][11] = 13;
	num[13][12] = 30;
	num[13][13] = 64;
	
	num[14][8] = 1;
	num[14][9] = 2;
	num[14][10] = 5;
	num[14][11] = 13;
	num[14][12] = 30;
	num[14][13] = 64;
	num[14][14] = 128;*/
	
	num[9][9] = 2; //n = 8
	
	num[10][9] = 1;
	num[10][10] = 4;
	
	num[9][11] = 1;
	num[10][11] = 3;
	num[11][11] = 8;
	
	num[12][9] = 1;
	num[12][10] = 2;
	num[12][11] = 6;
	num[12][12] = 16;
	
	std::vector<Graph> allGraphs = generate(n,2);
	
	//Graph C7({{0,1,1},{1,2,1},{2,3,1},{3,4,1},{4,5,1},{5,6,1},{6,0,1}},7,2);
	//Graph C5prism({{0,1,1},{1,2,1},{2,3,1},{3,4,1},{4,0,1},{5,6,1},{6,7,1},{7,8,1},{8,9,1},{9,5,1},{0,5,1},{1,6,1},{2,7,1},{3,8,1},{4,9,1}},10,2);
	Graph test({{0,2,1},{0,3,1},{1,2,1},{1,3,1},{2,3,1},{0,4,1},{1,5,1},{2,6,1},{3,7,1}},8,2);
	
	
	
	allGraphs = {test};
	
	std::cout << "Starting number of graphs: " << allGraphs.size() << std::endl;
	
	//Remove all graphs with too many edges (could just use complement)
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		if(allGraphs[i].getCanonLabel() < allGraphs[i].complement().getCanonLabel()) {
			allGraphs.erase(allGraphs.begin() +i);
			--i;
		}
	}
	
	std::cout << "After checking that they have few enough edges: " << allGraphs.size() << std::endl;
	
	//Remove all graphs with twins
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		if(allGraphs[i].containsTwins()) {
			allGraphs.erase(allGraphs.begin() + i);
			--i;
		}		
	}
	
	std::cout << "After removing all graphs with twins: " << allGraphs.size() << std::endl;
	
	//Remove all graphs which aren't connected
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		if(!allGraphs[i].connected()) {
			allGraphs.erase(allGraphs.begin()+i);
			--i;
		}
	}
	
	std::cout << "After removing all graphs who are not connected: " << allGraphs.size() << std::endl;
	
	//Remove all graphs whose complements aren't connected
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		if(!allGraphs[i].complement().connected()) {
			allGraphs.erase(allGraphs.begin()+i);
			--i;
		}
	}
	
	std::cout << "After removing all graphs whose complements are not connected: " << allGraphs.size() << std::endl;
	
	std::pair<int,int> noTwins;
	noTwins.first = -1;
	noTwins.second = -1;
	
	for(int i = n+1; i <= 12; ++i) { //CHANGE THIS
		std::cout << std::endl;
		for(int j = 0; j < allGraphs.size(); ++j) {
			bool show = true;
			std::cout << "(" << i << ", " << j << ") out of (12, " << allGraphs.size() << ")" << std::endl;
			std::vector<Graph> S = {allGraphs[j]};
			
			for(int k = n+1; k <= i; ++k) {
				S = expandGraphs(S,{});
				
				//Needs enough subgraphs
				for(int l = 0; l < S.size(); ++l) {
					if(k != i) {
						if(numSubgraphs(allGraphs[j], S[l]) < num[i][k]) {
							S.erase(S.begin() +l);
							--l;
						}
					
						else if(S[l].findTwins() != noTwins) {
							S.erase(S.begin()+l);
							--l;
						}
					}	
						
					else {
						int temp = numSubgraphs(allGraphs[j], S[l]);	
						if((temp > num[i][k]) || ((temp == num[i][k]) && (S[l].findTwins() == noTwins) )) {
							S[l].printEdges();
							allGraphs.erase(allGraphs.begin() + j);
							l = S.size();						
							--j;
							show = false;
						}
					}
				}
				
				if(show && (i == 12) && (k == 12)) {
					allGraphs[j].printEdges();
				}
			}
			
			S.clear();
		}
		
		for(int j = 0; j < allGraphs.size(); ++j) {
			std::cout << "Graph " << j << std::endl;
			allGraphs[j].printEdges();
			std::cout << std::endl << std::endl;
		}
	}
	
	std::cout << "Total number of potential fractalizers = " << allGraphs.size() << std::endl << std::endl;
	
	for(int i = 0; i < allGraphs.size(); ++i) {
		allGraphs[i].printEdges();
	}

	auto end=std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	std::cout << "Running time in seconds: " << duration.count() << std::endl << std::endl;
	
	return 0;
}

