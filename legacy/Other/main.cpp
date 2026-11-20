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

int oddPath(Graph &test, const std::vector<Graph> &allPaths, const std::unordered_set<std::string> allPathsCanon) {
	int n = test.getN();
	
	if(test.getNumEdges() == 0) {
		return 0;
	}
	
	if(allPathsCanon.find(test.getCanonLabel()) != allPathsCanon.end()) {
		return 1;
	}
	
	
	std::vector<Graph> candidates = {test};
	int counter = 0;
	
	//TODO Odd vertices - Can assume start or end is on odd vertex until out of out vertices - how to do with linear ordering?
	while(true) {
		++counter;
		//std::cout << counter << std::endl;
		std::vector<Graph> newCandidates; 
		std::unordered_set<std::string> canon;
		
		for(int i = 0; i < candidates.size(); ++i) {
			
			int odd = -1;
			
			for(int j = 0; j < n; ++j) {
				if((candidates[i].getDegree(j) % 2) == 1) {
					odd = j;
					j = n;
				}
			}
		
			for(int j = 0; j < allPaths.size(); ++j) {
				int endPt1 = -1;
				int endPt2 = -1;
				
				//Could do using coefficients
				for(int k = 0; k < n; ++k) {
					if(allPaths[j].getDegree(k) == 1) {
						if(endPt1 == -1) {
							endPt1 = k;
						}
					
						else{
							endPt2 = k;
							k = n;
						}
					}
				}
			
			
				if((odd == -1) || (endPt1 == odd) || (endPt2 == odd)) {
					std::vector<Edge> newEdges;
					
					for(int k = 0; k < n-1; ++k) {
						for(int l = k+1; l < n; ++l) {
							if(((candidates[i].getEdgeColor(k,l) == 0) && (allPaths[j].getEdgeColor(k,l) == 1)) || ((candidates[i].getEdgeColor(k,l) == 1) && (allPaths[j].getEdgeColor(k,l) == 0))) {
								newEdges.push_back({k,l,1});
							}
						} 
					}
					
					Graph temp({newEdges,n,2});
					temp.setCoefficient(Frac(j+1,1)); 
					
					if(allPathsCanon.find(temp.getCanonLabel()) != allPathsCanon.end()) {
						return counter + 1; 
					}
					
					else if(canon.find(temp.getCanonLabel()) == canon.end()) {
						Graph temp({newEdges,n,2});
						newCandidates.push_back(temp);
						canon.insert(temp.getCanonLabel());
					}
				}
			}
		}
		
		candidates = newCandidates;
	}

}

int main() {
	//omp_set_num_threads(1); //Use for debugging
	auto start=std::chrono::high_resolution_clock::now();
	bool conj = true;
	
	std::ofstream myFile1;
	std::ofstream myFile2;
	myFile1.open("oddPathDictionary1.txt");
	myFile2.open("oddPathDictionary2.txt");
		

	for(int i = 2; i < 9; ++i) {
		std::vector<Graph> allGraphs = generate(i,2,{});
		
		//Generate all paths
		std::vector<Graph> allPaths1;
		std::unordered_set<std::string> allPathsCanon1;
		for(int j = 2; j <= i; ++j) {
			std::vector<int> subset;
		
			for(int k = 0; k < j; ++k) {
				subset.push_back(k);
			}  
		
			do {
				int* subsetOrdered = &subset[0];
				do {
					if(subsetOrdered[0] > subsetOrdered[j-1]) {
						std::vector<Edge> edges;
						for(int k = 0; k < j-1; ++k) {
							edges.push_back({subsetOrdered[k],subsetOrdered[k+1],1});
						
						}
						Graph temp({edges,i,2});
						temp.setCoefficient(Frac(subsetOrdered[0],1));
						allPaths1.push_back(temp);
						allPathsCanon1.insert(temp.getCanonLabel());
					}
				} while(std::next_permutation(subsetOrdered,subsetOrdered+j));
			
			} while(nextSubset(subset,i,j));
		}
		
		//Generate all paths
		/*std::vector<Graph> allPaths2;
		std::unordered_set<std::string> allPathsCanon2;
		for(int j = 2; j <= i+1; ++j) {
			std::vector<int> subset;
		
			for(int k = 0; k < j; ++k) {
				subset.push_back(k);
			}  
		
			do {
				int* subsetOrdered = &subset[0];
				do {
					if(subsetOrdered[0] > subsetOrdered[j-1]) {
						std::vector<Edge> edges;
						for(int k = 0; k < j-1; ++k) {
							edges.push_back({subsetOrdered[k],subsetOrdered[k+1],1});
						
						}
						Graph temp({edges,i+1,2});
						allPaths2.push_back(temp);
						allPathsCanon2.insert(temp.getCanonLabel());
					}
				} while(std::next_permutation(subsetOrdered,subsetOrdered+j));
			
			} while(nextSubset(subset,i+1,j));
		}
		
		#pragma omp parallel
		{
			#pragma omp for schedule(dynamic)
			for(int j = 0; j < (int)allGraphs.size(); ++j) {
				//allGraphs.setCoefficient(Frac(0,1));
				std::cout << "(" << i << ", " << j << ") out of (7, " << allGraphs.size() << ") \n";
				allGraphs[j].printEdges();
				int temp1 = oddPath(allGraphs[j],allPaths1,allPathsCanon1);
				allGraphs[j].addVertex();
				int temp2 = oddPath(allGraphs[j],allPaths2,allPathsCanon2);
			
				if(temp1 != temp2) {
					std::cout << "CONJECTURE IS WRONG!!!!" << std::endl;
					conj = false;
					//return 0;
				}
			}
		}*/
		
		for(int j = 0; j < (int)allGraphs.size(); ++j) {
			int temp = oddPath(allGraphs[j],allPaths1,allPathsCanon1);
			std::string tempStr = allGraphs[j].getCanonLabel();
			
			if(!tempStr.empty()) {
				tempStr.pop_back();
			}
			
			myFile2 << "Vertices = " << i << " Edges = ";
			
			for(int k = 0; k < i-1; ++k) {
				for(int l = k+1; l < i; ++l) {
					if(allGraphs[j].getEdgeColor(k,l) == 1) {
						myFile2 << "(" << k << ", " << l << ") ";
					}
				}
			}
			
			myFile2 << " Odd Path Cover Number = " << temp << "\n";
			
			std::cout << "(" << i << ", " << j << ") out of (7, " << allGraphs.size() << ") \n";
			myFile1 << tempStr << " " << temp << "\n";
			//std::cout << tempStr << " " << temp << "\n";
		}	
	}
	

	//if(!conj) {
		//std::cout << "CONJECTURE IS WRONG!!!!" << std::endl;
	//}	

	
	auto end=std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	std::cout << "Running time in seconds: " << duration.count() << std::endl << std::endl;
	
	return 0;
}

