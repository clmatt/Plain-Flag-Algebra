//Nauty
//#define MAXN 1000    /* Define this before including nauty.h */
extern "C" {
	#include "nauty.h"   
	#include "naututil.h"
	#include "gtools.h"
}

const std::vector<long long int> numberOfGraphs = {1,1,2,4,11, 34, 156, 1044, 12346, 274668, 12005168, 1018997864, 165091172592}; //Up to 12 vertices

//---------------------
//-----Edge Struct-----
//---------------------

//This has some inherient ordering within it
//In general, we will assume a < b, unless we are in a directed graph
//If an edge isn't specified in Graph, it will assumed to be color 0
struct Edge {
	int a;
	int b;
	int color;
};

//---------------------
//---------------------
//-----Graph Class-----
//---------------------
//---------------------

//Vertices and colors are 0 indexed
//No loops or multiple edges
//Automatically gives canonLabel for every graph using Nauty 
class Graph {

	private:
		
		int n = -1; //Number of vertices
		int numColors = -1; //Number of colors
		Frac coefficient = Frac(1,1); //Helpful for averaging, can assume 1 if not using
		int sizeOfFlag = 0;
		std::vector<std::vector<int> > adjMat; //Adjacency Matrix
		std::vector<int> flag; //The flag can be determined by a subset of edges
		std::string canonLabel;
		std::vector<std::vector<int> > orbits; //orbits in automorphsim group (from nauty) orbits[i] is the ith orbit of vertices, if vertex is a flag it is in its own orbit (not iff)
		int numOrbits = 0;

	public:
		
		//-------------------------------
		//-----Constructor for Graph-----
		//-------------------------------
		
		//Need to specify NUMCOLORS, because it's possible that it differs
		//from actual number of colors
		Graph(const std::vector<Edge> &edgeList, const int N, const int NUMCOLORS = 2, const std::vector<int> &FLAG = {}) {
			numColors = NUMCOLORS;
			n = N;
			
			if(numColors < 0) {
				std::cout << "In constructor need a non-negative number of colors" << std::endl << std::endl;
				throw std::exception();
			}
			
			if(n < 0) {
				std::cout << "In constructor need a non-negative number of vertices" << std::endl << std::endl;
				throw std::exception();
			}
			
			//Intialize Adjacency Matrix
			//Colors of edges assumed to be 0 unless otherwise specified
			adjMat.resize(n);
			for(int i = 0; i < n; ++i) {
				adjMat[i].resize(n);
			}
			
			for(int i = 0; i < n; ++i) {
				for(int j = 0; j < n; ++j) {
					adjMat[i][j] = 0;
				}
			}
			
			//Insert edges for edgeList
			for (auto &edge: edgeList) {
				if(edge.color >= numColors) {
					std::cout << "In constructor too few colors specified, or range is wrong, note that non-edges count as a color, and colors are zero indexed"  << std::endl << std::endl;
					throw std::exception();
				}
				
				if(edge.a >= n) {
					std::cout << "In constructor an edge has a vertex outside of the range." << std::endl << std::endl;
					throw std::exception();
				}
				
				if(edge.b >= n) {
					std::cout << "In constructor an edge has a vertex outside of the range." << std::endl << std::endl;
					throw std::exception();
				}
			
            adjMat[edge.a][edge.b] = edge.color;
            adjMat[edge.b][edge.a] = edge.color;
  			}
  			
  			sizeOfFlag = 0;
		
			for(auto &i : FLAG) {
				++sizeOfFlag;
				
				if((i >= n) || (i < 0)) {
					std::cout << "In setFlag in constructor an element is outside of the vertex range." << std::endl << std::endl;
					throw std::exception();
				}	
			}
			
			for(int i = 0; i < (int)FLAG.size(); ++i) {
				for(int j = i+1; j < (int)FLAG.size(); ++j) {
					if(FLAG[i] == FLAG[j]) {
						std::cout << "In setFlag in constructor all elements of the std::vector must be different." << std::endl << std::endl;
						throw std::exception();
					}
				}
			}
		
			flag = FLAG;
  			
  			Frac frac(1,1);
  			coefficient = frac;
  			
  			if(n == 0) {
  				canonLabel = "\n";
  			}
  			
  			else {
  				canonRelabel();
  			}
		}

		
		//---------------
		//-----Clone-----
		//---------------
		
		Graph clone() const{
			std::vector<Edge> edges;
		
			for(int i = 0; i < n-1; ++i) {
				for(int j = i+1; j < n; ++j) {
					edges.push_back({i,j,adjMat[i][j]});
				}
			}
			
			Graph G(edges,n,numColors,flag);
			
			G.setCoefficient(coefficient);
			
			return G;
		}
		
		
		//-----------------------------
		//-----Restrict Vertex Set-----
		//-----------------------------
		
		//Can be a reorder if sigma is a bijection
		//Make sigma a size n vertex such that sigma[i] = j means vertex i in new graph is j
		//E.g. range is 0,1,..,k-1 and -1
		//Make everything else -1
		//Doesn't actually change the graph
		Graph restriction(const std::vector<int> &sigma) const{
			if(int(sigma.size()) != n) {
				std::cout << "In restriction, sigma size wrong." << std::endl;
				std::cout << "Sigma size = " << sigma.size() << std::endl;
				std::cout << "n = " << n << std::endl << std::endl;
				throw std::exception();
			}		
					
			//Check that it's a subset
			int k = 0;
					
			for(int i = 0; i < n; ++i) {
				if(sigma[i] != -1) {
				++k;
					for(int j = i+1; j < n; ++j) {
						if(sigma[i] == sigma[j]) {
							std::cout << "In restriction there is a duplicate in the input." << std::endl;
							std::cout << "Indices : " << i << ", " << j << std::endl;
							std::cout << "Value: " << sigma[i] << std::endl << std::endl;
							throw std::exception();
						}
					}
				}
			}
			
			//Check if sigma has elements outside of range
			for(int i = 0; i < n; ++i) {
				if((sigma[i] >= k) || (sigma[i] < -1)) {
					std::cout << "In restriction an element of sigma falls outside of range of vertices." << std::endl;
					std::cout << "Index: " << i << std::endl;
					std::cout << "Value: " << sigma[i] << std::endl << std::endl;
					throw std::exception();
				}
			}
			
			//Edges
			std::vector<Edge> newEdges;
			for(int i = 0; i < n-1; ++i) {
				for(int j = i+1; j < n; ++j) {
					if((adjMat[i][j] != 0) && (sigma[i] != -1) && (sigma[j] != -1)) {
						newEdges.push_back({sigma[i],sigma[j],adjMat[i][j]});
					}
				}
			}
			
			//Flag
			std::vector<int> newFlag;
			for(int j = 0; j < sizeOfFlag; ++j) {
				if(sigma[flag[j]] != -1) {
					newFlag.push_back(sigma[flag[j]]);
				}
			}
			
			Graph G(newEdges, k, numColors,newFlag);

			return G;
		}
		
		
		//--------------------------------
		//-----Print Adjacency Matrix-----
		//--------------------------------
		
		void printAdjMat() const{
			for(int i = 0; i < n; ++i) {
				for(int j = 0; j < n; ++j) {
					std::cout << adjMat[i][j] << " ";
				}
				std::cout << std::endl;
			}
		}
		
		//----------------------------------------
		//-----Print Adjacency Matrix to File-----
		//----------------------------------------
		
		void printAdjMatToFile(std::ofstream& myFile) const{
			for(int i = 0; i < n; ++i) {
				for(int j = 0; j < n; ++j) {
					myFile << adjMat[i][j] << " ";
				}
				myFile << std::endl;
			}
		}	
	
	
		//-----------------------------
		//-----Print Flag Vertices-----
		//-----------------------------
		
		void printFlagVertices() const{
			if (sizeOfFlag == 0) {
				std::cout << "There are no flag vertices." << std::endl;
			}	
			
			else {
				for(auto i : flag) {
					std::cout << i << " ";
				}
				std::cout << std::endl;
			} 
		}	
		
		
		//----------------------
		//-----Print Orbits-----
		//----------------------
		
		void printOrbits() const{
			for(int i = 0; i < int(orbits.size()); ++i) {
				std::cout << "Orbit " << i+1 << " contains vertices ";
				for(int j = 0; j < int(orbits[i].size()); ++j) {
					std::cout << orbits[i][j] << " "; 
				}
				std::cout << std::endl;
			}
		}
		
		
		//-------------------
		//-----Get Orbit-----
		//-------------------
		
		int getOrbit(int i, int j) const{
			return orbits[i][j];
		}
		
		
		//------------------------
		//-----Get Orbit Size-----
		//------------------------
		int getOrbitSize(int i) const{
			return (int)orbits[i].size();
		}
		
		
		//-----------------------
		//-----Get numOrbits-----
		//-----------------------
		
		int getNumOrbits() const {
			return numOrbits;
		}
		
		
		//---------------
		//-----Get n-----
		//---------------
		
		int getN() const{
			return n;
		} 
		
		//-----------------------
		//-----Get numColors-----
		//-----------------------
		
		int getNumColors() const{
			return numColors;
		}
		
		
		//-----------------------------
		//-----Get number of Edges-----
		//-----------------------------
		
		int getNumEdges() const{
			int numEdges = 0;
		
			for(int i = 0; i < n-1; ++i) {
				for(int j = i+1; j < n; ++j) {
					if(adjMat[i][j] != 0) {
						++numEdges;
					}
				}
			}
			
			return numEdges;
		}
		
		
		//-----------------------
		//-----Set numColors-----
		//-----------------------
		
		void setNumColors(int c) {
			for(int i = 0; i < n-1; ++i) {
				for(int j = 0; j < n; ++j) {
					if (c <= adjMat[i][j]) {
						std::cout << "In setNumColors, input too small." << std::endl << std::endl;
						throw std::exception();
					}
				}
			}
			
			numColors = c;
			canonRelabel();
		}
		
		//----------------------------
		//-----Get Edge Color------
		//----------------------------
		
		int getEdgeColor(int i, int j) const{
			if((i < 0) || (i >= n)) {
				std::cout << "First vertex out of range in getEdgeColor." << std::endl << std::endl;
				throw std::exception();
			}
			
			else if((j < 0) || (j >= n)) {
				std::cout << "Second vertex out of range in getEdgeColor." << std::endl << std::endl;
				throw std::exception();
			}
		
			return adjMat[i][j];
		}		
		
		
		//-------------------------
		//-----Get Canon Label-----
		//-------------------------
		
		std::string getCanonLabel() const{ 
			return canonLabel;
		}
	
		
		//---------------------
		//-----Flag Vertex-----
		//---------------------
		
		bool flagVertex(int i) const{
			for(auto j : flag) {
				if(j == i) {
					return true;
				}
			}
			
			return false;
		}
		
		//--------------------
		//-----Get Degree-----
		//--------------------
		
		//This one does degree in any color (except 0)
		int getDegree(int i) const{
			int deg = 0;
		
			for(int j = 0; j < n; ++j) {
				if((adjMat[i][j] != 0) && (i != j)) {
					++deg;
				}
			}
			
			return deg;
		}
		
		//c is the color of the edge
		int getDegree(int i, int c) const{
			int deg = 0;
		
			for(int j = 0; j < n; ++j) {
				if((adjMat[i][j] == c) && (i != j)) {
					++deg;
				}
			}
			
			return deg;
		}
		
	
		//---------------------
		//-----Print Edges-----
		//---------------------
		
		void printEdges() const{
			for(int i = 0; i < numColors; ++i) {
				bool comma = false;
				std::cout << "The edge list for color " << i << " is: {";
				for(int j = 0; j < n-1; ++j) {
					for(int k = j+1; k < n; ++k) {
						if(adjMat[j][k] == i) {
							if(comma) {
								std::cout << ",";
							}
							comma = true;
							std::cout << "{" << j << "," << k << "}";
						}
					}
				}
				std::cout << "}" << std::endl;
			}	
			std::cout << std::endl;
		}
		
		
		//---------------------
		//-----Change Edge-----
		//---------------------
		
		void changeEdge(Edge edge) {
			if((edge.a >= n) || (edge.b >= n)) {
				std::cout << "In changeEdge an edge has a vertex outside of the range." << std::endl << std::endl;
				throw std::exception();
			}
			
			if(edge.color >= numColors) { 
				std::cout << "In changeEdge the color is outside of the range." << std::endl << std::endl;
				throw std::exception();
			}
		
			adjMat[edge.a][edge.b] = edge.color;
			adjMat[edge.b][edge.a] = edge.color;
			
			canonRelabel();
		}
		

		//--------------------
		//-----Add Vertex-----
		//--------------------
		
		//Always adds a vertex with highest label
		//Default color is 0
		void addVertex() {
			std::vector<int> zeros;
			
			for(int i = 0; i < n; ++i) {
				adjMat[i].push_back(0);
				zeros.push_back(0);
			} 
			
			n = n + 1;
			
			zeros.push_back(0);
			adjMat.push_back(zeros);
			
			canonRelabel();
		}
		
		//-----------------------
		//-----Remove Vertex-----
		//-----------------------
		
		//Relabels all vertices to be from 0 to n-1
		void removeVertex(int i) {
			if(i >= n) {
				std::cout << "In removeVertex the input is outside of the range." << std::endl << std::endl;
				throw std::exception();
			}
			
			//Remove from flag
			int index = 0;
			int vertexPos = -1;
			
			
			for(auto j : flag) {
				if(j == i) {
					vertexPos = index;
					--sizeOfFlag;
				}

				else if(j > i) {
					flag[index] = flag[index]-1; 
				}
				
				++index;
			}
			
			if(vertexPos > -1) {
				flag.erase(flag.begin()+vertexPos);
			}
			
			//Remove from graph
			n = n-1;
			adjMat.erase(adjMat.begin()+i);
			
			for(int j = 0; j < n; ++j) {
				adjMat[j].erase(adjMat[j].begin()+i);
			}
			
			canonRelabel();
		}
		
		
		//-------------------------
		//-----Set Coefficient-----
		//-------------------------
		
		void setCoefficient(Frac val) {
			coefficient = val;
		}
		
		void setCoefficient(int val) {
			coefficient = convertFromInt(val);
		}
		
		
		//-------------------------
		//-----Get Coefficient-----
		//-------------------------
		
		Frac getCoefficient() const{
			return coefficient;
		}
		
		
		//-------------------
		//-----Set Flag------
		//-------------------
		
		//Note that flags have order, so order of &FLAG does matter
		void setFlag(std::vector<int> const &FLAG) {
			sizeOfFlag = 0;
		
			for(auto &i : FLAG) {
				++sizeOfFlag;
				
				if((i >= n) || (i < 0)) {
					std::cout << "In setFlag an element is outside of the vertex range." << std::endl << std::endl;
					throw std::exception();
				}	
			}
			
			for(int i = 0; i < (int)FLAG.size(); ++i) {
				for(int j = i+1; j < (int)FLAG.size(); ++j) {
					if(FLAG[i] == FLAG[j]) {
						std::cout << "In setFlag all elements of the std::vector must be different." << std::endl << std::endl;
						throw std::exception();
					}
				}
			}
		
			flag = FLAG;
			
			canonRelabel();
		}
		
		//-------------------------
		//-----Get Flag Vertex-----
		//-------------------------
		
		int getFlagVertex(int i) const{
			return flag[i];
		}
		
		
		//-----------------
		//-----Is flag-----
		//-----------------
		
		bool isFlag(int i) const{
			if((i < 0) || (i >= n)) {
				std::cout << "i is outside of range in isFlag." << std::endl;
				std::cout << "i = " << i << std::endl << std::endl;
				throw std::exception();
			}
			
			for(int j = 0; j < sizeOfFlag; ++j) {
				if(flag[j] == i) {
					return true;
				}
			}
			
			return false;
		}
		
		
		//------------------
		//-----Get Flag-----
		//------------------
		
		//The graph given relabels the vertices, preserving the order
		//Also maintains that it is still a flag
		Graph getFlag() const{
			
			//Create Edges
			int i = 0;
			int j = 0; 
			std::vector<Edge> edges;
			for(int k : flag) {
				for(int l : flag) {
					edges.push_back({i, j, adjMat[k][l]});
					++j;
				}
				j = 0;
				++i;
			}
			
			std::vector<int> flag;
			for(int i = 0; i < sizeOfFlag; ++i) {
				flag.push_back(i);
			}
			
			Graph G(edges, sizeOfFlag, numColors,flag);
			
			return G;
			
		}
		
		//-----------------------------
		//-----Print Flag Vertices-----
		//-----------------------------
		
		void printFlag() const{
			std::cout << "The vertices of the flag are (in order): ";
			for(auto j : flag) {
				std::cout << j << " ";
			}
			std::cout << std::endl << std::endl;
		}
		
		
		//--------------------------
		//-----Get Size of Flag-----
		//--------------------------
		
		int getSizeOfFlag() const{
			return sizeOfFlag;
		}
		
		
		//----------------------------------
		//-----Convert to layered graph-----
		//----------------------------------

		//Used in covertToNauty
		//Nauty only accepts vertex colored graphs, so this graph outputs a layered graph
		//Also edits the input in c to give the vertex coloring
		//Note that this removes flags but the info is encoded in c
		//This just returns adjMat as otherwise we get infinite loop when creating graph
		void convertToLayer(std::vector<std::vector<int> > &newAdjMat, std::vector<int> &c, int &np) {		
			int k = floor(log2(getNumColors()-1)+1); //number of layers
			np = n * k;
			
			newAdjMat.resize(np);
			for(int i = 0; i < np; ++i) {
				newAdjMat[i].resize(np);
			}
			
			for(int i = 0; i < np; ++i) {
				for(int j = 0; j < np; ++j) {
					newAdjMat[i][j] = 0;
				}
			}

			//Add paths between layers
			for(int i = 0; i < getN(); ++i) {
				for(int j = 0; j < k-1; ++j) {
					newAdjMat[i+j*getN()][i+(j+1)*getN()] = 1;
					newAdjMat[i+(j+1)*getN()][i+j*getN()] = 1;
				}
			}

			//Add edges between each layer
			for(int i = 0; i < getN(); ++i) {
				for(int j = i+1; j < getN(); ++j) {
					if(getEdgeColor(i,j) != 0) {
						for(int a = 0; a < k; ++a) {
							if((getEdgeColor(i,j) >> a) & 1) {
								newAdjMat[i+a*getN()][j+a*getN()] = 1;
								newAdjMat[j+a*getN()][i+a*getN()] = 1;
							}
						}
					}
				}
			}
			
			//Make c correct
			int l = getSizeOfFlag();
			int flagCounter = 0;
			c.resize(np);
			
			for(int a = 0; a < k; ++a) {
				for(int i = 0; i < getN(); ++i) {
					if(!flagVertex(i)) {
						c[i+a*getN()] = l;
					}
				}	
				
				for(int i = 0; i < getSizeOfFlag(); ++i) {
					c[getFlagVertex(i) + a*getN()] = flagCounter;
					++flagCounter;
				}
				
				flagCounter = l+1;
				l = l + getSizeOfFlag() + 1;
			}
		}
		
		//-------------------------------
		//-----Canoncial Relabelling-----
		//-------------------------------
		
		//This doesn't actually change the graph but gives the unique std::string from Nauty 
		void canonRelabel() {

			DYNALLSTAT(int,lab,lab_sz);
			DYNALLSTAT(int,ptn,ptn_sz);
			DYNALLSTAT(int,orbitsTEMP,orbitsTEMP_sz);
			static DEFAULTOPTIONS_SPARSEGRAPH(options);
			statsblk stats;
			
			options.getcanon = TRUE;
			options.digraph = FALSE;
			options.defaultptn = FALSE;

			sparsegraph* sg = (sparsegraph*)malloc(sizeof(sparsegraph));
			SG_INIT((*sg));	
			int np,mp;
			
			std::vector<int> c;
			std::vector<std::vector<int> > newAdjMat;
			convertToLayer(newAdjMat,c,np); //Nauty doesn't allow for edge colored graphs
			
			mp = SETWORDSNEEDED(np);
			nauty_check(WORDSIZE,mp,np,NAUTYVERSIONID);
			
			DYNALLOC1(int,lab,lab_sz,np,"malloc");
			DYNALLOC1(int,ptn,ptn_sz,np,"malloc");
			DYNALLOC1(int,orbitsTEMP,orbitsTEMP_sz,np,"malloc");
			
			int maxC = -1;
			for(int i = 0; i < np; ++i) {
				if(maxC < c[i]) {
					maxC = c[i];
				}
			}
			
			int temp = 0;
			for(int i = 0; i <= maxC; ++i) {
				for(int j = 0; j < np; ++j) {
					if(c[j] == i) {
						*(lab+temp) = j;
						*(ptn+temp) = 1;
						++temp;
					}
				}
				*(ptn+temp-1) = 0;
			}
			
			int numEdges = 0;
			for(int i = 0; i < np; ++i) {
				for(int j = 0; j < np; ++j) {
					if(newAdjMat[i][j] != 0) {
						++numEdges;
					}
				}
			}
			
			SG_ALLOC((*sg),np,numEdges,"malloc");
			
			sg->nv = np;	
			sg->nde = numEdges;
			
			temp = 0;
			int degree;
			for(int i = 0; i < np; ++i) {
				sg->v[i] = temp;
				degree = 0; 
				for(int j = 0; j < np; ++j) {
					if (newAdjMat[i][j] != 0) {
						++degree; 
						++temp;
					}
				}	
				sg->d[i] = degree;
			}
			
			for(int i = 0; i < np; ++i) {
				temp = 0;
				for(int j = 0; j < np; ++j) {
					if (newAdjMat[i][j] != 0) {
						sg->e[(sg->v[i]) + temp] = j;
						++temp;
					}
				}
			}
			
			sparsegraph* canong = (sparsegraph*)malloc(sizeof(sparsegraph));
			SG_INIT((*canong));
			
			sparsenauty(sg,lab,ptn,orbitsTEMP,&options,&stats,canong);
			sortlists_sg(canong);
			
			numOrbits = 0;
			orbits.clear();
			
			for(int i = 0; i < n; ++i) {
				if((*(orbitsTEMP+i)) == i) {
					orbits.push_back({i});
					
					for(int j = i+1; j < n; ++j) {
						if((*(orbitsTEMP+j)) == i) {
							orbits[numOrbits].push_back(j);
						}
					}
					
					++numOrbits;
				} 
			}
			
			canonLabel = std::string(sgtos6(canong));
			

			DYNFREE(lab,lab_sz);
			DYNFREE(ptn,ptn_sz);
			DYNFREE(orbitsTEMP,orbitsTEMP_sz);
			
			SG_FREE((*sg));
			SG_FREE((*canong));
			free(sg);
			free(canong);
		}
		
		
		//---------------------
		//-----Remove Flag-----
		//---------------------
		
		//Used in Average
		//Sigma are vertices kept
		void removeFlag() {
			sizeOfFlag = 0;
			flag = {};
			
			canonRelabel();
		}
		
		
		//---------------------
		//-----Average All-----
		//---------------------	
		
		void averageAll() {
			bool isomorphic(const Graph&, const Graph&);
			void returnSubgraphs(const Graph&, const Graph&, std::vector<std::vector<int> >&);
			int k = sizeOfFlag;
			
			Graph G = clone();
			Graph H = clone();
			Graph GFlag = G.getFlag();
			GFlag.removeFlag();
			
			H.removeFlag();
			removeFlag();
			
			//Determine change in coefficient
			long long int den = 0;
			long long int num = 0;
			
			if(k == 1) {
				for(int i = 0; i < H.getNumOrbits(); ++i) {
					for(int j = 0; j < H.getOrbitSize(i); ++j) {
						if(G.getFlagVertex(0) == H.getOrbit(i,j)) {
							num = H.getOrbitSize(i);
							j = H.getOrbitSize(i);
						}
					}
				}
			}
			
			else {
				//Add flag back in in anyway possible
				std::vector<int> subset;
				subset.resize(k);
				
				std::vector<std::vector<int> > subgraphs;
				returnSubgraphs(GFlag,H,subgraphs);
				
				for(auto X: subgraphs) {
					subset.clear();
					for(int i = 0; i < (int)X.size(); ++i) {
						if(X[i] != -1) {
							subset.push_back(i);
						}
					}
					
					//Iterate over all permutations of the subset
					do {
						H.setFlag(subset);	
														
						if(isomorphic(G,H)) {
							++num;
						}
						
					} while(next_permutation(subset.begin(),subset.end()));
					
					sort(subset.begin(),subset.end());
				}
			}
			
			den = choose(n,k) * factorial(k);
			
			coefficient = coefficient*Frac(num,den);
			
			canonRelabel();
		}
		
		//-----------------
		//-----Is twin-----
		//-----------------

		bool isTwin(const int u, const int v) const{
			if(sizeOfFlag != 0) {
				std::cout << "In isTwin, size of flag must be 0." << std::endl << std::endl;
				throw std::exception();	
			}
			
			for(int i = 0; i < n; ++i) {
				if((i != u) && (i != v)) {
					if(adjMat[i][u] != adjMat[i][v]) {
						return false;
					}
				}
			}
			
			return true;
		}

		//------------------------
		//-----Contains Twins-----
		//------------------------

		bool containsTwins() const{
			if(sizeOfFlag != 0) {
				std::cout << "In isTwin, size of flag must be 0." << std::endl << std::endl;
				throw std::exception();	
			}
			
			for(int i = 0; i < n-1; ++i) {
				for(int j = i+1; j < n; ++j) {
					if(isTwin(i,j)) {
						return true;
					}
				}
			}
			
			return false;
		}


		//-----------------------------
		//-----Find a set of twins-----
		//-----------------------------

		//Return (-1,-1) if no twins
		std::pair<int,int> findTwins() const{
			if(sizeOfFlag != 0) {
				std::cout << "In isTwin, size of flag must be 0." << std::endl << std::endl;
				throw std::exception();	
			}
			
			std::pair<int,int> output;
			
			for(int i = 0; i < n-1; ++i) {
				for(int j = i+1; j < n; ++j) {
					if(isTwin(i,j)) {
						output.first = i;
						output.second = j;
						return output;
					}
				}
			}
			
			output.first = -1;
			output.second = -1;
			
			return output;
		}
		
		//--------------------
		//-----Complement-----
		//--------------------

		Graph complement() const{
			if(numColors != 2) {
				std::cout << "Can only take a graph complement if there are two colors." << std::endl << std::endl;
				throw std::exception();
			}
			
			if(sizeOfFlag != 0) {
				std::cout << "In complement, size of flag must be 0." << std::endl << std::endl;
				throw std::exception();	
			}
			
			std::vector<Edge> newEdges;
			
			for(int i = 0; i < n-1; ++i) {
				for(int j = i+1; j < n; ++j) {
					if(adjMat[i][j] == 0) {
						newEdges.push_back({i,j,1});
					}
				}
			}
			
			Graph output(newEdges,n,2);
			
			return output;
		}


		//-------------------
		//-----Connected-----
		//-------------------

		//Assumes for more than two colors it is colorblind for >0
		//Not the fastest, but good enough for now
		bool connected() const{
			if(sizeOfFlag != 0) {
				std::cout << "In isTwin, size of flag must be 0." << std::endl << std::endl;
				throw std::exception();	
			}
			
			std::vector<bool> component(n, false);
			std::vector<bool> newComponent(n, false);
			
			newComponent[0] = true;
			bool cont = true;
			
			while(component != newComponent) {
				component = newComponent;
				
				for(int i = 0; i < n; ++i) {
					if(component[i]) {
						for(int j = 0; j < n; ++j) {
							if(adjMat[i][j] != 0) {
								newComponent[j] = true;
							}
						}
					}  
				}	
			}
			
			for(int i = 0; i < n; ++i) {
				if(!component[i]) {
					return false;
				}
			}
			
			return true;
		}
};

//------------------------------------
//------------------------------------
//-----Functions Involving Graphs-----
//------------------------------------
//------------------------------------


//------------------------
//-----Graph Equality-----
//------------------------

bool operator==(const Graph &G, const Graph &H) {
	if(G.getCanonLabel() != H.getCanonLabel()) {
		return false;
	}
	
	if(G.getNumColors() != H.getNumColors()) {
		return false;
	}
	
	return true;
}


//----------------------------------------------
//-----Convert from Nauty to my Graph Class-----
//----------------------------------------------

//Hopefully won't need very much b/c we can just get canonical labeling from Nauty then get rid of that graph
//Careful with this- when converting to Nauty we use the layerGraph so this isn't a direct inverse, it drops the flags and colors and lots of other things, mostly just used for debugging

Graph convertFromNauty(const sparsegraph &sg) {
	std::vector<Edge> edges;
	
	for(int i = 0; i < sg.nv; ++i) {
		for(int j = 0; j < sg.d[i]; ++j) {
			edges.push_back({i,sg.e[sg.v[i]+j],1});
		}
	}
	
	Graph G = Graph(edges,sg.nv,2);	
	
	return G;
}	


//-------------------------------
//-----Same Adjacency Matrix-----
//-------------------------------

//Technically supplemented by Nauty
bool sameAdjMat(const Graph &G, const Graph &H) {
	if(G.getN() != H.getN()) {
		return false;
	}
	
	int n = G.getN();
	
	for(int i = 0; i < n-1; ++i) {
		for(int j = i+1; j < n; ++j) {
			if(G.getEdgeColor(i,j) != H.getEdgeColor(i,j)) {
				return false;
			}
		}
	}
	
	return true;
}

//---------------------
//-----Isomorphic------
//---------------------

//Uses Nauty
bool isomorphic(const Graph &G, const Graph &H) {
	if(G.getNumColors() != H.getNumColors()) {
		return false;
	}
	
	if(G.getCanonLabel() == H.getCanonLabel()) {
		return true;
	}
	else {
		return false;
	}
}


//------------------------------------
//-----Returns Subgraphs No Flags-----
//------------------------------------

//Return of std::vectors that could be used in restriction
//Almost the same as numSubgraphs, just with different output
//Needs no flags (used in the version with flags)
void returnSubgraphsNoFlags(const Graph &H, const Graph &G, std::vector<std::vector<int> > &output) {
	//No flags
	if((H.getSizeOfFlag() != 0) || (G.getSizeOfFlag() != 0)) {
		std::cout << "In returnSubgraphsNoFlags H and G can't have flags." << std::endl << std::endl;
		throw std::exception();
	}
	
	if(H.getNumColors() != G.getNumColors()) {
		return;
	}

	if((G.getN() < H.getN())) {
		return;
	}

	if(G.getN() == H.getN()) {
		if(isomorphic(H,G)) {
			std::vector<int> outputEntry;
			for(int i = 0; i < G.getN(); ++i) {
				outputEntry.push_back(i);
			}
			
			output.push_back(outputEntry);
			return;
		}
		
		else {
			return;
		}
	}
	
	//Single Vertex
	if(H.getN() == 1) {
		for(int i = 0; i < G.getN(); ++i) {
			std::vector<int> temp(G.getN(),-1);
			temp[i] = 0;
			output.push_back(temp);
		}
		
		return;
	}
	
	//Fastest way to do edges
	if(H.getN() == 2) {
		std::vector<int> temp(G.getN(),-1);
			
		for(int i = 0; i < G.getN(); ++i) {
			for(int j = i+1; j < G.getN(); ++j) {
				if(G.getEdgeColor(i,j) == H.getEdgeColor(0,1)) {
					temp[i] = 0;
					temp[j] = 1;
					output.push_back(temp);
					temp[i] = -1;
					temp[j] = -1;
				}
			}
		}
		
		return;
	}
	
	//Guarantees the least number of calls to the function (helps prune quickly)
	int vertex = -1;
	long long int minVal = (1 << G.getN());
	
	for(int i = 0; i < G.getNumOrbits(); ++i) {
		long long int temp1 = 0;
	
		for(int j = 0; j < H.getNumOrbits(); ++j) {
			int temp2 = 1;
			for(int c = 0; c < H.getNumColors(); ++c) {
				temp2 = temp2*choose(G.getDegree(G.getOrbit(i,0),c),H.getDegree(H.getOrbit(j,0),c));
			}
			temp1 = temp1 + temp2;
		}
		
		if (temp1 < minVal) {
			minVal = temp1;
			vertex = G.getOrbit(i,0);
		}
	}
	
	if(minVal == 0) {
		Graph Gcopy = G;
		Gcopy.removeVertex(vertex);
		std::vector<std::vector<int> > output2;
		
		returnSubgraphsNoFlags(H,Gcopy,output2);
		
		for(auto X : output2) {
			std::vector<int> Y;
			Y.resize(G.getN(),-1);
		
			for(int i = 0; i <= (int)X.size(); ++i) {
				//Add vertex back in
				if(i < vertex) {
					Y[i] = X[i];
				}
				
				else if(i > vertex) {
					Y[i] = X[i-1];
				}
			}
			
			output.push_back(Y);
		}
		
		return;
	}
	
	for(int i = 0; i < H.getNumOrbits(); ++i) {
		bool val = true;
		std::vector< std::vector < std::vector < int > > > possible; //Subsets of vertices in each collor which could create a copy of H
		int Hvertex = H.getOrbit(i,0);	
		bool dominating = true; //If our vertex is dominating we don't have to make final isomorphism call

		for(int c = 0; c < H.getNumColors(); ++c) {
			possible.push_back({});
			possible[c].clear();
				
			int Gdegree = G.getDegree(vertex,c);					
			int Hdegree = H.getDegree(Hvertex,c);
			
			if((Hdegree != 0) && (Hdegree != H.getN()-1)) {
				dominating = false;
			}
				
			//Can we actually find H in the nbrhd of vertex?
			if((Gdegree >= Hdegree) && (Hdegree > 0)) {
				//Recursively call returnSubgraphsNoFlags for each nbrhd
				
				std::vector<int> Grestriction(G.getN(),-1);
				std::vector<int> Hrestriction(H.getN(),-1);
				
				int temp = 0;
				for(int j = 0; j < G.getN(); ++j) {
					if(j != vertex) {
						if(G.getEdgeColor(vertex,j) == c) {
							Grestriction[j] = temp;
							++temp;
						}
					}
				}
				
				temp = 0;
				for(int j = 0; j < H.getN(); ++j) {
					if(j != Hvertex) {
						if(H.getEdgeColor(Hvertex,j) == c) {
							Hrestriction[j] = temp;
							++temp;
						}
					}
				}
				
				std::vector< std::vector< int > > possibleTemp; 
					returnSubgraphsNoFlags(H.restriction(Hrestriction),G.restriction(Grestriction),possibleTemp);
				//Make possibleTemp actually correspond to indices in G		
				for(int j = 0; j < (int)possibleTemp.size(); ++j) {
					int index = 0;
					std::vector<int> tempVec(G.getN(),-1);
					
					for(int k = 0; k < G.getN(); ++k) {
						if((k != vertex) && G.getEdgeColor(vertex,k) == c) {
							tempVec[k] = possibleTemp[j][index];
							++index;
						}
					}
					
					possible[c].push_back(tempVec);
				}
			}
			
			if((possible[c].size() == 0) && (Hdegree > 0)) {
				val = false;
				c = H.getNumColors();
			}
		}
		
		//Go through all possibilities and see if any of them combine to give a subgraph
		if(val) {	
			std::vector<int> maxVals; //Use in next_list
			std::vector<int> list;
			list.resize(H.getNumColors(),0);
				
			for(int c = 0; c < H.getNumColors(); ++c) {
				if(possible[c].size() == 0) {
					maxVals.push_back(0);
				}
					
				else {
					maxVals.push_back(possible[c].size()-1);
				}
			}
				
			do {
				std::vector<int> Grestriction;
				Grestriction.resize(G.getN(),-1);
				int index = 0;
					
				for(int c = 0; c < H.getNumColors(); ++c) {
					if(possible[c].size() > 0) {
						for(int j = 0; j < G.getN(); ++j) {
							if(possible[c][list[c]][j] != -1) {
								Grestriction[j] = index;
								++index;
							}
						}
					}
				}
					
				Grestriction[vertex] = index;
				
				if(!dominating) {
					if(isomorphic(G.restriction(Grestriction),H)) { //Not a cheap call, try to avoid it
						output.push_back(Grestriction);
					}
				}
				
				//If we have a dominating vertex don't need to do a second isomorphism check
				else {
					output.push_back(Grestriction);
				}
			} while(nextList(list, maxVals));
		}
	}
	
	Graph Gcopy = G;
	Gcopy.removeVertex(vertex);
	std::vector<std::vector<int> > output2;
	
	returnSubgraphsNoFlags(H,Gcopy,output2);
	
	for(auto X : output2) {
		std::vector<int> Y;
		Y.resize(G.getN(),-1);
	
		for(int i = 0; i <= (int)X.size(); ++i) {
			//add vertex back in
			if(i < vertex) {
				Y[i] = X[i];
			}
			
			else if(i > vertex) {
				Y[i] = X[i-1];
			}
		}
		
		output.push_back(Y);
	}
	
	return;
}

//---------------------------
//-----Returns Subgraphs-----
//---------------------------

//Return of std::vectors that could be used in restriction
void returnSubgraphs(const Graph &H, const Graph &G, std::vector<std::vector<int> > &output) {
	int HFlagSize = H.getSizeOfFlag();
	int GFlagSize = G.getSizeOfFlag();
	
	if(HFlagSize != GFlagSize) {
		return;
	}
	
	if(HFlagSize == 0) {
		returnSubgraphsNoFlags(H,G,output);
		return;
	}

	if(!isomorphic(H.getFlag(),G.getFlag())) {
		return;
	}
	
	if(H.getNumColors() != G.getNumColors()) {
		return;
	}

	if((G.getN() < H.getN())) {
		return;
	}

	if(G.getN() == H.getN()) {
		if(isomorphic(H,G)) {
			std::vector<int> outputEntry;
			for(int i = 0; i < G.getN(); ++i) {
				outputEntry.push_back(i);
			}
			
			output.push_back(outputEntry);
			return;
		}
		
		else {
			return;
		}
	}
	
	if(H.getN() == HFlagSize) {
		std::vector<int> temp(G.getN(),-1);
		for(int i = 0; i < HFlagSize; ++i) {
			temp[G.getFlagVertex(0)] = 0;
		}
		output.push_back(temp);

		return;
	}
	
	//We must have that HFlagSize == 1 based on other constraints
	if(H.getN() == 2) {
		for(int i = 0; i < G.getN(); ++i) {
			if((i != G.getFlagVertex(0)) && (G.getEdgeColor(i,G.getFlagVertex(0)) == H.getEdgeColor(0,1))) {
				std::vector<int> temp(G.getN(),-1);
				temp[G.getFlagVertex(0)] = 0;
				temp[i] = 1;
				
				output.push_back(temp);
			}	 
		}
	
		return;
	}
	
	//Remove Flags from H and G and then use returnSubgraphsNoFlags
	
	Graph HNoFlags = H;
	Graph GNoFlags = G;
	
	do {
		HNoFlags.removeVertex(HNoFlags.getFlagVertex(0));
	} while(HNoFlags.getSizeOfFlag() != 0);
	
	do {
		GNoFlags.removeVertex(GNoFlags.getFlagVertex(0));
	} while(GNoFlags.getSizeOfFlag() != 0);
	
	std::vector< std::vector <int> > outputNoFlags;
	
	returnSubgraphsNoFlags(HNoFlags,GNoFlags,outputNoFlags);
	
	//TODO prune vertices not in flags that don't connect correctly
	
	
	for(int i = 0; i < (int)outputNoFlags.size(); ++i) {
		std::vector<int> possibleOutput(G.getN(),-1);
		for(int i = 0; i < GFlagSize; ++i) {
			possibleOutput[G.getFlagVertex(i)] = i;
		}
		
		for(int j = 0; j < GNoFlags.getN(); ++j) {
			if(outputNoFlags[i][j] != -1) {
				int index = -1;
				
				for(int k = 0; k < G.getN(); ++k) {
					if(!G.flagVertex(k)) {
						++index;
					}
					
					if(index == j) {
						possibleOutput[k] = outputNoFlags[i][j] + GFlagSize;
						k = G.getN();
					}
				}
			}
		}
		
		if(isomorphic(G.restriction(possibleOutput),H)) {
			output.push_back(possibleOutput);
		}
	}
	
	return;
}

//--------------------------------------
//-----Number of Subgraphs No Flags-----
//--------------------------------------

//Copied almost directly from returnSubgraphsNoFlags
int numSubgraphsNoFlags(const Graph &H, const Graph &G) {
	//No flags
	if((H.getSizeOfFlag() != 0) || (G.getSizeOfFlag() != 0)) {
		std::cout << "In returnSubgraphsNoFlags H and G can't have flags." << std::endl << std::endl;
		throw std::exception();
	}
	
	if(H.getNumColors() != G.getNumColors()) {
		return 0;
	}

	if((G.getN() < H.getN())) {
		return 0;
	}

	if(G.getN() == H.getN()) {
		if(isomorphic(H,G)) {
			return 1;
		}
		
		else {
			return 0;
		}
	}
	
	//Single Vertex
	if(H.getN() == 1) {
		return G.getN();
	}
	
	//Fastest way to do edges
	if(H.getN() == 2) {
		int output = 0;
			
		for(int i = 0; i < G.getN(); ++i) {
			for(int j = i+1; j < G.getN(); ++j) {
				if(G.getEdgeColor(i,j) == H.getEdgeColor(0,1)) {
					++output;
				}
			}
		}
		
		return output;
	}
	
	//Guarantees the least number of calls to the function (helps prune quickly)
	int vertex = -1;
	long long int minVal = (1 << G.getN());
	
	for(int i = 0; i < G.getNumOrbits(); ++i) {
		long long int temp1 = 0;
	
		for(int j = 0; j < H.getNumOrbits(); ++j) {
			int temp2 = 1;
			for(int c = 0; c < H.getNumColors(); ++c) {
				temp2 = temp2*choose(G.getDegree(G.getOrbit(i,0),c),H.getDegree(H.getOrbit(j,0),c));
			}
			temp1 = temp1 + temp2;
		}
		
		if (temp1 < minVal) {
			minVal = temp1;
			vertex = G.getOrbit(i,0);
		}
	}
	
	if(minVal == 0) {
		Graph Gcopy = G;
		Gcopy.removeVertex(vertex);
		
		return numSubgraphsNoFlags(H,Gcopy);
	}
	
	int output = 0;
	
	for(int i = 0; i < H.getNumOrbits(); ++i) {
		bool val = true;
		std::vector< std::vector < std::vector < int > > > possible; //Subsets of vertices in each collor which could create a copy of H
		int Hvertex = H.getOrbit(i,0);
		bool dominating = true;

		for(int c = 0; c < H.getNumColors(); ++c) {
			possible.push_back({});
			possible[c].clear();
				
			int Gdegree = G.getDegree(vertex,c);					
			int Hdegree = H.getDegree(Hvertex,c);
			
			if((Hdegree != 0) && (Hdegree != H.getN()-1)) {
				dominating = false;
			}
			
			//Can we actually find H in the nbrhd of vertex?
			if((Gdegree >= Hdegree) && (Hdegree > 0)) {
				//Recursively call returnSubgraphsNoFlags for each nbrhd
				
				std::vector<int> Grestriction(G.getN(),-1);
				std::vector<int> Hrestriction(H.getN(),-1);
				
				int temp = 0;
				for(int j = 0; j < G.getN(); ++j) {
					if(j != vertex) {
						if(G.getEdgeColor(vertex,j) == c) {
							Grestriction[j] = temp;
							++temp;
						}
					}
				}
				
				temp = 0;
				for(int j = 0; j < H.getN(); ++j) {
					if(j != Hvertex) {
						if(H.getEdgeColor(Hvertex,j) == c) {
							Hrestriction[j] = temp;
							++temp;
						}
					}
				}
				
				std::vector< std::vector< int > > possibleTemp; 
					
				
				//Make possibleTemp actually correspond to indices in G		
				if(!dominating) {						
					returnSubgraphsNoFlags(H.restriction(Hrestriction),G.restriction(Grestriction),possibleTemp);
					for(int j = 0; j < (int)possibleTemp.size(); ++j) {
						int index = 0;
						std::vector<int> tempVec(G.getN(),-1);
						
						for(int k = 0; k < G.getN(); ++k) {
							if((k != vertex) && G.getEdgeColor(vertex,k) == c) {
								tempVec[k] = possibleTemp[j][index];
								++index;
							}
						}
						
						possible[c].push_back(tempVec);
					}
				}
				
				else {
					output = output+numSubgraphsNoFlags(H.restriction(Hrestriction),G.restriction(Grestriction));
					c = H.getNumColors();
				}	
			}
			
			if((possible[c].size() == 0) && (Hdegree > 0)) {
				val = false;
				c = H.getNumColors();
			}
		}
		
		//Go through all possibilities and see if any of them combine to give a subgraph
		if(val && !dominating) {	
			std::vector<int> maxVals; //Use in next_list
			std::vector<int> list;
			list.resize(H.getNumColors(),0);
				
			for(int c = 0; c < H.getNumColors(); ++c) {
				if(possible[c].size() == 0) {
					maxVals.push_back(0);
				}
					
				else {
					maxVals.push_back(possible[c].size()-1);
				}
			}
				
			do {
				std::vector<int> Grestriction;
				Grestriction.resize(G.getN(),-1);
				int index = 0;
					
				for(int c = 0; c < H.getNumColors(); ++c) {
					if(possible[c].size() > 0) {
						for(int j = 0; j < G.getN(); ++j) {
							if(possible[c][list[c]][j] != -1) {
								Grestriction[j] = index;
								++index;
							}
						}
					}
				}
					
				Grestriction[vertex] = index;
				
				if(isomorphic(G.restriction(Grestriction),H)) { //Not a cheap call, try to avoid it
					++output;
				}
				
			} while(nextList(list, maxVals));
		}
	}
	
	Graph Gcopy = G;
	Gcopy.removeVertex(vertex);
	std::vector<std::vector<int> > output2;
	
	return (output + numSubgraphsNoFlags(H,Gcopy));
}

//-----------------------------
//-----Number of Subgraphs-----
//-----------------------------

//Return of std::vectors that could be used in restriction
int numSubgraphs(const Graph &H, const Graph &G) {
	int HFlagSize = H.getSizeOfFlag();
	int GFlagSize = G.getSizeOfFlag();
	
	if(HFlagSize != GFlagSize) {
		return 0;
	}
	
	if(HFlagSize == 0) {
		return numSubgraphsNoFlags(H,G);
	}

	if(!isomorphic(H.getFlag(),G.getFlag())) {
		return 0;
	}
	
	if(H.getNumColors() != G.getNumColors()) {
		return 0;
	}

	if((G.getN() < H.getN())) {
		return 0;
	}

	if(G.getN() == H.getN()) {
		if(isomorphic(H,G)) {
			return 1;
		}
		
		else {
			return 0;
		}
	}
	
	if(H.getN() == HFlagSize) {
		return 1;
	}
	
	//We must have that HFlagSize == 1 based on other constraints
	if(H.getN() == 2) {
		int output = 0;
		for(int i = 0; i < G.getN(); ++i) {
			if((i != G.getFlagVertex(0)) && (G.getEdgeColor(i,G.getFlagVertex(0)) == H.getEdgeColor(0,1))) {
				++output;
			}	 
		}
	
		return output;
	}
	
	//Remove Flags from H and G and then use returnSubgraphsNoFlags
	
	Graph HNoFlags = H;
	Graph GNoFlags = G;
	
	do {
		HNoFlags.removeVertex(HNoFlags.getFlagVertex(0));
	} while(HNoFlags.getSizeOfFlag() != 0);
	
	do {
		GNoFlags.removeVertex(GNoFlags.getFlagVertex(0));
	} while(GNoFlags.getSizeOfFlag() != 0);
	
	std::vector< std::vector <int> > outputNoFlags;
	
	returnSubgraphsNoFlags(HNoFlags,GNoFlags,outputNoFlags);
	
	//TODO prune vertices not in flags that don't connect correctly
	
	int output = 0;
	for(int i = 0; i < (int)outputNoFlags.size(); ++i) {
		std::vector<int> possibleOutput(G.getN(),-1);
		for(int i = 0; i < GFlagSize; ++i) {
			possibleOutput[G.getFlagVertex(i)] = i;
		}
		
		for(int j = 0; j < GNoFlags.getN(); ++j) {
			if(outputNoFlags[i][j] != -1) {
				int index = -1;
				
				for(int k = 0; k < G.getN(); ++k) {
					if(!G.flagVertex(k)) {
						++index;
					}
					
					if(index == j) {
						possibleOutput[k] = outputNoFlags[i][j] + GFlagSize;
						k = G.getN();
					}
				}
			}
		}
		
		if(isomorphic(G.restriction(possibleOutput),H)) {
			++output;
		}
	}
	
	return output;
}

//------------------
//-----Subgraph-----
//------------------

//Copied almost directly from returnSubgraphsNoFlags
//This assumes we won't find a subgraph, could assume we do find a subgraph instead.
bool subgraph(const Graph &HWithFlag, const Graph &GWithFlag) {
	int HFlagSize = HWithFlag.getSizeOfFlag();
	int GFlagSize = GWithFlag.getSizeOfFlag();
	
	if(!isomorphic(HWithFlag.getFlag(),GWithFlag.getFlag())) {
		return false;
	}
	
	if(HWithFlag.getNumColors() != GWithFlag.getNumColors()) {
		return false;
	}

	if((GWithFlag.getN() < HWithFlag.getN())) {
		return false;
	}

	if(GWithFlag.getN() == HWithFlag.getN()) {
		if(isomorphic(HWithFlag,GWithFlag)) {
			return true;
		}
		
		else {
			return false;
		}
	}
	
	if(HWithFlag.getN() == HFlagSize) {
		return true;
	}
	
	//Single Vertex
	if(HWithFlag.getN() == 1) {
		return true;
	}
	
	//We must have that HFlagSize == 0 or 1 based on other constraints
	if(HWithFlag.getN() == 2) {
		if(HFlagSize == 0) {
			for(int i = 0; i < GWithFlag.getN(); ++i) {
				for(int j = i+1; j < GWithFlag.getN(); ++j) {
					if(GWithFlag.getEdgeColor(i,j) == HWithFlag.getEdgeColor(0,1)) {
						return true;
					}
				}
			}
		}
	
		else if(HFlagSize == 1) {
			for(int i = 0; i < GWithFlag.getN(); ++i) {
				if((i != GWithFlag.getFlagVertex(0)) && (GWithFlag.getEdgeColor(i,GWithFlag.getFlagVertex(0)) == HWithFlag.getEdgeColor(0,1))) {
					return true;
				}	 
			}
		}
	
		return false;
	}
	
	//Remove Flags from H and G and then use returnSubgraphsNoFlags
	
	Graph H = HWithFlag;
	Graph G = GWithFlag;
	
	if(HFlagSize != 0) {
		do {
			H.removeVertex(H.getFlagVertex(0));
		} while(H.getSizeOfFlag() != 0);
		
		do {
			G.removeVertex(G.getFlagVertex(0));
		} while(G.getSizeOfFlag() != 0);
	}

	
	//Guarantees the least number of calls to the function (helps prune quickly)
	int vertex = -1;
	long long int minVal = (1 << G.getN());
	
	for(int i = 0; i < G.getNumOrbits(); ++i) {
		long long int temp1 = 0;
	
		for(int j = 0; j < H.getNumOrbits(); ++j) {
			int temp2 = 1;
			for(int c = 0; c < H.getNumColors(); ++c) {
				temp2 = temp2*choose(G.getDegree(G.getOrbit(i,0),c),H.getDegree(H.getOrbit(j,0),c));
			}
			temp1 = temp1 + temp2;
		}
		
		if (temp1 < minVal) {
			minVal = temp1;
			vertex = G.getOrbit(i,0);
		}
	}
	
	if(minVal == 0) {
		Graph GWithFlagCopy = GWithFlag;
		int index = -1;
		int newVertex = -1;
		for(int i = 0; i < GWithFlag.getN(); ++i) {
			if(!GWithFlag.flagVertex(i)) {
				++index;
			}
			
			if(index == vertex) {
				newVertex = i;
				i = GWithFlag.getN();
			}
		}
		GWithFlagCopy.removeVertex(newVertex);
		
		return subgraph(HWithFlag,GWithFlagCopy);
	}
	
	for(int i = 0; i < H.getNumOrbits(); ++i) {
		bool val = true;
		std::vector< std::vector < std::vector < int > > > possible; //Subsets of vertices in each collor which could create a copy of H
		int Hvertex = H.getOrbit(i,0);
		bool dominating = true; 

		for(int c = 0; c < H.getNumColors(); ++c) {
			possible.push_back({});
			possible[c].clear();
				
			int Gdegree = G.getDegree(vertex,c);					
			int Hdegree = H.getDegree(Hvertex,c);
			
			if((Hdegree != 0) && (Hdegree != H.getN()-1)) {
				dominating = false;
			}
				
			//Can we actually find H in the nbrhd of vertex?
			if((Gdegree >= Hdegree) && (Hdegree > 0)) {
				//Recursively call returnSubgraphsNoFlags for each nbrhd
				
				std::vector<int> Grestriction(G.getN(),-1);
				std::vector<int> Hrestriction(H.getN(),-1);
				
				int temp = 0;
				for(int j = 0; j < G.getN(); ++j) {
					if(j != vertex) {
						if(G.getEdgeColor(vertex,j) == c) {
							Grestriction[j] = temp;
							++temp;
						}
					}
				}
				
				temp = 0;
				for(int j = 0; j < H.getN(); ++j) {
					if(j != Hvertex) {
						if(H.getEdgeColor(Hvertex,j) == c) {
							Hrestriction[j] = temp;
							++temp;
						}
					}
				}
				
				std::vector< std::vector< int > > possibleTemp; 
					returnSubgraphsNoFlags(H.restriction(Hrestriction),G.restriction(Grestriction),possibleTemp);
					
				if((possibleTemp.size() != 0) && dominating) {
					return true;
				}
				
				//Make possibleTemp actually correspond to indices in G		
				for(int j = 0; j < (int)possibleTemp.size(); ++j) {
					int index = 0;
					std::vector<int> tempVec(G.getN(),-1);
					
					for(int k = 0; k < G.getN(); ++k) {
						if((k != vertex) && G.getEdgeColor(vertex,k) == c) {
							tempVec[k] = possibleTemp[j][index];
							++index;
						}
					}
					
					possible[c].push_back(tempVec);
				}
			}
			
			if((possible[c].size() == 0) && (Hdegree > 0)) {
				val = false;
				c = H.getNumColors();
			}
		}
		
		//Go through all possibilities and see if any of them combine to give a subgraph
		if(val) {	
			std::vector<int> maxVals; //Use in next_list
			std::vector<int> list;
			list.resize(H.getNumColors(),0);
				
			for(int c = 0; c < H.getNumColors(); ++c) {
				if(possible[c].size() == 0) {
					maxVals.push_back(0);
				}
					
				else {
					maxVals.push_back(possible[c].size()-1);
				}
			}
				
			do {
				std::vector<int> Grestriction;
				Grestriction.resize(G.getN(),-1);
				int index = 0;
					
				for(int c = 0; c < H.getNumColors(); ++c) {
					if(possible[c].size() > 0) {
						for(int j = 0; j < G.getN(); ++j) {
							if(possible[c][list[c]][j] != -1) {
								Grestriction[j] = index;
								++index;
							}
						}
					}
				}
					
				Grestriction[vertex] = index;
					
				if(isomorphic(G.restriction(Grestriction),H)) {
					if(HFlagSize == 0) {
						return true;
					}

					//Now check if isomorphic with flag added back in
					std::vector<int> possibleSubgraph(GWithFlag.getN(),-1);
					
					for(int j = 0; j < GFlagSize; ++j) {
						possibleSubgraph[GWithFlag.getFlagVertex(j)] = j;
					}
					
					for(int j = 0; j < GWithFlag.getN(); ++j) {
						if(Grestriction[j] != -1) {
							int index = -1;
							
							for(int k = 0; k < GWithFlag.getN(); ++k) {
								if(!GWithFlag.flagVertex(k)) {
									++index;
								}
								
								if(index == j) {
									possibleSubgraph[k] = Grestriction[j] + GFlagSize;
									k = GWithFlag.getN();
								}
							}
						}
					}
					
					if(isomorphic(GWithFlag.restriction(possibleSubgraph),HWithFlag)) {
						return true;
					}
				}
			} while(nextList(list, maxVals));
		}
	}
	
	Graph GWithFlagCopy = GWithFlag;
	int index = -1;
	int newVertex = -1;
	for(int i = 0; i < GWithFlag.getN(); ++i) {
		if(!GWithFlag.flagVertex(i)) {
			++index;
		}
			
		if(index == vertex) {
			newVertex = i;
			i = GWithFlag.getN();
		}
	}
	GWithFlagCopy.removeVertex(newVertex);
		
	return subgraph(HWithFlag,GWithFlagCopy);
}


//-----------------------
//-----Expand Graphs-----
//-----------------------

//Expands a set of graphs of the same size by one vertex
//If input has flags, zeros applies to the graphs with the flags removed
//Need all flag vertices to be first and in order
//Doesn't check for duplicates in input (does in output) so if needed check that ahead of time
//TODO add version which preserves coefficients (shouldn't be too bad)
//TODO replace with set
std::vector<Graph> expandGraphs(const std::vector<Graph> &input, const std::vector<Graph> &zeros) {
	int size = input.size();
	
	if(size == 0) {
		return {};
	}
	
	int n = input[0].getN();
	int numColors = input[0].getNumColors();
	Graph flag = input[0].getFlag();
	int sizeOfFlag = flag.getN();
	
	for(int i = 1; i < size; ++i) {
		if(input[i].getN() != n) {
			std::cout << "In expandGraphs, everything in input must be the same size." << std::endl << std::endl;
			throw std::exception();
		}
		
		if(input[i].getNumColors() != numColors) {
			std::cout << "In expandGraphs, everything in input must have the same number of colors." << std::endl << std::endl;
			throw std::exception();
		}
		
		if(input[i].getFlag().getCanonLabel() != flag.getCanonLabel()) {
			std::cout << "In expandGraphs, everything in input must have the same flag." << std::endl << std::endl;
			throw std::exception();
		}
		
		for(int j = 0; j < sizeOfFlag; ++j) {
			if(input[i].getFlagVertex(j) != j) {
				std::cout << "In expandGraphs, in input, flag vertices must be first in the correct order." << std::endl << std::endl;
				throw std::exception();
			}
		}
	}
	
	for(int i = 0; i < (int)zeros.size(); ++i) {
		if(zeros[i].getNumColors() != numColors) {
			std::cout << "In expandGraphs, at least one zero doesn't have the correct number of colors." << std::endl << std::endl;
			throw std::exception();
		}
		
		if(zeros[i].getFlag().getN() != 0) {
			std::cout << "In expandGraphs, everything in zeros must not have a flag." << std::endl << std::endl;
			throw std::exception();
		}
	}
	
	//Actual generation
	std::unordered_set<std::string> canonLabels;
	
	
	//Needed everytime since they are always the same
	std::vector<Edge> flagEdges;
	
	for(int i = 0; i < sizeOfFlag-1; ++i) {
		for(int j = i+1; j < sizeOfFlag; ++j) {
			if(flag.getEdgeColor(i,j) != 0) {
				flagEdges.push_back({i,j,flag.getEdgeColor(i,j)});
			}
		}
	}

	std::vector<Graph> output;
	
	#pragma omp parallel for
	for(auto G: input) {
		std::vector<Edge> GEdges = flagEdges;
		std::vector<int> permanentZeroDegree(n+1,0); //Value 0 for 0,1,...,sizeOfFlag
		
		for(int i = 0; i < n-1; ++i) {
			for(int j = std::max(sizeOfFlag,i+1); j < n; ++j) {
				if(G.getEdgeColor(i,j) != 0) {
					GEdges.push_back({i,j,G.getEdgeColor(i,j)});
				}
				
				else {
					if(i >= sizeOfFlag) {
						++permanentZeroDegree[i];
						++permanentZeroDegree[j];
					}
				}
			}
		}
		
		for(int i = 0; i < myPow(numColors,n); ++i) {		
		 	int temp = i;
		 	std::vector<Edge> edges = GEdges;
		 	std::vector<int> zeroDegree = permanentZeroDegree;
		 	std::vector<int> check; //Used in an isomorphism check
		 			
		 	for(int j = 0; j < n; ++j) {
		 		if((temp % numColors) != 0) {
			 		edges.push_back({j,n,temp % numColors});
			 	}
			 	
			 	else {
			 		if(j >= sizeOfFlag) {
			 			++zeroDegree[n];
			 			++zeroDegree[j];
					}
			 	}
			 	
			 	check.push_back(temp % numColors);
			 	temp = temp / numColors;
			}
			
			std::vector<int> GvFlag;
				
			for(int j = 0; j < sizeOfFlag; ++j) {
				GvFlag.push_back(j);
			}
				
			Graph Gv(edges,n+1,numColors,GvFlag);
																	
					
			//Removes zeros
			 for(int j = 0; j < (int)zeros.size(); ++j) {
			 	if(subgraph(zeros[j],Gv)) {
			 		cont = false;
			 		j = zeros.size();
			 	}
				}
				 	
			#pragma omp critical
			{			
				if(cont && (canonLabels.count(Gv.getCanonLabel()) == 0)) {
					output.push_back(Gv);
			 		canonLabels.insert(Gv.getCanonLabel());
			 	}
		 	}	
		}
	}
	
	return output;
}


//-------------------------
//-----Generate Graphs-----
//-------------------------

//Generates all graphs of size n
std::vector<Graph> generate(const int n, const int numColors, const std::vector<Graph> &zeros = {}) {
	if(numColors < 2) {
		std::cout << "You need at least two colors in generate." << std::endl << std::endl;
		throw std::exception();
	}
	
	//Check all zeros have the correct number of colors
	for(int i = 0; i < (int)zeros.size(); ++i) {
		if(zeros[i].getNumColors() != numColors) {
			std::cout << "In generate all zeros must have the correct number of colors." << std::endl << std::endl;
			throw std::exception();
		}
	}
	
	
	Graph G({},1,numColors);
	std::vector<Graph> output = {G};
	for(int k = 2; k <= n; ++k) {
		output = expandGraphs(output,zeros);
	}
	
	return output;
}


//-----------------------
//-----Add All Flags-----
//-----------------------

//Gives all possible ways to add 
//Used in plainFlagAlgebra to get v std::vectors
//flag isn't a flag but a graph
std::vector<Graph> addAllFlags(const Graph &G, const Graph &flag) {
	int n = G.getN();
	std::unordered_set<std::string> canonLabels;
	std::vector<Graph> output;
	std::vector<std::vector<int> > restrictions;
	Graph noFlag = flag;
	noFlag.removeFlag();
	
	returnSubgraphs(noFlag,G,restrictions);
	
	if(G.getSizeOfFlag() != 0) {
		std::cout << "In addAllFlag G must not have a flag." << std::endl << std::endl;
		throw std::exception();
	}
	
	if(flag.getSizeOfFlag() != flag.getN()) {
		std::cout << "In addAllFlag flag must have all vertices as flag vertices." << std::endl << std::endl;
		throw std::exception();
	}
		
	//Convert restrictions into input for adding flag
	std::vector<std::vector <int> > possibleFlags;
	possibleFlags.resize(restrictions.size(), std::vector<int>(flag.getN()));
	
	for(int i = 0; i < (int)restrictions.size(); ++i) {
		int temp = 0;
		for(int j = 0; j < n; ++j) {	
			if(restrictions[i][j] != -1) {
				possibleFlags[i][temp] = j;
				++temp;
			}
		}
	}


	//Add flag in all possible ways
	for(int i = 0; i < (int)possibleFlags.size(); ++i) {
		sort(possibleFlags[i].begin(),possibleFlags[i].end());
		do {
			Graph GFlag = G;
			GFlag.setFlag(possibleFlags[i]);
			
			if(isomorphic(GFlag.getFlag(),flag)) {
				if(canonLabels.count(GFlag.getCanonLabel()) == 0) {
					canonLabels.insert(GFlag.getCanonLabel());
					GFlag.setCoefficient(Frac(1,1));
					output.push_back(GFlag);
				}
			}
		} while(next_permutation(possibleFlags[i].begin(),possibleFlags[i].end()));
	}
	
	return output;
}

//---------------------------------
//-----Generate First Index Map-----
//---------------------------------

//Generate map for first index in A used in plain flag algebra
std::unordered_map<std::string, int> firstIndex(const int n, const int numColors, const std::vector<Graph> &zeros) {
	std::unordered_map<std::string, int> output;
	int index = 0;
	
	for(int i = n/2; i <= n-1; ++i) {
		std::cout << "In firstIndex, iteration " << i << " out of " << n-1 << std::endl;
		int sizeOfFlag = 2*i-n;
		
		if(sizeOfFlag > 0) {
			std::vector<Graph> allGraphs = generate(sizeOfFlag, numColors, zeros);
			
			for(auto G: allGraphs) {
				output.insert(make_pair(G.getCanonLabel(),index));
				++index;	
			}
		}
	}
	
	//output.insert(std::make_pair("\n",0));
	
	return output;
}


//---------------------------------
//-----Generate Last Index Map-----
//---------------------------------

//Generate map for last two indices in A used in plain flag algebra
//Doesn't need to contain any info about the first index  (even though it is implicitly there)
//String = flag + flag comp + totalGraph
std::unordered_map<std::string, int> lastIndex(const int n, const int numColors, const std::vector<Graph> &zeros) {
	std::unordered_map<std::string, int> output;
	std::unordered_map<std::string, int> flags; //Need to say which flag we are on

	for(int i = n/2; i < n; ++i) {
		std::cout << "In lastIndex, iteration " << i << " out of " << n-1 << std::endl;
		int sizeOfFlag = 2*i-n;
		
		if(sizeOfFlag > 0) {
			std::vector<Graph> allGraphs = generate(i, numColors,zeros);
			
			for(auto G : allGraphs) {
				std::vector<int> sigma(sizeOfFlag);
				
				for(int j = 0; j < sizeOfFlag; ++j) {
					sigma[j] = j;
				}
				
				do {
					std::vector<int> restriction1(i,-1); //Flag
					std::vector<int> restriction2(i,-1); //Complements of each other
					
					for(int j = 0; j < sizeOfFlag; ++j) {
						restriction1[sigma[j]] = j;
					}
					
					int temp = 0;
					for(int j = 0; j < i; ++j) {
						if(restriction1[j] == -1) {
							restriction2[j] = temp;
							++temp;
						}
					}
					
					Graph G1 = G.restriction(restriction1);
					Graph G2 = G.restriction(restriction2);
					
					std::string tempStr = G1.getCanonLabel() + G2.getCanonLabel() + G.getCanonLabel();
					
					if(output.find(tempStr) == output.end()) {
						if(flags.find(G1.getCanonLabel()) == flags.end()) {
							output.insert(make_pair(tempStr,0));
							flags.insert(make_pair(G1.getCanonLabel(),1));
						}
					
						else {
							output.insert(make_pair(tempStr,flags.at(G1.getCanonLabel())));
							flags.at(G1.getCanonLabel()) = flags.at(G1.getCanonLabel()) + 1;
						}
					}
					
				} while(nextSubset(sigma,i,sizeOfFlag));
			}
		}
	}
	
	return output;
}


//----------------------------
//-----Generate all Flags-----
//----------------------------

//Used in generateV
//TODO make faster? Use iterated approach?
std::vector<Graph> generateFlags(const int n, const int numColors, const std::vector<Graph> &zeros) {
	std::vector<Graph> allGraphs = generate(n,numColors,zeros);
	std::vector<Graph> output;
	#pragma omp parallel 
	{
		std::vector<Graph> privateOutput;
				
		#pragma omp for nowait schedule(static) ordered
		for(int i = 0; i < (int)allGraphs.size(); ++i) {
			std::vector<int> sigma;
			std::unordered_set<std::string> canonLabels;
			
			for(int j = 0; j < n; ++j) {
				sigma.push_back(j);
			}
			
			do {
				Graph G = allGraphs[i];
				G.setFlag(sigma);
				
				if(canonLabels.count(G.getCanonLabel()) == 0) {
					canonLabels.insert(G.getCanonLabel());
					privateOutput.push_back(G);
				}
			} while(next_permutation(sigma.begin(),sigma.end()));
		}
		#pragma omp for schedule(static) ordered
    	for(int i=0; i<omp_get_num_threads(); i++) {
        #pragma omp ordered
        output.insert(output.end(), privateOutput.begin(), privateOutput.end());
    	}
	}
	
	return output;
}


//--------------------
//-----Generate V-----
//--------------------

//Use in plain flag algebra
std::vector< std::vector<Graph> > generateV(const int n, const int numColors, const std::vector<Graph> &zeros) {
	std::vector< std::vector< Graph > > output;

	for(int i = n/2; i <= n-1; ++i) {
		int sizeOfFlag = 2*i-n;
		
		if(sizeOfFlag > 0) {
			std::vector<Graph> flags = generateFlags(sizeOfFlag, numColors, zeros);
			#pragma omp parallel 
			{
				std::vector < std::vector <Graph> > privateOutput;
				
				#pragma omp for nowait schedule(static) ordered
				for(int j = 0; j < (int)flags.size(); ++j) {
					std::cout << "In generate v (" << i << ", " << j << ") out of (" << n-1 << ", " << flags.size() << ")" << std::endl;
					std::vector<Graph> expanded = {flags[j]};
				
					while((expanded[0].getN()*2 - sizeOfFlag) != n) {
						expanded = expandGraphs(expanded,zeros);
					}
					
					privateOutput.push_back(expanded);
				}
				
				#pragma omp for schedule(static) ordered
    			for(int i=0; i<omp_get_num_threads(); i++) {
     			   #pragma omp ordered
       			 output.insert(output.end(), privateOutput.begin(), privateOutput.end());
    			}
			}
		}
	}
	
	return output;
}


//----------------------------------------
//-----Generate All Graphs With Flags-----
//----------------------------------------

//Use in plain flag algebra
std::vector< std::vector<Graph> > generateAllGraphsWithFlags(const int n, const int numColors, const std::vector<Graph> &zeros) {
	//Need generateV and this to have the same indexing
	std::unordered_map<std::string,int> fFlag; //Takes flag canonLabels and maps to index in v
	int index = 0;
	
	for(int i = n/2; i <= n-1; ++i) {
		int sizeOfFlag = 2*i-n;
		
		if(sizeOfFlag > 0) {
			std::vector<Graph> flags = generateFlags(sizeOfFlag, numColors, zeros);
			for(int j = 0; j < (int)flags.size(); ++j) {
				fFlag[flags[j].getCanonLabel()] = index;
				++index;
			}
		}
	}
	
	//index is number of flags
	std::vector<Graph> allGraphs = generate(n,numColors,zeros);
	std::vector < std::vector < Graph > > output(index); //Index is number of flags
	#pragma omp parallel
	{
		std::vector<int> index2(index);
		
		for(int i = n/2; i <= n-1; ++i) {
			int sizeOfFlag = 2*i-n;
			
			if(sizeOfFlag > 0) {

				std::vector < std::vector < Graph > > privateOutput(index); //Index is number of flags
				
				#pragma omp for nowait schedule(dynamic)
				for(int j = 0; j < (int)allGraphs.size(); ++j) {
					std::cout << "In generateALLGraphsWIthFlags, iteration (" << i << ", " << j << ") out of (" << n-1 << ", " << allGraphs.size() << ")" << std::endl;
					std::unordered_map<std::string,int> fGraph;
					int den = choose(n, sizeOfFlag);
				
					std::vector<int> X(sizeOfFlag);
					for(int k = 0; k < sizeOfFlag; ++k) {
						X[k] = k;
					}
					
					do {
						do {
							Graph G = allGraphs[j];
							G.setFlag(X);
							Graph flag = G.getFlag();
							
							if(fGraph.find(G.getCanonLabel()) == fGraph.end()) {
								fGraph[G.getCanonLabel()] = index2[fFlag[flag.getCanonLabel()]];
								++index2[fFlag[flag.getCanonLabel()]];
								G.setCoefficient(Frac(1,den));
								privateOutput[fFlag[flag.getCanonLabel()]].push_back(G);
							}
							
							else {
								privateOutput[fFlag[flag.getCanonLabel()]][fGraph[G.getCanonLabel()]].setCoefficient(privateOutput[fFlag[flag.getCanonLabel()]][fGraph[G.getCanonLabel()]].getCoefficient()+Frac(1,den));
							}
							
						} while(next_permutation(X.begin(), X.end()));
					} while(nextSubset(X,n,sizeOfFlag));
				}
				
				#pragma omp critical
     			for(int j = 0; j < index; ++j) {
       			output[j].insert(output[j].end(), privateOutput[j].begin(), privateOutput[j].end());
       		}
			}
		}
	}
	
	return output;
}


//----------------------
//-----Random Graph-----
//----------------------

//p[0] + p[1] + ... + p[numColors-1] = 1
Graph randomGraph(const int n, const int numColors, const std::vector<double> &p, const int flagSize) {
	if((int)p.size() != numColors) {
		std::cout << "In randomGraph, p must have size numColors." << std::endl;
		std::cout << "P size = " << p.size() << std::endl;
		std::cout << "numColors = " << numColors << std::endl << std::endl;
		throw std::exception();
	}
	
	if(flagSize > n) {
		std::cout << "In randomGraph, we can't have flagSize > n." << std::endl;
		std::cout << "flagSize = " << flagSize << std::endl;
		std::cout << "n = " << n << std::endl;
		throw std::exception();
	}
	
	if(numColors < 0) {
		std::cout << "In randomGraph, we can't have numColors < 0." << std::endl;
		std::cout << "numColors = " << numColors << std::endl << std::endl;
		throw std::exception();
	}
	
	if(n < 0) {
		std::cout << "In randomGraph, we can't have n < 0." << std::endl;
		std::cout << "n = " << n << std::endl << std::endl;
		throw std::exception();
	}
	
	if(flagSize < 0) {
		std::cout << "In randomGraph, we can't have flagSize < 0." << std::endl;
		std::cout << "flagSize = " << flagSize << std::endl << std::endl;
		throw std::exception();
	}

	double temp = 0;
	for(int i = 0; i < numColors; ++i) {
		temp = temp + p[i];
	}
	
	if(abs(temp - 1) > 0.000001) {
		std::cout << "In randomGraph, p[0] + ... + p[numColors-1] = 1." << std::endl;
		for(int i = 0; i < numColors; ++i) {
			std::cout << "p[" << i << "] = " << p[i] << std::endl;
		}
		std::cout << std::endl;
		throw std::exception();
	}
	
	std::vector<Edge> Edges;
	//srand (time(NULL));
	
	for(int i = 0; i < n-1; ++i) {
		for(int j = i+1; j < n; ++j) {
			double ran = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
			double prob = 0;
			
			for(int k = 0; k < numColors; ++k) {
				prob = prob + p[k];
				
				if(ran < prob) {
					Edges.push_back({i,j,k});
					k = numColors;
				}
			}
		}
	}
	
	
	
	if(flagSize > 0) {
		std::vector<int> flag;
		
		for(int i = 0; i < flagSize; ++i) {
			int ran = rand()%n;
			
			bool inFlag = false;
			for(int j = 0; j < (int)flag.size(); ++j) {
				if(flag[j] == ran) {
					inFlag = true;
					j = flag.size();
				}
			}
			
			if(!inFlag) {
				flag.push_back(ran);
			}
			
			else {
			 --i;
			}
		}
		
		
		Graph G(Edges,n,numColors,flag);
		return G;
	}
	
	else {
		Graph G(Edges,n,numColors);
		return G;
	}
}


//------------------------------
//-----Uniform Random Graph-----
//------------------------------

Graph uniformRandomGraph(const int n, const int numColors = 2, const int flagSize = 0) {
	std::vector<double> p(numColors,1./numColors);

	return randomGraph(n, numColors, p, flagSize);
}


//------------------------
//------------------------
//-----Equation Class-----
//------------------------
//------------------------

class Equation {

	private: 
	
		std::vector<Graph> variables;
		int numVariables;
		std::vector<Graph> zeros; //No Flag on zeros (remove flag from Graph to check subgraph)
		int numZeros;
		int numColors = -1;
		Frac ans = Frac(0,1);
		int type; //0 for =, 1 for <= (if want >= make everything negative)
		bool valid = true; //Way of dealing with equations with zero variables
		
	public:
	
		//---------------------
		//-----Constructor-----
		//---------------------
		
		//ans must be a fraction
		//TYPE 0 ==, TYPE 1 <=
		Equation(std::vector<Graph> const &VARIABLES, std::vector<Graph> const &ZEROS = {}, Frac ANS = Frac(1,1), int TYPE = 0) {
			variables = VARIABLES;
			zeros = ZEROS;
			ans = ANS;
			type = TYPE;
			numVariables = variables.size();
			numZeros = zeros.size();
			
			//Check type is -1,0,or 1
			if((type != 0) && (type != 1)) {
				std::cout << "In equation constructor, type must be 0, or 1." << std::endl << std::endl;
				throw std::exception();
			}
			
			//Check if all graphs have same number of colors
			if(numVariables != 0) {
				numColors = variables[0].getNumColors();
			}
			
			for(auto G: variables) {
				if(G.getNumColors() != numColors) {
					std::cout << "All graphs in Equation Constructor must have same number of colors." << std::endl << std::endl;
					throw std::exception();
				}
			}
			
			//Remove all variables with 0 coefficients
			//Coefficients don't matter for zeros (even if it is zero)
			int index = 0;

			while(index < numVariables) {
				if(variables[index].getCoefficient() == 0) {
					variables.erase(variables.begin()+index);
					--numVariables;
				}
				++index;
			}

			//Checks that all flags in variables are the same
			//Could probably make faster after canonically relabelling 

			if(numVariables != 0) {
				for(int i = 1; i < numVariables; ++i) {
					if(!isomorphic(variables[i].getFlag(), variables[0].getFlag())) {
						std::cout << "All flag have to be isomorphic in equation constructor." << std::endl << std::endl; 
						throw std::exception();
					}
				}
			}
			
			fixZeros(); 
			
			combine();
		}
		
		
		//-----------------
		//-----Combine-----
		//-----------------
		
		//Checks for isomorphisms and the adds coefficients
		//Private?
		void combine() {
			std::unordered_map<std::string, int> map;
			
			for(int i = 0; i < numVariables; ++i) {
				if(map.find(variables[i].getCanonLabel()) == map.end()) {
					map[variables[i].getCanonLabel()] = i;
				}
					
				else{ 
					int firstIndex = map[variables[i].getCanonLabel()];
					variables[firstIndex].setCoefficient(variables[firstIndex].getCoefficient() + variables[i].getCoefficient());
					variables.erase(variables.begin()+i);
					--i;
					--numVariables;
				}
			}

			//Do I want this?
			for(int i = 0; i < numVariables; ++i) {
				if((variables[i].getN() == 1) && (variables[i].getSizeOfFlag() == 0)) {
					ans = ans - variables[i].getCoefficient(); 
					variables.erase(variables.begin()+i);
					i = numVariables;
					--numVariables;
				}
			}
		}	
		
		//-------------------
		//-----Fix Zeros-----
		//-------------------
		
		//Private?
		void fixZeros() {
			//Check if all graphs in zeros have same number of colors
			if(numVariables > 0) {
				for(auto G: zeros) {
					if(G.getNumColors() != numColors) {
						std::cout << "All zeros in Equation Constructor must have same number of colors as variables." << std::endl << std::endl;
						throw std::exception();
					}
				}
			}
			
			if((numVariables == 0) && (numZeros != 0)) {
				numColors = zeros[0].getNumColors();
				
				for(auto G: zeros) {
					if(G.getNumColors() != numColors) {
						std::cout << "All zeros in Equation Constructor must have same number of colors." << std::endl << std::endl;
						throw std::exception();
					}
				}
			}
			
			//Check that all zeros don't have flag
			//I'm not entirely sure that this needs to be true, but I can't think of any reasonable examples were it isn't though
			for(auto G: zeros) {
				if(G.getSizeOfFlag() != 0) {
					std::cout << "Zeros can't have any flags." << std::endl << std::endl;
					throw std::exception();
				}
			}
		
			//Make the list of zeros self contained
			//Also checks for isomorphisms at this point
			if(numZeros != 0) {
				for(int i = 0; i < numZeros; ++i) {
					for(int j = 0; j < numZeros; ++j) {
						if(i != j) {
							if(subgraph(zeros[i], zeros[j])) {
								zeros.erase(zeros.begin()+j);
								--j;
								--numZeros;
							}
						}
					}
				}
			}
			
			//Remove variables with zeros as subgraphs
			for(int i = 0; i < numVariables; ++i) {
				Graph G = variables[i].clone(); //Don't want to remove flag from actual graph
				G.removeFlag();
				bool isZero = false;
				
				for(auto H : zeros) {
					if((!isZero) && subgraph(H,G)) {
						isZero = true;
					}
				}
				
				if(isZero) {
					variables.erase(variables.begin()+i);
					--numVariables;
					--i; 
				}
			}
		}
		
		
		//--------------------------
		//-----Get numVariables-----
		//--------------------------
		
		int getNumVariables() const {
			return numVariables;
		}
		
		
		//-----------------------
		//-----Get numColors-----
		//-----------------------
		
		int getNumColors() const {
			return numColors;
		}
		
		
		//----------------------
		//-----Get numZeros-----
		//----------------------
		
		int getNumZeros() const {
			return numZeros;
		}
		
		
		//------------------
		//-----Get Type-----
		//------------------
		
		int getType() const {
			return type;
		}
		
		
		//-----------------
		//-----setType-----
		//-----------------
		
		void setType(int i) {
			type = i;
			
			valid = checkValid();
			
			if(!valid) {
				std::cout << "Equation not valid - no variables, and it is false." << std::endl << std::endl;
				throw std::exception(); 
			}
		}
		
		
		//-----------------
		//-----Get Ans-----
		//-----------------
		
		Frac getAns() const {
			return ans;
		}
		
		
		//-----------------
		//-----Set Ans-----
		//-----------------
		
		void setAns(Frac ANS) {
			ans = ANS;
			
			valid = checkValid();
			
			if(!valid) {
				std::cout << "Equation not valid - no variables, and it is false." << std::endl << std::endl;
				throw std::exception(); 
			}
		}
		
		
		//-------------------------
		//-----Get Coefficient-----
		//-------------------------
		
		Frac getCoefficient(int i) const {
			return variables[i].getCoefficient();
		}
		
		
		//-------------------------
		//-----Set Coefficient-----
		//-------------------------
		
		void setCoefficient(int i, Frac COEFFICIENT) {
			variables[i].setCoefficient(COEFFICIENT);
		}
		
		void setCoefficient(int i, int COEFFICIENT) {
			variables[i].setCoefficient(COEFFICIENT);
		}
		
		
		//----------------------
		//-----Get Variable-----
		//----------------------
		
		Graph getVariable(int i) const {
			if((i < 0) || (i >= numVariables)) {
				std::cout << "Trying to get variable outside of range." << std::endl;
				std::cout << "Index = " << i << std::endl;
				std::cout << "Number of Variables = " << numVariables << std::endl;
				throw std::exception();
			}
			
			return variables[i];
		}
		
		
		//------------------
		//-----Get Zero-----
		//------------------
		
		Graph getZero(int i) const {
			if((i < 0) || (i >= numZeros)) {
				std::cout << "Trying to get zero outside of range: i = " << i << ", numZeros = " << numVariables << std::endl << std::endl;
				throw std::exception();
			}
			
			return zeros[i];
		}
		
		
		//------------------
		//-----Add Zero-----
		//------------------
		
		//Not allowed to remove a zero b/c we may have already deleted some flags with that as a subgraph
		void addZero(Graph G) {
			++numZeros;
			zeros.push_back(G);
			
			fixZeros(); //Slightly more than we need - could fix
			
			valid = checkValid();
			
			if(!valid) {
				std::cout << "Equation not valid - no variables, and it is false." << std::endl << std::endl;
				throw std::exception(); 
			}
		}
		
		
		//-------------------------
		//-----Print Variables-----
		//-------------------------
		
		void printVariables() const {
			if(numVariables == 0) {
				std::cout << "There are no variables." << std::endl << std::endl;
			}
			
			else {
				for(int i = 0; i < numVariables; ++i) {
					std::cout << "Variable " << i << " has adjacency matrix: " << std::endl;
					variables[i].printAdjMat();
					if(variables[i].getSizeOfFlag() != 0) {
						std::cout << "The flag vertices are: ";
						variables[i].printFlagVertices();
					}
					else {
						std::cout << "There are no flag vertices." << std::endl;
					}
					std::cout << "It has coefficient: " << variables[i].getCoefficient() << std::endl << std::endl;
					//Could add this back in, but it doesn't seem necessary unless debugging
					//std::cout << "Its canon label is: " << variables[i].getCanon() << std::endl << std::endl;
				}
			}
		}
		
		
		//----------------------
		//-----Prints Zeros-----
		//----------------------
		
		void printZeros() const {
			if(numZeros == 0) {
				std::cout << "There are no zeros." << std::endl << std::endl;
			}
			
			else {
				for(int i = 0; i < numZeros; ++i) {
					std::cout << "Zero " << i << " has adjacency matrix: " << std::endl;
					zeros[i].printAdjMat();
					std::cout << std::endl;
				}
			}
		}
		
		
		//---------------------
		//-----Check Valid-----
		//---------------------
		
		//Mostly to deal with equations with no variables
		bool checkValid() const {
			if(numVariables > 0) {
				return true;
			}
			
			else if(type == 1) {
				return (0 <= ans);
			}
			
			else {
				return (0 == ans);
			}
		}
		
				
		//-----------------
		//-----Average-----
		//-----------------
		
		void averageAll() {
			for(int i = 0; i < numVariables; ++i) {
				variables[i].averageAll();
			}
			
			//Combine variables that are now the same
			combine();
		}
		
		//----------------
		//-----Negate-----
		//----------------
		
		//Negates coefficients and answer, but doesn't change type
		void negate() {
			for(int i = 0; i < numVariables; ++i) {
				variables[i].setCoefficient(-variables[i].getCoefficient());
			}
			
			ans = -ans;
		}
};


//---------------------------------------
//---------------------------------------
//-----Functions Involving Equations----- 
//---------------------------------------
//---------------------------------------


//-------------------
//-----Equaltiy------
//-------------------

//Mostly used in debugging
//TODO use sets
bool operator==(const Equation &eq1, const Equation &eq2) {
	std::vector<Graph> variables;
	std::vector<Graph> zeros; //No Flag on zeros (remove flag from Graph to check subgraph)
	
	if(eq1.getType() != eq2.getType()) {
		return false;
	}	
	
	if(eq1.getAns() != eq2.getAns()) {
		return false;
	}
	
	if(eq1.getNumColors() != eq2.getNumColors()) {
		return false;
	}
	
	if(eq1.getNumZeros() != eq2.getNumZeros()) {
		return false;
	}
	
	if(eq1.getNumVariables() != eq2.getNumVariables()) {
		return false;
	}
	
	for(int i = 0; i < eq1.getNumZeros(); ++i) {
		bool check = false;
	
		for(int j = 0; j < eq2.getNumZeros(); ++j) {
			if(isomorphic(eq1.getZero(i),eq2.getZero(j))) {
				check = true;
				j = eq2.getNumZeros();
			}
		}
		
		if(!check) {
			return false;
		}
	}
	
	for(int i = 0; i < eq1.getNumVariables(); ++i) {
		bool check = false;
	
		for(int j = 0; j < eq2.getNumVariables(); ++j) {
			if(isomorphic(eq1.getVariable(i),eq2.getVariable(j)) && (eq1.getVariable(i).getCoefficient() == eq2.getVariable(i).getCoefficient())) {
				check = true;
				j = eq2.getNumVariables();
			}
		}
		
		if(!check) {
			return false;
		}
	}
	
	return true;
}


//------------------
//-----Addition-----
//------------------

Equation operator+(const Equation &eq1, const Equation &eq2) {
	//Deal with type
	int type;
	if((eq1.getType() == 1) || (eq2.getType() == 1)) {
		type = 1;
	}
	
	else {
		type = 0;
	}
	
	//Check that zeros are the same
	int numZeros = eq1.getNumZeros();
	
	if(numZeros != eq2.getNumZeros()) {
		std::cout << "When adding equations, they must have the same set of zeros." << std::endl << std::endl;
		throw std::exception();
	}
	
	std::vector<Graph> zeros;
	std::unordered_set<std::string> zeros1;
	std::unordered_set<std::string> zeros2;
	
	for(int i = 0; i < numZeros; ++i) {
		zeros.push_back(eq1.getZero(i));
		
		zeros1.insert(eq1.getZero(i).getCanonLabel());
		zeros2.insert(eq2.getZero(i).getCanonLabel());
	}
	
	if(zeros1 != zeros2) {
		std::cout << "When adding equations, they must have the same set of zeros." << std::endl << std::endl;
		throw std::exception();
	}
	
	//Constructor actually does the adding already
	std::vector<Graph> variables;
	
	for(int i = 0; i < eq1.getNumVariables(); ++i) {
		variables.push_back(eq1.getVariable(i));
	}
	
	for(int i = 0; i < eq2.getNumVariables(); ++i) {
		variables.push_back(eq2.getVariable(i));
	}
	
	return Equation(variables, zeros, eq1.getAns() + eq2.getAns(), type);
}

//If adding Graph to Equation, like adding G = 0 to it
Equation operator+(const Graph &G, const Equation &eq) {
	std::vector<Graph> variables = {G};
	
	std::vector<Graph> zeros;
	for(int i = 0; i < eq.getNumZeros(); ++i) {
		zeros.push_back(eq.getZero(i));
	}
	
	Frac ans(0,1);
	
	return (eq + Equation(variables, zeros, ans, 0));
}

Equation operator+(Equation eq, Graph G) {
	return(G + eq);
}


//------------------------
//-----Multiplication-----
//------------------------

//First Define for graphs - put into equation with = 0, and no zeros
//Really more of a std::vector (or set) of graphs, but in this form it allows us to use overridden +
//TODO add edge between parts in a smarter way?
Equation multiply(const Graph &G1, const Graph &H1, const std::vector<Graph> &zeros) {
	std::vector<int> HReorder;
	std::vector<int> GReorder;
	
	//Make sure all the flags are first
	HReorder.resize(H1.getN(),H1.getN());
	GReorder.resize(G1.getN(),G1.getN());
	
	int temp = 0;
	for(int i = 0; i < H1.getSizeOfFlag(); ++i) {
		HReorder[H1.getFlagVertex(i)] = temp;
		++temp;
	}
	
	for(int i = 0; i < H1.getN(); ++i) {
		if(!H1.isFlag(i)) {
			HReorder[i] = temp;
			++temp;
		}
	}
	
	temp = 0;
	for(int i = 0; i < G1.getSizeOfFlag(); ++i) {
		GReorder[G1.getFlagVertex(i)] = temp;
		++temp;
	}
	
	for(int i = 0; i < G1.getN(); ++i) {
		if(!G1.isFlag(i)) {
			GReorder[i] = temp;
			++temp;
		}
	}
	
	Graph H = H1.restriction(HReorder);
	Graph G = G1.restriction(GReorder);

	Frac GCoefficient = G1.getCoefficient();
	Frac HCoefficient = H1.getCoefficient();

	//Double up some checks as in actual operator but that's ok b/c this can be used on it's own
	
	//Flags must be the same	
	if(!isomorphic(G.getFlag(),H.getFlag())) {
		std::cout << "Flags must be the same in multiplication." << std::endl << std::endl;
		G.printEdges();
		G.printFlag();
		H.printEdges();
		H.printFlag();
		throw std::exception();
	}
	
	//Check that numColors are the same
	if(G.getNumColors() != H.getNumColors()) {
		std::cout << "Graphs in multiplication must have the same number of colors." << std::endl << std::endl;
		throw std::exception();
	}
	
	//Can't multiply empty graph
	if((G.getN() == 0) || (H.getN() == 0)) {
		std::cout << "Can't multiply empty graphs." << std::endl << std::endl;
		throw std::exception();
	}
	
	//Need to reorder vertices so flags are first
	int sizeOfFlag = G.getSizeOfFlag();
	int c = G.getNumColors();
	int nG = G.getN();
	int nH = H.getN();
	
	std::vector<Graph> variables;
	std::unordered_set<std::string> variablesSet;

	//Create all possible variables
	#pragma omp parallel for
	for(int i = 0; i < myPow(c,(nG - sizeOfFlag)*(nH-sizeOfFlag)); ++i) { //c-ary mask
		//Convert i to c-ary number
		std::vector<int> ary;
		int temp = i;
		
		for(int j = 0; j < (nG-sizeOfFlag)*(nH-sizeOfFlag); ++ j) {
			int temp2 = temp/myPow(c,(nG-sizeOfFlag)*(nH-sizeOfFlag)-j-1); //Need to be careful with rounding? 
		   ary.push_back(temp2);
		   temp = temp-temp2*myPow(c,(nG-sizeOfFlag)*(nH-sizeOfFlag)-j-1);
		   
		   if((temp2 < 0) || (temp2 >= c)) { //This really shouldn't happen I've checked, but better safe than sorry
		   	std::cout << "Something went wrong in graph multiplication." << std::endl << std::endl;
		   	throw std::exception();
		   }
		}
		
		//Create edge list for new Graph
		std::vector<Edge> edges;
		
		//First vertices G
		for(int j = 0; j < nG-1; ++j) {
			for(int k = j+1; k < nG; ++k) {
				if(G.getEdgeColor(j,k) != 0) {
					edges.push_back({j,k,G.getEdgeColor(j,k)});
				}
			}
		}

		//Flag Plus Last Vertices = H
		for(int j = 0; j < nH-1; ++j) {
			for(int k = j+1; k < nH; ++k) {
				if(H.getEdgeColor(j,k) != 0) {
					if((j < sizeOfFlag) && (k < sizeOfFlag)) {
						edges.push_back({j,k,H.getEdgeColor(j,k)});
					}
					
					else if((j < sizeOfFlag) && (k >= sizeOfFlag)) {
						edges.push_back({j,k+nG-sizeOfFlag,H.getEdgeColor(j,k)});
					}
					
					else {
						edges.push_back({j+nG-sizeOfFlag,k+nG-sizeOfFlag,H.getEdgeColor(j,k)});
					}
				}
			}
		}
		
		//Add edges from i
		//Easier to just iterate through than to make a function
		int index = 0;
		for(int j = nG; j < nG+nH-sizeOfFlag; ++j) {
			for(int k = sizeOfFlag; k < nG; ++k) {
				if(ary[index] != 0) {
					edges.push_back({j,k,ary[index]});
				}
				++index;
			}
		}
		
		//Create Flag
		std::vector<int> flag;
		
		//Can do this because they are in canonical form
		for(int j = 0; j < sizeOfFlag; ++j) {
			flag.push_back(G.getFlagVertex(j));
		}
		
		Graph GH(edges, nG+nH-sizeOfFlag, c,flag);
		
		#pragma omp critical
		{
			if(variablesSet.count(GH.getCanonLabel()) == 0) {
				variablesSet.insert(GH.getCanonLabel());
				variables.push_back(GH);
			}
		}
	}

	Equation eq(variables,zeros, Frac(0,1), 0); //When converting to equation, removes isomorphisms
	//Makes coefficients correct

	for(int i = 0; i < eq.getNumVariables(); ++i) {
		int num = 0;
		
		std::vector<std::vector<int> > subgraphs; 
		returnSubgraphs(G,eq.getVariable(i),subgraphs);
		
		for(auto X: subgraphs) {
			std::vector<int> restrictionComp;
			restrictionComp.resize(nG+nH-sizeOfFlag,-1);
			int index = 0;
			
			for(int i = 0; i < sizeOfFlag; ++i) {
				restrictionComp[i] = index;
				++index;
			}
			
			for(int i = sizeOfFlag; i < nG+nH-sizeOfFlag; ++i) {
				if(X[i] == -1) {
					restrictionComp[i] = index;
					++index;
				}
			}
			
			if(isomorphic(eq.getVariable(i).restriction(restrictionComp), H)) {
				++num;
			}
		}
	
		
		int den = choose(nG+nH-2*sizeOfFlag,nG-sizeOfFlag);
		eq.setCoefficient(i,Frac(num,den)*GCoefficient*HCoefficient);
	}

	return eq;
}

Equation operator*(const Equation &eq1, const Equation &eq2) {
	//Make sure types are right
	if((eq1.getType() == 1) && (eq2.getType() == 1)) {
		std::cout << "In multiplication, both types can't be 1, inequalities possibly not preserved." << std::endl << std::endl;
		throw std::exception();
	}
	
	int type;
	
	if((eq1.getType() == 0) || (eq2.getType() == 1)) {
		type = 1; //May need to make things negative (see end of function)
	}
	
	else {
		type = 0;
	}
	
	//Deal with empty equation(s)
	//Not sure about throwing an excpetion - may change
	if((eq1.getNumVariables() == 0) || (eq2.getNumVariables() == 0)) {
		std::cout << "Can't multiply two equations together when one is empty." << std::endl << std::endl;
		throw std::exception();
	}
	
	//Check the flags are the same
	if(!isomorphic(eq1.getVariable(0).getFlag(),eq2.getVariable(0).getFlag())) {
		std::cout << "Flags must be the same in multiplication." << std::endl << std::endl;
		throw std::exception();
	}
	
	//Check that numColors are the same
	if(eq1.getVariable(0).getNumColors() != eq2.getVariable(0).getNumColors()) {
		std::cout << "Graphs in multiplication must have the same number of colors." << std::endl << std::endl;
		throw std::exception();
	}
	
	//Check if Zeros are the same
	int numZeros = eq1.getNumZeros();
	
	if(numZeros != eq2.getNumZeros()) {
		std::cout << "When adding equations, they must have the same set of zeros." << std::endl << std::endl;
		throw std::exception();
	}
	
	std::vector<Graph> zeros;
	
	bool val;
	for(int i = 0; i < numZeros; ++i) {
		val = false;
		for(int j = 0; j < numZeros; ++j) {
			if(isomorphic(eq1.getZero(i),eq2.getZero(j))) {
				val = true;
				j = numZeros;
				zeros.push_back(eq1.getZero(i));
			}
		}
		
		if(!val) {
			std::cout << "When adding equations, they must have the same set of zeros." << std::endl << std::endl;
			throw std::exception();
		}
	}

	
	//Distributive
	Equation eq = multiply(eq1.getVariable(0), eq2.getVariable(0), zeros);

	//#pragma omp parallel for collapse(2) schedule(dynamic)
	for(int i = 0; i < eq1.getNumVariables(); ++i) {
		for(int j = 0; j < eq2.getNumVariables(); ++j) {
			std::cout << "In multiply (" << i << ", " << j << ") out of (" << eq1.getNumVariables() << ", " << eq2.getNumVariables() << ")." << std::endl;
			if((i != 0) || (j != 0)) {
				Equation toAdd = multiply(eq1.getVariable(i), eq2.getVariable(j), zeros);
				//#pragma omp critical
				eq = eq + toAdd;
			}
		}
	}
	
	eq.setType(type);
	eq.setAns(eq1.getAns() * eq2.getAns());
	
	
	//Switch coefficients if we have type 0 and 1, and a negative answer for the equation of type 0
	if((eq1.getType() == 0) && (eq2.getType() == 1) && (eq1.getAns() < 0)) {
		for(int i = 0; i < eq.getNumVariables(); ++i) {
			eq.setCoefficient(i,-eq.getCoefficient(i));
		}
	}
	
	if((eq1.getType() == 1) && (eq2.getType() == 0) && (eq2.getAns() < 0)) {
		for(int i = 0; i < eq.getNumVariables(); ++i) {
			eq.setCoefficient(i,-eq.getCoefficient(i));
		}
	}
	
	return eq; 
}

//Scale by constant
Equation operator*(Frac k, const Equation &eq) {
	Equation output = eq;

	for(int i = 0; i < eq.getNumVariables(); ++i) {
		Frac j = output.getVariable(i).getCoefficient();
		output.setCoefficient(i,j*k);
	}
	
	return output;
}

Equation operator*(int k, const Equation &eq) {
	Frac temp(k,1);
	return temp*eq;
}

//-------------------------
//-----Resize Equation-----
//-------------------------

//Currently only works for equations with no flags
//Make this a class function?
//Use generate to get generated
Equation resize(const Equation &eq, const std::vector<Graph> &generated) {
	if(eq.getNumVariables() == 0) {
		std::cout << "In resize, need at least one graph in the equation to resize." << std::endl << std::endl;
		throw std::exception();
	}
	
	for(int i = 0; i < eq.getNumVariables(); ++i) {
		if(eq.getVariable(i).getFlag().getN() != 0) {
			std::cout << "In resize, we don't allow flags (there is no reason for this I'm just lazy, so future me you can implement this)." << std::endl << std::endl;
			throw std::exception();
		}
		
		if(eq.getVariable(i).getN() > generated[0].getN()) {
			std::cout << "In resize, we can't have the variables having larger size than things in generated." << std::endl;
			throw std::exception();
		}
	}
	
	std::vector<Graph> zeros;
	for(int i = 0; i < eq.getNumZeros(); ++i) {
		zeros.push_back(eq.getZero(i));
	}	
	
	Frac ans = eq.getAns();
	
	
	Equation output({},zeros,eq.getAns(),eq.getType());
	
	#pragma omp parallel for
	for(int i = 0; i < eq.getNumVariables(); ++i) {
		std::vector<Graph> resized;
		for(int j = 0; j < int(generated.size()); ++j) {
			Graph H = generated[j];
			int num = numSubgraphs(eq.getVariable(i),generated[j]);
			int den = choose(generated[j].getN(), eq.getVariable(i).getN());
			
			
			if(num != 0) {
				H.setCoefficient(Frac(num,den)*eq.getVariable(i).getCoefficient());
				resized.push_back(H);
			}
		}
		
		//I'm not a huge fan of this implementation because it rechecks for zeros
		//Maybe make an "unsafe" mode for creating equations to not check for zeros?
		Equation temp(resized,zeros,eq.getAns(),eq.getType());
		
		#pragma omp critical
		{
			output = output+temp;
		}
	}
	
	output.setAns(ans);
	
	return output;
}


//----------------------------
//-----Plain Flag Algebra-----
//----------------------------

//Prints to either Mosek files or CSDP files
//f can be thought of as a linear combo of all graphs that we want to max/min
//Rather than taking a v this function takes the number of vertices to compute on (n)
//If using maximize, use 1-answer
void plainFlagAlgebra(std::vector<Graph> &f, int n, std::vector<Graph> &zeros, std::vector<Equation> &known, bool maximize = true) {
	std::cout << "Starting plainFlagAlgebra." << std::endl;
	
	int fSize = f.size();
	int zerosSize = zeros.size();
	int knownSize = known.size();
		
	if(n <= 1) {
		std::cout << "Plain flag algebra method not set up for graphs with fewer than two vertices." << std::endl << std::endl;
		throw std::exception();
	}

	if(fSize == 0) {
		std::cout << "Need at least one graph in f in plainFlagAlgebra." << std::endl << std::endl;
		throw std::exception();
	}
	
	//Check if every equation in known has at least one variable
	for(int i = 0; i < knownSize; ++i) {
		if(known[i].getNumVariables() == 0) {
			std::cout << "Equations in known must have at least one graph in plainFlagAlgebra." << std::endl << std::endl;
			throw std::exception();
		}
	}
	
	//Make every equation in known type 1 (<=)
	/*for(int i = 0; i < knownSize; ++i) {
		if(known[i].getType() != 1) {
			std::vector<Graph> variablesTemp;
			
			for(int j = 0; j < known[i].getNumVariables(); ++j) {
				Graph graphTemp = known[i].getVariable(j);
				graphTemp.setCoefficient(-graphTemp.getCoefficient());
				variablesTemp.push_back(graphTemp);
			}
			
			Equation knownTemp(variablesTemp,zeros,-known[i].getAns(),1);
			
			known[i].setType(1);
			known.insert(known.begin()+i, knownTemp);
			++knownSize;
		}
	}*/

	int fN = f[0].getN();
	
	//Graphs in f can't have flags
	for(int i = 0; i < fSize; ++i) {
		if(f[i].getFlag().getN() != 0) {
			std::cout << "All graphs in f in plainFlagAlgebra must not have any flags." << std::endl << std::endl;
			throw std::exception();
		}
	}

	//Graphs in known can't have flags
	//Maybe fix?
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			if(known[i].getVariable(j).getFlag().getN() != 0) {
				std::cout << "All graphs in known in plainFlagAlgebra must not have any flags." << std::endl << std::endl;
				throw std::exception();
			}
		}
	}
	
	//Make sure everything has correct number of colors
	int numColors = f[0].getNumColors();
	
	for(int i = 1; i < fSize; ++i) {
		if(f[i].getNumColors() != numColors) {
			std::cout << "Everything in plainFlagAlgebra in f must have the same number of colors." << std::endl << std::endl;
			throw std::exception();
		}
	}
	
	for(int i = 0; i < zerosSize; ++i) {
		if(zeros[i].getNumColors() != numColors) {
			std::cout << "Everything in plainFlagAlgebra in zeros must have the same number of colors." << std::endl << std::endl;
			throw std::exception();
		}
	}
	
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			if(known[i].getVariable(j).getNumColors() != numColors) {
				std::cout << "Everything in plainFlagAlgebra in known must have the same number of colors." << std::endl << std::endl;
				throw std::exception();
			}
		}
	}
	
	//Not sure if strictly necessary, but it would probably give trash bounds otherwise
	if(fN > n) {
		std::cout << "Make n large enough so it has vertices at least as many vertices when multiplied by itself as n." << std::endl << std::endl;
		throw std::exception();
	}
	
	std::cout << std::endl;
	
	Equation fEq(f,zeros,Frac(1,1),0); //Type doesn't matter
	std::vector<Edge> edges {};
	Graph H(edges,1,numColors);
	Equation eq1({H},zeros,Frac(1,1),0);
	Equation eq2({H},zeros,Frac(1,1),0);
	
	std::cout << "Generating graphs to be used in resize." << std::endl;
	std::vector<Graph> allGraphs = generate(n,numColors,zeros);
	
	std::cout << std::endl << "Resizing f." << std::endl;
	Equation fEqResized = resize(fEq,allGraphs);
	f.clear();
	for(int i = 0; i < fEqResized.getNumVariables(); ++i) {
		f.push_back(fEqResized.getVariable(i));
	}
	std::cout << std::endl;
	
	fSize = f.size();
	
	//Resize known
	std::cout << "Resizing known." << std::endl;
	//#pragma omp parallel for
	for(int i = 0; i < knownSize; ++i) {
		std::cout << i+1 << " out of " << knownSize << std::endl;
		known[i] = resize(known[i],allGraphs);
	}
	std::cout << std::endl;

	std::vector<Equation> knownEq;
	std::vector<Equation> knownLess;
	
	for(int i = 0; i < knownSize; ++i) {
		if(known[i].getType() == 1) {
			knownLess.push_back(known[i]);
		}
		
		else {
			knownEq.push_back(known[i]);
		}	
	}
	
	std::ofstream myFile;
	myFile.open("Duals.txt");
	
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		myFile << "Label: " << i+1 << std::endl;
		allGraphs[i].printAdjMatToFile(myFile);
		myFile << std::endl;
	}
	
	//Used to give indices of A
	//Rather than multiplying we look at all graphs and work backwards to get coefficients
	//This works well as we have already generated all graphs of size n and multiplication basically just has us partially doing that every time so it saves time on generation
	//The only slow down could come from looking up indices but std::unordered_maps are SO fast it doesn't matter
	std::unordered_map<std::string, int> allGraphsMap;
					
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		allGraphsMap.insert(std::pair<std::string, int>(allGraphs[i].getCanonLabel(),i));
	}
	
	//All graphs of size n with flags such that parity issues work out
	std::cout << "Generating v." << std::endl;
	std::vector < std::vector < Graph > > v = generateV(n,numColors,zeros);
	std::cout << std::endl;
	std::cout << "Generating allGraphsWithFlags." << std::endl;
	std::vector < std::vector < Graph > > allGraphsWithFlags = generateAllGraphsWithFlags(n,numColors,zeros);
	std::cout << std::endl;

	std::vector<Frac> B; //Gives numbers to be printed
	std::vector< std::vector<Frac> > C; //From Known
	std::vector<Frac> C1;
	
	Frac zeroFrac(0,1);
	//Calculating B
	std::cout << "Calculating B." << std::endl << std::endl;
	B.resize(allGraphs.size(),zeroFrac);
	
	//TODO use map
	#pragma omp parallel for
	for(int i = 0; i < fSize; ++i) {
		for(int j = 0; j < (int)allGraphs.size(); ++j) {
			if(isomorphic(f[i],allGraphs[j])) {
				B[j] = f[i].getCoefficient();
				j = (int)allGraphs.size();
			}
		}
	}
	
	
	//First has 0 for ==, 1 for <=
	//Next entry is bound
	//Finally, it's the std::vector of coefficients in known (after it's been resized)
	
	//Calculating C
	std::cout << "Calculating C." << std::endl << std::endl;
	C.resize(allGraphs.size());
	C1.resize(knownSize,zeroFrac);
	
	#pragma omp parallel for
	for(int i = 0; i < allGraphs.size(); ++i) {
		C[i].resize(knownSize,zeroFrac);
	}
	
	//Note I'm changing the order of indices of C from when I had the Python Script
	#pragma omp parallel for
	for(int i = 0; i < knownSize; ++i) {
		C1[i] = known[i].getAns();
		
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			for(int k = 0; k < (int)allGraphs.size(); ++k) {
				if(isomorphic(known[i].getVariable(j),allGraphs[k])) {
					C[k][i] = known[i].getVariable(j).getCoefficient();
					k = (int)allGraphs.size();
				}
			}
		}
	}
	
	//Stupid C++ doesn't have hash for int[2]
	//std::vector< std::vector< std::unordered_map<int[2],Frac,boost::hash<int[2]> > > > newA;
	std::vector< std::vector< std::unordered_map<std::pair<int,int>,Frac,boost::hash<std::pair<int,int> > > > > newA;
	newA.resize(allGraphs.size());
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		newA[i].resize(v.size()); 
	}
	
	int FCOORDcounter = 0;
	
	for(int i = 0; i < (int)v.size(); ++i) { //i is first index of A
		std::cout << "In A, iteration " << i+1 << " out of " << v.size() << std::endl;  
		
		//Create hash map for third and fourth indices of A
		std::unordered_map<std::string, int> f;
		for(int j = 0; j < (int)v[i].size(); ++j) {
			f.insert(std::pair<std::string, int>(v[i][j].getCanonLabel(),j));
		}	
		
		int sizeOfFlag = v[i][0].getSizeOfFlag();
		
		#pragma omp parallel
		{
			
			#pragma omp for nowait schedule(dynamic) 
			for(int index2 = 0; index2 < (int)allGraphsWithFlags[i].size(); ++index2) {
				Graph G = allGraphsWithFlags[i][index2];
				
				//Make sure all flag vertices in G are first
				//TODO make a function of this
				std::vector<int> reordering(n);
				
				for(int j = 0; j < sizeOfFlag; ++j) {
					reordering[G.getFlagVertex(j)] = j;
				}
				int index = sizeOfFlag;
				
				if(sizeOfFlag != G.getSizeOfFlag()) {
					std::cout << "AHHHHHHHHHHHH" << std::endl;
				}
				
				for(int j = 0; j < n; ++j) {
					if(!G.isFlag(j)) {
						reordering[j] = index;
						++index;
					}
				}
				
				Frac coeff = G.getCoefficient();
				G = G.restriction(reordering);
				G.setCoefficient(coeff);
				
				Graph Gcopy = G;
				Gcopy.removeFlag();
			
				int b = allGraphsMap[Gcopy.getCanonLabel()];	
				Frac e = G.getCoefficient() * Frac(1,choose(n-sizeOfFlag-1,(n-sizeOfFlag)/2 - 1));
				
				std::vector<int> X((n-sizeOfFlag)/2 - 1); //X is zero indexed but it should actually be sizeOfFlag+1 indexed
				//Wlog first non-flag vertex in X
				for(int j = 0; j < (n-sizeOfFlag)/2 - 1; ++j) {
					X[j] = j;
				}
				
				do {
					std::vector<int> restriction1(n,-1);
					std::vector<int> restriction2(n,-1);
					
					for(int j = 0; j < sizeOfFlag; ++j) {
						restriction1[j] = j;
						restriction2[j] = j;
					}
					
					int counter1 = sizeOfFlag+1;
					int counter2 = sizeOfFlag;
					int Xcounter = 0;
					
					//Know flag vertices of G are first
					restriction1[sizeOfFlag] = sizeOfFlag;
					for(int j = sizeOfFlag + 1; j < n; ++j) {
						if((X.size() != 0) && (X[Xcounter] == (j - sizeOfFlag - 1))) {
							++Xcounter;
							restriction1[j] = counter1;
							++counter1; 
						}
						
						else {
							restriction2[j] = counter2;
							++counter2;
						}
					}
					
					Graph G1 = G.restriction(restriction1);
					Graph G2 = G.restriction(restriction2);
					
					int G1Index = f[G1.getCanonLabel()];
					int G2Index = f[G2.getCanonLabel()];
					
					int c = std::min(G1Index,G2Index);
					int d = std::max(G1Index,G2Index);
					
					#pragma omp critical
					{					
						auto it = newA[b][i].find({c,d});
						
						if( it == newA[b][i].end() ) {
							newA[b][i].insert({{c,d},e});
							++FCOORDcounter;
						}
						
						else {
							it->second = it->second + e; 
						}	
					}
				
				} while(nextSubset(X,n-sizeOfFlag-1,(n-sizeOfFlag)/2 - 1));
			}
		}
	}
	std::cout << std::endl;
	
	//return;

	std::cout << "Writing to file." << std::endl;

	//Make my own mosek output
	std::ofstream mosekFile("sdp.cbf");
	
	mosekFile << "VER\n3\n\n";
	
	mosekFile << "OBJSENSE\nMAX\n\n";
	
	mosekFile << "PSDVAR\n" << v.size() << "\n";
	for(int i = 0; i < v.size(); ++i) {
		mosekFile << v[i].size() << "\n";
	} 
	mosekFile << "\n";
	
	mosekFile << "VAR\n";
	if(knownLess.size() != 0) {
		mosekFile << knownSize + 1 << " " << 2 << "\n";
		mosekFile << "F " << knownEq.size()+1 << "\n";
		mosekFile << "L+ " << knownLess.size() << "\n\n";
	}
	else {
		mosekFile << knownSize + 1 << " 1\n";
		mosekFile << "F " << knownSize + 1 << "\n\n";
	}
	
	mosekFile << "CON\n";
	mosekFile << allGraphs.size() << " 1\n";
	mosekFile << "L- " << allGraphs.size() << "\n\n";
	
	mosekFile << "OBJACOORD\n";
	int OBJACOORDcounter = 1;
	long long int mult = 1;
	for(int i = 0; i < (int)known.size(); ++i) {
		if(known[i].getAns().getNum() != 0) {
			mult = lcm(mult,known[i].getAns().getDen());
			++OBJACOORDcounter;
		}
	} 
	
	mosekFile << OBJACOORDcounter << "\n";
	mosekFile << "0 " << mult << "\n";
	for(int i = 0; i < (int)knownEq.size(); ++i) {
		if(knownEq[i].getAns().getNum() != 0) {
			mosekFile << i+1 << " " << -knownEq[i].getAns().getNum() * (mult/knownEq[i].getAns().getDen()) << "\n" ;
		}
	}
	for(int i = 0; i < (int)knownLess.size(); ++i) {
		if(knownLess[i].getAns().getNum() != 0) {
			mosekFile << i+1+knownEq.size() << " " << -knownLess[i].getAns().getNum() * (mult/knownLess[i].getAns().getDen()) << "\n" ;
		}
	}
	mosekFile << "\n";
	
	std::cout << "Making Integer." << std::endl;
	
	std::vector<long long int> constraintMult;
	int BCOORDcounter = 0;
	
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		constraintMult.push_back(B[i].getDen());
		
		for(int j = 0; j < v.size(); ++j) {
			for(auto k : newA[i][j]) {
				constraintMult[i] = lcm(constraintMult[i],k.second.getDen());
			}
		}
		
		for(int j = 0; j < (int)known.size(); ++j) {
			for(int k = 0; k < known[j].getNumVariables(); ++k) {
				constraintMult[i] = lcm(constraintMult[i], known[j].getVariable(k).getCoefficient().getDen());
			}
		}
	}
	
	mosekFile << "FCOORD\n";
	mosekFile << FCOORDcounter << "\n";
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		std::cout << "In writing FCOORD iteration " << i + 1  << " out of " << allGraphs.size() << std::endl;
		for(int j = 0; j < (int)v.size(); ++j) {
			//TODO Probably faster to use auto for loop, but easier to debug this way
			for(int k = 0; k < (int)v[j].size(); ++k) {
				for(int l = k; l < (int)v[j].size(); ++l) {
					auto it = newA[i][j].find({k,l});
				
					if(it != newA[i][j].end()) {	
						if(k != l) {
							mosekFile << i << " " << j << " " << l << " " << k << " " << it->second.getNum() *(constraintMult[i] / (2.*it->second.getDen())) << "\n";
						} //May give fractions, but at most 1/2 so we don't have precision issue
						
						else {
							mosekFile << i << " " << j << " " << l << " " << k << " " << it->second.getNum() *(constraintMult[i] / (it->second.getDen())) << "\n";
						}
					}
				}
			}
		}
	}
	mosekFile <<"\n";
	
	mosekFile << "ACOORD\n";
	
	std::cout << "Writing ACOORD." << std::endl;
	
	int ACOORDcounter = 0;
	std::vector< std::array<int, 3> > ACOORDvec;
	
	#pragma omp parallel for
	for(int i = 0; i < (int)knownEq.size(); ++i) {
		for(int j = 0; j < knownEq[i].getNumVariables(); ++j) {
			for(int k = 0; k < (int)allGraphs.size(); ++k) {
				if(isomorphic(knownEq[i].getVariable(j),allGraphs[k])) {		
					#pragma omp critical 
					{
						++ACOORDcounter; 
						ACOORDvec.push_back({k,i+1, static_cast<int> ( knownEq[i].getVariable(j).getCoefficient().getNum() * (constraintMult[k] / knownEq[i].getVariable(j).getCoefficient().getDen()) )});
					}
					k = (int)allGraphs.size();
				}
			}
		}
	}
	
	#pragma omp parallel for
	for(int i = 0; i < (int)knownLess.size(); ++i) {
		for(int j = 0; j < knownLess[i].getNumVariables(); ++j) {
			for(int k = 0; k < (int)allGraphs.size(); ++k) {
				if(isomorphic(knownLess[i].getVariable(j),allGraphs[k])) {		
					#pragma omp critical 
					{
						++ACOORDcounter; 
						ACOORDvec.push_back({k,i+1+(int)knownEq.size(), static_cast<int> ( knownLess[i].getVariable(j).getCoefficient().getNum() * (constraintMult[k] / knownLess[i].getVariable(j).getCoefficient().getDen()))});
					}
					k = (int)allGraphs.size();
				}
			}
		}
	}
	
	mosekFile << ACOORDcounter + allGraphs.size() << "\n"; 	
	for(int i = 0; i < allGraphs.size(); ++i) {
		mosekFile << i << " 0 " << constraintMult[i] << "\n";
	}
	for(int i = 0; i < (int)ACOORDvec.size(); ++i) {
		mosekFile << ACOORDvec[i][0] << " " << ACOORDvec[i][1] << " " << -ACOORDvec[i][2] << "\n";
	}
	mosekFile << "\n";
	
	std::cout << "Writing BCOORD." << std::endl;
	
	mosekFile << "BCOORD\n";
	mosekFile << allGraphs.size() << "\n";
	
	if(maximize) {
		for(int i = 0; i < allGraphs.size(); ++i) {
			mosekFile << i << " " << -constraintMult[i] + B[i].getNum() * (constraintMult[i]/B[i].getDen()) << "\n";
		}
	}
	
	else {
		for(int i = 0; i < allGraphs.size(); ++i) {
			mosekFile << i << " " << -B[i].getNum() * (constraintMult[i]/B[i].getDen()) << "\n";
		}
	}
	
	//Flush Buffer
	mosekFile << std::endl;
	/*
	std::ofstream multFile("multiplication.txt");
	multFile << mult << "\n";
	for(int i = 0; i < allGraphs.size(); ++i) {
		multFile << constraintMult[i] << "\n";
	}
	
	multFile << std::endl;
	std::cout << std::endl;
	
	std::ofstream CSDPfile("CSDP.txt");
	
	CSDPfile << allGraphs.size() << "\n";
	
	CSDPfile << 2+v.size() << "\n";
	
	CSDPfile << "-" << (1+known.size()) << " ";
	for(int i = 0; i < v.size(); ++i) {
		CSDPfile << v[i].size() << " ";
	}
	CSDPfile << "\n";
	
	for(int i = 0; i < allGraphs.size(); ++i) {
		if(maximize) {
			CSDPfile <<-constraintMult[i] + B[i].getNum() * (constraintMult[i]/B[i].getDen()) << " ";
		}
	
		else {
			CSDPfile << B[i].getNum() * (constraintMult[i]/B[i].getDen()) << " ";
		}
	}
	CSDPfile << "\n";
	
	CSDPfile << "0 1 1 1 " << mult << "\n";
	
	for(int i = 0; i < known.size(); ++i) {
		CSDPfile << "0 1 " << i+2 << " " << i+2 << " " << -known[i].getAns().getNum() * (mult/known[i].getAns().getDen()) << "\n";
	}
	
	for(int i = 0; i < allGraphs.size(); ++i) {
		std::cout << "In writing to CSDP, iteration " << i << " out of " << allGraphs.size() << std::endl;
	
		CSDPfile << i+1 << " 1 1 1 " << constraintMult[i] << "\n";
		for(int j = 0; j < known.size(); ++j) {
			for(int k = 0; k < known[j].getNumVariables(); ++k) {
				if(isomorphic(known[j].getVariable(k),allGraphs[i])) {
					CSDPfile << i+1 << " 1 " << j+2 << " " << j+2 << " " <<  known[j].getVariable(k).getCoefficient().getNum() * (constraintMult[i] / known[j].getVariable(k).getCoefficient().getDen()) << "\n";
				}
			}
		}
		
		for(int j = 0; j < v.size(); ++j) {
			for(int k = 0; k < v[j].size(); ++k) {
				for(int l = k; l < v[j].size(); ++l) {
					auto it = newA[i][j].find({k,l});
					
					if(it != newA[i][j].end()) {
						if(k == l) {
					//Maybe divide by 2?
							CSDPfile << i+1 << " " << j+2 << " " << k+1 << " " << l+1 << " " << it->second.getNum() * (constraintMult[i] / (it->second.getDen())) << "\n";
						}
						
						else {
							CSDPfile << i+1 << " " << j+2 << " " << k+1 << " " << l+1 << " " << it->second.getNum() * (constraintMult[i] / (2.*it->second.getDen())) << "\n";
						}
					}
				}
			}
		}
	}
	
	CSDPfile << std::endl;
	std::cout << std::endl;
	*/
	std::cout << "Divide answers by " << mult << std::endl <<std::endl;
	
	return;
}

//Prints to plainFlagAlgerba1.txt & plainFlagAlgebra2.txt necessary files for python SDP code 
//f can be thought of as a linear combo of all graphs that we want to max/min
//Rather than taking a v this function takes the number of vertices to compute on (n)
//If using maximize, use 1-answer
void fastPlainFlagAlgebra(std::vector<Graph> &f, int n, std::vector<Graph> &zeros, std::vector<Equation> &known, bool maximize = true) {
	std::cout << "Starting plainFlagAlgebra." << std::endl;
	
	int fSize = f.size();
	int zerosSize = zeros.size();
	int knownSize = known.size();
		
	if(n <= 1) {
		std::cout << "Plain flag algebra method not set up for graphs with fewer than two vertices." << std::endl << std::endl;
		throw std::exception();
	}

	if(fSize == 0) {
		std::cout << "Need at least one graph in f in plainFlagAlgebra." << std::endl << std::endl;
		throw std::exception();
	}
	
	//Check if every equation in known has at least one variable
	for(int i = 0; i < knownSize; ++i) {
		if(known[i].getNumVariables() == 0) {
			std::cout << "Equations in known must have at least one graph in plainFlagAlgebra." << std::endl << std::endl;
			throw std::exception();
		}
	}

	int fN = f[0].getN();
	
	//Graphs in f can't have flags
	for(int i = 0; i < fSize; ++i) {
		if(f[i].getFlag().getN() != 0) {
			std::cout << "All graphs in f in plainFlagAlgebra must not have any flags." << std::endl << std::endl;
			throw std::exception();
		}
	}

	//Graphs in known can't have flags
	//Maybe fix?
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			if(known[i].getVariable(j).getFlag().getN() != 0) {
				std::cout << "All graphs in known in plainFlagAlgebra must not have any flags." << std::endl << std::endl;
				throw std::exception();
			}
		}
	}
	
	//Make sure everything has correct number of colors
	int numColors = f[0].getNumColors();
	
	for(int i = 1; i < fSize; ++i) {
		if(f[i].getNumColors() != numColors) {
			std::cout << "Everything in plainFlagAlgebra in f must have the same number of colors." << std::endl << std::endl;
			throw std::exception();
		}
	}
	
	for(int i = 0; i < zerosSize; ++i) {
		if(zeros[i].getNumColors() != numColors) {
			std::cout << "Everything in plainFlagAlgebra in zeros must have the same number of colors." << std::endl << std::endl;
			throw std::exception();
		}
	}
	
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			if(known[i].getVariable(j).getNumColors() != numColors) {
				std::cout << "Everything in plainFlagAlgebra in known must have the same number of colors." << std::endl << std::endl;
				throw std::exception();
			}
		}
	}
	
	//Not sure if strictly necessary, but it would probably give trash bounds otherwise
	if(fN > n) {
		std::cout << "Make n large enough so it has vertices at least as many vertices when multiplied by itself as n." << std::endl << std::endl;
		throw std::exception();
	}
	
	std::cout << std::endl;
	
	Equation fEq(f,zeros,Frac(1,1),0); //Type doesn't matter
	std::vector<Edge> edges {};
	Graph H(edges,1,numColors);
	Equation eq1({H},zeros,Frac(1,1),0);
	Equation eq2({H},zeros,Frac(1,1),0);
	
	std::cout << "Generating graphs to be used in resize." << std::endl;
	std::vector<Graph> allGraphs = generate(n,numColors,zeros);
	
	std::cout << std::endl << "Resizing f." << std::endl;
	Equation fEqResized = resize(fEq,allGraphs);
	f.clear();
	for(int i = 0; i < fEqResized.getNumVariables(); ++i) {
		f.push_back(fEqResized.getVariable(i));
	}
	std::cout << std::endl;
	
	fSize = f.size();
	
	//Resize known
	std::cout << "Resizing known." << std::endl;
	//#pragma omp parallel for
	for(int i = 0; i < knownSize; ++i) {
		std::cout << i+1 << " out of " << knownSize << std::endl;
		known[i] = resize(known[i],allGraphs);
	}
	std::cout << std::endl;

	std::vector<Equation> knownEq;
	std::vector<Equation> knownLess;
	
	for(int i = 0; i < knownSize; ++i) {
		if(known[i].getType() == 1) {
			knownLess.push_back(known[i]);
		}
		
		else {
			knownEq.push_back(known[i]);
		}	
	}
	
	std::ofstream myFile;
	myFile.open("Duals.txt");
	
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		myFile << "Label: " << i+1 << std::endl;
		allGraphs[i].printAdjMatToFile(myFile);
		myFile << std::endl;
	}
	
	//Used to give indices of A
	//Rather than multiplying we look at all graphs and work backwards to get coefficients
	//This works well as we have already generated all graphs of size n and multiplication basically just has us partially doing that every time so it saves time on generation
	//The only slow down could come from looking up indices but std::unordered_maps are SO fast it doesn't matter
	
	std::cout << "Generating first index map." << std::endl;
	std::unordered_map<std::string, int> firstIndexMap = firstIndex(n,numColors,zeros);
	std::cout << std::endl;
	
	std::cout << "Generating second index map." << std::endl;
	std::unordered_map<std::string, int> secondIndexMap;
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		secondIndexMap.insert(std::pair<std::string, int>(allGraphs[i].getCanonLabel(),i));
	}
	std::cout << std::endl;
	
	std::cout << "Generating last index map." << std::endl;
	std::unordered_map<std::string, int> lastIndexMap = lastIndex(n,numColors,zeros);
	std::cout << std::endl;

	std::vector<Frac> B; //Gives numbers to be printed
	std::vector< std::vector<Frac> > C; //From Known
	std::vector<Frac> C1;
	
	Frac zeroFrac(0,1);
	//Calculating B
	std::cout << "Calculating B." << std::endl << std::endl;
	B.resize(allGraphs.size(),zeroFrac);
	
	#pragma omp parallel for
	for(int i = 0; i < fSize; ++i) {
		B[secondIndexMap.at(f[i].getCanonLabel())] = f[i].getCoefficient();
	}
	
	
	//First has 0 for ==, 1 for <=
	//Next entry is bound
	//Finally, it's the std::vector of coefficients in known (after it's been resized)
	
	//Calculating C
	std::cout << "Calculating C." << std::endl << std::endl;
	C.resize(allGraphs.size());
	C1.resize(knownSize,zeroFrac);
	
	#pragma omp parallel for
	for(int i = 0; i < allGraphs.size(); ++i) {
		C[i].resize(knownSize,zeroFrac);
	}
	
	//#pragma omp parallel for
	for(int i = 0; i < knownSize; ++i) {
		C1[i] = known[i].getAns();
		
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			C[secondIndexMap.at(known[i].getVariable(j).getCanonLabel())][i] = known[i].getVariable(j).getCoefficient();
		}
	}
	
	//Stupid C++ doesn't have hash for int[2]
	//Use this because otherwise we use an absurd amount of memory - only small loss in speed b/c unordered
	//std::vector< std::vector< std::unordered_map<int[2],Frac,boost::hash<int[2]> > > > newA;
	std::vector< std::vector< std::unordered_map<std::pair<int,int>,Frac,boost::hash<std::pair<int,int> > > > > newA;
	newA.resize(allGraphs.size());
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		newA[i].resize(secondIndexMap.size()); 
	}
	
	int FCOORDcounter = 0;
	
	std::vector<int> dim(firstIndexMap.size(),0); 
	std::vector< std::unordered_set<std::string> > dimSet(firstIndexMap.size());
	
	#pragma omp parallel
	{
		#pragma omp for nowait schedule(dynamic) 
		for(int i = 0; i < allGraphs.size(); ++i) {
			Graph G = allGraphs[i];
			std::cout << "While generating A, iteration " << i+1 << " out of " << allGraphs.size() << std::endl;
			
			int b = secondIndexMap.at(G.getCanonLabel()); //Throws exception if out of range!
			
			for(int j = n/2; j <= n-1; ++j) {
				
				int sizeOfFlag = 2*j-n;
				
				if(sizeOfFlag > 0) {
					std::vector<int> X(sizeOfFlag); //flag
					for(int k = 0; k < sizeOfFlag; ++k) {
						X[k] = k;
					}
					
					Frac e = Frac(1,choose(n,j)*choose(j,sizeOfFlag));	
					
					do {
						std::vector<int> restriction1(n,-1); 
						 
						for(int k = 0; k < sizeOfFlag; ++k) {
							restriction1[X[k]] = k;
						}
						
						Graph G1 = G.restriction(restriction1); 
					
						std::vector<int> Y(j - sizeOfFlag); //Every that isn't the flag
						for(int k = 0; k < j- sizeOfFlag; ++k) {
							Y[k] = k;
						}
						
						do {
							//Need to adjust Y to be on [n] not [n-sizeOfFlag]
							std::vector<int> realY(j - sizeOfFlag);
							std::vector<int> adjusted(n-sizeOfFlag); //List of all non-flag vertices
							int tempIndex = 0;
							for(int k = 0; k < n; ++k) {
								bool isFlag = false;
								
								for(int l = 0; l < sizeOfFlag; ++l) {
									if(X[l] == k) {
										isFlag = true;
										l = sizeOfFlag;
									}
								}
								
								if(!isFlag) {
									adjusted[tempIndex] = k;
									++tempIndex;
								}
							}
							
							for(int k = 0; k < j - sizeOfFlag; ++k) {
								realY[k] = adjusted[Y[k]];
							}
						
							std::vector<int> restriction2(n,-1);
							std::vector<int> restriction3(n,-1);
							
							std::vector<int> restriction2c(n,-1);
							std::vector<int> restriction3c(n,-1);
							
							for(int k = 0; k < sizeOfFlag; ++k) {
								restriction3[X[k]] = k;
							}
							
							for(int k = 0; k < j - sizeOfFlag; ++k) {
								restriction2[realY[k]] = k;
								restriction3[realY[k]] = k+sizeOfFlag;
							}
							
							int tempIndex2 = 0; 
							int tempIndex3 = 0;
							for(int k = 0; k < n; ++k) {
								bool inY = false;
								for(int l = 0; l < j - sizeOfFlag; ++l) {
									if(realY[l] == k) {
										inY = true;
										l = j - sizeOfFlag;
									}
								}
								
								if(!inY) {
									restriction3c[k] = tempIndex3; //Double Check
									++tempIndex3; 
									
									bool inFlag = false;
									for(int l = 0; l < sizeOfFlag; ++l) {
										if(X[l] == k) {
											inFlag = true;
											l = sizeOfFlag;
										}
									}
									
									if(!inFlag) {
										restriction2c[k] = tempIndex2;
										++tempIndex2;
									}
									
								}
							}

							Graph G2 = G.restriction(restriction2);
							Graph G3 = G.restriction(restriction3);
							Graph G2c = G.restriction(restriction2c);
							Graph G3c = G.restriction(restriction3c);
							
							std::string lastIndexString1 = G1.getCanonLabel() + G2.getCanonLabel() + G3.getCanonLabel();
							std::string lastIndexString2 = G1.getCanonLabel() + G2c.getCanonLabel() + G3c.getCanonLabel();
							
							int a = firstIndexMap.at(G1.getCanonLabel());
					
							int index1 = lastIndexMap.at(lastIndexString1);
							int index2 = lastIndexMap.at(lastIndexString2);
							int c = std::min(index1,index2);
							int d = std::max(index1,index2);
							
							#pragma omp critical
							{					
								auto it = newA[i][a].find({c,d});
								
								if( it == newA[i][a].end() ) {
									newA[i][a].insert({{c,d},e});
									++FCOORDcounter;
								}
								
								else {
									it->second = it->second + e; 
								}	
								
								if(dimSet[a].find(lastIndexString1) == dimSet[a].end()) {
									++dim[a];
									dimSet[a].insert(lastIndexString1);
								}
								
								if(dimSet[a].find(lastIndexString2) == dimSet[a].end()) {
									++dim[a];
									dimSet[a].insert(lastIndexString2);
								}
							}
							
						} while(nextSubset(Y,n-sizeOfFlag,j-sizeOfFlag));
					} while(nextSubset(X,n,sizeOfFlag));
					
				}	
			}
		}
	}
	std::cout << std::endl;
	
	for(int i = 0; i < allGraphs.size(); ++i) {
		allGraphs[i].printEdges();
	}
	
	/*for(int i = 0; i < secondIndexMap.size(); ++i) {
		v[i][0].printEdges();
		v[i][0].printFlag();
	}*/
	
	//return;

	std::cout << "Writing to file." << std::endl;

	//Make my own mosek output
	std::ofstream mosekFile("sdp.cbf");
	
	mosekFile << "VER\n3\n\n";
	
	mosekFile << "OBJSENSE\nMAX\n\n";
	
	mosekFile << "PSDVAR\n" << firstIndexMap.size() << "\n";
	for(int i = 0; i < firstIndexMap.size(); ++i) {
		mosekFile << dim[i] << "\n";
	} 
	mosekFile << "\n";
	
	mosekFile << "VAR\n";
	if(knownLess.size() != 0) {
		mosekFile << knownSize + 1 << " " << 2 << "\n";
		mosekFile << "F " << knownEq.size()+1 << "\n";
		mosekFile << "L+ " << knownLess.size() << "\n\n";
	}
	else {
		mosekFile << knownSize + 1 << " 1\n";
		mosekFile << "F " << knownSize + 1 << "\n\n";
	}
	
	mosekFile << "CON\n";
	mosekFile << allGraphs.size() << " 1\n";
	mosekFile << "L- " << allGraphs.size() << "\n\n";
	
	mosekFile << "OBJACOORD\n";
	int OBJACOORDcounter = 1;
	long long int mult = 1;
	for(int i = 0; i < (int)known.size(); ++i) {
		if(known[i].getAns().getNum() != 0) {
			mult = lcm(mult,known[i].getAns().getDen());
			++OBJACOORDcounter;
		}
	} 
	
	mosekFile << OBJACOORDcounter << "\n";
	mosekFile << "0 " << mult << "\n";
	for(int i = 0; i < (int)knownEq.size(); ++i) {
		if(knownEq[i].getAns().getNum() != 0) {
			mosekFile << i+1 << " " << -knownEq[i].getAns().getNum() * (mult/knownEq[i].getAns().getDen()) << "\n" ;
		}
	}
	for(int i = 0; i < (int)knownLess.size(); ++i) {
		if(knownLess[i].getAns().getNum() != 0) {
			mosekFile << i+1+knownEq.size() << " " << -knownLess[i].getAns().getNum() * (mult/knownLess[i].getAns().getDen()) << "\n" ;
		}
	}
	mosekFile << "\n";
	
	std::cout << "Making Integer." << std::endl;
	
	std::vector<long long int> constraintMult;
	int BCOORDcounter = 0;
	
	
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		constraintMult.push_back(B[i].getDen());
	}
	
	#pragma omp parallel for
	for(int i = 0; i < (int)allGraphs.size(); ++i) {	
		for(int j = 0; j < firstIndexMap.size(); ++j) {
			for(auto k : newA[i][j]) {
				constraintMult[i] = lcm(constraintMult[i],k.second.getDen());
			}
		}
		
		for(int j = 0; j < (int)known.size(); ++j) {
			for(int k = 0; k < known[j].getNumVariables(); ++k) {
				constraintMult[i] = lcm(constraintMult[i], known[j].getVariable(k).getCoefficient().getDen());
			}
		}
	}
	
	
	
	mosekFile << "FCOORD\n";
	mosekFile << FCOORDcounter << "\n";
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		std::cout << "In writing FCOORD iteration " << i + 1  << " out of " << allGraphs.size() << std::endl;
		for(int j = 0; j < (int)firstIndexMap.size(); ++j) {
			//TODO Probably faster to use auto for loop, but easier to debug this way
			for(int k = 0; k < (int)dim[j]; ++k) {
				for(int l = k; l < (int)dim[j]; ++l) {
					auto it = newA[i][j].find({k,l});
				
					if(it != newA[i][j].end()) {	
						if(k != l) {
							mosekFile << i << " " << j << " " << l << " " << k << " " << it->second.getNum() *(constraintMult[i] / (2.*it->second.getDen())) << "\n";
						} //May give fractions, but at most 1/2 so we don't have precision issue
						
						else {
							mosekFile << i << " " << j << " " << l << " " << k << " " << it->second.getNum() *(constraintMult[i] / (it->second.getDen())) << "\n";
						}
					}
				}
			}
		}
	}
	mosekFile <<"\n";
	
	mosekFile << "ACOORD\n";
	
	std::cout << "Writing ACOORD." << std::endl;
	
	int ACOORDcounter = 0;
	std::vector< std::array<int, 3> > ACOORDvec;
	
	#pragma omp parallel for
	for(int i = 0; i < (int)knownEq.size(); ++i) {
		for(int j = 0; j < knownEq[i].getNumVariables(); ++j) {
			for(int k = 0; k < (int)allGraphs.size(); ++k) {
				if(isomorphic(knownEq[i].getVariable(j),allGraphs[k])) {		
					#pragma omp critical 
					{
						++ACOORDcounter; 
						ACOORDvec.push_back({k,i+1, static_cast<int> ( knownEq[i].getVariable(j).getCoefficient().getNum() * (constraintMult[k] / knownEq[i].getVariable(j).getCoefficient().getDen()) )});
					}
					k = (int)allGraphs.size();
				}
			}
		}
	}
	
	#pragma omp parallel for
	for(int i = 0; i < (int)knownLess.size(); ++i) {
		for(int j = 0; j < knownLess[i].getNumVariables(); ++j) {
			for(int k = 0; k < (int)allGraphs.size(); ++k) {
				if(isomorphic(knownLess[i].getVariable(j),allGraphs[k])) {		
					#pragma omp critical 
					{
						++ACOORDcounter; 
						ACOORDvec.push_back({k,i+1+(int)knownEq.size(), static_cast<int> ( knownLess[i].getVariable(j).getCoefficient().getNum() * (constraintMult[k] / knownLess[i].getVariable(j).getCoefficient().getDen()))});
					}
					k = (int)allGraphs.size();
				}
			}
		}
	}
	
	mosekFile << ACOORDcounter + allGraphs.size() << "\n"; 	
	for(int i = 0; i < allGraphs.size(); ++i) {
		mosekFile << i << " 0 " << constraintMult[i] << "\n";
	}
	for(int i = 0; i < (int)ACOORDvec.size(); ++i) {
		mosekFile << ACOORDvec[i][0] << " " << ACOORDvec[i][1] << " " << -ACOORDvec[i][2] << "\n";
	}
	mosekFile << "\n";
	
	std::cout << "Writing BCOORD." << std::endl;
	
	mosekFile << "BCOORD\n";
	mosekFile << allGraphs.size() << "\n";
	
	if(maximize) {
		for(int i = 0; i < allGraphs.size(); ++i) {
			mosekFile << i << " " << -constraintMult[i] + B[i].getNum() * (constraintMult[i]/B[i].getDen()) << "\n";
		}
	}
	
	else {
		for(int i = 0; i < allGraphs.size(); ++i) {
			mosekFile << i << " " << -B[i].getNum() * (constraintMult[i]/B[i].getDen()) << "\n";
		}
	}
	
	//Flush Buffer
	mosekFile << std::endl;
	/*
	std::ofstream multFile("multiplication.txt");
	multFile << mult << "\n";
	for(int i = 0; i < allGraphs.size(); ++i) {
		multFile << constraintMult[i] << "\n";
	}
	
	multFile << std::endl;
	std::cout << std::endl;
	
	std::ofstream CSDPfile("CSDP.txt");
	
	CSDPfile << allGraphs.size() << "\n";
	
	CSDPfile << 2+firstIndexMap.size() << "\n";
	
	CSDPfile << "-" << (1+known.size()) << " ";
	for(int i = 0; i < firstIndexMap.size(); ++i) {
		CSDPfile << dim[i] << " ";
	}
	CSDPfile << "\n";
	
	for(int i = 0; i < allGraphs.size(); ++i) {
		if(maximize) {
			CSDPfile <<-constraintMult[i] + B[i].getNum() * (constraintMult[i]/B[i].getDen()) << " ";
		}
	
		else {
			CSDPfile << B[i].getNum() * (constraintMult[i]/B[i].getDen()) << " ";
		}
	}
	CSDPfile << "\n";
	
	CSDPfile << "0 1 1 1 " << mult << "\n";
	
	for(int i = 0; i < known.size(); ++i) {
		CSDPfile << "0 1 " << i+2 << " " << i+2 << " " << -known[i].getAns().getNum() * (mult/known[i].getAns().getDen()) << "\n";
	}
	
	for(int i = 0; i < allGraphs.size(); ++i) {
		std::cout << "In writing to CSDP, iteration " << i << " out of " << allGraphs.size() << std::endl;
	
		CSDPfile << i+1 << " 1 1 1 " << constraintMult[i] << "\n";
		for(int j = 0; j < known.size(); ++j) {
			for(int k = 0; k < known[j].getNumVariables(); ++k) {
				if(isomorphic(known[j].getVariable(k),allGraphs[i])) {
					CSDPfile << i+1 << " 1 " << j+2 << " " << j+2 << " " <<  known[j].getVariable(k).getCoefficient().getNum() * (constraintMult[i] / known[j].getVariable(k).getCoefficient().getDen()) << "\n";
				}
			}
		}
		
		for(int j = 0; j < firstIndexMap.size(); ++j) {
			for(int k = 0; k < dim[j]; ++k) {
				for(int l = k; l < dim[j]; ++l) {
					auto it = newA[i][j].find({k,l});
					
					if(it != newA[i][j].end()) {
						if(k == l) {
					//Maybe divide by 2?
							CSDPfile << i+1 << " " << j+2 << " " << k+1 << " " << l+1 << " " << it->second.getNum() * (constraintMult[i] / (it->second.getDen())) << "\n";
						}
						
						else {
							CSDPfile << i+1 << " " << j+2 << " " << k+1 << " " << l+1 << " " << it->second.getNum() * (constraintMult[i] / (2.*it->second.getDen())) << "\n";
						}
					}
				}
			}
		}
	}
	
	CSDPfile << std::endl;
	std::cout << std::endl;
	*/
	std::cout << "Divide answers by " << mult << std::endl <<std::endl;
	
	return;
	return;
}




















void plainFlagAlgebraApprox(std::vector<Graph> &f, int n, int r, std::vector<Graph> &zeros, std::vector<Equation> &known, bool maximize = true) {
	std::cout << "Starting plainFlagAlgebra." << std::endl;
	
	//The way the python script is setup we need at least one known, this just says that edge density is <= 1
	//if(known.size() == 0) {
		//Graph K2({{}},2,f[0].getNumColors());
		//Equation knownTemp({K2},zeros,Frac(1,1),1);
		//known.push_back(knownTemp);
	//}
	
	int fSize = f.size();
	int zerosSize = zeros.size();
	int knownSize = known.size();
		
	if(n <= 1) {
		std::cout << "Plain flag algebra method not set up for graphs with fewer than two vertices." << std::endl << std::endl;
		throw std::exception();
	}

	if(fSize == 0) {
		std::cout << "Need at least one graph in f in plainFlagAlgebra." << std::endl << std::endl;
		throw std::exception();
	}
	
	//Check if every equation in known has at least one variable
	for(int i = 0; i < knownSize; ++i) {
		if(known[i].getNumVariables() == 0) {
			std::cout << "Equations in known must have at least one graph in plainFlagAlgebra." << std::endl << std::endl;
			throw std::exception();
		}
	}
	
	//Make every equation in known type 1 (<=)
	/*for(int i = 0; i < knownSize; ++i) {
		if(known[i].getType() != 1) {
			std::vector<Graph> variablesTemp;
			
			for(int j = 0; j < known[i].getNumVariables(); ++j) {
				Graph graphTemp = known[i].getVariable(j);
				graphTemp.setCoefficient(-graphTemp.getCoefficient());
				variablesTemp.push_back(graphTemp);
			}
			
			Equation knownTemp(variablesTemp,zeros,-known[i].getAns(),1);
			
			known[i].setType(1);
			known.insert(known.begin()+i, knownTemp);
			++knownSize;
		}
	}*/

	int fN = f[0].getN();
	
	//Graphs in f can't have flags
	for(int i = 0; i < fSize; ++i) {
		if(f[i].getFlag().getN() != 0) {
			std::cout << "All graphs in f in plainFlagAlgebra must not have any flags." << std::endl << std::endl;
			throw std::exception();
		}
	}

	//Graphs in known can't have flags
	//Maybe fix?
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			if(known[i].getVariable(j).getFlag().getN() != 0) {
				std::cout << "All graphs in known in plainFlagAlgebra must not have any flags." << std::endl << std::endl;
				throw std::exception();
			}
		}
	}
	
	//Make sure everything has correct number of colors
	int numColors = f[0].getNumColors();
	
	for(int i = 1; i < fSize; ++i) {
		if(f[i].getNumColors() != numColors) {
			std::cout << "Everything in plainFlagAlgebra in f must have the same number of colors." << std::endl << std::endl;
			throw std::exception();
		}
	}
	
	for(int i = 0; i < zerosSize; ++i) {
		if(zeros[i].getNumColors() != numColors) {
			std::cout << "Everything in plainFlagAlgebra in zeros must have the same number of colors." << std::endl << std::endl;
			throw std::exception();
		}
	}
	
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			if(known[i].getVariable(j).getNumColors() != numColors) {
				std::cout << "Everything in plainFlagAlgebra in known must have the same number of colors." << std::endl << std::endl;
				throw std::exception();
			}
		}
	}
	
	//Not sure if strictly necessary, but it would probably give trash bounds otherwise
	if(fN > n) {
		std::cout << "Make n large enough so it has vertices at least as many vertices when multiplied by itself as n." << std::endl << std::endl;
		throw std::exception();
	}
	
	std::cout << std::endl;
	
	Equation fEq(f,zeros,Frac(1,1),0); //Type doesn't matter
	std::vector<Edge> edges {};
	Graph H(edges,1,numColors);
	Equation eq1({H},zeros,Frac(1,1),0);
	Equation eq2({H},zeros,Frac(1,1),0);
	
	std::cout << "Generating graphs to be used in resize." << std::endl;
	std::vector<Graph> allGraphs = generate(n,numColors,zeros);
	
	std::cout << std::endl << "Resizing f." << std::endl;
	Equation fEqResized = resize(fEq,allGraphs);
	f.clear();
	for(int i = 0; i < fEqResized.getNumVariables(); ++i) {
		f.push_back(fEqResized.getVariable(i));
	}
	std::cout << std::endl;
	
	fSize = f.size();
	
	//Resize known
	std::cout << "Resizing known." << std::endl;
	//#pragma omp parallel for
	for(int i = 0; i < knownSize; ++i) {
		std::cout << i+1 << " out of " << knownSize << std::endl;
		known[i] = resize(known[i],allGraphs);
	}
	std::cout << std::endl;

	std::vector<Equation> knownEq;
	std::vector<Equation> knownLess;
	
	for(int i = 0; i < knownSize; ++i) {
		if(known[i].getType() == 1) {
			knownLess.push_back(known[i]);
		}
		
		else {
			knownEq.push_back(known[i]);
		}	
	}
	
	std::ofstream myFile;
	myFile.open("Duals.txt");
	
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		myFile << "Label: " << i+1 << std::endl;
		allGraphs[i].printAdjMatToFile(myFile);
		myFile << std::endl;
	}
	
	//Used to give indices of A
	//Rather than multiplying we look at all graphs and work backwards to get coefficients
	//This works well as we have already generated all graphs of size n and multiplication basically just has us partially doing that every time so it saves time on generation
	//The only slow down could come from looking up indices but std::unordered_maps are SO fast it doesn't matter
	std::unordered_map<std::string, int> allGraphsMap;
					
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		allGraphsMap.insert(std::pair<std::string, int>(allGraphs[i].getCanonLabel(),i));
	}
	
	//All graphs of size n with flags such that parity issues work out
	std::cout << "Generating v." << std::endl;
	std::vector < std::vector < Graph > > v = generateV(n,numColors,zeros);
	std::cout << std::endl;
	std::cout << "Generating allGraphsWithFlags." << std::endl;
	std::vector < std::vector < Graph > > allGraphsWithFlags = generateAllGraphsWithFlags(n,numColors,zeros);
	std::cout << std::endl;

	std::vector<Frac> B; //Gives numbers to be printed
	std::vector< std::vector<Frac> > C; //From Known
	std::vector<Frac> C1;
	
	Frac zeroFrac(0,1);
	//Calculating B
	std::cout << "Calculating B." << std::endl << std::endl;
	B.resize(allGraphs.size(),zeroFrac);
	
	//TODO use map
	#pragma omp parallel for
	for(int i = 0; i < fSize; ++i) {
		for(int j = 0; j < (int)allGraphs.size(); ++j) {
			if(isomorphic(f[i],allGraphs[j])) {
				B[j] = f[i].getCoefficient();
				j = (int)allGraphs.size();
			}
		}
	}
	
	
	//First has 0 for ==, 1 for <=
	//Next entry is bound
	//Finally, it's the std::vector of coefficients in known (after it's been resized)
	
	//Calculating C
	std::cout << "Calculating C." << std::endl << std::endl;
	C.resize(allGraphs.size());
	C1.resize(knownSize,zeroFrac);
	
	#pragma omp parallel for
	for(int i = 0; i < allGraphs.size(); ++i) {
		C[i].resize(knownSize,zeroFrac);
	}
	
	//Note I'm changing the order of indices of C from when I had the Python Script
	#pragma omp parallel for
	for(int i = 0; i < knownSize; ++i) {
		C1[i] = known[i].getAns();
		
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			for(int k = 0; k < (int)allGraphs.size(); ++k) {
				if(isomorphic(known[i].getVariable(j),allGraphs[k])) {
					C[k][i] = known[i].getVariable(j).getCoefficient();
					k = (int)allGraphs.size();
				}
			}
		}
	}
	
	//Stupid C++ doesn't have hash for int[2]
	//std::vector< std::vector< std::unordered_map<int[2],Frac,boost::hash<int[2]> > > > newA;
	std::vector< std::vector< std::unordered_map<std::pair<int,int>,Frac,boost::hash<std::pair<int,int> > > > > newA;
	newA.resize(allGraphs.size());
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		newA[i].resize(v.size()); 
	}
	
	int FCOORDcounter = 0;
	
	for(int i = 0; i < (int)v.size(); ++i) { //i is first index of A
		std::cout << "In A, iteration " << i+1 << " out of " << v.size() << std::endl;  
		
		//Create hash map for third and fourth indices of A
		std::unordered_map<std::string, int> f;
		for(int j = 0; j < (int)v[i].size(); ++j) {
			f.insert(std::pair<std::string, int>(v[i][j].getCanonLabel(),j));
		}	
		
		int sizeOfFlag = v[i][0].getSizeOfFlag();
		
		#pragma omp parallel
		{
			
			#pragma omp for nowait schedule(dynamic) 
			for(int index2 = 0; index2 < (int)allGraphsWithFlags[i].size(); ++index2) {
				Graph G = allGraphsWithFlags[i][index2];
				
				//Make sure all flag vertices in G are first
				//TODO make a function of this
				std::vector<int> reordering(n);
				
				for(int j = 0; j < sizeOfFlag; ++j) {
					reordering[G.getFlagVertex(j)] = j;
				}
				int index = sizeOfFlag;
				
				for(int j = 0; j < n; ++j) {
					if(!G.isFlag(j)) {
						reordering[j] = index;
						++index;
					}
				}
				
				Frac coeff = G.getCoefficient();
				G = G.restriction(reordering);
				G.setCoefficient(coeff);
				
				Graph Gcopy = G;
				Gcopy.removeFlag();
			
				int b = allGraphsMap[Gcopy.getCanonLabel()];	
				Frac e = G.getCoefficient() * Frac(1,choose(n-sizeOfFlag-1,(n-sizeOfFlag)/2 - 1));
				
				std::vector<int> X((n-sizeOfFlag)/2 - 1); //X is zero indexed but it should actually be sizeOfFlag+1 indexed
				//Wlog first non-flag vertex in X
				for(int j = 0; j < (n-sizeOfFlag)/2 - 1; ++j) {
					X[j] = j;
				}
				
				do {
					std::vector<int> restriction1(n,-1);
					std::vector<int> restriction2(n,-1);
					
					for(int j = 0; j < sizeOfFlag; ++j) {
						restriction1[j] = j;
						restriction2[j] = j;
					}
					
					int counter1 = sizeOfFlag+1;
					int counter2 = sizeOfFlag;
					int Xcounter = 0;
					
					//Know flag vertices of G are first
					restriction1[sizeOfFlag] = sizeOfFlag;
					for(int j = sizeOfFlag + 1; j < n; ++j) {
						if((X.size() != 0) && (X[Xcounter] == (j - sizeOfFlag - 1))) {
							++Xcounter;
							restriction1[j] = counter1;
							++counter1; 
						}
						
						else {
							restriction2[j] = counter2;
							++counter2;
						}
					}
					
					Graph G1 = G.restriction(restriction1);
					Graph G2 = G.restriction(restriction2);
					
					int G1Index = f[G1.getCanonLabel()];
					int G2Index = f[G2.getCanonLabel()];
					
					int c = std::min(G1Index,G2Index);
					int d = std::max(G1Index,G2Index);
					
					#pragma omp critical
					{					
						auto it = newA[b][i].find({c,d});
						
						if( it == newA[b][i].end() ) {
							newA[b][i].insert({{c,d},e});
							++FCOORDcounter;
						}
						
						else {
							it->second = it->second + e; 
						}	
					}
				
				} while(nextSubset(X,n-sizeOfFlag-1,(n-sizeOfFlag)/2 - 1));
			}
		}
	}
	std::cout << std::endl;
	
	for(int i = 0; i < v.size(); ++i) {
		v[i][0].printEdges();
		v[i][0].printFlag();
	}
	
	//return;
	
	std::cout << "Generating phi." << std::endl << std::endl;
	std::vector< std::vector < std::vector<int> > > phi; //Almost certainly not the best data structure
	std::vector<int> phiInv1;
	std::vector<int> phiInv2;
	std::vector<int> phiInv3;
	int phiSize = 0;
	
	phi.resize(v.size());
	
	for(int j = 0; j < v.size(); ++j) {
		phi[j].resize(v[j].size());
		for(int k = 0; k < v[j].size(); ++k) {
			phi[j][k].resize(v[j].size());
			for(int l = k; l < v[j].size(); ++l) {
				phi[j][k][l] = phiSize;
				phiInv1.push_back(j);
				phiInv2.push_back(k);
				phiInv3.push_back(l);
				++phiSize;
			}
		}
	}

	//This has a bunhc of multiplying by 2 that the other one doesn't because of how Mosek deals with PSD matrices

	std::cout << "Writing to file." << std::endl;

	//Make my own mosek output
	std::ofstream mosekFile("sdp.cbf");
	
	mosekFile << "VER\n3\n\n";
	
	mosekFile << "OBJSENSE\nMAX\n\n";
	
	mosekFile << "VAR\n";
	if(knownLess.size() != 0) {
		mosekFile << knownSize + 1 + phiSize << " " << 2 << "\n";
		mosekFile << "F " << knownEq.size()+1+phiSize << "\n";
		mosekFile << "L+ " << knownLess.size() << "\n\n";
	}
	else {
		mosekFile << knownSize + 1 + phiSize << " 1\n";
		mosekFile << "F " << knownSize + 1 + phiSize << "\n\n";
	}
	
	mosekFile << "CON\n";
	//Generate number of new constraints
	int newConstraints = 0;
	for(int j = 0; j < v.size(); ++j) {
		newConstraints += choose(v[j].size()+r+1, v[j].size()-1);
	}
	mosekFile << allGraphs.size() + newConstraints << " 1\n";
	mosekFile << "L- " << allGraphs.size() + newConstraints << "\n\n";
	
	mosekFile << "OBJACOORD\n";
	int OBJACOORDcounter = 1;
	long long int mult = 1;
	for(int i = 0; i < (int)known.size(); ++i) {
		if(known[i].getAns().getNum() != 0) {
			mult = lcm(mult,known[i].getAns().getDen());
			++OBJACOORDcounter;
		}
	} 
	
	mosekFile << OBJACOORDcounter << "\n";
	mosekFile << "0 " << mult << "\n";
	for(int i = 0; i < (int)knownEq.size(); ++i) {
		if(knownEq[i].getAns().getNum() != 0) {
			mosekFile << i+1 << " " << -knownEq[i].getAns().getNum() * (mult/knownEq[i].getAns().getDen()) << "\n" ;
		}
	}
	for(int i = 0; i < (int)knownLess.size(); ++i) {
		if(knownLess[i].getAns().getNum() != 0) {
			mosekFile << i+1+knownEq.size() << " " << -knownLess[i].getAns().getNum() * (mult/knownLess[i].getAns().getDen()) << "\n" ;
		}
	}
	mosekFile << "\n";

	
	std::cout << "Making Integer." << std::endl;
	
	std::vector<long long int> constraintMult;
	int BCOORDcounter = 0;
	
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		constraintMult.push_back(B[i].getDen());
		
		for(int j = 0; j < v.size(); ++j) {
			for(auto k : newA[i][j]) {
				constraintMult[i] = lcm(constraintMult[i],k.second.getDen());
			}
		}
		
		for(int j = 0; j < (int)known.size(); ++j) {
			for(int k = 0; k < known[j].getNumVariables(); ++k) {
				constraintMult[i] = lcm(constraintMult[i], known[j].getVariable(k).getCoefficient().getDen());
			}
		}
	}
	
	/*mosekFile << "FCOORD\n";
	mosekFile << FCOORDcounter << "\n";
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		std::cout << "In writing FCOORD iteration " << i + 1  << " out of " << allGraphs.size() << std::endl;
		for(int j = 0; j < (int)v.size(); ++j) {
			//TODO Probably faster to use auto for loop, but easier to debug this way
			for(int k = 0; k < (int)v[j].size(); ++k) {
				for(int l = k; l < (int)v[j].size(); ++l) {
					auto it = newA[i][j].find({k,l});
				
					if(it != newA[i][j].end()) {	
						if(k != l) {
							mosekFile << i << " " << j << " " << l << " " << k << " " << it->second.getNum() *(constraintMult[i] / (2.*it->second.getDen())) << "\n";
						} //May give fractions, but at most 1/2 so we don't have precision issue
						
						else {
							mosekFile << i << " " << j << " " << l << " " << k << " " << it->second.getNum() *(constraintMult[i] / (it->second.getDen())) << "\n";
						}
					}
				}
			}
		}
	}
	mosekFile <<"\n";*/
	
	mosekFile << "ACOORD\n";
	
	std::cout << "Writing ACOORD." << std::endl;
	
	int ACOORDcounter = 0;
	std::vector< std::array<int, 3> > ACOORDvec;
	
	
	
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		for(int j = 0; j < phiSize; ++j) {
			auto it = newA[i][phiInv1[j]].find({phiInv2[j],phiInv3[j]});
			if(it != newA[i][phiInv1[j]].end()) {
				++ACOORDcounter;
				if(phiInv2[j] == phiInv3[j]) {
					ACOORDvec.push_back({i,1+j,static_cast<int>(it->second.getNum() *(constraintMult[i] / static_cast<int>(it->second.getDen())))});
				}
				else {
					ACOORDvec.push_back({i,1+j,static_cast<int>(it->second.getNum() *(constraintMult[i] / static_cast<int>(it->second.getDen())))});
				}
			}
		}
	}

	int firstIndex = (int)allGraphs.size();
	for(int j = 0; j < v.size(); ++j) {
		std::cout << "In writing ACOORD iteration " << j << " out of " << v.size() << std::endl;
		std::vector<int> sum;
		sum.push_back(r+2);
		
		for(int k = 1; k < v[j].size(); ++k) {
			sum.push_back(0);
		}
		
		do {
			
			for(int k = 0; k < v[j].size(); ++k) {
				for(int l = k; l < v[j].size(); ++l) {
					if((k != l) && (sum[k] > 0) && (sum[l] > 0)) {
						++ACOORDcounter;
						
						std::vector<int> multinomialInput = sum;
						--multinomialInput[k];
						--multinomialInput[l];
						ACOORDvec.push_back({firstIndex, 1+phi[j][k][l], 2*-(int)multinomial(r,multinomialInput)});
					}
				
					else if((k == l) && (sum[k] > 1)) {
						++ACOORDcounter;
							
						std::vector<int> multinomialInput = sum;
						--multinomialInput[k];
					--multinomialInput[k];
						ACOORDvec.push_back({firstIndex, 1+phi[j][k][k], -(int)multinomial(r,multinomialInput)});
					}
				}
			}
			++firstIndex;
		} while(nextOrderedSum(sum));
	} 
	
	#pragma omp parallel for
	for(int i = 0; i < (int)knownEq.size(); ++i) {
		for(int j = 0; j < knownEq[i].getNumVariables(); ++j) {
			for(int k = 0; k < (int)allGraphs.size(); ++k) {
				if(isomorphic(knownEq[i].getVariable(j),allGraphs[k])) {		
					#pragma omp critical 
					{
						++ACOORDcounter; 
						ACOORDvec.push_back({k,i+1+phiSize, -static_cast<int> ( knownEq[i].getVariable(j).getCoefficient().getNum() * (constraintMult[k] / knownEq[i].getVariable(j).getCoefficient().getDen()) )});
					}
					k = (int)allGraphs.size();
				}
			}
		}
	}
	
	#pragma omp parallel for
	for(int i = 0; i < (int)knownLess.size(); ++i) {
		for(int j = 0; j < knownLess[i].getNumVariables(); ++j) {
			for(int k = 0; k < (int)allGraphs.size(); ++k) {
				if(isomorphic(knownLess[i].getVariable(j),allGraphs[k])) {		
					#pragma omp critical 
					{
						++ACOORDcounter; 
						ACOORDvec.push_back({k,i+1+(int)knownEq.size()+phiSize, -static_cast<int> ( knownLess[i].getVariable(j).getCoefficient().getNum() * (constraintMult[k] / knownLess[i].getVariable(j).getCoefficient().getDen()))});
					}
					k = (int)allGraphs.size();
				}
			}
		}
	}
	
	mosekFile << ACOORDcounter + allGraphs.size() << "\n"; 	
	for(int i = 0; i < allGraphs.size(); ++i) {
		mosekFile << i << " 0 " << constraintMult[i] << "\n";
	}
	for(int i = 0; i < (int)ACOORDvec.size(); ++i) {
		mosekFile << ACOORDvec[i][0] << " " << ACOORDvec[i][1] << " " << ACOORDvec[i][2] << "\n";
	}
	mosekFile << "\n";
	
	std::cout << "Writing BCOORD." << std::endl;
	
	mosekFile << "BCOORD\n";
	mosekFile << allGraphs.size() + newConstraints << "\n";
	
	if(maximize) {
		for(int i = 0; i < allGraphs.size(); ++i) {
			mosekFile << i << " " << -constraintMult[i] + B[i].getNum() * (constraintMult[i]/B[i].getDen()) << "\n";
		}
	}
	
	else {
		for(int i = 0; i < allGraphs.size(); ++i) {
			mosekFile << i << " " << -B[i].getNum() * (constraintMult[i]/B[i].getDen()) << "\n";
		}
	}
	
	for(int i = 0; i < newConstraints; ++i) {
		mosekFile << i + (int)allGraphs.size() << " 0\n";
	}
	
	//Flush Buffer
	mosekFile << std::endl;
	
	std::ofstream multFile("multiplication.txt");
	multFile << mult << "\n";
	for(int i = 0; i < allGraphs.size(); ++i) {
		multFile << constraintMult[i] << "\n";
	}
	
	multFile << std::endl;
	std::cout << std::endl;
	
	std::cout << "Divide answers by " << mult << std::endl <<std::endl;
	
	return;
}


