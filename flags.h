//Nauty
#define MAXN 1000    /* Define this before including nauty.h */
extern "C" {
	#include "nauty.h"   
	#include "naututil.h"
	#include "gtools.h"
}

const vector<long long int> numberOfGraphs = {1,1,2,4,11, 34, 156, 1044, 12346, 274668, 12005168, 1018997864, 165091172592}; //Up to 12 vertices

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
		vector<vector<int> > adjMat; //Adjacency Matrix
		vector<int> flag; //The flag can be determined by a subset of edges
		string canonLabel;
		vector<vector<int> > orbits; //orbits in automorphsim group (from nauty) orbits[i] is the ith orbit of vertices, if vertex is a flag it is in its own orbit (not iff)
		int numOrbits = 0;

	public:
		
		//-------------------------------
		//-----Constructor for Graph-----
		//-------------------------------
		
		//Need to specify NUMCOLORS, because it's possible that it differs
		//from actual number of colors
		Graph(const vector<Edge> &edgeList, const int N, const int NUMCOLORS) {
			numColors = NUMCOLORS;
			n = N;
			
			if(numColors < 0) {
				cout << "In constructor need a non-negative number of colors" << endl << endl;
				throw exception();
			}
			
			if(n < 0) {
				cout << "In constructor need a non-negative number of vertices" << endl << endl;
				throw exception();
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
					cout << "In constructor too few colors specified, or range is wrong, note that non-edges count as a color, and colors are zero indexed"  << endl << endl;
					throw exception();
				}
				
				if(edge.a >= n) {
					cout << "In constructor an edge has a vertex outside of the range." << endl << endl;
					throw exception();
				}
				
				if(edge.b >= n) {
					cout << "In constructor an edge has a vertex outside of the range." << endl << endl;
					throw exception();
				}
			
            adjMat[edge.a][edge.b] = edge.color;
            adjMat[edge.b][edge.a] = edge.color;
  			}
  			
  			Frac frac(1,1);
  			coefficient = frac;
  			
  			canonRelabel();
		}

		
		//---------------
		//-----Clone-----
		//---------------
		
		Graph clone() const{
			vector<Edge> edges;
		
			for(int i = 0; i < n-1; ++i) {
				for(int j = i+1; j < n; ++j) {
					edges.push_back({i,j,adjMat[i][j]});
				}
			}
			
			Graph G(edges,n,numColors);
			
			G.setCoefficient(coefficient);
			G.setFlag(flag);
			
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
		Graph restriction(const vector<int> &sigma) const{
			if(int(sigma.size()) != n) {
				cout << "In restriction, sigma size wrong." << endl;
				cout << "Sigma size = " << sigma.size() << endl;
				cout << "n = " << n << endl << endl;
				throw exception();
			}		
					
			//Check that it's a subset
			int k = 0;
					
			for(int i = 0; i < n; ++i) {
				if(sigma[i] != -1) {
				++k;
					for(int j = i+1; j < n; ++j) {
						if(sigma[i] == sigma[j]) {
							cout << "In restriction there is a duplicate in the input." << endl;
							cout << "Indices : " << i << ", " << j << endl;
							cout << "Value: " << sigma[i] << endl << endl;
							throw exception();
						}
					}
				}
			}
			
			//Check if sigma has elements outside of range
			for(int i = 0; i < n; ++i) {
				if((sigma[i] >= k) || (sigma[i] < -1)) {
					cout << "In restriction an element of sigma falls outside of range of vertices." << endl;
					cout << "Index: " << i << endl;
					cout << "Value: " << sigma[i] << endl << endl;
					throw exception();
				}
			}
			
			//Edges
			vector<Edge> newEdges;
			for(int i = 0; i < n-1; ++i) {
				for(int j = i+1; j < n; ++j) {
					if((adjMat[i][j] != 0) && (sigma[i] != -1) && (sigma[j] != -1)) {
						newEdges.push_back({sigma[i],sigma[j],adjMat[i][j]});
					}
				}
			}
			
			Graph G(newEdges, k, numColors);
			
			//Flag
			/*vector<int> newFlag;
			for(int i = 0; i < n; ++i) {
				if(sigma[i] != -1) {
					for(int j = 0;  j < sizeOfFlag; ++j) {
						if(i == flag[j]) {
							newFlag.push_back(sigma[i]);
						}
					}
				}
			}*/
			
			//Flag
			vector<int> newFlag;
			for(int j = 0; j < sizeOfFlag; ++j) {
				if(sigma[flag[j]] != -1) {
					newFlag.push_back(sigma[flag[j]]);
				}
			}
			
			G.setFlag(newFlag);
			//G.canonRelabel();
			return G;
		}
		
		
		//--------------------------------
		//-----Print Adjacency Matrix-----
		//--------------------------------
		
		void printAdjMat() const{
			for(int i = 0; i < n; ++i) {
				for(int j = 0; j < n; ++j) {
					cout << adjMat[i][j] << " ";
				}
				cout << endl;
			}
		}
		
		//----------------------------------------
		//-----Print Adjacency Matrix to File-----
		//----------------------------------------
		
		void printAdjMatToFile(ofstream& myFile) const{
			for(int i = 0; i < n; ++i) {
				for(int j = 0; j < n; ++j) {
					myFile << adjMat[i][j] << " ";
				}
				myFile << endl;
			}
		}	
	
	
		//-----------------------------
		//-----Print Flag Vertices-----
		//-----------------------------
		
		void printFlagVertices() const{
			if (sizeOfFlag == 0) {
				cout << "There are no flag vertices." << endl;
			}	
			
			else {
				for(auto i : flag) {
					cout << i << " ";
				}
				cout << endl;
			} 
		}	
		
		
		//----------------------
		//-----Print Orbits-----
		//----------------------
		
		void printOrbits() const{
			for(int i = 0; i < int(orbits.size()); ++i) {
				cout << "Orbit " << i+1 << " contains vertices ";
				for(int j = 0; j < int(orbits[i].size()); ++j) {
					cout << orbits[i][j] << " "; 
				}
				cout << endl;
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
		
		int getNumEdges() {
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
						cout << "In setNumColors, input too small." << endl << endl;
						throw exception();
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
				cout << "First vertex out of range in getEdgeColor." << endl << endl;
				throw exception();
			}
			
			else if((j < 0) || (j >= n)) {
				cout << "Second vertex out of range in getEdgeColor." << endl << endl;
				throw exception();
			}
		
			return adjMat[i][j];
		}		
		
		
		//-------------------------
		//-----Get Canon Label-----
		//-------------------------
		
		string getCanonLabel() const{ 
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
				cout << "The edge list for color " << i << " is: {";
				for(int j = 0; j < n-1; ++j) {
					for(int k = j+1; k < n; ++k) {
						if(adjMat[j][k] == i) {
							if(comma) {
								cout << ",";
							}
							comma = true;
							cout << "{" << j << "," << k << "}";
						}
					}
				}
				cout << "}" << endl;
			}	
			cout << endl;
		}
		
		
		//---------------------
		//-----Change Edge-----
		//---------------------
		
		void changeEdge(Edge edge) {
			if((edge.a >= n) || (edge.b >= n)) {
				cout << "In changeEdge an edge has a vertex outside of the range." << endl << endl;
				throw exception();
			}
			
			if(edge.color >= numColors) { 
				cout << "In changeEdge the color is outside of the range." << endl << endl;
				throw exception();
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
			vector<int> zeros;
			
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
				cout << "In removeVertex the input is outside of the range." << endl << endl;
				throw exception();
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
		void setFlag(vector<int> const &FLAG) {
			sizeOfFlag = 0;
		
			for(auto &i : FLAG) {
				++sizeOfFlag;
				
				if((i >= n) || (i < 0)) {
					cout << "In setFlag an element is outside of the vertex range." << endl << endl;
					throw exception();
				}	
			}
			
			for(int i = 0; i < (int)FLAG.size(); ++i) {
				for(int j = i+1; j < (int)FLAG.size(); ++j) {
					if(FLAG[i] == FLAG[j]) {
						cout << "In setFlag all elements of the vector must be different." << endl << endl;
						throw exception();
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
				cout << "i is outside of range in isFlag." << endl;
				cout << "i = " << i << endl << endl;
				throw exception();
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
			vector<Edge> edges;
			for(int k : flag) {
				for(int l : flag) {
					edges.push_back({i, j, adjMat[k][l]});
					++j;
				}
				j = 0;
				++i;
			}
			
			Graph G(edges, sizeOfFlag, numColors);
			
			vector<int> flag;
			for(int i = 0; i < sizeOfFlag; ++i) {
				flag.push_back(i);
			}
			
			G.setFlag(flag);
			
			return G;
			
		}
		
		//-----------------------------
		//-----Print Flag Vertices-----
		//-----------------------------
		
		void printFlag() const{
			cout << "The vertices of the flag are (in order): ";
			for(auto j : flag) {
				cout << j << " ";
			}
			cout << endl << endl;
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

		void convertToLayer(vector<vector<int> > &newAdjMat, vector<int> &c, int &np) {		
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
		
		//This doesn't actually change the graph but gives the unique string from Nauty 
		
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
			
			vector<int> c;
			vector<vector<int> > newAdjMat;
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
			
			canonLabel = string(sgtos6(canong));
			

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
			void returnSubgraphs(const Graph&, const Graph&, vector<vector<int> >&);
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
				vector<int> subset;
				subset.resize(k);
				
				vector<vector<int> > subgraphs;
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
};

//------------------------------------
//------------------------------------
//-----Functions Involving Graphs-----
//------------------------------------
//------------------------------------

//----------------------------------------------
//-----Convert from Nauty to my Graph Class-----
//----------------------------------------------

//Hopefully won't need very much b/c we can just get canonical labeling from Nauty then get rid of that graph
//Careful with this- when converting to Nauty we use the layerGraph so this isn't a direct inverse, it drops the flags and colors and lots of other things, mostly just used for debugging

Graph convertFromNauty(const sparsegraph &sg) {
	vector<Edge> edges;
	
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

//Used in Isomorphic, when considering relabeling
//Only takes into account adjacency matrix, nothing else (no flags)
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

//Return of vectors that could be used in restriction
//Almost the same as numSubgraphs, just with different output
//Needs no flags (used in the version with flags)
void returnSubgraphsNoFlags(const Graph &H, const Graph &G, vector<vector<int> > &output) {
	//No flags
	if((H.getSizeOfFlag() != 0) || (G.getSizeOfFlag() != 0)) {
		cout << "In returnSubgraphsNoFlags H and G can't have flags." << endl << endl;
		throw exception();
	}
	
	if(H.getNumColors() != G.getNumColors()) {
		return;
	}

	if((G.getN() < H.getN())) {
		return;
	}

	if(G.getN() == H.getN()) {
		if(isomorphic(H,G)) {
			vector<int> outputEntry;
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
			vector<int> temp(G.getN(),-1);
			temp[i] = 0;
			output.push_back(temp);
		}
		
		return;
	}
	
	//Fastest way to do edges
	if(H.getN() == 2) {
		vector<int> temp(G.getN(),-1);
			
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
		vector<vector<int> > output2;
		
		returnSubgraphsNoFlags(H,Gcopy,output2);
		
		for(auto X : output2) {
			vector<int> Y;
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
		vector< vector < vector < int > > > possible; //Subsets of vertices in each collor which could create a copy of H
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
				
				vector<int> Grestriction(G.getN(),-1);
				vector<int> Hrestriction(H.getN(),-1);
				
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
				
				vector< vector< int > > possibleTemp; 
					returnSubgraphsNoFlags(H.restriction(Hrestriction),G.restriction(Grestriction),possibleTemp);
				//Make possibleTemp actually correspond to indices in G		
				for(int j = 0; j < (int)possibleTemp.size(); ++j) {
					int index = 0;
					vector<int> tempVec(G.getN(),-1);
					
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
			vector<int> maxVals; //Use in next_list
			vector<int> list;
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
				vector<int> Grestriction;
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
	vector<vector<int> > output2;
	
	returnSubgraphsNoFlags(H,Gcopy,output2);
	
	for(auto X : output2) {
		vector<int> Y;
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

//Return of vectors that could be used in restriction
void returnSubgraphs(const Graph &H, const Graph &G, vector<vector<int> > &output) {
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
			vector<int> outputEntry;
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
		vector<int> temp(G.getN(),-1);
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
				vector<int> temp(G.getN(),-1);
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
	
	vector< vector <int> > outputNoFlags;
	
	returnSubgraphsNoFlags(HNoFlags,GNoFlags,outputNoFlags);
	
	//TODO prune vertices not in flags that don't connect correctly
	
	
	for(int i = 0; i < (int)outputNoFlags.size(); ++i) {
		vector<int> possibleOutput(G.getN(),-1);
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
		cout << "In returnSubgraphsNoFlags H and G can't have flags." << endl << endl;
		throw exception();
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
		vector< vector < vector < int > > > possible; //Subsets of vertices in each collor which could create a copy of H
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
				
				vector<int> Grestriction(G.getN(),-1);
				vector<int> Hrestriction(H.getN(),-1);
				
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
				
				vector< vector< int > > possibleTemp; 
					
				
				//Make possibleTemp actually correspond to indices in G		
				if(!dominating) {						
					returnSubgraphsNoFlags(H.restriction(Hrestriction),G.restriction(Grestriction),possibleTemp);
					for(int j = 0; j < (int)possibleTemp.size(); ++j) {
						int index = 0;
						vector<int> tempVec(G.getN(),-1);
						
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
			vector<int> maxVals; //Use in next_list
			vector<int> list;
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
				vector<int> Grestriction;
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
	vector<vector<int> > output2;
	
	return (output + numSubgraphsNoFlags(H,Gcopy));
}

//-----------------------------
//-----Number of Subgraphs-----
//-----------------------------

//Return of vectors that could be used in restriction
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
	
	vector< vector <int> > outputNoFlags;
	
	returnSubgraphsNoFlags(HNoFlags,GNoFlags,outputNoFlags);
	
	//TODO prune vertices not in flags that don't connect correctly
	
	int output = 0;
	for(int i = 0; i < (int)outputNoFlags.size(); ++i) {
		vector<int> possibleOutput(G.getN(),-1);
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
		vector< vector < vector < int > > > possible; //Subsets of vertices in each collor which could create a copy of H
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
				
				vector<int> Grestriction(G.getN(),-1);
				vector<int> Hrestriction(H.getN(),-1);
				
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
				
				vector< vector< int > > possibleTemp; 
					returnSubgraphsNoFlags(H.restriction(Hrestriction),G.restriction(Grestriction),possibleTemp);
					
				if((possibleTemp.size() != 0) && dominating) {
					return true;
				}
				
				//Make possibleTemp actually correspond to indices in G		
				for(int j = 0; j < (int)possibleTemp.size(); ++j) {
					int index = 0;
					vector<int> tempVec(G.getN(),-1);
					
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
			vector<int> maxVals; //Use in next_list
			vector<int> list;
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
				vector<int> Grestriction;
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
					vector<int> possibleSubgraph(GWithFlag.getN(),-1);
					
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
vector<Graph> expandGraphs(const vector<Graph> &input, const vector<Graph> &zeros) {
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
			cout << "In expandGraphs, everything in input must be the same size." << endl << endl;
			throw exception();
		}
		
		if(input[i].getNumColors() != numColors) {
			cout << "In expandGraphs, everything in input must have the same number of colors." << endl << endl;
			throw exception();
		}
		
		if(input[i].getFlag().getCanonLabel() != flag.getCanonLabel()) {
			cout << "In expandGraphs, everything in input must have the same flag." << endl << endl;
			throw exception();
		}
		
		for(int j = 0; j < sizeOfFlag; ++j) {
			if(input[i].getFlagVertex(j) != j) {
				cout << "In expandGraphs, in input, flag vertices must be first in the correct order." << endl << endl;
				throw exception();
			}
		}
	}
	
	for(int i = 0; i < (int)zeros.size(); ++i) {
		if(zeros[i].getNumColors() != numColors) {
			cout << "In expandGraphs, at least one zero doesn't have the correct number of colors." << endl << endl;
			throw exception();
		}
		
		if(zeros[i].getFlag().getN() != 0) {
			cout << "In expandGraphs, everything in zeros must not have a flag." << endl << endl;
			throw exception();
		}
	}
	
	//Actual generation
	unordered_set<string> canonLabels;
	vector<Graph> output;
	
	//Needed everytime since they are always the same
	vector<Edge> flagEdges;
	
	for(int i = 0; i < sizeOfFlag-1; ++i) {
		for(int j = i+1; j < sizeOfFlag; ++j) {
			if(flag.getEdgeColor(i,j) != 0) {
				flagEdges.push_back({i,j,flag.getEdgeColor(i,j)});
			}
		}
	}
	
	for(auto G: input) {
		vector<Edge> GEdges = flagEdges;
		vector<int> permanentZeroDegree(n+1,0); //Value 0 for 0,1,...,sizeOfFlag
		
		for(int i = 0; i < n-1; ++i) {
			for(int j = max(sizeOfFlag,i+1); j < n; ++j) {
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
		 	vector<Edge> edges = GEdges;
		 	vector<int> zeroDegree = permanentZeroDegree;
		 	vector<int> check; //Used in an isomorphism check
		 			
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
			
			//Check if new vertex has highest degee in color 0
			bool cont = true;
			for(int j = sizeOfFlag; j < n; ++j) {
				if(zeroDegree[j] > zeroDegree[n]) {
					cont = false;
					j = n;
				}
			}
			
			if(cont) {
				//First non-flag vertex in some orbit has highest color to added vertex
				//Can only do one orbit b/c they don't play nicely together
				for(int j = 0; j < G.getNumOrbits(); ++j) {
					if(G.getOrbitSize(j) > 1) {
						for(int k = 1; k < G.getOrbitSize(j); ++k) {
							if(check[G.getOrbit(j,0)] < check[G.getOrbit(j,k)]) {
								cont = false;
								k = n;
							}
						}
						
						j = n;
					}
				}
			
				if(cont) {	
					Graph Gv(edges,n+1,numColors);
					vector<int> GvFlag;												
					
					//Removes zeros
			 		for(int j = 0; j < (int)zeros.size(); ++j) {
			 			if(subgraph(zeros[j],Gv)) {
			 				cont = false;
			 				j = zeros.size();
			 			}
			 		}
			 		
			 		if(cont) {
				 		for(int j = 0; j < sizeOfFlag; ++j) {
							GvFlag.push_back(j);
						}
				 		Gv.setFlag(GvFlag);
			 		}
				 				
				 	if(cont && (canonLabels.count(Gv.getCanonLabel()) == 0)) {
						output.push_back(Gv);
				 		canonLabels.insert(Gv.getCanonLabel());
				 	}	
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
vector<Graph> generate(const int n, const int numColors, const vector<Graph> &zeros) {
	if(numColors < 2) {
		cout << "You need at least two colors in generate." << endl << endl;
		throw exception();
	}
	
	//Check all zeros have the correct number of colors
	for(int i = 0; i < (int)zeros.size(); ++i) {
		if(zeros[i].getNumColors() != numColors) {
			cout << "In generate all zeros must have the correct number of colors." << endl << endl;
			throw exception();
		}
	}
	
	
	Graph G({},1,numColors);
	vector<Graph> output = {G};
	for(int k = 2; k <= n; ++k) {
		output = expandGraphs(output,zeros);
	}
	
	return output;
}


//-----------------------
//-----Add All Flags-----
//-----------------------

//Gives all possible ways to add 
//Used in plainFlagAlgebra to get v vectors
//flag isn't a flag but a graph
vector<Graph> addAllFlags(const Graph &G, const Graph &flag) {
	int n = G.getN();
	unordered_set<string> canonLabels;
	vector<Graph> output;
	vector<vector<int> > restrictions;
	Graph noFlag = flag;
	noFlag.removeFlag();
	
	returnSubgraphs(noFlag,G,restrictions);
	
	if(G.getSizeOfFlag() != 0) {
		cout << "In addAllFlag G must not have a flag." << endl << endl;
		throw exception();
	}
	
	if(flag.getSizeOfFlag() != flag.getN()) {
		cout << "In addAllFlag flag must have all vertices as flag vertices." << endl << endl;
		throw exception();
	}
		
	//Convert restrictions into input for adding flag
	vector<vector <int> > possibleFlags;
	possibleFlags.resize(restrictions.size(), vector<int>(flag.getN()));
	
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


//----------------------------
//-----Generate all Flags-----
//----------------------------

//Used in generateV
//TODO make faster? Use iterated approach?
vector<Graph> generateFlags(const int n, const int numColors, const vector<Graph> &zeros) {
	vector<Graph> allGraphs = generate(n,numColors,zeros);
	vector<Graph> output;
	
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		vector<int> sigma;
		unordered_set<string> canonLabels;
		
		for(int j = 0; j < n; ++j) {
			sigma.push_back(j);
		}
		
		do {
			Graph G = allGraphs[i];
			G.setFlag(sigma);
			
			if(canonLabels.count(G.getCanonLabel()) == 0) {
				canonLabels.insert(G.getCanonLabel());
				output.push_back(G);
			}
		} while(next_permutation(sigma.begin(),sigma.end()));
	}
	
	return output;
}


//--------------------
//-----Generate V-----
//--------------------

//Use in plain flag algebra
vector< vector<Graph> > generateV(const int n, const int numColors, const vector<Graph> &zeros) {
	vector< vector< Graph > > output;

	for(int i = n/2; i <= n-1; ++i) {
		int sizeOfFlag = 2*i-n;
		
		if(sizeOfFlag > 0) {
			vector<Graph> flags = generateFlags(sizeOfFlag, numColors, zeros);
			
			for(int j = 0; j < (int)flags.size(); ++j) {
				cout << "In generate v (" << i << ", " << j << ") out of (" << n-1 << ", " << flags.size() << ")" << endl;
				vector<Graph> expanded = {flags[j]};
			
				while((expanded[0].getN()*2 - sizeOfFlag) != n) {
					expanded = expandGraphs(expanded,zeros);
				}
				
				output.push_back(expanded);
			}
		}
	}
	
	return output;
}


//Use in plain flag algebra
vector< vector<Graph> > generateAllGraphsWithFlags(const int n, const int numColors, const vector<Graph> &zeros) {
	//Need generateV and this to have the same indexing
	unordered_map<string,int> fFlag; //Takes flag canonLabels and maps to index in v
	int index = 0;
	
	for(int i = n/2; i <= n-1; ++i) {
		int sizeOfFlag = 2*i-n;
		
		if(sizeOfFlag > 0) {
			vector<Graph> flags = generateFlags(sizeOfFlag, numColors, zeros);
			for(int j = 0; j < (int)flags.size(); ++j) {
				fFlag[flags[j].getCanonLabel()] = index;
				++index;
			}
		}
	}
	
	//index is number of flags
	vector<Graph> allGraphs = generate(n,numColors,zeros);
	vector < vector < Graph > > output(index);
	vector<int> index2(index);
	
	for(int i = n/2; i <= n-1; ++i) {
		int sizeOfFlag = 2*i-n;
		
		if(sizeOfFlag > 0) {
			for(int j = 0; j < (int)allGraphs.size(); ++j) {
				cout << "In generateALLGraphsWIthFlags, iteration (" << i << ", " << j << ") out of (" << n-1 << ", " << allGraphs.size() << ")" << endl;
				unordered_map<string,int> fGraph;
				int den = choose(n, sizeOfFlag);
			
				vector<int> X(sizeOfFlag);
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
							output[fFlag[flag.getCanonLabel()]].push_back(G);
						}
						
						else {
							output[fFlag[flag.getCanonLabel()]][fGraph[G.getCanonLabel()]].setCoefficient(output[fFlag[flag.getCanonLabel()]][fGraph[G.getCanonLabel()]].getCoefficient()+Frac(1,den));
						}
						
					} while(next_permutation(X.begin(), X.end()));
				} while(nextSubset(X,n,sizeOfFlag));
			}
		}
	}
	
	return output;
}


//----------------------
//-----Random Graph-----
//----------------------

//p[0] + p[1] + ... + p[numColors-1] = 1
Graph randomGraph(const int n, const int numColors, const vector<double> &p, const int flagSize) {
	if((int)p.size() != numColors) {
		cout << "In randomGraph, p must have size numColors." << endl;
		cout << "P size = " << p.size() << endl;
		cout << "numColors = " << numColors << endl << endl;
		throw exception();
	}
	
	if(flagSize > n) {
		cout << "In randomGraph, we can't have flagSize > n." << endl;
		cout << "flagSize = " << flagSize << endl;
		cout << "n = " << n << endl;
		throw exception();
	}
	
	if(numColors < 0) {
		cout << "In randomGraph, we can't have numColors < 0." << endl;
		cout << "numColors = " << numColors << endl << endl;
		throw exception();
	}
	
	if(n < 0) {
		cout << "In randomGraph, we can't have n < 0." << endl;
		cout << "n = " << n << endl << endl;
		throw exception();
	}
	
	if(flagSize < 0) {
		cout << "In randomGraph, we can't have flagSize < 0." << endl;
		cout << "flagSize = " << flagSize << endl << endl;
		throw exception();
	}

	double temp = 0;
	for(int i = 0; i < numColors; ++i) {
		temp = temp + p[i];
	}
	
	if(abs(temp - 1) > 0.000001) {
		cout << "In randomGraph, p[0] + ... + p[numColors-1] = 1." << endl;
		for(int i = 0; i < numColors; ++i) {
			cout << "p[" << i << "] = " << p[i] << endl;
		}
		cout << endl;
		throw exception();
	}
	
	vector<Edge> Edges;
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
	
	Graph G(Edges,n,numColors);
	
	if(flagSize > 0) {
		vector<int> flag;
		
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
		
		
		G.setFlag(flag);
	}
	
	return G;
}


//------------------------------
//-----Uniform Random Graph-----
//------------------------------

Graph uniformRandomGraph(const int n, const int numColors, const int flagSize) {
	vector<double> p(numColors,1./numColors);

	return randomGraph(n, numColors, p, flagSize);
}


//------------------------
//------------------------
//-----Equation Class-----
//------------------------
//------------------------

class Equation {

	private: 
	
		vector<Graph> variables;
		int numVariables;
		vector<Graph> zeros; //No Flag on zeros (remove flag from Graph to check subgraph)
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
		//unsafe doesn't do combine (should also fix combine probably)
		Equation(vector<Graph> const &VARIABLES, vector<Graph> const &ZEROS, Frac ANS, int TYPE, bool safe = true) {
			variables = VARIABLES;
			zeros = ZEROS;
			ans = ANS;
			type = TYPE;
			numVariables = variables.size();
			numZeros = zeros.size();
			
			//Check type is -1,0,or 1
			if((type != 0) && (type != 1)) {
				cout << "In equation constructor, type must be 0, or 1." << endl << endl;
				throw exception();
			}
			
			//Check if all graphs have same number of colors
			if(numVariables != 0) {
				numColors = variables[0].getNumColors();
			}
			
			for(auto G: variables) {
				if(G.getNumColors() != numColors) {
					cout << "All graphs in Equation Constructor must have same number of colors." << endl << endl;
					throw exception();
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
			//TODO Do we actually want the same adjacency matrices? 
			if(numVariables != 0) {
				for(int i = 1; i < numVariables; ++i) {
					if(!isomorphic(variables[i].getFlag(), variables[0].getFlag())) {
						cout << "All flag have to be isomorphic in equation constructor." << endl << endl; 
						throw exception();
					}
				}
			}
			
			fixZeros(); 
			
			valid = checkValid();
			
			//Not entirely sure I want to throw an exception here, I'll change it later if it doesn't work right
			/*if(!valid) {
				cout << "Equation not valid - no variables, and it is false." << endl << endl;
				throw exception(); 
			}*/
			
			//Check for isomorphism in Variables
			//If isomorphic add the coefficients

			combine();
			
		}
		
		
		//-----------------
		//-----Combine-----
		//-----------------
		
		//Checks for isomorphisms and the adds coefficients
		//Private?
		void combine() {
			for(int i = 0; i < numVariables-1; ++i) {
				for(int j = i+1; j < numVariables; ++j) {
					if(isomorphic(variables[i],variables[j])) {
						variables[i].setCoefficient(variables[i].getCoefficient() + variables[j].getCoefficient());
						variables.erase(variables.begin()+j);
						--j;
						--numVariables;
					}
				}
			}
			
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
						cout << "All zeros in Equation Constructor must have same number of colors as variables." << endl << endl;
						throw exception();
					}
				}
			}
			
			if((numVariables == 0) && (numZeros != 0)) {
				numColors = zeros[0].getNumColors();
				
				for(auto G: zeros) {
					if(G.getNumColors() != numColors) {
						cout << "All zeros in Equation Constructor must have same number of colors." << endl << endl;
						throw exception();
					}
				}
			}
			
			//Check that all zeros don't have flag
			//I'm not entirely sure that this needs to be true, but I can't think of any reasonable examples were it isn't though
			for(auto G: zeros) {
				if(G.getSizeOfFlag() != 0) {
					cout << "Zeros can't have any flags." << endl << endl;
					throw exception();
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
				cout << "Equation not valid - no variables, and it is false." << endl << endl;
				throw exception(); 
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
				cout << "Equation not valid - no variables, and it is false." << endl << endl;
				throw exception(); 
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
				cout << "Trying to get variable outside of range." << endl;
				cout << "Index = " << i << endl;
				cout << "Number of Variables = " << numVariables << endl;
				throw exception();
			}
			
			return variables[i];
		}
		
		
		//------------------
		//-----Get Zero-----
		//------------------
		
		Graph getZero(int i) const {
			if((i < 0) || (i >= numZeros)) {
				cout << "Trying to get zero outside of range: i = " << i << ", numZeros = " << numVariables << endl << endl;
				throw exception();
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
				cout << "Equation not valid - no variables, and it is false." << endl << endl;
				throw exception(); 
			}
		}
		
		
		//-------------------------
		//-----Print Variables-----
		//-------------------------
		
		void printVariables() const {
			if(numVariables == 0) {
				cout << "There are no variables." << endl << endl;
			}
			
			else {
				for(int i = 0; i < numVariables; ++i) {
					cout << "Variable " << i << " has adjacency matrix: " << endl;
					variables[i].printAdjMat();
					if(variables[i].getSizeOfFlag() != 0) {
						cout << "The flag vertices are: ";
						variables[i].printFlagVertices();
					}
					else {
						cout << "There are no flag vertices." << endl;
					}
					cout << "It has coefficient: " << variables[i].getCoefficient() << endl << endl;
					//Could add this back in, but it doesn't seem necessary unless debugging
					//cout << "Its canon label is: " << variables[i].getCanon() << endl << endl;
				}
			}
		}
		
		
		//----------------------
		//-----Prints Zeros-----
		//----------------------
		
		void printZeros() const {
			if(numZeros == 0) {
				cout << "There are no zeros." << endl << endl;
			}
			
			else {
				for(int i = 0; i < numZeros; ++i) {
					cout << "Zero " << i << " has adjacency matrix: " << endl;
					zeros[i].printAdjMat();
					cout << endl;
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
bool operator==(const Equation &eq1, const Equation &eq2) {
	vector<Graph> variables;
	vector<Graph> zeros; //No Flag on zeros (remove flag from Graph to check subgraph)
	
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
			if(isomorphic(eq1.getVariable(i),eq2.getVariable(j))) {
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
		cout << "When adding equations, they must have the same set of zeros." << endl << endl;
		throw exception();
	}
	
	vector<Graph> zeros;
	
	bool val;
	for(int i = 0; i < numZeros; ++i) {
		val = false;
		for(int j = 0; j < numZeros; ++j) {
			if(isomorphic(eq1.getZero(i),eq2.getZero(j))) {
				val = true;
				j = eq1.getNumZeros();
				zeros.push_back(eq1.getZero(i));
			}
		}
		
		if(!val) {
			cout << "When adding equations, they must have the same set of zeros." << endl << endl;
			throw exception();
		}
	}
	
	//Constructor actually does the adding already
	vector<Graph> variables;
	
	for(int i = 0; i < eq1.getNumVariables(); ++i) {
		variables.push_back(eq1.getVariable(i));
	}
	
	for(int i = 0; i < eq2.getNumVariables(); ++i) {
		variables.push_back(eq2.getVariable(i));
	}
	
	for(int i = 0; i < (int)variables.size(); ++i) {
		for(int j = i+1; j < (int)variables.size(); ++j) {
			if(variables[i].getCanonLabel() == variables[j].getCanonLabel()) {
				variables[i].setCoefficient(variables[i].getCoefficient() + variables[j].getCoefficient());
				variables.erase(variables.begin()+j);
			}
		}
	}
	
	return Equation(variables, zeros, eq1.getAns() + eq2.getAns(), type);
}

//If adding Graph to Equation, like adding G = 0 to it
Equation operator+(const Graph &G, const Equation &eq) {
	vector<Graph> variables = {G};
	
	vector<Graph> zeros;
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

//First Define for graphs - put into equation with = 0
//Really more of a vector of graphs, but in this form it allows us to use overridden +
//TODO add edge between parts in a smarter way?
Equation multiply(const Graph &G1, const Graph &H1, const vector<Graph> &zeros) {
	if(G1.getN() > H1.getN()) {
		return multiply(H1,G1,zeros);
	} 

	vector<int> HReorder;
	vector<int> GReorder;
	
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
		cout << "Flags must be the same in multiplication." << endl << endl;
		G.printEdges();
		G.printFlag();
		H.printEdges();
		H.printFlag();
		throw exception();
	}
	
	//Check that numColors are the same
	if(G.getNumColors() != H.getNumColors()) {
		cout << "Graphs in multiplication must have the same number of colors." << endl << endl;
		throw exception();
	}
	
	//Can't multiply empty graph
	if((G.getN() == 0) || (H.getN() == 0)) {
		cout << "Can't multiply empty graphs." << endl << endl;
		throw exception();
	}
	
	//Need to reorder vertices so flags are first
	int sizeOfFlag = G.getSizeOfFlag();
	int c = G.getNumColors();
	int nG = G.getN();
	int nH = H.getN();
	
	vector<Graph> variables;
	unordered_set<string> variablesSet;
	
	vector<Edge> permanentEdges;
	//First vertices G
	for(int j = 0; j < nG-1; ++j) {
		for(int k = j+1; k < nG; ++k) {
			if(G.getEdgeColor(j,k) != 0) {
				permanentEdges.push_back({j,k,G.getEdgeColor(j,k)});
			}
		}
	}

	//Flag Plus Last Vertices = H
	for(int j = 0; j < nH-1; ++j) {
		for(int k = j+1; k < nH; ++k) {
			if(H.getEdgeColor(j,k) != 0) {
				if((j < sizeOfFlag) && (k < sizeOfFlag)) {
					permanentEdges.push_back({j,k,H.getEdgeColor(j,k)});
				}
					
				else if((j < sizeOfFlag) && (k >= sizeOfFlag)) {
					permanentEdges.push_back({j,k+nG-sizeOfFlag,H.getEdgeColor(j,k)});
				}
					
				else {
					permanentEdges.push_back({j+nG-sizeOfFlag,k+nG-sizeOfFlag,H.getEdgeColor(j,k)});
				}
			}
		}
	}
	
	//Create Flag
	vector<int> flag;
		
	//Can do this because they are in canonical form
	for(int j = 0; j < sizeOfFlag; ++j) {
		flag.push_back(G.getFlagVertex(j));
	}
	
	if(nG == sizeOfFlag) {
		Equation eq({H},zeros, Frac(0,1), 0, false);
		return eq;
	}

	//Create all possible variables
	//In first vertex (index=sizeOfFlag) order of edges in each orbit doesn't matter
	vector< vector< vector<Edge> > > edgesTemp;
	
	int index1 = 0;
	
	for(int i = 0; i < H.getNumOrbits(); ++i) {
		if(H.getOrbit(i,0) >= sizeOfFlag) {
			edgesTemp.push_back(vector< vector<Edge> >());
		
			vector<int> colors(c,0);		
			colors[0] = H.getOrbitSize(i);
			
			int index2 = 0;
			do {
				edgesTemp[index1].push_back(vector<Edge>());
					
				int index3 = 0;
				for(int j = 0; j < c; ++j) {
					for(int k = 0; k < colors[j]; ++k) {
						edgesTemp[index1][index2].push_back({sizeOfFlag,H.getOrbit(i,index3)+nG-sizeOfFlag,j});
						++index3;
					}	
				}
				++index2;
			} while(nextOrderedSum(colors));
			
			++index1;
		}
	}
	
	//Combine info from each orbit from edgesTemp
	vector< vector<Edge> > firstVertexEdge;
	vector<int> maxVals;
	vector<int> list(H.getNumOrbits()-sizeOfFlag,0);
	
	for(int i = 0; i < H.getNumOrbits(); ++i) {
		if(H.getOrbit(i,0) >= sizeOfFlag) {
			maxVals.push_back(H.getOrbitSize(i));
		}
	}
	
	index1 = 0;
	do {
		firstVertexEdge.push_back(vector<Edge>());
		
		int index2 = 0;
		
		for(int i = 0; i < H.getNumOrbits(); ++i) {
			if(H.getOrbit(i,0) >= sizeOfFlag) {
				for(int j = 0; j < H.getOrbitSize(i); ++j) {
					firstVertexEdge[index1].push_back(edgesTemp[index2][list[index2]][j]);
				}
				
				++index2;
			}
		}
		
		++index1;
	
	} while(nextList(list,maxVals));
	
	for(auto partialEdges : firstVertexEdge) {
	for(int i = 0; i < myPow(c,(nG - sizeOfFlag-1)*(nH-sizeOfFlag)); ++i) { //c-ary mask
		//Convert i to c-ary number
		vector<int> ary;
		int temp = i;
		
		for(int j = 0; j < (nG-sizeOfFlag-1)*(nH-sizeOfFlag); ++ j) {
			int temp2 = temp/myPow(c,(nG-sizeOfFlag-1)*(nH-sizeOfFlag)-j-1); //Need to be careful with rounding? 
		   ary.push_back(temp2);
		   temp = temp-temp2*myPow(c,(nG-sizeOfFlag-1)*(nH-sizeOfFlag)-j-1);
		}
		
		vector<Edge> newEdges;
		
		//Add edges from i
		//Easier to just iterate through than to make a function
		int index = 0;
		for(int j = nG; j < nG+nH-sizeOfFlag; ++j) {
			for(int k = sizeOfFlag+1; k < nG; ++k) {
				if(ary[index] != 0) {
					newEdges.push_back({k,j,ary[index]});
				}
				++index;
			}
		}
		
		for(int j = 0; j < (int)partialEdges.size(); ++j) {
			newEdges.push_back(partialEdges[j]);
		}
		
		//May assume that within each orbit the degrees are decreasing
		bool cont = true;

		vector<int> degreeG(nG-sizeOfFlag,0); //0 indexed rather than sizeOfFlag indexed
		for(int j = 0; j < (int)newEdges.size(); ++j) {
			if(newEdges[j].color == 1) {
				++degreeG[newEdges[j].a-sizeOfFlag];
			}
		}
			
		for(int j = 0; j < G.getNumOrbits(); ++j) {	
			for(int k = 1; k < G.getOrbitSize(j); ++k) { //Makes it so we don't look at flag vertices
				if(degreeG[G.getOrbit(j,k-1)-sizeOfFlag] < degreeG[G.getOrbit(j,k)-sizeOfFlag]) {
					cont = false;
					k = G.getOrbitSize(j);
				}
			}	
				
			if(!cont) {
				j = G.getNumOrbits();
			}
		}
		
		if(cont) {
			vector<Edge> edges = permanentEdges;
			for(int j = 0; j < (int)newEdges.size(); ++j) {
				edges.push_back(newEdges[j]);
			}

			
			Graph GH(edges, nG+nH-sizeOfFlag, c);
			GH.setFlag(flag);
			
			if(variablesSet.count(GH.getCanonLabel()) == 0) {
				variablesSet.insert(GH.getCanonLabel());
				variables.push_back(GH);
			}
		}
	}
	}

	Equation eq(variables,zeros, Frac(0,1), 0, false); //When converting to equation, removes isomorphisms
	
	//Makes coefficients correct
	for(int i = 0; i < eq.getNumVariables(); ++i) {
		int num = 0;
		
		vector<vector<int> > subgraphs; 
		returnSubgraphs(G,eq.getVariable(i),subgraphs);
		
		for(auto X: subgraphs) {
			vector<int> restrictionComp;
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
		cout << "In multiplication, both types can't be 1, inequalities possibly not preserved." << endl << endl;
		throw exception();
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
		cout << "Can't multiply two equations together when one is empty." << endl << endl;
		throw exception();
	}
	
	//Check the flags are the same
	if(!isomorphic(eq1.getVariable(0).getFlag(),eq2.getVariable(0).getFlag())) {
		cout << "Flags must be the same in multiplication." << endl << endl;
		throw exception();
	}
	
	//Check that numColors are the same
	if(eq1.getVariable(0).getNumColors() != eq2.getVariable(0).getNumColors()) {
		cout << "Graphs in multiplication must have the same number of colors." << endl << endl;
		throw exception();
	}
	
	//Check if Zeros are the same
	int numZeros = eq1.getNumZeros();
	
	if(numZeros != eq2.getNumZeros()) {
		cout << "When adding equations, they must have the same set of zeros." << endl << endl;
		throw exception();
	}
	
	vector<Graph> zeros;
	
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
			cout << "When adding equations, they must have the same set of zeros." << endl << endl;
			throw exception();
		}
	}

	
	//Distributive
	Equation eq = multiply(eq1.getVariable(0), eq2.getVariable(0), zeros);
	
	for(int i = 0; i < eq1.getNumVariables(); ++i) {
		for(int j = 0; j < eq2.getNumVariables(); ++j) {
			if((i != 0) || (j != 0)) {
				eq = eq + multiply(eq1.getVariable(i), eq2.getVariable(j), zeros);
			}
		}
	}
	
	eq.setType(type);
	eq.setAns(eq1.getAns() * eq2.getAns());
	
	//for(int i = 0; i < (int)zeros.size(); ++i) {
		//eq.addZero(zeros[i]);
	//}
	
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

//-------------------------
//-----Resize Equation-----
//-------------------------

//Currently only works for equations with no flags
//Make this a class function?
//Use generate to get generated
Equation resize(const Equation &eq, const vector<Graph> &generated) {
	if(eq.getNumVariables() == 0) {
		cout << "In resize, need at least one graph in the equation to resize." << endl << endl;
		throw exception();
	}
	
	for(int i = 0; i < eq.getNumVariables(); ++i) {
		if(eq.getVariable(i).getFlag().getN() != 0) {
			cout << "In resize, we don't allow flags (there is no reason for this I'm just lazy, so future me you can implement this)." << endl << endl;
			throw exception();
		}
		
		if(eq.getVariable(i).getN() > generated[0].getN()) {
			cout << "In resize, we can't have the variables having larger size than things in generated." << endl;
			throw exception();
		}
	}
	
	vector<Graph> zeros;
	for(int i = 0; i < eq.getNumZeros(); ++i) {
		zeros.push_back(eq.getZero(i));
	}	
	
	Frac ans = eq.getAns();
	
	vector<Graph> resized;
	Equation output({},zeros,eq.getAns(),eq.getType());

	for(int i = 0; i < eq.getNumVariables(); ++i) {
		resized.clear();
		
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
		output = output+temp;
	}
	
	output.setAns(ans);
	
	return output;
}


//--------------------------
//-----Print Matrix Col-----
//--------------------------

//Used in Python code for solving equations
//Everything is how the graphs are indexed (so must include everything)
//Only used for == (right now, may change later)
//Last entry is the ans
void printMatrixCol(const Equation &everything, const Equation &eq) {
	ofstream myFile;
	myFile.open("matrix.txt",ios::out | ios::app);
	//Check if type of eq is right
	if(eq.getType() == 1) {
		cout << "Can only call augmented Matrix with a type of 0." << endl << endl; 
		throw exception();
	}
	
	//Make sure no Equations are empty
	if((everything.getNumVariables() == 0) || (eq.getNumVariables() == 0)) {
		cout << "Can't use augmentedMatrix if an equation has no variables." << endl << endl;
		throw exception();
	}
	
	//Check the flags are the same
	if(!isomorphic(eq.getVariable(0).getFlag(),everything.getVariable(0).getFlag())) {
		cout << "Flags must be the same in multiplication." << endl << endl;
		throw exception();
	}
	
	int numVariables = everything.getNumVariables();
	
	vector<int> var(numVariables+1,0);
	
	//Determine lcm of all denominators so everything is a fraction
	long long int den = lcm(eq.getVariable(0).getCoefficient().getDen(),eq.getAns().getDen());
	
	for(int i = 1; i < eq.getNumVariables(); ++i) {
		den = lcm(den,eq.getVariable(i).getCoefficient().getDen());
	}
	
	var[numVariables] = eq.getAns().getNum() * den /eq.getAns().getDen(); 
	
	bool val;
	
	for(int i = 0; i < eq.getNumVariables(); ++i) {
		val = false;
		
		for(int j = 0; j < numVariables; ++j) {
			if(isomorphic(eq.getVariable(i),everything.getVariable(j))) {
				var[j] = eq.getVariable(i).getCoefficient().getNum() * den / eq.getVariable(i).getCoefficient().getDen();
				val = true;
				j = numVariables;
			}
		}
		
		if(!val) {
			cout << "In printMatrixCol, eq has a variable not in everything." << endl << endl;
			throw exception();
		}		
	}
	
	for(int i = 0; i < numVariables; ++i) {
		myFile << var[i] << ",";
	}
	
	myFile << var[numVariables] << endl;
	
}

//-----------------------------------
//-----Generate Augmented Matrix-----
//-----------------------------------

//TODO Set would probably be better data structure than vector 
void augmentedMatrix(vector<Equation> &known, int n) {

	//Check n
	if(n < 2) {
		cout << "Need n of at least 2 for augmentedMatrix." << endl << endl;
		throw exception();
	}
	
	//No equation can have 0 varaibles
	for(auto eq: known) {
		if(eq.getNumVariables() == 0) {
			cout << "No empty equation in augmentedMatrix." << endl << endl;
			throw exception();
		}
	}
	
	//Check that numColors is the same in everything
	int c = known[0].getVariable(0).getNumColors();
	
	for(auto eq: known) {
		for(int i = 0; i < eq.getNumVariables(); ++i) {
			if(eq.getVariable(i).getNumColors() != c) {
				cout << "In augmentedMatrix, everything must all numColors must be the same." << endl << endl;
				throw exception();
			}
		}
	}	
	
	//Check that everything in known is type 0
	for(auto eq: known) {
		if(eq.getType() != 0) {
			cout << "In augmentedMatrix, everything in known must be type 0." << endl << endl;
			throw exception();
		}	
	}
	
	//Checks all zeros are the same
	for(auto eq1: known) {
		for(auto eq2: known) {
			if(eq1.getNumZeros() != eq2.getNumZeros()) {
				cout << "In augmentedMatrix, all known must have the same zeros (in the same order)." << endl << endl;
				throw exception();
			}
			
			for(int i = 0; i < eq1.getNumZeros(); ++i) {
				if(!isomorphic(eq1.getZero(i),eq2.getZero(i))) { 
					cout << "In augmentedMatrix, all known must have the same zeros (in the same order)." << endl << endl;
					throw exception();
				}
			}
		}
	}
	
	//Gets 0 graphs
	vector<Graph> zeros;
	
	for(int i = 0; i < known[0].getNumZeros(); ++i) {
		zeros.push_back(known[0].getZero(i));
	}
	
	//Checks all known have equations of the same size
	for(auto eq: known) {
		int number = eq.getVariable(0).getN();
		for(int i = 1; i < eq.getNumVariables(); ++i) {
			if(number != eq.getVariable(i).getN()) {
				cout << "All equations in known in augmentedMatrix must have same size." << endl << endl;
				throw exception();
			}
		}
	}
		
	//List of all flags
	vector<Graph> flags;
	bool val;
	
	for(auto eq: known) {
		for(int i = 0; i < eq.getNumVariables(); ++i) {
			Graph G = eq.getVariable(0);
			val = false;
			
			for(auto H: flags) {
				if(isomorphic(G.getFlag(),H)) {
					val = true;
				}
			}
			
			if(!val) {
				flags.push_back(G.getFlag());
			}
		}
	}
	
	vector<Graph> flagsPlusVertex;
	
	for(auto G: flags) {
		if(G.getFlag().getN() == 0) {
			Graph V({},1,c);
			known.push_back(Equation({V},zeros,Frac(1,1),0));
		}
		
		if(G.getFlag().getN() != 0) {
			vector<Graph> variables;		
					
			for(int i = 0; i < myPow(c,G.getFlag().getN()); ++i) {
				//Convert i to c-ary number
				vector<int> ary;
				int temp = i;
				
				for(int j = 0; j < G.getFlag().getN(); ++ j) {
					int temp2 = temp/myPow(c,G.getFlag().getN()-j-1); //Need to be careful with rounding? 
					ary.push_back(temp2);
					temp = temp-temp2*myPow(c,G.getFlag().getN()-j-1);
					
					if((temp2 < 0) || (temp2 >= c)) { //This really shouldn't happen I've checked, but better safe than sorry
						cout << "Something went wrong in augmentedMatrix." << endl << endl;
						throw exception();
					}
				}
							
				Graph H = G.clone();
				H.addVertex();
				for(int j = 0; j < G.getFlag().getN(); ++j) {
					H.changeEdge({j,G.getFlag().getN(),ary[j]});
				}		
				variables.push_back(H);
			}
			
			known.push_back(Equation(variables,zeros,Frac(1,1),0));
		}
	}
	
	//Generate new known equations
	vector<Equation> newKnown;
	
	for(int i = 2; i <= n; ++i) {
		newKnown.clear(); //use this so for loops don't change size
		
		for(int j = 0; j < (int)known.size(); ++j) {
			for(int k = j; k < (int)known.size(); ++k) {
				if(isomorphic(known[j].getVariable(0).getFlag(), known[k].getVariable(0).getFlag())) {
					if(known[j].getVariable(0).getN() + known[k].getVariable(0).getN() - known[j].getVariable(0).getFlag().getN() == i) {
						Equation eq = known[j]*known[k];
						newKnown.push_back(eq);
						cout << known.size() + newKnown.size() << endl;
					}
				}
			}
		}
		
		for(auto eq: newKnown) {
			known.push_back(eq);
		}
	}
	
	//Converts everything to be about graphs with no flags
	//May change this, not sure if it's best, but it good for current project
	for(int i = 0; i < (int)known.size(); ++i) {
		known[i].averageAll();
	}
	
	//Could fix to only start at size of smallest flag +1
	Graph V({},1,c);
	Equation everything({V},zeros,Frac(0,1),0);
	
	for(int i = 2; i <= n; ++i) {
		Equation temp1({V},zeros,Frac(0,1),0);
		Equation temp2({V},zeros,Frac(0,1),0);
		
		for(int j = 2; j <= i; ++j) {
			temp1 = temp1 * temp2;
		} 
		
		everything = everything + temp1;
	}
	
	//Prints the variables first
	ofstream myFile;
	myFile.open("graphs.txt");
	
	for(int i = 0; i < everything.getNumVariables(); ++i) {
		myFile << "Variable " << i << " has adjacency matrix: " << endl;
		for(int j = 0; j < everything.getVariable(i).getN(); ++j) {
			for(int k = 0; k < everything.getVariable(i).getN(); ++k) {
				myFile << everything.getVariable(i).getEdgeColor(j,k) << " ";
			}
			myFile << endl;
		}				
		myFile << "There are no flag vertices." << endl;
		myFile << "It has coefficient: " << everything.getVariable(i).getCoefficient() << endl;
		myFile << "Its canon label is: " << everything.getVariable(i).getCanonLabel() << endl << endl;
	}
	myFile.close();
	
	//Specific to this problem, can remove
	cout << "The variables with no empty edges of full size are: ";
	for(int i = 0; i < everything.getNumVariables(); ++i) {
		if(everything.getVariable(i).getN() == n) {
			bool temp = false;
				
			for(int j = 0; j < n-1; ++j) {
				for(int k = j+1; k < n; ++k) {
					if(everything.getVariable(i).getEdgeColor(j,k) == 0) {
						temp = true; 
						k = n;
						j = n;
					}
				}
			}
			
			if(!temp) {
				cout << i << " ";
			}
		}
	}
	cout << endl << endl;
	
	//Clears matrix.txt
	myFile.open("matrix.txt");
	myFile.close();
	
	int numEquations = 0;
	for(auto eq: known) {
		printMatrixCol(everything, eq);
		++numEquations;
	}
	
	cout << endl << "The number of variables is: " << everything.getNumVariables() << endl;
	cout << "The number of equations there are is: " << numEquations << endl;
}


//-------------------------------------------
//-----Plain Flag Algebra (Key Function)-----
//-------------------------------------------

//Prints to plainFlagAlgerba1.txt & plainFlagAlgebra2.txt necessary files for python SDP code 
//f can be thought of as a linear combo of all graphs that we want to max/min
//v is vector multiplied by matrix in SDP
void plainFlagAlgebra(vector<Graph> &f, vector<Graph> &v, vector<Graph> &zeros, vector<Equation> &known) {
	cout << "Starting plainFlagAlgebra." << endl;
	int vSize = v.size();
	int fSize = f.size();
	int zerosSize = zeros.size();
	int knownSize = known.size();
	
	if(vSize == 0) {
		cout << "Need at least one vector in v in plainFlagAlgebra." << endl << endl;
		throw exception();
	}

	if(fSize == 0) {
		cout << "Need at least one graph in f in plainFlagAlgebra." << endl << endl;
		throw exception();
	}
	
	//Check if every equation in known has at least one variable
	for(int i = 0; i < knownSize; ++i) {
		if(known[i].getNumVariables() == 0) {
			cout << "Equations in known must have at least one graph in plainFlagAlgebra." << endl << endl;
			throw exception();
		}
	}
	
	//Check if every equation in known is type 1
	for(int i = 0; i < knownSize; ++i) {
		if(known[i].getType() != 1) {
			cout << "Every equation in known must by type 1 in plainFlagAlgebra." << endl << endl;
			throw exception();
		}
	}
	
	int vN = v[0].getN();
	int vFlagSize = v[0].getFlag().getN();
	int fN = f[0].getN();

	//Check v has all same size
	for(int i = 1; i < vSize; ++i) {
		if(v[i].getN() != vN) {
			cout << "All graphs in v in plainFlagAlgebra must be the same size." << endl << endl;
			throw exception();
		}
	}
	
	//Check if v has all same flags
	Graph flag = v[0].getFlag();
	for(int i = 1; i < vSize; ++i) {
		if(!isomorphic(flag,v[i].getFlag())) {
			cout << "All graphs in v in plainFlagAlgebra must have the same flag." << endl << endl;
			throw exception();
		}
	}
	
	//Graphs in f can't have flags
	for(int i = 0; i < fSize; ++i) {
		if(f[i].getFlag().getN() != 0) {
			cout << "All graphs in f in plainFlagAlgebra must not have any flags." << endl << endl;
			throw exception();
		}
	}

	//Graphs in known can't have flags
	//Maybe fix?
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			if(known[i].getVariable(j).getFlag().getN() != 0) {
				cout << "All graphs in known in plainFlagAlgebra must not have any flags." << endl << endl;
				throw exception();
			}
		}
	}
	
	//Make sure everything has correct number of colors
	int numColors = f[0].getNumColors();
	
	for(int i = 1; i < fSize; ++i) {
		if(f[i].getNumColors() != numColors) {
			cout << "Everything in plainFlagAlgebra in f must have the same number of colors." << endl << endl;
			throw exception();
		}
	}
	
	for(int i = 0; i < vSize; ++i) {
		if(v[i].getNumColors() != numColors) {
			cout << "Everything in plainFlagAlgebra in v must have the same number of colors." << endl << endl;
			throw exception();
		}
	}
	
	for(int i = 0; i < zerosSize; ++i) {
		if(zeros[i].getNumColors() != numColors) {
			cout << "Everything in plainFlagAlgebra in zeros must have the same number of colors." << endl << endl;
			throw exception();
		}
	}
	
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			if(known[i].getVariable(j).getNumColors() != numColors) {
				cout << "Everything in plainFlagAlgebra in known must have the same number of colors." << endl << endl;
				throw exception();
			}
		}
	}
	
	int n = 2*vN-vFlagSize;
	
	//Not sure if strictly necessary, but it would probably give trash bounds otherwise
	if(fN > n) {
		cout << "Make v large enough so it has vertices at least as many vertices when multiplied by itself as n." << endl << endl;
		throw exception();
	}
	
	cout << endl;
	
	Equation fEq(f,zeros,Frac(1,1),0); //Type doesn't matter
	vector<Edge> edges {};
	Graph H(edges,1,numColors);
	Equation eq1({H},zeros,Frac(1,1),0);
	Equation eq2({H},zeros,Frac(1,1),0);
	
	cout << "Generating graphs to be used in resize." << endl;
	vector<Graph> allGraphs = generate(n,numColors,zeros);
	
	cout << endl << "Resizing f." << endl;
	Equation fEqResized = resize(fEq,allGraphs);
	f.clear();
	for(int i = 0; i < fEqResized.getNumVariables(); ++i) {
		f.push_back(fEqResized.getVariable(i));
	}
	cout << endl;
	
	fSize = f.size();
	
	//Resize known
	cout << "Resizing known." << endl;
	for(int i = 0; i < knownSize; ++i) {
		cout << i+1 << " out of " << knownSize << endl;
		known[i] = resize(known[i],allGraphs);
	}
	cout << endl;
	
	//Calculates everything
	Graph G({},1,numColors);
	Equation eq({G},zeros,Frac(1,1),0);
	
	if(n == 1) {
		cout << "Plain flag algebra method not set up for graphs with one vertex." << endl << endl;
		throw exception();
	}
	
	cout << "Creating everything equation." << endl;
	Equation everything (allGraphs,zeros,Frac(1,1),0);
	cout << endl;
	
	ofstream myFile;
	myFile.open("Duals.txt");
	
	for(int i = 0; i < everything.getNumVariables(); ++i) {
		myFile << "Label: " << i+1 << endl;
		everything.getVariable(i).printAdjMatToFile(myFile);
		myFile << endl;
	}

	vector< vector< vector<Frac> > > A; //Give adjacency matrix values to be printed
	vector<Frac> B; //Gives numbers to be printed
	vector< vector<Frac> > C; //From Known
	
	//Not sure if it will always initalize right with my own data structure
	Frac zeroFrac(0,1);
	
	A.resize(everything.getNumVariables());
	for(int i = 0; i < everything.getNumVariables(); ++i) {
		A[i].resize(vSize);
		for(int j = 0; j < vSize; ++j) {
			A[i][j].resize(vSize,zeroFrac);
		}
	}
	
	for(int i = 0; i < everything.getNumVariables(); ++i) {
		B.push_back(zeroFrac);
	}
	
	C.resize(knownSize);
	for(int i = 0; i < knownSize; ++i) {
		C[i].resize(everything.getNumVariables()+1,zeroFrac);
	}
	
	unordered_map<string, int> everythingMap;
			
	for(int i = 0; i < everything.getNumVariables(); ++i) {
		everythingMap.insert(pair<string, int>(everything.getVariable(i).getCanonLabel(),i));
	}
	
	for(int i = 0; i < vSize; ++i) {	
		for(int j = i; j < vSize; ++j) {
			cout << "In calculating A, (" << i+1 << ", " << vSize << "), (" << j+1 << ", " << vSize << ")" << endl;
			Equation temp = multiply(v[i],v[j],zeros);
			cout << "Finished multiplying." << endl;
			temp.averageAll();
			cout << "Finished averaging." << endl;
			
			for(int k = 0; k < temp.getNumVariables(); ++k) {
				if(i == j) {
					A[everythingMap[temp.getVariable(k).getCanonLabel()]][i][j] = temp.getVariable(k).getCoefficient();
				}
				
				else {
					A[everythingMap[temp.getVariable(k).getCanonLabel()]][i][j] = 2*temp.getVariable(k).getCoefficient();
				}
			}
		}
	}
	
	for(int i = 0; i < fSize; ++i) {
		for(int j = 0; j < everything.getNumVariables(); ++j) {
			if(isomorphic(f[i],everything.getVariable(j))) {
				B[j] = f[i].getCoefficient();
				j = everything.getNumVariables();
			}
		}
	}
	cout << endl;
	
	//First has 0 for ==, 1 for <=
	//Next entry is bound
	//Finally, it's the vector of coefficients in known (after it's been resized)
	for(int i = 0; i < knownSize; ++i) {
		C[i][0] = known[i].getAns();
		
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			for(int k = 0; k < everything.getNumVariables(); ++k) {
				if(isomorphic(known[i].getVariable(j),everything.getVariable(k))) {
					C[i][k+1] = known[i].getVariable(j).getCoefficient();
					k = everything.getNumVariables();
				}
			}
		}
	}
	
	ofstream myFile1;
	myFile1.open("plainFlagAlgebra1.txt");
	
	for(int i = 0; i < everything.getNumVariables(); ++i) {
		for(int j = 0; j < vSize; ++j) {
			for(int k = 0; k < vSize; ++k) {
				myFile1 << (double)A[i][j][k].getNum() / (double)A[i][j][k].getDen() << " ";
			}
			myFile1 << endl;
		}
		myFile1 << endl;
	}
	
	ofstream myFile2;
	myFile2.open("plainFlagAlgebra2.txt");
	
	for(int i = 0; i < everything.getNumVariables(); ++i) {
		myFile2 << (double)B[i].getNum() / (double)B[i].getDen() << " ";
	}
	
	ofstream myFile3;
	myFile3.open("plainFlagAlgebra3.txt");
	
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < everything.getNumVariables()+1; ++j) {
			myFile3 << (double)C[i][j].getNum() / (double)C[i][j].getDen() << " ";
		}
	}
	
	cout << "Finished plainFlagAlgebra." << endl << endl;
} 


//--------------------------------------
//-----Secondary Plain Flag Algebra-----
//--------------------------------------

//Prints to plainFlagAlgerba1.txt & plainFlagAlgebra2.txt necessary files for python SDP code 
//f can be thought of as a linear combo of all graphs that we want to max/min
//Rather than taking a v this function takes the number of vertices to compute on (n)
void plainFlagAlgebra(vector<Graph> &f, int n, vector<Graph> &zeros, vector<Equation> &known) {
	cout << "Starting plainFlagAlgebra." << endl;
	
	//The way the python script is setup we need at least one known
	if(known.size() == 0) {
		Graph K2({{}},2,f[0].getNumColors());
		Equation knownTemp({K2},zeros,Frac(1,1),1);
		known.push_back(knownTemp);
	}
	
	int fSize = f.size();
	int zerosSize = zeros.size();
	int knownSize = known.size();
		
	if(n <= 1) {
		cout << "Plain flag algebra method not set up for graphs with fewer than two vertices." << endl << endl;
		throw exception();
	}

	if(fSize == 0) {
		cout << "Need at least one graph in f in plainFlagAlgebra." << endl << endl;
		throw exception();
	}
	
	//Check if every equation in known has at least one variable
	for(int i = 0; i < knownSize; ++i) {
		if(known[i].getNumVariables() == 0) {
			cout << "Equations in known must have at least one graph in plainFlagAlgebra." << endl << endl;
			throw exception();
		}
	}
	
	//Make every equation in known type 1 (<=)
	for(int i = 0; i < knownSize; ++i) {
		if(known[i].getType() != 1) {
			vector<Graph> variablesTemp;
			
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
	}

	int fN = f[0].getN();
	
	//Graphs in f can't have flags
	for(int i = 0; i < fSize; ++i) {
		if(f[i].getFlag().getN() != 0) {
			cout << "All graphs in f in plainFlagAlgebra must not have any flags." << endl << endl;
			throw exception();
		}
	}

	//Graphs in known can't have flags
	//Maybe fix?
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			if(known[i].getVariable(j).getFlag().getN() != 0) {
				cout << "All graphs in known in plainFlagAlgebra must not have any flags." << endl << endl;
				throw exception();
			}
		}
	}
	
	//Make sure everything has correct number of colors
	int numColors = f[0].getNumColors();
	
	for(int i = 1; i < fSize; ++i) {
		if(f[i].getNumColors() != numColors) {
			cout << "Everything in plainFlagAlgebra in f must have the same number of colors." << endl << endl;
			throw exception();
		}
	}
	
	for(int i = 0; i < zerosSize; ++i) {
		if(zeros[i].getNumColors() != numColors) {
			cout << "Everything in plainFlagAlgebra in zeros must have the same number of colors." << endl << endl;
			throw exception();
		}
	}
	
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			if(known[i].getVariable(j).getNumColors() != numColors) {
				cout << "Everything in plainFlagAlgebra in known must have the same number of colors." << endl << endl;
				throw exception();
			}
		}
	}
	
	//Not sure if strictly necessary, but it would probably give trash bounds otherwise
	if(fN > n) {
		cout << "Make n large enough so it has vertices at least as many vertices when multiplied by itself as n." << endl << endl;
		throw exception();
	}
	
	cout << endl;
	
	Equation fEq(f,zeros,Frac(1,1),0); //Type doesn't matter
	vector<Edge> edges {};
	Graph H(edges,1,numColors);
	Equation eq1({H},zeros,Frac(1,1),0);
	Equation eq2({H},zeros,Frac(1,1),0);
	
	cout << "Generating graphs to be used in resize." << endl;
	vector<Graph> allGraphs = generate(n,numColors,zeros);
	
	cout << endl << "Resizing f." << endl;
	Equation fEqResized = resize(fEq,allGraphs);
	f.clear();
	for(int i = 0; i < fEqResized.getNumVariables(); ++i) {
		f.push_back(fEqResized.getVariable(i));
	}
	cout << endl;
	
	fSize = f.size();
	
	//Resize known
	cout << "Resizing known." << endl;
	for(int i = 0; i < knownSize; ++i) {
		cout << i+1 << " out of " << knownSize << endl;
		known[i] = resize(known[i],allGraphs);
	}
	cout << endl;
	
	//Calculates everything
	Graph G({},1,numColors);
	Equation eq({G},zeros,Frac(1,1),0);
	
	cout << "Creating everything equation." << endl;
	Equation everything (allGraphs,zeros,Frac(1,1),0);
	cout << endl;
	
	ofstream myFile;
	myFile.open("Duals.txt");
	
	for(int i = 0; i < everything.getNumVariables(); ++i) {
		myFile << "Label: " << i+1 << endl;
		everything.getVariable(i).printAdjMatToFile(myFile);
		myFile << endl;
	}

	
	cout << "Generating v." << endl;
	vector< vector< Graph> > v = generateV(n,numColors,zeros);
	cout << endl;

	queue<tuple<int,int,int,int,Frac> > A;

	vector<Frac> B; //Gives numbers to be printed
	vector< vector<Frac> > C; //From Known
	
	//Used to give indices of A
	unordered_map<string, int> everythingMap;
			
	for(int i = 0; i < everything.getNumVariables(); ++i) {
		everythingMap.insert(pair<string, int>(everything.getVariable(i).getCanonLabel(),i));
	}
	
	for(int i = 0; i < (int)v.size(); ++i) { //i is first index of A
		cout << "In A, iteration " << i+1 << " out of " << v.size() << endl;
		for(int j = 0; j < (int)v[i].size(); ++j) {
			for(int k = j; k < (int)v[i].size(); ++k) {
				
				Equation temp = multiply(v[i][j],v[i][k],zeros);
				temp.averageAll();
				
				for(int l = 0; l < temp.getNumVariables(); ++l) {
					if(temp.getVariable(l).getCoefficient() != 0) {
						if(j == k) {
							auto tempTup = make_tuple(i, everythingMap[temp.getVariable(l).getCanonLabel()], j ,k, temp.getVariable(l).getCoefficient());
							A.push(tempTup);
						}	
						
						else {
							auto tempTup = make_tuple(i, everythingMap[temp.getVariable(l).getCanonLabel()], j ,k, temp.getVariable(l).getCoefficient());
							A.push(tempTup);
						}
					}
				}
			}
		}
	}
	cout << endl;
	
	Frac zeroFrac(0,1);
	//Calculating B
	cout << "Calculating B." << endl << endl;
	for(int i = 0; i < everything.getNumVariables(); ++i) {
		B.push_back(zeroFrac);
	}
	
	for(int i = 0; i < fSize; ++i) {
		for(int j = 0; j < everything.getNumVariables(); ++j) {
			if(isomorphic(f[i],everything.getVariable(j))) {
				B[j] = f[i].getCoefficient();
				j = everything.getNumVariables();
			}
		}
	}
	
	//First has 0 for ==, 1 for <=
	//Next entry is bound
	//Finally, it's the vector of coefficients in known (after it's been resized)
	
	//Calculating C
	cout << "Calculating C." << endl << endl;
	C.resize(knownSize);
	for(int i = 0; i < knownSize; ++i) {
		C[i].resize(everything.getNumVariables()+1,zeroFrac);
	}
	
	for(int i = 0; i < knownSize; ++i) {
		C[i][0] = known[i].getAns();
		
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			for(int k = 0; k < everything.getNumVariables(); ++k) {
				if(isomorphic(known[i].getVariable(j),everything.getVariable(k))) {
					C[i][k+1] = known[i].getVariable(j).getCoefficient();
					k = everything.getNumVariables();
				}
			}
		}
	}
	
	//Non-integer version
	typedef std::numeric_limits<long double> ldbl;
	
	ofstream myFile1;
	myFile1.open("plainFlagAlgebra1.txt");
	myFile1.precision(ldbl::max_digits10);
	
	//Tells python the size of the matrix
	myFile1 << v.size() << " " << everything.getNumVariables() << " ";
	for(int i = 0; i < (int)v.size(); ++i) {
		myFile1 << v[i].size() << " ";
	}
	myFile1 << endl;
	
	while(!A.empty()) {	
		auto tempTup = A.front();
		A.pop();
		
		myFile1 << get<0>(tempTup) << " " << get<1>(tempTup)  << " " << get<2>(tempTup)  << " " << get<3>(tempTup)  << " " << (long double)get<4>(tempTup).getNum() / (long double)get<4>(tempTup).getDen() << endl;
	}
	
	ofstream myFile2;
	myFile2.open("plainFlagAlgebra2.txt");
	myFile2.precision(ldbl::max_digits10);
	
	for(int i = 0; i < everything.getNumVariables(); ++i) {
		myFile2 << (long double)B[i].getNum() / (long double)B[i].getDen() << " ";
	}

	ofstream myFile3;
	myFile3.open("plainFlagAlgebra3.txt");
	myFile3.precision(ldbl::max_digits10);
	
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < everything.getNumVariables()+1; ++j) {
			myFile3 << (long double)C[i][j].getNum() / (long double)C[i][j].getDen() << " ";
		}
		myFile3 << endl;
	}
	
	
	//Integer Version
	/*long long int mult = 1;
	
	queue<tuple<int,int,int,int,Frac> > Acopy = A;
	
	while(!Acopy.empty()) {
		auto tempTup = Acopy.front();
		Acopy.pop();

		mult = lcm(mult,get<4>(tempTup).getDen());
	}
	
	for(int i = 0; i < everything.getNumVariables(); ++i) {
		mult = lcm(mult,B[i].getDen());
	}
	
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < everything.getNumVariables()+1; ++j) {
			mult = lcm(mult,C[i][j].getDen());
		}
	}
	
	ofstream myFile1;
	myFile1.open("plainFlagAlgebra1.txt");
	
	//Tells python the size of the matrix
	myFile1 << v.size() << " " << everything.getNumVariables() << " ";
	for(int i = 0; i < (int)v.size(); ++i) {
		myFile1 << v[i].size() << " ";
	}
	myFile1 << endl;
	
	while(!A.empty()) {	
		auto tempTup = A.front();
		A.pop();
		
		myFile1 << get<0>(tempTup) << " " << get<1>(tempTup)  << " " << get<2>(tempTup)  << " " << get<3>(tempTup)  << " " << get<4>(tempTup).getNum()*mult/get<4>(tempTup).getDen() << endl;
	}
	
	ofstream myFile2;
	myFile2.open("plainFlagAlgebra2.txt");
	
	for(int i = 0; i < everything.getNumVariables(); ++i) {
		myFile2 << B[i].getNum()*mult/B[i].getDen() << " ";
	}

	ofstream myFile3;
	myFile3.open("plainFlagAlgebra3.txt");
	
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < everything.getNumVariables()+1; ++j) {
			myFile3 << C[i][j].getNum()*mult/C[i][j].getDen() << " ";
		}
		myFile3 << endl;
	}
	
	ofstream myFile4;
	myFile4.open("plainFlagAlgebra4.txt");
	myFile4 << mult;
	
	cout << "Finished plainFlagAlgebra." << endl << endl;*/
}





void NEWplainFlagAlgebra(vector<Graph> &f, int n, vector<Graph> &zeros, vector<Equation> &known) {
	cout << "Starting plainFlagAlgebra." << endl;
	
	//The way the python script is setup we need at least one known, this just says that edge density is <= 1
	if(known.size() == 0) {
		Graph K2({{}},2,f[0].getNumColors());
		Equation knownTemp({K2},zeros,Frac(1,1),1);
		known.push_back(knownTemp);
	}
	
	int fSize = f.size();
	int zerosSize = zeros.size();
	int knownSize = known.size();
		
	if(n <= 1) {
		cout << "Plain flag algebra method not set up for graphs with fewer than two vertices." << endl << endl;
		throw exception();
	}

	if(fSize == 0) {
		cout << "Need at least one graph in f in plainFlagAlgebra." << endl << endl;
		throw exception();
	}
	
	//Check if every equation in known has at least one variable
	for(int i = 0; i < knownSize; ++i) {
		if(known[i].getNumVariables() == 0) {
			cout << "Equations in known must have at least one graph in plainFlagAlgebra." << endl << endl;
			throw exception();
		}
	}
	
	//Make every equation in known type 1 (<=)
	for(int i = 0; i < knownSize; ++i) {
		if(known[i].getType() != 1) {
			vector<Graph> variablesTemp;
			
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
	}

	int fN = f[0].getN();
	
	//Graphs in f can't have flags
	for(int i = 0; i < fSize; ++i) {
		if(f[i].getFlag().getN() != 0) {
			cout << "All graphs in f in plainFlagAlgebra must not have any flags." << endl << endl;
			throw exception();
		}
	}

	//Graphs in known can't have flags
	//Maybe fix?
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			if(known[i].getVariable(j).getFlag().getN() != 0) {
				cout << "All graphs in known in plainFlagAlgebra must not have any flags." << endl << endl;
				throw exception();
			}
		}
	}
	
	//Make sure everything has correct number of colors
	int numColors = f[0].getNumColors();
	
	for(int i = 1; i < fSize; ++i) {
		if(f[i].getNumColors() != numColors) {
			cout << "Everything in plainFlagAlgebra in f must have the same number of colors." << endl << endl;
			throw exception();
		}
	}
	
	for(int i = 0; i < zerosSize; ++i) {
		if(zeros[i].getNumColors() != numColors) {
			cout << "Everything in plainFlagAlgebra in zeros must have the same number of colors." << endl << endl;
			throw exception();
		}
	}
	
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			if(known[i].getVariable(j).getNumColors() != numColors) {
				cout << "Everything in plainFlagAlgebra in known must have the same number of colors." << endl << endl;
				throw exception();
			}
		}
	}
	
	//Not sure if strictly necessary, but it would probably give trash bounds otherwise
	if(fN > n) {
		cout << "Make n large enough so it has vertices at least as many vertices when multiplied by itself as n." << endl << endl;
		throw exception();
	}
	
	cout << endl;
	
	Equation fEq(f,zeros,Frac(1,1),0); //Type doesn't matter
	vector<Edge> edges {};
	Graph H(edges,1,numColors);
	Equation eq1({H},zeros,Frac(1,1),0);
	Equation eq2({H},zeros,Frac(1,1),0);
	
	cout << "Generating graphs to be used in resize." << endl;
	vector<Graph> allGraphs = generate(n,numColors,zeros);
	
	cout << endl << "Resizing f." << endl;
	Equation fEqResized = resize(fEq,allGraphs);
	f.clear();
	for(int i = 0; i < fEqResized.getNumVariables(); ++i) {
		f.push_back(fEqResized.getVariable(i));
	}
	cout << endl;
	
	fSize = f.size();
	
	//Resize known
	cout << "Resizing known." << endl;
	for(int i = 0; i < knownSize; ++i) {
		cout << i+1 << " out of " << knownSize << endl;
		known[i] = resize(known[i],allGraphs);
	}
	cout << endl;

	
	ofstream myFile;
	myFile.open("Duals.txt");
	
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		myFile << "Label: " << i+1 << endl;
		allGraphs[i].printAdjMatToFile(myFile);
		myFile << endl;
	}

	queue<tuple<int,int,int,int,Frac> > A;
	vector<Frac> B; //Gives numbers to be printed
	vector< vector<Frac> > C; //From Known
	
	//Used to give indices of A
	//Rather than multiplying we look at all graphs and work backwards to get coefficients
	//This works well as we have already generated all graphs of size n and multiplication basically just has us partially doing that every time so it saves time on generation
	//The only slow down could come from looking up indices but unordered_maps are SO fast it doesn't matter
	unordered_map<string, int> allGraphsMap;
			
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		allGraphsMap.insert(pair<string, int>(allGraphs[i].getCanonLabel(),i));
	}
	
	//All graphs of size n with flags such that parity issues work out
	cout << "Generating v." << endl;
	vector < vector < Graph > > v = generateV(n,numColors,zeros);
	cout << endl;
	cout << "Generating allGraphsWithFlags." << endl;
	vector < vector < Graph > > allGraphsWithFlags = generateAllGraphsWithFlags(n,numColors,zeros);
	cout << endl;
	
	cout << "Calculating A." << endl << endl;
	
	for(int i = 0; i < (int)v.size(); ++i) { //i is first index of A
		cout << "In A, iteration " << i+1 << " out of " << v.size() << endl;  
		//Create hash map for third and fourth indices of A
		unordered_map<string, int> f;
		for(int j = 0; j < (int)v[i].size(); ++j) {
			f.insert(pair<string, int>(v[i][j].getCanonLabel(),j));
		}
		
		vector< vector < vector <Frac> > > partialA( (int)allGraphs.size(), vector< vector<Frac> > ((int)v[i].size(), vector<Frac>((int)v[i].size(), Frac(0,1))));
		
		int sizeOfFlag = v[i][0].getSizeOfFlag();
		
		for(int index2 = 0; index2 < (int)allGraphsWithFlags[i].size(); ++index2) {
			Graph G = allGraphsWithFlags[i][index2];
			
			//Make sure all flag vertices in G are first
			//TODO make a function of this
			vector<int> reordering(n);
			
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
			
			vector<int> X((n-sizeOfFlag)/2 - 1); //X is zero indexed but it should actually be sizeOfFlag+1 indexed
			//Wlog first non-flag vertex in X
			for(int j = 0; j < (n-sizeOfFlag)/2 - 1; ++j) {
				X[j] = j;
			}
			
			do {
				vector<int> restriction1(n,-1);
				vector<int> restriction2(n,-1);
				
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
				
				int c = min(G1Index,G2Index);
				int d = max(G1Index,G2Index);
				
				partialA[b][c][d] = partialA[b][c][d] +  e;
			
			} while(nextSubset(X,n-sizeOfFlag-1,(n-sizeOfFlag)/2 - 1));
		}
		
		//Combine partialA into A
		for(int index1 = 0; index1 < (int)partialA.size(); ++index1) {
			for(int index2 = 0; index2 < (int)partialA[index1].size(); ++index2) {
				for(int index3 = 0; index3 < (int)partialA[index1][index2].size(); ++index3) {
					if(partialA[index1][index2][index3] != 0) {
						if(index2 != index3) {
							auto tempTup = make_tuple(i,index1,index2,index3,partialA[index1][index2][index3]/2);
							A.push(tempTup);
						}
						
						else {
							auto tempTup = make_tuple(i,index1,index2,index3,partialA[index1][index2][index3]);
							A.push(tempTup);
						}
					}
				}
			}
		}
	}
	
	cout << endl;
	
	Frac zeroFrac(0,1);
	//Calculating B
	cout << "Calculating B." << endl << endl;
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		B.push_back(zeroFrac);
	}
	
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
	//Finally, it's the vector of coefficients in known (after it's been resized)
	
	//Calculating C
	cout << "Calculating C." << endl << endl;
	C.resize(knownSize);
	for(int i = 0; i < knownSize; ++i) {
		C[i].resize(allGraphs.size()+1,zeroFrac);
	}
	
	for(int i = 0; i < knownSize; ++i) {
		C[i][0] = known[i].getAns();
		
		for(int j = 0; j < known[i].getNumVariables(); ++j) {
			for(int k = 0; k < (int)allGraphs.size(); ++k) {
				if(isomorphic(known[i].getVariable(j),allGraphs[k])) {
					C[i][k+1] = known[i].getVariable(j).getCoefficient();
					k = (int)allGraphs.size();
				}
			}
		}
	}
	
	//Non-integer version
	typedef std::numeric_limits<long double> ldbl;
	
	ofstream myFile1;
	myFile1.open("plainFlagAlgebra1.txt");
	myFile1.precision(ldbl::max_digits10);
	
	//Tells python the size of the matrix
	myFile1 << v.size() << " " << allGraphs.size() << " ";
	for(int i = 0; i < (int)v.size(); ++i) {
		myFile1 << v[i].size() << " ";
	}
	myFile1 << endl;
	
	while(!A.empty()) {	
		auto tempTup = A.front();
		A.pop();
		
		myFile1 << get<0>(tempTup) << " " << get<1>(tempTup)  << " " << get<2>(tempTup)  << " " << get<3>(tempTup)  << " " << (long double)get<4>(tempTup).getNum() / (long double)get<4>(tempTup).getDen() << endl;
	}
	
	ofstream myFile2;
	myFile2.open("plainFlagAlgebra2.txt");
	myFile2.precision(ldbl::max_digits10);
	
	for(int i = 0; i < (int)allGraphs.size(); ++i) {
		myFile2 << (long double)B[i].getNum() / (long double)B[i].getDen() << " ";
	}

	ofstream myFile3;
	myFile3.open("plainFlagAlgebra3.txt");
	myFile3.precision(ldbl::max_digits10);
	
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < (int)allGraphs.size()+1; ++j) {
			myFile3 << (long double)C[i][j].getNum() / (long double)C[i][j].getDen() << " ";
		}
		myFile3 << endl;
	}
	
	
	//Integer Version
	/*long long int mult = 1;
	
	queue<tuple<int,int,int,int,Frac> > Acopy = A;
	
	while(!Acopy.empty()) {
		auto tempTup = Acopy.front();
		Acopy.pop();

		mult = lcm(mult,get<4>(tempTup).getDen());
	}
	
	for(int i = 0; i < everything.getNumVariables(); ++i) {
		mult = lcm(mult,B[i].getDen());
	}
	
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < everything.getNumVariables()+1; ++j) {
			mult = lcm(mult,C[i][j].getDen());
		}
	}
	
	ofstream myFile1;
	myFile1.open("plainFlagAlgebra1.txt");
	
	//Tells python the size of the matrix
	myFile1 << v.size() << " " << everything.getNumVariables() << " ";
	for(int i = 0; i < (int)v.size(); ++i) {
		myFile1 << v[i].size() << " ";
	}
	myFile1 << endl;
	
	while(!A.empty()) {	
		auto tempTup = A.front();
		A.pop();
		
		myFile1 << get<0>(tempTup) << " " << get<1>(tempTup)  << " " << get<2>(tempTup)  << " " << get<3>(tempTup)  << " " << get<4>(tempTup).getNum()*mult/get<4>(tempTup).getDen() << endl;
	}
	
	ofstream myFile2;
	myFile2.open("plainFlagAlgebra2.txt");
	
	for(int i = 0; i < everything.getNumVariables(); ++i) {
		myFile2 << B[i].getNum()*mult/B[i].getDen() << " ";
	}

	ofstream myFile3;
	myFile3.open("plainFlagAlgebra3.txt");
	
	for(int i = 0; i < knownSize; ++i) {
		for(int j = 0; j < everything.getNumVariables()+1; ++j) {
			myFile3 << C[i][j].getNum()*mult/C[i][j].getDen() << " ";
		}
		myFile3 << endl;
	}
	
	ofstream myFile4;
	myFile4.open("plainFlagAlgebra4.txt");
	myFile4 << mult;
	
	cout << "Finished plainFlagAlgebra." << endl << endl;*/
}


