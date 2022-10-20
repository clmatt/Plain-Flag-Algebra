//-----------------------
//-----------------------
//-----Simplex Class-----
//-----------------------
//-----------------------


class Simplex {

	private:
		
		int n = -1; //dimension - assume n points in n dimension (simplicial dimension is n-1)
		std::vector< std::vector <double> > vertices; //List of vertices

	public:
		
		//---------------------------------
		//-----Constructor for Simplex-----
		//---------------------------------
		
		//Need to specify NUMCOLORS, because it's possible that it differs
		//from actual number of colors
		Simplex(const std::vector< std::vector<double> >  &VERTICES) {
			n = VERTICES.size();
			
			for(int i = 0; i < n; ++i) {
				if(VERTICES[i].size() != n) {
					std::cout << "Coordinates in simplex aren't of correct size." << std::endl <<std::endl;
					throw std::exception();
				}
			}
			
			vertices = VERTICES;
			
			//TODO check if valid - need a linear algebra solver - don't need yet however
		}
		
		
		//-----------------------
		//-----Get Dimension-----
		//-----------------------
		
		int getDim() {
			return (int)vertices.size();
		}
		
		
		//----------------------
		//-----Get Vertices-----
		//----------------------
		
		std::vector< std::vector<double> > getVertices() {
			return vertices;
		}
		
		
		//--------------------
		//-----Get Vertex-----
		//--------------------
				
		std::vector<double> getVertex(int i) {
			if((i < 0) || (i >= n)) {
				std::cout << "Vertex out of range in getVertex." << std::endl << std::endl;
				throw std::exception();
			}
			
			return vertices[i];
		}
		
		
		//------------------------
		//-----Get Coordinate-----
		//------------------------
		
		double getCoordinate(int i, int j) {
			if((i < 0) || (i >= n) || (j < 0) || (j >= n)) {
				std::cout << "Vertex out of range in getCoordinate." << std::endl << std::endl;
				throw std::exception();
			}
			
			return vertices[i][j];
		}
		
		
		//----------------------
		//-----Longest Edge-----
		//----------------------
		
		//Uses l1 norm for ease of computation and b/c it suffices for algorithm
		//Used in Copositive program solver
		std::pair<int, int> longestEdge() {
			double length = 0.;
			std::pair<int, int> output;
			
			for(int i = 0; i < n-1; ++i) {
				for(int j = i+1; j < n; ++j) {
					double temp = 0.;
					
					for(int k = 0; k < n; ++k) {
						temp = temp + abs(vertices[i][k] - vertices[j][k]);
					}
					
					if(temp > length) {
						length = temp;
						output.first = i;
						output.second = j;
					}
				}
			}
			
			return output;
		}	
		
};


//--------------------------
//-----Standard Simplex-----
//--------------------------

Simplex standardSimplex(const int n) {
	std::vector< std::vector<double> > vertices;
	vertices.resize(n);
	
	for(int i = 0; i < n; ++i) {
		vertices[i].resize(n,0.);
		vertices[i][i] = 1.;
	}
	
	return Simplex(vertices);
}


//-----------------------------
//-----------------------------
//-----Simplicial Complex------
//-----------------------------
//-----------------------------

//Assumes all simplices are the same dimension
class SimplicialComplex {

	private:
		
		std::vector<Simplex> simplices; //List of simplices
		std::unordered_map< std::vector<double>, std::vector< std::vector<int> >, containerHash< std::vector<double> > > verticesMap; //Maps to which simplices contain the vertex
		std::unordered_map< std::vector<double>, std::vector< std::vector<int> >, containerHash< std::vector<double> > > edgesMap; //Concatenate two vertices so I don't have to deal with pairs for the has function
		//Could make a simplex edge struct, but that seems like overkill
		int n;

	public:
		
		//--------------------------------------------
		//-----Constructor for Simplicial Complex-----
		//--------------------------------------------
		
		//Need to specify NUMCOLORS, because it's possible that it differs
		//from actual number of colors
		SimplicialComplex(const std::vector<Simplex> SIMPLICES) {
			simplices = SIMPLICES;
		
			if(simplices.size() > 0) {
				n = simplices[0].getDim();
				
				for(int i = 1; i < (int)simplices.size(); ++i) {
					if(simplices[i].getDim() != n) {
						std::cout << "In simplicial complex, all simplices must be same dimension." << std::endl << std::endl;
						throw std::exception();
					}
				}
				
				//Add to vertices set
				for(int i = 0; i < (int)simplices.size(); ++i) {
					for(int j = 0; j < n; ++j) {
						if(verticesMap.find(simplices[i].getVertex(j)) == verticesMap.end()) {
							verticesMap[simplices[i].getVertex(j)] = {{i,j}};
						}
					
						else {
							verticesMap[simplices[i].getVertex(j)].push_back({i,j});
						}
					}
				}
				
				//Add to edges set
				for(int i = 0; i < (int)simplices.size(); ++i) {
					for(int j = 0; j < n; ++j) {
						for(int k = j+1; k < n; ++k) {
							//Concatenate two endpoint of edges
							std::vector<double> edge = simplices[i].getVertex(j);
							edge.insert(edge.end(), simplices[i].getVertex(k).begin(), simplices[i].getVertex(k).end());
							
							if(edgesMap.find(edge) == edgesMap.end()) {
								edgesMap[edge] = {{i,j,k}};
							}
							else {
								edgesMap[edge].push_back({i,j,k});
							}
							
							//Include both ways to encode edges
							edge = simplices[i].getVertex(k);
							edge.insert(edge.end(), simplices[i].getVertex(j).begin(), simplices[i].getVertex(j).end());
							
							if(edgesMap.find(edge) == edgesMap.end()) {
								edgesMap[edge] = {{i,k,j}};
							}
							else {
								edgesMap[edge].push_back({i,k,j});
							}
						}
					}
				}
			}
		}
		
		
		//---------------------
		//-----Get Simplex-----
		//---------------------
		
		Simplex getSimplex(int i) {
			if((i < 0) || (i > simplices.size())) {
				std::cout << "In getSimplex, index is out of range." << std::endl << std::endl;
				throw std::exception();
			}	
			
			return simplices[i];
		}
		
		
		//---------------------
		//-----Add simplex-----
		//---------------------
		
		void addSimplex(Simplex delta) {
			if(delta.getDim() != n) {
				std::cout << "In Simplicial Complex, all simplices must be same dimension." << std::endl << std::endl;
				throw std::exception();
			}
			
			simplices.push_back(delta);
		}
		
		
		//-----------------------
		//-----Get Dimension-----
		//-----------------------
		
		int getDim() {
			return n;
		}
		
		
		//--------------------------
		//-----Get Vertices Map-----
		//--------------------------
		
		std::unordered_map< std::vector<double>, std::vector< std::vector<int> >, containerHash< std::vector<double> > > getVerticesMap() {
			return verticesMap;
		}
		
		
		//--------------------------
		//-----Get Edges Map-----
		//--------------------------
		
		std::unordered_map< std::vector<double>, std::vector< std::vector<int> >, containerHash< std::vector<double> > > getEdgesMap() {
			return edgesMap;
		}
		
		
		//-------------------
		//-----Subdivide-----
		//-------------------
		
		//Returns the vertices which with the new vertex are in a new edge - used for copositive approximation
		std::unordered_set< std::vector<double>, containerHash<std::vector<double> > > subdivide(const std::vector<double> &edge) {
			std::unordered_set< std::vector<double>, containerHash<std::vector<double> > > output;
		
			if(edgesMap.find(edge) == edgesMap.end()) {
				std::cout << "In subdivide in simplicial complex, the edge you are subdividing must exist." << std::endl << std::endl;
				throw std::exception();
			}
			
			std::vector<double> newVertex(n);
			for(int i = 0; i < n; ++i) {
				newVertex[i] = (edge[i] + edge[i+n])/2.;
			}
			
			std::vector< std::vector<int> > indices = edgesMap.at(edge);
			
			for(int i = 0; i < (int)indices.size(); ++i) {
				std::vector< std::vector<double> > vertices1 = simplices[indices[i][0]].getVertices();
				std::vector< std::vector<double> > vertices2 = simplices[indices[i][0]].getVertices();
				
				vertices1[indices[i][1]] = newVertex;
				vertices2[indices[i][2]] = newVertex;
				
				simplices.push_back(Simplex(vertices1));
				simplices.push_back(Simplex(vertices2));
			}
			
			for(int i = 0; i < (int)indices.size(); ++i) {
				for(int j = 0; j < n; ++j) {
					output.insert(simplices[indices[i][0]].getVertex(j));
				}
			}
			
			for(int i = 0; i < (int)indices.size(); ++i) {
				simplices.erase(simplices.begin() + indices[i][0]-i); //Okay because we push_back and b/c indices[i][0] is ordered
			}
			
			return output;
		}
		
		
		//----------------------
		//-----Longest Edge-----
		//----------------------
		
		//L1 norm
		std::vector<double> longestEdge() {
			std::vector<double> output;
			double length = 0.;
		
			for (auto edge : edgesMap) {
				double tempLength = 0.;
				
				for(int i = 0; i < n; ++n) {
					tempLength += abs(edge.first[i] - edge.first[i+n]);
				}
				
				if(tempLength > length) {
					output = edge.first;
					length = tempLength;
				}
			}
			
			return output;
		}
		
};



