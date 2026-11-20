//-----------------------
//-----------------------
//-----Simplex Class-----
//-----------------------
//-----------------------


class Simplex {

	private:
		
		int n = -1; //dimension - assume n points in n dimension (simplicial dimension is n-1)
		std::vector< std::vector <double> > vertices; //List of vertices
		std::unordered_set< std::vector<double>, vectorHash< double > > vertexSet;

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
			
			for(int i = 0; i < n; ++i) {
				vertexSet.insert(vertices[i]); 
			}
			
			//TODO check if valid - need a linear algebra solver - don't need yet however
		}
		
		
		//-----------------------
		//-----Get Dimension-----
		//-----------------------
		
		int getDim() const{
			return (int)vertices.size();
		}
		
		
		//----------------------
		//-----Get Vertices-----
		//----------------------
		
		std::vector< std::vector<double> > getVertices() const{
			return vertices;
		}
		
		
		//--------------------
		//-----Get Vertex-----
		//--------------------
				
		std::vector<double> getVertex(int i) const{
			if((i < 0) || (i >= n)) {
				std::cout << "Vertex out of range in getVertex." << std::endl << std::endl;
				throw std::exception();
			}
			
			return vertices[i];
		}
		
		
		//------------------------
		//-----Get Coordinate-----
		//------------------------
		
		double getCoordinate(int i, int j) const{
			if((i < 0) || (i >= n) || (j < 0) || (j >= n)) {
				std::cout << "Vertex out of range in getCoordinate." << std::endl << std::endl;
				throw std::exception();
			}
			
			return vertices[i][j];
		}
		
		
		//------------------------
		//-----Get Vertex Set-----
		//------------------------
		
		std::unordered_set< std::vector<double>, vectorHash< double > > getVertexSet() {
			return vertexSet;
		}
		
		
		//-----------------------
		//-----Contains Edge-----
		//-----------------------
		
		bool containsEdge(const std::vector<double> myEdge) {
			std::vector<double> vertex1(myEdge.begin(), myEdge.begin()+n);
			std::vector<double> vertex2(myEdge.begin()+n, myEdge.begin()+2*n);
			
			if(vertex1 == vertex2) {
				return false;
			}
			
			if((vertexSet.find(vertex1) != vertexSet.end()) && (vertexSet.find(vertex2) != vertexSet.end())) {
				return true;
			}
			
			return false;
		}
		
		
		//----------------------
		//-----Longest Edge-----
		//----------------------
		
		//Uses l1 norm for ease of computation and b/c it suffices for algorithm
		//Can be used in Copositive program solver
		std::pair<int, int> longestEdge() const{
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
		//Could make a simplex edge struct, but that seems like overkill
		int n;

	public:
		
		//--------------------------------------------
		//-----Constructor for Simplicial Complex-----
		//--------------------------------------------
		
		//Need to specify NUMCOLORS, because it's possible that it differs
		//from actual number of colors
		SimplicialComplex(const std::vector<Simplex> &SIMPLICES) {
			simplices = SIMPLICES;
		
			if(simplices.size() > 0) {
				n = simplices[0].getDim();
				
				for(int i = 1; i < (int)simplices.size(); ++i) {
					if(simplices[i].getDim() != n) {
						std::cout << "In simplicial complex, all simplices must be same dimension." << std::endl << std::endl;
						throw std::exception();
					}
				}
			}
		}
		
		
		//-----------------------
		//-----Get Simplices-----
		//-----------------------
		
		std::vector<Simplex> getSimplices() const{
			return simplices;
		}
		
		
		//---------------------
		//-----Get Simplex-----
		//---------------------
		
		Simplex getSimplex(int i) const{
			if((i < 0) || (i > simplices.size())) {
				std::cout << "In getSimplex, index is out of range." << std::endl << std::endl;
				throw std::exception();
			}	
			
			return simplices[i];
		}
		
		
		//-----------------------------
		//-----Number of Simplices-----
		//-----------------------------
		
		int numSimplices() const{
			return simplices.size();
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
		
		int getDim() const{
			return n;
		}
		
		
		//----------------------
		//-----Longest Edge-----
		//----------------------
		
		//L1 norm
		std::vector<double> longestEdge() {
			std::vector<double> output(2*n);
			double length = -1.;
			
			
			for(int i = 0; i < simplices.size(); ++i) {
				for(int j = 0; j < n; ++j) {
					for(int k = j+1; k < n; ++k) {
						double tempLength = 0.;
					
						for(int l = 0; l < n; ++l) {
							tempLength += abs(simplices[i].getCoordinate(j,l) - simplices[i].getCoordinate(k,l));
						}
						
						if(tempLength > length) {
							length = tempLength;
							for(int l = 0; l < n; ++l) {
								output[l] = simplices[i].getCoordinate(j,l);
								output[l+k] = simplices[i].getCoordinate(k,l);
							}
						}
					}
				}
			}
			
			return output;
		}
		
		
		//-------------------
		//-----Subdivide-----
		//-------------------
				
		//Output is edges in subdivided simplex - used in copositive 
		//Maybe not the fastest b/c need to recompute maps, but I don't think this is a bottleneck
		//Using a map is a bit hacky, but it's fine
		void subdivide(const std::vector<double> &myEdge, std::unordered_map<std::vector<double>, std::vector<double>, vectorHash<double> >  &output) {	

			std::vector<int> indices; //Gives which simplices the edge is in
			
			for(int i = 0; i < simplices.size(); ++i) {
				if(simplices[i].containsEdge(myEdge)) {
					indices.push_back(i);
				}
			}

			if(indices.size() == 0) {
				std::cout << "In subdivide in simplicial complex, the edge you are subdividing must exist." << std::endl << std::endl;
				throw std::exception();
			}
					
			std::vector<double> newVertex(n);
			std::vector<double> oldVertex1(n);
			std::vector<double> oldVertex2(n);
			
			for(int i = 0; i < n; ++i) {
				newVertex[i] = (myEdge[i] + myEdge[i+n])/2.;
				oldVertex1[i] = myEdge[i];
				oldVertex2[i] = myEdge[i+n];
			}
			
			for(int i = 0; i < (int)indices.size(); ++i) {
				std::vector< std::vector<double> > vertices = simplices[indices[i]].getVertices();
				int index1 = -1;
				int index2 = -1;
				
				for(int j = 0; j < n; ++j) {
					if((vertices[j] != oldVertex1) && (vertices[j] != oldVertex2)) {
						std::vector<double> myEdge1(2*n);
						std::vector<double> myEdge2(2*n);
						std::vector<double> myEdge3(2*n);
						
						for(int k = 0; k < n; ++k) {
							myEdge1[k] = vertices[j][k];
							myEdge2[k] = vertices[j][k];
							myEdge3[k] = vertices[j][k];
							myEdge1[k+n] = oldVertex1[k];
							myEdge2[k+n] = oldVertex2[k];
							myEdge3[k+n] = newVertex[k];
						}
						
						output[myEdge1] = myEdge3;
						output[myEdge2] = myEdge3;
					} 
					
					else {
						if(index1 == -1) {
							index1 = j;
						}
						
						else {
							index2 = j;
						}
					}
				}
				
				std::vector< std::vector<double> > vertices1 = vertices;
				std::vector< std::vector<double> > vertices2 = vertices;		
						
				vertices1[index1] = newVertex;
				vertices2[index2] = newVertex;
						
				simplices.push_back(Simplex(vertices1));
				simplices.push_back(Simplex(vertices2));
			}
					
			for(int i = 0; i < (int)indices.size(); ++i) {
				simplices.erase(simplices.begin() + indices[i] - i); //Okay because we push_back and b/c indices[i][0] is ordered
			}
			
			
		}
};






