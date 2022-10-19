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

