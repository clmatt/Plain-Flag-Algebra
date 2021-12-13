//-------------------
//-------------------
//----Frac Class-----
//-------------------
//-------------------

//Numerator and denominator always positive
class Frac {
	
	private:
	
		long long int num;
		long long int den;
		
	public:
		
		//------------------------------
		//-----Constructor for Frac-----
		//------------------------------
		
		Frac(long long int a, long long int b) {
			if( b == 0) {
				cout << "Error in Frac constructor, no fractions with 0 denominator" << endl << endl;
				throw exception();
			}
			
			//Copied from reduce (need here because reduce works on already intialized fractions)			
			if(a == 0) {
				num = 0;
				den = 1;
			}
			
			else {
				long long c = gcd(a,b);
				
				if((a < 0) && (b < 0)) {
					a = -a;
					b = -b;
				}
				
				else if((a > 0) && (b < 0)) {
					a = -a;
					b = -b;
				}
				
				num = a/c;
				den = b/c;
			}
			
		}
		
		
		//-----------------
		//-----Reduce------
		//-----------------
		
		void reduce() {
			if(den == 0) {
				cout << "Error in reduce (shouldn't throw this error, so something is broken) 0 denominator" << endl << endl;
				throw exception();
			}
			
			if(num == 0) {
				den = 1;
			}
			
			else {
				int c = gcd(num,den);
				
				if((num < 0) && (den < 0)) {
					num = -num;
					den = -den;
				}
				
				else if((num > 0) && (den < 0)) {
					num = -num;
					den = -den;
				}
				
				num = num/c;
				den = den/c;
			}
		}
		
		
		//-----------------
		//-----Get Num-----
		//-----------------
		
		long long int getNum() {
			return num;
		}
		
		
		//-----------------
		//-----Set Num-----
		//-----------------
		
		void setNum(long long int a) {
			num = a;
			reduce();
		}
		
		//-----------------
		//-----Get Den-----
		//-----------------
		
		long long int getDen() {
			return den;
		}
		
		
		//-----------------
		//-----Set Den-----
		//-----------------
		
		void setDen(long long int b) {
			if(b == 0) {
				cout << "Can't set denominator to 0 in Frac." << endl << endl;
			}
			
			den = b;
			
			reduce();
		}
		
};

//---------------------------------------
//---------------------------------------
//-----Functions Involving Fractions-----
//---------------------------------------
//---------------------------------------


//--------------------
//-----Is Integer-----
//--------------------

bool isInteger(Frac frac) {
	if(frac.getDen() == 1) {
		return true;
	}
	
	return false;
}


//------------------------
//-----Print Fraction-----
//------------------------

ostream& operator<<(ostream& os, Frac frac) {
	if(isInteger(frac)) {
		os << frac.getNum();
	}

	else {
		os << frac.getNum() << "/" << frac.getDen();
	}
	
	return os;
}


//----------------------------
//-----Convert to integer-----
//----------------------------

long long int convertToInt(Frac frac) {
	if(isInteger(frac)) {
		return frac.getNum();
	}
	
	else {
		cout << "Can't convert fraction to integer if not an integer." << endl << endl;
		throw exception();
	}
}


//------------------------------
//-----Convert from Integer-----
//------------------------------

Frac convertFromInt(int i) {
	return Frac(i,1);
}


//------------------------
//-----Multiplication-----
//------------------------

Frac operator*(Frac frac1, Frac frac2) {
	long long int a = frac1.getNum() * frac2.getNum();
	long long int b = frac1.getDen() * frac2.getDen();
		
	return Frac(a,b);
}


//----------------------------------
//-----Multiplication by number-----
//----------------------------------

Frac operator*(long long int i, Frac const &frac1) {	
	return(frac1*convertFromInt(i));
}

Frac operator*(Frac const &frac1, long long int i) {
	return(frac1*convertFromInt(i));
}



//-----------------
//-----Inverse-----
//-----------------

Frac inv(Frac frac) {
	if(frac.getNum() == 0) {
		cout << "Error in inv, can't take inverse of 0." << endl << endl;
		throw exception();
	}
	
	return Frac(frac.getDen(),frac.getNum());
}


//------------------
//-----Division-----
//------------------

Frac operator/(Frac frac1, Frac frac2) {
	return(frac1 * (inv(frac2)));
}


//------------------------------
//-----Division by a number-----
//------------------------------

Frac operator/(Frac frac, int i) {
	if(i == 0) {
		cout << "Can't divide by 0" << endl << endl;
		throw exception();
	}
	
	return(frac/convertFromInt(i));
}

Frac operator/(int i, Frac frac) {
	return(convertFromInt(i)/frac);
}


//---------------
//-----Power-----
//---------------

Frac operator^(Frac frac, int k) {
	if(k == 0) {
		return Frac(1,1);
	}
	
	Frac ans = frac;
	
	if (k < 0) {
		ans = inv(frac);
		
		for(int i = 0; i < -k+1; ++i) {
			ans = ans*inv(frac);
		}
		
	}
	
	else{
		for(int i = 0; i < k-1; ++i) {
			ans = ans*frac;
		}
	}
	
	return ans;
}


//------------------
//-----Addition-----
//------------------

Frac operator+(Frac frac1, Frac frac2) {
	long long int a = frac1.getNum() * frac2.getDen() + frac1.getDen() * frac2.getNum();
	long long int b = frac1.getDen() * frac2.getDen();
	
	return Frac(a,b);
}


//------------------------------
//-----Addition with number-----
//------------------------------

Frac operator+(Frac frac, int i) {
	return frac+convertFromInt(i);
}

Frac operator+(int i, Frac frac) {
	return frac+convertFromInt(i);
}


//---------------------
//-----Subtraction-----
//---------------------

Frac operator-(Frac frac1, Frac frac2) {
	return(frac1 + (-1 * frac2));
}


//---------------------------------
//-----Subtraction with Number-----
//---------------------------------

Frac operator-(Frac frac, int i) {
	return frac-convertFromInt(i);
}

Frac operator-(int i, Frac frac) {
	return convertFromInt(i)-frac;
}


//-------------------
//-----Negation------
//-------------------

Frac operator-(Frac frac) {
	return( (-1) * frac);
}


//-----------------
//-----Add one-----
//-----------------

Frac operator++(Frac frac) {
	return frac+1;
}


//----------------------
//-----Subtract one-----
//----------------------

Frac operator--(Frac frac) {
	return frac-1;
}


//----------------
//-----Equal------
//----------------

bool operator==(Frac frac1, Frac frac2) {
	//All fractions are always in reduced form because of constructor
	if((frac1.getNum() == frac2.getNum()) && (frac1.getDen() == frac2.getDen())) {
		return true;
	}
	
	return false;
}

bool operator==(Frac frac1, int val) {
	Frac frac2 = convertFromInt(val);
	
	return(frac1 == frac2);
}

bool operator==(int val, Frac frac1) {
	Frac frac2 = convertFromInt(val);
	
	return(frac1 == frac2);
}


//--------------------
//-----Not Equal------
//--------------------

bool operator!=(Frac frac1, Frac frac2) {
	return !(frac1 == frac2);
}

bool operator!=(Frac frac1, int val) {
	Frac frac2 = convertFromInt(val);
	
	return (frac1 != frac2);
}

bool operator!=(int val, Frac frac1) {
	Frac frac2 = convertFromInt(val);
	
	return (frac1 != frac2);
}

//----------------------
//-----Greater than-----
//----------------------

bool operator>(Frac frac1, Frac frac2) {
	Frac frac3 = frac1 - frac2;
	
	if (frac3.getNum() > 0) {
		return true;
	}
	
	return false;
}

bool operator>(Frac frac1, int val) {
	Frac frac2 = convertFromInt(val);
	
	return (frac1 > frac2);
}

bool operator>(int val, Frac frac1) {
	Frac frac2 = convertFromInt(val);
	
	return (frac1 > frac2);
}

//----------------------------------
//-----Greater than or equal to-----
//----------------------------------

bool operator>=(Frac frac1, Frac frac2) {
	if((frac1 > frac2) || (frac1 == frac2)) {
		return true;
	}
	
	return false;
}

bool operator>=(Frac frac1, int val) {
	Frac frac2 = convertFromInt(val);
	
	return (frac1 >= frac2);
}

bool operator>=(int val, Frac frac1) {
	Frac frac2 = convertFromInt(val);
	
	return (frac1 >= frac2);
}

//-------------------
//-----Less than-----
//-------------------

bool operator<(Frac frac1, Frac frac2) {
	return(!(frac1 >= frac2));
}

bool operator<(Frac frac1, int val) {
	Frac frac2 = convertFromInt(val);
	
	return (frac1 < frac2);
}

bool operator<(int val, Frac frac1) {
	Frac frac2 = convertFromInt(val);
	
	return (frac1 < frac2);
}

//-------------------------------
//-----Less than or equal to-----
//-------------------------------

bool operator<=(Frac frac1, Frac frac2) {
	return(!(frac1 > frac2)); 
}

bool operator<=(Frac frac1, int val) {
	Frac frac2 = convertFromInt(val);
	
	return (frac1 <= frac2);
}

bool operator<=(int val, Frac frac1) {
	Frac frac2 = convertFromInt(val);
	
	return (frac1 <= frac2);
}
