// Header file with Vector class

#ifndef vec_h
#define vec_h

#include <math.h>

using namespace std;

class dVector
{
private:
	int no_row;
	double *v;

public:
	// Constructors & Destructor
	dVector(void);
	dVector(const dVector & b);
	dVector(int nr);
	~dVector(void);

	// Operators
	dVector& operator=(const dVector & b);
	dVector& operator=(double b);
	dVector& operator=(int b);
	dVector operator+=(const dVector & b);
	dVector operator-=(const dVector & b);
	dVector operator+(const dVector & b);
	dVector operator-(const dVector & b);
	dVector operator*(double c);
	dVector operator*(int c);
	double& operator[](int i); // Index operator
	double getTheta(dVector a, dVector b); // angle between 2 2D vectors

	// Member functions
	int size();
	double length(void);

	// Friend functions
	friend dVector operator*(double c, dVector b);
	friend dVector operator*(int c, dVector b);
	friend double length(const dVector & b);
	friend double dot(const dVector & a, const dVector & b);
	friend dVector cross(const dVector & a, const dVector & b);
	friend ostream & operator << (ostream & os, const dVector & b);
	friend istream & operator >> (istream & is, const dVector & b);
};
dVector::dVector(void) // Default constructor
{
	no_row = 0;
	v = 0;
}
dVector::dVector(const dVector & b) // Copy constructor
{
	v = new double [no_row = b.no_row];
	v --;
	for (int i=1; i<=no_row; i++)
		v[i] = b.v[i];
}
dVector::dVector(int nr) // constructor
{
	v = new double [no_row = nr];
	v --;
	for (int i=1; i<=no_row; i++)
		v[i] = 0;
}
dVector::~dVector() // Destructor
{
	if (++v) delete[] v;
}
int dVector::size()
{
	return no_row;
}
double dVector::length() // Member method
{
	double l = 0;
	for (int i=1; i<=no_row; i++)
		l += v[i]*v[i];
	return sqrt(l);
}
double length (const dVector & b) // Friend method
{
	double l = 0;
	for (int i=1; i<=b.no_row; i++)
		l += b.v[i]*b.v[i];
	return sqrt(l);
}
dVector dVector::operator*(double c)
{
	for (int i=1; i<=no_row; i++)
		v[i] *= c;
	return *this;
}
dVector operator*(double c, dVector b)
{
	return b*c;
}
dVector dVector::operator*(int c)
{
	for (int i=1; i<=no_row; i++)
		v[i] *= c;
	return *this;
}
dVector operator*(int c, dVector b)
{
	return b*c;
}

dVector& dVector::operator=(const dVector & b) // Assignment operator
{
	if (v == b.v)
		return *this;
	v ++;
//	delete v;
	v = new double [no_row = b.no_row];
	v --;
	for (int i=1; i<=b.no_row; i++)
		v[i] = b.v[i];
	return *this;
}
dVector& dVector::operator=(double b) // Assignment operator
{
	for (int i=1; i<=no_row; i++)
		v[i] = b;
	return *this;
}
dVector& dVector::operator=(int b) // Assignment operator
{
	for (int i=1; i<=no_row; i++)
		v[i] = b;
	return *this;
}
double & dVector::operator[](int i) // Index operator
{
	if (v && i>=1 && i<=no_row)
		return v[i];
	else 
	{
		cout << " error in index in Vector: i, Vi are :   " << i << "  " << v[i];// ***   check these lines later
		return v[i];
	}
}
dVector dVector::operator+=(const dVector & b) // Addition with assign operator
{
	if (no_row == b.no_row)
		for (int i=1; i<=no_row; i++)
			v[i] += b.v[i];
	return *this;
}
dVector dVector::operator-=(const dVector & b) // Addition with assign operator
{
	if (no_row == b.no_row)
		for (int i=1; i<=no_row; i++)
			v[i] -= b.v[i];
	return *this;
}
dVector dVector::operator+(const dVector & b) // Addition operator
{
	dVector a(*this);
	a += b;
	return a;
}
dVector dVector::operator-(const dVector & b) // Addition operator
{
	dVector a(*this);
	a -= b;
	return a;
}
double dot(const dVector & a, const dVector & b)
{
	double prod = 0;
	if (a.no_row != b.no_row)
		cout << "\n*** error in sizes of dVectors a and b  ***";
	for (int i=1; i<=a.no_row; i++)
		prod += a.v[i] * b.v[i];
	return prod;
}
double getTheta(dVector a, dVector b)
{
	double theta;
	a = a*(1/length(a));
	b = b*(1/length(b));
	theta = asin(a[1]*b[2] - a[2]*b[1]);
	return theta;
}
dVector cross(const dVector & a, const dVector & b)
{
	dVector c(a.no_row);
	c.v[1] = a.v[2]*b.v[3] - a.v[3]*b.v[2];
	c.v[2] = a.v[3]*b.v[1] - a.v[1]*b.v[3];
	c.v[3] = a.v[1]*b.v[2] - a.v[2]*b.v[1];
	return c;
}
ostream& operator<<(ostream & os, const dVector & b)
{
	os.precision(4);
	os.setf(ios::showpoint);
	os << endl;
	for (int i=1; i<=b.no_row; i++)
		os << i << "   " << b.v[i] << endl;
	os << endl;
	return os;
}
istream& operator>>(istream& is, const dVector& b)
{
	for (int i=1; i<=b.no_row; i++)
		is >> b.v[i];
	return is;
}

class iVector
{
private:
	int no_row;
	int *v;

public:
	// Constructors & Destructor
	iVector(void);
	iVector(const iVector & b);
	iVector(int nr);
	~iVector(void);

	// Operators
	iVector& operator=(const iVector & b);
	iVector operator+=(const iVector & b);
	iVector operator+(const iVector & b);
	int & iVector::operator[](int i); // Index operator

	// Member functions
	int size();
	double length(void);

	// Friend functions
	friend double length(const iVector & b);
	friend int dot(const iVector & a, const iVector & b);
	friend iVector cross(const iVector & a, const iVector & b);
	friend ostream& operator<<(ostream & os, const iVector & b);
	friend istream& operator>>(istream & is, const iVector & b);
};
iVector::iVector(void) // Default constructor
{
	no_row = 0;
	v = 0;
}
iVector::iVector(const iVector & b) // Copy constructor
{
	v = new int [no_row = b.no_row];
	v --;
	for (int i=1; i<=no_row; i++)
		v[i] = b.v[i];
}
iVector::iVector(int nr) // constructor
{
	v = new int [no_row = nr];
	v --;
	for (int i=1; i<=no_row; i++)
		v[i] = 0;
}
iVector::~iVector() // Destructor
{
//	if (++v) delete[] v;
}
double iVector::length() // Member method
{
	double l = 0;
	for (int i=1; i<=no_row; i++)
		l += v[i]*v[i];
	return sqrt(l);
}
int iVector::size() // Member method
{
	return no_row;
}

double length (const iVector & b) // Friend method
{
	double l = 0;
	for (int i=1; i<=b.no_row; i++)
		l += b.v[i]*b.v[i];
	return sqrt(l);
}
iVector& iVector::operator=(const iVector & b) // Assignment operator
{
	if (v == b.v)
		return *this;
	v ++;
	//delete [] v;
	v = new int [no_row = b.no_row];
	v --;
	for (int i=1; i<=b.no_row; i++)
		v[i] = b.v[i];
	return *this;
}
int & iVector::operator[](int i) // Index operator
{
	if (v && i>=1 && i<=no_row)
		return v[i];
	else 
	{
		cout << " error in index ???";// ***   check these lines later
		return v[i];
	}
}
iVector iVector::operator+=(const iVector & b) // Addition with assign operator
{
	if (no_row == b.no_row)
		for (int i=1; i<=no_row; i++)
			v[i] += b.v[i];
	return *this;
}
iVector iVector::operator+(const iVector & b) // Addition operator
{
	iVector a(*this);
	a += b;
	return a;
}
int dot(const iVector & a, const iVector & b)
{
	int prod = 0;
	if (a.no_row != b.no_row)
		cout << "\n*** error in sizes of iVectors a and b  ***";
	for (int i=1; i<=a.no_row; i++)
		prod += a.v[i] * b.v[i];
	return prod;
}
iVector cross(const iVector & a, const iVector & b)
{
	iVector c(a.no_row);
	c.v[1] = a.v[2]*b.v[3] - a.v[3]*b.v[2];
	c.v[2] = a.v[3]*b.v[1] - a.v[1]*b.v[3];
	c.v[3] = a.v[1]*b.v[2] - a.v[2]*b.v[1];
	return c;
}
ostream& operator<<(ostream & os, const iVector & b)
{
	os.precision(4);
	os.setf(ios::showpoint);
	os << endl;
	for (int i=1; i<=b.no_row; i++)
		os << b.v[i] << endl;
	os << endl;
	return os;
}
istream& operator>>(istream& is, const iVector& b)
{
	for (int i=1; i<=b.no_row; i++)
		is >> b.v[i];
	return is;
}
#endif