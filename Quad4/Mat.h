// Header file with Matrix class
#ifndef mat_h
#define mat_h

#include <fstream>
#include <iostream>
#include <math.h>
#include "Vec.h"

using namespace std;

double const PI = 3.1415926535897932384626433832795;
int min(int a, int b);
int i, j, k;


class dMatrix // double Matrix
{
private:
	int no_row;
	int no_col;
	double **m;

public:
	// Constructors & Destructor
	dMatrix(void);
	dMatrix(const dMatrix& b);
	dMatrix(int nr, int nc);
	~dMatrix(void);

	dVector solve(dVector b); // gauss elimination 
							  // (fully populated symmetric matrix)
	dMatrix invert(); // invert square symmetric matrix
	
	// Operators
	dMatrix operator=(const dMatrix& b);
	dMatrix operator=(double b);
	dMatrix operator=(int b);
	dMatrix operator+=(const dMatrix& b);
	dMatrix operator+(const dMatrix& b);
	dVector operator*(dVector b); // multiplication by a vector
	dVector operator^(dVector b); // gauss elimination 
								  // with banded storage

	double& operator()(int i, int j); // Index operator
	dMatrix operator*(double b);
	dMatrix operator*(int b);
	dMatrix operator*(dMatrix b); // multiplication by a matrix
	dMatrix transpose();

	// Friend functions
	friend dMatrix operator*(double b, dMatrix& a);
	friend dMatrix operator*(int b, dMatrix& a);
	friend ostream& operator<<(ostream & os, const dMatrix & b);
	friend istream& operator>>(istream & is, const dMatrix & b);
};
dMatrix::dMatrix(void) // Default constructor
{
	no_row = 0;
	no_col = 0;
	m = 0;
}
dMatrix::dMatrix(const dMatrix & b) // Copy constructor
{
	m = new double *[no_row = b.no_row];
	m --;
	for (i=1; i<=no_row; i++)
	{
		m[i] = new double [no_col = b.no_col];
		m[i] --;
	}
	for (i=1; i<=no_row; i++)
		for (j=1; j<=no_col; j++)
			m[i][j] = b.m[i][j];
}
dMatrix::dMatrix(int nr, int nc) // constructor
{
	m = new double *[no_row = nr];
	m --;
	for (i=1; i<=no_row; i++)
	{
		m[i] = new double [no_col = nc];
		m[i] --;
	}
	for (i=1; i<=no_row; i++)
		for (j=1; j<=no_col; j++)
			m[i][j] = 0.0;
}
dMatrix::~dMatrix() // Destructor
{
	for (i=1; i<=no_row; i++)
		if (++m[i])
			delete m[i];
}
dMatrix dMatrix::operator=(const dMatrix & b) // Assignment operator
{
	if (m == b.m)
		return *this;
	m ++;
	for (i=1; i<=no_row; i++)
		m[i] --;
//	delete [][]m;
	m = new double *[no_row = b.no_row];
	m --;
	for (i=1; i<=no_row; i++)
	{
		m[i] = new double [no_col = b.no_col];
		m[i] --;
	}
	for (i=1; i<=no_row; i++)
		for (j=1; j<=no_col; j++)
			m[i][j] = b.m[i][j];
	return *this;
}

dMatrix dMatrix::operator=(double b) // Assignment operator
{
	for (i=1; i<=no_row; i++)
		for (j=1; j<=no_col; j++)
			m[i][j] = b;
	return *this;
}
dMatrix dMatrix::operator=(int b) // Assignment operator
{
	for (i=1; i<=no_row; i++)
		for (j=1; j<=no_col; j++)
			m[i][j] = b;
	return *this;
}

dMatrix dMatrix::transpose()
{
	int p = no_col;
	int q = no_row;
	dMatrix A(p,q);

	for (i=1; i<=no_row; i++)
		for (j=1; j<=no_col; j++)
			A.m[j][i] = m[i][j];

	return A;
}


double& dMatrix::operator()(int i, int j) // Index operator
{
	if (m && i>=1 && i<=no_row && j>=1 && j<=no_col)
		return m[i][j];
	else 
	{
		cout << " error in index";
		return m[i][j];
	}
}
dMatrix dMatrix::operator+=(const dMatrix & b) // Addition with assign operator
{
	for (i=1; i<=no_row; i++)
		for (j=1; j<=no_col; j++)
			m[i][j] += b.m[i][j];
	return *this;
}
dMatrix dMatrix::operator+(const dMatrix & b) // Addition operator
{
	dMatrix a(*this);
	a += b;
	return a;
}
dMatrix dMatrix::operator*(double b)
{
	dMatrix c(*this);
	for (i=1; i<=no_row; i++)
		for (j=1; j<=no_col; j++)
			c.m[i][j] = m[i][j]*b;
	return c;
}
dMatrix operator*(double b, dMatrix& a) // multiplication by a constant
{
	return a*b;
}
dMatrix dMatrix::operator*(int b)
{
	for (i=1; i<=no_row; i++)
		for (j=1; j<=no_col; j++)
			m[i][j] *= b;
	return *this;
}
dMatrix operator*(int b, dMatrix& a) // multiplication by a constant
{
	return a*b;
}
dMatrix dMatrix::operator*(dMatrix b) // multiplication by a constant
{
	dMatrix c(no_row,b.no_col);
	for (i=1; i<=c.no_row; i++)
		for (j=1; j<=c.no_col; j++)
		{
			c.m[i][j] = 0;
			for (k=1; k<=no_col; k++)
				c.m[i][j] += m[i][k]*b.m[k][j];
		}
	return c;
}
dVector dMatrix::operator*(dVector b)
{
	dVector c(b);
	for (i=1; i<=no_row; i++)
	{
		c[i] = 0;
		for (j=1; j<=no_col; j++)
			c[i] += m[i][j]*b[j];
	}
	return c;
}
dMatrix dMatrix::invert()
{
	double r;
	int i, j, k, p;
	int neq = no_row;
	dMatrix b(neq,neq);
	for (i=1; i<=neq; i++)
		b(i,i) = 1;

// Forward reduction
	for (k=2; k<=neq; ++k)
	{
		for (i=k; i<=neq; ++i)
		{
			r = m[i][k-1]/m[k-1][k-1];
			for (p=1; p<=neq; p++)
				b(p,i) -= r*b(p,k-1);
			for (j=k; j<=neq; ++j)
				m[i][j] -= r*m[k-1][j];
		}
	}
// Backsubstitution phase
	for (p=1; p<=neq; p++)
		b(p,neq) = b(p,neq)/m[neq][neq];
	for (k=neq-1; k>=1; --k)
	{
		for (j=k+1; j<=neq; ++j)
			for (p=1; p<=neq; p++)
				b(p,k) -= m[k][j] * b(p,j);
		for (p=1; p<=neq; p++)
			b(p,k) /= m[k][k];
	}
	return b;
}

dVector dMatrix::solve(dVector b)
{
	double r;
	int neq = no_row;
// Forward reduction
	for (k=2; k<=neq; ++k)
	{
		for (i=k; i<=neq; ++i)
		{
			r = m[i][k-1]/m[k-1][k-1];
			b[i] -= r*b[k-1];
			for (j=k; j<=neq; ++j)
				m[i][j] -= r*m[k-1][j];
		}
	}
// Backsubstitution phase
	b[neq] = b[neq]/m[neq][neq];
	for (k=neq-1; k>=1; --k)
	{
		for (j=k+1; j<=neq; ++j)
			b[k] -= m[k][j] * b[j];
		b[k] /= m[k][k];
	}
	return b;
}


dVector dMatrix::operator^(dVector b)
{
	int N = no_row;
	int S = no_col;
	if (N != b.size()) cout << "\n ***  Matrices not compatible  *** \n";
//From R D Cook et al (p595)
	long double dum;
	int i, j, k, l, n;
	int lim;
// Treat the case of one or more independent equations
	if (S <= 1)
	{
		for (i=1; i<=N; ++i)
			b[i] = b[i] / m[N][1];
		return b;
	}
// Forward reduction of stiffness matrix
	for (n=1; n<=N-1; ++n)
	{
		lim = min(S, N + 1 - n);
		for (l=2; l<=lim; ++l)
		{
			dum = m[n][l]/m[n][1];
			i = n + l - 1;
			j = 0;
			for (k=l; k<=lim; ++k)
			{
				j ++;
				m[i][j] -= dum*m[n][k];
			}
			m[n][l] = dum;
		}
	}
// Forward reduction of load vector
	for (n=1; n<=N-1; ++n)
	{
		lim = min(S, N + 1 - n);
		for (l=2; l<=lim; ++l)
		{
			i = n + l - 1;
			b[i] -= m[n][l]*b[n];
		}
		b[n] = b[n]/m[n][1];
	}
	b[N] = b[N]/m[N][1];
// Backsubstitution phase
	for (n=N-1; n>=1; --n)
	{
		lim = min(S, N + 1 - n);
		for (l=2; l<=lim; ++l)
		{
			k = n + l - 1;
			b[n] -= m[n][l] * b[k];
		}
	}
	return b;
}

ostream& operator<<(ostream & os, const dMatrix & b)
{
	os.precision(4);
	os.setf(ios::showpoint);
	for (i=1; i<=b.no_row; i++)
	{
		os << endl;
		for (j=1; j<=b.no_col; j++)
			os << setw(20) << b.m[i][j];
	}
	os << endl << endl;
	return os;
}
istream& operator>>(istream& is, const dMatrix& b)
{
	for (i=1; i<=b.no_row; i++)
		for (j=1; j<=b.no_col; j++)
			is >> b.m[i][j];
	return is;
}
class iMatrix // double Matrix
{
private:
	int no_row;
	int no_col;
	int **m;

public:
	// Constructors & Destructor
	iMatrix(void);
	iMatrix(const iMatrix& b);
	iMatrix(int nr, int nc);
	~iMatrix(void);

	// Operators
	iMatrix operator=(const iMatrix& b);
	iMatrix operator+=(const iMatrix& b);
	iMatrix operator+(const iMatrix& b);
//	iVector operator*(iVector b); // multiplication by a vector
	int& iMatrix::operator()(int i, int j); // Index operator
	iMatrix operator*(int b);
	iMatrix operator*(iMatrix b); // multiplication by a matrix

	// Friend functions
	friend iMatrix operator*(int b, iMatrix& a);
	friend ostream& operator<<(ostream & os, const iMatrix & b);
	friend istream& operator>>(istream & is, const iMatrix & b);
};
iMatrix::iMatrix(void) // Default constructor
{
	no_row = 0;
	no_col = 0;
	m = 0;
}
iMatrix::iMatrix(const iMatrix & b) // Copy constructor
{
	m = new int *[no_row = b.no_row];
	m --;
	for (i=1; i<=no_row; i++)
	{
		m[i] = new int [no_col = b.no_col];
		m[i] --;
	}
	for (i=1; i<=no_row; i++)
		for (int j=1; j<=no_col; j++)
			m[i][j] = b.m[i][j];
}
iMatrix::iMatrix(int nr, int nc) // constructor
{
	m = new int *[no_row = nr];
	m --;
	for (i=1; i<=no_row; i++)
	{
		m[i] = new int [no_col = nc];
		m[i] --;
	}
	for (i=1; i<=no_row; i++)
		for (int j=1; j<=no_col; j++)
			m[i][j] = 0;
}
iMatrix::~iMatrix() // Destructor
{
/*	if (++m)
	{
		for (i=1; i<=no_row; i++)
			m[i] ++;
		delete (m);
	}*/
}
iMatrix iMatrix::operator=(const iMatrix & b) // Assignment operator
{
	if (m == b.m)
		return *this;
	m ++;
	for (i=1; i<=no_row; i++)
		m[i] --;
//	delete [][]m;
	m = new int *[no_row = b.no_row];
	m --;
	for (i=1; i<=no_row; i++)
	{
		m[i] = new int [no_col = b.no_col];
		m[i] --;
	}
	for (i=1; i<=no_row; i++)
		for (int j=1; j<=no_col; j++)
			m[i][j] = b.m[i][j];
	return *this;
}
int& iMatrix::operator()(int i, int j) // Index operator
{
	if (m && i>=1 && i<=no_row && j>=1 && j<=no_col)
		return m[i][j];
	else 
	{
		cout << " error in index";
		return m[i][j];
	}
}
iMatrix iMatrix::operator+=(const iMatrix & b) // Addition with assign operator
{
	for (i=1; i<=no_row; i++)
		for (j=1; j<=no_col; j++)
			m[i][j] += b.m[i][j];
	return *this;
}
iMatrix iMatrix::operator+(const iMatrix & b) // Addition operator
{
	iMatrix a(*this);
	a += b;
	return a;
}
iMatrix iMatrix::operator*(int b)
{
	for (i=1; i<=no_row; i++)
		for (j=1; j<=no_col; j++)
			m[i][j] *= b;
	return *this;
}
iMatrix operator*(int b, iMatrix& a) // multiplication by a constant
{
	return a*b;
}
iMatrix iMatrix::operator*(iMatrix b) // multiplication by a constant
{
	iMatrix c(no_row,b.no_col);
	for (i=1; i<=c.no_row; i++)
		for (j=1; j<=c.no_col; j++)
		{
			c.m[i][j] = 0;
			for (k=1; k<=no_col; k++)
				c.m[i][j] += m[i][k]*b.m[k][j];
		}
	return c;
}
/*iVector iMatrix::operator*(iVector b)
{
	iVector c(b);
	for (i=1; i<=no_row; i++)
	{
		c[i] = 0;
		for (j=1; j<=no_col; j++)
			c[i] += m[i][j]*b[j];
	}
	return c;
}*/

ostream& operator<<(ostream & os, const iMatrix & b)
{
	os.precision(4);
	os.setf(ios::showpoint);
	os << endl;
	for (i=1; i<=b.no_row; i++)
	{
		os << endl;
		for (j=1; j<=b.no_col; j++)
			os << setw(20) << b.m[i][j];
	}
	os << endl << endl;
	return os;
}
istream& operator>>(istream& is, const iMatrix& b)
{
	for (i=1; i<=b.no_row; i++)
		for (j=1; j<=b.no_col; j++)
			is >> b.m[i][j];
	return is;
}

int min(int a, int b)
{
	if (a < b) 
		return (a);
	else
		return (b);
}


#endif