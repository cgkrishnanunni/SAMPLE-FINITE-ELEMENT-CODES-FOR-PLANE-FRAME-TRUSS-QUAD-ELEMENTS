// Mesh generator for a cantilever beam
// Uses three-noded triangular elements

#include <fstream>
using namespace std;

void main(void)
{
	int  nDivSpan, nDivDepth;
	int nNodes, nElems, probType;
	int i, j, k, elem, node;
	double x[10000], y[10000];
	int elemConn[10][10000];
	double span, depth, thick, E, po, load;

	ifstream fin("genTriangle.inp");
	ofstream fout("Tri.inp");
	
	fout.precision(20);
	
	fin >> probType;
	fin >> span >> depth >> thick >> nDivSpan >> nDivDepth;
	fin >> E >> po;
	fin >> load;

	nElems = 2*nDivDepth*nDivSpan;
	nNodes = (nDivDepth + 1)*(nDivSpan + 1);

// Generating nodal coordinate data
	node = 0;
	double xx;
	for (j=1; j<=nDivSpan + 1; ++j)
	{
		xx = span/nDivSpan*(j - 1);
		for (i=1; i<=nDivDepth + 1; ++i)
		{
			node ++;
			x[node] = xx;
			y[node] = depth/nDivDepth*(i - 1);
		}
	}

// Generate element connectivity

	// First element
	elemConn[1][1] = 1;
	elemConn[2][1] = nDivDepth + 2;
	elemConn[3][1] = 2;

	// Second element
	elemConn[1][2] = 2;
	elemConn[2][2] = nDivDepth + 2;
	elemConn[3][2] = nDivDepth + 3;
	
	// First layer along left strip
	elem = 2;
	for (i=2; i<=nDivDepth; ++i)
	{
		for (k=1; k<=2; k++) // elements like first/second elements 
		{
			elem ++;
			for (j=1; j<=3; ++j)
				elemConn[j][elem] = elemConn[j][elem-2] + 1;
		}
	}

	// Subsequent layers
	elem = 2*nDivDepth;
	for (i=2; i<=nDivSpan; ++i)
	{
		for (j=1; j<=nDivDepth; ++j)
		{
			for (int l=1; l<=2; l++) // elements like first/second elements 
			{
				elem ++;
				for (k=1; k<=3; ++k)
					elemConn[k][elem] = elemConn[k][elem-2*nDivDepth] + (nDivDepth + 1);
			}
		}
	}

	fout.precision(12);
// Printing results
	fout << "tri.m\n";
	fout << probType << "  " << nNodes << "  " << nElems << "\n";

// Printing nodal coordinates
	for (i=1; i<=nNodes; ++i)
		fout << "\n " << i << "   " << x[i] << "    " << y[i];

//	Printing element Connectivity
	for (i=1; i<=nElems; ++i)
	{
		fout << "\n"<< i << "  ";
		for (j=1; j<=3; ++j)
			fout << elemConn[j][i] << "  ";
		fout << "  " << thick << "  " << E << "  " << po;
	}

// Displacement boundary condition data
	fout << "\n" << nDivDepth + 1;
	fout << "\n1  1  0  1  0";
	for (j=2; j<=nDivDepth + 1; ++j)
		fout << "\n" << j << "  1  0  0  0";

// Tractions at the tip
	node = (nDivDepth + 1)*nDivSpan;
	int nLoad = nDivDepth + 1;
	fout << "\n" << nLoad;
	double lenth = depth/nDivDepth;
	for (i=1; i<=nLoad; ++i)
	{
		double eq_load;
		node ++;
		if (i == 1 ||  i == nLoad) // first or last node on the side
			eq_load = load*lenth/2.0;
		else
			eq_load = load*lenth;
		fout << "\n" << node << "  0  " << eq_load;
	}
}
