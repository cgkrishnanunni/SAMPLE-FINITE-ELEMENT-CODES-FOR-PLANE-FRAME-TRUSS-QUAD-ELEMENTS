// Mesh generator for a latticed beam made of truss elements

#include <fstream.h>
#include <iomanip.h>
#include "mat.h"
#include "vec.h"

void main(void)
{
	int  nDivSpan, nDivDepth;
	int nNodes, nElems;
	int i, j, k, elem, node;
	double span, depth, E, load, areaH, areaV, areaD; // area of horz, vert, and diag. elements

	ifstream fin("genTruss.inp");
	ofstream fout("genTruss.out");
	
	fout.precision(20);

	double magDisp;

	fin >> magDisp;
	fin >> span >> depth >> nDivSpan >> nDivDepth;
	fin >> areaH >> areaV >> areaD;
	fin >> E;
	fin >> load;

	nElems = 4*nDivDepth*nDivSpan + nDivDepth + nDivSpan;
	nNodes = (nDivDepth + 1)*(nDivSpan + 1);

	dVector x(nNodes), y(nNodes), area(nElems);
	iMatrix elemConn(2, nElems);

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

	// First unit at bottom left corner (one vert, one horz and two diags)
	// vert
	elemConn(1,1) = 1;
	elemConn(2,1) = 2;
	area[1] = areaV;
	
	// horz.
	elem = nDivDepth + 1;
	elemConn(1,elem) = 1;
	elemConn(2,elem) = nDivDepth + 2;
	area[elem] = areaH;

	// diag1
	elem ++;
	elemConn(1,elem) = 1;
	elemConn(2,elem) = nDivDepth + 3;
	area[elem] = areaD;

	// diag2
	elem ++;
	elemConn(1,elem) = 2;
	elemConn(2,elem) = nDivDepth + 2;
	area[elem] = areaD;

	// First layer along left strip
	for (i=2; i<=nDivDepth; ++i)
	{
		// vert
		elem = i;
		elemConn(1,elem) = elemConn(1,elem - 1) + 1;
		elemConn(2,elem) = elemConn(2,elem - 1) + 1;
		area[elem] = areaV;
		
		// horz.
		elem = nDivDepth + (i - 1)*4;
		elemConn(1,elem) = elemConn(1,elem - 3) + 1;
		elemConn(2,elem) = elemConn(2,elem - 3) + 1;
		area[elem] = areaH;

		// diag1
		elem ++;
		elemConn(1,elem) = elemConn(1,elem - 3) + 1;
		elemConn(2,elem) = elemConn(2,elem - 3) + 1;
		area[elem] = areaD;

		// diag2
		elem ++;
		elemConn(1,elem) = elemConn(1,elem - 3) + 1;
		elemConn(2,elem) = elemConn(2,elem - 3) + 1;
		area[elem] = areaD;
	}

	// top one horz.
	elem ++;
	elemConn(1,elem) = elemConn(1,elem - 3) + 1;
	elemConn(2,elem) = elemConn(2,elem - 3) + 1;
	area[elem] = areaH;

	// Subsequent layers
	for (i=2; i<=nDivSpan; ++i)
	{
		for (j=1; j<=nDivDepth; ++j)
		{
			for (k=1; k<=4; ++k)
			{
				elem ++;
				int elemIncr = 4*nDivDepth + 1;
				int nodeIncr = nDivDepth + 1;
				elemConn(1,elem) = elemConn(1,elem - elemIncr) + nodeIncr;
				elemConn(2,elem) = elemConn(2,elem - elemIncr) + nodeIncr;
				area[elem] = area[elem - elemIncr];
			}
		}
		// top one horz.
		elem ++;
		elemConn(1,elem) = elemConn(1,elem - 3) + 1;
		elemConn(2,elem) = elemConn(2,elem - 3) + 1;
		area[elem] = areaH;
	}

	// right edge vert. members
	for (i=1; i<=nDivDepth; ++i)
	{
		elem ++;
		elemConn(1,elem) = nNodes - nDivDepth + i - 1;
		elemConn(2,elem) = nNodes - nDivDepth + i;
		area[elem] = areaV;
	}

	fout.precision(12);
// Printing results
	fout << magDisp;
	fout << "\n" << nNodes << "   " << nElems;

// Printing nodal coordinates
	for (i=1; i<=nNodes; ++i)
		fout << "\n " << i << "   " << x[i] << "    " << y[i];

	fout << "\n\n" << E << endl;

//	Printing element connectivity and area
	for (i=1; i<=nElems; ++i)
	{
		fout << "\n "<< i << setw(10);
		for (j=1; j<=2; ++j)
			fout << elemConn(j,i) << setw(12);
		fout << area[i];
	}

// Displacement boundary condition data
	fout << "\n" << 2;
	fout << "\n1  1  1";
	fout << "\n" << nNodes - nDivDepth << "  0  1";

// Load
	fout << "\n\n1";
	fout << "\n" << nNodes/2 << "  0  " << load; 
}
