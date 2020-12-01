// finite element analysis - using truss elements
// Static Analysis

#include <fstream>
#include <iomanip>
#include<stdio.h>
#include<stdlib.h>

#include "Mat.h"
#include "Vec.h"

void trussElem();
void plot(int nNodes, int nElems, dVector xp, dVector yp, iMatrix nod);
void erase();

dMatrix eStiff(4,4);
double elemArea, E, L, cost, sint;
char color[3], file[10];

using namespace std;

int main(void)
{
	ifstream fin ("truss.inp");
//	ifstream fin ("gentruss.out");
	ofstream fout("truss.out");
	
    fin >> file;
	double magDisp;
	fin >> magDisp;

	int i, j, k;
	int nElems, nNodes;
	
	fin >> nNodes >> nElems;

	dVector x(nNodes), y(nNodes);
	
	// input coordinate data:
	for (i=1; i<=nNodes; i++)
	{
		fin >> k;
		fin >> x[k] >> y[k];
	}
	fout << "\n*********************************\n";
	fout << "\n FINITE ELEMENT ANALYSIS PROGRAM\n";
	fout << "\n*********************************\n\n";
	fout << "\nProblem Type: Truss Analysis\n";
	fout << "\nNumber of nodes    =  " << nNodes;
	fout << "\nNumber of elements =  " << nElems;
	fout << "\n\nNodal Coordinates";
	fout << "\n~~~~~~~~~~~~~~~~~\n";
	fout << "\nNode    X-coord.     Y-coord.\n";
	fout <<"=============================";

	// print nodal coordinates:
	for (i=1; i<=nNodes; i++)
		fout << "\n" << setw(2) << i << setw(12) << x[i] << setw(12) << y[i];

	iMatrix nod(2,nElems);
	dVector area(nElems);

	fin >> E; // Young's modulus
	
	// input nodal connectivity and element c/s area:
	for (i=1; i<=nElems; i++)
	{
		fin >> k;
		fin >> nod(1,k) >> nod(2,k) >> area[k];
	}
	fout << "\n\nElement Connectivity";
	fout << "\n~~~~~~~~~~~~~~~~~~~~\n";
	fout <<"\n  Elem      Nod1     Nod2      Area \n";
	fout <<"===================================\n";
	// print nodal connectivity:
	for (i=1; i<=nElems; i++)
		fout << "\n" << setw(4) << i << setw(10) << nod(1,i) << setw(10) << nod(2,i) << setw(10) << area[i];
	
	int nDispBC;//, nLoadedNodes;
	// input boundary constraint list:
	fin >> nDispBC;
	iMatrix id(2,nNodes);
	for (i=1; i<=nDispBC; i++)
	{
		fin >> k;
		fin >> id(1,k) >> id(2,k);
	}

	int neq = 0;
	// update id array:
	for (j=1; j<=nNodes; ++j)
		for (i=1; i<=2; ++i)
		{
			if (id(i, j) == 0)
			{
				neq++;
				id(i, j) = neq;
				continue;
			}
			else
				id(i, j) = 0;
		}
	// print id array:
	fout << "\n\nThe Destination Array:";
	fout << "\n~~~~~~~~~~~~~~~~~~~~~\n";
	fout << "\n Node      X-dof     Y-dof\n";
	fout <<"==========================";
	for (i=1; i<=nNodes; i++)
		fout << "\n" << setw(4) << i << setw(10) << id(1,i) << setw(10) << id(2,i);
	
	// input load data:
	int nLoadedNodes;
	fin >> nLoadedNodes;
	dVector gLoad(neq), gDisp(neq);
	for (i=1; i<=nLoadedNodes; i++)
	{
		double loadX, loadY;
		fin >> k >> loadX >> loadY;
		int dof1, dof2;
		dof1 = id(1,k);
		dof2 = id(2,k);
		gLoad[dof1] = loadX;
		gLoad[dof2] = loadY;
	}

	// Plot the undeformed configuration:
	color[0] = 'r'; color[1] = ':';
	erase();
	plot(nNodes, nElems, x, y, nod);

// Bandwidth calculation
	int small, large, diff;
	int n, kk[5];
	int semiband = 0;
	for (i=1; i<=nElems; ++i)  // scan over each element
	{
		small = neq;
		large = 1;
		for (j=1; j<=2; ++j)  // scan over each node
		{
			int node = nod(j, i);
			for (k=1; k<=2; ++k)  // scan over each dof
			{
				int dof = id(k, node);
				if (dof == 0) continue;
				if (dof < small) small = dof;
				if (dof > large) large = dof;
				
			}
		}
		diff = large - small;// diff gives the semi-bandwidth
		if (diff > semiband) // pick up the largest of diff
			semiband = diff;
	}
	semiband = semiband + 1;
	fout << "\n\nSemi-band width = " << semiband;
	dMatrix gStiff(neq, semiband);

	// assemble element stiffness 
	for (n=1; n<=nElems; n++)
	{
		elemArea = area[n];
		int nod1, nod2;
		nod1 = nod(1,n);
		nod2 = nod(2,n);
		double xl, yl;
		xl = x[nod2] - x[nod1];
		yl = y[nod2] - y[nod1];
		L = sqrt(xl*xl + yl*yl);
		cost = xl/L;
		sint = yl/L;
		
		trussElem();
		int dof = 0;
		for (i=1; i<=2; ++i)
		{
			int node = nod(i, n);
			for (j=1; j<=2; ++j)
				{
					dof ++;
					kk[dof] = id(j,node);
				}
		}
		for (i=1; i<=4; ++i)
		{
			if (kk[i] <= 0)
				continue;
			k = kk[i];
			for (j=1; j<=4; ++j)
			{
				if (kk[j] < k)
					continue;
				int l = kk[j] - k + 1;
				gStiff(k, l) += eStiff(i,j);
			}
		}
	}

	fout << "\n[K] = " << gStiff;
	fout << "\n{R} = " << gLoad;
	
	gDisp = gStiff ^ gLoad;
	
	fout <<"\n{gDisp}" << gDisp;
	
	fout << "\n\nGlobal Displacement Vector";
	fout <<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	fout << "\n\n Node                Disp-X                 Disp-Y";
	fout << "\n=======================================================";
	
	for (j=1; j<=nNodes; ++j)
	{
		fout << "\n" << setw(4) << j << setw(15);
		for (i=1; i<=2; ++i)
		{
			int dof = id(i, j);
			if (dof != 0)
				fout << setw(25) << gDisp[dof];
			else
				fout << setw(25) << "0.0000";
		}
	}	
// Plotting the deformed configuration:
	dVector xp(nNodes), yp(nNodes);
	for (j=1; j<=nNodes; j++)
	{
		int dof1, dof2;
		dof1 = id(1,j);
		dof2 = id(2,j);
		if (dof1 != 0)
			xp[j] = x[j] + magDisp*gDisp[dof1];
		else 
			xp[j] = x[j];
			if (dof2 != 0)
			yp[j] = y[j] + magDisp*gDisp[dof2];
		else 
			yp[j] = y[j];
	}
	color[0] = 'k'; color[1] = ' ';
	plot(nNodes, nElems, xp, yp, nod);

	// Internal forces in each element
	dVector barForce(nElems);
	for (n=1; n<=nElems; n++)
	{
		elemArea = area[n];
		int nod1, nod2;
		nod1 = nod(1,n);
		nod2 = nod(2,n);
		double xl, yl;
		xl = x[nod2] - x[nod1];
		yl = y[nod2] - y[nod1];
		L = sqrt(xl*xl + yl*yl);
		cost = xl/L;
		sint = yl/L;
		
		trussElem();
		int dof = 0;
		for (i=1; i<=2; ++i)
		{
			int node = nod(i, n);
			for (j=1; j<=2; ++j)
				{
					dof ++;
					kk[dof] = id(j,node);
				}
		}

		dVector eDisp(4), eLoad(4);
		for (i=1; i<=4; ++i)
		{
			if (kk[i] <= 0)
				continue;
			
			k = kk[i];
			eDisp[i] += gDisp[k];
		}
		fout << "\nElement No: " << n;
		fout << "\n{u} = " << eDisp;
		fout << "\n[k] = " << eStiff;
		eLoad = eStiff * eDisp;
		fout << "\n{r} = " << eLoad;
		barForce[n] = -(eLoad[1]*cost + eLoad[2]*sint);
	}
	fout << "\nMember forces are:" << barForce;
    return 0;
}

void trussElem()
{
	eStiff = 0;
	eStiff(1,1) = cost*cost;
	eStiff(1,2) = cost*sint;
	eStiff(1,3) = -cost*cost;
	eStiff(1,4) = -cost*sint;

	eStiff(2,2) = sint*sint;
	eStiff(2,3) = -cost*sint;
	eStiff(2,4) = -sint*sint;

	eStiff(3,3) = cost*cost;
	eStiff(3,4) = cost*sint;
	 
	eStiff(4,4) = sint*sint;

	for (int i=2; i<=4; i++)
		for (int j=1; j<i; j++)
			eStiff(i,j) = eStiff(j,i);

	double aeL = elemArea*E/L;
	eStiff = aeL* eStiff;
}

void plot(int nNodes, int nElems, dVector xp, dVector yp, iMatrix nod)
{
	ofstream fout(file, ios::app);
	for (int i=1; i<=nElems; ++i)
	{
		fout << "\nx = [" << xp[nod(1,i)] << " " << xp[nod(2,i)] << "];";
		fout << "\ny = [" << yp[nod(1,i)] << " " << yp[nod(2,i)] << "];";
		fout << "\nplot(x,y,'" << color << "','Linewidth',2)";
		fout << "\nhold on";
	}
	for (int i=1; i<=nNodes; ++i)
	{
		fout << "\nx = " << xp[i] << ";\ny = " << yp[i] << ";";
		fout <<"\nplot (x,y,'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',4)";
		fout << "\nhold on";
	}
	fout << "\naxis off";
}
void erase()
{
		ofstream fout(file);
}
