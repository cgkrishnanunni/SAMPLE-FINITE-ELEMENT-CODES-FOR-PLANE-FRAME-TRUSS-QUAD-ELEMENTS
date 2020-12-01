// Program for Finite Element Analysis of Plane Frames
// Non-zero displacement dof can be prescribed; uses dynamic memory allocation

#include <fstream>
#include <iostream>
#include <iomanip>

#include <math.h>

#include "Mat.h"
#include "Vec.h"

void frame();
void plotElem(char file[10]);
void erase(void);

dMatrix eStiff(6,6);
dVector  eLoad(6), xl(2), yl(2);
int neq;
double Area, EE, Izz;

char color[3];

using namespace std;

int main(void)
{
    char file[10] = "frame.m";
	ifstream fin ("planeFrame.inp");
	ofstream fout("planeFrame.out");
//	fout.setf(ios::showpoint);
//	fout.setf(ios::floatfield, ios::fixed);
	//fout.precision(4);
	int nElems, nNodes, nRestNodes, nLoadedNodes;
	int  i, j, k;
	
	double magna_disp;
	fin >> magna_disp;

	// Input data::
	fin >> nNodes >> nElems;
	
	dMatrix X(2,nNodes);
	for (i=1; i<=nNodes; ++i)
	{
		fin >> k;
		for (j=1; j<=2; ++j)
			fin >> X(j, k);
	}
	
	fout << "\n*********************************\n";
	fout << "\n FINITE ELEMENT ANALYSIS PROGRAM\n";
	fout << "\n*********************************\n\n";
	fout << "Problem Tyle: Plane Frame\n";
	fout << "\nNumber of nodes    =  " << nNodes;
	fout << "\nNumber of elements =  " << nElems;
	fout << "\n\nNodal Coordinates";
	fout << "\n~~~~~~~~~~~~~~~~~\n";
	fout << "\nNode    X-coord.     Y-coord.\n";
	fout <<"=============================";
	for (i=1; i<=nNodes; ++i)
	{
		fout << "\n" << setw(4) << i << setw(12);
		for (j=1; j<=2; ++j)
			fout << X(j, i) << setw(12);
	}
	iMatrix elemConn(2, nElems);
	fout << "\n\nElement Connectivity";
	fout << "\n~~~~~~~~~~~~~~~~~~~~\n";
	fout <<"\n  Elem Nod1  Nod2    A       E        I\n";
	fout <<"=================================================\n";
	dVector area(nElems), Iz(nElems), E(nElems);
	for (i=1; i<=nElems; ++i)
	{
		fin >> k;
		for (j=1; j<=2; ++j)
			fin >> elemConn(j, k);
		fin >> area[k] >> E[k] >> Iz[k];
	}
	for (i=1; i<=nElems; ++i)
	{
		fout << "\n" << setw(4) << i << setw(5);
		for (j=1; j<=2; ++j)
			fout << elemConn(j, i) << setw(7);
		fout << setw(15) << area[i] << setw(15) << E[i] << setw(15) << Iz[i];
	}
	
	dVector U(3*nNodes);
	iMatrix destn(3, nNodes);
	
	int dof;
	fin >> nRestNodes; // Note: Fix all the dof of the k-th node (node 3) of the beam element
	fout << "\n\nNumber of nodes at which displacement is prescribed  =  " << nRestNodes;
	for (i=1; i<=nRestNodes; ++i)
	{
		fin >> k;
		for (j=1; j<=3;++j)
		{
			dof = 3*(k-1) + j;
			fin >> destn(j,k) >> U[dof];
		}
	}
	neq = 0;
	for (j=1; j<=nNodes; ++j)
	{
		for (i=1; i<=3; ++i)
		{
			if (destn(i, j) == 0)
			{
				neq++;
				destn(i, j) = neq;
				continue;
			}
			else
				destn(i, j) = 0;
		}
	}
	fout << "\n\nThe Destination Array:";
	fout << "\n~~~~~~~~~~~~~~~~~~~~~\n";
	fout << "\nNode    X-dof     Y-dof      Z-dof\n";
	fout <<"==========================================";
	for (i=1; i<=nNodes; ++i)
	{
		fout << "\n" << setw(4) << i << setw(9);
		for (j=1; j<=3; ++j)
			fout << destn(j,i) << setw(9);
	}
	fout << "\n\nNo. of Degrees of Freedom =  " << neq;
	dVector gLoad(neq);
	dVector gDisp(neq);
	
	double load;
	fin >> nLoadedNodes;
	for (i=1; i<=nLoadedNodes; ++i)
	{
		fin >> k;
		for (j=1; j<=3; ++j)
		{
			dof = destn(j, k);
			fin >> load;
			if (dof != 0)
				gLoad[dof] += load;
		}
	}
	fout << "\ngLoad = " << gLoad;

// Bandwidth calculation
	int small, large, diff;
	int semiband = 0;
	for (i=1; i<=nElems; ++i)  // scan over each element
	{
		small = neq;
		large = 1;
		for (j=1; j<=2; ++j)  // scan over each node
		{
			int node = elemConn(j, i);
			for (k=1; k<=3; ++k)  // scan over each dof
			{
				int dof = destn(k, node);
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
	
	// Plot the undeformed mesh:
	color[0] = 'b';color[1] = ':';
	
	erase ();
	int node;
	for (i=1; i<=nElems; i++)
	{
		for (j=1; j<=2; ++j)
		{
			node = elemConn(j, i);
			xl[j] = X(1, node);
			yl[j] = X(2, node);
		}
		plotElem(file);
	}
	
	// Assembly of element matrices
	iVector kk(6);
	for (int n=1; n<=nElems; ++n)
	{
		dVector presDisp(6), eLoadDisp(6);
		for (i=1; i<=2; i++)
		{
			node  = elemConn(i,n);
			xl[i] = X(1,node);
			yl[i] = X(2,node);
		}
		Area = area[n];
		EE   = E[n];
		Izz  = Iz[n];

		frame ();
	
		dof = 0;
		for (i=1; i<=2; ++i)
		{
			for (j=1; j<=3; ++j)
			{
				dof ++;
				node = elemConn(i, n);
				kk[dof] = destn(j, node);
			}
		}

		int dof1 = 0;
		for (i=1; i<=2; ++i)
		{
			node = elemConn(i, n);
			for (j=1; j<=3; ++j)
			{
				dof1 ++;
				dof = 3*node - 3 + j;
				presDisp[dof1] = U[dof];
			}
		}
		eLoadDisp = eStiff * presDisp;
		fout << "\n" << n << "\n{ r_e } = " << eLoadDisp;
		for (i=1; i<=6; ++i)
		{
			if (kk[i] <= 0)
				continue;
			
			k = kk[i];
			gLoad[k] += eLoad[i] - eLoadDisp[i];
			for (j=1; j<=6; ++j)
			{
				if (kk[j]  < k)
					continue;
			
				int l = kk[j] - k + 1;
				gStiff(k, l) += eStiff(i,j);
			}
		}
	}
	
	fout << "\ngStiff = " << gStiff;
	fout << "\ngLoad = " << gLoad;
	fout <<"\n";
	
	gDisp = gStiff ^ gLoad;

	for (j=1; j<=nNodes; ++j)
	{
		for (i=1; i<=3; ++i)
		{
			dof = 3*(j-1) + i;
			if (destn(i, j) != 0)
				U[dof] = gDisp[destn(i, j)];
		}
	}
	fout << "\n\nGlobal Displacement Vector";
	fout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	fout << "\n\n Node        Disp-X         Disp-Y";
	fout << "\n====================================";
	for (j=1; j<=nNodes; ++j)
	{
		fout << "\n" << setw(4) << j << setw(15); 
		for (i=1; i<=3; ++i)
		{
			dof = 3*(j-1) + i;
			fout << U[dof] << setw(15);
		}
	}

/*	fout << "\ngDisp = " << gDisp;
	
	fout << "\n\nGlobal Displacement Vector";
	fout <<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	fout << "\n\n Node        Disp-X         Disp-Y          rotn";
	fout << "\n=======================================================";
	
	for (j=1; j<=nNodes; ++j)
	{
		fout << "\n" << setw(4) << j << setw(15);
		for (i=1; i<=3; ++i)
		{
			dof = destn(i, j);
			if (dof != 0)
				fout << setw(25) << gDisp[dof];
			else
				fout << setw(25) << "0.0000";
		}
	}
	*/
	color[0] = 'k';color[1] = ' ';
	for (i=1; i<=nElems; i++)
	{
		for (j=1; j<=2; ++j)
		{
			node = elemConn(j, i);
			int dof1 = destn(1,node);
			int dof2 = destn(2,node);

			if (dof1 == 0)
				xl[j] = X(1, node);
			else
				xl[j] = X(1, node) + U[3*node-2] * magna_disp;

			if (dof2 == 0)
				yl[j] = X(2, node);
			else
				yl[j] = X(2, node) + U[3*node-1] * magna_disp;
		}
		plotElem(file);
	}
	// Internal forces in each element
    fout << "\n\n";
	for (int n=1; n<=nElems; ++n)
	{
		for (i=1; i<=2; i++)
		{
			node  = elemConn(i,n);
			xl[i] = X(1,node);
			yl[i] = X(2,node);
		}
		Area = area[n];
		EE   = E[n];
		Izz  = Iz[n];

		frame();
	
		dof = 0;
		for (i=1; i<=2; ++i)
		{
			node = elemConn(i, n);
			for (j=1; j<=3; ++j)
			{
				dof ++;
				kk[dof] = destn(j, node);
			}
		}

		dVector eDisp(6);

		int dof1 = 0;
		for (i=1; i<=2; ++i)
		{
			node = elemConn(i, n);
			for (j=1; j<=3; ++j)
			{
				dof1 ++;
				dof = 3*(node - 1) + j;
				eDisp[dof1] = U[dof];
			}
		}

/*		
		for (i=1; i<=6; ++i)
		{
			if (kk[i] <= 0)
				continue;
			
			k = kk[i];
			eDisp[i] = gDisp[k];
		}*/
		fout << "\nElement No: " << n;
		fout << "\n{u} = " << eDisp;
		fout << "\n[k] = " << eStiff;
		eLoad = eStiff * eDisp;
		fout << "\n{r} = " << eLoad;
	}
    return 0;
}
void erase(void)
{
	ofstream fout("frame.m");
}
void frame()
{
//Refer RD Cook et al. for the following notations;
	double Z, A, B, K, M, F, G, P, H, Q, s, c, xL, yL, length;
	int i, j;
	xL = xl[2] - xl[1];
	yL = yl[2] - yl[1];
	length = sqrt(xL*xL + yL*yL);
	c = xL / length;
	s = yL / length;
	Z = Area * EE / length;
	A = 4.0 * EE * Izz / length;
	B = A / 2.;
	M = 1.5 * A / length;
	K = 2. * M / length;
	F = Z * c*c + K * s*s;
	G = (Z - K) * c*s;
	P = Z * s*s + K * c*c;
	H = -M *s;
	Q = M *c;
	
	eLoad  = 0;
	eStiff = 0;
	
	eStiff(1,1) = F;
	eStiff(1,2) = G;
	eStiff(1,3) = H;
	eStiff(1,4) = -F;
	eStiff(1,5) = -G;
	eStiff(1,6) = H;
	eStiff(2,2) = P;
	eStiff(2,3) = Q;
	eStiff(2,4) = -G;
	eStiff(2,5) = -P;
	eStiff(2,6) = Q;
	eStiff(3,3) = A;
	eStiff(3,4) = -H;
	eStiff(3,5) = -Q;
	eStiff(3,6) = B;
	eStiff(4,4) = F;
	eStiff(4,5) = G;
	eStiff(4,6) = -H;
	eStiff(5,5) = P;
	eStiff(5,6) = -Q;
	eStiff(6,6) = A;

	for (i=1; i<=5; ++i)
		for (j=i; j<=6; ++j)
			eStiff(j,i) = eStiff(i,j);
}

void plotElem(char file[10])
{
	ofstream fout(file, ios::app);
	fout << "\nx = [" << xl[1] << " " << xl[2] << "];";
	fout << "\ny = [" << yl[1] << " " << yl[2] << "];";
	fout << "\nplot(x,y,'" << color <<"','Linewidth',2)";
	fout << "\nhold on";
	fout << "\naxis off";
}

