// Program for Finite Element Analysis of Elastostatic Problems
// Using 3-noded Constant Strain Triangle Elements
// Non-zero displacement dof can be prescribed; uses dynamic memory allocation

#include <fstream>
#include <iostream>
#include <iomanip>

#include <math.h>
#include "Mat.h"
#include "Vec.h"

using namespace std;

void triang3(void);
void triang3stress(void);
void plotElem(void);
void erase(void);

dMatrix eStiff(6,6), B(3,6), Jac(2,2), D(3,3);
dVector eLoad(6), xl(3), yl(3), stress(3), strain(3), u(6);
int neq, semiband;
double detJac;
double thickness;
char color[3], file[10];

int main(void)
{
	ifstream fin ("femTriang3.inp");
//	ifstream fin ("Tri.inp");
	ofstream fout("femTriang3.out");
	fout.setf(ios::showpoint);
	fout.setf(ios::floatfield, ios::fixed);
	fout.precision(6);
	int probType, nElems, nNodes, nRestNodes, nLoadedNodes;
	int  i, j, k;

	// Input data::	
	fin >> file;
	fin >> probType;
	fin >> nNodes >> nElems;
	 
	// input nodal coordinates
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
	if (probType == 0)
		fout << "Problem Type: Plane Stress\n";
	else
		fout << "Problem Type: Plane Strain\n";
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

	iMatrix elemConn(6, nElems);
	fout << "\n\nElement Connectivity";
	fout << "\n~~~~~~~~~~~~~~~~~~~~\n";
	fout <<"\n  Elem Nod1  Nod2   Nod3   Thickness         E            nu\n";
	fout <<"==================================================================\n";
	
	dVector thick(nElems), nu(nElems), E(nElems);
	for (i=1; i<=nElems; ++i)
	{
		fin >> k;
		for (j=1; j<=3; ++j)
		fin >> elemConn(j, k);
		fin >> thick[k] >> E[k] >> nu[k];
	}

	for (i=1; i<=nElems; ++i)
	{
		fout << "\n" << setw(4) << i << setw(5);
		for (j=1; j<=3; ++j)
			fout << elemConn(j, i) << setw(7);
		fout << setw(10) << thick[i] << setw(15) << E[i] << setw(15) << nu[i];
	}

	dVector U(2*nNodes);
	iMatrix destn(2, nNodes);

	int dof;
	fin >> nRestNodes;
	fout << "\n\nNumber of nodes at which displacement is prescribed  =  " << nRestNodes;
	for (i=1; i<=nRestNodes; ++i)
	{
		fin >> k;
		for (j=1; j<=2;++j)
		{
			dof = 2*(k-1) + j;
			fin >> destn(j,k) >> U[dof];
		}
	}
	neq = 0;
	for (j=1; j<=nNodes; ++j)
		for (i=1; i<=2; ++i)
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
	fout << "\n\nThe Destination Array:";
	fout << "\n~~~~~~~~~~~~~~~~~~~~~\n";
	fout << "\nNode    X-dof     Y-dof\n";
	fout <<"========================";
	for (i=1; i<=nNodes; ++i)
	{
		fout << "\n" << setw(4) << i << setw(9);
		for (j=1; j<=2; ++j)
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
		for (j=1; j<=2; ++j)
		{
			dof = destn(j, k);
			fin >> load;
			if (dof != 0)
				gLoad[dof] += load;
		}
	}
	int node;
	
	// Plot the undeformed mesh:
	color[0] = 'b';color[1] = ':';
	erase ();
	for (i=1; i<=nElems; i++)
	{
		for (j=1; j<=3; ++j)
		{
			int node = elemConn(j, i);
			xl[j] = X(1, node);
			yl[j] = X(2, node);
		}
		plotElem();
	}
	

// Bandwidth calculation
	int small, large, diff;
	iVector kk(6);
	int n, l;
	semiband = 0;
	for (i=1; i<=nElems; ++i)  // scan over each element
	{
		small = neq;
		large = 1;
		for (j=1; j<=3; ++j)  // scan over each node
		{
			node = elemConn(j, i);
			for (k=1; k<=2; ++k)  // scan over each dof
			{
				dof = destn(k, node);
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
	fout << "\nSemi-band width = " << semiband;
	int dof1;
	dMatrix gStiff(neq, semiband);

	fout << "\n\nGlobal Load Vector:";
	fout << "\n~~~~~~~~~~~~~~~~~~~~~\n";
	fout << "\n Node        Load-X         Load-Y";
	fout <<"\n==================================";
	int dof2;
	for (i=1; i<=nNodes; ++i)
	{
		dof1 = destn(1, i);
		dof2 = destn(2, i);
		double zero = 0.0;
		fout << "\n" << setw(4) << i;
		if (dof1 != 0)
			fout << setw(15) << gLoad[dof1];
		else
			fout << setw(15) << zero;
		if (dof2 != 0)
			fout << setw(15) << gLoad[dof2];
		else
			fout << setw(15) << zero;
	}

// Assembly of element matrices	
	dVector presDisp(6), eLoadDisp(6);
	for (n=1; n<=nElems; ++n)
	{
		for (i=1; i<=3; i++)
		{
			node = elemConn(i,n);
			xl[i] = X(1,node);
			yl[i] = X(2,node);
		}
		if (probType == 0) 
		{
			// plane stress
			D(1,1) = D(2,2) = E[n]/(1 - nu[n]*nu[n]);
			D(1,2) = D(2,1) = D(1,1) * nu[n];
			D(3,3) = D(1,1) * (1 - nu[n])/2.0;
		}
		else
		{
			// plane strain
			double c;
			c = E[n]/(1 + nu[n])/(1 - 2*nu[n]);
			D(1,1) = D(2,2) = c*(1 - nu[n]);
			D(1,2) = D(2,1)  = c * nu[n];
			D(3,3)= c * (1 - 2*nu[n])/2.0;
		}
		D(1,3) = D(2,3) = D(3,1) = D(3,2)= 0;

		thickness = thick[n];
		triang3();
		dof = 0;
		for (i=1; i<=3; ++i)
		{
			node = elemConn(i, n);
			for (j=1; j<=2; ++j)
			{
				dof ++;
				kk[dof] = destn(j, node);
			}
		}

		dof1 = 0;
		for (i=1; i<=3; ++i)
		{
			node = elemConn(i, n);
			for (j=1; j<=2; ++j)
			{
				dof1 ++;
				dof = 2*(node - 1) + j;
				presDisp[dof1] = U[dof];
			}
		}

		eLoadDisp = eStiff * presDisp;
//		fout << "\neLoadDisp = " << eLoadDisp;
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
				l = kk[j] - k + 1;
				gStiff(k, l) += eStiff(i,j);
			}
		}
	}

	fout << "\n[ K ] = " << gStiff;
	fout << "\n[ R ] = " << gLoad;
	
	gDisp = gStiff^gLoad;
	
	fout << "\n\nGlobal Displacement Vector";
	fout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	fout << "\n\n Node        Disp-X         Disp-Y";
	fout << "\n====================================";
	for (j=1; j<=nNodes; ++j)
	{
		for (i=1; i<=2; ++i)
		{
			dof = 2*(j-1) + i;
			if (destn(i, j) != 0)
				U[dof] = gDisp[destn(i, j)];
		}
	}

	for (j=1; j<=nNodes; ++j)
	{
		fout << "\n" << setw(4) << j << setw(15); 
		for (i=1; i<=2; ++i)
		{
			dof = 2*(j-1) + i;
			fout << U[dof] << setw(15);
		}
	}

	// Update nodal coordinates for plotting the deformed mesh
	double big = 0.;
	double magna_disp;
	for (i=1; i<=neq; i++)
	{
		if (big <= fabs(gDisp[i]))
			big = fabs(gDisp[i]);
	}
	if (big <= 1.e-25)
		magna_disp = 0.;
	else
		magna_disp = 0.01*X(1, nNodes)/big;
	
	color[0] = 'k';color[1] = ' ';
	for (i=1; i<=nElems; i++)
	{
		for (j=1; j<=3; ++j)
		{
			node = elemConn(j, i);
			xl[j] = X(1, node) + U[2*node-1] * magna_disp;
			yl[j] = X(2, node) + U[2*node  ] * magna_disp;
		}
		plotElem();
	}
	
// Compute stresses
	fout << "\n\nConstant Stresses and Strains in Element";
	fout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	fout << "\n\n  Elem       Sig_x         Sig_y          Sig_xy ";
	fout << "\n================================================";
	for (n=1; n<=nElems; ++n)
	{
		thickness = thick[n];
		if (probType == 0) 
		{
			// plane stress
			D(1,1) = D(2,2) = E[n]/(1 - nu[n]*nu[n]);
			D(1,2) = D(2,1) = D(1,1) * nu[n];
			D(3,3) = D(1,1) * (1 - nu[n])/2.0;
		}
		else
		{
			// plane strain
			double c;
			c = E[n]/(1 + nu[n])/(1 - 2*nu[n]);
			D(1,1) = D(2,2) = c*(1 - nu[n]);
			D(1,2) = D(2,1)  = c * nu[n];
			D(3,3)= c * (1 - 2*nu[n])/2.0;
		}
		D(1,3) = D(2,3) = D(3,1) = D(3,2)= 0;
		for (int ii=1; ii<=3; ++ii)
		{
			node = elemConn(ii, n);
			xl[ii] = X(1, node);
			yl[ii] = X(2, node);
		}

		// Get element displacement vector
		int dof, dof1, node;
		dof1 = 0;
		for  (int ii=1; ii<=3; ++ii)
		{
			node = elemConn(ii, n);
			for (j=1; j<=2; ++j)
			{
				dof1 ++;
				dof = 2*(node - 1) + j;
				u[dof1] = U[dof];
			}
		}
		triang3stress();
		fout << "\n" << setw(4) << n ;
		for (i=1; i<=3; i++)
		{
			fout << setw(15) << stress[i];
		    
		}
    }
    return 0;
}

void erase(void)
{
	ofstream fout(file);
}
void triang3(void)
{
// generate element stiffness coefficients of 3-noded triangles
	double b1, b2, b3, c1, c2, c3, area;

	ofstream fout1("tri.out");
	B = 0;
	b1 = yl[2] - yl[3];
	b2 = yl[3] - yl[1];
	b3 = yl[1] - yl[2];
	
	c1 = xl[3] - xl[2];
	c2 = xl[1] - xl[3];
	c3 = xl[2] - xl[1];
	
	area = 0.5*(c3*b2 - c2*b3);
	
	B(1,1) = b1;
	B(1,3) = b2;
	B(1,5) = b3;

	B(2,2) = c1;
	B(2,4) = c2;
	B(2,6) = c3;

	B(3,1) = c1;
	B(3,3) = c2;
	B(3,5) = c3;
	B(3,2) = b1;
	B(3,4) = b2;
	B(3,6) = b3;

	eStiff = B.transpose() * D * B * thickness / area / 4.0;
	fout1 << "\n[ B ] = " << B;
}

void triang3stress(void)
{
// stress computation
	double b1, b2, b3, c1, c2, c3, area;

	B = 0;
	b1 = yl[2] - yl[3];
	b2 = yl[3] - yl[1];
	b3 = yl[1] - yl[2];
	
	c1 = xl[3] - xl[2];
	c2 = xl[1] - xl[3];
	c3 = xl[2] - xl[1];
	
	area = 0.5*(c3*b2 - c2*b3);
	
	B(1,1) = b1;
	B(1,3) = b2;
	B(1,5) = b3;

	B(2,2) = c1;
	B(2,4) = c2;
	B(2,6) = c3;

	B(3,1) = c1;
	B(3,3) = c2;
	B(3,5) = c3;
	B(3,2) = b1;
	B(3,4) = b2;
	B(3,6) = b3;

	stress = D * (B / area / 2.0) * u;
}
void plotElem()
{
	ofstream fout(file, ios::app);
	fout << "\nx = [" << xl[1] << " " << xl[2] << "];";
	fout << "\ny = [" << yl[1] << " " << yl[2] << "];";
	fout << "\nplot(x,y,'" << color <<"','Linewidth',2)";
	fout << "\nhold on";
	fout << "\nx = [" << xl[2] << " " << xl[3] << "];";
	fout << "\ny = [" << yl[2] << " " << yl[3] << "];";
	fout << "\nplot(x,y,'" << color <<"','Linewidth',2)";
	fout << "\nhold on";
	fout << "\nx = [" << xl[3] << " " << xl[1] << "];";
	fout << "\ny = [" << yl[3] << " " << yl[1] << "];";
	fout << "\nplot(x,y,'" << color <<"','Linewidth',2)";
	fout << "\nhold on";
	fout << "\naxis off";
}


