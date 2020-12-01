// Program for Finite Element Analysis of Elastostatic Problems
// Using Isoparametric 4-noded Quadrilateral Element
// Non-zero displacement dof can be prescribed; uses dynamic memory allocation

#include <iostream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <math.h>
#include "Mat.h"
#include "Vec.h"

void quad4(int, double);
void shapeFunc(double, double);
void quad4Stress(int, double, double, double);
void plotElem();
void erase();

dMatrix eStiff(8,8), B(3,8), Jac(3,3), D(3,3);
dVector eLoad(8), stress(3), strain(3), xl(4), yl(4), u(8);
int nGauss, neq, semiband;
double detJac;
char color[3];
char file[10];

using namespace std;

// Gaussian coordinates and weights
const long double place[5][5] = 
{
	{0., 0., 0.,                 0.,	             0.},
	{0,  0, -0.577350269189626, -0.774596669241483, -0.861136311594053},
	{0,  0,  0.577350269189626,  0,                 -0.339981043584856},
	{0,  0,  0,                  0.774596669241483,  0.339981043584856},
	{0,  0,  0,                  0,                  0.861136311594053}        
};
const long double wgt[5][5] =
{
	{0, 0,  0,  0,                   0},
	{0, 2., 1., 0.555555555555555,   0.347854845137454},
	{0, 0,  1., 0.888888888888888,   0.652145154862546},
	{0, 0,  0,  0.555555555555555,   0.652145154862546},
	{0, 0,  0,  0,                   0.347854845137454}
};
void main(void)
{
	ifstream fin ("femQuad4.inp");
	ofstream fout("femQuad4.out");
	fout.setf(ios::showpoint);
	fout.setf(ios::floatfield, ios::fixed);
	fout.precision(4);
	int probType, nElems, nNodes, nRestNodes, nLoadedNodes;
	int  i, j, k, ii;

	// Input data::
	fin >> file;
	fin >> nGauss >> probType;
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
	iMatrix elemConn(4, nElems);
	fout << "\n\nElement Connectivity";
	fout << "\n~~~~~~~~~~~~~~~~~~~~\n";
	fout <<"\n  Elem Nod1  Nod2   Nod3   Nod4  Thickness       E           Po\n";
	fout <<"====================================================================";
	dVector thick(nElems), nu(nElems), E(nElems);
	for (i=1; i<=nElems; ++i)
	{
		fin >> k;
		for (j=1; j<=4; ++j)
			fin >> elemConn(j, k);
		fin >> thick[k] >> E[k] >> nu[k];
	}
	for (i=1; i<=nElems; ++i)
	{
		fout << "\n" << setw(4) << i << setw(5);
		for (j=1; j<=4; ++j)
			fout << elemConn(j, i) << setw(7);
		fout << setw(10) << thick[i] << setw(15) << E[i] << setw(15) << E[i];
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
		for (j=1; j<=4; ++j)
		{
			int node = elemConn(j, i);
			xl[j] = X(1, node);
			yl[j] = X(2, node);
		}
		plotElem();
	}

// Bandwidth calculation
	int small, large, diff;
	int n, l, kk[9];
	semiband = 0;
	for (i=1; i<=nElems; ++i)  // scan over each element
	{
		small = neq;
		large = 1;
		for (j=1; j<=4; ++j)  // scan over each node
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
	int dof1, tot_dof;
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
	double thickness;
	dVector presDisp(8), eLoadDisp(8);
	for (n=1; n<=nElems; ++n)
	{
		for (i=1; i<=4; i++)
		{
			node = elemConn(i,n);
			xl[i] = X(1,node);
			yl[i] = X(2,node);
		}
		if (probType == 0) 
		{
			// plane stress
			D(1,1) = D(2,2) =  E[n]/(1 - nu[n]*nu[n]);
			D(1,2) = D(2,1) = D(1,1) * nu[n];
			D(3,3) = D(1,1) * (1 - nu[n])/2.0;
		}
		else
		{
			// plane strain
			double c;
			c = E[n]/(1 + nu[n])/(1 - 2*nu[n]);
			D(1,1) = D(2,2) = c*(1 - nu[n]);
			D(1,2) = D(2,1) = c * nu[n];
			D(3,3) = c * (1 - 2*nu[n])/2.0;
		}
		D(1,3) = D(2,3) = D(3,1) = D(3,2) = 0;

		thickness = thick[n];
		quad4(n, thickness);
		dof = 0;
		for (i=1; i<=4; ++i)
		{
			for (j=1; j<=2; ++j)
			{
				dof ++;
				node = elemConn(i, n);
				kk[dof] = destn(j, node);
			}
		}
		tot_dof = 8;
		dof1 = 0;
		for (i=1; i<=4; ++i)
		{
			node = elemConn(i, n);
			for (j=1; j<=2; ++j)
			{
				dof1 ++;
				dof = 2*(node - 1) + j;
				presDisp[dof1] = U[dof];
			}
		}
		eLoadDisp += eStiff * presDisp;
		for (i=1; i<=tot_dof; ++i)
		{
			if (kk[i] <= 0)
				continue;
			k = kk[i];
			gLoad[k] += eLoad[i] - eLoadDisp[i];
			for (j=1; j<=tot_dof; ++j)
			{
				if (kk[j]  < k)
					continue;
				l = kk[j] - k + 1;
				gStiff(k, l) += eStiff(i,j);
			}
		}
	}
	gDisp = gStiff ^ gLoad;
	for (j=1; j<=nNodes; ++j)
	{
		for (i=1; i<=2; ++i)
		{
			dof = 2*(j-1) + i;
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
		magna_disp = 0.05*X(1, nNodes)/big;
	
	color[0] = 'k';color[1] = ' ';
	for (i=1; i<=nElems; i++)
	{
		for (j=1; j<=4; ++j)
		{
			node = elemConn(j, i);
			xl[j] = X(1, node) + U[2*node-1] * magna_disp;
			yl[j] = X(2, node) + U[2*node  ] * magna_disp;
		}
		plotElem();
	}
	
// Compute stresses at element guassian points
	fout << "\n\nAverage Stresses and Strains at Element Centroid";
	fout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	fout << "\n\n  Elem       Sig_x         Sig_y          Sig_xy ";
	fout << "\n================================================";
	double pxi, pet;
	for (n=1; n<=nElems; ++n)
	{
		thickness = thick[n];
		if (probType == 0) 
		{
			// plane stress
			D(1,1) = D(2,2) =  E[n]/(1 - nu[n]*nu[n]);
			D(1,2) = D(2,1) = D(1,1) * nu[n];
			D(3,3) = D(1,1) * (1 - nu[n])/2.0;
		}
		else
		{
			// plane strain
			double c;
			c = E[n]/(1 + nu[n])/(1 - 2*nu[n]);
			D(1,1) = D(2,2) = c*(1 - nu[n]);
			D(1,2) = D(2,1) = c * nu[n];
			D(3,3) = c * (1 - 2*nu[n])/2.0;
		}
		D(1,3) = D(2,3) = D(3,1) = D(3,2) = 0;
		for (ii=1; ii<=4; ++ii)
		{
			node = elemConn(ii, n);
			xl[ii] = X(1, node);
			yl[ii] = X(2, node);
		}
		pxi = pet = 0; // = place[1][1]

		// Get element displacement vector
		int dof, dof1, node;
		dof1 = 0;
		for  (ii=1; ii<=4; ++ii)
		{
			node = elemConn(ii, n);
			for (j=1; j<=2; ++j)
			{
				dof1 ++;
				dof = 2*(node - 1) + j;
				u[dof1] = U[dof];
			}
		}
		quad4Stress(n, thickness, pxi, pet);
		fout << "\n" << setw(4) << n ;
		for (i=1; i<=3; i++)
			fout << setw(15) << stress[i];
	}
}

void erase(void)
{
	ofstream fout(file);
}
void quad4(int elem, double thick)
{
//Refer RD Cook et al. 
	int na, nb; 
	double dv, pxi, pet;

// Clear load vector and upper triangle of stiffness matrix
	eLoad  = 0;
	eStiff = 0;

// Start Gauss quadrature loop. Use nGauss by nGauss rule
	for (na=1; na<=nGauss; ++na)
	{
		pxi = place[na][nGauss];
		for (nb=1; nb<=nGauss; ++nb)
		{
			pet = place[nb][nGauss];
			shapeFunc(pxi,pet);
			dv = wgt[na][nGauss] * wgt[nb][nGauss] * thick * detJac;
			eStiff += B.transpose() * D * B * dv;
		}
	}
}

void shapeFunc(double pxi, double pet)
{
	dVector N(4), Nxi(4), Net(4);
	double dum;
	int i, l, j, k; 
	int xii[5] = {0, -1, 1, 1, -1};
	int eti[5] = {0, -1, -1, 1, 1};
// Shape functions and their derivatives
	for (i=1; i<=4; ++i)
	{
		double pX = 0.25 * (1 + pxi*xii[i]);
		double pT = 0.25 * (1 + pet*eti[i]);
		N[i]   = 4 * pX * pT;
		Nxi[i] = xii[i] * pT;
		Net[i] = eti[i] * pX;
	}
	
	// Clear array Jac and B
	Jac = 0;
	B   = 0;
	
	// Find Jacobian and its determinant. 
	Jac(1,1) += dot(Nxi,xl);
	Jac(1,2) += dot(Nxi,yl);
	Jac(2,1) += dot(Net,xl);
	Jac(2,2) += dot(Net,yl);

	detJac = Jac(1,1)*Jac(2,2) - Jac(1,2)*Jac(2,1);
	
	// Replace Jac by its inverse.
	dum = Jac(1,1)/detJac;
	Jac(1,1) = Jac(2,2)/detJac;
	Jac(1,2) = - Jac(1,2)/detJac;
	Jac(2,1) = - Jac(2,1)/detJac;
	Jac(2,2) = dum;
	
	// Form [B] matrix (zero entries are already set)
	for (j=1; j<=4; ++j)
	{
		l = 2 * j;
		k = l - 1;
		B(1,k) = Jac(1,1)*Nxi[j] + Jac(1,2)*Net[j];
		B(2,l) = Jac(2,1)*Nxi[j] + Jac(2,2)*Net[j];
		B(3,k) = B(2,l);
		B(3,l) = B(1,k);
	}
}

void quad4Stress(int elem, double thick, double pxi, double pet)
{
// Get stresses and strains
	shapeFunc(pxi, pet);

// Compute stress as DB * u and strains as B * u
	strain = B * u;
	stress = D * strain;
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
	fout << "\nx = [" << xl[3] << " " << xl[4] << "];";
	fout << "\ny = [" << yl[3] << " " << yl[4] << "];";
	fout << "\nplot(x,y,'" << color <<"','Linewidth',2)";
	fout << "\nhold on";
	fout << "\nx = [" << xl[4] << " " << xl[1] << "];";
	fout << "\ny = [" << yl[4] << " " << yl[1] << "];";
	fout << "\nplot(x,y,'" << color <<"','Linewidth',2)";
	fout << "\nhold on";
	fout << "\naxis off";
}
