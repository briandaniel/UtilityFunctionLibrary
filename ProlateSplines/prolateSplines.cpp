/*
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Copyright (C) 2017-2019 Brian D. Hong
 *
*/

/*
 * prolateSplines.cpp
 *
 *  Created on: Sep 28, 2017
 *      Author: brian
 *
 *   This implementation gives functions for 2D Splines in prolate spheroidal coordinates.
 *
 *   1. The splines are functions of (nu,phi)
 *   2. This was orginally used to define surfaces/deformation functions for the prolate spheroidal left ventricle geometry
 *
 */

#include "prolateSplines.hpp"

ProlateSplines::ProlateSplines( int NnuSplinesIn, int NphiSplinesIn, double nuMinIn, double muConstIn, double * s ){

	setSplineSize( NnuSplinesIn, NphiSplinesIn, nuMinIn, muConstIn );
	computeSplineParameters( s );
}

void ProlateSplines::computeConstantSpline( double f0 )
{

	// this function just sets the values to zero except the constant values which are set to f
	double * s = new double[NSplineParams];

	setConstParamValues( f0, s, NnuSplines, NSplineParams, NSplines );

	computeSplineParameters( s );

	delete [] s;
}

void ProlateSplines::createConstantProlateSplines( double f0, int NnuSplinesIn, int NphiSplinesIn, double nuMinIn, double muConstIn ){

	setSplineSize( NnuSplinesIn, NphiSplinesIn, nuMinIn, muConstIn );
	computeConstantSpline(f0);

}

void ProlateSplines::createRegularSplines( double * s, int NnuSplinesIn, int NphiSplinesIn, double nuMinIn, double muConstIn ){


	setSplineSize( NnuSplinesIn, NphiSplinesIn, nuMinIn, muConstIn );
	computeSplineParameters( s );

}

double ProlateSplines::evaluateSplines( double nu0, double phi0 )
{
	double pVal = 0.0;

	pVal = evaluateSpheroidSplines( nu0, phi0, pp, nubnd, phibnd, NSplines);

	return pVal;
}

void ProlateSplines::setSplineSize( int NnuSplinesIn, int NphiSplinesIn, double nuMinIn, double muConstIn )
{
	nuMin = nuMinIn;
	muConst = muConstIn;

	NnuSplines = NnuSplinesIn;
	NphiSplines = NphiSplinesIn;

	NSplines = NnuSplines*NphiSplines;

	NSplineParams = 3 + 4*NnuSplines*NphiSplines;

	pp.assign(NSplines, vector<double>(16,0) );
	nubnd.assign(NSplines, vector<double>(2,0) );
	phibnd.assign(NSplines, vector<double>(2,0) );
}


void ProlateSplines::computeSplineParameters( double * s ){

	double * nuVecSplines, * phiVecSplines;

	// Compute the initial surface spline coefficients
	nuVecSplines = new double[NnuSplines + 1];
	phiVecSplines = new double[NphiSplines + 1];
	computeSplineNuPhiVectors( muConst, nuMin, NnuSplines, NphiSplines, nuVecSplines, phiVecSplines);

	computeSpheroidSplineParams( s, NnuSplines, NphiSplines, nuVecSplines, phiVecSplines, pp, nubnd, phibnd);

	delete [] nuVecSplines;
	delete [] phiVecSplines;

}

void ProlateSplines::evaluateSplineDerivatives( double nu0, double phi0, double &dfdnu0, double &dfdphi0 ){


	evaluateSpheroidSplineDerivatives( nu0, phi0, pp, nubnd, phibnd, NSplines, dfdnu0, dfdphi0);
}

double ProlateSplines::evaluateSplineCrossDerivative( double nu0, double phi0  ){

	double dfdnu0dphi0 = evaluateSpheroidSplineCrossDerivative( nu0, phi0, pp, nubnd, phibnd, NSplines );

	return dfdnu0dphi0;

}

void ProlateSplines::evaluateSplineSecondDerivatives( double nu0, double phi0, double & d2fdnu2, double & d2fdnudphi, double & d2fdphi2 ){

	 evaluateSpheroidSplineSecondDerivatives( nu0, phi0, pp, nubnd, phibnd, NSplines, d2fdnu2, d2fdnudphi, d2fdphi2 );

}

double ProlateSplines::evaluateLaplaceSqInt(){

	double integral;
	integral = evaluateSpheroidSplineLaplacianSquaredIntegral( pp, nubnd, phibnd, NnuSplines, NphiSplines );

	return integral;

}

void computeSplineNuPhiVectors( double mu_const, double nuUp0, int NnuSplines, int NphiSplines, double * nuVec, double * phiVec)
{
	double localIntegral;
	double integrand, integrandPrev;
	double integral;

	int Nnu = NnuSplines+1;

	int Nnu_integral = 100;

	double * nu = new double [Nnu_integral];
	double dnu = (PI - nuUp0)/(Nnu_integral-1);
	double * d = new double [Nnu_integral];

	d[0] = 0; // first value is zero
	nu[0] = nuUp0; // first value of nu is nuUp0
	integrand = sqrt( pow(sinh(mu_const),2) + pow( sin(nu[0]) ,2) );

	integral = 0;
	for(int k = 1; k < Nnu_integral; k++)
	{
		integrandPrev = integrand;

		nu[k] = nu[k-1] + dnu;

		integrand = sqrt( pow(sinh(mu_const),2) + pow( sin(nu[k]) ,2) );

		localIntegral = dnu*( integrand + integrandPrev )/2.0;

		integral = integral + localIntegral;

		d[k] = integral;

	}


	double * dNode = new double [Nnu];
	linspace( 0.0, d[Nnu_integral-1], Nnu, dNode );


	// compute nuVec values
	for(int k = 0; k < Nnu-1; k ++ )
	{
		nuVec[k] = linearInterp( d, nu, Nnu_integral, dNode[k] );
	}
	nuVec[Nnu-1] = PI; // last value is pi, so this is to prevent roundoff error.

	// compute phiVec values
	linspace( 0, 2*PI, NphiSplines+1, phiVec );

	delete [] nu;
	delete [] d;
	delete [] dNode;

}



/* the parameters are read from "s" as follows:
 *
 * 1. The first three parameters determine (i) the value at the nu = pi and (ii) the slopes at nu = pi
 * 	  these are independent of the grid.
 *
 * 2. The next NnuSplines*NphiSplines parameters determine the values of the function at the nodes
 *	  the ordering starts at the smallest value of nu (nuUp0), traversing the phi direction first.
 *	  The value then shifts to the next nu value, and runs through the phi variable again.
 *	  (the periodic value of phi is excluded from s as it would be redundant)
 *
 * 3. The next NnuSplines*NphiSplines parameters are the dFdnu derivatives in the same order
 *
 * 4. The next NnuSplines*NphiSplines parameters are the dFdphi derivatives in the same order
 *
 * 5. The next NnuSplines*NphiSplines parameters are the d2Fdnudphi derivatives in the same order
 *
 *
 */
void computeSpheroidSplineParams( double * s, int NnuSplines, int NphiSplines, double * nuVec, double * phiVec,
		vector< vector<double> > & pp, vector< vector<double> > & nu0bnd, vector< vector<double> > & phi0bnd)
{
	double** nuMesh, ** phiMesh, ** F, ** dFdnu, ** dFdphi, ** d2Fdnudphi;
	int i, j, k;

	// Grid size is one more than the number of splines
	int Nnu = NnuSplines + 1;
	int Nphi = NphiSplines + 1;


	nuMesh = new double * [ Nnu ];
	phiMesh = new double * [ Nnu ];
	F = new double * [ Nnu ];
	dFdnu= new double * [ Nnu ];
	dFdphi  = new double * [ Nnu ];
	d2Fdnudphi = new double * [ Nnu ];

	for (int k = 0; k < Nnu; k++)
	{
		nuMesh[k] = new double [ Nphi ];
		phiMesh[k] = new double [Nphi];
		F[k] = new double [ Nphi ];
		dFdnu[k] = new double [ Nphi ];
		dFdphi[k] = new double [ Nphi ];
		d2Fdnudphi[k] = new double [ Nphi ];

	}


	for ( k = 0; k < Nnu; k++)
	{
		for ( j = 0; j < Nphi; j++)
		{
			nuMesh[k][j] = nuVec[k];
			phiMesh[k][j] = phiVec[j];

		}
	}

	// apex computations
	i = Nnu-1;
	for( j = 0; j < Nphi; j++)
	{
		F[i][j] = s[0];
		dFdnu[i][j] = s[1]*cos(phiVec[j]) + s[2]*sin(phiVec[j]);
		dFdphi[i][j] = 0;
		d2Fdnudphi[i][j] = s[2]*cos(phiVec[j]) - s[1]*sin(phiVec[j]);
	}


	// Regular values
	for( i = 0; i < Nnu-1; i++)
	{
		for( j = 0; j < Nphi - 1; j++)
		{
			F[i][j] = s[3 + j + (Nphi-1)*i];
			dFdnu[i][j] = s[3 + j + (Nphi-1)*i + NnuSplines*NphiSplines];
			dFdphi[i][j] = s[3 + j + (Nphi-1)*i + 2*NnuSplines*NphiSplines];
			d2Fdnudphi[i][j] = s[3 + j + (Nphi-1)*i + 3*NnuSplines*NphiSplines];
		}
	}

	// periodic extension
	j = Nphi - 1;
	for ( k = 0; k < Nnu; k++)
	{
		F[k][j] = F[k][0];
		dFdphi[k][j] = dFdphi[k][0];
		dFdnu[k][j] = dFdnu[k][0];
		d2Fdnudphi[k][j] = d2Fdnudphi[k][0];
	}

	computeSpheroidSplineCoef( nuMesh, phiMesh, F, dFdnu, dFdphi, d2Fdnudphi, pp, nu0bnd, phi0bnd, NnuSplines, NphiSplines);

	// Clean up
	for ( k = 0; k < Nnu; k++)
	{
		delete [] nuMesh[k];
		delete [] phiMesh[k];
		delete [] F[k];
		delete [] dFdnu[k];
		delete [] dFdphi[k];
		delete [] d2Fdnudphi[k];
	}

	delete [] nuMesh;
	delete [] phiMesh;
	delete [] F;
	delete [] dFdnu;
	delete [] dFdphi;
	delete [] d2Fdnudphi;

}

void computeSpheroidSplineCoef( double** nuMesh, double** phiMesh, double** F, double** dFdnu, double** dFdphi, double ** d2Fdnudphi,
		vector< vector<double> > & pp, vector< vector<double> > & nu0bnd, vector< vector<double> > & phi0bnd, int NnuSplines, int NphiSplines )
{

	// Define the static matrix Minv
	double**Minv = new double * [16];
	for(int i = 0; i <16; i++)
		Minv[i] = new double[16];

	double MinvStatic[16][16] = {   { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
									{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
									{-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
									{ 2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
									{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
									{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
									{ 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0},
									{ 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0},
									{-3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0},
									{ 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0},
									{ 9,-9,-9, 9, 6, 3,-6,-3, 6,-6, 3,-3, 4, 2, 2, 1},
									{-6, 6, 6,-6,-3,-3, 3, 3,-4, 4,-2, 2,-2,-2,-1,-1},
									{ 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
									{ 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0},
									{-6, 6, 6,-6,-4,-2, 4, 2,-3, 3,-3, 3,-2,-1,-2,-1},
									{ 4,-4,-4, 4, 2, 2,-2,-2, 2,-2, 2,-2, 1, 1, 1, 1} 	};
	for(int i = 0; i <16; i++)
	{
		for(int j = 0; j<16; j++)
		{
			Minv[i][j] = MinvStatic[i][j];
		}
	}


	int k;
	double nu0stretch, phi0stretch;
	double * fk = new double[16]; // fixed array for rhs


	for (int i = 0; i < NnuSplines; i++ )
	{
		for (int j = 0; j < NphiSplines; j++ )
		{

			k = i*NphiSplines + j; // Polynomial index

			// Boundaries of the current (k) polynomial
			nu0bnd[k][0] = nuMesh[i][j];
			nu0bnd[k][1] = nuMesh[i+1][j];
			phi0bnd[k][0] = phiMesh[i][j];
			phi0bnd[k][1] = phiMesh[i][j+1];

			// The derivatives are normalized so these provide the correct stretch factor
			nu0stretch = nu0bnd[k][1] - nu0bnd[k][0];
			phi0stretch = phi0bnd[k][1] - phi0bnd[k][0];

			// Compute the RHS
			fk[0] = F[i][j];
			fk[1] = F[i+1][j];
			fk[2] = F[i][j+1];
			fk[3] = F[i+1][j+1];
			fk[4] = dFdnu[i][j]*nu0stretch;
			fk[5] = dFdnu[i+1][j]*nu0stretch;
			fk[6] = dFdnu[i][j+1]*nu0stretch;
			fk[7] = dFdnu[i+1][j+1]*nu0stretch;
			fk[8] = dFdphi[i][j]*phi0stretch;
			fk[9] = dFdphi[i+1][j]*phi0stretch;
			fk[10] = dFdphi[i][j+1]*phi0stretch;
			fk[11] = dFdphi[i+1][j+1]*phi0stretch;
			fk[12] = d2Fdnudphi[i][j]*phi0stretch*nu0stretch;
			fk[13] = d2Fdnudphi[i+1][j]*phi0stretch*nu0stretch;
			fk[14] = d2Fdnudphi[i][j+1]*phi0stretch*nu0stretch;
			fk[15] = d2Fdnudphi[i+1][j+1]*phi0stretch*nu0stretch;

			// Compute the coefficients
			matrixVectorMultiply( Minv, fk, 16, 16, pp[k].data() );

		}
	}

	// clean up
	delete [] fk;
	for(int i = 0; i <16; i++)
		delete [] Minv[i];
	delete [] Minv;

}

double evaluateSpheroidSplines( double nu0, double phi0, vector< vector<double> > & pp,
		vector< vector<double> > & nu0bnd, vector< vector<double> > & phi0bnd, int Nspline)
{
	double x,y;
	double f = 0;
	int k = 0; // which spline to use

	// domain safeties
	phi0 = mod(phi0, 2*PI);
	if(nu0 > PI){ nu0 = PI; }
	if(nu0 < nu0bnd[0][0]){  nu0 = nu0bnd[0][0]; };


	// Determine which spline to evaluate (this could be computed directly because the spacing is fixed)
	for(int i = 0; i < Nspline; i++)
	{
		if( nu0 >= nu0bnd[i][0] && nu0 <= nu0bnd[i][1] )
		{
			if( phi0 >= phi0bnd[i][0] && phi0 <= phi0bnd[i][1] )
			{
				k = i;
			}
		}
	}

	// Evaluate the spline
	x = (nu0 - nu0bnd[k][0])/(nu0bnd[k][1] - nu0bnd[k][0]);
	y = (phi0 - phi0bnd[k][0])/(phi0bnd[k][1] - phi0bnd[k][0]);
	f =             ( pp[k][0]  + pp[k][1]*x  + pp[k][2]*pow(x,2)  + pp[k][3]*pow(x,3)  )
		 +        y*( pp[k][4]  + pp[k][5]*x  + pp[k][6]*pow(x,2)  + pp[k][7]*pow(x,3)  )
		 + pow(y,2)*( pp[k][8]  + pp[k][9]*x  + pp[k][10]*pow(x,2) + pp[k][11]*pow(x,3) )
		 + pow(y,3)*( pp[k][12] + pp[k][13]*x + pp[k][14]*pow(x,2) + pp[k][15]*pow(x,3) );


	return f;

}

void evaluateSpheroidSplineDerivatives( double nu0, double phi0, vector< vector<double> > & pp,
		vector< vector<double> > & nu0bnd, 	vector< vector<double> > & phi0bnd, int Nspline, double &dfdnu0, double &dfdphi0)
{
	double x,y;
	double dfdunuk = 0;
	double dfduphik = 0;
	int k = 0; // which spline to use


	// domain safeties
	phi0 = mod(phi0, 2*PI);
	if(nu0 > PI){ nu0 = PI; }



	// Determine which spline to evaluate (this could be computed directly because the spacing is fixed)
	for(int i = 0; i < Nspline; i++)
	{
		if( nu0 >= nu0bnd[i][0] && nu0 <= nu0bnd[i][1] )
		{
			if( phi0 >= phi0bnd[i][0] && phi0 <= phi0bnd[i][1] )
			{
				k = i;
			}
		}
	}

	// Evaluate the spline
	x = (nu0 - nu0bnd[k][0])/(nu0bnd[k][1] - nu0bnd[k][0]);
	y = (phi0 - phi0bnd[k][0])/(phi0bnd[k][1] - phi0bnd[k][0]);
	dfdunuk = pp[k][1] + 2*pp[k][2]*x + 3*pp[k][3]*pow(x,2)
			 + y*( pp[k][5] + 2*pp[k][6]*x + 3*pp[k][7]*pow(x,2) )
			 + pow(y,2)*( pp[k][9] + 2*pp[k][10]*x + 3*pp[k][11]*pow(x,2) )
			 + pow(y,3)*( pp[k][13] + 2*pp[k][14]*x + 3*pp[k][15]*pow(x,2) );

	dfduphik =  pp[k][4] + pp[k][5]*x + pp[k][6]*pow(x,2) + pp[k][7]*pow(x,3)
				 + 2*y*( pp[k][8] + pp[k][9]*x+ pp[k][10]*pow(x,2) + pp[k][11]*pow(x,3) )
				 + 3*pow(y,2)* (pp[k][12] + pp[k][13]*x + pp[k][14]*pow(x,2) + pp[k][15]*pow(x,3) );

	dfdnu0 = dfdunuk/(nu0bnd[k][1] - nu0bnd[k][0]);
	dfdphi0 = dfduphik/(phi0bnd[k][1] - phi0bnd[k][0]);

}


double evaluateSpheroidSplineCrossDerivative( double nu0, double phi0, vector< vector<double> > & pp,
		vector< vector<double> > & nu0bnd, vector< vector<double> > & phi0bnd, int Nspline )
{
	double dfdnu0dphi0;

	double x,y;
	double dfdunuk = 0;
	double dfdunuduphik = 0;
	int k = 0; // which spline to use

	// domain safeties
	phi0 = mod(phi0, 2*PI);
	if(nu0 > PI){ nu0 = PI; }

	// Determine which spline to evaluate (this could be computed directly because the spacing is fixed)
	for(int i = 0; i < Nspline; i++)
	{
		if( nu0 >= nu0bnd[i][0] && nu0 <= nu0bnd[i][1] )
		{
			if( phi0 >= phi0bnd[i][0] && phi0 <= phi0bnd[i][1] )
			{
				k = i;
			}
		}
	}

	// Evaluate the spline
	x = (nu0 - nu0bnd[k][0])/(nu0bnd[k][1] - nu0bnd[k][0]);
	y = (phi0 - phi0bnd[k][0])/(phi0bnd[k][1] - phi0bnd[k][0]);

	dfdunuduphik = ( pp[k][5] + 2*pp[k][6]*x + 3*pp[k][7]*pow(x,2) )
					 + 2*y*( pp[k][9] + 2*pp[k][10]*x + 3*pp[k][11]*pow(x,2) )
					 + 3*pow(y,2)*( pp[k][13] + 2*pp[k][14]*x + 3*pp[k][15]*pow(x,2) );

	dfdnu0dphi0 = dfdunuduphik/( (nu0bnd[k][1] - nu0bnd[k][0])*(phi0bnd[k][1] - phi0bnd[k][0]) );

	return dfdnu0dphi0;

}

void evaluateSpheroidSplineSecondDerivatives( double nu0, double phi0, vector< vector<double> > & pp, vector< vector<double> > & nu0bnd,
		vector< vector<double> > & phi0bnd, int Nspline, double & d2fdnu2, double & d2fdnudphi, double & d2fdphi2 )
{

	double x,y;
	double d2fdunu2 = 0;
	double d2fdunuduphi = 0;
	double d2fduphi2 = 0;

	double nuDiff;
	double phiDiff;

	int k = 0; // which spline to use


	// domain safeties
	phi0 = mod(phi0, 2*PI);
	if(nu0 > PI){ nu0 = PI; }


	// Determine which spline to evaluate (this could be computed directly because the spacing is fixed)
	for(int i = 0; i < Nspline; i++)
	{
		if( nu0 >= nu0bnd[i][0] && nu0 <= nu0bnd[i][1] )
		{
			if( phi0 >= phi0bnd[i][0] && phi0 <= phi0bnd[i][1] )
			{
				k = i;
			}
		}
	}

	// Evaluate the spline
	nuDiff = nu0bnd[k][1] - nu0bnd[k][0];
	phiDiff = phi0bnd[k][1] - phi0bnd[k][0];
	x = (nu0 - nu0bnd[k][0])/(nuDiff);
	y = (phi0 - phi0bnd[k][0])/(phiDiff);

	d2fdunu2 = 2*pp[k][2] + 6*pp[k][3]*x
			 + y*( 2*pp[k][6] + 6*pp[k][7]*x )
			 + pow(y,2)*( 2*pp[k][10] + 6*pp[k][11]*x )
			 + pow(y,3)*( 2*pp[k][14] + 6*pp[k][15]*x );

	d2fdunuduphi = ( pp[k][5] + 2*pp[k][6]*x + 3*pp[k][7]*pow(x,2) )
					 + 2*y*( pp[k][9] + 2*pp[k][10]*x + 3*pp[k][11]*pow(x,2) )
					 + 3*pow(y,2)*( pp[k][13] + 2*pp[k][14]*x + 3*pp[k][15]*pow(x,2) );

	d2fduphi2= 2*( pp[k][8] + pp[k][9]*x+ pp[k][10]*pow(x,2) + pp[k][11]*pow(x,3) )
			   + 6*y* (pp[k][12] + pp[k][13]*x + pp[k][14]*pow(x,2) + pp[k][15]*pow(x,3) );


	d2fdnu2 = d2fdunu2/( pow( nuDiff, 2 ) );
	d2fdnudphi= d2fdunuduphi/( nuDiff*phiDiff );
	d2fdphi2 = d2fduphi2/( pow( phiDiff, 2 ) );

}

double evaluateSpheroidSplineLaplacianSquaredIntegral( vector< vector<double> > & pp,
		vector< vector<double> > & nu0bnd,	vector< vector<double> > & phi0bnd, int NnuSplines, int NphiSplines )
{


	double p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16;
	double a, b;
	int k;
	double Istraight, Icross;

	double integral = 0;

	for(int i = 0; i < NnuSplines; i++)
	{
		for(int j = 0; j < NphiSplines; j++)
		{

			k = i*NphiSplines + j; // Polynomial index

			p1 = pp[k][0]; p2 = pp[k][1]; p3 = pp[k][2]; p4 = pp[k][3];
			p5 = pp[k][4]; p6 = pp[k][5]; p7 = pp[k][6]; p8 = pp[k][7];
			p9 = pp[k][8]; p10 = pp[k][9]; p11 = pp[k][10]; p12 = pp[k][11];
			p13 = pp[k][12]; p14 = pp[k][13]; p15 = pp[k][14]; p16 = pp[k][15];


			// Scale factors
			// The derivatives are normalized so these provide the correct stretch factor
			a = 1/(nu0bnd[k][1] - nu0bnd[k][0]);
			b = 1/(phi0bnd[k][1] - phi0bnd[k][0]);

			// This is the integral of the second partials squared
	        Istraight = 1/105*pow(a,4)*( 84*pow(p11,2) + 252*pow(p12,2) + 60*pow(p15,2) + 180*p15*p16 +
	                					 180*pow(p16,2) + 210*p15*p3 + 315*p16*p3 + 420*pow(p3,2) + 315*p15*p4 +
										 630*p16*p4 + 1260*p3*p4 + 1260*pow(p4,2) + 168*p15*p7 + 252*p16*p7 +
										 420*p3*p7 + 630*p4*p7 + 140*pow(p7,2) +
										 42*(6*p15 + 12*p16 + 15*p3 + 30*p4 + 10*p7)*p8 + 420*pow(p8,2) +
										 105*p12*(2*p15 + 4*p16 + 4*p3 + 8*p4 + 3*p7 + 6*p8) +
										 7*p11*(36*p12 + 5*(4*p15 + 6*p16 + 8*p3 + 12*p4 + 6*p7 + 9*p8))) +
										 pow(b,4)*((4*pow(p10,2))/3 + (4*pow(p11,2))/5 + (4*pow(p12,2))/7 + 3*p12*p13 +
										 12*pow(p13,2) + (12*p12*p14)/5 + 12*p13*p14 + 4*pow(p14,2) + 2*p12*p15 +
										 8*p13*p15 + 6*p14*p15 + (12*pow(p15,2))/5 + (12*p12*p16)/7 +
										 6*p13*p16 + (24*p14*p16)/5 + 4*p15*p16 + (12*pow(p16,2))/7 +
										 2*p12*p9 + 12*p13*p9 + 6*p14*p9 + 4*p15*p9 + 3*p16*p9 + 4*pow(p9,2) +
										 p11*((4*p12)/3 + 4*p13 + 3*p14 + (12*p15)/5 + 2*p16 + (8*p9)/3) +
										 p10*(2*p11 + (8*p12)/5 + 6*p13 + 4*p14 + 3*p15 + (12*p16)/5 + 4*p9));


	        // This is the integral of 2* the product of the second partials
	        Icross = 1/225*pow(a,2)*pow(b,2)*(200*pow(p11,2) + 75*p10*(4*p11 + 8*p12 +
												3*(p15 + 2*(p16 + 2*p3 + 4*p4 + p7 + 2*p8))) +
												75*p11*(8*p12 + 18*p13 + 9*p14 + 8*p15 + 9*p16 + 8*p3 + 18*p4 +
												4*p7 + 9*p8 + 8*p9) +
												3*(120*pow(p12,2) + 180*p14*p15 + 120*pow(p15,2) + 360*p14*p16 +
												360*p15*p16 + 216*pow(p16,2) + 450*p14*p3 + 300*p15*p3 +
												225*p16*p3 + 900*p14*p4 + 675*p15*p4 + 540*p16*p4 +
												300*p14*p7 + 200*p15*p7 + 150*p16*p7 + 600*p14*p8 +
												450*p15*p8 + 360*p16*p8 +
												30*p13*(12*p15 + 18*p16 + 5*(6*p3 + 9*p4 + 4*p7 + 6*p8)) +
												75*(2*p15 + 3*p16 + 8*p3 + 12*p4 + 4*p7 + 6*p8)*p9 +
												15*p12*(45*p13 + 30*p14 + 25*p15 + 24*p16 + 10*p3 + 24*p4 +
												5*p7 + 12*p8 + 20*p9)));

	        integral = integral + (Istraight + Icross)/(a*b);
		}
	}
	return integral;
}

void setConstParamValues( double f0, double * s, int NnuSplines, int NSplineParams, int NSplines )
{

	// this function just sets the values to zero except the constant values which are set to f0
	for(int k = 0; k < NSplineParams; k++)
	{
		s[k] = 0.0;
	}

	// set constant
	s[0] = f0;

	for(int k = 3; k < 3+NSplines; k ++ ){
		s[k] = f0;
	}

}


