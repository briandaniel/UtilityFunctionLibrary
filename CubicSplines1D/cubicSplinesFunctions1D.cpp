/*
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Copyright (C) 2017-2019 Brian D. Hong
 *
*/

/*
 * cubicSplinesFunctions1D.cpp
 *
 *  Created on: Nov 22, 2017
 *      Author: brian
 *
 * A simple implementation of 1D cubic splines
 *	1. Natural spline conditions
 *	2. Periodic conditions
 *	3. Specified derivatives
 *	4. Uneven spacing
 *
 */


#include "cubicSplinesFunctions1D.hpp"

// This sets the size of the splines, even if they have already been set previously
void Spline1D::resetSplineSize( int NnodesIn)
{
	// only reset the size if they are different
	Nnodes = NnodesIn;
	Nsplines = Nnodes - 1;

	xbnd.resize(Nnodes);
	a.resize(Nsplines);
	b.resize(Nsplines);
	c.resize(Nsplines);
	d.resize(Nsplines);

}



// this computes spline coefficients for an uneven grid
// with zero second derivative boundary condition
void Spline1D::computeNaturalSplineCoef( int NnodesIn, double * x, double * y )
{
	// reset the size of the splines
	resetSplineSize(NnodesIn);

	// set spline boundary positions
	for(int i = 0; i < Nnodes; i++)
		xbnd[i] = x[i];

	int i;
	double * dydx = new double [Nnodes];
	double * rhs = new double [Nnodes];
	double * dx = new double [Nnodes-1];
	double ** M = new double * [Nnodes];
	for( i = 0; i < Nnodes; i++ )
		M[i] = new double[Nnodes];


	for ( i = 0; i < Nnodes-1; i++ )
	{
		dx[i] = xbnd[i+1] - xbnd[i];
	}

	for( i = 0; i < Nnodes; i++)
	{
		for(int j = 0; j < Nnodes ; j++)
		{
			M[i][j] = 0.0;
		}
	}


	// regular domain
	for( i = 1; i < Nnodes - 1; i++)
	{
		M[i][i-1] = dx[i-1];
		M[i][i] = 2*dx[i-1] + 2*dx[i];
		M[i][i+1] = dx[i];

		rhs[i] = 3.0*( y[i+1] - y[i-1] );
	}

	// left endpoint
	i = 0;
	M[i][i] = 2*dx[i];
	M[i][i+1] = dx[i];
	rhs[i] = 3.0*(y[i+1] - y[i]);

	// right endpoint
	i = Nnodes-1;
	M[i][i-1] = dx[i-1];
	M[i][i] = 2*dx[i-1];
	rhs[i] = 3.0*( y[i] - y[i-1] );

	// print2DArray( M, Nnodes, Nnodes, 3);
	// print1DArray( rhs, Nnodes, 3);


	// solves the system Ax = b by lu factorization
	luSolve( M, rhs, Nnodes, dydx );

	computeSplineCoef( Nnodes, x, y, dydx );


	// clean up
	for( int i = 0; i < Nnodes; i++ )
		delete [] M[i];
	delete [] M;

	delete [] dx;
	delete [] rhs;
	delete [] dydx;


}


// this computes spline coefficients for an uneven grid
// this routine ASSUMES that the LAST POINT IN Y IS THE SAME as the FIRST POINT
// in other words that the periodic point is INCLUDED in the input
void Spline1D::computeNaturalPeriodicSplineCoef( int NnodesIn, double * x, double * y )
{
	// reset the size of the splines
	resetSplineSize(NnodesIn);


	// the actual number of nodes in the fitted periodic points is one less
	// to avoid double counting the periodic point
	int NnodesPeriodic = Nnodes-1;

	// set spline boundary positions
	for(int i = 0; i < Nnodes; i++)
		xbnd[i] = x[i];

	int i;
	double * dydx = new double [Nnodes];

	double * rhs = new double [NnodesPeriodic];
	double * dx = new double [Nsplines];
	double * dydxPeriodic = new double [NnodesPeriodic];
	double ** M = new double * [NnodesPeriodic];
	for( i = 0; i < NnodesPeriodic; i++ )
		M[i] = new double[NnodesPeriodic];


	for ( i = 0; i < Nsplines; i++ )
	{
		dx[i] = xbnd[i+1] - xbnd[i];
	}

	for( i = 0; i < NnodesPeriodic; i++)
	{
		for(int j = 0; j < NnodesPeriodic ; j++)
		{
			M[i][j] = 0.0;
		}
	}


	// regular domain
	for( i = 1; i < NnodesPeriodic-1; i++)
	{
		M[i][i-1] = dx[i-1];
		M[i][i] = 2*dx[i-1] + 2*dx[i];
		M[i][i+1] = dx[i];

		rhs[i] = 3.0*( y[i+1] - y[i-1] );
	}

	// periodic extensions
	i = 0;
	M[i][NnodesPeriodic-1] = dx[Nsplines-1];
	M[i][i] = 2*dx[Nsplines-1] + 2*dx[i];
	M[i][i+1] = dx[i];
	rhs[i] = 3.0*(y[i+1] - y[NnodesPeriodic-1]);

	i = NnodesPeriodic-1;
	M[i][i-1] = dx[Nsplines-1];
	M[i][i] = 2*dx[Nsplines-1] + 2*dx[0];
	M[i][0] = dx[0];
	rhs[i] = 3.0*( y[0] - y[i-1] );



	// solves the system Ax = b by lu factorization
	luSolve( M, rhs, NnodesPeriodic, dydxPeriodic );

	// extend dydx back to the original grid
	for(int k = 0; k < NnodesPeriodic; k++)
		dydx[k] = dydxPeriodic[k];

	dydx[Nnodes-1] = dydxPeriodic[0]; // periodic extension

	computeSplineCoef( Nnodes, x, y, dydx );

	// clean up
	for( int i = 0; i < NnodesPeriodic; i++ )
		delete [] M[i];
	delete [] M;

	delete [] dx;
	delete [] rhs;
	delete [] dydx;
	delete [] dydxPeriodic;

}

// this computes spline coefficients for an uneven grid
// of nodes x with function values y and known derivatives dydx
void Spline1D::computeSplineCoef( int NnodesIn, double * x, double * y, double * dydx )
{
	// reset the size of the splines
	resetSplineSize(NnodesIn);

	// set spline boundary positions
	for(int i = 0; i < NnodesIn; i++)
		xbnd[i] = x[i];

	double dx;
	for (int i = 0; i < Nsplines; i++ )
	{


		dx = xbnd[i+1] - xbnd[i];

		a[i] = y[i];
		b[i] = dx*dydx[i];
		c[i] = 3*( y[i+1] - y[i] ) - dx*( dydx[i+1] + 2*dydx[i] );
		d[i] = dx*( dydx[i+1] + dydx[i] ) + 2*( y[i] - y[i+1] );

	}

}


// This evaluates the splines at the point x.
// Because the grid is (possibly) uneven, a search is used.
// if x is outside the domain p = 0
double Spline1D::evaluateSpline( double x ){

	// default value outside the domain is zero
	double t;
	double p = 0;

	for( int i = 0; i < Nsplines; i++ ){

		if ( x >= xbnd[i] && x <= xbnd[i+1] )
		{
			t = ( x - xbnd[i] )/( xbnd[i+1] - xbnd[i] );

			p = a[i] + b[i]*t + c[i]*t*t + d[i]*t*t*t;

			break;
		}
	}

	return p;

}


// This evaluates the derivative of the spline at point x
// if x is outside the domain dpdx = 0
double Spline1D::evaluateSplineDerivative( double x ){

	// default value outside the domain is zero
	double t;
	double dpdx = 0;
	double dpdt = 0;
	double dx;

	for( int i = 0; i < Nsplines; i++ ){

		if ( x >= xbnd[i] && x <= xbnd[i+1] )
		{
			dx = xbnd[i+1] - xbnd[i];
			t = ( x - xbnd[i] )/dx;

			dpdt = b[i] + 2.0*c[i]*t + 3.0*d[i]*t*t;

			dpdx = dpdt/dx;

			break;
		}
	}

	return dpdx;

}



// This evaluates the energy norm
// \int  ( d^2p/dx^2 )^2 dx
double Spline1D::evaluateSplineEnergy(){

	// default value outside the domain is zero
	double dx;

	double integralSum = 0.0;

	for( int i = 0; i < Nsplines; i++ ){

		dx = xbnd[i+1] - xbnd[i];

		double integralLocal = ( 4.0*pow(c[i],2) + 12.0*c[i]*d[i] + 12.0*pow(d[i],2) ) / pow(dx,3);

		integralSum = integralSum + integralLocal;

	}

	return integralSum;

}



// this finds the minimum value in the spline
void Spline1D::findSplineMinimumAndMaximum( double & xMin, double & xMax, double & pMin, double & pMax)
{
	pMin = 1e100;
	pMax = -1e100;

	double tvec [4];
	double pvec [4];

	for( int i = 0; i < Nsplines; i++ )
	{
		double x0 = xbnd[i];
		double dx = xbnd[i+1] - xbnd[i];

		// check locations
		tvec[0] = 0;
		tvec[1] = 1;
		tvec[2] = ( -2*c[i] + sqrt( 4.0*pow(c[i],2)- 12.0*b[i]*d[i] ) )/( 6.0*d[i] );
		if ( (tvec[2] < 0  || tvec[2] > 1 ) || tvec[2] != tvec[2]){ tvec[2] = 0; } // make sure its in the domain, if not just set it to 0

		tvec[3] = ( -2*c[i] - sqrt( 4.0*pow(c[i],2)- 12.0*b[i]*d[i] ) )/( 6.0*d[i] );
		if ( (tvec[3] < 0  || tvec[3] > 1 ) || tvec[3] != tvec[3]){ tvec[3] = 0; } // make sure its in the domain, if not just set it to 0


  		for(int j = 0; j < 4; j++ )
		{
			double x = x0 + tvec[j]*dx;
			double p = evaluateSpline(x);

			// check if max
			if( p > pMax )
			{
				pMax = p;
				xMax = x;
			}

			// check if min
			if( p < pMin )
			{
				pMin = p;
				xMin = x;
			}

		}
	}


}























