/*
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Copyright (C) 2017-2019 Brian D. Hong
 *
*/

/*
 * cubicSplinesFunctions1D.hpp
 *
 *  Created on: Nov 22, 2017
 *      Author: brian
 */

#ifndef CUBICSPLINESFUNCTIONS1D_HPP_
#define CUBICSPLINESFUNCTIONS1D_HPP_

// standard includes
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "../UtilityFunctions/utilityFunctions.hpp"

#ifndef PI
#define PI 3.141592653589793
#endif


class Spline1D {

  public:


	// This sets the size of the splines, even if they have already been set previously
	void resetSplineSize( int NnodesIn );

	// this computes spline coefficients for an uneven grid
	// with zero second derivative boundary condition
	void computeNaturalSplineCoef( int NnodesIn, double * x, double * y );

	// this computes spline coefficients for an uneven grid
	// this routine ASSUMES that the LAST POINT IN Y IS THE SAME as the FIRST POINT
	// in other words that the periodic point is INCLUDED in the input
	void computeNaturalPeriodicSplineCoef( int NnodesIn, double * x, double * y );

	// this computes spline coefficients for an uneven grid
	// of nodes x with function values y and known derivatives dydx
	void computeSplineCoef( int NnodesIn, double * x, double * y, double * dydx );

	// This evaluates the splines at the point x.
	// Because the grid is (possibly) uneven, a search is used.
	// if x is outside the domain p = 0
	double evaluateSpline( double x );

	// This evaluates the derivative of the spline at point x
	// if x is outside the domain dpdx = 0
	double evaluateSplineDerivative( double x );

	// This evaluates the energy norm
	// \int  ( d^2p/dx^2 )^2 dx
	double evaluateSplineEnergy();

	// this finds the minimum value in the spline
	void findSplineMinimumAndMaximum( double & xMin, double & xMax, double & pMin, double & pMax);

	// default constructor
	Spline1D()
	{
		Nnodes = -1;
		Nsplines = -1;
	};

	~Spline1D(){};


  private:

	int Nnodes, Nsplines;
	vector <double> a;
	vector <double> b;
	vector <double> c;
	vector <double> d;
	vector <double> xbnd;

};



#endif /* CUBICSPLINESFUNCTIONS1D_HPP_ */
