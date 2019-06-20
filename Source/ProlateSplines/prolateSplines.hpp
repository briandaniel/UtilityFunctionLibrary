/*
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Copyright (C) 2017-2019 Brian D. Hong
 *
*/

/*
 * prolateSplines.hpp
 *
 *  Created on: Sep 28, 2017
 *      Author: brian
 */

#ifndef SOURCE_PROLATESPLINES_HPP_
#define SOURCE_PROLATESPLINES_HPP_

#ifndef PI
#define PI 3.141592653589793
#endif

// required includes
#include <math.h>

// local includes
#include "../UtilityFunctions/utilityFunctions.hpp"


using namespace std;

// functions outside the class, used by some of the class functions
// this original implementation was un-packaged
void computeSplineNuPhiVectors( double mu_const, double nuUp0, int NnuSplines, int NphiSplines, double * nuVec, double * phiVec);
void computeSpheroidSplineParams( double * s, int NnuSplines, int NphiSplines, double * nuVec, double * phiVec,
		vector< vector<double> > & pp, vector< vector<double> > & nu0bnd, vector< vector<double> > & phi0bnd);
void computeSpheroidSplineCoef( double** nuMesh, double** phiMesh, double** F, double** dFdnu, double** dFdphi, double ** d2Fdnudphi,
		vector< vector<double> > & pp, vector< vector<double> > & nu0bnd, vector< vector<double> > & phi0bnd, int NnuSplines, int NphiSplines );
double evaluateSpheroidSplines( double nu0, double phi0, vector< vector<double> > & pp,
		vector< vector<double> > & nu0bnd, vector< vector<double> > & phi0bnd, int Nspline);
void evaluateSpheroidSplineDerivatives( double nu0, double phi0, vector< vector<double> > & pp,
		vector< vector<double> > & nu0bnd, 	vector< vector<double> > & phi0bnd, int Nspline, double &dfdnu0, double &dfdphi0);
double evaluateSpheroidSplineCrossDerivative( double nu0, double phi0, vector< vector<double> > & pp,
		vector< vector<double> > & nu0bnd, vector< vector<double> > & phi0bnd, int Nspline );
double evaluateSpheroidSplineLaplacianSquaredIntegral( vector< vector<double> > & pp,
		vector< vector<double> > & nu0bnd,	vector< vector<double> > & phi0bnd, int NnuSplines, int NphiSplines );
void evaluateSpheroidSplineSecondDerivatives( double nu0, double phi0, vector< vector<double> > & pp, vector< vector<double> > & nu0bnd,
		vector< vector<double> > & phi0bnd, int Nspline, double & d2fdnu2, double & d2fdnudphi, double & d2fdphi2 );
void setConstParamValues( double f0, double * s, int NnuSplines, int NSplineParams, int NSplines );


// prolate splines class if you want them packaged
class ProlateSplines{

  public:

	int NnuSplines;
	int NphiSplines;
	int NSplines;
	int NSplineParams;

	double nuMin;
	double muConst;

	// spline coefficients, etc.
	vector< vector<double> > pp;
	vector< vector<double> > nubnd;
	vector< vector<double> > phibnd;

	ProlateSplines(){
		NnuSplines = -1;
		NphiSplines = -1;
		NSplines = -1;
		NSplineParams = -1;
		nuMin = -1;
		muConst = -1;
	};
	~ProlateSplines(){ };

	ProlateSplines( int NnuSplinesIn, int NphiSplinesIn, double nuMinIn, double muConstIn, double * s );
	void createConstantProlateSplines( double f0, int NnuSplinesIn, int NphiSplinesIn, double nuMinIn, double muConstIn );
	void createRegularSplines( double * s, int NnuSplinesIn, int NphiSplinesIn, double nuMinIn, double muConstIn );
	void computeConstantSpline( double f0 );
	void setSplineSize( int NnuSplinesIn, int NphiSplinesIn, double nuMinIn, double muConstIn );

	void computeSplineParameters( double * s );

	// main output function
	double evaluateSplines( double nu0, double phi0 );

	// other outputs
	void evaluateSplineDerivatives( double nu0, double phi0, double &dfdnu0, double &dfdphi0 );
	double evaluateSplineCrossDerivative( double nu0, double phi0  );
	void evaluateSplineSecondDerivatives( double nu0, double phi0, double & d2fdnu2, double & d2fdnudphi, double & d2fdphi2 );

	double evaluateLaplaceSqInt();

};




#endif /* SOURCE_PROLATESPLINES_HPP_ */
