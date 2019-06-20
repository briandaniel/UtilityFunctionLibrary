/*
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Copyright (C) 2017-2019 Brian D. Hong
 *
*/

/*
 * mathematicalOperators.cpp
 *
 *  Created on: Jul 4, 2018
 *      Author: brian
 *
 * This file describes mathematical operators that are too specific too be considered "utility"
 *
 */


#include "mathematicalOperators.hpp"

// dTdx has size 3x3x3 where the indices define
// row x column x derivative
void computeCartesianTensorDivergence( double *** dTdx, double * divT )
{
	for(int i = 0; i < 3; i++)
	{
		divT[i] = dTdx[i][0][0] + dTdx[i][1][1] + dTdx[i][2][2];
	}
}

// Prolate tensor divergence operator. Ew.
void computeProlateTensorDivergence(  double mu, double nu, double phi, double a, double ** S, double *** dSdi, double * divS_prolate )
{
	double *** gamma = new double ** [3];
	double ** D = new double * [3];
	double * h = new double [3];

	for(int i = 0; i < 3; i++)
	{
		D[i] = new double [3];
		gamma[i] = new double *[3];
		for(int j = 0; j < 3; j++)
		{
			gamma[i][j] = new double [3];
		}
	}

	// compute the coefficients for the prolate spheroidal coordinate system
	prolateSpheroidalChristoffelSymbols( mu, nu, phi, a, gamma, D, h );

	// compute the general divergence of a tensor (in particular, here, using prolate spheroid coordinates,
	//  but this could be reapplied to other curvilinear coordinate systems)
	for(int i = 0; i < 3; i++)
	{
		double Ai = 0;
		double Bi = 0;
		double Ci = 0;

		for(int j = 0; j < 3; j++)
		{
			Ai = Ai + 1/h[j]*dSdi[i][j][j] + S[i][j]*D[i][j];

			for(int m = 0; m < 3; m++)
			{
				Bi = Bi + 1/h[m]*gamma[m][j][j]*S[i][m];

				Ci = Ci - h[m]/(h[i]*h[j])*gamma[i][j][m]*S[m][j];
			}
		}

		divS_prolate[i] = Ai + Bi + Ci;
	}


	// clean up
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			delete [] gamma[i][j];
		}
		delete [] D[i];
		delete [] gamma[i];
	}

	delete [] h;
	delete [] D;
	delete [] gamma;

}

void prolateSpheroidalChristoffelSymbols( double mu, double nu, double phi, double a, double *** gamma, double ** D, double * h )
{

	double snu, cnu, shmu, chmu, sumsq;
	snu = sin(nu);
	cnu = cos(nu);
	shmu = sinh(mu);
	chmu = cosh(mu);
	sumsq = pow(shmu,2) + pow(snu,2);

	h[0] = a*sqrt( sumsq );
	h[1] = h[0];
	h[2] = a*shmu*snu;

	// mostly zeros
	for(int i = 0; i < 3; i++ )
	{
		for(int j = 0; j < 3; j++ )
		{
			D[i][j] = 0.0;
			for(int k = 0; k < 3; k++ )
			{
				gamma[i][j][k] = 0.0;
			}
		}
	}

	// nonzero entries (except symmetric entries)
	gamma[0][0][0] = shmu*chmu/sumsq;
	gamma[1][1][0] = -shmu*chmu/sumsq;
	gamma[2][2][0] = -shmu*chmu*pow(snu,2)/sumsq;
	gamma[0][1][0] = snu*cnu/sumsq;
	gamma[1][1][1] = snu*cnu/sumsq;
	gamma[0][0][1] = -snu*cnu/sumsq;
	gamma[2][2][1] = -pow(shmu,2)*snu*cnu/sumsq;
	gamma[1][0][1] = shmu*chmu/sumsq;
	gamma[2][0][2] = chmu/shmu;
	gamma[2][1][2] = cnu/snu;

	// symmetry entries
	gamma[1][0][0] = gamma[0][1][0];
	gamma[0][1][1] = gamma[1][0][1];
	gamma[0][2][2] = gamma[2][0][2];
	gamma[1][2][2] = gamma[2][1][2];

	// Two nonzero D entries
	D[2][0] = 1/h[2]*( chmu*pow(snu,3)/sqrt( pow(sumsq,3) ) );
	D[2][1] = 1/h[2]*( cnu*pow(shmu,3)/sqrt( pow(sumsq,3) ) );

}


























