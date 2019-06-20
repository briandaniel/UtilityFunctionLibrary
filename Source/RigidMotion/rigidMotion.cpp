/*
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Copyright (C) 2017-2019 Brian D. Hong
 *
*/

/*
 * rigidMotion.cpp
 *
 *  Created on: Sep 30, 2017
 *      Author: brian
 *
 *	This code uses 6 parameters to determine the rigid body motion of a mathematical model.
 *
 *   1. This file uses the Z-X-Z Euler angles to adjust the rigid orientation of a model.
 *   2. The rigid translation is also applied in cartesian coordinates.
 *   3. This code was originally built as part of a 3D echo deformable image registration
 *
 */


#include "rigidMotion.hpp"


// rigid motion from model coordinates to echo coordinates
void RigidMotion::computeRigid( double xModel, double yModel, double zModel,
		           	   	   	    double & xEcho, double & yEcho, double & zEcho ){

	double xr, yr, zr;

	// Forward rotation
	xr = R11*xModel + R12*yModel + R13*zModel;
	yr = R21*xModel + R22*yModel + R23*zModel;
	zr = R31*xModel + R32*yModel + R33*zModel;

	// shift to echo coordinates
	xEcho = xr - sx;
	yEcho = yr - sy;
	zEcho = zr - sz;

}


// inverse rigid motion: from echo coordinates back to model coordinates
void RigidMotion::computeRigidInv( double xEcho, double yEcho, double zEcho,
								   double & xModel, double & yModel, double & zModel ){


	double xsinv, ysinv, zsinv;

	// shift back
	xsinv = xEcho + sx;
	ysinv = yEcho + sy;
	zsinv = zEcho + sz;

	// Inverse rotation
	xModel = Rinv11*xsinv + Rinv12*ysinv + Rinv13*zsinv;
	yModel = Rinv21*xsinv + Rinv22*ysinv + Rinv23*zsinv;
	zModel = Rinv31*xsinv + Rinv32*ysinv + Rinv33*zsinv;


}





void RigidMotion::computeRigidPrm( double * theta, double * shift){

	// forward shift from the vector
	sx = shift[0];
	sy = shift[1];
	sz = shift[2];


	// matrices
	double ** Rz1 = new double * [3];
	double **Rx = new double * [3];
	double ** Rz2 = new double * [3];
	double **R = new double * [3];
	double **Rinv = new double * [3];
	double ** temp = new double * [3];
	for(int i = 0; i < 3; i++)
	{
		Rz1[i] = new double [3];
		Rx[i] = new double [3];
		Rz2[i] = new double [3];
		R[i] = new double [3];
		Rinv[i] = new double [3];
		temp[i] = new double [3];

		for(int j = 0; j < 3; j++)
		{
			Rz1[i][j] = 0.0;
			Rx[i][j] = 0.0;
			Rz2[i][j] = 0.0;
			R[i][j] = 0.0;
			Rinv[i][j] = 0.0;
			temp[i][j] = 0.0;
		}
	}



	//------------- Forward rotation matrix ----------------//
	// Compute rotation
	Rz1[0][0] = cos(theta[0]);
	Rz1[0][1] = -sin(theta[0]);
	Rz1[1][0] = sin(theta[0]);
	Rz1[1][1] = cos(theta[0]);
	Rz1[2][2] = 1;

	Rx[0][0] = 1;
	Rx[1][1] = cos(theta[1]);
	Rx[1][2] = -sin(theta[1]);
	Rx[2][1] = sin(theta[1]);
	Rx[2][2] = cos(theta[1]);

	Rz2[0][0] = cos(theta[2]);
	Rz2[0][1] = -sin(theta[2]);
	Rz2[1][0] = sin(theta[2]);
	Rz2[1][1] = cos(theta[2]);
	Rz2[2][2] = 1;


	// multiply to get overall rotation matrix
	matrixMultiply( Rz2, Rx, temp, 3, 3, 3);
	matrixMultiply( temp, Rz1, R, 3, 3, 3);


	// Copy R to the doubles that are part of the class
	R11 = R[0][0];
	R12 = R[0][1];
	R13 = R[0][2];

	R21 = R[1][0];
	R22 = R[1][1];
	R23 = R[1][2];

	R31 = R[2][0];
	R32 = R[2][1];
	R33 = R[2][2];


	//------------- Inverse rotation matrix ----------------//
	// Compute rotation
	Rz1[0][0] = cos(-theta[0]);
	Rz1[0][1] = -sin(-theta[0]);
	Rz1[1][0] = sin(-theta[0]);
	Rz1[1][1] = cos(-theta[0]);
	Rz1[2][2] = 1;

	Rx[0][0] = 1;
	Rx[1][1] = cos(-theta[1]);
	Rx[1][2] = -sin(-theta[1]);
	Rx[2][1] = sin(-theta[1]);
	Rx[2][2] = cos(-theta[1]);

	Rz2[0][0] = cos(-theta[2]);
	Rz2[0][1] = -sin(-theta[2]);
	Rz2[1][0] = sin(-theta[2]);
	Rz2[1][1] = cos(-theta[2]);
	Rz2[2][2] = 1;


	// multiply to get overall rotation matrix
	matrixMultiply( Rz1, Rx, temp, 3, 3, 3);
	matrixMultiply( temp, Rz2, Rinv, 3, 3, 3);


	// Copy R to the doubles that are part of the class
	Rinv11 = Rinv[0][0];
	Rinv12 = Rinv[0][1];
	Rinv13 = Rinv[0][2];

	Rinv21 = Rinv[1][0];
	Rinv22 = Rinv[1][1];
	Rinv23 = Rinv[1][2];

	Rinv31 = Rinv[2][0];
	Rinv32 = Rinv[2][1];
	Rinv33 = Rinv[2][2];


	// clean up
	for(int k = 0; k < 3; k++)
	{
		delete [] Rz1[k];
		delete [] Rx[k];
		delete [] Rz2[k];
		delete [] R[k];
		delete [] Rinv[k];
		delete [] temp[k];
	}
	delete [] Rz1;
	delete [] Rx;
	delete [] Rz2;
	delete [] R;
	delete [] Rinv;
	delete [] temp;



}
