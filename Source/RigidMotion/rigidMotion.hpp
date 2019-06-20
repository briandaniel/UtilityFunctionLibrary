/*
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Copyright (C) 2017-2019 Brian D. Hong
 *
*/

/*
 *
 * rigidMotion.hpp
 *
 *  Created on: Sep 30, 2017
 *      Author: brian
 *
 *
 */


#ifndef SOURCE_RIGIDMOTION_HPP_
#define SOURCE_RIGIDMOTION_HPP_


#include <math.h>
#include "../UtilityFunctions/utilityFunctions.hpp"


// define a rigid motion class which precomputes the rotation matrix
class RigidMotion {

  public:

	// For simplicity just use doubles because we only have two matrices
	double R11, R12, R13, R21, R22, R23, R31, R32, R33;
	double Rinv11, Rinv12, Rinv13, Rinv21, Rinv22, Rinv23, Rinv31, Rinv32, Rinv33;

	// shift values
	double sx, sy, sz;

	// precompute the parameters for rigid motion
	void computeRigidPrm( double * theta, double * shift);

	// rigid motion from model coordinates to echo coordinates
	void computeRigid( double xModel, double yModel, double zModel, double & xEcho, double & yEcho, double & zEcho );

	// inverse rigid motion: from echo coordinates back to model coordinates
	void computeRigidInv( double xEcho, double yEcho, double zEcho, double & xModel, double & yModel, double & zModel );

	// Empty constructor, must call computeRigidPrm before use
	RigidMotion(){};


};

#endif /* SOURCE_RIGIDMOTION_HPP_ */
