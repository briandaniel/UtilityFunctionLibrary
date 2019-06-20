/*
 * mathematicalOperators.hpp
 *
 *  Created on: Jul 4, 2018
 *      Author: brian
 */

#ifndef MATHEMATICALOPERATORS_HPP_
#define MATHEMATICALOPERATORS_HPP_


// Standard constant definitions
#define PI 3.141592653589793

// standard includes
#include <math.h>
#include <complex>
#include <iostream>

// local functions
void computeCartesianTensorDivergence( double *** dTdx, double * divT );
void computeProlateTensorDivergence(  double mu, double nu, double phi, double a, double ** S, double *** dSdi, double * divS_prolate );
void prolateSpheroidalChristoffelSymbols( double mu, double nu, double phi, double a, double *** gamma, double ** D, double * h );



#endif /* MATHEMATICALOPERATORS_HPP_ */
