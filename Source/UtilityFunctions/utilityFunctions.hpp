/*
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Copyright (C) 2017-2019 Brian D. Hong
 *
*/

/*
 * utilityFunctions.h
 *
 *  Created on: Sep 9, 2017
 *      Author: brian
 */

#ifndef UTILITYFUNCTIONS_HPP_
#define UTILITYFUNCTIONS_HPP_

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <ctime>
#include <random>
#include <mpi.h>

#ifndef PI
#define PI 3.141592653589793
#endif

using namespace std;

// Include all other utility function headers
#include "utilityFunctionsVector.hpp"


// Computes the modulus of two numbers
double mod (double x, double y);

// compute mean of a vector
double mean( double * v, int N);

// Compute the value in a B spline wavelet
double bSplineWavelet( double x );

// Compute periodic Bsplien wavelets
double periodicSplineWavelets( double phi0, int k );

// Computes y = A*x where A has size Ni x Nj, so that x is necessarily of size Nj and y has size Ni
void matrixVectorMultiply( double** A, double *x, int Ni, int Nj, double * y  );

// Computes the index of the minimum of an array x with size N
int arrayMin( double* x, int N);

// Converts prolate spheroid coordinates to cartesian coordinates
void prolate2xyz( double mu, double nu, double phi, double a, double &x, double &y, double &z);

// Converts cartesian to prolate spheroid coordinates
void xyz2prolate( double x, double y, double z, double a, double &mu, double &nu, double &phi );

// Compute the matrix determinant of A, where A is 3x3 only
// This is a poor description/implementation and should be replaced
double computeMatrixDeterminant( double ** A );

// Compute the inverse of a 3x3 matrix
void computeMatrixInverse33( double ** A, double ** Ainv );

// multiplies matrices so that C = A*B, with A = NixNk, B = NkxNj
void matrixMultiply( double ** A, double ** B, double ** C, int Ni, int Nj, int Nk);

// Legacy, these types of implementations are poor, use vectors.
void instantiate33Matrix( double ** A );

// Legacy, these types of implementations are poor, use vectors.
void free33Matrix( double ** A );

// Legacy, these types of implementations are poor, use vectors.
void freeVector( double * v );

// Compute the k-norm of a double array with size N
double vectorNorm ( double * v, int N, double k );

// Normalize a double array with length N according to the Euclidean norm
void normalizeVector2norm( double * v, int N );

// computes cross product of length 3 vectors
void crossProduct3( double * a, double * b, double * acrossb );

// Dot product of vectors of length N
double dotProd( double* a, double* b, int N );

// computes matrix transpose
void matrixTranspose( double ** A, double ** AT, int Ni, int Nj);

// Print 3x3 matrix A
void print33Matrix( double** A );

// Integral coefficients for an evenly spaced grid. Note that this
// still requires the delta x multiplication, since it is not
// included in this coefficient
double computeIntegralCoef ( int i, int N );

// Integral coefficients for a periodic function
// assumes that the function does not include the point phi = 2*pi
// uses a trapezoid rule
double computePeriodicIntegralCoef ( int i, int N );

// Computes the double dot of the two matrices (sum of pairwise multiplication)
double matrixDoubleDot( double ** A, double ** B, int Ni, int Nj );

// Print the values in a 2D array to the console
void print2DArray( double** A, int Ni, int Nj, int precision );

// Print the values in a 1D array to the console
void print1DArray( double* v, int Ni, int precision );

// Print the values in a 1D array to the console
void print1DArrayLine( double* v, int Ni, int precision, string name );

// Compute LU factorization
void lu(double **A, double **L, double **U, int n);

// Compute forward substitution (part of LU matrix solution)
void forwardSubstitution(double**L, double*y, double*b, int n);

// Compute backward substitution (part of LU matrix solution)
void backSubstitution(double**U, double*x, double*y, int n);

// solves the system Ax = b by lu factorization
void luSolve( double ** A, double * b, int n, double *x );

// x is assumed to be linearly increasing
// but equal spacing is not assumed
// so this could be made much faster by assuming equal spacing
// right now it just looks through all the x values
// should make another one called "linearInterpEqual" or something
double linearInterp( double * x, double * y, int N, double xx );

// Compute equal spacing with N points and place in vector x
// (Same as Matlab linspace)
void linspace( double a, double b, int N, double*x );

// Rotation matrix such that v_cart = Q*v_prolate
void computeProlateToCartesianRotationMatrix( double mu, double nu, double phi, double **  Q );

// Jacobian of the transformation (x,y,z)->(mu,nu,phi)
// with elements T_ij = dxi/dthetaj
// where x = [x,y,z], theta = [mu,nu,phi]
void computeProlateToCartesianTransformationJacobian( double mu, double nu, double phi, double a, double ** T );

// Find minimum value of a double array
void vectorMin (double * v, int Nv, double & minv, int & indexMin);

// Find maximum value of a double array
void vectorMax (double * v, int Nv, double & maxv, int & indexMax);

// Generates a random number between 0 and 1 using the time stamp
double timeRand();

// Returns the sign of a number
int sign( double x );

// Generates a random number between 0 and 1 using hardware
double hardRand();

// Determine the process ID
int getProcID();

// use the in if statement to execute only if on ROOT
bool isRootProc();


// Legacy, these types of implementations are poor, use vectors.
#define SAFE_DELETE( pPtr ) { if ( pPtr ) { delete [] pPtr; pPtr = NULL; } }
#define SAFE_DELETE_2D( pPtr, N1 ){ if ( pPtr ) { for(int i = 0; i < N1; i++) { delete [] pPtr[i]; } delete [] pPtr; pPtr = NULL; } }

#endif /* UTILITYFUNCTIONS_HPP_ */
















