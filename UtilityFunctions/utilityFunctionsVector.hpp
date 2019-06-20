/*
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Copyright (C) 2017-2019 Brian D. Hong
 *
*/

/*
 * utilityFunctionsVector.hpp
 *
 *  Created on: Nov 6, 2018
 *      Author: brian
 */

#ifndef UTILITYFUNCTIONS_UTILITYFUNCTIONSVECTOR_HPP_
#define UTILITYFUNCTIONS_UTILITYFUNCTIONSVECTOR_HPP_

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <ctime>

#ifndef PI
#define PI 3.141592653589793
#endif

using namespace std;

// compute sum of a integer vector
int sum( vector<int> & v );

// compute sum of a double vector
double sum( vector<double> & v );

// compute mean of a vector
double mean( vector<double> & v );

// compute mean of absolute value of entries of a vector
double absMean( vector<double> & v );

// Computes y = A*x where A has size Ni x Nj, so that x is necessarily of size Nj and y has size Ni
void matrixVectorMultiply(vector<vector<double>> &A, vector<double> &x , vector<double> &y  );

// multiplies matrices so that C = A*B, with A = NixNk, B = NkxNj
void matrixMultiply( vector<vector<double>> &A, vector<vector<double>> & B, vector<vector<double>> &C );

// compute the Euclidean norm of a vector
double vector2Norm ( vector<double> v );

// Dot product of vectors of length N
double dotProd( vector<double> & a, vector <double> & b );

// Compute matrix transpose
void matrixTranspose( vector<vector<double>> & A, vector<vector<double>>  & AT);

// Compute LU factorization
void lu( vector<vector<double>> & A, vector<vector<double>> & L, vector<vector<double>> & U);

// Compute forward substitution (part of LU matrix solution)
void forwardSubstitution( vector<vector<double>> & L,  vector<double> &y,  vector<double> &b);

// Compute backward substitution (part of LU matrix solution)
void backSubstitution(vector<vector<double>> & U, vector<double> &x, vector<double> &y );

// solves the system Ax = b by lu factorization
void luSolve( vector<vector<double>> & A, vector<double> & b, vector<double> & x );

// computes inverse of A(nxn) using LU factorization
void matrixInverse( vector<vector<double>> & A,  vector<vector<double>> & Ainv );

// Find minimum value of a vector
void vectorMin ( vector <double> & v, int Nv, double & minv, int & indexMin);

// find maximum value of a vector
void vectorMax ( vector <double> & v, int Nv, double & maxv, int & indexMax);

// Find minimum value of a vector (without indexing)
double vectorMin ( vector <double> & v );

// Find max of a vector (without index)
double vectorMax ( vector <double> & v );

// Print the values in a <double> vector in the console, standard use:
// cout << "<vector name>"; print1DVector(v);
void print1DVector( vector <double> & v );

// Print the values in a <int> vector in the console, standard use:
// cout << "<vector name>"; print1DVector(v);
void print1DVector( vector <int> & v );

// Print the values in a <bool> vector in the console, standard use:
// cout << "<vector name>"; print1DVector(v);
void print1DVector( vector <bool> & v );

// Print the values in a 2D vector to the console
void print2DVector ( vector< vector<double> > & A );

// Print the values in a 3D vector to the console
void print3DVector ( vector< vector < vector<double> > > & A );

// Compute equal spacing with N point and place in vector x
// (Same as Matlab linspace)
void linspace( double a, double b, int N, vector <double> & x );

// A set to identity (size of A unchanged)
void setIdentity( vector< vector<double> > & A );

// Linear interpolation on an uneven grid for evaluating a single value
/*
 * x is assumed to be linearly increasing but equal spacing is not assumed
 * so this could be made much faster by assuming equal spacing right now
 * it just looks through all the x values should make another one called
 * "linearInterpEqual" or something
 *
 */
double linearInterp( vector <double> & x, vector <double> & y, double xx );

// Linear interpolation on an uneven grid for evaluating multiple values
void linearInterpMulti( vector <double> & x, vector <double> & y, vector <double> & xx, vector <double> & yy );

// assumes that the two inputs already have the same size
void copyVector( vector<double> & v_original, vector<double> & v_copy );

// Compute the matrix determinant of A, where A is 3x3 only
// This is a poor description/implementation and should be replaced
double computeMatrixDeterminant( vector<vector<double>> & A );

// Compute the inverse of a 3x3 matrix
void computeMatrixInverse33( vector<vector<double>> & A, vector<vector<double>> & Ainv );

// Rotation matrix such that v_cart = Q*v_prolate
void computeProlateToCartesianRotationMatrix( double mu, double nu, double phi, vector< vector<double> > &  Q );

// Computes the double dot of the two matrices (sum of pairwise multiplication)
double matrixDoubleDot( vector<vector<double>> & A, vector<vector<double>> & B );

// computes cross product of length 3 vectors
void crossProduct3( vector<double> & a, vector<double> & b, vector<double> & acrossb );

// Find the maximum of the absolute difference between all elements of v
double maxElementAbsDiff( vector <double> & v);






#endif /* UTILITYFUNCTIONS_UTILITYFUNCTIONSVECTOR_HPP_ */
