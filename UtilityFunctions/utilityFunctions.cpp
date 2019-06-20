/*
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Copyright (C) 2017-2019 Brian D. Hong
 *
*/

/*
 * utilityFunctions.cpp
 *
 *  Created on: Sep 9, 2017
 *      Author: brian
 *
 * This file describes some standard light-weight mathematical functionality for use in computational codes.
 *
 */

#include "utilityFunctions.hpp"


// Computes the modulus of two numbers
double mod(double x, double y)
{
	int n = (int) (x/y);
	double z = x - n*y;

	// In the case where the remainder is less than zero
	// Add one period to return to the positive interval
	if(z < 0)
	{
		z = z+y;
	}

	return z;
}

// compute mean of a vector
double mean( double * v, int N){

	double sum = 0.0;
	for(int k = 0; k < N; k++)
	{
		sum += v[k];
	}
	double meanVal = sum/N;

	return meanVal;

}

// Compute the value in a B spline wavelet
double bSplineWavelet( double x )
{
	double u, N1, N2, N3, N4;
	double B = 0;

	u = 4*x;

	if ( u>=0 && u<1 )
		B = 1.0/6.0*pow(u,3);
	else if ( u >= 1 && u < 2)
		B = -0.5*pow(u,3) + 2*pow(u,2) - 2*u + 2.0/3.0;
	else if ( u >= 2 && u < 3)
		B  = 0.5*pow(u,3) - 4*pow(u,2) + 10*u - 22.0/3.0;
	else if ( u >= 3 && u <= 4)
		B = -1.0/6.0*pow(u,3) + 2*pow(u,2) - 8*u + 32.0/3.0;

	return B;

}

// Compute periodic Bsplien wavelets
double periodicSplineWavelets( double phi0, int k )
{
	double twoPi = 2*PI;
	phi0 = mod(phi0, twoPi);
	double Bphi0 = 0;

	if( k == 0 )
		Bphi0 = bSplineWavelet( phi0/twoPi );
	else if (k == 1 )
		Bphi0 = bSplineWavelet( ( phi0 - PI/2 )/ twoPi ) +  bSplineWavelet( ( phi0 + 3*PI/2 )/ twoPi );
	else if (k == 2 )
		Bphi0 = bSplineWavelet( ( phi0 - PI )/ twoPi ) +  bSplineWavelet( ( phi0 + PI )/ twoPi );
	else if (k == 3 )
		Bphi0 = bSplineWavelet( ( phi0 - 3*PI/2 )/ twoPi ) +  bSplineWavelet( ( phi0 + PI/2 )/ twoPi );

	return Bphi0;
}

// Computes y = A*x where A has size Ni x Nj, so that x is necessarily of size Nj and y has size Ni
void matrixVectorMultiply( double** A, double *x, int Ni, int Nj, double * y  )
{

	for (int i = 0; i < Ni; i++)
	{
		y[i] = 0;
		for (int j = 0; j < Nj; j++)
		{
			y[i] = y[i] + A[i][j]*x[j];
		}
	}

}

// Computes the index of the minimum of an array x with size N
int arrayMin( double* x, int N)
{
	int iMin = 0;
	for (int i = 0; i < N; i++)
	{
		if ( x[i] < x[iMin] )
		{
			iMin = i;
		}

	}
	return iMin;
}

// Converts prolate spheroid coordinates to cartesian coordinates
void prolate2xyz( double mu, double nu, double phi, double a, double &x, double &y, double &z)
{

	x = a*sinh(mu)*sin(nu)*cos(phi);
	y = a*sinh(mu)*sin(nu)*sin(phi);
	z = a*cosh(mu)*cos(nu);

}

// Converts cartesian to prolate spheroid coordinates
void xyz2prolate( double x, double y, double z, double a, double &mu, double &nu, double &phi )
{

	double dist1, dist2, coef, sigma, tau;

	dist1 = sqrt( pow(x,2) + pow(y,2) + pow(z+a,2) );
	dist2 = sqrt( pow(x,2) + pow(y,2) + pow(z-a,2) );
	coef = 1.0/(2*a);

	sigma = coef*(dist1 + dist2);
	tau = coef*( dist1 - dist2 );

	// Round off error can cause tau outside [-1,1] (even though this can't actually occur with valid numbers for x, y, z, and a)
	if( tau > 1 )
	{
		tau = 1;
	}else if( tau < -1 )
	{
		tau = -1;
	}

	mu = acosh(sigma);
	nu = acos(tau);

	phi = atan2( y, x );

	// put phi in [0,2pi)
	// although technically it really shouldnt matter
	if ( phi < 0 )
	{
		phi = phi + 2*PI;
	}
}

// Compute the matrix determinant of A, where A is 3x3 only
// This is a poor description/implementation and should be replaced
double computeMatrixDeterminant( double ** A )
{
	double detA = A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1]
				- A[0][2]*A[1][1]*A[2][0] - A[0][1]*A[1][0]*A[2][2] - A[0][0]*A[1][2]*A[2][1];

	return detA;
}

// Compute the inverse of a 3x3 matrix
void computeMatrixInverse33( double ** A, double ** Ainv )
{

	double detA = computeMatrixDeterminant( A );
	Ainv[0][0] = ( A[1][1]*A[2][2] - A[1][2]*A[2][1] )/ detA;
	Ainv[0][1] = ( A[0][2]*A[2][1] - A[0][1]*A[2][2] )/ detA;
	Ainv[0][2] = ( A[0][1]*A[1][2] - A[0][2]*A[1][1] )/ detA;

	Ainv[1][0] = ( A[1][2]*A[2][0] - A[1][0]*A[2][2] )/ detA;
	Ainv[1][1] = ( A[0][0]*A[2][2] - A[0][2]*A[2][0] )/ detA;
	Ainv[1][2] = ( A[0][2]*A[1][0] - A[0][0]*A[1][2] )/ detA;

	Ainv[2][0] = ( A[1][0]*A[2][1] - A[1][1]*A[2][0] )/ detA;
	Ainv[2][1] = ( A[0][1]*A[2][0] - A[0][0]*A[2][1] )/ detA;
	Ainv[2][2] = ( A[0][0]*A[1][1] - A[0][1]*A[1][0] )/ detA;

}

// multiplies matrices so that C = A*B, with A = NixNk, B = NkxNj
void matrixMultiply( double ** A, double ** B, double ** C, int Ni, int Nj, int Nk)
{
	for (int i = 0; i < Ni; i++)
	{
		for (int j = 0; j < Nj; j++)
		{
			C[i][j] = 0;
			for (int k = 0; k < Nk; k++)
			{
				C[i][j] = C[i][j] + A[i][k]*B[k][j];
			}
		}
	}
}

// Legacy, these types of implementations are poor, use vectors.
void instantiate33Matrix( double ** A )
{
	if (A == NULL)
	{
		A = new double *[3];
		for(int k = 0; k < 3; k++)
		{
			A[k] = new double[3];
		}

		for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j<3; j++)
			{
				A[i][j] = 0;
			}
		}
	}
}

// Legacy, these types of implementations are poor, use vectors.
void free33Matrix( double ** A )
{
	if (A != NULL)
	{
		for(int i = 0; i < 3; i++)
		{
			if(A[i] != NULL)
			{
				delete[] A[i];
			}
		}
		delete [] A;
		A = NULL;
	}
}

// Legacy, these types of implementations are poor, use vectors.
void freeVector( double * v )
{

	if ( v )
	{
		delete [] v;
	}
	v = NULL;
}

// Compute the k-norm of a double array with size N
double vectorNorm ( double * v, int N, double k )
{
	double normv, sumvsq;

	sumvsq = 0.0;
	for(int i = 0; i < N; i++)
		sumvsq = sumvsq + pow(v[i],k);

	normv = pow( sumvsq, 1.0/k );

	return normv;
}

// Normalize a double array with length N according to the Euclidean norm
void normalizeVector2norm( double * v, int N )
{
	double normv;
	double power = 2.0;

	// 2 norm
	normv = vectorNorm ( v, N, power );
	if(normv == 0)
	{
		normv = 1;
	}

	for(int k = 0; k < N; k++)
		v[k] = v[k] / normv;
}

// computes cross product of length 3 vectors
void crossProduct3( double * a, double * b, double * acrossb )
{
	acrossb[0] = a[1]*b[2] - a[2]*b[1];
	acrossb[1] = a[2]*b[0] - a[0]*b[2];
	acrossb[2] = a[0]*b[1] - a[1]*b[0];
}

// Dot product of vectors of length N
double dotProd( double* a, double* b, int N )
{
	double sum = 0;
	for(int k = 0; k<N; k++)
		sum = sum + a[k]*b[k];

	return sum;
}

// computes matrix transpose
void matrixTranspose( double ** A, double ** AT, int Ni, int Nj)
{
	for(int i = 0; i < Ni; i++)
	{
		for(int j = 0; j < Nj; j++)
		{
			AT[j][i] = A[i][j];
		}
	}
}

// Print 3x3 matrix A
void print33Matrix( double** A )
{
	for(int i = 0; i < 3; i++)
	{

		if (i == 1)
		{
			cout << "A = |";
		}
		else
		{
			cout << "    |";
		}
		for(int j = 0; j < 3; j++)
		{
			std::cout << std::setw(16) << A[i][j] << " ";
		}
		cout << "|" << endl;
	}
}

// Integral coefficients for an evenly spaced grid. Note that this
// still requires the delta x multiplication, since it is not
// included in this coefficient
double computeIntegralCoef ( int i, int N )
{
	double c = 0.0;

	if ( N >= 9 ) // extended simpsons
	{

		if( i == 0 || i == N-1 ){ c = 17.0 / 48.0; }
		else if( i == 1 || i == N-2 ){ c = 59.0 / 48.0; }
		else if( i == 2 || i == N-3 ){ c = 43.0 / 48.0; }
		else if( i == 3 || i == N-4 ){ c = 49.0 / 48.0; }
		else{ c = 1.0; }

	}else{ // just use trapezoid

		if( i == 0 || i == N-1 ){ c = 1.0/2.0; }
		else{ c = 1.0; }
	}

	//c = 1;
	return c;
}


// Integral coefficients for a periodic function
// assumes that the function does not include the point phi = 2*pi
// uses a trapezoid rule
double computePeriodicIntegralCoef ( int i, int N )
{
	double c;

	c = 1.0;

	return c;
}

// Computes the double dot of the two matrices (sum of pairwise multiplication)
double matrixDoubleDot( double ** A, double ** B, int Ni, int Nj )
{

	double sum = 0.0;

	for(int i = 0; i < Ni; i++)
	{
		for(int j = 0; j < Nj; j++)
		{
			sum = sum + A[i][j]*B[i][j];
		}
	}

	return sum;
}

// Print the values in a 2D array to the console
void print2DArray( double** A, int Ni, int Nj, int precision )
{

	cout << endl;
	for(int i = 0; i < Ni; i++)
	{

		if (i == floor(Ni/2.0))
		{
			cout << "A = |";
		}
		else
		{
			cout << "    |";
		}
		for(int j = 0; j < Nj; j++)
		{
			std::cout << std::setprecision( precision) << std::setw( precision + 3) << A[i][j] << " ";
		}
		cout << "|" << endl;
	}
	cout << endl;
}

// Print the values in a 1D array to the console
void print1DArray( double* v, int Ni, int precision )
{

	cout << endl;
	for(int i = 0; i < Ni; i++)
	{

		if (i == floor(Ni/2.0))
		{
			cout << "A = |";
		}
		else
		{
			cout << "    |";
		}

		std::cout << std::setprecision( precision) << std::setw( precision + 10) << v[i] << " ";
		cout << "|" << endl;
	}
	cout << endl;
}

// Print the values in a 1D array to the console
void print1DArrayLine( double* v, int Ni, int precision, string name )
{

	cout << std::left;
	cout << std::setw(20)  << name <<  " = [ ";
	for(int i = 0; i < Ni-1; i++)
	{

		std::cout << std::setprecision( precision) << std::setw( precision + 3) << v[i] << ", ";

	}
	std::cout << std::setprecision( precision) << std::setw( precision + 3) << v[Ni-1] << " ]" << endl;

}

// Compute LU factorization
void lu(double **A, double **L, double **U, int n)
{
	for (int k = 0; k < n; k++){
		for (int j = 0; j < n; j++){
			L[k][j] = 0.0;
			U[k][j] = 0.0;
		}
	}

	for (int k = 0; k < n; k++){
		L[k][k] = 1.0;

		for (int i = k + 1; i < n; i++){

			L[i][k] = A[i][k] / A[k][k];
			for (int j = k + 1; j < n; j++)
			{
				A[i][j] = A[i][j] - L[i][k] * A[k][j];
			}
		}
		for (int j = k; j < n; j++)
		{
			U[k][j] = A[k][j];
		}
	}

}

// Compute forward substitution (part of LU matrix solution)
void forwardSubstitution(double**L, double*y, double*b, int n)
{
	// Solves lower triangular system Ly = b
	double alpha;
	for (int i = 0; i < n; i++)
	{
		alpha = b[i];
		for (int j = 0; j < i ; j++){
			alpha = alpha - L[i][j] * y[j];
		}
		y[i] = alpha / L[i][i];
	}

}

// Compute backward substitution (part of LU matrix solution)
void backSubstitution(double**U, double*x, double*y, int n)
{
	// Solves upper triangular system Ux = y
	for (int i = n - 1; i >= 0; i--)
	{
		x[i] = y[i];
		for (int j = i + 1; j < n; j++){
			x[i] = x[i] - U[i][j] * x[j];
		}
		x[i] = x[i] / U[i][i];
	}

}

// solves the system Ax = b by lu factorization
void luSolve( double ** A, double * b, int n, double *x )
{
	double * v = new double [n];
	double ** L = new double * [n];
	double ** U = new double * [n];
	for(int k = 0; k < n; k++)
	{
		L[k] = new double [n];
		U[k] = new double [n];
	}

	// compute LU factorization
	lu(A, L, U, n);

	// solve Lv = b for v, where v = Ux
	forwardSubstitution(L, v, b, n);

	// solve Ux = v for x
	backSubstitution(U, x, v, n);

	// clean up
	delete [] v;

	for(int k = 0; k < n; k++)
	{
		delete [] L[k];
		delete [] U[k];
	}

	delete [] L;
	delete [] U;

}

// x is assumed to be linearly increasing
// but equal spacing is not assumed
// so this could be made much faster by assuming equal spacing
// right now it just looks through all the x values
// should make another one called "linearInterpEqual" or something
double linearInterp( double * x, double * y, int N, double xx )
{
	double dx1, dx2;
	double yy;
	int success = 0;

	for(int k = 0; k < N-1; k++ )
	{
		if( xx >= x[k] && xx <= x[k+1] )
		{
			dx1 = xx-x[k];
			dx2 = x[k+1] - xx;
			yy = ( dx2*y[k] + dx1*y[k+1] )/( dx1 + dx2 );

			success = 1;
			break;
		}
	}

	if (success == 0)
	{
		cout << "Value " << xx << " for linear interpolation is not inside the given range of x values: "<< x[0] << " to " << x[N-1] << endl;
		cout << "So the value is automatically set the output to 0." << endl;
		yy = 0;
	}

	return yy;
}

// Compute equal spacing with N points and place in vector x
// (Same as Matlab linspace)
void linspace( double a, double b, int N, double*x )
{
	double dx = (b - a)/(N-1);
	for (int k = 0; k < N; k++)
	{
		x[k] = a + dx*k;
	}

}

// Rotation matrix such that v_cart = Q*v_prolate
void computeProlateToCartesianRotationMatrix( double mu, double nu, double phi, double **  Q )
{

    double shmu = sinh(mu);
    double chmu = cosh(mu);
    double snu = sin(nu);
    double cnu = cos(nu);
    double sphi = sin(phi);
    double cphi = cos(phi);

    double b = 1/sqrt( pow(shmu,2) + pow(snu,2) );

    Q[0][0] = b*chmu*snu*cphi;
    Q[0][1] = b*shmu*cnu*cphi;
	Q[0][2] = -sphi;

	Q[1][0] = b*chmu*snu*sphi;
	Q[1][1] = b*shmu*cnu*sphi;
	Q[1][2] = cphi;

	Q[2][0] = b*shmu*cnu;
	Q[2][1] = -b*chmu*snu;
	Q[2][2] = 0;

}


// Jacobian of the transformation (x,y,z)->(mu,nu,phi)
// with elements T_ij = dxi/dthetaj
// where x = [x,y,z], theta = [mu,nu,phi]
void computeProlateToCartesianTransformationJacobian( double mu, double nu, double phi, double a, double ** T )
{

    double shmu = sinh(mu);
    double chmu = cosh(mu);
    double snu = sin(nu);
    double cnu = cos(nu);
    double sphi = sin(phi);
    double cphi = cos(phi);

    T[0][0] = a*chmu*snu*cphi;
    T[0][1] = a*shmu*cnu*cphi;
	T[0][2] = -a*shmu*snu*sphi;

	T[1][0] = a*chmu*snu*sphi;
	T[1][1] = a*shmu*cnu*sphi;
	T[1][2] = a*shmu*snu*cphi;

	T[2][0] = a*shmu*cnu;
	T[2][1] = -a*chmu*snu;
	T[2][2] = 0;

}

// Find minimum value of a double array
void vectorMin (double * v, int Nv, double & minv, int & indexMin)
{
	indexMin = 0;
	minv = v[0]*2;

	for( int k = 0; k < Nv; k++ )
	{
		if(v[k] < minv)
		{
			minv = v[k];
			indexMin = k;
		}
	}

}

// Find maximum value of a double array
void vectorMax (double * v, int Nv, double & maxv, int & indexMax)
{
	indexMax = 0;
	maxv = v[0]/2;

	for( int k = 0; k < Nv; k++ )
	{
		if(v[k] > maxv)
		{
			maxv = v[k];
			indexMax= k;
		}
	}

}

// Generates a random number between 0 and 1 using the time stamp
double timeRand()
{
	int aNum = rand();
	double num = (double)aNum / RAND_MAX;

	return num;
}

// Returns the sign of a number
int sign( double x )
{
	double value = 0;
	if(x == 0)
		value = 0;
	else
		value = ( x/fabs(x) );

	int valueInt = (int) value;
	return valueInt;
}

// Generates a random number between 0 and 1 using hardware
double hardRand()
{
	std::mt19937_64 randGen(std::random_device{}());
	uniform_real_distribution<double> unif;
	return unif(randGen);
}

// Determine the process ID
int getProcID()
{
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);
	return procID;
}

// use the in if statement to execute only if on ROOT
bool isRootProc()
{

	if( getProcID() == 0 )
	{
		return true;
	}

	return false;
}


