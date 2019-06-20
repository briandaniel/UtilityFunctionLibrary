/*
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Copyright (C) 2017-2019 Brian D. Hong
 *
*/

/*
 * utilityFunctionsVector.cpp
 *
 *  Created on: Nov 6, 2018
 *      Author: brian
 *
 * This file describes some standard light-weight mathematical functionality for use in computational codes
 * but using the std::vector
 *
 */


#include "utilityFunctionsVector.hpp"


// compute sum of a integer vector
int sum( vector<int> & v )
{

	int sum = 0.0;
	for(int k = 0; k < v.size(); k++)
	{
		sum += v[k];
	}
	return sum;

}


// compute sum of a double vector
double sum( vector<double> & v )
{

	double sum = 0.0;
	for(int k = 0; k < v.size(); k++)
	{
		sum += v[k];
	}
	return sum;

}


// compute mean of a vector
double mean( vector<double> & v )
{

	double sum = 0.0;
	for(int k = 0; k < v.size(); k++)
	{
		sum += v[k];
	}
	double meanVal = sum/v.size();

	return meanVal;

}

// compute mean of absolute value of entries of a vector
double absMean( vector<double> & v )
{

	double sum = 0.0;
	for(int k = 0; k < v.size(); k++)
	{
		sum += fabs(v[k]);
	}
	double meanVal = sum/v.size();

	return meanVal;

}



// Computes y = A*x where A has size Ni x Nj, so that x is necessarily of size Nj and y has size Ni
void matrixVectorMultiply(vector<vector<double>> &A, vector<double> &x , vector<double> &y  )
{

	for (int i = 0; i < A.size(); i++)
	{
		y[i] = 0;
		for (int j = 0; j < A[i].size(); j++)
		{
			y[i] = y[i] + A[i][j]*x[j];
		}
	}

}


// multiplies matrices so that C = A*B, with A = NixNk, B = NkxNj
void matrixMultiply( vector<vector<double>>& A, vector<vector<double>> &B, vector<vector<double>> &C )
{

	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < B[0].size(); j++)
		{
			C[i][j] = 0;
			for (int k = 0; k < B.size(); k++)
			{
				C[i][j] = C[i][j] + A[i][k]*B[k][j];

			}

		}
	}

}


// compute the Euclidean norm of a vector
double vector2Norm ( vector<double> v )
{
	double normv, sumvsq;

	sumvsq = 0.0;
	for(int i = 0; i < v.size(); i++)
		sumvsq = sumvsq + pow(v[i],2);

	normv = pow( sumvsq, 1.0/2.0);

	return normv;
}


// Dot product of vectors of length N
double dotProd( vector<double> & a, vector <double> & b )
{

	double sum = 0;
	for(int k = 0; k<a.size(); k++)
		sum = sum + a[k]*b[k];

	return sum;

}


// Compute matrix transpose
void matrixTranspose( vector<vector<double>> & A, vector<vector<double>>  & AT)
{
	for(int i = 0; i < A.size(); i++)
	{
		for(int j = 0; j < A[i].size(); j++)
		{
			AT[j][i] = A[i][j];
		}
	}

}


// Compute LU factorization
void lu( vector<vector<double>> & A, vector<vector<double>> & L, vector<vector<double>> & U)
{
	int n = A.size();
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
void forwardSubstitution( vector<vector<double>> & L,  vector<double> &y,  vector<double> &b)
{
	// Solves lower triangular system Ly = b
	int n = b.size();
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
void backSubstitution(vector<vector<double>> & U, vector<double> &x, vector<double> &y )
{
	int n = y.size();

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
void luSolve( vector<vector<double>> & A, vector<double> & b, vector<double> & x )
{
	int n = b.size();

	vector<double> y(n,0);
	vector<double> v(n,0);

	vector<vector<double> > L(n, vector<double>(n));
	vector<vector<double> > U(n, vector<double>(n));

	// compute LU factorization
	lu(A, L, U);

	// solve Lv = b for v, where v = Ux
	forwardSubstitution(L, v, b);

	// solve Ux = v for x
	backSubstitution(U, x, v );


}


// computes inverse of A(nxn) using LU factorization
void matrixInverse( vector<vector<double>> & A,  vector<vector<double>> & Ainv )
{
	int n = A.size();
	vector<double> ei(n,0);
	vector<double> xi(n,0);
	vector<double> v(n,0);

	vector<vector<double> > L(n, vector<double>(n));
	vector<vector<double> > U(n, vector<double>(n));

	// compute LU factorization
	lu( A, L, U);


	// Solve for ainv one column at a time
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
			ei[j] = 0;
		ei[i] = 1;

		forwardSubstitution(L, v, ei);

		backSubstitution(U, xi, v );

		for(int j = 0; j < n; j++)
		{
			Ainv[i][j] = xi[j];
		}


	}

}


// Find minimum value of a vector
void vectorMin ( vector <double> & v, int Nv, double & minv, int & indexMin)
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


// find maximum value of a vector
void vectorMax ( vector <double> & v, int Nv, double & maxv, int & indexMax)
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


// Find minimum value of a vector (without indexing)
double vectorMin ( vector <double> & v )
{
	double minv = v[0]*2;

	for( int k = 0; k < v.size(); k++ )
	{
		if(v[k] < minv)
		{
			minv = v[k];
		}
	}

	return minv;
}


// Find max of a vector (without index)
double vectorMax ( vector <double> & v )
{
	double maxv = v[0]/2;

	for( int k = 0; k < v.size(); k++ )
	{
		if(v[k] > maxv)
		{
			maxv = v[k];
		}
	}

	return maxv;
}


// Print the values in a <double> vector in the console, standard use:
// cout << "<vector name>"; print1DVector(v);
void print1DVector( vector <double> & v )
{
	int precision = 6;
	cout << " [ ";
	for(int i = 0; i < v.size()-1; i++)
	{

		std::cout << std::setprecision( precision) << std::setw( precision + 6) << v[i] << ", ";

	}
	std::cout << std::setprecision( precision) << std::setw( precision + 6) << v[v.size()-1] << " ]" << endl;

}


// Print the values in a <int> vector in the console, standard use:
// cout << "<vector name>"; print1DVector(v);
void print1DVector( vector <int> & v )
{
	int precision = 6;
	cout << " [ ";
	for(int i = 0; i < v.size()-1; i++)
	{

		std::cout << std::setprecision( precision) << std::setw( precision + 6) << v[i] << ", ";

	}
	std::cout << std::setprecision( precision) << std::setw( precision + 6) << v[v.size()-1] << " ]" << endl;

}


// Print the values in a <bool> vector in the console, standard use:
// cout << "<vector name>"; print1DVector(v);
void print1DVector( vector <bool> & v )
{
	int precision = 3;
	cout << " [ ";
	for(int i = 0; i < v.size()-1; i++)
	{

		std::cout << std::setprecision( precision) << std::setw( precision + 6) << v[i] << ", ";

	}
	std::cout << std::setprecision( precision) << std::setw( precision + 6) << v[v.size()-1] << " ]" << endl;

}


// Print the values in a 2D vector to the console
void print2DVector ( vector< vector<double> > & A )
{
	int precision = 6;

	cout << endl;
	for(int i = 0; i < A.size(); i++)
	{
		cout << "  |";
		for(int j = 0; j < A[i].size(); j++)
		{
			std::cout << std::setprecision( precision) << std::setw( precision + 2) << A[i][j] << "     ";
		}
		cout << "|" << endl;
	}
	cout << endl;
}


// Print the values in a 3D vector to the console
void print3DVector ( vector< vector < vector<double> > > & A )
{
	int precision = 6;

	cout << endl;

	for(int i = 0; i < A.size(); i++)
	{

		for(int j = 0; j < A[i].size(); j++)
		{
			if(j == floor(A[i].size()/2))
			{
				cout <<  " A [" << i << ",:,:] = " << " |";
			}
			else
			{
				cout << std::setw( 15) << " |";
			}
			for(int k = 0; k < A[i][j].size(); k++)
			{
				std::cout << std::setprecision( precision) << std::setw( precision + 2) << A[i][j][k] << "     ";
			}
			cout << "|" << endl;
		}
		cout << endl << endl;
	}
	cout << endl;
}


// Compute equal spacing with N points and place in vector x
// (Same as Matlab linspace)
void linspace( double a, double b, int N, vector <double> & x )
{
	x.resize(N);
	double dx = (b - a)/(N-1);
	for (int k = 0; k < N; k++)
	{
		x[k] = a + dx*k;
	}

}


// A set to identity (size of A unchanged)
void setIdentity( vector< vector<double> > & A )
{
	int precision = 6;

	for(int i = 0; i < A.size(); i++)
	{
		for(int j = 0; j < A[i].size(); j++)
		{
			if(i == j)
			{
				A[i][j] = 1;
			}
			else
			{
				A[i][j] = 0.0;
			}
		}
	}
}


// Linear interpolation on an uneven grid for evaluating a single value
/*
 * x is assumed to be linearly increasing but equal spacing is not assumed
 * so this could be made much faster by assuming equal spacing right now
 * it just looks through all the x values should make another one called
 * "linearInterpEqual" or something
 *
 */
double linearInterp( vector <double> & x, vector <double> & y, double xx )
{
	int N = x.size();
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


// Linear interpolation on an uneven grid for evaluating multiple values
void linearInterpMulti( vector <double> & x, vector <double> & y, vector <double> & xx, vector <double> & yy )
{
	yy.resize(xx.size());
	for(int i = 0; i < xx.size(); i++)
	{
		yy[i] = linearInterp( x, y, xx[i]);
	}

}


// assumes that the two inputs already have the same size
void copyVector( vector<double> & v_original, vector<double> & v_copy )
{
	for(int i = 0; i < v_original.size(); i++ )
	{
		v_copy[i] = v_original[i];
	}
}


// Compute the matrix determinant of A, where A is 3x3 only
// This is a poor description/implementation and should be replaced
double computeMatrixDeterminant( vector<vector<double>> & A )
{
	double detA = A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1]
				- A[0][2]*A[1][1]*A[2][0] - A[0][1]*A[1][0]*A[2][2] - A[0][0]*A[1][2]*A[2][1];

	return detA;
}


// Compute the inverse of a 3x3 matrix
void computeMatrixInverse33( vector<vector<double>> & A, vector<vector<double>> & Ainv )
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


// Rotation matrix such that v_cart = Q*v_prolate
void computeProlateToCartesianRotationMatrix( double mu, double nu, double phi, vector< vector<double> > &  Q )
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


// Computes the double dot of the two matrices (sum of pairwise multiplication)
double matrixDoubleDot( vector<vector<double>> & A, vector<vector<double>> & B )
{

	double sum = 0.0;

	for(int i = 0; i < A.size(); i++)
	{
		for(int j = 0; j < A[i].size(); j++)
		{
			sum = sum + A[i][j]*B[i][j];
		}
	}

	return sum;
}


// computes cross product of length 3 vectors
void crossProduct3( vector<double> & a, vector<double> & b, vector<double> & acrossb )
{
	acrossb[0] = a[1]*b[2] - a[2]*b[1];
	acrossb[1] = a[2]*b[0] - a[0]*b[2];
	acrossb[2] = a[0]*b[1] - a[1]*b[0];
}


// Find the maximum of the absolute difference between all elements of v
double maxElementAbsDiff( vector <double> & v)
{

	double maxElAbsDiff = 0;
	for(int i = 0; i < v.size(); i++)
	{
		for(int j = 0; j < v.size(); j++)
		{
			if( fabs( v[i] - v[j] ) > maxElAbsDiff )
			{
				maxElAbsDiff = fabs( v[i] - v[j] );
			}
		}

	}

	return maxElAbsDiff;

}
























