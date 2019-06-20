/*
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Copyright (C) 2017-2019 Brian D. Hong
 *
*/

/*
 * importExport.h
 *
 *  Created on: Sep 23, 2017
 *      Author: brian
 */

#ifndef SOURCE_IMPORTEXPORT_HPP_
#define SOURCE_IMPORTEXPORT_HPP_


// Includes
#include <iostream>
#include <fstream>
#include <string>
#include <inttypes.h>
#include <sstream>
#include <cstring>
#include <iomanip>
#include <mpi.h>
#include <math.h>
#include <vector>

// local includes
#include "../UtilityFunctions/utilityFunctions.hpp"

using namespace std;

#ifndef ROOT_ID
#define ROOT_ID 0
#endif

// Looks through the parameter file and finds keys
void findKeys( ifstream &fileID, vector<string> & keys, int verbose );

// determines whether or not the keys are 1D or 2D and their size
void findKeySize( ifstream &fileID, string key, int & dim, int & N1, int & N2, bool & stringFound, int verbose );

// This finds a single value from a line that contains the key "string"
// Doesnt read lines that contain "//"
double findValue( ifstream &fileID, string key, int verbose  );

// This finds a vector of values from a line that contains the key "string"
void findArray( ifstream &fileID, string key, int N, double * values, int verbose );

// This is a poor re-implementation of the above findArray
// --> really this should just call findArray of doubles then convert them to integers
void findArray( ifstream &fileID, string key, int N, int * values, int verbose  );

// This finds a vector of values from a line that contains the key "string"
void find2DArray( ifstream &fileID, string key, int Nj, int Nk, double ** values, int verbose );

// This finds a string that follows the key
string findString( ifstream &fileID, string key, int verbose  );

// Prints a matlab formatted value (note: this annoying version also tells matlab to print the value when it runs)
void printMatlabVariable( ofstream &fileID, string key, double value  );

// Prints a matlab formatted value (note: this version doesn't print when loaded into matlab)
void printMatlabVariableSimple( ofstream &fileID, string key, double value );

// Prints a matlab formatted 1D array (note: this annoying version also tells matlab to print the value when it runs)
void printMatlabArray( ofstream &fileID, string key, double * values, int N );

// Prints a matlab formatted array ((note: this version doesn't print when loaded into matlab)
void printMatlabArraySimple( ofstream &fileID, string key, double * values, int N );

// Prints a matlab formatted 2D array (note: this annoying version also tells matlab to print the value when it runs)
void printMatlab2DArray( ofstream &fileID, string key, double ** values, int N1, int N2 );

// Prints a matlab formatted 2D array (note: this version doesn't print when loaded into matlab)
void printMatlab2DArraySimple( ofstream &fileID, string key, double ** values, int N1, int N2  );

// Prints a matlab formatted 3D array (note: this version doesn't print when loaded into matlab)
void printMatlab3DArraySimple( ofstream &fileID, string key, double *** values, int N1, int N2, int N3  );

// Prints a matlab formatted 3D array (note: this version doesn't print when loaded into matlab)
void printMatlab4DArraySimple( ofstream &fileID, string key, double **** values, int N1, int N2, int N3, int N4 );

// Prints a matlab formatted 1D array from a 1D vector of doubles (doesnt print)
void printMatlab1DArray( ofstream &fileID, string key, vector<double> & values);


#endif /* SOURCE_IMPORTEXPORT_HPP_ */
