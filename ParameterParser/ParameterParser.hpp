/*
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Copyright (C) 2017-2019 Brian D. Hong
 *
*/

/*
 * ParameterParser.hpp
 *
 *  Created on: Nov 20, 2018
 *      Author: brian
 */

#ifndef PARAMETERPARSER_PARAMETERPARSER_HPP_
#define PARAMETERPARSER_PARAMETERPARSER_HPP_

// Standard includes
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

// local dependencies
#include "../ImportExport/importExport.hpp"
#include "../DataContainer/dataContainer.hpp"

#ifndef ROOT_ID
#define ROOT_ID 0
#endif

// local functions
void readPrmVector( string paramName, string paramPrefix, DataContainer & prmData, vector <double> & prmOut );
void readPrmVector( string paramName, string paramPrefix, DataContainer & prmData, vector <int> & prmOut );
void readPrmVector( string paramName, string paramPrefix, DataContainer & prmData, double * prmOut, int N );
double readPrmVectorValue( int idx, string paramName, string paramPrefix, DataContainer & prmData );
double readPrmValue( string paramName, string paramPrefix, DataContainer & prmData );
string readPrmString( string paramName, string paramPrefix, DataContainer & prmData );

void readParseParameterFiles( string paramFileNamesLocation, DataContainer & prmData, int verbose );



#endif /* PARAMETERPARSER_PARAMETERPARSER_HPP_ */
