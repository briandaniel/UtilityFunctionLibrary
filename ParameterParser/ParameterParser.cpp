/*
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Copyright (C) 2017-2019 Brian D. Hong
 *
*/

/*
 * ParameterParser.cpp
 *
 *  Created on: Nov 20, 2018
 *      Author: brian
 *
 * The goal of this is to incorporate multiple simple parameter files into one
 * easy to use container called "prmData"
 *
 * 	1. Each parameter is tagged as "<paramPrefix>_<prmName>"
 * 	2. Parameters can only be read if they are 1D or values by this function,
 * 	   however, the data container can still be read for higher dimensional sets
 *  3. The param prefix is the name given by the parameter file names defined by the
 * 	   file that was used to give the parameter file names
 * 	4. This is still a very early version, so although it serves its purpose,
 * 	   it would require substantial work for broad usefulness
 *
 */

#include "ParameterParser.hpp"


void readPrmVector( string paramName, string paramPrefix, DataContainer & prmData, vector <double> & prmOut )
{
	vector <int> sizes;
	int Ndim;
	string tag = paramPrefix  + "_" + paramName;

	// return data information
	prmData.getDataSize( tag, Ndim, sizes );
	int N = 0;
	if( sizes.size() > 0)
		N = sizes[0];

	prmOut.resize(N);

	for(int i = 0; i < N; i++)
	{
		prmOut[i] = prmData.readDataValue( tag, i );
	}
}

void readPrmVector( string paramName, string paramPrefix, DataContainer & prmData, vector <int> & prmOut )
{
	vector <int> sizes;
	int Ndim;
	string tag = paramPrefix  + "_" + paramName;

	// return data information
	prmData.getDataSize( tag, Ndim, sizes );
	int N = 0;
	if( sizes.size() > 0)
		N = sizes[0];

	prmOut.resize(N);

	for(int i = 0; i < N; i++)
	{
		prmOut[i] = (int) prmData.readDataValue( tag, i );
	}
}


void readPrmVector( string paramName, string paramPrefix, DataContainer & prmData, double * prmOut, int N )
{
	for(int i = 0; i < N; i++)
	{
		prmOut[i] = prmData.readDataValue(  paramPrefix  + "_" + paramName, i );
	}
}

double readPrmVectorValue( int idx, string paramName, string paramPrefix, DataContainer & prmData )
{
	double value = prmData.readDataValue(  paramPrefix  + "_" + paramName, idx ); // read single value

	return value;
}

double readPrmValue( string paramName, string paramPrefix, DataContainer & prmData )
{
	double value = prmData.readDataValue(  paramPrefix  + "_" + paramName, 0 ); // read single value

	return value;
}


string readPrmString( string paramName, string paramPrefix, DataContainer & prmData )
{
	string out = prmData.getDataSetStringData(  paramPrefix  + "_" + paramName ); // read single value

	return out;
}


void readParseParameterFiles( string paramFileNamesLocation, DataContainer & prmData, int verbose )
{

	DataContainer fileLocations;
	fileLocations.readDataFileMPI(paramFileNamesLocation, verbose);
	int NprmFiles = fileLocations.getNDataSets();

	int N1 = 0;
	int N2 = 0;

	for(int k = 0; k < NprmFiles; k++)
	{
		string fileTag = fileLocations.getDataSetTag(k);

		string fileName = fileLocations.getDataSetStringData(k);

		prmData.readDataFileMPI( fileName, verbose );

		N2 = prmData.getNDataSets() ;

		for(int i = N1; i < N2; i++)
		{
			prmData.renameDataSet(i, fileTag + '_' + prmData.getDataSetTag(i) );
		}

		N1 = N2;

	}

}


