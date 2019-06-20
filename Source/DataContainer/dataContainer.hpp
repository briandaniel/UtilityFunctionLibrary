/*
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Copyright (C) 2017-2019 Brian D. Hong
 *
*/

/*
 * dataContainer.hpp
 *
 *  Created on: Dec 8, 2017
 *      Author: brian
 */

#ifndef DATACONTAINER_HPP_
#define DATACONTAINER_HPP_

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


#ifndef ROOT_ID
#define ROOT_ID 0
#endif



// class stores data up to 3 dimensions
class DataMember{

  public:
	DataMember(); // default constructor (empty)
	DataMember( double * x, int N1 ); // 1D array constructor
	DataMember( double ** x, int N1, int N2 ); // 2D array constructor
	DataMember( double *** x, int N1, int N2, int N3 ); // 3D array constructor

	~DataMember(); // default destructor clears the data
	void clearData(); // clears the data

	// ----------- read functions ------------ //
	bool checkTag( string tag ); // check if this member has the right tag
	string getTag(); // return the stag

	// pull single indexed data point from 1D array
	double read( int i );
	// pull single indexed data point from 2D array
	// --> these functions require updating for error catching
	double read( int i, int j );
	// pull single indexed data point from 3D array
	// --> these functions require updating for error catching
	double read( int i, int j, int k );

	// pull values to an array
	// --> these functions require updating for sizes and error catching
	void exportValues( double * x );
	// pull values to a 2D array
	// --> these functions require updating for sizes and error catching
	void exportValues( double ** x );
	// pull values to a 3D array
	// --> these functions require updating for sizes and error catching
	void exportValues( double *** x );

	// pull the size of the array
	void getSize( int & NdimOut, vector <int> & sizeOut );

	// ----------- write functions ---------- //
	void addDataPoint( double value ); // add a data point

	void replaceDataPoint( double value, int i ); // replace data point (1D)
	void replaceDataPoint( double value, int i, int j ); // replace data point (2D)
	void replaceDataPoint( double value, int i, int j, int k ); // replace data point (3D)

	void setTag( string tag ); // set the tag of the data member

	void setSize( int N ); // set data size 1D array
	void setSize( int N1, int N2 ); // set data size 2D array
	void setSize( int N1, int N2, int N3); // set data size 3D array

	void importArray( double * x, int N1 ); // import a 1D array
	void importArray( double ** x, int N1, int N2 ); // import a 2D array
	void importArray( double *** x, int N1, int N2, int N3 ); // import a 3D array

	void printData( int precision ); // print data


	void setStringData( string inputString ){ dataString = inputString; } // Set the data string
	string getStringData( ){ return dataString; } // get the string

	void setVerbosity( int verboseIn ){ verbose = verboseIn; } // set verbosity

  // ------------ Local storage ------------ //
  private:

	string dataTag;
	int Ndim;
	int Ndata;
	int verbose;

	// Vector of dimensions
	vector <int> size;

	// Vector that contains the data
	vector <double> data;

    // In case the data is a string
	string dataString;
	bool stringData; // indicates if text data was found

};




class DataContainer{

  public:

	DataContainer(); // default constructor sets size to zero
	~DataContainer(); // clear all data when the container goes out of scope

	void clearData(); // remove all data sets

	void addDataSet( double * x, string tag, int N1 ); // add a 1D data set
	void addDataSet( double ** x, string tag, int N1, int N2 ); // add a 2D data set
	void addDataSet( double *** x, string tag, int N1, int N2, int N3 ); // add a 3D data set

	void removeDataSet( string tag ); // remove data set

	void printDataSetList(); // print all the tags currently in the container

	void printDataSet( string tag, int precision ); // Print the requested data set

	void printDataSet( int dataSetIndex, int precision ); // Print the data set by index

	void printDataSetMPI( int dataSetIndex, int precision ); // Print the data set only on the root proc

	void readDataFile( string dataFileName, int verboseIn ); // Read a data file and store all the values

	/* wrapper for reading in parallel to avoid errors from accessing the input file
	 * each processor reads the file in serial
	 * --> this could be improved by reading only on root and the sending the data
	 *     to the other processors, but that would require some new code
	 */
	void readDataFileMPI( string dataFileName, int verboseIn );

	void newDataVectorFromFile( ifstream &fileID, string key, int N  ); // Pull 1D data set from file

	void readMatrixFromFile( ifstream &fileID, string key, int N1, int N2 ); // Pull 2D data set from file

	void getDataSize( string tag, int & Ndim, vector <int> & sizes ); // Get the data size

	int getDataLength( string tag ); // Gets the first dimension length of the data

	void readDataSet( string tag, double * x ); // Read 1D data set
	void readDataSet( string tag, double ** x ); // Read 2D data set
	void readDataSet( string tag, double *** x ); // Read 3D data set

	double readDataValue( string tag, int i ); // Read single value from 1D data set
	double readDataValue( string tag, int i, int j ); // Read single value from 2D data set
	double readDataValue( string tag, int i, int j, int k ); // Read single value from 3D data set

	int findDataSetIndex( string tag  ); // Find the index of a data set using its tag

	string getDataSetStringData( string tag  ); // Read single value from 1D data set

	void replaceDataSet( string tag, vector <double> & newDataVector ); // Replace the values in a (1D) data set
	void replaceDataSet( string tag, double newData ); // Replace single value data

	void copy1DDataSet( string tag, vector<double> &v ); // Pull 1D data set into a vector

	int getNDataSets(){ return NDataSets; };

	string getDataSetTag( int dataSetIndex ){ return data[dataSetIndex].getTag(); }

	string getDataSetStringData( int dataSetIndex ){ return data[dataSetIndex].getStringData(); }

	void renameDataSet( int dataSetIndex, string newName ){ data[dataSetIndex].setTag(newName); }

	void setVerbosity( int verboseIn ){ verbose = verboseIn; }

  private:

	int NDataSets;
	vector <DataMember> data;
	int verbose;

};




#endif /* DATACONTAINER_HPP_ */
