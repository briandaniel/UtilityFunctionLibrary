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
 * dataContainer.cpp
 *
 *  Created on: Dec 8, 2017
 *      Author: brian
 *
 * These classes are meant to store information from input files or other data sources.
 * 	1. The implementation avoids memory issues and stores all data in each member in a single vector
 * 	2. Each data member is called based on the tag given to it
 * 	3. The data container is just a vector of DataMembers and functions to find and pull desired values
 *	4. The functions are generally pretty self-explanatory
 *
 */

#include "dataContainer.hpp"


//--------------------------------------------------------------------------///
// DataMember class functions

// Default constructor
DataMember::DataMember(){
	setSize(0);
	verbose = 0;
}

// 1D array constructor
DataMember::DataMember( double * x, int N1 ){

	importArray( x, N1 );

}

// 2D array constructor
DataMember::DataMember( double ** x, int N1, int N2 ){

	importArray( x, N1, N2 );

}

// 3D array constructor
DataMember::DataMember( double *** x, int N1, int N2, int N3 ){

	importArray( x, N1, N2, N3 );

}

// default destructor clears the data
DataMember::~DataMember(){
	clearData();
}

// clears the data
void DataMember::clearData(){
	Ndim = 0;
	Ndata = 0;
	size.clear();
	data.clear();
}

// check if this member has the right tag
bool DataMember::checkTag( string tag ){

	// if the compare == 0, then they are the same
	if( !dataTag.compare(tag) )
	{
		return true;
	}

	return false;
}

// return the string tag
string DataMember::getTag(){
	return dataTag;
}

// pull single indexed data point from 1D array
double DataMember::read( int i )
{
	int index = i;
	double value = 0;

	if( data.size() > i )
		value = data[index];
	else
	{
		// If value is outside the range, throw an error
		if(getProcID() == 0)
			cout << "DataContainer error = Data value does not exist!! Returning 0 instead." << endl;
	}
	return value;
}

// pull single indexed data point from 2D array
// --> these functions require updating for error catching
double DataMember::read( int i, int j )
{

	int index = i + size[0]*j;
	double value = data[index];
	return value;

}

// pull single indexed data point from 3D array
// --> these functions require updating for error catching
double DataMember::read( int i, int j, int k )
{
	int index = i + size[0]*j + size[0]*size[1]*k;
	double value = data[index];
	return value;
}

// pull values to an array
// --> these functions require updating for sizes and error catching
void DataMember::exportValues( double * x )
{
	for(int i = 0; i < size[0]; i ++ )
	{
		x[i] = data[i];
	}

}

// pull values to a 2D array
// --> these functions require updating for sizes and error catching
void DataMember::exportValues( double ** x )
{
	for(int i = 0; i < size[0]; i ++ )
	{
		for(int j = 0; j < size[1]; j ++ )
		{
			int index = i + size[0]*j;
			x[i][j] = data[index];
		}
	}

}

// pull values to a 3D array
// --> these functions require updating for sizes and error catching
void DataMember::exportValues( double *** x )
{
	for( int i = 0; i < size[0]; i ++ )
	{
		for( int j = 0; j < size[1]; j ++ )
		{
 			for( int k = 0; k < size[2]; k++ )
 			{
 				int index = i + size[0]*j + size[0]*size[1]*k;
 				x[i][j][k] = data[index];
 			}
		}
	}
}


// pull the size of the array
void DataMember::getSize( int & NdimOut, vector <int> & sizeOut ){

	sizeOut.clear();

	NdimOut = Ndim;
	for(int dim = 0; dim < Ndim; dim++)
	{
		sizeOut.push_back( size[dim] );
	}
}


// ----------- write functions ---------- //
// add a data point
void DataMember::addDataPoint( double value )
{
	data.push_back(value);
}

// replace data point (1D)
void DataMember::replaceDataPoint( double value, int i ) // replace value in 1D array
{
	int index = i;
	data[index] = value;
}

// replace data point (2D)
void DataMember::replaceDataPoint( double value, int i, int j ) // replace value in 2D array
{
	int index = i + size[0]*j;
	data[index] = value;

}

// replace data point (3D)
void DataMember::replaceDataPoint( double value, int i, int j, int k ) // replace value in 3D array
{
	int index = i + size[0]*j + size[0]*size[1]*k;
	data[index] = value;
}

// set the tag of the data member
void DataMember::setTag( string tag )
{
	dataTag = tag;
}

// set data size 1D array
void DataMember::setSize( int N )
{
	clearData();

	Ndim = 1;
	Ndata = N;
	size.push_back(N);
	data.resize(Ndata);

}

// set data size 2D array
void DataMember::setSize( int N1, int N2 )
{
	clearData();

	Ndim = 2;
	Ndata = N1*N2;
	size.push_back(N1);
	size.push_back(N2);
	data.resize(Ndata);

}

// set data size 3D array
void DataMember::setSize( int N1, int N2, int N3)
{
	clearData();

	Ndim = 3;
	Ndata = N1*N2*N3;
	size.push_back(N1);
	size.push_back(N2);
	size.push_back(N3);
	data.resize(Ndata);
}

// import 1D array
void DataMember::importArray( double * x, int N1 )
{
	setSize(N1);

	for(int i = 0; i < size[0]; i ++ )
	{
		data[i] = x[i];
	}

}

// import 2D array
void DataMember::importArray( double ** x, int N1, int N2 )
{
	setSize(N1,N2);

	for(int i = 0; i < size[0]; i ++ )
	{
		for(int j = 0; j < size[1]; j ++ )
		{
			int index = i + size[0]*j;
			data[index] = x[i][j];
		}
	}

}

// import 3D array
void DataMember::importArray( double *** x, int N1, int N2, int N3 )
{

	setSize(N1,N2,N3);

	for( int i = 0; i < size[0]; i ++ )
	{
		for( int j = 0; j < size[1]; j ++ )
		{
 			for( int k = 0; k < size[2]; k++ )
 			{

 				int index = i + size[0]*j + size[0]*size[1]*k;
 				data[index] = x[i][j][k];
 			}
		}
	}

}


// print data
void DataMember::printData( int precision ){

	int success = 0;


	if(Ndim == 0 || data.size() == 0){

		cout << "No data to print!" << endl;
		success = 1;

	} else if (Ndim == 1){

		cout <<  "Vector " << dataTag << " with length " << size[0] << ": "<< endl;
		cout << dataTag << " = [ ";
		for( int i = 0; i < size[0]-1; i ++ )
		{
			std::cout << std::setprecision( precision) << std::setw( precision + 3) << data[i] << ", ";

		}
		std::cout << std::setprecision( precision) << std::setw( precision + 3) << data[size[0]-1] << " ]" << endl;
		cout << endl;

		success = 1;
	} else if (Ndim == 2){

		cout <<  "Matrix " << dataTag << " with size " << size[0] << " x " << size[1] << ": "<< endl;
		for(int i = 0; i < size[0]; i++)
		{
			string temp = dataTag + " = | ";

			if (i == floor(size[0]/2.0))
			{
				cout <<  temp;
			}
			else
			{
				cout << std::setw( temp.size() ) << " | ";
			}
			for(int j = 0; j < size[1]; j++)
			{
				int index = i + size[0]*j;


				std::cout << std::setprecision( precision) << std::setw( precision + 3) << data[index]  << " ";
			}
			cout << "|" << endl;
		}
		cout << endl;

		success = 1;
	} else if (Ndim == 3){

		cout <<  "3D vector " << dataTag << " with size " << size[0] << " x " << size[1] << " x " << size[2] <<  ": "<< endl;
		for(int i = 0; i < size[0]; i++)
		{

			for(int j = 0; j < size[1]; j++)
			{
				string temp = dataTag + "_i="+ to_string(i) + " = | ";
				if (j == floor(size[1]/2.0))
				{
					cout << temp;
				}
				else
				{
					cout << std::setw( temp.size() ) << " | ";
				}
				for(int k = 0; k < size[2]; k++)
				{
	 				int index = i + size[0]*j + size[0]*size[1]*k;
					std::cout << std::setprecision( precision) << std::setw( precision + 4) << data[index]  << " ";
				}
				cout << "|" << endl;
			}
			cout << endl;


		}
		success = 1;
	}


	if(success == 0)
	{
		cout << "Data was not printed.." << endl;
	}



}


//--------------------------------------------------------------------------///
// DataContainer class functions

// default constructor
DataContainer::DataContainer()
{
	NDataSets = 0;
	verbose = 0;
}

// clear all data when the container goes out of scope
DataContainer::~DataContainer()
{
	clearData();
}

// remove all data sets
void DataContainer::clearData(){

	// clear the data containers
	for(int k = 0; k < NDataSets; k++)
	{
		data[k].clearData();
	}

	// now clear the vector of containers
	NDataSets = 0;
	data.clear();
}

// add a 1D data set
void DataContainer::addDataSet( double * x, string tag, int N1 ){

	// temporary data member
	DataMember dmt(x,N1);
	dmt.setTag(tag);

	data.push_back(dmt);
	NDataSets = data.size();
}

// add a 2D data set
void DataContainer::addDataSet( double ** x, string tag, int N1, int N2 ){

	// temporary data member
	DataMember dmt(x,N1,N2);
	dmt.setTag(tag);

	data.push_back(dmt);
	NDataSets = data.size();

}

// add a 3D data set
void DataContainer::addDataSet( double *** x, string tag, int N1, int N2, int N3 ){

	// temporary data member
	DataMember dmt(x,N1,N2,N3);
	dmt.setTag(tag);

	data.push_back(dmt);
	NDataSets = data.size();
}


// remove data set
void DataContainer::removeDataSet( string tag ){

	int success = 0;
	for(int i = 0; i < data.size(); i++)
	{
		if ( data[i].checkTag(tag) )
		{
			data.erase( data.begin() + i );
			success = 1;
		}
	}
	if(success == 0)
		cout << "Data set " << tag << "cannot be erased because it does not exist!" << endl;

}

// print all the tags currently in the container
void DataContainer::printDataSetList()
{

	if(getProcID() == ROOT_ID)
	{
		cout << " Data container has: " << endl;
		for(int i = 0; i < data.size(); i++)
		{
			cout << "  set " << i << ": " << data[i].getTag() << endl;
		}
		cout << endl;
	}

}

// Print the requested data set
void DataContainer::printDataSet( string tag, int precision ){

	int success = 0;

	for(int i = 0; i < data.size(); i++)
	{
		if ( data[i].checkTag(tag) )
		{
			data[i].printData( precision );
			success = 1;
		}
	}
	if(success == 0)
		cout << "Data set " << tag << " cannot be printed because it does not exist!" << endl;

}

// Print the data set by index
void DataContainer::printDataSet( int dataSetIndex, int precision )
{
	int success = 0;

	if( dataSetIndex < data.size() )
	{
		data[dataSetIndex].printData( precision );
		success = 1;
	}
	if(success == 0)
		cout << "Data set with index " << dataSetIndex << " cannot be printed because it does not exist!" << endl;

}

// Print the data set only on the root proc
void DataContainer::printDataSetMPI( int dataSetIndex, int precision )
{
	// Get the number of processes
	int Nprocs, procID;
	MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	if(procID == ROOT_ID)
	{
		int success = 0;

		if( dataSetIndex < data.size() )
		{
			data[dataSetIndex].printData( precision );
			success = 1;
		}
		if(success == 0)
			cout << "Data set with index " << dataSetIndex << " cannot be printed because it does not exist!" << endl;
	}
}

// Read a data file and store all the values
void DataContainer::readDataFile( string dataFileName, int verboseIn ){

	verbose = verboseIn;

	int Nprocs, procID;
	MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);


	ifstream fileID( dataFileName.c_str() , ios::in  );
	if( !fileID.is_open() )
	{
		cout << "File \"" << dataFileName  << "\" requested does not exist..." << endl;
	}
	else
	{

		// Look for keys in the data file
		vector<string> keys;
		findKeys( fileID, keys, verbose );

		int dim, N1, N2;
		bool stringFound;

		for(int i = 0; i < keys.size(); i++ )
		{
			// Determine the size of the data labeled by the i^{th} key
			// Also determine if it has a string
			findKeySize( fileID, keys[i], dim, N1, N2, stringFound, verbose  );

			// Now read the data file
			if(dim == 1)
			{
				newDataVectorFromFile( fileID, keys[i], N2  );
			}
			else if(dim == 2)
			{
				readMatrixFromFile( fileID, keys[i], N1, N2);
			}
			else
			{
				cout << "While reading data from file " << dataFileName << " key " << keys[i] << " failed!" << endl;
			}

			// Save the string if one was found
			if( stringFound )
			{
				int idx = data.size()-1;
				data[idx].setStringData( findString( fileID, keys[i], verbose  ) );
			}
		}
	}
	fileID.close();

}

/* wrapper for reading in parallel to avoid errors from accessing the input file
 * each processor reads the file in serial
 * --> this could be improved by reading only on root and the sending the data
 *     to the other processors, but that would require some new code
 */
void DataContainer::readDataFileMPI( string dataFileName, int verboseIn ){

	if(getProcID() == ROOT_ID)
	{
		cout << "Data container reading \"" << dataFileName << "\" (MPI safe). " << endl;
	}

	// Get the number of processes
	int Nprocs, procID;
	MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	int file_free = 0;
	MPI_Status status;

	// master gets permission to read parameter file first, others wait in line
	if ( procID == ROOT_ID )
		file_free = 1;
	else
		MPI_Recv(&file_free, 1, MPI_INT, procID-1, 1, MPI_COMM_WORLD, &status);


	// this process calls the routine to read and process input file if it is free
	if (file_free == 1)
	{
		readDataFile( dataFileName, verboseIn );
	}

	// give read file permission to the next process
	if (procID != Nprocs-1)
		MPI_Send(&file_free, 1, MPI_INT, procID+1, 1, MPI_COMM_WORLD);

}

// Pull 1D data set from file
void DataContainer::newDataVectorFromFile( ifstream &fileID, string key, int N  ){

	double * x = new double [N];

	findArray( fileID, key, N, x, verbose );

	addDataSet( x, key, N );

	// clean up
	delete [] x;
}

// Pull 2D data set from file
void DataContainer::readMatrixFromFile( ifstream &fileID, string key, int N1, int N2 ){

	double ** M = new double * [N1];
	for(int i = 0; i < N1; i++ ){
		M[i] = new double [N2];
	}

	find2DArray( fileID, key, N1, N2, M , verbose );
	addDataSet( M, key, N1, N2 );

	// clean up
	for(int i = 0; i < N1; i++ ){ delete [] M[i]; }
	delete [] M;

}

// Get the data size
void DataContainer::getDataSize( string tag, int & Ndim, vector <int> & sizes )
{
	int success = 0;
	for(int i = 0; i < data.size(); i++)
	{
		if ( data[i].checkTag(tag) )
		{
			data[i].getSize( Ndim, sizes );
			success = 1;
			break;
		}
	}

	if(success == 0)
		cout << "Data set " << tag << " cannot be accessed because it does not exist!" << endl;
}

// Gets the first dimension length of the data
int DataContainer::getDataLength( string tag )
{
	int Ndim;
	vector <int> sizes;
	getDataSize( tag, Ndim, sizes );
	if(Ndim > 1)
	{
		cout << "Warning! Length of array with more than one dimension accessed." << endl;
	}
	int length = sizes[0];

	return length;
}

// Read 1D data set
void DataContainer::readDataSet( string tag, double * x )
{
	int success = 0;
	for(int i = 0; i < data.size(); i++)
	{
		if ( data[i].checkTag(tag) )
		{
			data[i].exportValues(x);
			success = 1;
			break;
		}
	}

	if(success == 0)
		cout << "Data set " << tag << " cannot be accessed because it does not exist!" << endl;
}

// Read 2D data set
void DataContainer::readDataSet( string tag, double ** x ){
	int success = 0;
	for(int i = 0; i < data.size(); i++)
	{
		if ( data[i].checkTag(tag) )
		{
			data[i].exportValues(x);
			success = 1;
			break;
		}
	}

	if(success == 0)
		cout << "Data set " << tag << " cannot be accessed because it does not exist!" << endl;
}


// Read 3D data set
void DataContainer::readDataSet( string tag, double *** x ){
	int success = 0;
	for(int i = 0; i < data.size(); i++)
	{
		if ( data[i].checkTag(tag) )
		{
			data[i].exportValues(x);
			success = 1;
			break;
		}
	}

	if(success == 0)
		cout << "Data set " << tag << " cannot be accessed because it does not exist!" << endl;
}

// Read single value from 1D data set
double DataContainer::readDataValue( string tag, int i ){

	double value = 0;

	int success = 0;
	for(int n = 0; n < data.size(); n++)
	{
		if ( data[n].checkTag(tag) )
		{
			value = data[n].read(i);
			success = 1;
			break;
		}
	}

	if(success == 0)
		cout << "Data value " << tag << " cannot be accessed because it does not exist!" << endl;

	return value;
}

// Read single value from 2D data set
double DataContainer::readDataValue( string tag, int i, int j ){

	double value = 0;

	int success = 0;
	for(int n = 0; n < data.size(); n++)
	{
		if ( data[n].checkTag(tag) )
		{
			value = data[n].read(i,j);
			success = 1;
			break;
		}
	}

	if(success == 0)
		cout << "Data value " << tag << " cannot be accessed because it does not exist!" << endl;

	return value;

}

// Read single value from 3D data set
double DataContainer::readDataValue( string tag, int i, int j, int k ){

	double value = 0;

	int success = 0;
	for(int n = 0; n < data.size(); n++)
	{
		if ( data[n].checkTag(tag) )
		{
			value = data[n].read(i,j,k);
			success = 1;
			break;
		}
	}

	if(success == 0)
		cout << "Data set " << tag << " cannot be accessed because it does not exist!" << endl;

	return value;

}

// Find the index of a data set using its tag
int DataContainer::findDataSetIndex( string tag  ){

	int idx = -1;
	int success = 0;
	for(int n = 0; n < data.size(); n++)
	{
		if ( data[n].checkTag(tag) )
		{
			idx = n;
			success = 1;
			break;
		}
	}

	if(success == 0)
		cout << "Data index " << tag << " cannot be accessed because it does not exist!" << endl;

	return idx;
}

// Read single value from 1D data set
string DataContainer::getDataSetStringData( string tag  )
{
	int idx = findDataSetIndex( tag  );
	return data[idx].getStringData();
}

// Replace the values in a (1D) data set
void DataContainer::replaceDataSet( string tag, vector <double> & newDataVector )
{

	int idx = -1;
	int success = 0;
	for(int n = 0; n < data.size(); n++)
	{
		if ( data[n].checkTag(tag) )
		{
			idx = n;
			success = 1;
			break;
		}
	}


	if(success == 0)
	{
		cout << "Data index " << tag << " cannot be replaced because it does not exist!" << endl;
	}
	else
	{
		// now import the data
		data[idx].importArray( newDataVector.data(), newDataVector.size() ); // import a 1D array
	}
}

// Replace single value data
void DataContainer::replaceDataSet( string tag, double newData )
{
	vector <double> newDataVector;
	newDataVector.push_back(newData);
	replaceDataSet( tag, newDataVector );

}

// Pull 1D data set into a vector
void DataContainer::copy1DDataSet( string tag, vector<double> &v )
{
	int Ndim;
	vector<int> sizes;
	getDataSize(tag,Ndim,sizes);
	int dataIdx = findDataSetIndex( tag  );
	v.resize(sizes[0]);

	for(int k = 0; k < sizes[0]; k++)
		v[k] = data[dataIdx].read(k);
}


