/*
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Copyright (C) 2017-2019 Brian D. Hong
 *
*/

/*
 * importExport.cpp
 *
 *  Created on: Sep 23, 2017
 *      Author: brian
 *
 * A simple and lightweight I/O system for reading parameters and strings
 *
 *	1. The goal of this implementation is to avoid "overthinking" the I/O system
 *	   the input format is a simple text list of terms like
 *
 *	   bx = 20e2
 *	   fileName = dataFile.txt
 *	   v = [1, 3, 2.8]
 *
 *	   these are used by the data container to store parm values as needed
 *
 *	2. This code does allow you to make errors, but tries to avoid the most common ones
 *	3. The input will read normal values (eg. 1.3 ) or scientific values (e.g. 1e-5 )
 *	4. The code will automatically also save a string, so that it can be pulled in your code
 *	5. 2D arrays are also supported
 *	6. The output formats are meant to be read by matlab, but only using the simple
 *	   standard .m file, not a storage format, such as .mat
 *	7. As a result of (6), large output arrays tend to be bulky and read slowly by matlab
 *  8. Any line in the input file with // in it anywere is automatically ignored
 *
 *
 */

#include "importExport.hpp"


// Looks through the parameter file and finds keys
void findKeys( ifstream &fileID, vector<string> & keys, int verbose )
{

	// clear the string
	keys.clear();

	// Get the process
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	string outString;
	bool success = 0;

	string str;
	size_t foundKey;

	stringstream valueStr;

	fileID.clear();
	fileID.seekg(0, ios::beg);

	while (std::getline(fileID, str))
	{

		char * checkChar = new char[ str.length()+1];
		std::strcpy(checkChar, str.c_str());

		// skip if found "//"
		size_t foundEscape = str.find("//");

		// if the first value is alphabetical and not an escape, give a key
		if ( foundEscape == string::npos && isalpha( checkChar[0] ) )
		{

			// once the key is found, search for the number
			char * strChar = new char[ str.length()+1];
			std::strcpy(strChar, str.c_str());

			int i = 0;
			char c;

			valueStr.str("");

			while (strChar[i])
			{
				c=strChar[i];

				// skip all nonalphabetical characters

				if( c != ' ' && c != '=' )
				{

					valueStr << c;
				}
				else if ( c == '=' )
				{
					break;
				}


				i++;
			}

			outString = valueStr.str();
			keys.push_back(outString);
			outString.clear();
			success = 1;

			delete [] strChar;

		}
		delete [] checkChar;

	}

	if ( success == 1 ){
		if( procID == ROOT_ID && verbose >= 2  ){
			cout << "Successfully found keys:";
			for(int k = 0; k < keys.size(); k++)
			{
				cout << " " << keys[k] ;
			}
			cout << "." << endl << endl;
		}

	}else {
		if( procID == ROOT_ID && verbose >= 0)  cout << "!!! No keys found. String is left empty. !!!" << endl;
	}

}




// determines whether or not the keys are 1D or 2D and their size
void findKeySize( ifstream &fileID, string key, int & dim, int & N1, int & N2, bool & stringFound, int verbose )
{

	stringFound = false;
	dim = 0;
	N1 = 0;
	N2 = 0;


	// Get the process
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);


	double success = 0;
	double value;

	string str;
	size_t foundKey;
	stringstream numStr;

	fileID.clear();
	fileID.seekg(0, ios::beg);

	int k = 0;
	int i;

	bool condition1 = false;
	bool condition2 = false;

	while (std::getline(fileID, str))
	{

		// look through file for the key
		foundKey = str.find(key);

		char * checkChar = new char[ str.length()+1];
		std::strcpy(checkChar, str.c_str());


		size_t foundEscape = str.find("//");

		if ( foundKey != string::npos
			 && ( foundKey == 0 || isspace(checkChar[foundKey-1]) )
			 && ( isspace(checkChar[foundKey+key.length()]) ||( checkChar[foundKey+key.length()] == '=') )
			 && foundEscape == string::npos )
		{

			N2 = 0;
			N1 = 1;
			condition1 = true;

			// erase the key, as it might contain a digit
			str.erase(foundKey, key.length());

			// add something on the end so it doesnt miss the last number
			// but dont add a letter or a number!
			str.append(" | ");

			// once the key is found, search for the number
			char * strChar = new char[ str.length()+1];
			std::strcpy(strChar, str.c_str());

			success = 0.5;

			// once the key is found, search for the values
			char c;
			char cprev = 'e';
			i = 0;

			while (strChar[i])
			{
				c=strChar[i];

				if( isdigit(c) || c == '.' || ( c == 'e' && isdigit(cprev) ) || c == '+' || c == '-' )
				{
					numStr << c;
				} else if ( numStr.str().length() > 0) {
					N2++;

					numStr.str( std::string() );
					numStr.clear();

				} else if ( isalpha(c)  ) { // if alphabetic indicate that a string was found
					stringFound = true;
				}
				else
				{
					// if its not a number or a . and no number is stored then just skip

				}

				cprev = c;
				i++;
			}

			// after the first line has been found, check to see if there is another (i.e., is it a 2D array?)
			while (std::getline(fileID, str))
			{

				char * checkChar2 = new char[ str.length()+1];
				std::strcpy(checkChar2, str.c_str());

				bool condition2Local = false;

				i = 0;
				// if a number is found, set to true
				while (checkChar2[i])
				{
					if( isdigit(checkChar2[i]) ){
						condition2Local = true;
					}
					i++;
				}


				i = 0;
				// if an equal is found, set back to false
				while (checkChar2[i])
				{
					if( (isalpha(checkChar2[i]) && checkChar2[i]!='e') || checkChar2[i] == '=' ){
						condition2Local = false;
					}
					i++;
				}

				if(condition2Local == true)
				{
					N1 ++;
				}
				else
				{
					break;
				}

				delete [] checkChar2;
			}
			delete [] strChar;
		}
		delete [] checkChar;
	}

	if(N1 > 1)
		condition2 = true;

	if(condition2 == true )
	{
		dim = 2;
	}else if (condition1 == true){
		dim = 1;
	}else{
		dim = 0;
	}

	// success if k got to n
	if( procID == ROOT_ID && verbose >= 2) cout << "Key " << key << " is " << dim << " dimensional with size " << N1 << " x "  << N2 << endl;

}


// This finds a single value from a line that contains the key "string"
// Doesnt read lines that contain "//"
double findValue( ifstream &fileID, string key, int verbose  )
{

	// Get the process
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	double value;
	bool success = 0;

	string str;
	size_t foundKey;

	stringstream numStr;

	fileID.clear();
	fileID.seekg(0, ios::beg);

	while (std::getline(fileID, str))
	{

		// look through file for the key
		foundKey = str.find(key);


		char * checkChar = new char[ str.length()+1];
		std::strcpy(checkChar, str.c_str());

		size_t foundEscape = str.find("//");


		if ( foundKey != string::npos
			 && ( foundKey == 0 || isspace(checkChar[foundKey-1]) )
			 && ( isspace(checkChar[foundKey+key.length()]) ||( checkChar[foundKey+key.length()] == '=') )
			 && foundEscape == string::npos )
		{

			// erase the key, as it might contain a digit
			str.erase(foundKey, key.length());

			// once the key is found, search for the number
			char * strChar = new char[ str.length()+1];
			std::strcpy(strChar, str.c_str());

			int i = 0;
			char c;
			char cprev = 'e';


			while (strChar[i])
			{
				c=strChar[i];

				if( isdigit(c) || c == '.' || ( c == 'e' && isdigit(cprev) ) || c == '+' || c == '-' )
				{
					numStr << c;
				}
				else if ( numStr.str().length() > 0)
				{
					break;

				}

				cprev = c;

				i++;
			}
			numStr >> value;
			success = 1;

			delete [] strChar;

		}
		delete [] checkChar;
	}

	if ( success == 1 ){
		if( procID == ROOT_ID && verbose >= 1)  cout << "Successfully set " << key << " = " << value << endl;
	}else {
		if( procID == ROOT_ID &&  verbose >= 0)  cout << "!!! Parameter " << key << " not found. Parameter arbitrarily set to 0. !!!" << endl;
		value = 0;
	}

	return value;
}



// This finds a vector of values from a line that contains the key "string"
void findArray( ifstream &fileID, string key, int N, double * values, int verbose ){

	// Get the process
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	double success = 0;
	double value;

	string str;
	size_t foundKey;
	stringstream numStr;

	fileID.clear();
	fileID.seekg(0, ios::beg);

	int k = 0;
	while (std::getline(fileID, str))
	{

		// look through file for the key
		foundKey = str.find(key);

		char * checkChar = new char[ str.length()+1];
		std::strcpy(checkChar, str.c_str());


		size_t foundEscape = str.find("//");

		if ( foundKey != string::npos
			 && ( foundKey == 0 || isspace(checkChar[foundKey-1]) )
			 && ( isspace(checkChar[foundKey+key.length()]) ||( checkChar[foundKey+key.length()] == '=') )
			 && foundEscape == string::npos )
		{

			// erase the key, as it might contain a digit
			str.erase(foundKey, key.length());

			// add something on the end so it doesnt miss the last number
			str.append(" | ");

			// once the key is found, search for the number
			char * strChar = new char[ str.length()+1];
			std::strcpy(strChar, str.c_str());



			success = 0.5;

			// once the key is found, search for the values
			char c;
			char cprev = 'e';
			int i = 0;

			while (strChar[i])
			{
				c=strChar[i];

				if( isdigit(c) || c == '.' || ( c == 'e' && isdigit(cprev) ) || c == '+' || c == '-' )
				{
					numStr << c;
				} else if ( numStr.str().length() > 0) {
					numStr >> values[k];

					k = k+1;
					if (k > N )
					{
						success = -1;
						break;
					}
					numStr.str( std::string() );
					numStr.clear();

				} else {
					// if its not a number or a . and no number is stored then just skip
				}

				cprev = c;
				i++;
			}

			delete [] strChar;
		}

		delete [] checkChar;
	}

	// success if k got to n
	if(k == N && success != -1)
	{
		if( procID == ROOT_ID && verbose >= 1) {
			cout << "   "; print1DArrayLine(values,N,4,key);
		}
		success = 1;
	}

	if ( success == 0 ){
		if( procID == ROOT_ID && verbose >= 0)  cout << "!!! Vector " << key << " not found. Vector arbitrarily set to 0. !!!" << endl;
	}else if (success == 0.5){
		if( procID == ROOT_ID && verbose >= 0)  cout << "!!! Vector " << key << " did not have enough elements! Remaining elements arbitrarily set to 0. !!!" << endl;
	}else if (success == -1){
		if( procID == ROOT_ID && verbose >= 0)  cout << "!!! Vector " << key << " had too many elements! Returning vector with only "<< N << " elements !!!"<< endl;
	}

}



// This is a poor re-implementation of the above findArray
// --> really this should just call findArray of doubles then convert them to integers
void findArray( ifstream &fileID, string key, int N, int * values, int verbose  )
{

	// Get the process
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	double success = 0;
	double value;

	string str;
	size_t foundKey;
	stringstream numStr;

	fileID.clear();
	fileID.seekg(0, ios::beg);

	int k = 0;
	while (std::getline(fileID, str))
	{

		// look through file for the key
		foundKey = str.find(key);

		char * checkChar = new char[ str.length()+1];
		std::strcpy(checkChar, str.c_str());


		size_t foundEscape = str.find("//");

		if ( foundKey != string::npos
			 && ( foundKey == 0 || isspace(checkChar[foundKey-1]) )
			 && ( isspace(checkChar[foundKey+key.length()]) ||( checkChar[foundKey+key.length()] == '=') )
			 && foundEscape == string::npos )
		{

			// erase the key, as it might contain a digit
			str.erase(foundKey, key.length());

			// add something on the end so it doesnt miss the last number
			str.append(" | ");

			// once the key is found, search for the number
			char * strChar = new char[ str.length()+1];
			std::strcpy(strChar, str.c_str());

			success = 0.5;

			// once the key is found, search for the values
			char c;
			char cprev = 'e';
			int i = 0;

			while (strChar[i])
			{
				c=strChar[i];

				if( isdigit(c) || c == '.' || ( c == 'e' && isdigit(cprev) ) || c == '+' || c == '-' )
				{
					numStr << c;
				} else if ( numStr.str().length() > 0) {
					double valueTemp;
					numStr >> valueTemp;

					values[k] = round(valueTemp);

					k = k+1;
					if (k > N )
					{
						success = -1;
						break;
					}
					numStr.str( std::string() );
					numStr.clear();

				} else {
					// if its not a number or a . and no number is stored then just skip
				}

				cprev = c;
				i++;
			}
			delete [] strChar;
		}
		delete [] checkChar;

	}

	// success if k got to n
	if(k == N && success != -1)
	{
		if( procID == ROOT_ID && verbose >= 1) cout << "Vector " << key << " successfully copied with " << N << " elements."  << endl;
		success = 1;
	}

	if ( success == 0 ){
		if( procID == ROOT_ID &&  verbose >= 0)  cout << "!!! Vector " << key << " not found. Vector arbitrarily set to 0. !!!" << endl;
	}else if (success == 0.5){
		if( procID == ROOT_ID &&  verbose >= 0)  cout << "!!! Vector " << key << " did not have enough elements! Remaining elements arbitrarily set to 0. !!!" << endl;
	}else if (success == -1){
		if( procID == ROOT_ID &&  verbose >= 0)  cout << "!!! Vector " << key << " had too many elements! Returning vector with only "<< N << " elements !!!"<< endl;
	}

}

// This finds a vector of values from a line that contains the key "string"
void find2DArray( ifstream &fileID, string key, int Nj, int Nk, double ** values, int verbose )
{

	for(int j = 0; j < Nj; j++)
	{
		for(int k = 0; k < Nk; k++)
		{
			values[j][k] = 0.0;
		}
	}

	// Get the process
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	double success = 0;
	double value;

	string str;
	size_t foundKey;
	stringstream numStr;

	fileID.clear();
	fileID.seekg(0, ios::beg);

	int k = 0;
	int j = 0;
	while (std::getline(fileID, str))
	{

		// look through file for the key
		foundKey = str.find(key);

		char * checkChar = new char[ str.length()+1];
		std::strcpy(checkChar, str.c_str());


		size_t foundEscape = str.find("//");

		if ( foundKey != string::npos
			 && ( foundKey == 0 || isspace(checkChar[foundKey-1]) )
			 && ( isspace(checkChar[foundKey+key.length()]) ||( checkChar[foundKey+key.length()] == '=') )
			 && foundEscape == string::npos )
		{

			// erase the key, as it might contain a digit
			str.erase(foundKey, key.length());

			// add something on the end so it doesnt miss the last number
			str.append(" | ");

			// once the key is found, search for the number
			char * strChar = new char[ str.length()+1];
			std::strcpy(strChar, str.c_str());

			success = 0.5;

			// once the key is found, search for the values
			char c;
			char cprev = 'e';
			int i = 0;

			while (strChar[i])
			{
				c=strChar[i];

				if( isdigit(c) || c == '.' || ( c == 'e' && isdigit(cprev) ) || c == '+' || c == '-' )
				{
					numStr << c;
				} else if ( numStr.str().length() > 0) {
					numStr >> values[j][k];

					k = k+1;
					if (k > Nk )
					{
						success = -1;
						break;
					}
					numStr.str( std::string() );
					numStr.clear();

				} else {
					// if its not a number or a . and no number is stored then just skip
				}

				cprev = c;
				i++;
			}

			delete [] strChar;
			j++; k = 0;

			while( j < Nj )
			{

				std::getline(fileID, str);

				// add something on the end so it doesnt miss the last number
				str.append(" end line ");


				// copy to cstring
				char * strChar2 = new char[ str.length()+1];
				std::strcpy(strChar2, str.c_str());

				success = 0.5;

				// search for the values
				char c;
				char cprev = 'e';
				int i = 0;

				while (strChar2[i])
				{
					c=strChar2[i];

					if( isdigit(c) || c == '.' || ( c == 'e' && isdigit(cprev) ) || c == '+' || c == '-' )
					{
						numStr << c;
					} else if ( numStr.str().length() > 0) {
						numStr >> values[j][k];


						k = k+1;
						if (k > Nk )
						{
							success = -1;
							break;
						}
						numStr.str( std::string() );
						numStr.clear();

					} else {
						// if its not a number or a . and no number is stored then just skip
					}

					cprev = c;
					i++;
				}

				delete [] strChar2;
				j++; k = 0;
			}
		}

		delete [] checkChar;
	}

	// success if k got to n
	if(j == Nj && success != -1)
	{
		if( procID == ROOT_ID && verbose >= 1 ) cout << "Matrix " << key << " successfully copied with " << Nj << " x " << Nk << "  elements."  << endl;
		success = 1;
	}

	if ( success == 0 ){
		if( procID == ROOT_ID && verbose >= 0)  cout << "!!! Vector " << key << " not found. Vector arbitrarily set to 0. !!!" << endl;
	}else if (success == 0.5){
		if( procID == ROOT_ID && verbose >= 0 )  cout << "!!! Vector " << key << " did not have enough elements! Remaining elements arbitrarily set to 0. !!!" << endl;
	}else if (success == -1){
		if( procID == ROOT_ID && verbose >= 0)  cout << "!!! Vector " << key << " had too many elements! Returning vector with only "<< Nj << " x " << Nk << "  elements."  << endl;
	}

}

// This finds a string that follows the key
string findString( ifstream &fileID, string key, int verbose  ){

	// Get the process
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	string outString;
	bool success = 0;

	string str;
	size_t foundKey;

	stringstream valueStr;

	fileID.clear();
	fileID.seekg(0, ios::beg);

	while (std::getline(fileID, str))
	{

		// look through file for the key
		foundKey = str.find(key);


		char * checkChar = new char[ str.length()+1];
		std::strcpy(checkChar, str.c_str());

		size_t foundEscape = str.find("//");


		if ( foundKey != string::npos
			 && ( foundKey == 0 || isspace(checkChar[foundKey-1]) )
			 && ( isspace(checkChar[foundKey+key.length()]) ||( checkChar[foundKey+key.length()] == '=') )
			 && foundEscape == string::npos )
		{

			// erase the key, as it might contain a digit
			str.erase(foundKey, key.length());

			// once the key is found, search for the number
			char * strChar = new char[ str.length()+1];
			std::strcpy(strChar, str.c_str());

			int i = 0;
			char c;
			char cprev = 'e';


			while (strChar[i])
			{
				c=strChar[i];

				// skip '=' or spaces
				if( c != ' ' && c!= '=' && c != '\r' )
				{
					valueStr << c;
				}
				else if ( valueStr.str().length() > 0)
				{
					break;
				}
				cprev = c;
				i++;
			}
			outString = valueStr.str();
			success = 1;

			delete [] strChar;

		}
		delete [] checkChar;
	}

	if ( success == 1 ){
		if( procID == ROOT_ID && verbose >= 1)
		{
			cout << std::left;
			cout << std::setw(20)  << "   "  +key << "    = " << outString << endl;

		}
	}else {
		if( procID == ROOT_ID &&  verbose >= 0)  cout << "!!! String" << key << " not found. String is left empty. !!!" << endl;
		outString = "";
	}

	return outString;
}

// Prints a matlab formatted value (note: this annoying version also tells matlab to print the value when it runs)
void printMatlabVariable( ofstream &fileID, string key, double value  )
{
	// Get the process
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	int precision = 16;

	if (fileID.is_open()) {
		fileID << std::setprecision( precision );

	    fileID << key << " = " << value << ";" << endl;
	    fileID << "display( [ \'"<< key << " = \', num2str(" << key <<  ") ] ) " << ";" << endl << endl;

	}else{
		if( procID == ROOT_ID )  cout << "Matlab variable file was not open. Cannot print output." << endl;
	}
}

// Prints a matlab formatted value (note: this version doesn't print when loaded into matlab)
void printMatlabVariableSimple( ofstream &fileID, string key, double value )
{
	// Get the process
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	int precision = 16;

	if (fileID.is_open()) {
		fileID << std::setprecision( precision );

	    fileID << key << " = " << value << ";" << endl << endl;
	}else{
		if( procID == ROOT_ID )  cout << "Matlab variable file was not open. Cannot print output." << endl;
	}
}

// Prints a matlab formatted 1D array (note: this annoying version also tells matlab to print the value when it runs)
void printMatlabArray( ofstream &fileID, string key, double * values, int N )
{
	// Get the process
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	int precision = 16;

	if (fileID.is_open()) {
		fileID << std::setprecision( precision );

	    fileID << key << " = [ ";
	    for(int k = 0; k < N; k ++ )
	    {
	    	if (k < N-1)
	    		fileID << values[k] << ", ";
			else
				fileID << values[k] << " ];" << endl;
	    }
	    fileID << "display( [ \'"<< key << " = \', num2str(" << key <<  ") ] ) " << ";" << endl << endl;

	}else{
		if( procID == ROOT_ID )  cout << "Matlab variable file was not open. Cannot print output." << endl;
	}
}

// Prints a matlab formatted array ((note: this version doesn't print when loaded into matlab)
void printMatlabArraySimple( ofstream &fileID, string key, double * values, int N )
{
	// Get the process
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	int precision = 16;

	if (fileID.is_open()) {
		fileID << std::setprecision( precision );

	    fileID << key << " = [ ";
	    for(int k = 0; k < N; k ++ )
	    {
	    	if (k < N-1)
	    		fileID << values[k] << ", ";
			else
				fileID << values[k] << " ];" << endl;
	    }
	    fileID  << endl;
	}else{
		if( procID == ROOT_ID)   cout << "Matlab variable file was not open. Cannot print output." << endl;
	}
}

// Prints a matlab formatted 2D array (note: this annoying version also tells matlab to print the value when it runs)
void printMatlab2DArray( ofstream &fileID, string key, double ** values, int N1, int N2 )
{
	// Get the process
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	int precision = 16;

	if (fileID.is_open()) {
		fileID << std::setprecision( precision );

	    fileID << key << " = [ ";
	    for(int k = 0; k < N1; k++ )
	    {
	    	for( int j = 0; j < N2; j++ )
	    	{

		    	if (j < N2-1)
		    		fileID << values[k][j] << ", ";
				else
					fileID << values[k][j] << " ;";
	    	}
	    	if (k < N1 - 1)
	    		fileID << endl << "            ";
	    	else
	    		fileID << "];" << endl ;


	    }
	    fileID << "display( \'"<< key << " = \');" << endl;
	    fileID << "display( num2str(" << key <<  " ) )" << ";" << endl << endl;

	}else{
		if( procID == ROOT_ID )  cout << "Matlab variable file was not open. Cannot print output." << endl;
	}
}

// Prints a matlab formatted 2D array (note: this version doesn't print when loaded into matlab)
void printMatlab2DArraySimple( ofstream &fileID, string key, double ** values, int N1, int N2  )
{
	// Get the process
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	int precision = 16;

	if (fileID.is_open()) {
		fileID << std::setprecision( precision );

	    fileID << key << " = [ ";
	    for(int k = 0; k < N1; k++ )
	    {
	    	for( int j = 0; j < N2; j++ )
	    	{

		    	if (j < N2-1)
		    		fileID << values[k][j] << ", ";
				else
					fileID << values[k][j] << " ;";
	    	}
	    	if (k < N1 - 1)
	    		fileID << endl << "            ";
	    	else
	    		fileID << "];" << endl ;


	    }

	    fileID << endl;

	}else{
		if( procID == ROOT_ID )  cout << "Matlab variable file was not open. Cannot print output." << endl;
	}
}

// Prints a matlab formatted 3D array (note: this version doesn't print when loaded into matlab)
void printMatlab3DArraySimple( ofstream &fileID, string key, double *** values, int N1, int N2, int N3  )
{
	// Get the process
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	int precision = 16;


	if( procID == ROOT_ID )
	{
		if (fileID.is_open()) {
			fileID << std::setprecision( precision );


			fileID << key << " = zeros( " << N1 << ", " << N2 << ", " << N3 << ");" << endl;
			for (int i = 0; i < N3; i++)
			{
				fileID << key << "(:,:," << i+1 << ") = [ ";
				for(int k = 0; k < N1; k++ )
				{
					for( int j = 0; j < N2; j++ )
					{

						if (j < N2-1)
							fileID << values[k][j][i] << ", ";
						else
							fileID << values[k][j][i] << " ;";
					}
					if (k < N1 - 1)
						fileID << endl << "            ";
					else
						fileID << "];" << endl ;
				}
			}

			fileID << endl;

		}else{
			if( procID == ROOT_ID )  cout << "Matlab variable file was not open. Cannot print output." << endl;
		}
	}
}

// Prints a matlab formatted 3D array (note: this version doesn't print when loaded into matlab)
void printMatlab4DArraySimple( ofstream &fileID, string key, double **** values, int N1, int N2, int N3, int N4 )
{
	// Get the process
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	int precision = 16;

	if( procID == ROOT_ID )
	{
		if (fileID.is_open()) {
			fileID << std::setprecision( precision );


			fileID << key << " = zeros( " << N1 << ", " << N2 << ", " << N3 <<  ", " << N4 << ");" << endl;
			for (int m = 0; m < N4; m++)
			{
				for (int i = 0; i < N3; i++)
				{
					fileID << key << "(:,:," << i+1 << "," << m+1 << ") = [ ";
					for(int k = 0; k < N1; k++ )
					{
						for( int j = 0; j < N2; j++ )
						{

							if (j < N2-1)
								fileID << values[k][j][i][m] << ", ";
							else
								fileID << values[k][j][i][m] << " ;";
						}
						if (k < N1 - 1)
							fileID << endl << "            ";
						else
							fileID << "];" << endl ;
					}
				}
			}
			fileID << endl;

		}else{
			if( procID == ROOT_ID )  cout << "Matlab variable file was not open. Cannot print output." << endl;
		}
	}
}

// Prints a matlab formatted 2D array from a 2D vector of doubles (doesnt print)
void printMatlab2DArray( ofstream &fileID, string key, vector<vector<double>> & values  )
{
	int N1 = values.size();

	// Get the process
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	int precision = 16;

	if (fileID.is_open()) {
		fileID << std::setprecision( precision );

	    fileID << key << " = [ ";
	    for(int k = 0; k < N1; k++ )
	    {
	    	int N2 = values[k].size();
	    	for( int j = 0; j < N2; j++ )
	    	{

		    	if (j < N2-1)
		    		fileID << values[k][j] << ", ";
				else
					fileID << values[k][j] << " ;";
	    	}
	    	if (k < N1 - 1)
	    		fileID << endl << "            ";
	    	else
	    		fileID << "];" << endl ;


	    }

	    fileID << endl;

	}else{
		if( procID == ROOT_ID )  cout << "Matlab variable file was not open. Cannot print output." << endl;
	}
}

// Prints a matlab formatted 1D array from a 1D vector of doubles (doesnt print)
void printMatlab1DArray( ofstream &fileID, string key, vector<double> & values)
{
	int N = values.size();

	// Get the process
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	int precision = 16;

	if (fileID.is_open()) {
		fileID << std::setprecision( precision );

	    fileID << key << " = [ ";
	    for(int k = 0; k < N; k ++ )
	    {
	    	if (k < N-1)
	    		fileID << values[k] << ", ";
			else
				fileID << values[k] << " ];" << endl;
	    }
	    fileID  << endl;
	}else{
		if( procID == ROOT_ID)   cout << "Matlab variable file was not open. Cannot print output." << endl;
	}
}









