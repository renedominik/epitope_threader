/////////////////////////////////////////////////////////////////////////
//
//  This file is part of the EpitopeThreader program
//
//  Copyright (C) 2013 by Rene Staritzbichler
//  renedominik@yahoo.de
//
//
//
//
// date:   Dez 16, 2013
/////////////////////////////////////////////////////////////////////////


#ifndef STRING_FUNCTIONS
#define STRING_FUNCTIONS

#include <vector>
#include <sstream>
#include <iostream>
#include <cassert>
#include <limits>

// translate a char array into a string, removing eventual non alphanumerical endings
void BufferToFileName( const char BUFFER[], std::string &STR, const int &SIZE);


// trims spaces from beginning and end of a copied string
std::string Trim( const std::string &STRING);


// test whether string is numerical (double, double, size_t, int) with '.', '-', leading and tailing spaces are allowed
bool IsNumerical( const std::string &STRING);


std::vector< std::string> Split( const std::string &STRING, const std::string &SPLITTER = " ");

std::vector< std::string> Split( const std::string &STRING, const std::string &SPLITTER_A, const std::string &SPLITTER_B);


// converts std::string to object of template class T1
template< class T1>
T1 StringToValue( const std::string &STRING)
{
	// remove tabs, spaces, new lines, etc.
	const std::string string( Trim( STRING));

	// if string is empty an undefined is returned
	if( string.empty())
	{ return std::numeric_limits< T1>::max();}

	assert( IsNumerical( string));
	T1 new_t;
	std::stringstream iss( string);
	iss >> new_t;
	return new_t;
}


template<class T>
std::string ValueToString(const T& VALUE)
{
	std::stringstream strm;
	strm << VALUE;
	return strm.str();
}


template<class T>
std::string ValueToString(const T& VALUE, const int & TOTAL, const int &PRECISION)
{
	std::string str( TOTAL, ' ');

	char
		*nr = new char[4],
		*prec = new char[4],
		*char_ptr = new char[ TOTAL];

	sprintf( prec, "%d", PRECISION);
	sprintf( nr, "%d", TOTAL);

	std::string
		format = "%";
	format.append( nr);
	format +=  ".";
	format.append( prec);
	format += "f";

	//      std::cout << "format: " << format << std::endl;

	int n = sprintf( char_ptr, format.c_str(), VALUE);

	// the resulting char array might be longer than TOTAL
	// adjust the precision in that case
	if( n > TOTAL)
	{
		//	  std::cout << "long stuff" << std::endl;
		int new_size = std::max( 0, PRECISION - ( n - TOTAL));
		sprintf( prec, "%d", new_size);
		format = "%";
		format.append( nr);
		format +=  ".";
		format.append( prec);
		format += "f";
		//	  std::cout << "new format: " << format << std::endl;

		sprintf( char_ptr, format.c_str(), VALUE);
	}

	// pass char array to correct sized string
	std::copy( char_ptr, char_ptr + TOTAL, str.begin());

	return str;
}


// transforms given string into lower case
std::string &AllToLowerCase( std::string &STRING);

// returns a lower case copy of string
std::string ToLower( const std::string &STRING);

// transforms given string to upper case
std::string &AllToUpperCase( std::string &STRING);

// returns an upper case copy of string
std::string ToUpper( const std::string &STRING);

// translates strings up to size 9 to size_t; needed for example for sending messages among processes -msgsnd, msgrcv cannot handle strings correctly on all OS
size_t StringToSizeT( const std::string &STRING);

// translates size_t to string up to size 9
std::string SizeTToString( const size_t &VALUE);

// searches for spaces and removes them from the string
inline
std::string RemoveSpacesFromString( const std::string &STRING)
{
	std::string cleaned_string;
	for( size_t i = 0 ;i < STRING.size(); i++ )
	{
		if( STRING[i] != ' ') cleaned_string.push_back( STRING[i]);
	}
	return cleaned_string;
  }

inline
bool
IsCapitolLetter( const char &CHAR)
{
	size_t nr( CHAR);
	return nr >= 65 && nr <= 90;
}

inline
bool
IsSmallLetter( const char &CHAR)
{
	size_t nr( CHAR);
	return nr >= 97 && nr <= 122;
}

inline
bool
IsLetter( const char &CHAR)
{
	size_t nr( CHAR);
	return ( nr >= 65 && nr <= 90) || ( nr >= 97 && nr <= 122);
}

inline
bool
IsNumber( const char &CHAR)
{
	size_t nr( CHAR);
	return nr >= 48 && nr <= 57;
}

inline
std::string
SortAndGlue( const std::string &S1, const std::string &S2)
{
	if( S1 < S2)
	{
		return S1 + S2;
	}
	return S2 + S1;
}

inline
std::string
ReplaceAll( const std::string &ORIG, const std::string &TO_REPLACE, const std::string &REPLACE_WITH)
{
	std::string copy = ORIG;
	size_t pos = 0;
	while( (pos = ORIG.find( TO_REPLACE, pos)) != std::string::npos)
	{
		copy.replace( pos, REPLACE_WITH.size(), REPLACE_WITH);
		pos += REPLACE_WITH.length();
	}
	return copy;
}

#endif
