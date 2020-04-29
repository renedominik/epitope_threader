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

#include "../include/io.h"

#include <stdlib.h>

void Open( std::ifstream &STREAM, const std::string &FILE)
{
	STREAM.open( FILE.c_str());
	if(!STREAM)
	{
		std::cerr << "ERROR: " << FILE << " could not be opened for input" << std::endl;
		exit(1);
	}
}

void Close( std::ifstream &STREAM)
{
	STREAM.close();
	STREAM.clear();
}

void Open( std::ofstream &STREAM, const std::string &FILE)
{
	STREAM.open( FILE.c_str());
	if(!STREAM)
	{
		std::cerr << "ERROR: " << FILE << " could not be opened for output" << std::endl;
		exit(1);
	}
}

void Close( std::ofstream &STREAM)
{
	STREAM.close();
	STREAM.clear();
}

	std::istream &IgnoreComments( std::istream &STREAM, const char &IDENTIFIER)
	{
	    char c;
	    std::string tmp;

	    while( (c = STREAM.get()) == IDENTIFIER || c == '\n' || c == ' ')
	    {
	        if( c == IDENTIFIER)
	        {
	            std::getline( STREAM, tmp);
	        }
	    }
	    if( !STREAM.eof())
	    {
	        STREAM.putback( c);
	    }
	    return STREAM;
	}


