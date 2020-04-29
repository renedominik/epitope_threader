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



#ifndef IO_H_
#define IO_H_

#include <iostream>
#include <fstream>

void Open( std::ifstream &STREAM, const std::string &FILE);

void Close( std::ifstream &STREAM);

void Open( std::ofstream &STREAM, const std::string &FILE);

void Close( std::ofstream &STREAM);

std::istream &IgnoreComments( std::istream &STREAM, const char &IDENTIFIER = '#');


#endif /* IO_H_ */
