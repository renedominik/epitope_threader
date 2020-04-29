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



#ifndef MANUAL_H_
#define MANUAL_H_

#include <iostream>
#include <fstream>

std::ostream &WriteHeader( std::ostream& OUT = std::cout);

std::ostream &WriteHelp( std::string TYPE = "", std::ostream& OUT = std::cout);

std::ostream &WriteManual( std::ostream& OUT = std::cout);

#endif /* MANUAL_H_ */
