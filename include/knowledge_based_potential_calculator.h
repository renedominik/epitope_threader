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



#ifndef KNOWLEDGE_BASED_POTENTIAL_CALCULATOR_H_
#define KNOWLEDGE_BASED_POTENTIAL_CALCULATOR_H_

#include <iostream>
#include <vector>

void
KB_Potential_Backbone( std::istream &PDB_LIST, const std::string &OUTPATH);

void
KB_Potential_AAPair( std::istream &PDB_LIST, const std::string &OUTPATH);

void
KB_Potential_AASolvation( std::istream &PDB_LIST, const std::string &OUTPATH);

std::vector< double>
DistributionToPotential( const std::vector< double> &DISTRIBUTION);


#endif /* KNOWLEDGE_BASED_POTENTIAL_CALCULATOR_H_ */
