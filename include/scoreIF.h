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

#ifndef SCOREIF_H_
#define SCOREIF_H_

#include <vector>

class Molecule;

class ScoreIF
{
public:
	virtual ~ScoreIF(){}

	virtual double Score( const std::vector< Molecule> &MOLS) const = 0;

	virtual double Score( const std::vector< Molecule> &MOLS, const std::vector< std::vector< size_t> > &IDS) const = 0;

	virtual const double &GetCurrentScore() const = 0;
};


#endif /* SCOREIF_H_ */
