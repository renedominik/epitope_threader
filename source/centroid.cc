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

#include "../include/centroid.h"

#include "../include/molecule.h"
#include "../include/std_functions.t.h"
#include "../include/global.h"

void CbetaCentroids( std::istream &IN, std::ostream &OUT)
{
	DebugFct;
	std::vector< Molecule>
		mols = ReadVector<Molecule>( IN);

	CbetaCentroids( mols);

	WriteVector( mols, OUT);
}

void CalphaCentroids( std::istream &IN, std::ostream &OUT)
{
	DebugFct;
	std::vector< Molecule>
		mols = ReadVector<Molecule>( IN);

	CalphaCentroids( mols);

	WriteVector( mols, OUT);
}



