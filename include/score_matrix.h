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


#ifndef SCORE_MATRIX_H_
#define SCORE_MATRIX_H_

#include <map>

#include "molecule.h"
#include "forcefield.h"



class ScoreMatrix
{
protected:
	std::vector< std::map< char, double> >  m_Scores;    //< position specific score
	std::string                             m_Who;       //<  who am i ?? ...
	// instead of map one could exploit int( char) - offset
public:
	ScoreMatrix(){}


	virtual ~ScoreMatrix(){}

	const std::string & Who() const
	{
		return m_Who;
	}

	int Size() const
	{
		return m_Scores.size();
	}

	ScoreMatrix Read( std::istream &IN);

	virtual double operator() ( const std::string &EPITOPE) const;

	virtual double operator() ( const std::string &EPITOPE, const std::vector< double> &WEIGHTS) const;

};


void
CalcScoreMatrix( std::istream &IN_PDB, char EPITOPE_CHAIN, std::istream &IN_AA, std::ostream&OUT);

void
CalculateScoreMatrix
(
		const std::vector< Molecule> &MOLS,
		size_t CHAIN_ID,
		const ForceField &FORCEFIELD,
		const Molecule &AMINO_ACIDS,
		std::ostream&OUT
);



ScoreMatrix
ReadScoreMatrix( std::istream &IN);


std::vector< ScoreMatrix>
ReadScoreMatrices( std::istream &IN);

std::vector< std::pair< double, ScoreMatrix> >
ReadWeightedScoreMatrix( std::istream &IN);



//	object
//		parameters = io::ReadParameterFile( PARAM_FILE);
//
//	object
//		topology = io::ReadTopologyFile( TOP_FILE);



#endif /* SCORE_MATRIX_H_ */
