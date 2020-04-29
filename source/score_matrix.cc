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

#include "../include/score_matrix.h"

#include "../include/io.h"
#include "../include/std_functions.t.h"
#include "../include/string_functions.h"


double ScoreMatrix::operator() ( const std::string &EPITOPE) const
{
	double score = 0.0;
	std::string::const_iterator chr = EPITOPE.begin();
	std::vector< std::map< char, double> >::const_iterator scr = ScoreMatrix::m_Scores.begin();
	for( ; chr != EPITOPE.end(); ++chr, ++scr)
	{
		score += SafeGet(*scr, *chr);
	}
	return score;
}


double ScoreMatrix::operator() ( const std::string &EPITOPE, const std::vector< double> &WEIGHTS) const
{
	std::cout << "### ";
	double score , sum = 0.0;
	std::vector< double>::const_iterator wtr = WEIGHTS.begin();
	std::string::const_iterator chr = EPITOPE.begin();
	std::vector< std::map< char, double> >::const_iterator scr = ScoreMatrix::m_Scores.begin();
	for( ; chr != EPITOPE.end(); ++chr, ++scr, ++wtr)
	{
		score = SafeGet(*scr, *chr);
		std::cout << score << " ";
		sum += *wtr * score;
	}
	std::cout << std::endl;
	return sum;
}


ScoreMatrix ScoreMatrix::Read( std::istream &IN)
{
	DebugFct;
	char c;
	double f;
	int epitope_length, nr_amino_acids, count;

//	IN >> m_Who;
	IN >> epitope_length >> nr_amino_acids;

//	std::cout << m_Who << " " << epitope_length << " " << nr_amino_acids << std::endl;

	m_Scores = std::vector< std::map< char, double> >( epitope_length);

	for( int i = 0; i < epitope_length; ++i)
		for( int j = 0; j < nr_amino_acids; ++j)
		{
			IN >> count >> c >> f;
			m_Scores[count][c] = f;
		}
	return *this;
}



void CalcScoreMatrix( std::istream &IN_MOLS, char EPITOPE_CHAIN, std::istream &IN_AA, std::ostream&OUT)
{
	DebugFct;
	std::vector< Molecule>
		mols = ReadMolecules( IN_MOLS);

	Molecule
		amino_acids = ReadMolecules( IN_AA)[0];

	size_t
		chain_id = Chain2ID( EPITOPE_CHAIN, mols);

	// knowledge-based, electrostatic
	ForceField
		forcefield;

	// single conf, shake, rotamers, md
	CalculateScoreMatrix( mols, chain_id, forcefield, amino_acids, OUT);
}


void CalculateScoreMatrix
(
		const std::vector< Molecule> &MOLS,
		size_t CHAIN_ID,
		const ForceField &FORCEFIELD,
		const Molecule &AMINO_ACIDS,
		std::ostream&OUT
)
{
	DebugFct;

	std::cout << __FUNCTION__ << " chain: " << CHAIN_ID << std::endl;

	size_t
		size;
	double
		sum;

	Residue
		amino_acid;
	const Molecule
		*ref,
		*mol;

	OUT << /*MOLS.size() << " " <<*/ MOLS[ CHAIN_ID].Residues().size() << " " << AMINO_ACIDS.Residues().size() << std::endl;

//	for( std::vector< Residue>::const_iterator itr = AMINO_ACIDS.Residues().begin(); itr != AMINO_ACIDS.Residues().end(); ++itr)
//	{
//		OUT.width( 14);
//		OUT << std::left << itr->Name() << " ";
//	}
//	OUT << std::endl;

	for( size_t i = 0; i < MOLS[CHAIN_ID].Residues().size(); ++i)     // position in epitope
	{
		ref = &MOLS[CHAIN_ID];
		for( size_t j = 0; j < AMINO_ACIDS.Residues().size(); ++j)        // all amino acids
		{
			sum = 0.0;

			amino_acid = AMINO_ACIDS.Residues()[j];  // get from general amino acid definition
			amino_acid.CentroidAt( ref->Residues()[i].Atoms()[0].Position());  // place all atoms at same coordinates // reference residue is centroid

			for( size_t k = 0; k < MOLS.size(); ++k)        // go through all other chains
				if( k != CHAIN_ID)
				{
					mol = &MOLS[k];
					size = mol->Residues().size();
					for( size_t l = 0; l < size; ++l)   // sum interaction with all residues
					{
						sum += FORCEFIELD.Score( mol->Residues()[l], amino_acid);
					}
				}
			OUT.width(4);
			OUT << std::left << i << " " << Converter().SingleLetter( amino_acid.Name()) << "    " << sum << std::endl;
		}
		OUT << std::endl;
	}
	OUT << std::endl;
}




ScoreMatrix
ReadScoreMatrix( std::istream &IN)
{
	ScoreMatrix score;
	score.Read( IN);
	return score;
}


std::vector< ScoreMatrix>
ReadScoreMatrices( std::istream &IN)
{
	// file based version
	std::vector< ScoreMatrix>
		scores;
	std::string
		line;
	std::ifstream
		inmate;
	while( IN >> line)
	{
		if( line.length() > 0)
		{
			Open( inmate, line);
			scores.push_back( ReadScoreMatrix( inmate));
			Close( inmate);
		}
	}
	DebugFct;
	std::cout << __FUNCTION__ << scores.size() << " scoring matrices, peptide size: " << scores[0].Size() << std::endl;
	return scores;
}


std::vector< std::pair< double, ScoreMatrix> >
ReadWeightedScoreMatrix( std::istream &IN)
{
	// file based version
	std::vector< std::pair< double,ScoreMatrix> >
		scores;
	std::vector< std::string>
		cols;
	double
		weight;
	std::string
		line;
	std::ifstream
		inmate;
	while( IN)
	{
		line.clear();
		IN >> line;
		if( line.length() > 0)
		{
			cols = Split( line);
			Open( inmate, cols[1]);
			weight = StringToValue<double>( cols[0]);
			scores.push_back( std::make_pair( weight, ReadScoreMatrix( inmate)));
			Close( inmate);
		}
	}
	DebugFct;
	std::cout << scores.size() << "->" << scores[0].second.Size() << std::endl;
	return scores;
}
