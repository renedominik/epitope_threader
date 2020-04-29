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



#ifndef SCORE_MANAGER_H_
#define SCORE_MANAGER_H_

#include "scoreIF.h"

#include "score_kb.h"
#include "constraint.h"
#include "density.h"

#include "io.h"

class ScoreManager
{
private:
	std::vector< ScoreIF*>  m_Scores;   //smart pointer?
	std::vector< double>	m_Weights;
	std::vector< std::string>  m_Names;

public:
	ScoreManager( const std::vector< ScoreIF*> &SCORES = std::vector< ScoreIF*>(), const std::vector< double> &WEIGHTS = std::vector< double>())
	: m_Scores(),
	m_Weights(),
	m_Names()
	{}

	virtual ~ScoreManager(){}

	virtual double Score( const std::vector< Molecule> &MOLS) const
	{
		double score = 0.0;
		std::vector<double>::const_iterator wtr = m_Weights.begin();
		for( std::vector< ScoreIF*>::const_iterator itr = m_Scores.begin(); itr != m_Scores.end(); ++itr, ++wtr)
		{
			score += *wtr * (*itr)->Score( MOLS);
		}
		return score;
	}

	virtual double Score( const std::vector< Molecule> &MOLS, const std::vector< std::vector< size_t> > &IDS) const
	{
		double
			score = 0.0;
		std::vector<double>::const_iterator
			wtr = m_Weights.begin();
//		std::vector< std::string>::const_iterator
//			ntr = m_Names.begin();
		for( std::vector< ScoreIF*>::const_iterator itr = m_Scores.begin(); itr != m_Scores.end(); ++itr, ++wtr)
		{
			score += *wtr * (*itr)->Score( MOLS, IDS);
//			std::cout << *ntr++ << " " <<  *wtr * (*itr)->Score( MOLS, IDS) << "  " << score << std::endl;
		}
		return score;
	}

	virtual void WriteScores( const std::vector< Molecule> &MOLS, const std::vector< std::vector< size_t> > &IDS, std::ostream &OUT) const
	{
		std::vector<double>::const_iterator
			wtr = m_Weights.begin();
		std::vector< std::string>::const_iterator
			ntr = m_Names.begin();
		double
			score,
			total = 0.0;
		for( std::vector< ScoreIF*>::const_iterator itr = m_Scores.begin(); itr != m_Scores.end(); ++itr, ++ntr, ++wtr)
		{
			score = (*itr)->GetCurrentScore();
			OUT << "# " << *ntr << " weight: " << *wtr << " raw: " << score << " weighted: " << *wtr * score << std::endl; // no recalculation!
			total += *wtr * score;
		}
		OUT << "# total: " << total << std::endl;
	}


	void Add( ScoreIF *SCORE, const double &WEIGHT = 1.0, const std::string &NAME = "")
	{
		m_Weights.push_back( WEIGHT);
		m_Names.push_back( NAME);
		m_Scores.push_back( SCORE);
	}


//	math::Vector3N Force() const;  // force from a pyramid of scores around thisposition? sum of vector forces between edges // compare to wuerfel // use center point as fifth point and all forces to it from edges (3,4,7,..)
//	std::pair< double, math::Vector3N> ScoreAndForce() const;

};

void
ReadScores( std::istream &IN, ScoreManager &MANAGER)
{
	std::cout << __FUNCTION__ << std::endl;
	std::string
		file,
		type;
	double
		weight;
	std::ifstream
		in;
	while( IN >> type >> weight)
	{
		std::cout << type << " " << weight << " ";
		if( type == "distc")
		{
			std::cout << std::endl;
			MANAGER.Add( ReadConstraint( IN), weight, type);
		}
		else if( type == "ff")
		{
			std::cout << std::endl;
			MANAGER.Add( &ForceField(), weight, type);
			flush( std::cout );
		}
		else if( type == "density")
		{
			IN >> file;
			std::cout << file << std::endl;
			MANAGER.Add( (ScoreIF *) new Density( file), weight, type);
			flush( std::cout );
		}
		else if( type == "aapair")
		{
			IN >> file;
			std::cout << file << std::endl;
			Open( in, file);
			MANAGER.Add( (ScoreIF *) new ScoreAAPair( in), weight, type);
			Close(in);
		}
		else if( type == "solv")
		{
			IN >> file;
			std::cout << file << std::endl;
			Open( in, file);
			MANAGER.Add( (ScoreIF *) new ScoreAASolvation( in), weight, type);
			Close(in);
		}
		else if( type == "bb")
		{
			IN >> file;
			std::cout << file << "  " << std::endl;
			Open( in, file);
			MANAGER.Add( (ScoreIF *) new ScorePhi( in), weight, type + "_phi");
			flush( std::cout );
			Close(in);

			IN >> file;
			std::cout << file << std::endl;
			Open( in, file);
			MANAGER.Add( (ScoreIF *) new ScorePsi( in), weight, type + "_psi");
			Close(in);
		}
		else
		{
			std::cerr << "ERROR: " << __FUNCTION__ << " no valid score: " << type << std::endl;
			exit(1);
		}
//		std::cout << std::endl;
	}

}


#endif /* SCORE_MANAGER_H_ */
