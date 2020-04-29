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



#ifndef SCORE_KB_H_
#define SCORE_KB_H_

#include "distribution.h"
#include "../include/std_functions.t.h"
#include "../include/string_functions.h"


class ScoreKB
: public ScoreIF
{
protected:
	std::map< std::string, Distribution>  m_Score;
	mutable double                        m_Current;


public:
	ScoreKB( std::istream &IN)
	: m_Score(),
	m_Current()
	{
		std::cout << __FUNCTION__ << " construct from stream" << std::endl;
		Read(IN);
	}

	ScoreKB()
	: m_Score()
	{}


	virtual ~ScoreKB(){}


	virtual void
	Read(std::istream &IN)
	{
		std::cout << __FUNCTION__ << std::endl;
		IN >> m_Score;
		std::cout << m_Score.size() << " entries" << std::endl;
		Patch();
	}


	void Patch()
	{
		std::cout << __FUNCTION__ << std::endl;

		std::vector< std::pair< std::string, Distribution> >
			pool;
		std::string
			newstr;
		for( std::map< std::string, Distribution>::iterator itr = m_Score.begin(); itr != m_Score.end(); ++itr)
		{
			if( itr->first.find( "HIS") != std::string::npos)
			{
				pool.push_back( *itr);
			}
		}

		for( std::vector< std::pair< std::string, Distribution> >::const_iterator itr = pool.begin(); itr != pool.end(); ++itr)
		{
			newstr = ReplaceAll( itr->first, "HIS", "HSE");
//			newstr.replace( newstr.find( "HIS"), 3 ,"HSE");
			std::cout << __FUNCTION__ << " " << itr->first << " => " << newstr << std::endl;
			m_Score[ newstr] = itr->second;
		}
		std::cout << m_Score.size() << " entries" << std::endl;
	}

	virtual const double &GetCurrentScore() const
	{
		return m_Current;
	}


	virtual double
	Score( const Residue &R1, const Residue &R2) const
	{
		std::cerr << "ERROR: base class function should not be called " << __PRETTY_FUNCTION__ << std::endl;
		exit(1);
		return std::numeric_limits<double>::max();
	}


	virtual double
	Score( const std::vector< Molecule> &MOLS) const // backbone definition for phi,psi
	{
		m_Current = 0.0;
		for( std::vector<Molecule>::const_iterator mol = MOLS.begin(); mol != MOLS.end(); ++mol)
			for( std::vector< Residue>::const_iterator res = mol->Residues().begin(); res + 1 != mol->Residues().end(); ++res)
			{
				m_Current += Score( *res, *(res+1));
			}
		return m_Current;
	}


	virtual double
	Score( const std::vector< Molecule> &MOLS, const std::vector< std::vector< size_t> > &IDS) const // backbone definition for phi,psi // TODO cleanup
	{
		m_Current = 0.0;
		int i = 0;
		for( std::vector<Molecule>::const_iterator mol = MOLS.begin(); mol != MOLS.end(); ++mol, ++i)
			if( IDS[i].size() > 0)
				for( std::vector< Residue>::const_iterator res = mol->Residues().begin() + IDS[i][0]; res + 1 != mol->Residues().end(); ++res)
				{
					m_Current += Score( *res, *(res+1));
				}
		return m_Current;
	}

};


class ScorePhi
: public ScoreKB
{
public:
	ScorePhi( std::istream &IN)
	: ScoreKB( IN)
	{}

	virtual ~ScorePhi(){}

	virtual double Score( const Residue &R1, const Residue &R2) const
	{
//		std::cout << "ScorePhi.Score(R1,R2): dihedral: " << math::Dihedral( R1.Atoms()[2].Position(), R2.Atoms()[0].Position(), R2.Atoms()[1].Position(), R2.Atoms()[2].Position()) << std::endl;
		return SafeGet( m_Score, SortAndGlue( R1.Name(), R2.Name())).InterpolatedValue( math::Dihedral( R1.Atoms()[2].Position(), R2.Atoms()[0].Position(), R2.Atoms()[1].Position(), R2.Atoms()[2].Position()));
	}
};


class ScorePsi
: public ScoreKB
{
public:
	ScorePsi( std::istream &IN)
	: ScoreKB( IN)
	{}

	virtual ~ScorePsi(){}

	virtual double Score( const Residue &R1, const Residue &R2) const
	{
//		std::cout << "ScorePsi.Score(R1,R2): dihedral: " <<  math::Dihedral( R1.Atoms()[0].Position(), R1.Atoms()[1].Position(), R1.Atoms()[2].Position(), R2.Atoms()[0].Position()) << std::endl;
		return SafeGet( m_Score, SortAndGlue( R1.Name(), R2.Name())).InterpolatedValue( math::Dihedral( R1.Atoms()[0].Position(), R1.Atoms()[1].Position(), R1.Atoms()[2].Position(), R2.Atoms()[0].Position()));
	}

};


double
ScoreBackbone( const Molecule &MOL, const size_t &ID, const ScorePhi &SCOREPHI, const ScorePsi &SCOREPSI)
{
	double
		energy = 0.0;

	for( size_t i = std::max( size_t(1), ID); i < MOL.Residues().size(); ++i)
	{
		energy += SCOREPHI.Score( MOL.Residues()[i-1], MOL.Residues()[i]);
	}
	for( size_t i = ID; i < MOL.Residues().size() - 1; ++i)
	{
		energy += SCOREPSI.Score( MOL.Residues()[i], MOL.Residues()[i+1]);
	}
	return energy;
}



class ScoreAAPair
: public ScoreKB
{
public:
	ScoreAAPair( std::istream &IN)
	: ScoreKB( IN)
	{}

	virtual ~ScoreAAPair(){}

	virtual double Score( const Residue &R1, const Residue &R2) const
	{
//		std::cout << "ScoreAAPair.Score(R1,R2) dist: " << math::Distance( R1.Atoms()[1].Position(), R2.Atoms()[1].Position())  << " " << R1.Name() << " " << R2.Name() << std::endl;
		return SafeGet( m_Score, SortAndGlue( R1.Name(), R2.Name())).InterpolatedValue( math::Distance( R1.Atoms()[1].Position(), R2.Atoms()[1].Position()));
	}


	virtual double
	Score( const std::vector< Molecule> &MOLS) const
	{
		m_Current = 0.0;
		for( std::vector<Molecule>::const_iterator mol1 = MOLS.begin(); mol1 != MOLS.end(); ++mol1)
			for( std::vector< Residue>::const_iterator res1 = mol1->Residues().begin(); res1 != mol1->Residues().end(); ++res1)
				for( std::vector<Molecule>::const_iterator mol2 = MOLS.begin(); mol2 != MOLS.end(); ++mol2)
					for( std::vector< Residue>::const_iterator res2 = mol2->Residues().begin(); res2 != mol2->Residues().end(); ++res2)
						if( ! ( mol1 == mol2 && abs( res1 - res2) < 3))
						{
							m_Current += Score( *res1, *res2);
						}
		return m_Current;
	}


	virtual double
	Score( const std::vector< Molecule> &MOLS, const std::vector< std::vector< size_t> > &IDS) const // backbone definition for phi,psi
	{
		m_Current = 0.0;
		for( size_t i = 0; i < IDS.size(); ++i)
			for( size_t j = 0; j < IDS[i].size(); ++j)
				for( size_t k = 0; k < MOLS.size(); ++k)
					for( size_t l = 0; l < MOLS[k].Residues().size(); ++l)
						if( ! ( i == k && abs( j - l) < 3))
						{
							m_Current += Score( MOLS[i].Residues()[j], MOLS[k].Residues()[l]);
						}
		return m_Current;
	}

};


class ScoreAASolvation
: public ScoreKB
{
public:
	ScoreAASolvation( std::istream &IN)
	: ScoreKB( IN)
	{}

	virtual ~ScoreAASolvation(){}


	virtual double Score( const Residue &R1, const Residue &R2) const
	{
		std::cout << "ERROR: ScoreAASolvation.Score(R1,R2) makes no sense!" << std::endl;
		exit(10);
	}

	virtual double Score( const std::vector< Molecule> &MOLS) const
	{
//		std::cout << "ScoreAASolvation.Score(MOLS)" << std::endl;
		std::vector< Triplet<char,double,char> >
			neighbor_counts = NeighborCount( MOLS, "all", std::vector< std::string>());
		m_Current = 0.0;
		for( std::vector< Triplet<char,double,char> >::const_iterator itr = neighbor_counts.begin(); itr != neighbor_counts.end(); ++itr)
		{
			m_Current += SafeGet( m_Score, Converter().ThreeLetter( itr->first)).Value( itr->second);
		}

		return m_Current;
	}

	virtual double
	Score( const std::vector< Molecule> &MOLS, const std::vector< std::vector< size_t> > &IDS) const // backbone definition for phi,psi
	{
//		std::cout << "ScoreAASolvation.Score(MOLS,IDS)" << std::endl;
		return m_Current  = Score( MOLS);
	}


};



#endif /* SCORE_KB_H_ */


