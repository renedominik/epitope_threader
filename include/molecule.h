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



#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <map>

#include "vector3N.h"
#include "math_functions.h"


class Atom
{
private:
	math::Vector3N      m_Pos;
	std::string         m_Name;
	std::string         m_Type;
//	std::string         m_ResName;
	double               m_Charge;
	double               m_VdwEpsilon;
	double               m_VdwRadius;
	double               m_Mass;
public:

	Atom()
	: m_Pos(),
	m_Name( "undef"),
	m_Type( "undef"),
	m_Charge(),
	m_VdwEpsilon(),
	m_VdwRadius(),
	m_Mass()
	{}

	Atom
	(
			const math::Vector3N &POS,
			const std::string &ATOM_NAME,
//			const std::string &RES_NAME,
			const std::string &ATOM_TYPE = "",
			double CHARGE = 0.0,
			double EPSILON = 0.0,
			double RADIUS = 0.0,
			double MASS = 0.0
	)
	: m_Pos( POS),
	m_Name( ATOM_NAME),
	m_Type( ATOM_TYPE),
//	m_ResName( RES_NAME),
	m_Charge( CHARGE),
	m_VdwEpsilon( EPSILON),
	m_VdwRadius( RADIUS),
	m_Mass( MASS)
	{}

	Atom( const Atom& ATOM)
	: m_Pos( ATOM.m_Pos),
	  m_Name( ATOM.m_Name),
	  m_Type( ATOM.m_Type),
	  m_Charge( ATOM.m_Charge),
	  m_VdwEpsilon( ATOM.m_VdwEpsilon),
	  m_VdwRadius( ATOM.m_VdwRadius),
	  m_Mass( ATOM.m_Mass)
	{}


	~Atom(){}

	math::Vector3N &Position(){ return m_Pos;}
	const math::Vector3N &Position() const{ return m_Pos;}

	std::string &Name(){ return m_Name;}
	const std::string &Name() const{ return m_Name;}

	std::string &Type(){ return m_Type;}

	double &Charge(){ return m_Charge;}
	const double &Charge() const{ return m_Charge;}

	double &Mass(){ return m_Mass;}
	const double &Mass() const{ return m_Mass;}

	double &VdwEpsilon(){ return m_VdwEpsilon;}
	const double &VdwEpsilon() const{ return m_VdwEpsilon;}

	double &VdwRadius(){ return m_VdwRadius;}
	const double &VdwRadius() const{ return m_VdwRadius;}

	void ReadFromPdbLine( const std::string &LINE);

	std::ostream &Write( std::ostream &OUT = std::cout) const;

	std::istream &Read( std::istream &IN);

};



class Residue
{
protected:
	std::string         m_Name;
	int					m_ID;
	std::vector< Atom>  m_Atoms;
//	double               m_Exposure;

public:
	Residue
	(
			const std::string &NAME = "",
			const std::vector< Atom> &ATOMS = std::vector< Atom>()
	)
	:m_Name( NAME),
	 m_ID(),
	m_Atoms( ATOMS)
//	m_Exposure( 0.0)
	{}

	Residue
	(
			const std::string &NAME,
			int RESID
	)
	:m_Name( NAME),
	 m_ID( RESID),
	m_Atoms()
//	m_Exposure( 0.0)
	{}

	int &ID()
	{
		return m_ID;
	}

	const int &ID() const
	{
		return m_ID;
	}


	std::string & Name()
	{
		return m_Name;
	}

	const std::string & Name() const
	{
		return m_Name;
	}

	std::vector< Atom>& Atoms()
	{
		return m_Atoms;
	}

	const std::vector< Atom>& Atoms() const
	{
		return m_Atoms;
	}

//	double &Exposure()
//	{
//		return m_Exposure;
//	}

//	const double &Exposure() const
//	{
//		return m_Exposure;
//	}

	std::ostream &Write( std::ostream &OUT = std::cout) const;

	std::istream &Read( std::istream &IN);

	void FuseToCalphaCentroid();

	void FuseToCbetaCentroid();

	void PatchNTerminus();

	void PatchCTerminus();

	void CentroidAt( const math::Vector3N &POS);

	Atom*
	FindCBeta() const;

	Atom*
	FindAtomByName( const std::string &NAME);

	const Atom*
	FindAtomByName( const std::string &NAME) const;

	Atom&
	FindAtomByType( const std::string &NAME);

	void
	BuildRandomSidechain();

	void
	ReplaceRandomRotamer( math::Vector3N *FIRST);

	void
	SortBackbone();

//	void Clear()
//	{
//		m_Name = "";
//		m_Atoms.clear();
//	}
};





class Molecule
{
protected:
	std::string              m_Name;
	char                     m_Chain;
	std::vector< Residue>    m_Residues;  //  amino acids
	math::Vector3N			 m_CMS;
public:
	Molecule( const std::string &NAME = "", char CHAIN = 'x', const std::vector< Residue> &RES = std::vector< Residue>())
	: m_Name( NAME),
	  m_Chain( CHAIN),
	  m_Residues( RES),
	  m_CMS()
	{
		if( m_Chain == ' ')
		{
			m_Chain = 'x';
		}
	}

	std::vector< Residue> & Residues(){ return m_Residues;}
	const std::vector< Residue> & Residues() const { return m_Residues;}

	std::string & Name(){ return m_Name;}

	char &Chain(){ return m_Chain;}
	const char &Chain() const{ return m_Chain;}

	size_t NrAtoms() const;

	void FuseToCalphaCentroid();

	void FuseToCbetaCentroid();

	std::ostream &Write( std::ostream &OUT = std::cout) const;

	std::istream &Read( std::istream &IN);

	math::Vector3N &
	CalcCMS();

	const math::Vector3N &
	GetCMS() const;

};


std::ostream &
Write( const std::vector< Molecule> &MOS, std::ostream &OUT = std::cout);

std::istream&
operator >>( std::istream &IN, Molecule &MOL);

std::ostream &
operator << ( std::ostream &OUT, const Molecule &MOL);


std::istream&
operator >>( std::istream &IN, Atom &A);

std::ostream &
operator << ( std::ostream &OUT, const Atom &A);

void
CalphaCentroids( std::vector< Molecule> &MOLS);

void
CbetaCentroids( std::vector< Molecule> &MOLS);

std::vector< Molecule>
ReadMolecules( std::istream &IN);

size_t
Chain2ID( const char CHAIN, const std::vector< Molecule> &MOLS);

size_t
NrAtoms( const std::vector< Molecule> &MOLS);


/*
std::vector< Molecule>
ReadMols( std::istream &IN);
*/


std::vector< Molecule>
ReadPdb( std::istream &IN, const std::string &MOL_NAME = "undefined");

void
WritePdb( const std::vector< Molecule> &MOLS, std::ostream &OUT);


void
WritePdb( const std::vector< std::vector< Atom> > &MOLS, std::ostream &OUT);

void
WritePdb( const Atom &ATOM, std::ostream &OUT, size_t ATOMID, size_t RESID, const std::string &RESNAME, char CHAIN);

void
SortResidueBackbone( std::vector< Molecule> &MOLS);

void
Exposure( std::istream &IN, const std::string &MODE, std::ostream &OUT, const std::vector< std::string> &CHAINS);

std::vector< Triplet<char,double,char> >
NeighborCount( const std::vector< Molecule> &MOLS, const std::string &MODE, const std::vector< std::string> &CHAINS);

double
ContactFunc( const Atom& A1, const Atom& A2, const std::string &MODE = "SINUS", double CONST_DIST = 7.0, double TRANSITION = 3.3);

void
WriteBFactor( std::istream &IN_PDB, std::ostream &OUT_PDB, std::istream &IN_FILE, size_t COL);

void
SSEDefinitions( const std::string &INPDB, const std::string &OUTPATH);

void
Profile( std::istream &IN_PDB, std::istream &IN_SCALE, std::ostream &OUT, const std::string &WINDOW_TYPE, const size_t &WINDOW_SIZE);


std::vector< char>
BuildAA();

bool
IsNaturalAminoAcid( const std::string &NAME);

std::vector< std::vector< Atom> >
Atoms( const std::vector< Molecule> &MOLS);

//double
//ForAllAtoms( const std::vector< Molecule> &MOLS, double (* FCT)(Atom));
//
//double
//ForAllAtoms( const std::vector< std::vector< Atom> > &MOLS, double (* FCT)(Atom));



class Converter
{
private:
	static std::map< char, std::string> m_1To3;
	static std::map< std::string, char> m_3To1;

public:
	Converter(){}
	~Converter(){}

	char SingleLetter( const std::string &AA);

	std::string ThreeLetter( char AA);

	bool Contains( char AA);
	bool Contains( const std::string &AA);

private:
	std::map< char, std::string> Build1To3();
	std::map< std::string, char> Build3To1();
};


void
SortPDBChainsByEpitopeLength( const std::string &PDB, size_t LENGTH, const std::string &OUTPATH);

math::Vector3N
CMS( const std::vector< Molecule> &MOLS);


double
InternalCADistanceDiff( const Molecule &MOL1, const Molecule &MOL2);


double
InternalCADistanceDiff( const std::vector< Molecule> &MOLS1, const std::vector< Molecule> &MOLS2);


/*
class
InternalDistanceDiff
{
private:
	std::vector< std::pair< std::pair< math::Vector3N*, math::Vector3N*>, std::pair< math::Vector3N*, math::Vector3N*> > >     m_Pos;
public:
	InternalDistanceDiff( const std::vector< Molecule> &MOLS1, const std::vector< Molecule> &MOLS2)
	{
		SetPositionPointers( MOLS1, MOLS2);
	}
	void SetPositionPointers( const std::vector< Molecule> &MOLS1, const std::vector< Molecule> &MOLS2)
	{
	}
	double operator()()
	{
		Histogram< double> hist( 0.0, 0.5, 100);
		for( std::vector< std::pair< std::pair< math::Vector3N*, math::Vector3N*>, std::pair< math::Vector3N*, math::Vector3N*> > >::const_iterator itr = m_Pos.begin(); itr != m_Pos.end(); ++itr)
		{
			hist( std::abs( math::Distance( itr->first.first, itr->first.second) - math::Distance( itr->second.first, itr->second.second)));
		}
	}
};
*/


inline
double
kT()
{
	return 0.593;
}


enum SSEType {e_Sheet, e_Helix, e_Random, e_Extended};

//std::vector< Residue>
//BuildFromTopologyFile( std::istream &IN);

#endif /* MOLECULE_H_ */
