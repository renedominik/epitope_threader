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



#ifndef FORCEFIELD_H_
#define FORCEFIELD_H_

//#include <set>
#include <map>

#include "molecule.h"
#include "scoreIF.h"


// class ForceFieldParameters stores all force field relevant data
// read from charmm topology file:
// residue:atomtype -> charge , this requires CorrectTermini()
// atomname -> atomtype
// read from charmm parameter file:
// atomtype -> mass, radii, epsilon

class ForceFieldParameters
{
private:
	std::map< std::string, double>   		m_TypeChargeMap;
	std::map< std::string, double>   		m_TypeEpsilonMap;
	std::map< std::string, double>   		m_TypeRadiusMap;
	std::map< std::string, double>   		m_TypeMassMap;
	std::map< std::string, std::string>  	m_NameTypeMap;

public:
	ForceFieldParameters()
	:m_TypeChargeMap(),
	 m_TypeEpsilonMap(),
	 m_TypeRadiusMap(),
	 m_TypeMassMap(),
	 m_NameTypeMap()
	{}

	~ForceFieldParameters(){}


	const double &Charge( const std::string &TYPE) const;  // needed?

	const double &Epsilon( const std::string &TYPE) const;

	const double &Radius( const std::string &TYPE) const;

	const double &Mass( const std::string &TYPE) const;

	void ReadPartialChargeMapFromCharmmTopologyFile( std::istream &STREAM);

	void ReadEpsilonAndRadiusFromCharmmParameterFile( std::istream &STREAM);

	void ReadMassFromCharmmTopologyFile( std::istream &STREAM);

	void ReadAtomTypeFromCharmmTopologyFile( std::istream& STREAM);

	Molecule &WriteAtomTypeIntoMolecule( Molecule &MOL) const;

	Molecule &WritePartialChargeIntoMolecule( Molecule &MOL) const;

	Molecule &WriteEpsilonAndRadiusIntoMolecule( Molecule &MOL) const;

	Molecule &WriteMassIntoMolecule( Molecule &MOL) const;

	Molecule &WriteAllIntoMolecule( Molecule &MOL) const;

	std::ostream &Write( std::ostream& OUT = std::cout) const;

//	const std::string &NameToType( const std::string &NAME) const;
};


ForceFieldParameters ReadForceFieldParameters( std::istream &IN_PAR, std::istream &IN_TOP);

void WriteForceFieldParametersIntoMolecules( const ForceFieldParameters &PAR, std::vector< Molecule> &MOLS);



template< typename T>
T CorrectTermini
(
		const std::map< std::string, T> &MAP,
		const std::string &RESIDUE,
		const std::string &ATOM
)
{
	std::string atom = ATOM;
	typename std::map< std::string, T>::const_iterator itr = MAP.find( RESIDUE + ":" + ATOM);
    if( itr == MAP.end())
    {
    	if( ATOM == "H1")
    	{
    		atom = "HT1";
    	}
    	else if( ATOM == "H2")
    	{
    		atom = "HT2";
    	}
    	else if( ATOM == "H3")
    	{
    		atom = "HT3";
    	}
    	else if( ATOM == "HN1")
    	{
    		atom = "HT1";
    	}
    	else if( ATOM == "HN2")
    	{
    		atom = "HT2";
    	}
    	else if( ATOM == "HN3")
    	{
    		atom = "HT3";
    	}

        itr = MAP.find( "NTER:" + atom);
        if( itr == MAP.end())
        {
            itr = MAP.find( "CTER:" + atom);
            if( itr == MAP.end())
            {
                std::cout << "ERROR: " << __FUNCTION__ << ": "<< RESIDUE << ":" << ATOM << "(" << atom << ") NOT found in map!" << std::endl;
                exit(1);
            }
        }
    }
    return itr->second;
}


class ForceField
: public ScoreIF
{
private:
	mutable double   m_Current;
public:
	ForceField()
	: m_Current()
	{}

	virtual ~ForceField(){}


	virtual const double &GetCurrentScore() const
	{
		return m_Current;
	}


	virtual
	double Score( const std::vector< Molecule> &MOLS) const;

	virtual
	double Score( const std::vector< Molecule> &MOLS, const std::vector< std::vector< size_t> > &SORTED_IDS) const;


	double Score( const std::vector< std::vector< Atom> > &MOLS) const;

	double Score( const Molecule &M1, const Molecule &M2) const;

	double Score( const Molecule &M) const;

	double Score( const Molecule &MOL, const size_t &ID) const;

	double Score( const Residue & RES, const std::vector< Molecule> &MOLS) const;

	double Score( const Residue &R1, const Residue &R2) const;

	double Score( const Atom &A1, const Atom &A2) const;

	double ScoreNextNeighbor( const Molecule &MOL, const size_t &ID) const;

};

double Coulomb( const Atom &A1, const Atom &A2);


double VanDerWaals( const Atom &A1, const Atom &A2);


double Coulomb( const Atom &A1, const Atom &A2, const double &DIST);


double VanDerWaals( const Atom &A1, const Atom &A2, const double &DIST);


void Pdb2mol
(
		std::istream &IN_PAR,
		std::istream &IN_TOP,
		std::istream &IN_PDB,
		const std::string &MOL_NAME,
		std::ostream &OUT
);


void DefaultAminoAcidsFromCharmm
(
		std::istream &IN_PAR,
		std::istream &IN_TOP,
		std::ostream &OUT
);


#endif /* FORCEFIELD_H_ */
