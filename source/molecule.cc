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

#include "../include/molecule.h"

#include <fstream>
#include <algorithm>
#include <stdio.h>

#include "../include/string_functions.h"
#include "../include/global.h"
#include "../include/std_functions.t.h"
#include "../include/math_functions.h"
#include "../include/io.h"
#include "../include/histogram.h"

void Atom::ReadFromPdbLine( const std::string &LINE)
{
	m_Type = "";
	m_Charge = 0.0;
	m_VdwEpsilon = 0.0;
	m_VdwRadius = 0.0;
	m_Mass = 0.0;
	m_Pos(0) = StringToValue<double>( LINE.substr( 30, 8));
	m_Pos(1) = StringToValue<double>( LINE.substr( 38, 8));
	m_Pos(2) = StringToValue<double>( LINE.substr( 46, 8));
	m_Name = Trim( LINE.substr( 12, 4));
}


std::ostream &Atom::Write( std::ostream &OUT) const
{
	OUT << m_Name << "   ";
	if( m_Type == "")
	{
		OUT << "x   ";
	}
	else
	{
		OUT << m_Type << "   ";
	}
	OUT << m_Pos[0] << "   " << m_Pos[1] << "   " << m_Pos[2];
	OUT << "   " << m_Charge << "   " << m_VdwEpsilon << "   " << m_VdwRadius << "   " << m_Mass << std::endl;

	return OUT;
}


std::istream &Atom::Read( std::istream &IN)
{
	IN >> m_Name >> m_Type >> m_Pos[0] >> m_Pos[1] >> m_Pos[2];
	IN >> m_Charge >> m_VdwEpsilon >> m_VdwRadius >> m_Mass;
	return IN;
}


std::ostream &Residue::Write( std::ostream &OUT) const
{
	OUT << m_Name << "  ";
	OUT << m_Atoms.size() << std::endl;
	for( unsigned int i = 0; i < m_Atoms.size(); ++i)
	{
		m_Atoms[i].Write( OUT);
	}
	return OUT;
}


std::istream &Residue::Read( std::istream &IN)
{
	unsigned int size;
	IN >> m_Name >> size;
	m_Atoms = std::vector< Atom>( size, Atom( math::Vector3N(), "undef"));
	for( size_t i = 0; i < size; ++i)
	{
		m_Atoms[i].Read( IN);
	}
	return IN;
}


void Residue::FuseToCalphaCentroid()
{
	math::Vector3N pos;
	std::vector< Atom>::iterator itr = m_Atoms.begin();
	bool found = false;
	while( itr != m_Atoms.end() && !found)
	{
		if(itr->Name() == "CA")
		{
			pos = itr->Position();
			found = true;
		}
		++itr;
	}
	for( itr = m_Atoms.begin(); itr != m_Atoms.end(); ++itr)
	{
		itr->Position() = pos;
	}
}


void Residue::FuseToCbetaCentroid()
{
	math::Vector3N pos;
	std::vector< Atom>::iterator itr = m_Atoms.begin();
	bool found = false;
	while( itr != m_Atoms.end() && !found)
	{
		if(m_Name != "GLY" && itr->Name() == "CB")
		{
			pos = itr->Position();
			found = true;
		}
		// place GLY centroid at center of HA1 and HA2
		else if( m_Name == "GLY" && (itr->Name() == "HA1" || itr->Name() == "HA2"))
		{
			pos += 0.5 * itr->Position();
		}
		++itr;
	}
	for( itr = m_Atoms.begin(); itr != m_Atoms.end(); ++itr)
	{
		itr->Position() = pos;
	}
}


void Residue::CentroidAt( const math::Vector3N &POS)
{
	for( std::vector< Atom>::iterator itr = m_Atoms.begin(); itr != m_Atoms.end(); ++itr)
	{
		itr->Position() = POS;
	}
}

Atom *
Residue::FindCBeta() const
{
	if( m_Name == "GLY")
	{
		for( std::vector< Atom>::const_iterator itr = m_Atoms.begin(); itr != m_Atoms.end(); ++itr)
			if( itr->Name() == "H1")
			{
//				std::cout << __FUNCTION__ << " h1 found" << std::endl;
				return (Atom*) &*itr;
			}
	}
	else
	{
		for( std::vector< Atom>::const_iterator itr = m_Atoms.begin(); itr != m_Atoms.end(); ++itr)
			if( itr->Name() == "CB")
			{
//				std::cout << __FUNCTION__ << " cb found" << std::endl;
				return (Atom*) &*itr;
			}
	}
		//	std::cerr << "CB not found in " << __FUNCTION__ << std::endl;
	return NULL;
//	return *std::find( m_Atoms.begin(), m_Atoms.end(), std::bind2nd( std::mem_fun_ref( &Atom::Name));
}



Atom *Residue::FindAtomByName( const std::string &NAME)
{
	for( std::vector< Atom>::iterator itr = m_Atoms.begin(); itr != m_Atoms.end(); ++itr)
		if( itr->Name() == NAME)
		{
			return (Atom*) &*itr;
		}
	std::cerr << "WARNING: " << NAME << " not found in " << __FUNCTION__ << " " << m_Name << " atoms " << m_Atoms << std::endl;
	return NULL;
//	return *std::find( m_Atoms.begin(), m_Atoms.end(), std::bind2nd( std::mem_fun_ref( &Atom::Name));
}


const Atom *Residue::FindAtomByName( const std::string &NAME) const
{
	for( std::vector< Atom>::const_iterator itr = m_Atoms.begin(); itr != m_Atoms.end(); ++itr)
		if( itr->Name() == NAME)
		{
			return (const Atom*) &*itr;
		}
	std::cerr << "WARNING: " << NAME << " not found in " << __FUNCTION__  << " " << m_Name << " atoms " << m_Atoms << std::endl;
	return NULL;
//	return *std::find( m_Atoms.begin(), m_Atoms.end(), std::bind2nd( std::mem_fun_ref( &Atom::Name));
}


Atom &Residue::FindAtomByType( const std::string &TYPE)
{
	for( std::vector< Atom>::iterator itr = m_Atoms.begin(); itr != m_Atoms.end(); ++itr)
		if( itr->Type() == TYPE)
		{
			return *itr;
		}
	std::cerr << TYPE << " not found in " << __FUNCTION__ << std::endl;
	exit(1);
	return m_Atoms[0];
//	return *std::find( m_Atoms.begin(), m_Atoms.end(), std::bind2nd( std::mem_fun_ref( &Atom::Name));
}



void
Residue::BuildRandomSidechain()
{

}


void
Residue::ReplaceRandomRotamer( math::Vector3N *FIRST)
{
	for( std::vector< Atom>::iterator itr = m_Atoms.begin(); itr != m_Atoms.end(); ++itr)
	{
		itr->Position() = *FIRST++;
	}
}

void
Residue::SortBackbone()
{
	std::vector< Atom>
		resorted(3);
	int
		count = 0;
	bool
		ca = false,c=false,n=false;
	for( std::vector< Atom>::const_iterator itr = m_Atoms.begin(); itr != m_Atoms.end(); ++itr)
	{
		if( itr->Name() == "N")
		{
			resorted[0] = *itr;
			n = true;
			++count;
		}
		else if( itr->Name() == "CA")
		{
			resorted[1] = *itr;
			ca = true;
			++count;
		}
		else if( itr->Name() == "C")
		{
			resorted[2] = *itr;
			c = true;
			++count;
		}
		else
		{
			resorted.push_back( *itr);
		}
	}
	if( m_Atoms.size() != resorted.size())
	{
		std::cout << "ERROR: " << __FUNCTION__ << " before:\n" << m_Atoms << "after:\n" << resorted;
		exit(1);
	}
	if( !ca || !c || !n)
	{
		std::cout << "ERROR: " << __FUNCTION__ << " not all backbone atoms found!" << std::endl;
		Write( std::cout);
		exit(1);
	}
	if( count != 3)
	{
		std::cout << "ERROR: " << __FUNCTION__ << " incorrect number of backbone atoms found " << count << std::endl;
		Write();
		exit(1);
	}
	m_Atoms = resorted;
	assert( m_Atoms[0].Name() == "N" && m_Atoms[1].Name() == "CA" && m_Atoms[2].Name() == "C");
}


std::ostream &Molecule::Write( std::ostream &OUT) const
{
	OUT << m_Name << "  " << m_Chain << "  ";
	OUT << m_Residues.size() << std::endl;
	for( size_t i = 0; i < m_Residues.size(); ++i)
	{
		m_Residues[i].Write( OUT);
	}
	return OUT;
}


size_t Molecule::NrAtoms() const
{
	size_t nr = 0;
	for( std::vector< Residue>::const_iterator itr = m_Residues.begin(); itr != m_Residues.end(); ++itr)
	{
		nr += itr->Atoms().size();
	}
	return nr;
}


void Molecule::FuseToCalphaCentroid()
{
	DebugFct;
	std::for_each( m_Residues.begin(), m_Residues.end(), std::mem_fun_ref( &Residue::FuseToCalphaCentroid));
}


void Molecule::FuseToCbetaCentroid()
{
	DebugFct;
	std::for_each( m_Residues.begin(), m_Residues.end(), std::mem_fun_ref( &Residue::FuseToCbetaCentroid));
}

math::Vector3N &
Molecule::CalcCMS()
{
	math::Vector3N
		cms;
	double
		mass = 0;
	for( std::vector< Residue>::const_iterator ritr = m_Residues.begin(); ritr != m_Residues.end(); ++ritr)
		for( std::vector< Atom>::const_iterator aitr = ritr->Atoms().begin(); aitr != ritr->Atoms().end(); ++aitr)
		{
			cms += aitr->Mass() * aitr->Position();
			mass += aitr->Mass();
		}
	m_CMS = cms / mass;
	return m_CMS;
}

const math::Vector3N &
Molecule::GetCMS() const
{
	return m_CMS;
}



std::ostream &Write( const std::vector< Molecule> &MOLS, std::ostream &OUT)
{
	DebugFct;
	OUT << MOLS.size() << std::endl;
	for( size_t i = 0; i < MOLS.size(); ++i)
	{
		MOLS[i].Write( OUT);
	}
	return OUT;
}


std::istream &Molecule::Read( std::istream &IN)
{
	int size;
	IN >> m_Name >> m_Chain >> size;
	m_Residues = std::vector< Residue>( size);
	for( int i = 0; i < size; ++i)
	{
		m_Residues[i].Read( IN);
	}
	return IN;
}


std::istream& operator >>( std::istream &IN, Molecule &MOL)
{
	return MOL.Read(IN);
}

std::ostream & operator << ( std::ostream &OUT, const Molecule &MOL)
{
	return MOL.Write( OUT);
}

std::istream&
operator >>( std::istream &IN, Atom &A)
{
	return A.Read( IN);
}

std::ostream &
operator << ( std::ostream &OUT, const Atom &A)
{
	return A.Write( OUT);
}



void CalphaCentroids( std::vector< Molecule> &MOLS)
{
	DebugFct;
	std::for_each( MOLS.begin(), MOLS.end(), std::mem_fun_ref( &Molecule::FuseToCalphaCentroid));
}

void CbetaCentroids( std::vector< Molecule> &MOLS)
{
	DebugFct;
	std::for_each( MOLS.begin(), MOLS.end(), std::mem_fun_ref( &Molecule::FuseToCbetaCentroid));
}



std::vector< Molecule> ReadMolecules( std::istream &IN)
{
	DebugFct;
	size_t size;
	IN >> size;
	std::vector< Molecule> mols( size);
	for( size_t i = 0; i < size; ++i)
	{
		mols[i].Read(IN);
	}
	return mols;
}


size_t Chain2ID( const char CHAIN, const std::vector< Molecule> &MOLS)
{
	size_t count = 0;
	while( count < MOLS.size())
	{
		if( MOLS[count].Chain() == CHAIN)
		{
			break;
		}
		++count;
	}
	if( count >= MOLS.size())
	{
		std::cerr << "ERROR: epitope chain <" << CHAIN << "> not found" << std::endl;
		exit(1);
	}
	return count;
}


size_t NrAtoms( const std::vector< Molecule> &MOLS)
{
	size_t nr = 0;
	for( std::vector< Molecule>::const_iterator itr = MOLS.begin(); itr != MOLS.end(); ++itr)
	{
		nr += itr->NrAtoms();
	}
	return nr;
}


/*
std::vector< Molecule>
ReadMols( std::istream &IN)
{
	int size;
	IN >> size;
	std::vector< Molecule> mols( size);
	for( int i = 0; i < size; ++i)
	{
		mols[i].Read( IN);
	}
	return mols;
}
*/


std::vector< Molecule>
ReadPdb( std::istream &IN, const std::string &MOL_NAME)
{
	DebugFct;
	std::string
		line,
		residue_type,
		previous_type;
	char
		alternative,
		previous_chain = 'x',
		chain;
	int
		residue_id,
		previous_id( std::numeric_limits< int>::max());

	Molecule
		molecule;
	Residue
		residue;
	Atom
		atom;
	std::vector< Molecule>
		all_molecules;

	// collecting and fitting limits
	while( IN)
	{
		line.clear();
		std::getline( IN, line);
		if( line.length() > 6 && line.substr( 0, 4) == "ATOM")
		{
			alternative = line[16];
			if( alternative != ' ' && alternative != 'A')
			{
//				std::cout << "WARNING: line is ignored, alternative location is not empty nor 'A': " << line << std::endl;
				continue;
			}

			atom.ReadFromPdbLine( line);

			chain = line[ 21];
			residue_id = StringToValue<int>( line.substr( 22, 4));
			residue_type = Trim( line.substr( 17, 3));

			// ignore not natural AAs
			if( residue_type.size() < 3)
			{
				continue;
			}

			// new residue ?
			if( residue_id != previous_id || residue_type != previous_type)
			{
				previous_id = residue_id;
				previous_type = residue_type;

				if( residue.Atoms().size() > 0)
				{
					molecule.Residues().push_back( residue);
				}
				residue = Residue( residue_type, residue_id);
			}

			// new molecule = new chain ?
			if( chain != previous_chain)
			{
//				std::cout << "new chain" << std::endl;
				previous_chain = chain;
				if( molecule.Residues().size() > 0)
				{
					all_molecules.push_back( molecule);
				}
				molecule = Molecule( MOL_NAME, chain);
			}
			// add atom to residue
			residue.Atoms().push_back( atom);
		}
		else if( line.length() > 12 && line.substr( 0, 5) == "MODEL")
		{
			std::vector< std::string> cols = Split( line);
			if( cols[1] != "1")
			{
				std::cout << __FUNCTION__ << " break reading at: " << line << std::endl;
				break;
			}
		}
	}
	if( residue.Atoms().size() > 0)
	{
		molecule.Residues().push_back( residue);
	}
	if( molecule.Residues().size() > 0)
	{
		all_molecules.push_back( molecule);
	}
	return all_molecules;
}


void
WritePdb( const std::vector< Molecule> &MOLS, std::ostream &OUT)
{
	size_t
		atom_count = 1,
		residue_count = 1;
	std::string
		resname;
	char
		chain;

	for( std::vector< Molecule>::const_iterator mol = MOLS.begin(); mol != MOLS.end(); ++mol)
	{
		residue_count = 1;
		chain = mol->Chain();
		for( std::vector< Residue>::const_iterator aa = mol->Residues().begin(); aa != mol->Residues().end(); ++aa, ++residue_count)
		{
			resname = aa->Name();
			for( std::vector< Atom>::const_iterator atr = aa->Atoms().begin(); atr != aa->Atoms().end(); ++atr, ++atom_count)
			{
				WritePdb( *atr, OUT, atom_count, residue_count, resname, chain);
			}
		}
	}
}

void
WritePdb( const std::vector< std::vector< Atom> > &MOLS, std::ostream &OUT)
{
	size_t
		residue_count,
		mol_count = 0,
		atom_count = 1;

	for( std::vector< std::vector< Atom> >::const_iterator mol = MOLS.begin(); mol != MOLS.end(); ++mol, ++mol_count)
	{
		residue_count = 0;
		for( std::vector< Atom>::const_iterator atr = mol->begin(); atr != mol->end(); ++atr, ++atom_count)
		{
			if( atr->Name() == "C" || atr->Name() == "CA" || atr->Name() == "N")
			{
				if( atr->Name() == "N")
				{
					residue_count++;
				}
				WritePdb( *atr, OUT, atom_count, residue_count, "GLY", char( 97 + mol_count));
			}
		}
	}
}




void
WritePdb( const Atom &ATOM, std::ostream &OUT, size_t ATOMID, size_t RESID, const std::string &RESNAME, char CHAIN)
{
	char buffer[90];

	math::Vector3N pos = ATOM.Position();
	std::string name;
	double bfact = 0;
	double charge = 0.0;

	sprintf (buffer, "ATOM  %5u %4s %3s %1c%4u    %8.3f%8.3f%8.3f      %6.2f          %4f\n", (unsigned int) ATOMID % 100000, ATOM.Name().c_str(), RESNAME.c_str(), CHAIN, (unsigned int) RESID % 10000, pos[0], pos[1], pos[2], bfact, charge);

	OUT << buffer;
}


void
SortResidueBackbone( std::vector< Molecule> &MOLS)
{
	for( std::vector< Molecule>::iterator mol = MOLS.begin(); mol != MOLS.end(); ++mol)
		for( std::vector< Residue>::iterator aa = mol->Residues().begin(); aa != mol->Residues().end(); ++aa)
		{
			aa->SortBackbone();
		}
}



void
Exposure( std::istream &IN, const std::string &MODE, std::ostream &OUT, const std::vector< std::string> &CHAINS)
{
	std::vector< Triplet<char,double,char> >
		vec = NeighborCount( ReadPdb( IN), MODE, CHAINS);
	for( std::vector< Triplet<char,double,char> >::const_iterator itr = vec.begin(); itr != vec.end(); ++itr)
	{
		OUT << itr->first << "  " << itr->second << "  " << itr->third << std::endl;
	}
//	WriteVector( NeighborCount( ReadPdb( IN), MODE), OUT, '\n');
}


std::vector< Triplet<char,double,char> >
NeighborCount( const std::vector< Molecule> &MOLS, const std::string &MODE, const std::vector< std::string> &CHAINS)
{
	std::vector< Triplet<char,double,char> > counts;
	Atom *ca1, *ca2;
	double contact;

	std::vector< char>
		first_block_of_chains,
		second_block_of_chains;
	// if chain selections are defined
	if( CHAINS.size() == 2)
	{
		for( size_t i = 0; i < CHAINS[0].size(); i+=2)
		{
			first_block_of_chains.push_back( CHAINS[0][i]);
		}
		for( size_t i = 0; i < CHAINS[1].size(); i+=2)
		{
			second_block_of_chains.push_back( CHAINS[1][i]);
		}
		std::cout << __FUNCTION__ << " chains in each block: \n" << first_block_of_chains << "\n" << second_block_of_chains << std::endl;
	}

	for( std::vector< Molecule>::const_iterator m1 = MOLS.begin(); m1 != MOLS.end(); ++m1)
	{
		for( std::vector< Residue>::const_iterator r1 = m1->Residues().begin(); r1 != m1->Residues().end(); ++r1)
		{
			double count = 0.0;
			ca1 = (Atom*) r1->FindAtomByName( "CA");
			if( ca1 == NULL)
			{
				continue;
			}
			if( MODE == "all")
			{
				for( std::vector< Molecule>::const_iterator m2 = MOLS.begin(); m2 != MOLS.end(); ++m2)
				{
					for( std::vector< Residue>::const_iterator r2 = m2->Residues().begin(); r2 != m2->Residues().end(); ++r2)
					{
						if( r1 == r2) continue;

						ca2 = (Atom*) r2->FindAtomByName( "CA"); // TODO: take this out of the loop (pre run)
						if( ca2 != NULL)
						{
							contact = ContactFunc( *ca1, *ca2, "SINUS", 4.0, 7.4);  // see Durham, ..,Staritzbichler, Meiler 2009 J Mol Model  // restype independent cutoffs
//						((Residue) *r1).Exposure() += contact;
							count += contact;
//						((Residue) *r2).Exposure() += contact;
						}
					}
				}
			}
			else if( MODE == "inter")
			{
				for( std::vector< Molecule>::const_iterator m2 = MOLS.begin(); m2 != MOLS.end(); ++m2)
				{
					if( m1 == m2
							|| ( CHAINS.size() == 2
									&& ! ( ( std::find( first_block_of_chains.begin(), first_block_of_chains.end(), m1->Chain()) != first_block_of_chains.end()
													&& std::find( second_block_of_chains.begin(), second_block_of_chains.end(), m2->Chain()) != second_block_of_chains.end())
											|| 	( std::find( first_block_of_chains.begin(), first_block_of_chains.end(), m2->Chain()) != first_block_of_chains.end()
													&& std::find( second_block_of_chains.begin(), second_block_of_chains.end(), m1->Chain()) != second_block_of_chains.end()))))
					{ /*std::cout << "ignore " << m1->Chain() << " " << m2->Chain() << std::endl;*/ continue;}
//					else
//					{ std::cout << "go " << m1->Chain() << " " << m2->Chain() << std::endl;}

					for( std::vector< Residue>::const_iterator r2 = m2->Residues().begin(); r2 != m2->Residues().end(); ++r2)
					{
						ca2 = (Atom*) r2->FindAtomByName( "CA"); // TODO: take this out of the loop (pre run)
						if( ca2 != NULL)
						{
							count += ContactFunc( *ca1, *ca2, "SINUS", 4.0, 7.4);  // see Durham, ..,Staritzbichler, Meiler 2009 J Mol Model  // restype independent cutoffs
//						((Residue) *r1).Exposure() += contact;
//						((Residue) *r2).Exposure() += contact;
						}
					}
				}
			}
			else if( MODE == "intra")
			{
				for( std::vector< Residue>::const_iterator r2 = m1->Residues().begin(); r2 != m1->Residues().end(); ++r2)
				{
					if( r1 == r2) continue;

					ca2 = (Atom*) r2->FindAtomByName( "CA"); // TODO: take this out of the loop (pre run)
					if( ca2 != NULL)
					{
						count +=  ContactFunc( *ca1, *ca2, "SINUS", 4.0, 7.4);  // see Durham, ..,Staritzbichler, Meiler 2009 J Mol Model  // restype independent cutoffs
//					((Residue) *r1).Exposure() += contact;
//						((Residue) *r2).Exposure() += contact;
					}
				}
			}
			char c = Converter().SingleLetter( r1->Name());
			counts.push_back( Triplet<char,double,char>( c, count, m1->Chain()));
		}
	}

//	for( std::vector< Molecule>::const_iterator m1 = MOLS.begin(); m1 != MOLS.end(); ++m1)
//	{
//		for( std::vector< Residue>::const_iterator r1 = m1->Residues().begin(); r1 != m1->Residues().end(); ++r1)
//		{
//			counts.push_back( ((Residue) *r1).Exposure());
//		}
//	}

	return counts;
}



void
WriteBFactor( std::istream &IN_PDB, std::ostream &OUT_PDB, std::istream &IN_FILE, size_t COL)
{
	std::string
		line,
		val_str,
		resname,
		prev_resname = "zzzzzzz";
	std::vector< std::string>
		cols,
		vals;
	size_t
		resid, prev_resid = 9999999999;
	int
		i = -1;

	std::vector< char>
		names;

	while( IN_FILE)
	{
		line.clear();
		std::getline( IN_FILE, line);
		line = Trim( line);
		cols = Split( line);
		if( cols.size() < COL + 1)
		{
			std::cout << "WARNING: line does not contain columns and is ignored <" << line << ">" << std::endl;
		}
		else
		{
			val_str = ValueToString( math::Round( StringToValue< double>( cols[ COL]), 2));
			if( val_str.size() < 6)
			{
//				val_str.append( 6 - val_str.size(), ' ');
				val_str = std::string( 6 - val_str.size(), ' ') + val_str;
			}
			else if( val_str.size() > 6)
			{
				std::cout << "WARNING: value is too large to be placed into PDB file: " << val_str << std::endl;
			}
			vals.push_back( val_str);
			names.push_back( cols[0][0]);
		}
	}

	std::cout << "loaded " << vals.size() << " values" << std::endl;

	while( IN_PDB)
	{
		line.clear();
		std::getline( IN_PDB, line);

//		if( line.size() == 0 || (line.substr(0, 4) != "ATOM" && line.substr( 0, 6) != "HETATM")) continue;
		if( line.size() == 0 || line.substr(0, 4) != "ATOM") continue;
//		std::cout << line << "|" << line.size() << std::endl;
		resid = StringToValue< size_t>( line.substr( 22, 4)); // actually not necessary, can be checked as strings
		resname = line.substr( 17,3);

		if( resid != prev_resid || resname != prev_resname)
		{
			++i;
			prev_resid = resid;
			prev_resname = resname;
			val_str = vals[i];

			if( names[i] != Converter().SingleLetter( line.substr( 17, 3)))
			{
				std::cout << "WARNING: " <<  Converter().ThreeLetter( names[i]) << " >> " << line << std::endl;
			}
		}
//		std::cout << line << "|" << line.size() << std::endl;
		if (line.size() < 80)
		{
			line.append( 80 - line.size(), ' '); // fill up line
		}
//		std::cout << line << "|" << line.size() << std::endl;
		line.replace ( 60, 6, val_str);
//		std::cout << line << "|" << line.size() << std::endl;
		OUT_PDB << line << std::endl;
	}


}





double ContactFunc( const Atom& A1, const Atom& A2, const std::string &MODE, double CONST_DIST, double TRANSITION) // TODO: do separate functions instead of mode-if for speed up
{
	double dist = Distance( A1.Position(), A2.Position());
	if( dist <= CONST_DIST)
	{
		return 1.0;
	}
	if( dist < CONST_DIST + TRANSITION)
	{
		if( MODE == "SINUS")
		{
			return 0.5 + 0.5 * cos( M_PI / TRANSITION * ( dist - CONST_DIST));
		}
		if( MODE == "LINEAR")
		{
			return 1.0 - (dist - CONST_DIST) / TRANSITION;
		}
	}
	return 0.0;
}

struct SSE
{
	char        m_chain;
	int         m_first_id;
	std::string m_first_resname;
	int         m_last_id;
	std::string m_last_resname;

	SSE( char CHAIN, int FIRSTID, const std::string &FIRST_RESNAME, int SECONDID, const std::string &SECOND_RESNAME)
	: m_chain( CHAIN),
	  m_first_id( FIRSTID),
	  m_first_resname( FIRST_RESNAME),
	  m_last_id( SECONDID),
	  m_last_resname( SECOND_RESNAME)
	{}
};


bool
IsWithinSSE( char CHAIN, int RESID, const std::string &NAME, const std::vector< SSE> &V)
{
	bool found = false;

	for( std::vector< SSE>::const_iterator itr = V.begin(); itr != V.end() ; ++itr)
	{
//		std::cout << CHAIN << " " << itr->m_chain << " " << RESID << " "  << itr->m_first_id << " " << itr->m_last_id << " " << NAME << " " << itr->m_first_resname << " " << itr->m_last_resname << std::endl;
		if( itr->m_chain == CHAIN && RESID >= itr->m_first_id && RESID <= itr->m_last_id)
		{
//			std::cout << "found " << std::endl;
			found = true;
			// security
			if( RESID == itr->m_first_id && NAME != itr->m_first_resname)
			{
				std::cout << "WARNING: residue names of first aa in sse do not match in " << __FUNCTION__ << ": " << CHAIN << " " << RESID << ":" << NAME << " != " << itr->m_first_resname << std::endl;
			}
			if( RESID == itr->m_last_id && NAME != itr->m_last_resname)
			{
				std::cout << "WARNING: residue names of first aa in sse do not match in " << __FUNCTION__ << ": " << CHAIN << " " << RESID << ":" << NAME << " != " << itr->m_last_resname << std::endl;
			}

		}

	}

	return found;
}


void
SSEDefinitions( const std::string &INPDB, const std::string &OUTPATH)
{
	std::string
		line;
	std::vector< SSE>
		sheets,
		helices;
	std::vector<Molecule>
		mol;
	char
		chain;
	bool
		is_helix,
		is_sheet;

	size_t
		first = INPDB.find_last_of( '/'),
		first_id,
		last_id;
	if( first == std::string::npos)
	{
		first = 0;
	}
	else
	{
		++first;
	}

	std::string
		base = INPDB.substr( first, INPDB.size() - first - 4);

	std::ifstream
		in;
	std::ofstream
		out;
	Open( in, INPDB);
	Open( out, OUTPATH + "/" + base + "_sse.txt");



	while( in)
	{
		line.clear();
		std::getline( in, line);
		if( line.size() < 6) continue;
		if( line.substr(0, 5) == "HELIX")
		{
			first_id = StringToValue<int>( line.substr( 21,4));
			last_id = StringToValue<size_t>( line.substr( 33,4));
			chain = line[19];
			if( chain != line[31])
			{
				std::cout << "WARNING: in " << __FUNCTION__ << ": chainID of initial and last residue in HELIX does not match!: " << chain << " " << line[31] << std::endl;
			}
			helices.push_back( SSE( chain, first_id, line.substr( 15,3), last_id, line.substr( 27,3)));
		}
		else if (line.substr(0, 5) == "SHEET")
		{
			first_id = StringToValue<size_t>( line.substr( 22,4));
			last_id = StringToValue<size_t>( line.substr( 33,4));
			chain = line[21];
			if( chain != line[32])
			{
				std::cout << "WARNING: in " << __FUNCTION__ << ": chainID of initial and last residue in SHEET does not match!: " << line << std::endl;
			}
			sheets.push_back( SSE( chain, first_id, line.substr( 17,3), last_id, line.substr( 28,3)));
		}
	}

	std::cout << "found " << helices.size() << " helices and " << sheets.size() << " sheets" << std::endl;

	in.clear();
	in.seekg(0, std::ios::beg);
	// read sequence information
	mol = ReadPdb( in);
	Close( in);

//	// clean up helices and sheets
//	CleanSSE( helices);
//	CleanSSE( sheets);


	// write helices and sheets
	for( std::vector< Molecule>::const_iterator mtr = mol.begin(); mtr != mol.end(); ++mtr)
		for( std::vector< Residue>::const_iterator rtr = mtr->Residues().begin(); rtr != mtr->Residues().end(); ++rtr)
		{
			is_helix = IsWithinSSE( mtr->Chain(), rtr->ID(), rtr->Name(), helices);
			is_sheet = IsWithinSSE( mtr->Chain(), rtr->ID(), rtr->Name(), sheets);
			out << Converter().SingleLetter( rtr->Name()) << "   " << is_helix << "  " << is_sheet << "  " << mtr->Chain() << std::endl;
		}


	Close( out);
	Close( out);
}


std::map< char, double>
ReadScale( std::istream &IN)
{
	std::string
		key;
	double
		value;
	std::map< char, double>
		scale;
	while (IN >> key >> value)
	{
		key = key.substr(0,3);
		std::transform( key.begin(), key.end(), key.begin(), ::toupper);
		scale[ Converter().SingleLetter( key)] = value;
	}
	return scale;
}


void
Profile( std::istream &IN_PDB, std::istream &IN_SCALE, std::ostream &OUT, const std::string &WINDOW_TYPE, const size_t &WINDOW_SIZE)
{
	std::vector< Molecule>
		mol = ReadPdb( IN_PDB);
	std::map< char, double>
		scale = ReadScale( IN_SCALE);
	char
		name;

	if( WINDOW_TYPE == "none")
	{
		for( std::vector< Molecule>::const_iterator mtr = mol.begin(); mtr != mol.end(); ++mtr)
			for( std::vector< Residue>::const_iterator rtr = mtr->Residues().begin(); rtr != mtr->Residues().end(); ++rtr)
			{
				name = Converter().SingleLetter( rtr->Name());
				OUT << name << "  " << SafeGet( scale, name) << std::endl;
			}
		OUT << std::endl;
	}
	else if( WINDOW_TYPE == "basic")
	{
		size_t
			min,
			max,
			count,
			k,
			width = WINDOW_SIZE / 2;
		double
			sum;

		for( size_t i = 0; i < mol.size(); ++i)
			for( size_t j = 0; j < mol[i].Residues().size(); ++j)
			{
				sum = 0.0;
				count = 0;
				min = std::max( (size_t) 0, j - width);
				max = std::min( mol[i].Residues().size() - 1, j + width);
				for( k = min; k <= max; ++k, ++count)  // could be improved
				{
					name = Converter().SingleLetter( mol[i].Residues()[k].Name());
					sum += SafeGet( scale, name);
//					std::cout << i << " " << j << " " << min << " " << max << " " << k << " " << name << " " << sum << std::endl;
				}
				sum /= (double) count;
				name = Converter().SingleLetter( mol[i].Residues()[j].Name());
//				std::cout << name << " " << sum << std::endl;
				OUT << name << "  " << sum << std::endl;
			}
	}
	else if( WINDOW_TYPE == "triangular")
	{
		std::cout << "triangular window is NOT IMPLEMENTED YET" << std::endl;
	}
}





std::vector< char>
BuildAA()
{
	DebugFct;
	std::vector< char>
		aa( 24);

	aa[0] = 'A';
	aa[1] = 'R';
	aa[2] = 'N';
	aa[3] = 'D';
	aa[4] = 'C';
	aa[5] = 'Q';
	aa[6] = 'E';
	aa[7] = 'G';
	aa[8] = 'H';
	aa[9] = 'I';
	aa[10] = 'L';
	aa[11] = 'K';
	aa[12] = 'M';
	aa[13] = 'F';
	aa[14] = 'P';
	aa[15] = 'S';
	aa[16] = 'T';
	aa[17] = 'W';
	aa[18] = 'Y';
	aa[19] = 'V';
	aa[20] = 'B';
	aa[21] = 'Z';
	aa[22] = '-';
	aa[23] = 'U';

	return aa;
}



std::map< char, std::string> Converter::m_1To3 = Converter().Build1To3();
std::map< std::string, char> Converter::m_3To1 = Converter().Build3To1();

char Converter::SingleLetter( const std::string &AA)
{
	return SafeGet( m_3To1, AA);
}

std::string Converter::ThreeLetter( char AA)
{
	return SafeGet( m_1To3, AA);
}

std::map< char, std::string> Converter::Build1To3()
{
	DebugFct;
	std::map< char, std::string> map;

	map['A'] = ToUpper("Ala");
	map['R'] = ToUpper("Arg");
	map['N'] = ToUpper("Asn");
	map['D'] = ToUpper("Asp");
	map['C'] = ToUpper("Cys");
	map['E'] = ToUpper("Glu");
	map['Q'] = ToUpper("Gln");
	map['G'] = ToUpper("Gly");
	map['H'] = ToUpper("His");
	map['I'] = ToUpper("Ile");
	map['L'] = ToUpper("Leu");
	map['K'] = ToUpper("Lys");
	map['M'] = ToUpper("Met");
	map['F'] = ToUpper("Phe");
	map['P'] = ToUpper("Pro");
	map['S'] = ToUpper("Ser");
	map['T'] = ToUpper("Thr");
	map['W'] = ToUpper("Trp");
	map['Y'] = ToUpper("Tyr");
	map['V'] = ToUpper("Val");
	map['U'] = ToUpper("Unk");
	map['Z'] = ToUpper("Glx");
	map['B'] = ToUpper("Asx");

	return map;
}


std::map< std::string, char> Converter::Build3To1()
{
	DebugFct;
	std::map< std::string, char> map;

	map[ToUpper("Ala")] = 'A';
	map[ToUpper("Arg")] = 'R';
	map[ToUpper("Asn")] = 'N';
	map[ToUpper("Asp")] = 'D';
	map[ToUpper("Cys")] = 'C';
	map[ToUpper("Glu")] = 'E';
	map[ToUpper("Gln")] = 'Q';
	map[ToUpper("Gly")] = 'G';
	map[ToUpper("His")] = 'H';
	map[ToUpper("HSE")] = 'H';  // TODO: DOUBLE CHECK!!
	map[ToUpper("Ile")] = 'I';
	map[ToUpper("Leu")] = 'L';
	map[ToUpper("Lys")] = 'K';
	map[ToUpper("Met")] = 'M';
	map[ToUpper("Phe")] = 'F';
	map[ToUpper("Pro")] = 'P';
	map[ToUpper("Ser")] = 'S';
	map[ToUpper("Thr")] = 'T';
	map[ToUpper("Trp")] = 'W';
	map[ToUpper("Tyr")] = 'Y';
	map[ToUpper("Val")] = 'V';
	map[ToUpper("Unk")] = 'U';
	map[ToUpper("Glx")] = 'Z';
	map[ToUpper("Asx")] = 'B';

	return map;
}

bool Converter::Contains( char AA)
{
//	DebugFct;
	return m_1To3.find( AA) != m_1To3.end();
}
bool Converter::Contains( const std::string &AA)
{
//	DebugFct;
//	std::cout << AA << " " << ( m_3To1.find( AA) != m_3To1.end() ? "TRUE" : "FALSE") << "  " << m_3To1.size() <<  std::endl;
	return m_3To1.find( AA) != m_3To1.end();
}


bool IsNaturalAminoAcid( const std::string &NAME)
{
//	DebugFct;
	return Converter().Contains( NAME);
}



void SortPDBChainsByEpitopeLength( const std::string &PDB, size_t LENGTH, const std::string &OUTPATH)
{
	std::vector< std::vector< std::string> >
		chains;
	std::vector< std::string>
		block;
	std::vector< size_t>
		counters;
	std::string
		line,
		pdb_id = PDB;
	char
		previous_chain = 'x',
		chain;
	size_t
		residue_count = 0,
		resid,
		hits = 0,
		hit = std::numeric_limits< size_t>::max(),
		i_letter = 66; // = 'B'
	int
		prev_resid = -99999;
	std::ifstream
		in;
	std::ofstream
		out;

	size_t found = PDB.find( '/');
	if( found != std::string::npos)
	{
		pdb_id = PDB.substr( found + 1);
	}

	if( pdb_id.find( '.') != std::string::npos)
	{
		pdb_id = pdb_id.substr(0, pdb_id.find('.'));
	}


	Open( in, PDB);

	while( in)
	{
		line.clear();
		getline( in, line);
		if( line.size() > 50 && line.substr(0,4) == "ATOM")
		{
			chain = line[21];
			if( chain != previous_chain)
			{
				if( block.size() > 0)
				{
					chains.push_back( block);
					counters.push_back( residue_count);
					residue_count = 0;
					block.clear();
					prev_resid = -99999;
				}
				previous_chain = chain;
			}

			resid = StringToValue< size_t>( line.substr( 22, 4));
			if( (signed) resid != prev_resid)
			{
				++residue_count;
				prev_resid = resid;
			}

			block.push_back( line);
		}
	}
	if( block.size() > 0)
	{
		chains.push_back( block);
		counters.push_back( residue_count);
	}
	Close(in);

	std::cout << "found " << chains.size() << " chains" << std::endl << "sizes: " << std::endl;
	for( size_t i  = 0; i < counters.size(); ++i)
	{
		std::cout << counters[i] << std::endl;
	}

	for( unsigned int i = 0; i < counters.size(); ++i)
		if( LENGTH == counters[i])
		{
			++hits;
			hit = i;
		}

	if( hits == 0)
	{
		std::cerr << "ERROR: no chain has expected size!!" << std::endl;
		exit(1);
	}

	if( hits > 1)
	{
		std::cout << "WARNING: more than one chain has expected size, check output, last matching chain found will be considered epitope and assigned chain A!!!" << std::endl;
	}


	for( unsigned int i = 0; i < chains.size(); ++i)
	{
		if( i == hit)
		{
			chain = 'A';
		}
		else
		{
			chain = (char) i_letter++;
		}
		Open( out, OUTPATH + pdb_id + "_" + chain + ".pdb");

		for( std::vector< std::string>::iterator itr = chains[i].begin(); itr != chains[i].end(); ++itr)
		{
			itr->replace( 21,1, &chain,1);
			out << *itr << std::endl;
		}


		Close( out);
	}


}

math::Vector3N
CMS( const std::vector< Molecule> &MOLS)
{
	math::Vector3N
		cms;
	double
		mass = 0;
	for( std::vector< Molecule>::const_iterator mitr = MOLS.begin(); mitr != MOLS.end(); ++mitr)
		for( std::vector< Residue>::const_iterator ritr = mitr->Residues().begin(); ritr != mitr->Residues().end(); ++ritr)
			for( std::vector< Atom>::const_iterator aitr = ritr->Atoms().begin(); aitr != ritr->Atoms().end(); ++aitr)
			{
				cms += aitr->Mass() * aitr->Position();
				mass += aitr->Mass();
			}
	return cms / mass;
}


double
InternalCADistanceDiff( const Molecule &MOL1, const Molecule &MOL2)
{
	Histogram< double> hist( 0.0, 0.5, 100);
	for( std::vector< Residue>::const_iterator r1 = MOL1.Residues().begin(), r2 = MOL2.Residues().begin(); r1 + 1 != MOL1.Residues().end() && r2 + 1 != MOL2.Residues().end(); ++r1, ++r2)
		for( std::vector< Residue>::const_iterator e1 = r1 + 1, e2 = r2 + 1; e1 != MOL1.Residues().end() && e2 != MOL2.Residues().end(); ++e1, ++e2)
		{
			hist( std::abs( math::Distance( r1->Atoms()[1].Position(), e1->Atoms()[1].Position()) - math::Distance( r2->Atoms()[1].Position(), e2->Atoms()[1].Position())));
		}

	// auc like sum
	int
		total = 0.0,
		sum = 0;
	for( std::vector< size_t>::const_iterator itr = hist.GetData().begin(); itr != hist.GetData().end(); ++itr)
	{
		total += (sum += *itr);
	}
	return double( total) / double( hist.GetData().size() * sum);
}


double
InternalCADistanceDiff( const std::vector< Molecule> &MOLS1, const std::vector< Molecule> &MOLS2)
{
	double
		score = 0.0;
	for( std::vector< Molecule>::const_iterator m1 = MOLS1.begin(), m2 = MOLS2.begin(); m1 != MOLS1.end(), m2 != MOLS2.end(); ++m1, ++m2)
	{
		score += InternalCADistanceDiff( *m1, *m2);
	}
	return score;
}



std::vector< std::vector< Atom> >
Atoms( const std::vector< Molecule> &MOLS)
{
	std::vector< std::vector< Atom> >
		mols;

	std::vector< Atom>
		atoms;

	for( std::vector< Molecule>::const_iterator mol = MOLS.begin(); mol != MOLS.end(); ++mol)
	{
		atoms.clear();
		for( std::vector< Residue>::const_iterator aa = mol->Residues().begin(); aa != mol->Residues().end(); ++aa)
			for( std::vector< Atom>::const_iterator atom = aa->Atoms().begin(); atom != aa->Atoms().end(); ++atom)
			{
				atoms.push_back( *atom);
			}
		mols.push_back( atoms);
	}
	return mols;
}




//double
//ForAllAtoms( const std::vector< Molecule> &MOLS, double (* FCT)(Atom))
//{
//	double sum = 0.0;
//	for( std::vector< Molecule>::const_iterator mol = MOLS.begin(); mol != MOLS.end(); ++mol)
//		for( std::vector< Residue>::const_iterator aa = mol->Residues().begin(); aa != mol->Residues().end(); ++aa)
//			for( std::vector< Atom>::const_iterator atom = aa->Atoms().begin(); atom != aa->Atoms().end(); ++atom)
//			{
//				sum += (*FCT)( *atom);
//			}
//	return sum;
//}
// one needs to pass both a pointer to the object as well as a pointer to the member function, then object and function can be linked



//std::vector< Residue>
//BuildFromTopologyFile( std::istream &IN)
//{
//	std::vector<Residue>
//		resiudes;
//
//
//	return residues;
//}


