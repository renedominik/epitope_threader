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


#include "../include/forcefield.h"
#include "../include/string_functions.h"
#include "../include/std_functions.t.h"
#include "../include/io.h"
#include "../include/global.h"

const double &ForceFieldParameters::Charge( const std::string &TYPE) const
{
	std::map< std::string, double>::const_iterator itr = ForceFieldParameters::m_TypeChargeMap.find( TYPE);
	if( itr != ForceFieldParameters::m_TypeChargeMap.end())
	{
		return itr->second;
	}
	std::cout << "WARNING: " << TYPE << " not found in m_TypeChargeMap in class ForceFieldParameters" << std::endl;
	exit(1);
}


const double &ForceFieldParameters::Epsilon( const std::string &TYPE) const
{
	std::map< std::string, double>::const_iterator itr = ForceFieldParameters::m_TypeEpsilonMap.find( TYPE);
	if( itr != ForceFieldParameters::m_TypeEpsilonMap.end())
	{
		return itr->second;
	}
	std::cout << "WARNING: " << TYPE << " not found in m_TypeEpsilonMap in class ForceFieldParameters" << std::endl;
	exit(1);
}


const double &ForceFieldParameters::Radius( const std::string &TYPE) const
{
	std::map< std::string, double>::const_iterator itr = ForceFieldParameters::m_TypeRadiusMap.find( TYPE);
	if( itr != ForceFieldParameters::m_TypeRadiusMap.end())
	{
		return itr->second;
	}
	std::cout << "WARNING: " << TYPE << " not found in m_TypeRadiusMap in class ForceFieldParameters" << std::endl;
	exit(1);
}

const double &ForceFieldParameters::Mass( const std::string &TYPE) const
{
	return SafeGet( m_TypeMassMap, TYPE);
}



void ForceFieldParameters::ReadPartialChargeMapFromCharmmTopologyFile( std::istream &STREAM)
{
	DebugFct;
    std::string line, residue, atom, block, key;
    double charge;
    while( STREAM)
    {
        line.clear();
        std::getline( STREAM, line);
        if( line.size() > 5)
        {
            std::vector< std::string> vec = Split( line);
            if( vec[0] == "RESI" || vec[0] == "PRES")
            {
                block = vec[0];
                residue = vec[1];
                if( residue == "HSD") // repair
                {
                    residue = "HIS";
                }
            }
            else if( vec[0] == "ATOM" && (   block == "RESI"
                    || ( block == "PRES" && residue == "CTER")
                    || ( block == "PRES" && residue == "NTER")))
            {
                atom = vec[1];
                charge = StringToValue< double>( vec[3]);
                key = residue + ":" + atom;
                UniqueInsert( m_TypeChargeMap, key, charge);
            }
        }
    }
}


void ForceFieldParameters::ReadMassFromCharmmTopologyFile( std::istream &STREAM)
{
	DebugFct;
    std::string line, atomtype;
    double mass;
    while( STREAM)
    {
    	line.clear();
        std::getline( STREAM, line);
        if( line.size() > 0)
        {
            std::vector< std::string> vec( Split( line));
            if( vec[0] == "MASS")
            {
            	atomtype = vec[2];
                mass = StringToValue< double>( vec[3]);
                UniqueInsert( m_TypeMassMap, atomtype, mass);
            }
            if( vec[0] == "DECL")
            { break;}
        }
    }
}



void ForceFieldParameters::ReadEpsilonAndRadiusFromCharmmParameterFile( std::istream &STREAM)
{
	DebugFct;
    std::string line, residue, atomtype;
    bool not_found( true);
    while( STREAM && not_found)
    {
        std::getline( STREAM, line);
        if( line.find( "epsilon") != std::string::npos)
        {
            not_found = false;
        }
    }
    while( STREAM)
    {
        IgnoreComments( STREAM, '!');
        std::getline( STREAM, line);
        if( line.size() > 0)
        {
            std::vector< std::string> vec( Split( line));
            if( vec.size() > 2 && vec[0].length() <= 4 && IsCapitolLetter( vec[0][0]) && IsNumerical( vec[2]))
            {
            	atomtype = vec[0];
                UniqueInsert( m_TypeEpsilonMap, atomtype, StringToValue<double>( vec[2]));
                UniqueInsert( m_TypeRadiusMap, atomtype, StringToValue<double>( vec[3]));
            }
        }
    }
}



void ForceFieldParameters::ReadAtomTypeFromCharmmTopologyFile( std::istream &STREAM)
{
	DebugFct;
    std::string line, residue, block, ct3 = "CT3", cd1 = "ILE:CD1";
    std::vector< std::string> vec;
    while( STREAM)
    {
        line.clear();
        std::getline( STREAM, line);
        if( line.size() > 5)
        {
            vec = Split( line);
            if( vec[0] == "RESI" || vec[0] == "PRES")
            {
                block = vec[0];
                residue = vec[1];
                if( residue == "HSD") // repair
                {
                    residue = "HIS";
                }
            }
            else if( vec[0] == "ATOM" && (   block == "RESI"
                    || ( block == "PRES" && residue == "CTER")
                    || ( block == "PRES" && residue == "NTER")))
            {

            	UniqueInsert( m_NameTypeMap, residue + ":" + vec[1], vec[2]);
            }
        }
    }
   	UniqueInsert( m_NameTypeMap, cd1, ct3);
}



Molecule & ForceFieldParameters::WritePartialChargeIntoMolecule( Molecule &MOL) const
{
	DebugFct;
    std::string residue, atomtype;
    for( std::vector< Residue>::iterator rtr = MOL.Residues().begin(); rtr != MOL.Residues().end(); ++rtr)
    {
		residue = rtr->Name();
    	for( std::vector< Atom>::iterator atr = rtr->Atoms().begin(); atr != rtr->Atoms().end(); ++atr)
		{
    		atomtype = atr->Name();
			if( residue == "ILE" && atomtype == "CD1") // repair
			{
				atomtype = "CD";
			}
			else if( residue == "LYS" && atomtype == "OXT") // repair
			{
				atomtype = "O";
			}
			atr->Charge() = CorrectTermini( m_TypeChargeMap, residue, atomtype);
//			atr->Charge() = SafeGet( m_TypeChargeMap, CorrectTermini( m_NameTypeMap, residue, atomtype));
		}
    }
    return MOL;
}


Molecule &ForceFieldParameters::WriteMassIntoMolecule( Molecule &MOL) const
{
	DebugFct;
    std::string atomtype, residue;
    for( std::vector< Residue>::iterator rtr = MOL.Residues().begin(); rtr != MOL.Residues().end(); ++rtr)
    {
    	residue = rtr->Name();
    	for( std::vector< Atom>::iterator atr = rtr->Atoms().begin(); atr != rtr->Atoms().end(); ++atr)
    	{
    		atomtype = atr->Type();
//                if( residue == "ILE" && atom == "CD1") // repair
//                {
//                    atom = "CD";
//                }
//                else if( residue == "LYS" && atom == "OXT") // repair
//                {
//                    atom = "O";
//                }
    		atr->Mass() = SafeGet( m_TypeMassMap, atomtype);
//    		atr->Mass() = CorrectTermini( m_TypeMassMap, residue, atomtype); // WHICH TO USE???
    	}
    }
//    std::cout  << __FUNCTION__ <<  " done" <<  std::endl;
    return MOL;
}



Molecule &ForceFieldParameters::WriteEpsilonAndRadiusIntoMolecule( Molecule &MOL) const
{
	DebugFct;
    std::string
    	residue,
    	atomtype;

    for( std::vector< Residue>::iterator rtr = MOL.Residues().begin(); rtr != MOL.Residues().end(); ++rtr)
    {
    	residue = rtr->Name();
    	for( std::vector< Atom>::iterator atr = rtr->Atoms().begin(); atr != rtr->Atoms().end(); ++atr)
    	{
    		atomtype = atr->Type();

			if( atomtype == "OXT")
			{
				atomtype = "OC";
			}
			atr->VdwEpsilon() = SafeGet( m_TypeEpsilonMap, atomtype);
			atr->VdwRadius() = SafeGet( m_TypeRadiusMap, atomtype);
    	}
    }
    return MOL;
}


Molecule &ForceFieldParameters::WriteAtomTypeIntoMolecule( Molecule &MOL) const
{
	DebugFct;
	std::string resname;
    for( std::vector< Residue>::iterator rtr = MOL.Residues().begin(); rtr != MOL.Residues().end(); ++rtr)
    {
    	resname = rtr->Name();
    	for( std::vector< Atom>::iterator atr = rtr->Atoms().begin(); atr != rtr->Atoms().end(); ++atr)
    	{
    		atr->Type() = CorrectTermini( m_NameTypeMap, resname, atr->Name());
    	}
    }
    return MOL;
}



Molecule &ForceFieldParameters::WriteAllIntoMolecule( Molecule &MOL) const
{
	DebugFct;
	WriteAtomTypeIntoMolecule( MOL);
	WriteMassIntoMolecule( MOL);
	WritePartialChargeIntoMolecule( MOL);
	WriteEpsilonAndRadiusIntoMolecule(MOL);
	return MOL;
}


std::ostream &ForceFieldParameters::Write( std::ostream& OUT) const
{
	DebugFct;
	std::cout << "m_TypeChargeMap" << std::endl;
	WriteMap<std::string, double>( m_TypeChargeMap, OUT);
	std::cout << "m_TypeEpsilonMap" << std::endl;
	WriteMap( m_TypeEpsilonMap, OUT);
	std::cout << "m_TypeRadiusMap" << std::endl;
	WriteMap( m_TypeRadiusMap, OUT);
	std::cout << "m_TypeMassMap" << std::endl;
	WriteMap( m_TypeMassMap, OUT);
	std::cout << "m_NameTypeMap" << std::endl;
	WriteMap( m_NameTypeMap, OUT);
	return OUT;
}




ForceFieldParameters ReadForceFieldParameters( std::istream &IN_PAR, std::istream &IN_TOP)
{
	DebugFct;
	ForceFieldParameters parameters;
	parameters.ReadAtomTypeFromCharmmTopologyFile( IN_TOP);
	// mv stream back to beginning of input file
	IN_TOP.clear();
	IN_TOP.seekg( 0, std::ios::beg);
	parameters.ReadPartialChargeMapFromCharmmTopologyFile( IN_TOP);
	IN_TOP.clear();
	IN_TOP.seekg( 0, std::ios::beg);
	parameters.ReadMassFromCharmmTopologyFile( IN_TOP);

	parameters.ReadEpsilonAndRadiusFromCharmmParameterFile( IN_PAR);
#ifdef DEBUG
	parameters.Write();
#endif
	return parameters;
}


void WriteForceFieldParametersIntoMolecules( const ForceFieldParameters &PAR, std::vector< Molecule> &MOLS)
{
	DebugFct;
	for( std::vector< Molecule>::iterator itr = MOLS.begin(); itr != MOLS.end(); ++itr)
	{
		PAR.WriteAllIntoMolecule( *itr);
	}
}


double ForceField::Score( const std::vector< Molecule> &MOLS, const std::vector< std::vector< size_t> > &SORTED_IDS) const
{
//	std::cout << __FUNCTION__ << " mols ids" << std::endl;
	double
		energy = 0.0;
	size_t
		id,
		first;
	std::vector< Molecule>
		first_mols,
		second_mols;
	std::vector< size_t>::const_iterator
		itr;

	// all mols score , e.g. exposure, solvation
//	if( !sm_AASolvationScore.empty())
//	{
//		energy += ScoreAASolvation( MOLS);
//	}

	for( size_t chain = 0; chain < SORTED_IDS.size(); ++chain)
	{
//		std::cout << "iteration: " << chain << std::endl;
		itr = SORTED_IDS[chain].begin(); // first mutated residue
		if( itr != SORTED_IDS[chain].end())
		{
			first = *itr;
		}
		else
		{
			continue;
		}
		id = first;
		first_mols = std::vector< Molecule>( MOLS.begin(), MOLS.begin() + chain); // potential shit, too much copying?
		second_mols = std::vector< Molecule>( MOLS.begin() + chain + 1, MOLS.begin() + MOLS.size());
//		std::cout << "inter chain " << first_mols.size()
//			<< " " << second_mols.size() << std::endl;

		for( std::vector< Residue>::const_iterator res = MOLS[chain].Residues().begin() + first; res < MOLS[chain].Residues().begin() + MOLS[chain].Residues().size(); ++res, ++id)
		{
			// inner chain
			energy += Score( MOLS[chain], id);

			// inter chain
			energy += Score( *res, first_mols);
			energy += Score( *res, second_mols);
		}
	}
	m_Current = energy;
	return energy;
}




double ForceField::Score( const Molecule &MOL, const size_t &RESID) const
{
//	std::cout << __FUNCTION__ << " mol id " << RESID << std::endl;
	double
		energy = ScoreNextNeighbor( MOL, RESID);

//	energy += ScoreBackbone( MOL, RESID);

	const Residue
		&r1 = MOL.Residues()[RESID];

	for( std::vector< Residue>::const_iterator res = MOL.Residues().begin(); res < MOL.Residues().begin() + RESID - 1; ++res)
	{
		energy += Score( r1, *res);
//		std::cout << int(res - MOL.Residues().begin()) << " ";
	}
	for( std::vector< Residue>::const_iterator res = MOL.Residues().begin() + RESID + 2; res < MOL.Residues().begin() + MOL.Residues().size(); ++res)
	{
		energy += Score( r1, *res);
//		std::cout << int(res - MOL.Residues().begin()) << " ";
	}
//	std::cout << std::endl;
//	std::cout << __FUNCTION__ << ": " << i << " " << energy << std::endl;
	return energy;
}



double ForceField::ScoreNextNeighbor( const Molecule &MOL, const size_t &RESID) const
{
//	std::cout << __FUNCTION__ << " mol id " << RESID << " of " << MOL.Residues().size() <<  std::endl;
	if( RESID + 1 >= MOL.Residues().size() || RESID - 1 < 0)
	{
//		std::cout << __FUNCTION__ << " mol id " << RESID << " return 0 " <<  std::endl;
		return 0.0;
	}

	double
		energy = 0.0;

	const Residue
		*r2,
		*r1 = &MOL.Residues()[RESID];
	for( size_t i = RESID - 1; i <= RESID + 1; i+=2)
	{
		r2 = &MOL.Residues()[i];
		for( std::vector< Atom>::const_iterator a1 = r1->Atoms().begin() + 3; a1 != r1->Atoms().end(); ++a1)
				for( std::vector< Atom>::const_iterator a2 = r2->Atoms().begin() + 3; a2 != r2->Atoms().end(); ++a2)
//					if( a1->Name() != "N" && a1->Name() != "CA" && a1->Name() != "C" && a2->Name() != "N" && a2->Name() != "CA" && a2->Name() != "C") //  mv away, rather check consitency of input data once at start!!
					{
						energy += Score( *a1, *a2);
					}
	}
//	std::cout << __FUNCTION__ << " mol id " << RESID << " return " <<  energy << std::endl;
	return energy;
}





double ForceField::Score( const Residue & RES, const std::vector< Molecule> &MOLS) const
{
//	std::cout << __FUNCTION__ << " mols " << MOLS.size() << std::endl;
	double
		energy = 0.0;

	for( std::vector< Molecule>::const_iterator mol = MOLS.begin(); mol != MOLS.end(); ++mol)
		for( std::vector< Residue>::const_iterator res = mol->Residues().begin(); res != mol->Residues().end(); ++res)
		{
			// inner chain
			energy += Score( RES, *res);
		}
	return energy;
}





double ForceField::Score( const std::vector< Molecule> &MOLS) const
{
//	std::cout << __FUNCTION__ << " mols" << std::endl;
	double energy = 0.0;

	// intramolecular
	for( std::vector< Molecule>::const_iterator m = MOLS.begin(); m != MOLS.end(); ++m)
	{
		energy += Score( *m);
	}

	// inter molecular
	for( std::vector< Molecule>::const_iterator m1 = MOLS.begin(); m1 != MOLS.end(); ++m1)
		for( std::vector< Molecule>::const_iterator m2 = m1+1; m2 != MOLS.end(); ++m2)
			{
				energy += Score( *m1, *m2);
			}

	m_Current = energy;
	return energy;
}


double ForceField::Score( const Molecule &M) const
{
//	std::cout << __FUNCTION__ << " mol " << std::endl;

	double energy = 0.0;

	// inner-residue ????

	// next neighbor: ignore backbone atoms
	for( std::vector< Residue>::const_iterator r1 = M.Residues().begin(), r2 = r1 + 1; r1 + 1 != M.Residues().end(); ++r1, ++r2)
		for( std::vector< Atom>::const_iterator a1 = r1->Atoms().begin() + 3; a1 < r1->Atoms().end(); ++a1)
			for( std::vector< Atom>::const_iterator a2 = r2->Atoms().begin() + 3; a2 < r2->Atoms().end(); ++a2)
//				if( a1->Name() != "N" && a1->Name() != "CA" && a1->Name() != "C" && a2->Name() != "N" && a2->Name() != "CA" && a2->Name() != "C")
				{
					energy += Score( *a1, *a2);
				}

	// non-next neighbors (min resdist 2)
	for( std::vector< Residue>::const_iterator r1 = M.Residues().begin(); r1 + 2 != M.Residues().end(); ++r1)
		for( std::vector< Residue>::const_iterator r2 = r1 + 2; r2 < M.Residues().end(); ++r2)
		{
			energy += Score( *r1, *r2);
		}
	return energy;
}


double ForceField::Score( const Molecule &M1, const Molecule &M2) const
{
//	std::cout << __FUNCTION__ << " mol mol" << std::endl;
	double energy = 0.0;
	for( std::vector< Residue>::const_iterator r1 = M1.Residues().begin(); r1 != M1.Residues().end(); ++r1)
			for( std::vector< Residue>::const_iterator r2 = M2.Residues().begin(); r2 != M2.Residues().end(); ++r2)
			{
				energy += Score( *r1, *r2);
			}
	return energy;
}


double ForceField::Score( const Residue &R1, const Residue &R2) const
{
//	std::cout << __FUNCTION__ << " res res" << std::endl;
	double energy = 0.0;
	for( std::vector< Atom>::const_iterator a1 = R1.Atoms().begin(); a1 != R1.Atoms().end(); ++a1)
			for( std::vector< Atom>::const_iterator a2 = R2.Atoms().begin(); a2 != R2.Atoms().end(); ++a2)
			{
				energy += Score( *a1, *a2);
			}
	return energy;
}


double ForceField::Score( const Atom &A1, const Atom &A2) const
{
	double dist = math::Distance( A1.Position(), A2.Position());
	return Coulomb( A1, A2, dist) + VanDerWaals( A1, A2, dist);
}


double ForceField::Score( const std::vector< std::vector< Atom> > &MOLS) const
{
	std::cout << __FUNCTION__ << std::endl;
	double score = 0.0;

	size_t a2;
	for( size_t m1 = 0; m1 < MOLS.size(); ++m1)
		for( size_t a1 = 0; a1 < MOLS[m1].size(); ++a1)
			for( size_t m2 = m1; m2 < MOLS.size(); ++m2)
			{
				if( m1 == m2)
				{

					if( a1 + 15 < MOLS[m1].size())
					{
						a2 = a1 + 15;
					}
					else
					{
						a2 = std::numeric_limits< size_t>::max();
					}
				}
				else
				{
					a2 = 0;
				}

				for( ; a2 < MOLS[m2].size(); ++a2)
				{

					try
					{
						score += Coulomb( MOLS[m1][a1], MOLS[m2][a2]) + VanDerWaals( MOLS[m1][a1], MOLS[m2][a2]);
						if( std::isnan(score))
						{
							std::cout << "isnan: coul: " << Coulomb( MOLS[m1][a1], MOLS[m2][a2]) << " vdw: " <<  VanDerWaals( MOLS[m1][a1], MOLS[m2][a2])	<< " dist: " << math::Distance( MOLS[m1][a1].Position(), MOLS[m2][a2].Position()) << std::endl;
							MOLS[m1][a1].Write();
							MOLS[m2][a2].Write();
							std::cout << "sizes: " << MOLS[m1].size() << " " << MOLS[m2].size() << std::endl;
//							<< " names: " << MOLS[m1][a1].Name() << " " << MOLS[m1][a1].VdwRadius() << " " << MOLS[m1][a1].Charge() << " " << MOLS[m1][a1].Mass()
//								<< " " << MOLS[m2][a2].Name()  << " " << MOLS[m2][a2].VdwRadius() << " " << MOLS[m2][a2].Charge() << " " << MOLS[m2][a2].Mass()<< std::endl;

							std::cerr << m1 << " " << a1 << " " << m2 <<  "  " << a2 << std::endl;
							exit(1);
						}
					}
					catch(...)
					{
						std::cerr << m1 << " " << a1 << " " << m2 <<  "  " << a2 << std::endl;
						exit(1);
					}
				}
			}
	m_Current = score;
	std::cout << __FUNCTION__ << ": " << score << std::endl;
	return score;
}


double Coulomb( const Atom &A1, const Atom &A2, const double &DIST)
{
	return 332.063711 * A1.Charge() * A2.Charge() / ( DIST * DIST);
}

double Coulomb( const Atom &A1, const Atom &A2)
{
	double factor = 332.063711; // Coulomb's constant , see NAMD user guide 2.9 page 62 // in the source code it is 332.0636 ...
	double dist = math::Distance( A1.Position(), A2.Position());
//	std::cout << __FUNCTION__ << " dist: " << dist << std::endl;
//	A1.Write();
//	A2.Write();
	return factor * A1.Charge() * A2.Charge() / ( dist * dist);
}

double VanDerWaals( const Atom &A1, const Atom &A2, const double &DIST)
{
    double
    	rmindist = (A1.VdwRadius() + A2.VdwRadius()) / DIST;
    return sqrt( A1.VdwEpsilon() * A2.VdwEpsilon()) * ( std::pow( rmindist, 12) - 2.0 * std::pow( rmindist, 6));
}

double VanDerWaals( const Atom &A1, const Atom &A2)
{
/*
    !V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
    !
    !epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
    !Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
 */

    double
    	rmindist = (A1.VdwRadius() + A2.VdwRadius()) / math::Distance( A1.Position(), A2.Position());

 //   std::cout << __FUNCTION__ << " rmindist: " << rmindist << std::endl;

    return sqrt( A1.VdwEpsilon() * A2.VdwEpsilon()) * ( std::pow( rmindist, 12) - 2.0 * std::pow( rmindist, 6));
}


double SumFunction( const Atom &A1, const Atom &A2)
{
	static std::vector < double (*)( const Atom &A1, const Atom &A2)> sg_fct_ptr;
	double sum = 0.0;
	for( std::vector< double (*)( const Atom &A1, const Atom &A2)>::const_iterator itr = sg_fct_ptr.begin(); itr != sg_fct_ptr.end(); ++itr)
	{
		sum += (**itr)( A1, A2);
	}
	return sum;
}

void Pdb2mol
(
		std::istream &IN_PAR,
		std::istream &IN_TOP,
		std::istream &IN_PDB,
		const std::string &MOL_NAME,
		std::ostream &OUT
)
{
	DebugFct;
	ForceFieldParameters
		parameters = ReadForceFieldParameters( IN_PAR, IN_TOP);

	std::vector< Molecule>
		mols = ReadPdb( IN_PDB, MOL_NAME);

	WriteForceFieldParametersIntoMolecules( parameters, mols);

	Write( mols, OUT);
}




void DefaultAminoAcidsFromCharmm
(
		std::istream &IN_PAR,
		std::istream &IN_TOP,
		std::ostream &OUT
)
{
	DebugFct;
	Molecule mol( "default_amino_acids");
	Residue residue;
	double charge;
    std::string line, resname, block, name, type;
    std::vector< std::string> vec;
    math::Vector3N zero;

    ForceFieldParameters parameters;

    parameters.ReadEpsilonAndRadiusFromCharmmParameterFile( IN_PAR);
    parameters.ReadMassFromCharmmTopologyFile( IN_TOP);

    IN_TOP.clear();
	IN_TOP.seekg( 0, std::ios::beg);

    while( IN_TOP)
    {
        line.clear();
        std::getline( IN_TOP, line);
        if( line.size() > 5)
        {
            vec = Split( line);
            if( vec[0] == "RESI")
            {
                block = vec[0];
                resname = vec[1];
                if( resname == "HSD") // repair
                {
                	resname = "HIS";
                }
            	if( residue.Atoms().size() > 0 && IsNaturalAminoAcid( residue.Name()))
            	{
            		mol.Residues().push_back( residue);
            	}
        		residue = Residue( resname);
            }
            else if( vec[0] == "ATOM" &&  block == "RESI")
            {
            	name = vec[1];
            	type = vec[2];
            	charge = StringToValue< double>( vec[3]);
            	residue.Atoms().push_back( Atom( zero, name, type, charge, parameters.Epsilon( type), parameters.Radius( type), parameters.Mass( type)));
            }
        }
    }
	if( residue.Atoms().size() > 0 && IsNaturalAminoAcid( residue.Name()))
	{
		mol.Residues().push_back( residue);
	}
	OUT << "1" << std::endl;
	mol.Write(OUT);
}


