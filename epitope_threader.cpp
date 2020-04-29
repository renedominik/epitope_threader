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

#include <vector>

#include <cassert>

#include "include/io.h"
#include "include/centroid.h"
#include "include/threader.h"
#include "include/string_functions.h"
#include "include/manual.h"
#include "include/knowledge_based_potential_calculator.h"


int main( int ARGC, const char **ARGV)
{
	WriteHeader();

	if( ARGC == 1)
	{
		std::cout << "nothing to do" << std::endl;
		return 0;
	}

	std::vector< std::string>
		args( ARGC, "");
	std::ofstream
		out;
	std::ifstream
		list,
		in;

	for( int i = 0; i < ARGC - 1; ++i)
	{
		args[i] = ARGV[i+1];
	}
	if( args[0] == "-cb_centroid")
	{
		if( ARGC != 4)
		{
			return 0;
		}

		Open( in, args[1]);
		Open( out, args[2]);
		CbetaCentroids( in, out);
		Close( in);
		Close( out);
	}
	else if( args[0] == "-ca_centroid")
	{
		if( ARGC != 4)
		{
			return 0;
		}

		Open( in, args[1]);
		Open( out, args[2]);
		CalphaCentroids( in, out);
		Close( in);
		Close( out);
	}
	else if( args[0] == "-score")
	{
		if( ARGC != 6)
		{
			return 0;
		}

		char
			chain = args[2][0];
		std::ifstream
			in_aa;

		Open( in, args[1]);
		Open( in_aa, args[3]);
		Open( out, args[4]);
		CalcScoreMatrix( in, chain, in_aa, out);
		Close( in);
		Close( out);
	}
	else if( args[0] == "-scan")
	{
		if( ARGC != 5 )
		{
			return 0;
		}
		Open( in, args[1]);
		Open( list, args[2]);
		Open( out, args[3]);
		ETScan( in, list, out);
		Close( in);
		Close( out);
		Close( list);
	}
	else if( args[0] == "-pdb2mol")
	{
		if( ARGC != 6)
		{
			return 0;
		}
		std::ifstream
			in_par,
			in_top,
			in_pdb;
		Open( in_par, args[1]);  // param
		Open( in_top, args[2]); // topol
		Open( in_pdb, args[3]); // topol
		Open( out, args[4]); // output
		Pdb2mol( in_par, in_top, in_pdb, args[3], out);
		Close( in_par);
		Close( in_top);
		Close( in_pdb);
		Close( out);
	}
	else if( args[0] == "-default_amino_acids")
	{
		if( ARGC != 5)
		{
			return 0;
		}
		std::ifstream
			in_par,
			in_top;
		Open( in_par, args[1]);  // param
		Open( in_top, args[2]); // topol
		Open( out, args[3]); // output
		DefaultAminoAcidsFromCharmm( in_par, in_top, out);
		Close( in_par);
		Close( in_top);
		Close( out);
	}
	else if( args[0] == "-sort_chains")
	{
		size_t
			length = StringToValue< size_t>( args[2]);
		std::string
			pdb = args[1],
			outpath = args[3];
		SortPDBChainsByEpitopeLength( pdb, length, outpath);
	}
	else if( args[0] == "-kb:pot:dihedral")
	{
		std::ifstream
			pdb_list;
		Open( pdb_list, args[1]);
		KB_Potential_Backbone( pdb_list, /*outpath*/ args[2]);
		Close( pdb_list);
	}
	else if( args[0] == "-kb:pot:aapair")
	{
		std::ifstream
			pdb_list;
		Open( pdb_list, args[1]);
		KB_Potential_AAPair( pdb_list, /*outpath*/ args[2]);
		Close( pdb_list);
	}
	else if( args[0] == "-kb:pot:solvation")
	{
		std::ifstream
			pdb_list;
		Open( pdb_list, args[1]);
		KB_Potential_AASolvation( pdb_list, /*outpath*/ args[2]);
		Close( pdb_list);
	}
	else
	{
		std::cerr << "ERROR: UNDEFINED MODE (" << args[0] << ")" << std::endl << std::endl;
	}
}




