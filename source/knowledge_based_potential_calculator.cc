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


#include <float.h>

#include "../include/knowledge_based_potential_calculator.h"
#include "../include/distribution.h"
#include "../include/math_functions.h"
#include "../include/molecule.h"
#include "../include/std_functions.t.h"
#include "../include/string_functions.h"
#include "../include/io.h"


void
KB_Potential_Backbone( std::istream &PDB_LIST, const std::string &OUTPATH)
{
	size_t
		nr_bins = 36,
		rough_bins = 27;
	std::string
		line;
	std::ifstream
		in_pdb;
	std::ofstream
		ff_phi,
		ff_psi,
		out;
	std::map< std::string, Histogram<double> >
		phi_histograms,
		psi_histograms,
		phi_pair_histograms,
		psi_pair_histograms,
		phi_triplet_histograms,
		psi_triplet_histograms;
	std::map< std::string, Histogram<double> >::iterator
		hist_itr;
	Histogram<double>
		default_hist( -M_PI, 2.0 * M_PI / (float) nr_bins, nr_bins),
		rough_hist( -M_PI, 2.0 * M_PI / (float) rough_bins, rough_bins),
		all_hist( default_hist);
	std::vector< Atom>::const_iterator
		atom;
	double
		phi,
		psi;

	Open( ff_phi, OUTPATH + "/kbdihedral_" + ValueToString<size_t>( nr_bins) + "_phi.txt");
	Open( ff_psi, OUTPATH + "/kbdihedral_" + ValueToString<size_t>( nr_bins) + "_psi.txt");

	// filling histograms
	while( PDB_LIST >> line)
	{
		std::cout <<line << std::endl;
		Open( in_pdb, line);
		std::vector< Molecule>
			mols = ReadPdb( in_pdb);
		Close( in_pdb);

		for( std::vector< Molecule>::const_iterator mol = mols.begin(); mol != mols.end(); ++mol)
			for( std::vector< Residue>::const_iterator res = mol->Residues().begin() + 1; res < mol->Residues().begin() + mol->Residues().size() - 1; ++res) // ignore first and last residue
			{
				atom = res->Atoms().begin();
				if( (res-1)->Atoms().size() < 3 || res->Atoms().size() < 3 || (res+1)->Atoms().size() < 1 || (res-1)->Atoms()[2].Name() != "C" || atom->Name() != "N" || (atom+1)->Name() != "CA" || (atom+2)->Name() != "C" || (res+1)->Atoms()[0].Name() != "N")
				{
					std::cout << "WARNING: incorrect backbone atoms!:  " << res->Name() << " " << int( res - mol->Residues().begin()) << " of " << mol->Residues().size() << ":  ";
					if( (res-1)->Atoms().size() < 3 || res->Atoms().size() < 3 || (res+1)->Atoms().size() < 1)
					{
						std::cout << "too few atoms in residues: prev: " << (res-1)->Atoms().size() << " this: " << res->Atoms().size() << " next: " <<  (res+1)->Atoms().size() << std::endl;
					}
					else
					{
						std::cout << "wrong atom types: "  << (res-1)->Atoms()[2].Name() << "  "<< res->Atoms()[0].Name()<< "  " << res->Atoms()[1].Name()  << "  " << res->Atoms()[2].Name() << " " << (res+1)->Atoms()[0].Name() <<  std::endl;
					}
					continue;
				}
				phi = math::Dihedral( (res-1)->Atoms()[2].Position(), (atom)->Position(), (atom+1)->Position(), (atom+2)->Position()); // check names in fct
				//atom = res->Atoms().begin();
				psi = Dihedral( atom->Position(), (atom+1)->Position(), (atom+2)->Position(), (res+1)->Atoms()[0].Position());

				SafeAdd( phi_histograms, res->Name(), phi, default_hist);
				SafeAdd( psi_histograms, res->Name(), psi, default_hist);

				SafeAdd( phi_pair_histograms, (res-1)->Name() + res->Name(), phi, default_hist);
				SafeAdd( psi_pair_histograms, res->Name() + (res+1)->Name(), psi, default_hist);

				SafeAdd( phi_triplet_histograms, (res-1)->Name() + res->Name() + (res+1)->Name(), phi, rough_hist);
				SafeAdd( psi_triplet_histograms, (res-1)->Name() + res->Name() + (res+1)->Name(), psi, rough_hist);
			}
	}

	std::cout << "write histograms	 " << std::endl;
	// writing histrograms
	for( hist_itr = phi_histograms.begin(); hist_itr != phi_histograms.end(); ++hist_itr)
	{
		Open( out, OUTPATH + "/hist_" + ValueToString<size_t>( nr_bins) + "_phi_" + hist_itr->first + ".txt");
		hist_itr->second.Write(out);
		Close( out);

		Distribution
			normalized( hist_itr->second);
		normalized.ScaleData( (float) nr_bins / normalized.Sum());
		Open( out, OUTPATH + "/norm_" + ValueToString<size_t>( nr_bins) + "_phi_" + hist_itr->first + ".txt");
		normalized.Write( out);
		Close( out);

		Open( out, OUTPATH + "/kbpot_" + ValueToString<size_t>( nr_bins) + "_phi_" + hist_itr->first + ".txt");
		TransformToPotential( normalized).Write( out);
		Close( out);
	}
	for( hist_itr = psi_histograms.begin(); hist_itr != psi_histograms.end(); ++hist_itr)
	{
		Open( out, OUTPATH + "/hist_" + ValueToString<size_t>( nr_bins) + "_psi_" + hist_itr->first + ".txt");
		hist_itr->second.Write(out);
		Close( out);

		Distribution
			normalized( hist_itr->second);
		normalized.ScaleData( (float) nr_bins / normalized.Sum());
		Open( out, OUTPATH + "/norm_" + ValueToString<size_t>( nr_bins) + "_psi_" + hist_itr->first + ".txt");
		normalized.Write( out);
		Close( out);

		Open( out, OUTPATH + "/kbpot_" + ValueToString<size_t>( nr_bins) + "_psi_" + hist_itr->first + ".txt");
		TransformToPotential( normalized).Write( out);
		Close( out);
	}
	for( hist_itr = phi_pair_histograms.begin(); hist_itr != phi_pair_histograms.end(); ++hist_itr)
	{
		Open( out, OUTPATH + "/hist_" + ValueToString<size_t>( nr_bins) + "_phi_" + hist_itr->first + ".txt");
		hist_itr->second.Write(out);
		Close( out);

		Distribution
			normalized( hist_itr->second);
		normalized.ScaleData( (float) nr_bins / normalized.Sum());
		Open( out, OUTPATH + "/norm_" + ValueToString<size_t>( nr_bins) + "_phi_" + hist_itr->first + ".txt");
		normalized.Write( out);
		Close( out);

		Open( out, OUTPATH + "/kbpot_" + ValueToString<size_t>( nr_bins) + "_phi_" + hist_itr->first + ".txt");
		Distribution
			potential = TransformToPotential( normalized);
		potential.Write( out);
		Close( out);

		if( hist_itr->first.find( "UNK") == std::string::npos)
		{
			ff_phi << hist_itr->first << std::endl << potential;
		}
	}
	for( hist_itr = psi_pair_histograms.begin(); hist_itr != psi_pair_histograms.end(); ++hist_itr)
	{
		Open( out, OUTPATH + "/hist_" + ValueToString<size_t>( nr_bins) + "_psi_" + hist_itr->first + ".txt");
		hist_itr->second.Write(out);
		Close( out);

		Distribution
			normalized( hist_itr->second);
		normalized.ScaleData( (float) nr_bins / normalized.Sum());
		Open( out, OUTPATH + "/norm_" + ValueToString<size_t>( nr_bins) + "_psi_" + hist_itr->first + ".txt");
		normalized.Write( out);
		Close( out);

		Open( out, OUTPATH + "/kbpot_" + ValueToString<size_t>( nr_bins) + "_psi_" + hist_itr->first + ".txt");
		Distribution
			potential = TransformToPotential( normalized);
		potential.Write( out);
		Close( out);

		if( hist_itr->first.find( "UNK") == std::string::npos)
		{
			ff_psi << hist_itr->first << std::endl << potential;
		}
	}
	for( hist_itr = phi_triplet_histograms.begin(); hist_itr != phi_triplet_histograms.end(); ++hist_itr)
	{
		Open( out, OUTPATH + "/hist_" + ValueToString<size_t>( nr_bins) + "_phi_" + hist_itr->first + ".txt");
		hist_itr->second.Write(out);
		Close( out);

		Distribution
			normalized( hist_itr->second);
		normalized.ScaleData( (float) nr_bins / normalized.Sum());
		Open( out, OUTPATH + "/norm_" + ValueToString<size_t>( nr_bins) + "_phi_" + hist_itr->first + ".txt");
		normalized.Write( out);
		Close( out);

		Open( out, OUTPATH + "/kbpot_" + ValueToString<size_t>( nr_bins) + "_phi_" + hist_itr->first + ".txt");
		TransformToPotential( normalized).Write( out);
		Close( out);
	}
	for( hist_itr = psi_triplet_histograms.begin(); hist_itr != psi_triplet_histograms.end(); ++hist_itr)
	{
		Open( out, OUTPATH + "/hist_" + ValueToString<size_t>( nr_bins) + "_psi_" + hist_itr->first + ".txt");
		hist_itr->second.Write(out);
		Close( out);

		Distribution
			normalized( hist_itr->second);
		normalized.ScaleData( (float) nr_bins / normalized.Sum());
		Open( out, OUTPATH + "/norm_" + ValueToString<size_t>( nr_bins) + "_psi_" + hist_itr->first + ".txt");
		normalized.Write( out);
		Close( out);

		Open( out, OUTPATH + "/kbpot_" + ValueToString<size_t>( nr_bins) + "_psi_" + hist_itr->first + ".txt");
		TransformToPotential( normalized).Write( out);
		Close( out);
	}
	Close( ff_phi);
	Close( ff_psi);
}






void
KB_Potential_AAPair( std::istream &PDB_LIST, const std::string &OUTPATH)
{
	size_t
		inner_chain_aa_dist = 12, // as in bcl::score
		nr_bins = 200;
	std::string
		bin_str = ValueToString<size_t>( nr_bins),
		res_str,
		line;
	std::ifstream
		in_pdb;
	std::ofstream
		ff_ca,
//		ff_cb,
		out;
	std::map< std::string, Histogram<double> >
		ca_histograms;
//		cb_histograms;
	std::map< std::string, Histogram<double> >::iterator
		hist_itr;
	Histogram<double>
		default_hist( 0.0, 0.25, nr_bins),
		all_hist( default_hist);
	std::vector< Atom>::const_iterator
		atom1,
		atom2;
	Atom
		*a1,
		*a2;
	double
		dist;

	Open( ff_ca, OUTPATH + "/kbaadist_ca_" + bin_str + ".txt");
//	Open( ff_cb, OUTPATH + "/kbaadist_cb_" + bin_str + ".txt");

	// filling histograms
	while( PDB_LIST >> line)
	{
		std::cout <<line << std::endl;
		Open( in_pdb, line);
		std::vector< Molecule>
			mols = ReadPdb( in_pdb);
		Close( in_pdb);

		// inner chain
		for( std::vector< Molecule>::const_iterator mol1 = mols.begin(); mol1 != mols.end(); ++mol1)
			for( std::vector< Residue>::const_iterator res1 = mol1->Residues().begin(); res1 < mol1->Residues().begin() + mol1->Residues().size(); ++res1)
				for( std::vector< Molecule>::const_iterator mol2 = mol1; mol2 < mols.begin() + mols.size(); ++mol2)
					for( std::vector< Residue>::const_iterator res2 = mol2->Residues().begin(); res2 < mol2->Residues().begin() + mol2->Residues().size(); ++res2)
					{
						if( (mol1 == mol2 && res2 - res1 < (int) inner_chain_aa_dist) || res1->Atoms().size() < 2 || res2->Atoms().size() < 2)
						{
							continue;
						}
//						std::cout << "NOW:  " << res1->Name() << " " << int( res1 - mol1->Residues().begin()) << " of " << mol1->Residues().size() << " and "<< res2->Name() << " " << int( res2 - mol2->Residues().begin()) << " of " << mol2->Residues().size() << std::endl;

						atom1 = res1->Atoms().begin() + 1;
						atom2 = res2->Atoms().begin() + 1;

//						std::cout << "atoms : "  << atom1->Name()  << " " << atom2->Name() <<  std::endl;

						if( atom1->Name() != "CA" || atom2->Name() != "CA")
						{
							a1 = (Atom *) res1->FindAtomByName( "CA");
							a2 = (Atom *) res2->FindAtomByName( "CA");
							if( a1 == NULL || a2 == NULL)
							{
								std::cout << "WARNING: incorrect backbone atoms!:  " << std::endl;
								std::cout << res1->Name() << " " << int( res1 - mol1->Residues().begin()) << " of " << mol1->Residues().size() << std::endl;
								std::cout << res2->Name() << " " << int( res2 - mol2->Residues().begin()) << " of " << mol2->Residues().size() << std::endl;
								std::cout << " wrong atom types: "  << atom1->Name()  << " " << atom2->Name() <<  std::endl;
								continue;
							}
							atom1 = std::vector< Atom>::const_iterator( a1);
							atom2 = std::vector< Atom>::const_iterator( a2);
						}

						res_str = SortAndGlue( res1->Name(), res2->Name());

						dist = math::Distance( atom1->Position(), atom2->Position());
						if( dist < 1.5)
						{
							std::cout << "WARNING: very close CALPHAs: " << dist;
							std::cout << res1->Name() << " " << int( res1 - mol1->Residues().begin()) << " of " << mol1->Residues().size() << std::endl;
							std::cout << res2->Name() << " " << int( res2 - mol2->Residues().begin()) << " of " << mol2->Residues().size() << std::endl;
							std::cout  << *atom1  << " \n" << *atom2 <<  std::endl;
						}

						all_hist( dist);

						SafeAdd( ca_histograms, res_str, dist, default_hist);

//						hist_itr = ca_histograms.find( res_str);
//						if( hist_itr == ca_histograms.end())
//						{
//							std::cout << "new " << res_str << std::endl;
//							ca_histograms.insert( std::make_pair( res_str, default_hist));
//							cb_histograms.insert( std::make_pair( res_str, default_hist));
//							hist_itr = ca_histograms.find( res_str);
//						}
//						hist_itr->second( dist);

////						std::cout << "find CB" << std::endl;
//						a1 = res1->FindCBeta();
//						a2 = res2->FindCBeta();
//						if( a1 != NULL && a2 != NULL)
//						{
//							dist = math::Distance( a1->Position(), a2->Position());
//							cb_histograms.find( res_str)->second( dist);
//						}
//						else
//						{
////							std::cout << "WARNING: cbeta not found " << res1->Name() << " " << int( res1 - mol1->Residues().begin()) << " of " << mol1->Residues().size() << " and "<< res2->Name() << " " << int( res2 - mol2->Residues().begin()) << " of " << mol2->Residues().size() << std::endl;
//						}
					}
	}


	Open( out, OUTPATH + "/hist_ca_all.txt");
	all_hist.Write( out);
	Close( out);

	std::cout << "normalize all dist histograms	 " << std::endl;
	Distribution
		all_dist( all_hist);
	all_dist.Normalize();

	Open( out, OUTPATH + "/norm_ca_all.txt");
	all_dist.Write( out);
	Close( out);


	std::cout << "write histograms	 " << std::endl;
	// writing histrograms
	for( hist_itr = ca_histograms.begin(); hist_itr != ca_histograms.end(); ++hist_itr)
	{
		Open( out, OUTPATH + "/hist_ca_" + hist_itr->first + ".txt");
		hist_itr->second.Write(out);
		Close( out);

		Distribution
			ca_dist( hist_itr->second),
			normalized;
		ca_dist.Normalize();

		normalized = ca_dist / all_dist;
		Open( out, OUTPATH + "/ratio_ca_" + hist_itr->first + ".txt");
		normalized.Write( out);
		Close( out);

//		normalized.ScaleData( 100.0);
//		Open( out, OUTPATH + "/norm_ca_" + hist_itr->first + ".txt");
//		normalized.Write( out);
//		Close( out);

		normalized = TransformToPotential( normalized);
		Open( out, OUTPATH + "/kbpot_ca_" + hist_itr->first + ".txt");
		normalized.Write( out);
		Close( out);

		ff_ca << hist_itr->first << std::endl << normalized << std::endl;
	}
//	for( hist_itr = cb_histograms.begin(); hist_itr != cb_histograms.end(); ++hist_itr)
//	{
//		Open( out, OUTPATH + "/hist_cb_" + hist_itr->first + ".txt");
//		hist_itr->second.Write(out);
//		Close( out);
//
//		Distribution
//			cb_dist( hist_itr->second),
//			normalized;
//		cb_dist.Normalize();
//
//		normalized = Devide( cb_dist , all_dist);
//		Open( out, OUTPATH + "/ratio_cb_" + hist_itr->first + ".txt");
//		normalized.Write( out);
//		Close( out);
//
////		normalized.ScaleData( 100.0);
////		Open( out, OUTPATH + "/norm_cb_" + hist_itr->first + ".txt");
////		normalized.Write( out);
////		Close( out);
//
//		normalized = TransformToPotential( normalized);
//		Open( out, OUTPATH + "/kbpot_cb_" + hist_itr->first + ".txt");
//		normalized.Write( out);
//		Close( out);
//
//		ff_cb << hist_itr->first << std::endl << normalized << std::endl;
//	}

	Close( ff_ca);
//	Close( ff_cb);
}






void
KB_Potential_AASolvation( std::istream &PDB_LIST, const std::string &OUTPATH)
{
	size_t
//		inner_chain_aa_dist = 12, // as in bcl::score
		nr_bins = 100;
	std::string
		bin_str = ValueToString<size_t>( nr_bins),
		res_str,
		line;
	std::ifstream
		in_pdb;
	std::ofstream
		ff,
		out;
	std::map< std::string, Histogram<double> >
		aa_histograms;
	std::map< std::string, Histogram<double> >::iterator
		hist_itr;
	Histogram<double>
		default_hist( 0.0, 0.25, nr_bins),
		all_hist( default_hist);

	Open( ff, OUTPATH + "/kbaasolv_" + bin_str + ".txt");


	// filling histograms
	while( PDB_LIST >> line)
	{
		std::cout <<line << std::endl;
		Open( in_pdb, line);
		std::vector< Molecule>
			mols = ReadPdb( in_pdb);
		Close( in_pdb);

		std::vector< Triplet<char,double,char> >
			neighbor_counts = NeighborCount( mols, "all", std::vector< std::string>());

		for( std::vector< Triplet< char,double, char> >::const_iterator itr = neighbor_counts.begin(); itr != neighbor_counts.end(); ++itr)
		{

			if( itr->second <= 2.0)
			{
				std::cout << "WARNING " << __FUNCTION__ << " very low neighbor count " << itr->second << " " << Converter().ThreeLetter( itr->first) << "  " << itr->third << std::endl;
			}

			all_hist( itr->second);

			SafeAdd( aa_histograms, Converter().ThreeLetter( itr->first), itr->second, default_hist);
		}
	}

	Open( out, OUTPATH + "/hist_solvation_all.txt");
	all_hist.Write( out);
	Close( out);

	std::cout << "normalize all dist histograms	 " << std::endl;
	Distribution
		all_dist( all_hist);
	all_dist.Normalize();

	Open( out, OUTPATH + "/norm_solvation_all.txt");
	all_dist.Write( out);
	Close( out);

	std::cout << "write histograms	 " << std::endl;
	// writing histrograms
	for( hist_itr = aa_histograms.begin(); hist_itr != aa_histograms.end(); ++hist_itr)
	{
		Open( out, OUTPATH + "/hist_solvation_" + hist_itr->first + ".txt");
		hist_itr->second.Write(out);
		Close( out);

		Distribution
			aa_dist( hist_itr->second),
			normalized;
		aa_dist.Normalize();
		normalized = aa_dist / all_dist;

		Open( out, OUTPATH + "/ratio_solvation_" + hist_itr->first + ".txt");
		normalized.Write( out);
		Close( out);

//		normalized.ScaleData( 100.0);
//		Open( out, OUTPATH + "/norm_solvation_" + hist_itr->first + ".txt");
//		normalized.Write( out);
//		Close( out);

		normalized = TransformToPotential( normalized);
		Open( out, OUTPATH + "/kbpot_solvation_" + hist_itr->first + ".txt");
		normalized.Write( out);
		Close( out);

		ff << hist_itr->first << std::endl << normalized;
	}
	Close( ff);

}





std::vector< double>
DistributionToPotential( const std::vector< double> &DISTRIBUTION)
{
	std::vector< double>
		v( DISTRIBUTION.size());
	std::vector< double>::iterator vtr = v.begin();
	for( std::vector< double>::const_iterator dtr = DISTRIBUTION.begin(); dtr != DISTRIBUTION.end(); ++dtr, ++vtr)
	{
		*vtr = -kT() * log10( *dtr + DBL_EPSILON);
	}
	return v;
}
