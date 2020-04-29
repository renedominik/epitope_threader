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

#include "../include/threader.h"

#include "../include/std_functions.t.h"
#include "../include/string_functions.h"
#include "../include/io.h"

//#define DistMatrix std::pair< std::string, std::vector< std::vector< std::pair< std::string, double> > > >


void
ETScan
(
		std::istream & SCORE_MATRIX_LIST,
		std::istream & EPITOPE_LIST,
		std::ostream & OUTFILE
)
{
	DebugFct;
	std::vector< ScoreMatrix>
		score_matrices = ReadScoreMatrices( SCORE_MATRIX_LIST);  // HookUp: read from file or connect to SQL
	std::cout << __FUNCTION__ << " " << score_matrices.size() << std::endl;
	std::vector< std::string>
		epitopes = ReadVector< std::string>( EPITOPE_LIST);

	std::cout << score_matrices.size() << " scores, " << epitopes.size() << " epitopes loaded" << std::endl;
	ScanEpitopes( OUTFILE, epitopes, score_matrices);
}


void
ScanEpitopes
(
		std::ostream &OUT,
		const std::vector< std::string> &EPITOPES,
		const std::vector< ScoreMatrix> &SCORES
)
{
	DebugFct;
	double
		value,
		inv = 1.0 / (double) SCORES.size();
	// output:  mhc  tcr  epitope  score
	// distmatrix => mhc,tcr
	// could turn in daemon loop waiting for message
	OUT << "#";
	for( std::vector< ScoreMatrix>::const_iterator scr = SCORES.begin(); scr != SCORES.end(); ++scr)
	{
		OUT  << scr->Who();
		if( scr + 1 != SCORES.end()){ OUT << "::";}
	}
	OUT << std::endl;
	for( std::vector< std::string>::const_iterator etr = EPITOPES.begin(); etr != EPITOPES.end(); ++etr)
	{
		value = 0.0;
		for( std::vector< ScoreMatrix>::const_iterator scr = SCORES.begin(); scr != SCORES.end(); ++scr)
		{
			value += (*scr)( *etr);
		}
		OUT << *etr << ":   " << value * inv << std::endl;
	}
//	for( std::vector< ScoreMatrix>::const_iterator scr = SCORES.begin(); scr != SCORES.end(); ++scr)
//	{
//		for( std::vector< std::string>::const_iterator etr = EPITOPES.begin(); etr != EPITOPES.end(); ++etr)
//		{
//			OUT  << scr->Who() << ":" << *etr << ":   " << (*scr)( *etr) << std::endl;
//		}
//	}
}

void
ScanEpitopes
(
		std::ostream &OUT,
		const std::vector< std::string> &EPITOPES,
		const std::vector< std::pair< double, ScoreMatrix> > &SCORES
)
{
	DebugFct;
	double
		sum,
		value,
		mean,
		minimum = std::numeric_limits< double>::max(),
		inv = 1.0 / (double) SCORES.size();
	std::vector< double>
		values( SCORES.size());
	std::vector< double>::iterator
		valitr;
	// output:  mhc  tcr  epitope  score
	// distmatrix => mhc,tcr
	// could turn in daemon loop waiting for message
	OUT << "#";
	for( std::vector< std::pair< double, ScoreMatrix> >::const_iterator scr = SCORES.begin(); scr != SCORES.end(); ++scr)
	{
		OUT  << scr->first << ":" << scr->second.Who();
		if( scr + 1 != SCORES.end()){ OUT << "::";}
	}
	OUT << std::endl;
	for( std::vector< std::string>::const_iterator etr = EPITOPES.begin(); etr != EPITOPES.end(); ++etr)
	{
		sum = 0.0;
		OUT << *etr << " ";
		valitr = values.begin();
		for( std::vector< std::pair< double, ScoreMatrix> >::const_iterator scr = SCORES.begin(); scr != SCORES.end(); ++scr)
		{
			value = scr->first * scr->second( *etr);
			OUT << value << " ";
			sum += value;
			minimum = std::min( minimum, value);
			*valitr++ = value;
		}
		mean = sum * inv;
		OUT << " min: " << minimum << " mean: " << mean << " dev: " << math::StandardDeviation( values, mean) << std::endl;
	}
}



std::multimap< double, std::string>
ScanEpitopes
(
		const std::vector< std::string> &EPITOPES,
		const std::vector< ScoreMatrix> &SCORES,
		const std::vector< double> &WEIGHTS
)
{
	DebugFct;
	std::multimap< double, std::string>
		sorted;
	double
		sum,
		inv = 1.0 / (double) SCORES.size();
	// output:  mhc  tcr  epitope  score
	// distmatrix => mhc,tcr
	// could turn in daemon loop waiting for message
	for( std::vector< std::string>::const_iterator etr = EPITOPES.begin(); etr != EPITOPES.end(); ++etr)
	{
		sum = 0.0;
		for( std::vector< ScoreMatrix>::const_iterator scr = SCORES.begin(); scr != SCORES.end(); ++scr)
		{
			sum += (*scr)( *etr, WEIGHTS);
		}
		sorted.insert( std::make_pair( sum * inv, *etr));
	}
	return sorted;
}



double
AUC( const std::multimap< double, std::string> &SCORED, const std::map< std::string, int> &REFERENCE)
{
	double
		auc = 0.0;
	size_t
		val,
		cc = 0,
		cy = 0;
	for( std::multimap< double, std::string>::const_iterator scr = SCORED.begin(); scr != SCORED.end(); ++scr)
	{
		val =  REFERENCE.find( scr->second)->second;
		if( val == 1)
		{
			++cy;
		}
		else if( val == 0)
		{
			auc += cy;
			++cc;
		}
	}
	return auc / double( cc * cy);
}
