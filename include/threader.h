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



#ifndef THREADER_H_
#define THREADER_H_

#include "score_matrix.h"



void
ETScan
(
		std::istream & SCORE_MATRIX_LIST,
		std::istream & EPITOPE_LIST,
		std::ostream & OUTFILE
);


void
ScanEpitopes
(
		std::ostream &OUT,
		const std::vector< std::string> &EPITOPES,
		const std::vector< ScoreMatrix> &SCORES
);


void
ScanEpitopes
(
		std::ostream &OUT,
		const std::vector< std::string> &EPITOPES,
		const std::vector< std::pair< double,ScoreMatrix> > &SCORES
);



double
AUC
(
		const std::multimap< double, std::string> &SCORED,
		const std::map< std::string, int> &REFERENCE
);


std::multimap< double, std::string>
ScanEpitopes
(
		const std::vector< std::string> &EPITOPES,
		const std::vector< ScoreMatrix> &SCORE,
		const std::vector< double> &WEIGHTS
);

#endif /* THREADER_H_ */
