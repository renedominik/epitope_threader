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


#include "../include/distribution.h"

#include <float.h>
#include <math.h>


Distribution TransformToPotential( const Distribution &D)
{
	Distribution d( D);
	for( std::vector< double>::iterator mtr = d.m_Data.begin(); mtr != d.m_Data.end(); ++mtr)
	{
		*mtr = -0.593 * log10( *mtr + DBL_EPSILON);
	}
	return d;
}


Distribution
Devide( const Distribution &D1, const Distribution &D2)
{
	Distribution
		dist( D1);

	std::vector< double>::const_iterator d2 = D2.GetData().begin();
	for( std::vector< double>::iterator d1 = dist.GetData().begin(); d1 != dist.GetData().end(); ++d1, ++d2)
	{

		*d1 /= *d2 + DBL_EPSILON;  // DBL_EPSILON to avoid nan when count equal zero
	}

	return dist;
}


Distribution
operator / ( const Distribution &D1, const Distribution &D2)
{
	return Devide( D1, D2);
}

