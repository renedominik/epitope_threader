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


#ifndef DISTRIBUTION_H_
#define DISTRIBUTION_H_

#include <cassert>
#include "histogram.h"

class Distribution
{
protected:
    ///////////////
    //    std::vector< size_t>     //
    ///////////////
    double                 m_Min;
    double                 m_Delta;
    size_t                 m_NrBins;  // = m_Data.size() !!
    double                 m_InverseDelta;
    double                 m_BelowValue;
    double                 m_AboveValue;

private:
	std::vector< double>   m_Data;

public:
	Distribution()
	: m_Min(),
	  m_Delta(),
	  m_NrBins(),
	  m_InverseDelta(),
	  m_BelowValue( 10.0),
	  m_AboveValue( 0.0),
	  m_Data()
	{}
	Distribution( const double &MIN, const double &DELTA, const size_t &NR_BINS)
	: m_Min( MIN),
	  m_Delta( DELTA),
	  m_NrBins( NR_BINS),
	  m_InverseDelta( 1.0 / DELTA),
	  m_BelowValue( 10.0),
	  m_AboveValue( 0.0),
	  m_Data( NR_BINS)
	{}
	Distribution( const double &MIN, const double &DELTA, const size_t &NR_BINS, const std::vector< double> &V)
	: m_Min( MIN),
	  m_Delta( DELTA),
	  m_NrBins( NR_BINS),
	  m_InverseDelta( 1.0 / DELTA),
	  m_BelowValue( 10.0),
	  m_AboveValue( 0.0),
	  m_Data( V)
	{assert( V.size() == NR_BINS);}

	// copy constructor
	Distribution( const Distribution &ORIGINAL)
	: m_Min( ORIGINAL.m_Min),
	  m_Delta( ORIGINAL.m_Delta),
	  m_NrBins( ORIGINAL.m_NrBins),
	  m_InverseDelta( ORIGINAL.m_InverseDelta),
	  m_BelowValue( 10.0),
	  m_AboveValue( 0.0),
	  m_Data( ORIGINAL.m_Data)
	{
//		std::cout << "Distribution copy constructor" << std::endl;
	}

	// construct from histogram
	Distribution( const Histogram<double> &HIST)
	: m_Min( HIST.GetMinimum()),
	  m_Delta( HIST.GetDelta()),
	  m_NrBins( HIST.GetNrBins()),
	  m_InverseDelta( 1.0 / HIST.GetDelta()),
	  m_BelowValue( 10.0),
	  m_AboveValue( 0.0),
	  m_Data( HIST.GetNrBins())
	{
//		std::cout << "Construct Distribution from Histogram" << std::endl;
		std::vector< size_t>::const_iterator itr = HIST.GetData().begin();
		for( std::vector< double>::iterator dtr = m_Data.begin(); dtr != m_Data.end(); ++itr, ++dtr)
		{
			*dtr = double( *itr);
		}
	}

	// virtual destructor
	virtual ~Distribution(){}

	// virtual copy constructor
	virtual Distribution *Clone() const{ return new Distribution( *this);}


	double Value( const double &X) const
	{
		int id = ValueToBin(X);
		if( id < 0 )
		{
			return m_BelowValue;
		}
		else if( id >= (signed) m_Data.size())
		{
			return m_AboveValue;
		}
		return m_Data[ id];
	}

	double InterpolatedValue( const double &X) const
	{
		int id = ValueToBin( X);
		if( id < 0 )
		{
			return m_BelowValue;
		}
		else if( id >= (signed) m_Data.size())
		{
			return m_AboveValue;
		}
		double rest = m_InverseDelta * ( X - m_Min) - (double) id;
//		std::cout << X << " "<< id << " " << rest << std::endl;
		if( rest >= 0.5 && rest <= 1)
		{
			if( id == (signed) m_Data.size()-1)
			{
				return m_Data[id];
			}
			return (1.5 - rest) * m_Data[id] + (rest - 0.5) * m_Data[id+1];
		}
		else if( rest < 0.5 && rest >= 0)
		{
			if( id == 0)
			{
				return m_Data[id];
			}
			return (0.5 + rest) * m_Data[id] + (0.5-rest) * m_Data[id-1];
		}
		else
		{
			std::cout << "WARNING: " << __FUNCTION__ << " should be between 0 and 1: " << rest << " x: " << X << " (x-min)/delta: " << m_InverseDelta * ( X - m_Min) << " id: " << id << " min: " << m_Min << " invdelta: " << m_InverseDelta << " delta: " << m_Delta <<  std::endl;
			return 0.0;
		}
		return m_Data[ id];
	}


	std::vector< double> & GetData()
	{
		return m_Data;
	}


	const std::vector< double> & GetData() const
	{
		return m_Data;
	}


    int ValueToBin( const double &VALUE) const
     {
         if( VALUE < m_Min)
         {
         	return -1;
         }
         else if( VALUE > m_Min + m_NrBins * m_Delta)
         {
         	return std::numeric_limits< int>::max();
         }
         return size_t( ( VALUE - m_Min) * m_InverseDelta);
     }


    void
    ScaleData( const double &FACTOR)
    {
    	for( std::vector< double>::iterator mtr = m_Data.begin(); mtr != m_Data.end(); ++mtr)
    	{
    		*mtr = *mtr * FACTOR;
    	}
    }

    void
    Normalize()
    {
    	double
    		factor = 1.0 / Sum();
    	for( std::vector< double>::iterator mtr = m_Data.begin(); mtr != m_Data.end(); ++mtr)
    	{
    		*mtr = *mtr * factor;
    	}
    }

    double
    Sum() const
    {
    	double sum = 0.0;
    	for( std::vector< double>::const_iterator itr = m_Data.begin(); itr != m_Data.end(); ++itr)
    	{
    		sum += *itr;
    	}
    	return sum;
    }

    virtual std::ostream &
    Write( std::ostream &STREAM) const
    {
    	double x = m_Min + 0.5 * m_Delta;
    	STREAM << m_Min << std::endl;
    	STREAM << m_Delta << std::endl;
    	STREAM << m_NrBins << std::endl;
    	for( std::vector< double>::const_iterator itr = m_Data.begin(); itr != m_Data.end(); ++itr, x += m_Delta)
        {
            STREAM << x << "  " << *itr << std::endl;
        }
        return STREAM;
    }

	std::istream &
	Read( std::istream &STREAM)
	{
		double x;
		STREAM >> m_Min >> m_Delta >> m_NrBins;
		m_InverseDelta = 1.0 / m_Delta;
		m_BelowValue = 10.0;
	    m_AboveValue = 0.0;
		m_Data = std::vector< double>( m_NrBins);
		for( std::vector< double>::iterator itr = m_Data.begin(); itr != m_Data.end(); ++itr)
		{
			STREAM >> x >> *itr;
		}
		return STREAM;
	}



    friend Distribution TransformToPotential( const Distribution &D);
};


Distribution
Devide( const Distribution &D1, const Distribution &D2);

Distribution
operator / ( const Distribution &D1, const Distribution &D2);



Distribution
TransformToPotential( const Distribution &D);

inline
std::istream &
operator >> ( std::istream &STREAM, Distribution & D)
{
	return D.Read( STREAM);
}

inline
std::ostream &
operator << ( std::ostream &STREAM, const Distribution & D)
{
	return D.Write( STREAM);
}

#endif /* DISTRIBUTION_H_ */
