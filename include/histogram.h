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


#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include <vector>
#include <iterator>
#include <iostream>
#include <limits>

//namespace math
//{
    template< typename t_INPUT>
    class Histogram
    {
    protected:
        ///////////////
        //    std::vector< size_t>     //
        ///////////////
        t_INPUT                 m_Min;
        t_INPUT                 m_Delta;
        size_t                  m_NrBins;  // = m_Data.size() !!
        t_INPUT                 m_InverseDelta;

    private:
        std::vector< size_t>    m_Data;
        mutable size_t          m_BelowLowestBin;
        mutable size_t          m_AboveHighestBin;
        std::vector< t_INPUT>   m_Above;
        std::vector< t_INPUT>   m_Below;


    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        // construct from std::vector< size_t>
        Histogram( const t_INPUT &MIN, const t_INPUT &DELTA, const size_t &NR_BINS)
        :
        m_Min( MIN),
        m_Delta( DELTA),
        m_NrBins( NR_BINS),
	    m_InverseDelta( t_INPUT(1) / DELTA),
        m_Data( NR_BINS),
        m_BelowLowestBin( 0),
        m_AboveHighestBin( 0),
	    m_Above(),
	    m_Below()
        {}

        // copy constructor
        Histogram( const Histogram &ORIGINAL)
        :
        m_Min( ORIGINAL.m_Min),
        m_Delta( ORIGINAL.m_Delta),
        m_NrBins( ORIGINAL.m_NrBins),
	    m_InverseDelta( ORIGINAL.m_InverseDelta),
        m_Data( ORIGINAL.m_Data),
        m_BelowLowestBin( ORIGINAL.m_BelowLowestBin),
        m_AboveHighestBin( ORIGINAL.m_AboveHighestBin),
	    m_Above( ORIGINAL.m_Above),
	    m_Below( ORIGINAL.m_Below)
        {
//        	std::cout << "Histogram copy constructor" << std::endl;
        }

        // virtual destructor
        virtual ~Histogram(){}

        // virtual copy constructor
        virtual Histogram *Clone() const{ return new Histogram( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////
        virtual size_t operator() ( const t_INPUT &INPUT)
        {
            size_t
				id = ValueToBin( INPUT);
            if( id != std::numeric_limits< size_t>::max() && id >= 0)
            {
            	return ++m_Data[ id];
            }
            return std::numeric_limits<size_t>::max();
        }

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        virtual const t_INPUT &GetMinimum() const
        {
            return m_Min;
        }

        virtual const t_INPUT &GetDelta() const
        {
            return m_Delta;
        }

        virtual const size_t &GetNrBins() const
        {
            return m_NrBins;
        }

        const std::vector< size_t> &GetData() const
        {
            return m_Data;
        }


        const t_INPUT &GetInverseDelta() const
        {
            return m_InverseDelta;
        }

        const size_t &GetAboveHighestBin() const
        {
            return m_AboveHighestBin;
        }

        const size_t &BelowLowestBin() const
        {
            return m_BelowLowestBin;
        }


        size_t ValueToBin( const t_INPUT &VALUE) const
        {
            if( VALUE < m_Min)
            {
//            	std::cout <<__FUNCTION__ << " " << VALUE << " is below " << m_Min << std::endl;
            	++m_BelowLowestBin;
            	return std::numeric_limits< size_t>::max();
            }
            else if( VALUE > m_Min + m_NrBins * m_Delta)
            {
//            	std::cout <<__FUNCTION__ << " " << VALUE << " is above " << m_Min + m_NrBins * m_Delta << std::endl;
            	++m_AboveHighestBin;
            	return std::numeric_limits< size_t>::max();
            }
            return size_t( ( VALUE - m_Min) * m_InverseDelta);
        }

        size_t Sum() const
        {
            size_t sum( 0);
            for( std::vector< size_t>::const_iterator itr = m_Data.begin(); itr != m_Data.end(); ++itr)
            {
                sum += *itr;
            }
            return sum;
        }

        void Clear()
        {
        	std::fill( m_Data.begin(), m_Data.end(), 0);
        }

        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM)
        {
            return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
//            std::cout << "events_below_bins: " << m_BelowLowestBin << std::endl;
            t_INPUT lower( m_Min), upper( m_Min + m_Delta);
            STREAM << "#sum of values found below lowest bin: " << m_BelowLowestBin << std::endl;
            for( std::vector< size_t>::const_iterator itr = m_Data.begin(); itr != m_Data.end(); ++itr)
            {
                STREAM << lower + 0.5 * m_Delta << "  " << *itr << std::endl;
                lower += m_Delta;
                upper += m_Delta;
            }
            STREAM << "#sum of values found above hightest bin: " << m_AboveHighestBin << std::endl << "# ";
	    std::copy( m_Above.begin(), m_Above.end(), std::ostream_iterator< t_INPUT>( STREAM, " "));
	    STREAM << std::endl << "# ";
	    std::copy( m_Below.begin(), m_Below.end(), std::ostream_iterator< t_INPUT>( STREAM, " "));
	    STREAM << std::endl;
            return STREAM;
        }

    }; // end class Histogram
//} // end namespace math




#endif /* HISTOGRAM_T_H_ */
