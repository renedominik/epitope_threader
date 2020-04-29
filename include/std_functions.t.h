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



#ifndef STD_FUNCTIONS_H_
#define STD_FUNCTIONS_H_


#include "global.h"

#include <vector>
#include <map>
#include <set>
#include <stdlib.h>

template< typename T>
std::ostream &
WriteVector( const std::vector< T> &V, std::ostream &OUT = std::cout, const char &c = ' ')
{
	DebugFct;
	OUT << V.size() << std::endl;
	for( unsigned int i = 0; i < V.size(); ++i)
	{
		OUT << V[i] << c;
//		V[i].Write( OUT);
	}
	return OUT;
}

template< typename T>
std::ostream&
operator << ( std::ostream &OUT, const std::vector< T> &V)
{
	return WriteVector( V, OUT);
}

template< typename T>
std::ostream&
WriteVecVec( const std::vector< std::vector< T> > &V, std::ostream &OUT)
{
	OUT << V.size() << std::endl << std::endl;
	for( size_t i = 0; i < V.size(); ++i)
	{
		WriteVector( V[i], OUT);
		OUT << std::endl;
	}
	return OUT;
}

template< typename T>
std::ostream&
operator << ( std::ostream &OUT, const std::vector< std::vector< T> > &V)
{
	return WriteVecVec( V, OUT);
}


template< typename TA, typename TB>
std::ostream &
operator << ( std::ostream &OUT, const std::pair< TA, TB> &P)
{
	OUT << P.first << " " << P.second << std::endl;
	return OUT;
}

template< typename T>
std::vector< T>
ReadVector( std::istream &IN)
{
	DebugFct;
	int size;
	IN >> size;
	std::vector< T> v( size);
	for( int i = 0; i < size; ++i)
	{
		IN >> v[i];
	}
	return v;
}


template< class t_KEY, class t_VALUE>
const t_VALUE &
SafeGet( const std::map< t_KEY, t_VALUE> &MAP, const t_KEY &KEY)
{
	typename std::map< t_KEY, t_VALUE>::const_iterator itr = MAP.find( KEY);
	if( itr == MAP.end())
	{
		std::cerr << "ERROR: " << KEY << " not found in " << __FUNCTION__ << "\n";
		exit(1);
	}
	return itr->second;
}

template < class t_KEY, class t_VALUE>
bool
CheckedGet( const std::map<t_KEY,t_VALUE> &MAP, const t_KEY &KEY, t_VALUE &RETURNVAL)
{
    //    std::cout << __FUNCTION__ << " " << KEY << std::endl;
	typename std::map< t_KEY, t_VALUE>::const_iterator itr = MAP.find( KEY);
	if( itr == MAP.end())
	{
	  //		std::cerr << "ERROR: " << KEY << " not found in " << __FUNCTION__ << "\n";
		return false;
	}
	//    std::cout << __FUNCTION__ << " found: " << KEY << " " << itr->second << std::endl;
	
	RETURNVAL = itr->second;
	return true;
}
    


// forward declaration
template< typename t_VALUE>
class Histogram;

template< class t_KEY, class t_VALUE>
void
SafeAdd( std::map< t_KEY, Histogram< t_VALUE> > &MAP, const t_KEY &KEY, const t_VALUE &VALUE, const Histogram< t_VALUE> &DEFAULT)
{
	typename std::map< t_KEY, Histogram< t_VALUE> >::iterator itr = MAP.find( KEY);
	if( itr == MAP.end())
	{
		MAP.insert( std::make_pair( KEY, DEFAULT));
		itr = MAP.find( KEY);
	}
	itr->second(VALUE);
}

template< class t_KEY, class t_VALUE>
void
SafeAdd( std::map< t_KEY, std::vector< t_VALUE> > &MAP, const t_KEY &KEY, const t_VALUE &VALUE)
{
	typename std::map< t_KEY, std::vector< t_VALUE> >::iterator itr = MAP.find( KEY);
	if( itr == MAP.end())
	{
		MAP[ KEY] = std::vector< t_VALUE>(1, VALUE);
	}
	else
	{
		itr->second.push_back( VALUE);
	}
}

template< class t_KEY, class t_VALUE>
std::ostream &
WriteMap( const std::map< t_KEY, t_VALUE> &MAP, std::ostream &OUT)
{
	DebugFct;
	typename std::map< t_KEY, t_VALUE>::const_iterator itr = MAP.begin();
	for( ;itr != MAP.end(); ++itr)
	{
		OUT << itr->first << "   " << itr->second << "\n";
	}
	return OUT;
}

template< class t_KEY, class t_VALUE>
std::map< t_KEY, t_VALUE>
ReadMap( std::istream &IN)
{
	DebugFct;
	std::map< t_KEY, t_VALUE> map;
	t_KEY key;
	t_VALUE value;
	while( IN >> key >> value)
	{
//		std::cout << __FUNCTION__ <<" " << key << " " << value << std::endl;
		map[key] = value;
	}
	return map;
}

template< class t_KEY, class t_VALUE>
std::istream &
ReadMap( std::map< t_KEY, t_VALUE> &MAP, std::istream &IN)
{
	DebugFct;
	t_KEY key;
	t_VALUE value;
	while( IN >> key >> value)
	{
//		std::cout << __FUNCTION__ << key << " " << value << std::endl;
		MAP[key] = value;
	}
	return IN;
}

template< class t_KEY, class t_VALUE>
std::ostream &
WriteMultimap( const std::multimap< t_KEY, t_VALUE> &MAP, std::ostream &OUT)
{
	DebugFct;
	typename std::multimap< t_KEY, t_VALUE>::const_iterator itr = MAP.begin();
	for( ;itr != MAP.end(); ++itr)
	{
		OUT << itr->first << ":   " << itr->second << "\n";
	}
	return OUT;
}


template< class t_KEY, class t_VALUE>
const void
UniqueInsert(  std::map< t_KEY, t_VALUE> &MAP, const t_KEY &KEY, const t_VALUE &VAL)
{
	typename std::map< t_KEY, t_VALUE>::const_iterator itr = MAP.find( KEY);
	if( itr == MAP.end())
	{
		MAP[KEY] = VAL;
	}
	else
	{
		std::cerr << "WARNING: key (" << KEY << ") already exists (" << KEY << "," << itr->second << "), new value (" << VAL << ") skipped, message from " << __FUNCTION__ << "\n";
	}
}


template< typename t_KEY, typename t_VALUE>
std::istream&
operator >> ( std::istream &IN, std::map< t_KEY, t_VALUE> &M)
{
	return ReadMap( M, IN);
}


template< typename t_KEY, typename t_VALUE>
std::ostream&
operator << ( std::ostream &OUT, const std::map< t_KEY, t_VALUE> &V)
{
	return WriteMap( V, OUT);
}


template< typename t_KEY, typename t_VALUE>
std::ostream&
operator << ( std::ostream &OUT, const std::multimap< t_KEY, t_VALUE> &V)
{
	return WriteMultimap( V, OUT);
}


template< typename t_VALUE>
t_VALUE RandomElement( const std::vector< t_VALUE> &V)
{
	return V[ Random( V.size())];
}


template< typename t_VALUE>
bool HaveCommonElements( const std::vector< t_VALUE> &V1, const std::vector< t_VALUE> &V2)
{
	for( typename std::vector< t_VALUE>::const_iterator itr = V1.begin(); itr != V1.end(); ++itr)
		if( std::find( V2.begin(), V2.end(), *itr) != V2.end())
		{
			return true;
		}
	return false;
}

template< typename T>
int
RemoveDuplicatesUnsorted(std::vector<T>& V)
{
	int b = V.size();
	std::set<T> seen;

	typename std::vector<T>::iterator itr = V.begin();
	while(itr != V.end())
    {
		if( seen.find(*itr) != seen.end())
		{
			itr = V.erase(itr);
		}
		else
        {
			seen.insert(*itr);
			itr++;
        }
    }
	return b - seen.size();
}

template< typename T>
std::vector<T>& Add( std::vector<T> &V, int ID, int ADD)
{
	V[ID] += ADD;
	return V;
}



#endif /* STD_FUNCTIONS_H_ */
