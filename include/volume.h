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



#ifndef VOLUME_H_
#define VOLUME_H_

#include<set>
//#include <iterator>
#include "std_functions.t.h"

typedef int PointID;


class VolumeID
{
private:
	std::vector< PointID>	  m_Points;
	std::vector< PointID>	  m_Surf;
//	std::vector< PointID>	  m_Cover;  // points outside, touching surface
//	PointID                   m_Min;     // limits of cube around volume
//	PointID                   m_Max;

public:

	VolumeID()
	: m_Points(),
	m_Surf()
	{}

	VolumeID( const VolumeID &V)
	: m_Points( V.m_Points),
	m_Surf( V.m_Surf)
	{/* static int count = 0; std::cout << __FUNCTION__ << " copy constructor " << count << std::endl; ++count;*/}


	VolumeID&
	Add( const PointID &P)
	{
		m_Points.push_back( P);
		return *this;
	}


	VolumeID&
	AddSurf( const PointID &P)
	{
		m_Surf.push_back( P);
		return *this;
	}


	VolumeID&
	MoveLastToSurf()
	{
		m_Surf.push_back( m_Points.back());
		m_Points.pop_back();
		return *this;
	}

	int
	NrBins()
	{
		return m_Points.size() + m_Surf.size();
	}

	std::vector< PointID>&
	GetPoints()
	{
		return m_Points;
	}

	const std::vector< PointID>&
	GetPoints() const
	{
		return m_Points;
	}

	std::vector< PointID>&
	GetSurf()
	{
		return m_Surf;
	}

	const std::vector< PointID>&
	GetSurf() const
	{
		return m_Surf;
	}

	std::vector< PointID>::iterator
	Find( const PointID &P)
	{
		return std::find( m_Points.begin(), m_Points.end(), P);
	}

	std::vector< PointID>::iterator
	FindSurf( const PointID &P)
	{
		return std::find( m_Surf.begin(), m_Surf.end(), P);
	}


	void
	RemoveSurf( const PointID &P)
	{
//		std::cout << __FUNCTION__ << "  " << m_Surf.size() << std::endl;
//		std::cout << " p " << P << std::endl;
		std::vector< PointID>::iterator itr = std::find( m_Surf.begin(), m_Surf.end(), P);
//		std::cout << __FUNCTION__ << " found " << itr - m_Surf.begin() << " " << m_Surf.end() - m_Surf.begin() << std::endl;
		if( itr != m_Surf.end())
		{
			std::cout << __FUNCTION__ << " erase  " << std::endl;
			m_Surf.erase( itr);
		}
//		std::cout << __FUNCTION__ << "  done " << m_Surf.size() << std::endl;
	}

	void
	RemoveSurf( std::vector< PointID>::iterator PTR)
	{
//		std::cout << __PRETTY_FUNCTION__ << "  " << m_Surf.size() << std::endl;
//		std::cout << " p " << PTR->Data() << std::endl;
//		std::cout << " p " << *PTR << std::endl;
		if( PTR != m_Surf.end())
		{
			std::cout << __FUNCTION__ << " erase  " << std::endl;
			m_Surf.erase( PTR);
		}
//		std::cout << __FUNCTION__ << "  done " << m_Surf.size() << std::endl;
	}


	void
	Remove( const PointID &P)
	{
		std::vector< PointID>::iterator itr = std::find( m_Points.begin(), m_Points.end(), P);
		if( itr != m_Points.end())
		{
			m_Points.erase( itr);
		}
	}


	bool
	IsVol( const PointID &P)
	{
		return  std::find( m_Points.begin(), m_Points.end(), P) != m_Points.end();
	}

	bool
	IsSurf( const PointID &P)
	{
		return  std::find( m_Surf.begin(), m_Surf.end(), P) != m_Surf.end();
	}

	PointID&
	Last()
	{
		return m_Points.back();
	}



	std::ostream &Write( std::ostream &OUT) const
	{
		OUT << "points:\n" <<m_Points;
		OUT << "surf:\n" <<m_Surf;
		return OUT;
	}

	VolumeID&
	Insert( const VolumeID &VOL)
	{
		m_Points.insert( m_Points.end(), VOL.m_Points.begin(), VOL.m_Points.end());
		return *this;
	}

};  // end class VolumeID


std::ostream& operator << ( std::ostream &OUT, const VolumeID &P)
{
	return P.Write( OUT);
}


bool
AreOverlapping( const VolumeID &V1, const VolumeID &V2)
{
	for( std::vector< PointID>::const_iterator i1 = V1.GetPoints().begin(); i1 != V1.GetPoints().end(); ++i1)
		for( std::vector< PointID>::const_iterator i2 = V2.GetPoints().begin(); i2 != V2.GetPoints().end(); ++i2)
			if( *i1 == *i2)
			{
				return true;
			}
	return false;
}

bool
MergeIfOverlapping( VolumeID &V1, const VolumeID &V2)
{
		int
			a,b,c,size;

//			std::vector< PointID> intersection;
//			std::set_intersection( atr->GetPoints().begin(), atr->GetPoints().end(), btr->GetPoints().end(), btr->GetPoints().end(), std::back_inserter( intersection));  // wasting time here!!!

		if( HaveCommonElements( V1.GetPoints(), V2.GetPoints()))
		{
			std::cout << __FUNCTION__ << " have common elements" << std::endl;
			a = V1.GetPoints().size();
			b = V2.GetPoints().size();
			V1.GetPoints().insert( V1.GetPoints().end(), V2.GetPoints().begin(), V2.GetPoints().end());
			c = V1.GetPoints().size();
			size = RemoveDuplicatesUnsorted( V1.GetPoints());
			if( size > 0)
			{
				std::cout << " removed points " <<size  << " (" << a << ":" << b << "=" << c << " vs " << V1.GetPoints().size() << ")\n";
				flush( std::cout);
			}
			std::cout << std::endl;
			return true;
		}
		return false;
}

// returns a list of ids of points surrounding a given volume (not its surface, but its cover, points are outside of volume)
std::vector< PointID>
FindCoverPointIDs( const VolumeID &VOL, const std::vector< int>& OFFSET)
{
//	std::cout << __FUNCTION__ << std::endl;
	std::set< PointID>
		neighbors;

	for( std::vector< PointID>::const_iterator itr = VOL.GetPoints().begin(); itr != VOL.GetPoints().end(); ++itr)
		for( std::vector< int>::const_iterator off = OFFSET.begin(); off != OFFSET.end(); ++off)
			if( std::find( VOL.GetPoints().begin(), VOL.GetPoints().end(), *itr + *off) == VOL.GetPoints().end())
			{
				neighbors.insert( *itr + *off);
			}
//	std::cout << __FUNCTION__ <<  " " << neighbors.size() << std::endl;

	return std::vector< PointID>(neighbors.begin(), neighbors.end());
}



bool
AreNeighbors( const std::vector< int> &COVER, const VolumeID &VOL, int MIN)
{
//	std::cout << __FUNCTION__ << " " << NEIGHBORS.size() << " " << VOL.GetPoints().size() << std::endl;
	int
		contacts = 0;

	for( std::vector< PointID>::const_iterator itr = COVER.begin(); itr != COVER.end(); ++itr)
		if( std::find( VOL.GetPoints().begin(), VOL.GetPoints().end(), *itr) != VOL.GetPoints().end() && ++contacts >= MIN) // iterates from beginning to end each time, actually unnecessary
		{
			std::cout << __FUNCTION__ << " true" << std::endl;
			return true;
		}

//	std::cout << __FUNCTION__ << " false, #contacts: " << contacts << std::endl;
	return false;
}




class VolumeIDArray
{
private:
	std::vector< VolumeID>      m_Volumes;
public:

	VolumeID&
	Last()
	{
		return m_Volumes.back();
	}

	VolumeIDArray&
	Add( const VolumeID &V)
	{
		if( V.GetPoints().size() == 0 && V.GetSurf().size() == 0)
		{
			std::cout << "WARNING: " << __FUNCTION__ << " empty" << std::endl;
			return *this;
		}

//		std::cout << __FUNCTION__ << " vol" << std::endl;
		m_Volumes.push_back( V);
		return *this;
	}

	std::vector< VolumeID>&
	Data()
	{
		return m_Volumes;
	}

	const std::vector< VolumeID>&
	Data() const
	{
		return m_Volumes;
	}

	std::vector< VolumeID>::iterator
	FindVolumeContainingPoint( const PointID &P)
	{
		for( std::vector< VolumeID>::iterator vtr = m_Volumes.begin(); vtr != m_Volumes.end(); ++vtr) // reverse?
		{
			if( vtr->Find( P) != vtr->GetPoints().end())  // find point but return vol ... somewhat weird ...
			{
				return vtr;
			}
		}
		return m_Volumes.end();
	}



	void
	FuseLast( const PointID &POINT)
	{
		VolumeID
			last = m_Volumes.back();
		std::vector< VolumeID>::iterator
			vol = FindVolumeContainingPoint( POINT);

		if( vol == m_Volumes.end())
		{
			return;
		}
		if( &*vol == &last)
		{
			return;
		}
		last.GetPoints().insert( last.GetPoints().end(), vol->GetPoints().begin(), vol->GetPoints().end());

		RemoveDuplicatesUnsorted(last.GetPoints());

		m_Volumes.back() = last;

		m_Volumes.erase( vol);
	}



	void
	RemoveDuplicates()
	{
		unsigned int
			size;
		std::cout << __FUNCTION__ ;
		for( std::vector< VolumeID>::iterator atr = m_Volumes.begin(); atr != m_Volumes.end(); ++atr)
		{
			size = atr->GetPoints().size();
			RemoveDuplicatesUnsorted(atr->GetPoints());
			if( size != atr->GetPoints().size())
			{
				std::cout << " removed duplicates (" << size << ":" << atr->GetPoints().size() << ") ";
			}
		}
		std::cout << std::endl;
	}


	void
	FuseOverlapping()
	{
		std::cout << __FUNCTION__ << " before: " << m_Volumes.size();
		flush( std::cout);
		for( std::vector< VolumeID>::iterator atr = m_Volumes.begin(); atr + 1 < m_Volumes.end(); ++atr)
			for( std::vector< VolumeID>::iterator btr = atr + 1; btr < m_Volumes.end(); ++btr)
			{
				if( MergeIfOverlapping( *atr, *btr))
				{
					btr = m_Volumes.erase( btr);
					--btr;
				}
			}
		std::cout << " after: " << m_Volumes.size() << std::endl;
	}

	void
	FuseNeighbors( int MIN, const std::vector< int>& OFFSET)
	{
		std::cout << __FUNCTION__ << " #vols: " << m_Volumes.size()  << " min contacts: " << MIN << std::endl;
		int before, size;
		for( std::vector< VolumeID>::iterator atr = m_Volumes.begin(); atr + 1 < m_Volumes.end(); ++atr)
		{
			std::cout << m_Volumes.size() - (atr - m_Volumes.begin()) << ":" << atr->GetPoints().size() << ":" ;
			flush( std::cout);
			std::vector<int>
				neighbors = FindCoverPointIDs( *atr, OFFSET);
			std::cout << neighbors.size() << "  ";
			flush( std::cout);
			for( std::vector< VolumeID>::iterator btr = atr + 1; btr < m_Volumes.end(); ++btr)
				if( AreNeighbors( neighbors, *btr, MIN))
				{
					before = atr->GetPoints().size();
					atr->GetPoints().insert( atr->GetPoints().end(), btr->GetPoints().begin(), btr->GetPoints().end());
					size = RemoveDuplicatesUnsorted( atr->GetPoints());
					btr = m_Volumes.erase( btr);
					--btr;
				}
		}
		std::cout << "\n" << std::endl;
	}



	void
	Sort()
	{
		for( std::vector< VolumeID>::iterator atr = m_Volumes.begin(); atr + 1 < m_Volumes.end(); ++atr)
		{
			std::sort( atr->GetPoints().begin(), atr->GetPoints().end());
		}
	}

};  // end class VolumeIDArray


void
Insert( std::vector< std::set< const VolumeID*> > &GROUPS, const VolumeID *FIRST, const VolumeID *SECOND)
{
	std::cout << __FUNCTION__  << std::endl;
	bool
		brand_new = true;
	for( std::vector< std::set< const VolumeID*> >::iterator itr = GROUPS.begin(); itr != GROUPS.end() && brand_new; ++itr)
		if( itr->find( FIRST) != itr->end() || itr->find( SECOND) != itr->end())
		{
			itr->insert( FIRST);
			itr->insert( SECOND);
			brand_new = false;
		}
	if( brand_new)
	{
		std::set< const VolumeID*>
			set;
		set.insert( FIRST);
		set.insert( SECOND);
		GROUPS.push_back( set);
	}
}


// keep all HIGHs, merge with LOWS
VolumeIDArray
MergeOverlappingPairs( const VolumeIDArray &HIGH, const VolumeIDArray &LOW)
{
	std::cout << __FUNCTION__  << std::endl;
	VolumeIDArray
		merged;
	std::vector< std::set< const VolumeID*> >
		groups;

	std::vector< VolumeID>::iterator
		htr;
	bool
		has_overlap;
	// only one point of high needs to be within low
	for( std::vector< VolumeID>::const_iterator high = HIGH.Data().begin(); high != HIGH.Data().end(); ++high)
	{
		std::cout << __FUNCTION__ << " check high for overlaps with low"  << std::endl;
		has_overlap = false;

		for( std::vector< VolumeID>::const_iterator low = LOW.Data().begin(); low != LOW.Data().end(); ++low)
			if( AreOverlapping( *high, *low))
			{
				std::cout << __FUNCTION__ << " overlap"  << std::endl;
				has_overlap = true;
				Insert( groups, high.base(), low.base());
			}

		if( !has_overlap)
		{
			std::cout << __FUNCTION__ << " new"  << std::endl;
			std::set< const VolumeID*>
				set;
			set.insert( high.base());
			groups.push_back( set);
		}
	}

	// add groups to merged
	for( std::vector< std::set< const VolumeID*> >::const_iterator itr = groups.begin(); itr != groups.end(); ++itr)
	{
		VolumeID
			vol; // = **itr->begin();
		for( std::set< const VolumeID*>::const_iterator vtr = itr->begin(); vtr != itr->end(); ++vtr)
		{
			vol.Insert( **vtr);
		}
		merged.Add( vol);
	}

	merged.RemoveDuplicates();
	std::cout << __FUNCTION__ << " done"  << std::endl;
	return merged;
}

#endif /* VOLUME_H_ */
