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



#include <cassert>
#include "../include/vector_functions.h"
#include "../include/math_functions.h"

namespace math
{
  Vector
  NormalizedVector( const Vector &V)
  {
    Vector v( V);
    return v.Normalize();
  }

  double
  ScalarProduct( const Vector &V1, const Vector &V2)
  {
    assert( V1.size() == V2.size());
    double scalar(0.0);
    for( std::vector< double>::const_iterator itr1( V1.begin()), itr2( V2.begin()); itr1 != V1.end() && itr2 != V2.end(); ++itr1, ++itr2)
      { scalar += ( *itr1) * ( *itr2);}
    return scalar;
  }

  std::ostream &
  operator << ( std::ostream & STREAM, const Vector& V)
  { return V.Write( STREAM);}

  std::istream &
  operator >> ( std::istream &STREAM, Vector &V)
  { return V.Read( STREAM);}

  Vector
  operator - ( const Vector &V1, const Vector &V2)
    {
      assert( V1.size() == V2.size());
      Vector difference( V1.size());
      std::vector< double>::iterator itr_diff( difference.begin());
      for( std::vector<double>::const_iterator itr_v1( V1.begin()), itr_v2( V2.begin()); itr_v1 != V1.end(); ++itr_v1, ++itr_v2, ++itr_diff)
      { *itr_diff = *itr_v1 - *itr_v2;}
      return difference;
  }

  Vector
  operator + ( const Vector &V1, const Vector &V2)
    {
      assert( V1.size() == V2.size());
      Vector sum( V1.size());
      std::vector< double>::iterator itr_sum( sum.begin());
      for( std::vector<double>::const_iterator itr_v1( V1.begin()), itr_v2( V2.begin()); itr_v1 != V1.end(); ++itr_v1, ++itr_v2, ++itr_sum)
      { *itr_sum = *itr_v1 + *itr_v2;}
      return sum;
  }

  double
  operator * ( const Vector &V1, const Vector &V2)
    {
      assert( V1.size() == V2.size());
      double scalar( 0.0);
      for( std::vector< double>::const_iterator itr_v1( V1.begin()), itr_v2( V2.begin()); itr_v1 != V1.end(); ++itr_v1, ++itr_v2)
      { scalar += *itr_v1 * *itr_v2;}
      return scalar;
    }

  // multiply vector with factor
  Vector
  operator * ( const double &VALUE, const Vector &V)
  {
      Vector vec( V);
      for( std::vector< double>::iterator itr = vec.begin(); itr != vec.end(); ++itr)
      {
          *itr *= VALUE;
      }
      return vec;
  }

  // multiply vector with factor
  Vector
  operator * ( const Vector &V, const double &VALUE)
{
      Vector vec( V);
      for( std::vector< double>::iterator itr = vec.begin(); itr != vec.end(); ++itr)
      {
          *itr *= VALUE;
      }
      return vec;
}

  Vector
  operator / ( const Vector &V, const double &VALUE)
  {
	  Vector vec( V);
      for( std::vector< double>::iterator itr = vec.begin(); itr != vec.end(); ++itr)
      {
          *itr /= VALUE;
      }
      return vec;
  }


  Vector
  operator / ( const Vector &V1, const Vector &V2)
  {
      assert( V1.size() == V2.size());
      Vector vec( V1.size());
      std::vector< double>::iterator itr = vec.begin();
      for( std::vector< double>::const_iterator itr_v1( V1.begin()), itr_v2( V2.begin()); itr_v1 != V1.end(); ++itr_v1, ++itr_v2, ++itr)
      { *itr = *itr_v1 / *itr_v2;}
      return vec;
  }

  double
  SquaredDistance( const Vector &POSITION_A, const Vector &POSITION_B)
  {
    assert( POSITION_A.size() == POSITION_B.size());
    double dist( 0.0);
    for( std::vector< double>::const_iterator itr_a( POSITION_A.begin()), itr_b( POSITION_B.begin()); itr_a != POSITION_A.end() && itr_b != POSITION_B.end(); ++itr_a, ++itr_b)
      { dist += ( *itr_b - *itr_a) * ( *itr_b - *itr_a);}
    return dist;
  }

  double
  Distance( const Vector &POSITION_A, const Vector &POSITION_B)
  {
#ifdef DEBUG
	  if( POSITION_A.size() != POSITION_B.size())
	  {
    	std::cout << "ERROR: Distance Vector sizes do not match, v1: " << POSITION_A << std::flush << " v2: " << POSITION_B << std::endl;
    	throw 1;
//    	exit(1);
      }
#endif

    double dist( 0.0);
    for( std::vector< double>::const_iterator itr_a( POSITION_A.begin()), itr_b( POSITION_B.begin()); itr_a != POSITION_A.end() && itr_b != POSITION_B.end(); ++itr_a, ++itr_b)
      { dist += ( *itr_b - *itr_a) * ( *itr_b - *itr_a);}

#ifdef DEBUG
    if( dist == 0.0)
    {
    	std::cout << "ERROR: dist 0" << std::endl;
    	throw 2;
    }
#endif

    return sqrt( dist);
  }

  double
  Angle( const Vector &DIRECTION_A, const Vector &DIRECTION_B)
  {
      double 
	  a_sqr =  DIRECTION_A.SquaredLength(),
	  b_sqr =  DIRECTION_B.SquaredLength();

      if( a_sqr == 0 || b_sqr == 0)
      {
          std::cout << "====> angle of vectors of size zero requested!" << std::endl;
          return std::numeric_limits< double>::max(); // TODO: return undefined
      }
      double scalar( ScalarProduct( DIRECTION_A, DIRECTION_B) / sqrt( a_sqr * b_sqr));
      if( scalar > 1.0 && scalar < 1.00001)
      {
    	  scalar = 1.0;
      }
      else if( scalar < -1.0 && scalar > -1.00001)
      {
          scalar = -1.0;
      }
      else if( scalar > 1.0 || scalar < -1.0)
      {
          std::cout << "===> undefined values in scalar product: " << scalar << std::endl;
          exit( -1);
      }
    return acos( scalar);
  }


  bool
  IsDistanceSmallerThan( const Vector &POSITION_A, const Vector &POSITION_B, const double &THRESHOLD)
  {
        assert( POSITION_A.size() == POSITION_B.size());
        double squared_distance( 0.0);
        for( std::vector< double>::const_iterator itr_a( POSITION_A.begin()), itr_b( POSITION_B.begin()); itr_a != POSITION_A.end() && itr_b != POSITION_B.end(); ++itr_a, ++itr_b)
          { squared_distance += ( *itr_b - *itr_a) * ( *itr_b - *itr_a);}
        return squared_distance < THRESHOLD * THRESHOLD;
  }

  // the length of the projection of vector A on vector B
  double
  ProjectionOnVector( const math::Vector &A, const math::Vector &B)
  {
       return ScalarProduct( A, B) / B.Length();
  }

  // returns elementwise minimum between two vectors A and B
   Vector
   Min( const Vector &A, const Vector &B)
   {
	      assert( A.size() == B.size());
	      Vector vec( A.size());
	      std::vector< double>::iterator itr = vec.begin();
	      for( std::vector< double>::const_iterator itr_v1( A.begin()), itr_v2( B.begin()); itr_v1 != A.end(); ++itr_v1, ++itr_v2, ++itr)
	      { *itr = std::min( *itr_v1, *itr_v2);}
	      return vec;
   }


    // returns elementwise maximum between two vectors A and B
    Vector
    Max( const Vector &A, const Vector &B)
    {
 	      assert( A.size() == B.size());
 	      Vector vec( A.size());
 	      std::vector< double>::iterator itr = vec.begin();
 	      for( std::vector< double>::const_iterator itr_v1( A.begin()), itr_v2( B.begin()); itr_v1 != A.end(); ++itr_v1, ++itr_v2, ++itr)
 	      { *itr = std::max( *itr_v1, *itr_v2);}
 	      return vec;
    }


	Vector
	Randomize( const Vector &V, const double &MAX)
	{
		Vector
			v = V;
		for( std::vector< double>::iterator itr = v.begin(); itr != v.end(); ++itr)
		{
			itr += Random( -MAX, MAX);
		}
		return v;
	}

	Vector
	RandomVector( int SIZE)
	{
		Vector v( SIZE);
		for( std::vector< double>::iterator itr = v.begin(); itr != v.end(); ++itr)
		{
			*itr = Random( -1.0, 1.0);
		}
		return v;
	}

	Vector
	RandomVector( const double &LENGTH, int SIZE)
	{
		Vector v( SIZE);
		while( (v = RandomVector(SIZE)).Length() > 1.0)
		{}
		return v.SetToLength( LENGTH);
	}


} // end namespace math

