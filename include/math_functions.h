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


#ifndef MATH_DOUBLE_FUNCTIONS_H
#define MATH_DOUBLE_FUNCTIONS_H

// access to M_PI from math.h
#define _USE_MATH_DEFINES

#include <math.h>
#include <cstdlib>

template <typename T1, typename T2, typename T3>
struct Triplet
{
	T1 first;
	T2 second;
	T3 third;

	Triplet( const T1 &FIRST, const T2 &SECOND, const T3 &THIRD)
	: first( FIRST), second( SECOND), third( THIRD)
	{}

};



namespace math
{


  inline double Unsign( const double &X)
  {
    if( X < 0){ return -X;}
    return X;
  }

  // rounds positive (!) double values, returns an size_t
  template< typename t_TYPE>
  inline t_TYPE Round( const double &X)
  {
    if( X < 0.0)
    { return t_TYPE( int( X - 0.5));}
    return t_TYPE( X + 0.5);
  }

  // rounds double values, returns a double
  inline double Round( const double &X, size_t N = 0)
  {
	  double n = pow( 10.0, N);
	  if( X < 0.0)
	  { return double( int( n * X - 0.5)) / n;}
	  return double( int( n * X + 0.5)) / n;
  }


  // truncates a double value, returns an double
  inline double Truncate( const double &X)
  { return double( int( X));}

  // random [MIN,MAX)
  inline double Random( const double &MIN, const double &MAX)
  {
      return double( rand() / ( double( RAND_MAX) + 1.0) * ( MAX - MIN) + MIN);
  }

  // random (0, 1)
  inline double Random()
  {
	  return double( rand()) / double( RAND_MAX);
  }

  // random (0, max)
  inline double Random( const double &MAX)
  {
	  return double( rand()) / double( RAND_MAX) * MAX;
  }

  inline
  double Square( const double &VALUE)
    { return VALUE * VALUE;}

  inline
  double AngleFromLawOfCosinus( const double &OPPOSITE_SIDE_LENGTH, const double &FIRST_NEXT_SIDE_LENGTH, const double &SECOND_NEXT_SIDE_LENGTH)
  { return acos( ( FIRST_NEXT_SIDE_LENGTH * FIRST_NEXT_SIDE_LENGTH + SECOND_NEXT_SIDE_LENGTH * SECOND_NEXT_SIDE_LENGTH - OPPOSITE_SIDE_LENGTH * OPPOSITE_SIDE_LENGTH) / ( 2.0 * FIRST_NEXT_SIDE_LENGTH * SECOND_NEXT_SIDE_LENGTH));}

  inline
  double CosinusAngleFromLawOfCosinus( const double &OPPOSITE_SIDE_LENGTH, const double &FIRST_NEXT_SIDE_LENGTH, const double &SECOND_NEXT_SIDE_LENGTH)
  { return ( FIRST_NEXT_SIDE_LENGTH * FIRST_NEXT_SIDE_LENGTH + SECOND_NEXT_SIDE_LENGTH * SECOND_NEXT_SIDE_LENGTH - OPPOSITE_SIDE_LENGTH * OPPOSITE_SIDE_LENGTH) / ( 2.0 * FIRST_NEXT_SIDE_LENGTH * SECOND_NEXT_SIDE_LENGTH);}

  inline
  double LawOfCosinus( const double &AN_1, const double &AN_2, const double &ANGLE)
  { return ( AN_1 * AN_1 + AN_2 * AN_2 - 2.0 * AN_1 * AN_2 * cos( ANGLE));}

  inline
  double ProjectedLawOfCosinus( const double &FIRST_SIDE_LENGTH, const double &FIRST_ANGLE, const double &SECOND_SIDE_LENGTH, const double &SECOND_ANGLE)
  { return FIRST_SIDE_LENGTH * cos( SECOND_ANGLE) + SECOND_SIDE_LENGTH * cos( FIRST_ANGLE);}

  inline
  double LengthFromLawOfSinus( const double &OPPOSITE_ANGLE, const double &OTHER_SIDE_LENGTH, const double &OTHER_OPPOSITE_ANGLE)
  { return OTHER_SIDE_LENGTH * sin( OPPOSITE_ANGLE) / sin( OTHER_OPPOSITE_ANGLE);}

  inline
  double AngleFromLawOfSinus( const double & OPPOSITE_SIDE_LENGTH, const double &OTHER_SIDE_LENGTH, const double &OTHER_ANGLE)
  { return asin( OPPOSITE_SIDE_LENGTH / OTHER_SIDE_LENGTH * sin( OTHER_ANGLE));}

  inline
  double RadiansToDegrees( const double &RAD)
  { return RAD * 180.0 / M_PI;}

  inline
  double DegreesToRadians( const double &DEG)
  { return DEG * M_PI / 180.0;}

  inline
  bool IsLarger( const double& X, const double &Y)
  { return X > Y;}

  inline
  bool IsSmaller( const double& X, const double &Y)
  { return X < Y;}

  inline
  bool IsEqual( const double& X, const double &Y)
  { return X == Y;}

  inline
  bool IsEqualWithinThreshold( const double& X, const double &Y, const double &THRESHOLD)
  { return fabs( X - Y) < THRESHOLD;}

  inline
  bool IsNonEqual( const double& X, const double &Y)
  { return X != Y;}

  inline
  double
  Modulo( const double &VALUE, const double &NOMINATOR)
  {
    return VALUE - double( int( VALUE / NOMINATOR)) * NOMINATOR;
  }

  inline
  double
  GaussDistribution( const double &X, const double &DELTA)
  {
	  return  1.0 / ( DELTA * sqrt( 2.0 * M_PI)) * exp( -0.5 * ( X * X / ( DELTA * DELTA)));
  }

  inline
  double
  StandardDeviation( const std::vector< double> &V, double X)
  {
	  double dev = 0;
	  for( std::vector< double>::const_iterator itr = V.begin(); itr != V.end(); ++itr)
	  {
		  dev += Square( X - *itr);
	  }
	  return sqrt( dev / (double) V.size());
  }

  inline
  double
  Angle2Radian( double ANGLE)
  {
	  return M_PI * ANGLE / 180.0;
  }

  inline
  double
  Radian2Angle( double RAD)
  {
	  return RAD * 180.0 / M_PI;
  }

} // end namespace math


#endif
