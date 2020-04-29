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


#ifndef VECTOR_FUNCTIONS_H
#define VECTOR_FUNCTIONS_H


#include "vector.h"

namespace math
{
  // returns normalized vector, conserves argument vector
  Vector NormalizedVector( const Vector &V);


  // returns the scalar product of two vectors V1 and V2
  double ScalarProduct( const Vector &V1, const Vector &V2);

  // write vector into std::ostream
  std::ostream &operator << ( std::ostream & STREAM, const Vector& V);

  // read vector from std::istream
  std::istream &operator >> ( std::istream &STREAM, Vector &V);

  // difference of two vectors
  Vector operator - ( const Vector &V1, const Vector &V2);

  // sum of two vectors
  Vector operator + ( const Vector &V1, const Vector &V2);

  // scalar product of two vectors
  double operator * ( const Vector &V1, const Vector &V2);

  // multiply vector with factor
  Vector operator * ( const double &VALUE, const Vector &V);

  // multiply vector with factor
  Vector operator * ( const Vector &V, const double &VALUE);

  // devide vector by factor
  Vector operator / ( const Vector &V, const double &VALUE);

  // devide vector by vector
  Vector operator / ( const Vector &V1, const Vector &V2);

  // the squared distance of two points
  double SquaredDistance( const Vector &POS_A, const Vector &POS_B);

  // the distance of two points
   double Distance( const Vector &POS_A, const Vector &POS_B);

   // the angle between two vectors
   double Angle( const Vector &DIRECTION_A, const Vector &DIRECTION_B);

   // boolean function checking whether the distance of two points is smaller than threshold
   bool IsDistanceSmallerThan( const Vector &POSITION_A, const Vector &POSITION_B, const double &THRESHOLD);

   // the length of the projection of vector A on vector B
   double ProjectionOnVector( const math::Vector &A, const math::Vector &B);

   // returns elementwise minimum between two vectors A and B
   Vector Min( const Vector &A, const Vector &B);

   // returns elementwise maximum between two vectors A and B
   Vector Max( const Vector &A, const Vector &B);

	Vector
	Randomize( const Vector &V, const double &MAX);

	Vector
	RandomVector( int SIZE);

	Vector
	RandomVector( const double &LENGTH = 1.0, int SIZE = 3);

	template< typename T>
	double Norm( const std::vector< T> &V)
	{
		double sum = 0.0;
		for( typename std::vector<T>::const_iterator itr = V.begin(); itr != V.end(); ++itr)
		{
			sum += (double) *itr * *itr;
		}
		return sqrt(sum);
	}


} // end namespace math


#endif
