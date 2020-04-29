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
#include "../include/vector3N.h"


namespace math
{

  Vector3N::Vector3N()
    : Vector(3, 0.0)
  {}

  Vector3N::Vector3N( const Vector3N &V)
    : Vector( V)
  { assert( V.size() == 3);}

  Vector3N::Vector3N( const Vector &V)
    : Vector( V)
  { assert( V.size() == 3);}

  Vector3N::Vector3N( const std::vector< double> &V)
    : Vector( V)
  { assert( V.size() == 3);}

  Vector3N::Vector3N( const double &X)
    : Vector( 3)
  {
    ( *this)[0] = X;
    ( *this)[1] = X;
    ( *this)[2] = X;
  }

  Vector3N::Vector3N( const double &X, const double &Y, const double &Z)
    : Vector( 3)
  {
    ( *this)[0] = X;
    ( *this)[1] = Y;
    ( *this)[2] = Z;
  }

  Vector3N::~Vector3N(){}

  Vector3N *Vector3N::Clone() const
  { return new Vector3N( *this);}

  double &Vector3N::operator()( const size_t &ID)
  {
//    assert( ID < 3);
    return ( *this)[ ID];
  }

  const double& Vector3N::operator() ( const size_t &ID) const
  {
//    assert( ID < 3);
    return ( *this)[ ID];
  }

  Vector3N &Vector3N::operator = ( const std::vector< double> &VEC)
  {
      assert( VEC.size() == 3);
        ( *this)[0] = VEC[0];
        ( *this)[1] = VEC[1];
        ( *this)[2] = VEC[2];
        return *this;
  }


  std::istream &Vector3N::Read( std::istream &STREAM)
  { return Vector::Read( STREAM);}


  std::ostream &Vector3N::Write( std::ostream &STREAM) const
  { return Vector::Write( STREAM);}


  double Dihedral( const Vector &V1, const Vector &V2, const Vector &V3)
  {
	  math::Vector v = CrossProduct( V2, V3);
	  return atan2( V2.Length() * ScalarProduct( V1, v), ScalarProduct( CrossProduct( V1, V2), v));
  }

  double Dihedral( const Vector &P1, const Vector &P2, const Vector &P3, const Vector &P4)
  {
//	  std::cout << __FUNCTION__ << P1 << " " << P2 << " " << P3 << " " << P4 << std::endl;
	  return Dihedral( P2-P1, P3-P2, P4-P3);
  }



} // end namespace math
