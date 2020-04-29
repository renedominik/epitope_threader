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


#ifndef MATH_VECTOR3N_H
#define MATH_VECTOR3N_H

#include "vector.h"
#include "vector_functions.h"

namespace math
{

  class Vector3N
    : public Vector
  {
  public:
    Vector3N();

    Vector3N( const Vector3N &V);

    Vector3N( const Vector &V);
    Vector3N( const std::vector< double> &V);

    Vector3N( const double &X);

    Vector3N( const double &X, const double &Y, const double &Z);

    virtual ~Vector3N();

    virtual Vector3N *Clone() const;

    virtual const double& operator() ( const size_t &ID) const;

    virtual double& operator() ( const size_t &ID);

    Vector3N &operator = ( const std::vector< double> &VEC);

    virtual std::ostream &Write( std::ostream &STREAM) const;

    virtual std::istream &Read( std::istream &STREAM);

    // insert all elements of Vector in this
    virtual Vector3N &Merge( const Vector3N &V)
    {
        std::cout << "calling this function of a vector of fixed size does not make any sense!" << std::endl;
        exit( -1);
        return *this;
    }

    // remove multiple copies
    virtual Vector3N &Unique()
    {
        std::cout << "calling this function of a vector of fixed size does not make any sense!" << std::endl;
        exit( -1);
        return *this;
    }



  }; // class Vector3N

    inline
    Vector3N CrossProduct( const Vector3N &X, const Vector3N &Y)
    {
        return Vector3N
        (
                X( 1) * Y( 2) - X( 2) * Y( 1),
             X( 2) * Y( 0) - X( 0) * Y( 2),
             X( 0) * Y( 1) - X( 1) * Y( 0)
        );
    }

    // Spatprodukt - the volume created from the three vectors
    inline
    double TripleProduct( const Vector3N &X, const Vector3N &Y, const Vector3N &Z)
    {
        return ScalarProduct( CrossProduct( X, Y), Z);
    }

    // dihedral of three vectors
    double Dihedral( const Vector &V1, const Vector &V2, const Vector &V3);

    // dihedral of three points
    double Dihedral( const Vector &P1, const Vector &P2, const Vector &P3, const Vector &P4);



} // namespace math

#endif
