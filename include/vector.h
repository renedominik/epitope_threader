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


#ifndef MATH_VECTOR_H
#define MATH_VECTOR_H

#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <algorithm>


namespace math
{

  ///////////////////////////////////////////////
  // Wrapper class around std::vector< double>
  // with secure element access.
  //
  // author Rene Staritzbichler
  // date 19.11.2008
  // @example "../../example/math/vector.cpp", "../../example/math/vector3N.cpp","../../example/math/vector_functions.cpp"
  ///////////////////////////////////////////////

  class Vector
    : public std::vector< double>
  {
  public:
    // default constructor
    Vector()
    : std::vector<double>()
    {}

    // construction from size
    Vector( const size_t &SIZE)
    : std::vector< double>( SIZE)
    { /*SetAll( 0.0);*/}

    // construction from size and value
    Vector( const size_t &SIZE, const double &VALUE)
    : std::vector< double>( SIZE, VALUE)
    { /*SetAll( 0.0);*/}

      Vector( const std::vector< size_t> &V)
        : std::vector< double>( V.size())
        {
          std::vector< size_t>::const_iterator v_itr( V.begin());
          for( std::vector< double>::iterator this_itr( this->begin()); this_itr != this->end() && v_itr != V.end(); ++v_itr, ++this_itr)
        { *this_itr = double( *v_itr);}
        }

      Vector( const std::vector< double> &V)
        : std::vector< double>( V)
        {}


      // copy constructor
        Vector( const Vector &V)
          : std::vector< double>( V)
          {}


      // virtual destructor
      virtual ~Vector(){}

      // mutable element access
      virtual double& operator()( const size_t &ID);

      // unmutable element access
      virtual const double& operator()( const size_t &ID) const;

      // add another vector
      virtual Vector& operator += ( const Vector &V);

      // substract another vector
      virtual Vector& operator -= ( const Vector &V);

//      // set this vector to the values of another vector
//      virtual Vector operator = ( const Vector &V);

      // multiply each element of this vector with VALUE
      virtual Vector& operator *= ( const double &VALUE);
      // multiply each element of this vector with VALUE
 //     virtual Vector operator *= ( const double &VALUE) const;
      // devide each element of this vector by VALUE
      virtual Vector& operator /= ( const double &VALUE);
      // substract VALUE from each element of this vector
      virtual Vector& operator -= ( const double &VALUE);
      // add VALUE to each element of this vector
      virtual Vector& operator += ( const double &VALUE);
      // set each element of this vector to VALUE
      virtual Vector &operator = ( const double &VALUE);
//      // check whether this vector is equal to V
//      virtual bool operator == ( const Vector &V);

      // checks whether condition '<' is fullilled for all elements
      virtual bool operator < ( const Vector &V) const;

      // checks whether condition '<=' is fullilled for all elements
      virtual bool operator <= ( const Vector &V) const;

      // checks whether condition ">' is fullilled for all elements
      virtual bool operator > ( const Vector &V) const;

      // checks whether condition '>=' is fullilled for all elements
      virtual bool operator >= ( const Vector &V) const;

      virtual double SumOfElements() const;

      virtual double SubSumOfElements( const size_t &BEGIN, const size_t &NR_ELEMENTS) const;

      virtual double ProductOfElements() const;

      virtual double SubProductOfElements( const size_t &BEGIN, const size_t &NR_ELEMENTS) const;

      virtual Vector SubVector( const size_t &BEGIN, const size_t &NR_ELEMENTS) const;

      virtual double MeanOfElements() const;

      virtual std::pair< double, double> RmsdOfElements() const;

      virtual double MaxElement() const;

      virtual double MinElement() const;

      virtual double DeltaOfElements() const;

      virtual void SetAll( const double &VALUE);

      // check whether length of vector is larger zero (faster than explicit calculation of the length)
      virtual bool IsLengthLargerZero() const;

      // returns the length of the vector
      virtual double Length() const;

      // returns the squared length of the vector
      virtual double SquaredLength() const;

      // normalize the vector
      virtual Vector &Normalize();

      // return normalized copy of vector
      virtual Vector NormalizedCopy() const;

      // returns length and normalizes vector
      virtual double GetLengthAndNormalize();

      // write vector to ostream
      virtual std::ostream& Write( std::ostream &STREAM) const;

      // read vector from istream
      virtual std::istream& Read( std::istream &STREAM);

      // randomize all elements within limits
      virtual Vector &Randomize( const double &MIN, const double &MAX);

      // randomize all elements within limits
      virtual Vector &Randomize( const Vector &MIN, const Vector &MAX);

      // set vector length to MAX if larger MAX
      virtual Vector &SetToLength( const double &MAX);

      // sets all negative elements to -1, all positive to +1, 0 remains 0
      virtual Vector &SetNonZeroElementsToMinusOrPlusOne();

      // devide each element of this vector by element of VECTOR
      virtual Vector &operator /= ( const Vector &VECTOR);

      // insert all elements of Vector in this
      virtual Vector &Merge( const Vector &V);

      // checks whether all element in the vector are equal VALUE
      virtual bool AreAllElements( const double &VALUE) const;


      // transforms vector into new coordinate system
      virtual Vector CooTrafo( const std::vector< Vector> &COO);


//      // remove multiple copies
//      virtual Vector &RemoveDuplicates();



  }; // end class Vector

} // end namespace math

#endif
