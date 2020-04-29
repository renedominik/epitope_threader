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


#include "../include/vector.h"
#include "../include/math_functions.h"
#include <limits>

namespace math
{

  double& Vector::operator()( const size_t &ID)
  {
    assert( ID < this->size());
    return this->operator[]( ID);
  }


  const double& Vector::operator()( const size_t &ID) const
  {
    assert( ID < this->size());
    return this->operator[]( ID);
  }

  Vector &Vector::operator += ( const Vector &V)
  {
      assert( this->size() == V.size());
      std::vector< double>::const_iterator v_itr( V.begin());
      for( std::vector< double>::iterator itr = this->begin(); itr != this->end(); ++itr, ++v_itr)
      {
          *itr += *v_itr;
      }
      return *this;
  }


  Vector &Vector::operator -= ( const Vector &V)
  {
    assert( this->size() == V.size());
    std::vector< double>::const_iterator v_itr( V.begin());
    for( std::vector< double>::iterator itr = this->begin(); itr != this->end(); ++itr, ++v_itr)
    {
    	*itr -= *v_itr;
    }
    return *this;
  }


  Vector &Vector::operator *= ( const double &VALUE)
  {
      for( std::vector< double>::iterator itr = this->begin(); itr != this->end(); ++itr)
      {
          *itr *= VALUE;
      }
      return *this;
  }

//  // multiply each element of this vector with VALUE
//  Vector Vector::operator *= ( const double &VALUE) const
//  {
//      Vector vec( *this);
//      for( std::vector< double>::iterator itr = vec.begin(); itr != vec.end(); ++itr)
//      {
//          *itr *= VALUE;
//      }
//      return vec;
//  }

  Vector &Vector::operator /= ( const double &VALUE)
  {
      double inverse( 1.0 / VALUE);
      for( std::vector< double>::iterator itr = this->begin(); itr != this->end(); ++itr)
      {
          *itr *= inverse;
      }
      return *this;
  }

  Vector &Vector::operator /= ( const Vector &VEC)
  {
      assert( this->size() == VEC.size());
      std::vector< double>::const_iterator vec_itr = VEC.begin();
      for( std::vector< double>::iterator this_itr = this->begin(); this_itr != this->end(); ++this_itr, ++vec_itr)
      {
          *this_itr /= *vec_itr;
      }
      return *this;
  }

  Vector &Vector::operator -= ( const double &VALUE)
  {
      for( std::vector< double>::iterator itr = this->begin(); itr != this->end(); ++itr)
      {
          *itr -= VALUE;
      }
      return *this;
  }

  Vector &Vector::operator += ( const double &VALUE)
  {
       for( std::vector< double>::iterator itr = this->begin(); itr != this->end(); ++itr)
       {
           *itr += VALUE;
       }
       return *this;
  }

  Vector &Vector::operator = ( const double &VALUE)
  {
    for( std::vector< double>::iterator itr = this->begin(); itr != this->end(); ++itr)
      {
    *itr = VALUE;
      }
    return *this;
  }

  bool Vector::operator < ( const Vector &V) const
  {
      assert( this->size() == V.size());
      for( std::vector< double>::const_iterator this_itr( this->begin()), v_itr( V.begin()); this_itr != this->end(); ++this_itr, ++v_itr)
          if( *this_itr >= *v_itr)
          { return false;}
      return true;
  }

  bool Vector::operator <= ( const Vector &V) const
  {
      assert( this->size() == V.size());
      for( std::vector< double>::const_iterator this_itr( this->begin()), v_itr( V.begin()); this_itr != this->end(); ++this_itr, ++v_itr)
          if( *this_itr > *v_itr)
          { return false;}
      return true;
  }

  bool Vector::operator > ( const Vector &V) const
  {
      assert( this->size() == V.size());
      for( std::vector< double>::const_iterator this_itr( this->begin()), v_itr( V.begin()); this_itr != this->end(); ++this_itr, ++v_itr)
          if( *this_itr <= *v_itr)
          { return false;}
      return true;
  }

  bool Vector::operator >= ( const Vector &V) const
  {
      assert( this->size() == V.size());
      for( std::vector< double>::const_iterator this_itr( this->begin()), v_itr( V.begin()); this_itr != this->end(); ++this_itr, ++v_itr)
          if( *this_itr < *v_itr)
          { return false;}
      return true;
  }


  double Vector::SumOfElements() const
  {
    double sum( 0.0);
    for( std::vector< double>::const_iterator itr( this->begin()); itr != this->end(); ++itr)
      { sum += *itr;}
    return sum;
  }

  double Vector::SubSumOfElements( const size_t &BEGIN, const size_t &NR_ELEMENTS) const
  {
    size_t count( 0);
    double sum( 0.0);
    for( std::vector< double>::const_iterator itr( this->begin() + BEGIN); itr != this->end() && count < NR_ELEMENTS; ++itr, ++count)
      { sum += *itr;}
    return sum;
  }


  double Vector::ProductOfElements() const
  {
    double product( 1.0);
    for( std::vector< double>::const_iterator itr( this->begin()); itr != this->end(); ++itr)
      { product *= *itr;}
    return product;
  }

  double Vector::SubProductOfElements( const size_t &BEGIN, const size_t &NR_ELEMENTS) const
  {
    size_t count( 0);
    double product( 1.0);
    for( std::vector< double>::const_iterator itr( this->begin() + BEGIN); itr != this->end() && count < NR_ELEMENTS; ++itr, ++count)
      { product *= *itr;}
    return product;
  }

  Vector Vector::SubVector( const size_t &BEGIN, const size_t &NR_ELEMENTS) const
  {
    Vector sub( NR_ELEMENTS);
    std::vector< double>::iterator new_itr( sub.begin());
    for( std::vector< double>::const_iterator this_itr( this->begin() + BEGIN); this_itr != this->end() && new_itr != sub.end(); ++this_itr, ++new_itr)
      { *new_itr = *this_itr;}
    return sub;
  }

  double Vector::MeanOfElements() const
  { return SumOfElements() / this->size();}


  std::pair< double, double> Vector::RmsdOfElements() const
  { // could be optimized
    double rmsd( 0.0);
    double mean( MeanOfElements());
    for( std::vector< double>::const_iterator itr = this->begin(); itr != this->end(); ++itr)
      { rmsd += math::Square( ( *itr) - mean);}
    return std::pair< double, double>( mean, sqrt( rmsd / double( this->size() - 1)));
  }

  double Vector::MaxElement() const
  {
    double max( *this->begin());
    for( std::vector< double>::const_iterator itr = this->begin() + 1; itr != this->end(); ++itr)
      if( *itr > max)
    { max = *itr;}
    return max;
  }

  double Vector::MinElement() const
  {
    double min( *this->begin());
    for( std::vector< double>::const_iterator itr = this->begin() + 1; itr != this->end(); ++itr)
      if( *itr < min)
    { min = *itr;}
    return min;
  }

  double Vector::DeltaOfElements() const
  { return MaxElement() - MinElement();}



  std::ostream& Vector::Write( std::ostream &STREAM) const
  {
//	  STREAM << mystr::GetClassName( std::string( __PRETTY_FUNCTION__)) << std::endl;
	  STREAM << this->size() << std::endl;
	  for( std::vector< double>::const_iterator itr = this->begin(); itr != this->end(); ++itr)
      {
		  STREAM.width( 16);
		  STREAM << *itr << "   ";
      }
	  STREAM << std::endl;
	  return STREAM;
  }


  std::istream& Vector::Read( std::istream &STREAM)
  {
//    std::string str;
//    STREAM >> str;
//    if( str != mystr::GetClassName( std::string( __PRETTY_FUNCTION__)))
//    {
//    	std::cout << "===> not the correct id <" << str << "> in " << __PRETTY_FUNCTION__ << std::endl;
//    	STREAM >> str;
//    	std::cout << "try next: " << str << std::endl;
//    	exit( -1);
//    }
    size_t size;
    STREAM >> size;
    *this = Vector( size);
    for( std::vector< double>::iterator itr = this->begin(); itr != this->end(); ++itr)
    {
    	STREAM >> *itr;
    }
    return STREAM;
  }


  bool Vector::IsLengthLargerZero() const
  {
	  for( std::vector< double>::const_iterator itr( this->begin()); itr != this->end(); ++itr)
		  if( *itr != 0)
		  { return true;}
	  return false;
  }


  double Vector::Length() const
  {
	  double sum = 0.0;
	  for( std::vector< double>::const_iterator itr = this->begin(); itr != this->end(); ++itr)
	  {
		  sum += *itr * *itr;
	  }
	  return sqrt( sum);
  }

  double Vector::SquaredLength() const
  {
	  double sum = 0.0;
	  for( std::vector< double>::const_iterator itr = this->begin(); itr != this->end(); ++itr)
	  {
		  sum += *itr * *itr;
	  }
	  return sum;
  }


  Vector &Vector::Normalize()
  {
    *this *= 1.0 / Length();
    return *this;
  }

  Vector Vector::NormalizedCopy() const
  {
      Vector copy( *this);
      copy *= 1.0 / Length();
      return copy;
  }

  double Vector::GetLengthAndNormalize()
  {
	  double length = Length();
      *this *= 1.0 / length;
      return length;
  }

  void Vector::SetAll( const double& X)
  {
      for( std::vector< double>::iterator itr = this->begin(); itr != this->end(); ++itr)
      { *itr = X;}
  }

  Vector &Vector::Randomize( const double &MIN, const double &MAX)
  {
      for( std::vector< double>::iterator itr = this->begin(); itr != this->end(); ++itr)
      {  *itr = math::Random( MIN, MAX);}
      return *this;
  }

  // randomize all elements within limits
  Vector &Vector::Randomize( const Vector &MIN, const Vector &MAX)
  {
      std::vector< double>::const_iterator min_itr = MIN.begin(), max_itr = MAX.begin();
      for( std::vector< double>::iterator itr = this->begin(); itr != this->end(); ++itr, ++max_itr, ++min_itr)
      {
          *itr = math::Random( *min_itr, *max_itr);
      }
      return *this;
  }


  Vector &Vector::SetToLength( const double &MAX)
  {
      double length( Length());
//      DebugWrite( "Vector::SetToLength(): length_before: " << length);
      if( length > 0 && MAX >= 0.0 && MAX < std::numeric_limits< double>::max())
      {
    	  *this *= ( MAX / length);
#ifdef DEBUG
    	  std::cout << "Vector::SetToLength(): length_after: " << Length() << std::endl;
      }
      else
      {
    	  std::cout << "==> nothing done in " << __FUNCTION__ << ": useless values: length: " << length << " value: " << MAX << std::endl;
#endif
      }
      return *this;
  }


  Vector &Vector::SetNonZeroElementsToMinusOrPlusOne()
  {
      for( std::vector< double>::iterator itr = this->begin(); itr != this->end(); ++itr)
      {
          if( *itr > 0)
          {
              *itr = 1;
          }
          else if( *itr < 0)
          {
              *itr = -1;
          }
      }
      return *this;
  }

  Vector &Vector::Merge( const Vector &V)
  {
      this->insert( this->end(), V.begin(), V.end());
      return *this;
  }

  bool Vector::AreAllElements( const double &VALUE) const
  {
      for( std::vector< double>::const_iterator itr = this->begin(); itr != this->end(); ++itr)
    	  if( *itr != VALUE)
		  {
    		  return false;
		  }
      return true;
  }


  // transforms vector into new coordinate system
   Vector Vector::CooTrafo( const std::vector< Vector> &COO)
   {
	   size_t size = this->size();

	   Vector v( size);

	   double a[ size*size], b[size];
	   for( size_t i = 0; i < COO[0].size(); ++i)
		   for( size_t j = 0; j < COO.size(); ++j)
		   {
			   a[i] = COO[j][i];
		   }
	   for( size_t i = 0; i < size; ++i)
	   {
		    b[i] = this->operator[](i);
	   }

	   exit(1);
//	  gsl_matrix_view m
//		 = gsl_matrix_view_array (a, size, size);
//
//	   gsl_vector_view c
//		 = gsl_vector_view_array (b, size);
//
//	   gsl_vector *x = gsl_vector_alloc (size);
//
//	   int s;
//
//	   gsl_permutation * p = gsl_permutation_alloc (size);
//
//	   gsl_linalg_LU_decomp (&m.matrix, p, &s);
//
//	   gsl_linalg_LU_solve (&m.matrix, p, &c.vector, x);
//
//	   for( size_t i = 0 ; i < size ; ++i )
//	     v[i] = gsl_vector_get( x , i);
//
//	   gsl_permutation_free (p);
//	   gsl_vector_free (x);

	   return v;
   }




//  Vector &Vector::RemoveDuplicates()
//  {
//      std::vector< double> tmp;
//      std::vector< double>::iterator new_end = std::unique( this->begin(), this->end());
//      tmp.insert( tmp.end(), this->begin(), new_end);
//      *this = tmp;
//      return *this;
//  }

} // end namespace math
