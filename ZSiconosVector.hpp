/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*! \file SiconosVector.hpp
 */

#ifndef __ZSiconosVector__
#define __ZSiconosVector__

#include "SiconosAlgebraTypeDef.hpp"
#include "SiconosVectorFriends.hpp"
#include "SiconosSerialization.hpp" // For ACCEPT_SERIALIZATION

#include <memory_resource>

DEFINE_SPTR(ZSiconosVector)

struct SiconosVectorIterator;
struct SiconosVectorConstIterator;

using ZDenseVect = ublas::vector<double, std::vector<double, std::pmr::polymorphic_allocator<double>  >>;

/** Union to gather all types of ublas vectors used in Siconos */
union ZVECTOR_UBLAS_TYPE
{
  ZDenseVect *Dense; // num = 1
  SparseVect *Sparse; // num = 4
};


/**
   Vectors of double. (Interface to various types of Boost-Ublas vectors).
   
   Two possible types: Siconos::DENSE (default) and Siconos:SPARSE.
   
*/
class ZSiconosVector : public std::enable_shared_from_this<ZSiconosVector>
{
protected:
  ACCEPT_SERIALIZATION(ZSiconosVector);

  bool _dense = true;

  /**
   * Union of pointers to the ublas vector type (dense or sparse)
   */
  ZVECTOR_UBLAS_TYPE vect;

public:

  /***************************** CONSTRUCTORS ****************************/

  /** Creates a zero-size vector. */
  ZSiconosVector();

  /** creates a vector, all components set to zero.
   *
   *  \param row the size of the vector
   *  \param type the type of vector (dense or sparse)
   */
  ZSiconosVector(unsigned row, Siconos::UBLAS_TYPE type = Siconos::DENSE);

  /** creates a vector and initializes its content with a single value
   *
   *  \param row size of the new vector
   *  \param val value to initialize its content
   *  \param type type of vector (dense or sparse)
   */
  ZSiconosVector(unsigned row, double val, Siconos::UBLAS_TYPE type = Siconos::DENSE);

  /** creates a dense vector from a copy of a stl vector.
   *
   *  \param vec vector to be copied
   *  \param type of the vector (dense or sparse)
   */
  ZSiconosVector(const std::vector<double>& vec, Siconos::UBLAS_TYPE type = Siconos::DENSE);

  /** copy constructor
   *
   *  \param v source vector to be copied
   */
  ZSiconosVector(const ZSiconosVector& v);

  /** creates a dense vector, with a copy.
   *
   *  \param v source vector (ublas dense)
   */
  ZSiconosVector(const ZDenseVect& v);

  /** creates a sparse vector, with a copy.
   *
   *  \param v source vector (ublas sparse)
   */
  ZSiconosVector(const SparseVect& v);

  /** creates a vector from data in a file
   *
   *  \param filename file name (possibly with path)
   *  \param is_ascii file format (true if ascii, false if binary)
   */
  ZSiconosVector(const std::string& filename, bool is_ascii);

  /** constructor from the concatenation of two vectors
   *
   *  \param v1 the first vector
   *  \param v2 the second vector
   */
  ZSiconosVector(const ZSiconosVector& v1, const ZSiconosVector& v2);

  /** constructor from a BlockVector.
   *  explicit to forbid implicit conversion/conversion constructor.
   *
   *  \param input source vector
   */
  explicit ZSiconosVector(const BlockVector& input);//, bool = false);

  /** destructor
   */
  ~ZSiconosVector();

  /** get the vector size, ie the total number of (double) elements in the vector
   *
   *  \return unsigned int
   */
  unsigned int size() const;

  /** Get the type number of the current vector.
   *
   *  \return an unsigned int
   */
  Siconos::UBLAS_TYPE num() const
  {
    if (_dense) return Siconos::DENSE;
    else return Siconos::SPARSE;
  }

  /** get a pointer to the ublas embedded vector if it's type is Dense
   *
   *  \return a ZDenseVect*
   */
  inline ZDenseVect* dense() const
  {
    return vect.Dense;
  };

  /** get a pointer to the ublas embedded vector if it's type is Sparse
   *
   *  \return a SparseVect*
   */
  SparseVect* sparse() const;

  /** \return the array of double values of the vector
   */
  double* getArray() const;

  /** sets all the values of the vector to 0.0 */
  void zero();

  /** Resize the vector. The existing elements may be preseved if specified.
   *
   *  \param size new size of the vector
   *  \param preserve true if the content of the vector must be preserved.
   */
  void resize(unsigned int size, bool preserve= true);

  /** \return the infinite norm of the vector */
  double normInf()const;

  /** \return the Euclidian norm of the vector */
  double norm2() const ;

  /** \return the sum of all elements of the vector */
  double vector_sum() const;

  /** display vector content */
  void display(void) const;

  /** set all values of the vector to input value.
   *
   *  \param a input value
   */
  void fill(double a);

  /** \return the content of the vector as a string */
  std::string toString() const;

  /** for iterator interface */
  typedef SiconosVectorIterator iterator;

  /** for iterator interface */
  typedef SiconosVectorConstIterator const_iterator;

  /** \return an iterator pointing to the first element in the vector. */
  iterator begin();

  /** \return an iterator pointing to the first element in the vector. */
  const_iterator begin() const;

  /**  \return an iterator referring to the past-the-end element in the vector container. */
  iterator end();

  /**  \return an iterator referring to the past-the-end element in the vector container. */
  const_iterator end() const;

  /** cast a ZSiconosVector into a std::vector<double> (performs copy) */
  operator std::vector<double>();

  //************************** VECTORS HANDLING AND OPERATORS *******************************

  /** Get a component of the vector
   *
   *  \param i index of the required component
   *  \return the component value
   */
  double getValue(unsigned int i) const ;

  /** set a component of the vector
   *
   *  \param i index of the required component
   *  \param value of the component
   */
  void setValue(unsigned int i, double value);

  /** get a component of the vector
   *
   *  \param i index of the required component
   *  \return value of the component
   */
  double& operator()(unsigned int i);

  /** get a component of the vector
   *
   *  \param i index of the required component
   *  \return value of the component
   */
  double operator()(unsigned int i) const;

  /** set a sub-block of the current vector
      
      \param i the beginning of the destination range
      \param v vector to be copied
  */
  void setBlock(unsigned int i, const ZSiconosVector& v);

  /** 
      copy a part of the vector into another 
    
      \param vOut destination vector
      \param sizeB number of the elements to copy
      \param startIn the beginning of the range of elements to copy from
      \param startOut the beginning of the destination range
  */
  void toBlock(ZSiconosVector& vOut, unsigned int sizeB,
               unsigned int startIn, unsigned int startOut) const;

  /** 
      add the input vector to a sub-block of the current vector
      
      \param i the beginning of the destination range
      \param v the source vector to be added
   */
  void addBlock(unsigned int i, const ZSiconosVector& v);

  /** 
      subtract the input vector to a sub-block of the current vector
      
      \param i the beginning of the destination range
      \param v the source vector to be added
   */
  void subBlock(unsigned int i, const ZSiconosVector& v);

  /** copy the vector into an array
   *
   *  \param data the memory where to copy the data
   *  \return the number of element written (size of the vector)
   */
  unsigned copyData(double* data) const;

  /** operator =
   *
   *  \param v the vector to be copied
   *  \return  ZSiconosVector&
   */
  ZSiconosVector& operator = (const ZSiconosVector& v);

  /** operator =
   *
   *  \param b the vector to be copied
   *  \return  ZSiconosVector&
   */
  ZSiconosVector& operator = (const BlockVector& b);

  /** operator =
   *
   *  \param v the vector to be copied
   *  \return  ZSiconosVector&
   */
  ZSiconosVector& operator = (const ZDenseVect& v);

  /** operator =
   *
   *  \param sp the vector to be copied
   *  \return  ZSiconosVector&
   */
  ZSiconosVector& operator = (const SparseVect& sp);

  /** operator =
   *
   *  \param d data to put the in vector
   *  \return  ZSiconosVector&
   */
  ZSiconosVector& operator = (const double* d);

  /** operator +=
   *
   *  \param v the vector to add
   *  \return  ZSiconosVector&
   */
  ZSiconosVector& operator +=(const ZSiconosVector& v);

  /** operator +=
   *
   *  \param v the vector to add
   *  \return  ZSiconosVector&
   */
  ZSiconosVector& operator +=(const BlockVector& v);

  /** operator -=
   *
   *  \param  v the vector to subtract
   *  \return  ZSiconosVector&
   */
  ZSiconosVector& operator -=(const ZSiconosVector& v);
  /** operator -=
   *
   *  \param  v the vector to subtract
   *  \return  ZSiconosVector&
   */
  ZSiconosVector& operator -=(const BlockVector& v);

  /** \defgroup ZSiconosVectorFriends

      List of friend functions of the ZSiconosVector class

      @{
  */

  /** send data of the vector to an ostream
   *
   *  \param os An output stream
   *  \param sv a ZSiconosVector
   *  \return The same output stream
   */
  friend std::ostream& operator<<(std::ostream& os, const ZSiconosVector& sv);

  friend ZSiconosVector& operator *= (ZSiconosVector& v, const double& s);

  friend ZSiconosVector& operator /= (ZSiconosVector& v, const double& s);

  friend bool operator ==(const ZSiconosVector&, const ZSiconosVector&);

  friend ZSiconosVector operator * (double, const ZSiconosVector&);

  friend ZSiconosVector operator * (const ZSiconosVector&, double);

  friend ZSiconosVector operator / (const ZSiconosVector&, double);

  friend ZSiconosVector operator + (const ZSiconosVector&, const ZSiconosVector&);

  friend void add(const ZSiconosVector&, const ZSiconosVector&, ZSiconosVector&);

  friend ZSiconosVector operator - (const ZSiconosVector&, const ZSiconosVector&);

  friend void sub(const ZSiconosVector&, const ZSiconosVector&, ZSiconosVector&);

  friend void axpby(double, const ZSiconosVector&, double, ZSiconosVector&);

  friend void axpy(double, const ZSiconosVector&, ZSiconosVector&);

  friend double inner_prod(const ZSiconosVector&, const ZSiconosVector&);

  friend SimpleMatrix outer_prod(const ZSiconosVector&, const ZSiconosVector&);

  friend void scal(double, const ZSiconosVector&, ZSiconosVector&, bool);

  friend void subscal(double, const ZSiconosVector&, ZSiconosVector&, const Index&, bool);

  friend void cross_product(const ZSiconosVector&, const ZSiconosVector&, ZSiconosVector&);

  friend void abs_wise(const ZSiconosVector&, ZSiconosVector&);

  friend void getMax(const ZSiconosVector&, double &, unsigned int &);

  friend void  getMin(const ZSiconosVector&, double &, unsigned int &);

  friend struct IsDense;

  friend struct IsSparse;

  friend struct IsBlock;

  friend class TestDense;

  /** End of Friend functions group @} */

  //  temporary workaround, the visitor has to be removed or rework -- xhub
//  ACCEPT_NONVIRTUAL_VISITORS();

};


#include "Question.hpp"

struct IsDense : public Question<bool>
{
  using SiconosVisitor::visit;

  void visit(const ZSiconosVector& v)
  {
    answer = v._dense;
  }

  void visit(const BlockVector& v)
  {
    answer = false;
  }
};

struct IsSparse : public Question<bool>
{

  using SiconosVisitor::visit;

  void visit(const ZSiconosVector& v)
  {
    answer = !v._dense;
  }

  void visit(const BlockVector& v)
  {
    answer = false;
  }
};

struct IsBlock : public Question<bool>
{
  using SiconosVisitor::visit;

  void visit(const ZSiconosVector& v)
  {
    answer = false;
  }

  void visit(const BlockVector& v)
  {
    answer = true;
  }
};


#endif
