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

#include "SiconosConfig.h"

#include "SiconosAlgebraTypeDef.hpp"
#include "SiconosVectorIterator.hpp"
#include "Tools.hpp"

#include <boost/numeric/ublas/io.hpp>            // for >> 
//#include <boost/numeric/ublas/vector_proxy.hpp>  // for project
#include <boost/numeric/ublas/vector_sparse.hpp>


#include "boost/numeric/bindings/ublas/vector_proxy.hpp"
#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/std/vector.hpp>
namespace siconosBindings = boost::numeric::bindings::blas;

#include "SimpleMatrix.hpp"
#include "BlockVector.hpp"
#include "ioVector.hpp"
#include "XSiconosVector.hpp"
#include "SiconosAlgebra.hpp"
#include <cmath>        // std::exp(double)
#include <algorithm>    // std::transform
#include "SiconosException.hpp"
//#define DEBUG_MESSAGES
#include "siconos_debug.h"

// Do not document
/// @cond


/// @endcond


// =================================================
//                CONSTRUCTORS
// =================================================

// Default
XSiconosVector::XSiconosVector()
{
  _dense = true;
  vect.Dense = new XDenseVect(ublas::zero_vector<double>());
}

// parameters: dimension and type.
XSiconosVector::XSiconosVector(unsigned row, Siconos::UBLAS_TYPE type)
{
  if(type == Siconos::SPARSE)
  {
    _dense = false;
    vect.Sparse = new SparseVect(ublas::zero_vector<double>(row));
  }
  else if(type == Siconos::DENSE)
  {
    _dense = true;
    vect.Dense = new XDenseVect(ublas::zero_vector<double>(row));
  }
  else
  {
    THROW_EXCEPTION("invalid type");
  }
}

// parameters: dimension, default value for all components and type.
XSiconosVector::XSiconosVector(unsigned row, double val, Siconos::UBLAS_TYPE type)
{
  if(type == Siconos::SPARSE)
  {
    _dense = false;
    vect.Sparse = new SparseVect(row);
    fill(val);
  }
  else if(type == Siconos::DENSE)
  {
    _dense = true;
    vect.Dense = new XDenseVect(ublas::scalar_vector<double>(row, val));
  }
  else
  {
    THROW_EXCEPTION("invalid type");
  }
}

// parameters: a vector (stl) of double and the type.
XSiconosVector::XSiconosVector(const std::vector<double>& v, Siconos::UBLAS_TYPE typ)
{
  if(typ != Siconos::DENSE)
    THROW_EXCEPTION("invalid type");

  _dense = true;
  vect.Dense = new XDenseVect(v.size());
  std::copy(v.begin(), v.end(), (vect.Dense)->begin());
}

// Copy
XSiconosVector::XSiconosVector(const XSiconosVector &svect) : std::enable_shared_from_this<XSiconosVector>()
{
  if(true)  // dense
  {
    _dense = true;
    vect.Dense = new XDenseVect(svect.size());
    noalias(*vect.Dense) = (*svect.dense());
    // std::copy((vect.Dense)->begin(), (vect.Dense)->end(), (svect.dense())->begin());
  }
  else //sparse
  {
    _dense = false;
    vect.Sparse = new SparseVect(svect.size());
    noalias(*vect.Sparse) = (*svect.sparse());
    //std::copy((vect.Sparse)->begin(), (vect.Sparse)->end(), (svect.sparse())->begin());
  }
  // Note FP: using constructor + noalias = (or std::copy) seems to be more
  // efficient than a call to ublas::vector copy constructor, this for
  // large or small vectors.
}


XSiconosVector::XSiconosVector(const XDenseVect& m)
{
  _dense = true;
  vect.Dense = new XDenseVect(m.size());
  noalias(*vect.Dense) = m;

}

XSiconosVector::XSiconosVector(const SparseVect& m)
{
  _dense = false;
  vect.Sparse = new SparseVect(m.size());
  noalias(*vect.Sparse) = m;
}

XSiconosVector::XSiconosVector(const std::string &file, bool ascii)
{
}

XSiconosVector::XSiconosVector(const XSiconosVector& v1, const XSiconosVector& v2)
{
}

// Copy a block vector into a XSiconosVector
// This is mostly used to handle contiguous memory.
XSiconosVector::XSiconosVector(const BlockVector & input)
{

}


XSiconosVector::~XSiconosVector()
{
  if(_dense)
    delete(vect.Dense);
  else
    delete(vect.Sparse);
}


// =================================================
//        get Ublas component (dense or sparse)
// =================================================

SparseVect* XSiconosVector::sparse()const
{

  if(_dense)
    THROW_EXCEPTION("the current vector is not a Sparse vector");

  return vect.Sparse;
}

double* XSiconosVector::getArray() const
{
  assert(vect.Dense && "XSiconosVector::getArray() : not yet implemented for sparse vector.");

  return &(((*vect.Dense).data())[0]);
}

// ===========================
//       fill vector
// ===========================

void XSiconosVector::zero()
{
  if(_dense)
    siconosBindings::scal(0.0, *vect.Dense);

  else
  {
    assert(vect.Sparse);
    *vect.Sparse *= 0.0;
  }

}

void XSiconosVector::fill(double value)
{
  if(!_dense)
  {
    for(unsigned int i = 0; i < (vect.Sparse)->size(); ++i)
      (vect.Sparse)->push_back(i, value);
  }
  else
    siconosBindings::set(value, *vect.Dense);


}

//=======================
// set vector dimension
//=======================

void XSiconosVector::resize(unsigned int n, bool preserve)
{
  if(_dense)
    (vect.Dense)->resize(n, preserve);
  else
    (vect.Sparse)->resize(n, preserve);
}

//=======================
//       get norm
//=======================

double XSiconosVector::normInf() const
{
  if(_dense)
    return norm_inf(*vect.Dense);
  else //if(num==Siconos::SPARSE)
    return norm_inf(*vect.Sparse);
}

double XSiconosVector::norm2() const
{
  if(_dense)
    return ublas::norm_2(*vect.Dense);
  else //if(num==Siconos::SPARSE)
    return ublas::norm_2(*vect.Sparse);
}
//======================================
// get sum of all elements of the vector
//=====================================
double XSiconosVector::vector_sum() const
{
  if(_dense)
    return ublas::sum(*vect.Dense);
  else
    return ublas::sum(*vect.Sparse);
}

//=====================
// screen display
//=====================

void XSiconosVector::display()const
{
  std::cout.setf(std::ios::scientific);
  std::cout.precision(6);
  if(_dense)
    std::cout << *vect.Dense << std::endl;
  else if(vect.Sparse)
    std::cout << *vect.Sparse << std::endl;
}

//============================
// Convert vector to a std::string
//============================

std::string XSiconosVector::toString() const
{
  return ::toString(*this);
}

//=====================
// convert to an ostream
//=====================

std::ostream& operator<<(std::ostream& os, const XSiconosVector& sv)
{
  if(sv._dense)
    os << *sv.vect.Dense;
  else
    os << *sv.vect.Sparse;
  return os;
}

//=============================
// Elements access (get or set)
//=============================

double XSiconosVector::getValue(unsigned int row) const
{
  assert(row < size() && "XSiconosVector::getValue(index) : Index out of range");

  if(_dense)
    return (*vect.Dense)(row);
  else
    return (*vect.Sparse)(row);
}

void XSiconosVector::setValue(unsigned int row, double value)
{
  assert(row < size() && "XSiconosVector::setValue(index, value) : Index out of range");
  if(_dense)
    (*vect.Dense)(row) = value ;
  else
    (*vect.Sparse)(row) = value;
}

double& XSiconosVector::operator()(unsigned int row)
{
  assert(row < size() && "XSiconosVector::operator ( index ): Index out of range");

  if(_dense)
    return (*vect.Dense)(row);
  else
    return (*vect.Sparse)(row).ref();
}

double XSiconosVector::operator()(unsigned int row) const
{
  assert(row < size() && "XSiconosVector::operator ( index ): Index out of range");

  if(_dense)
    return (*vect.Dense)(row);
  else
    return ((*vect.Sparse)(row)).ref();
}

//============================================
// Access (get or set) to blocks of elements
//============================================

void XSiconosVector::setBlock(unsigned int index, const XSiconosVector& vIn)
{
  // Set current vector elements, starting from position "index", to the values of vector vIn

  // Exceptions ...
  assert(&vIn != this && "XSiconosVector::this->setBlock(pos,vIn): vIn = this.");

  assert(index < size() && "XSiconosVector::setBlock : invalid ranges");

  auto end = vIn.size() + index;
  assert(end <= size() && "XSiconosVector::setBlock : invalid ranges");

  assert(vIn.num() == num() && "XSiconosVector::setBlock: inconsistent types.");

  if(_dense)
    noalias(ublas::subrange(*vect.Dense, index, end)) = *vIn.dense();
  else
    noalias(ublas::subrange(*vect.Sparse, index, end)) = *vIn.sparse();
}

void XSiconosVector::toBlock(XSiconosVector& vOut, unsigned int sizeB, unsigned int startIn, unsigned int startOut) const
{
  // To copy a subBlock of the vector (from position startIn to startIn+sizeB) into vOut (from pos. startOut to startOut+sizeB).
  // Check dim ...
  assert(startIn < size() && "vector toBlock(v1,v2,...): start position in input vector is out of range.");

  assert(startOut < vOut.size() && "vector toBlock(v1,v2,...): start position in output vector is out of range.");

  assert(startIn + sizeB <= size() && "vector toBlock(v1,v2,...): end position in input vector is out of range.");
  assert(startOut + sizeB <= vOut.size() && "vector toBlock(v1,v2,...): end position in output vector is out of range.");

  unsigned int endOut = startOut + sizeB;
  Siconos::UBLAS_TYPE numIn = num();
  Siconos::UBLAS_TYPE numOut = vOut.num();

  if(numIn == numOut)
  {
    if(numIn == Siconos::DENSE)
      noalias(ublas::subrange(*vOut.dense(), startOut, endOut)) = ublas::subrange(*vect.Dense, startIn, startIn + sizeB);
    else // if(numIn == Siconos::SPARSE)// vIn / vOut are Sparse
      noalias(ublas::subrange(*vOut.sparse(), startOut, endOut)) = ublas::subrange(*vect.Sparse, startIn, startIn + sizeB);
  }
  else // vIn and vout of different types ...
  {
    if(numIn == Siconos::DENSE)  // vIn Dense
      noalias(ublas::subrange(*vOut.sparse(), startOut, endOut)) = ublas::subrange(*vect.Dense, startIn, startIn + sizeB);
    else // if(numIn == Siconos::SPARSE)// vIn Sparse
      noalias(ublas::subrange(*vOut.dense(), startOut, endOut)) = ublas::subrange(*vect.Sparse, startIn, startIn + sizeB);
  }
}

void XSiconosVector::addBlock(unsigned int index, const XSiconosVector& vIn)
{
  // Add vIn to the current vector, starting from position "index".

  if(&vIn == this)
    THROW_EXCEPTION("try to add a vector to itself.");

  unsigned int end = vIn.size();
  if((index + end) > size())
    THROW_EXCEPTION("invalid ranges");

  Siconos::UBLAS_TYPE numVin = vIn.num();

  if(numVin != num())
    THROW_EXCEPTION("inconsistent types.");

  if(_dense)
    noalias(ublas::subrange(*vect.Dense, index, index + end)) += *vIn.dense();
  else
    noalias(ublas::subrange(*vect.Sparse, index, index + end)) += *vIn.sparse();
}

void XSiconosVector::subBlock(unsigned int index, const XSiconosVector& vIn)
{
  // Add vIn from the current vector, starting from position "index".

  unsigned int end = vIn.size();
  if((index + end) > size())
    THROW_EXCEPTION("invalid ranges");

  Siconos::UBLAS_TYPE numVin = vIn.num();
  if(numVin != num())
    THROW_EXCEPTION("inconsistent types.");

  if(_dense)
    noalias(ublas::subrange(*vect.Dense, index, index + end)) -= *vIn.dense();
  else
    noalias(ublas::subrange(*vect.Sparse, index, index + end)) -= *vIn.sparse();
}

//===============
//  Assignment
//===============

XSiconosVector& XSiconosVector::operator = (const XSiconosVector& vIn)
{
  if(&vIn == this) return *this;  // auto-assignment.

  assert(size() == vIn.size() && "XSiconosVector::operator = failed: inconsistent sizes.");

  Siconos::UBLAS_TYPE vInNum = vIn.num();
  {
    switch(num())
    {
    case Siconos::DENSE:
      switch(vInNum)
      {
      case Siconos::DENSE:
        //siconosBindings::copy(*vIn.dense(),*vect.Dense);
        noalias(*vect.Dense) = *vIn.dense();
        break;
      case Siconos::SPARSE:
        noalias(*vect.Dense) = *vIn.sparse();
        break;
      default:
        THROW_EXCEPTION("invalid type");
        break;
      }
      break;
    case Siconos::SPARSE:
      if (vInNum == Siconos::DENSE)
        noalias(*vect.Sparse) = *vIn.dense();
      else if(vInNum == Siconos::SPARSE)
        noalias(*vect.Sparse) = *vIn.sparse();
      else
        THROW_EXCEPTION("invalid type");
      break;
    default:
        THROW_EXCEPTION("invalid type");
      break;
    }
  }
  return *this;
}

XSiconosVector& XSiconosVector::operator = (const BlockVector& vIn)
{
}


XSiconosVector& XSiconosVector::operator = (const XDenseVect& d)
{
  if(!_dense)
    THROW_EXCEPTION("the current vector is not dense.");
  if(d.size() != size())
    THROW_EXCEPTION("inconsistent size.");

  siconosBindings::copy(d, *vect.Dense);
  return *this;
}

XSiconosVector& XSiconosVector::operator = (const SparseVect& sp)
{
  if(_dense)
    THROW_EXCEPTION("current vector is not sparse.");
  if(sp.size() != size())
    THROW_EXCEPTION("inconsistent size.");

  noalias(*vect.Sparse) = sp;

  return *this;
}

XSiconosVector& XSiconosVector::operator = (const double* d)
{
  assert(_dense && "XSiconosVector::operator = double* : forbidden: the current vector is not dense.");

  siconosBindings::detail::copy(vect.Dense->size(), d, 1, getArray(), 1);
  return *this;
}

unsigned XSiconosVector::copyData(double* data) const
{
  assert(_dense && "XSiconosVector::copyData : forbidden: the current vector is not dense.");

  unsigned size = vect.Dense->size();
  siconosBindings::detail::copy(vect.Dense->size(), getArray(), 1, data, 1);
  return size;
}


//=================================
// Op. and assignment (+=, -= ... )
//=================================

XSiconosVector& XSiconosVector::operator += (const XSiconosVector& vIn)
{
  if(&vIn == this)  // alias
  {
    // Note: using this *= 2.0 is much more time-consuming.
    switch(num())
    {
    case Siconos::DENSE:
      *vect.Dense += *vect.Dense;
      break;
    case Siconos::SPARSE:
      *vect.Sparse += *vect.Sparse;
      break;
    default:
      THROW_EXCEPTION("invalid type.");
      break;
    }
    return *this;
  }

  Siconos::UBLAS_TYPE vInNum = vIn.num();
  {
    switch(num())
    {
    case Siconos::DENSE:
      switch(vInNum)
      {
      case Siconos::DENSE:
        noalias(*vect.Dense) += *vIn.dense();
        break;
      case Siconos::SPARSE:
        noalias(*vect.Dense) += *vIn.sparse();
        break;
      default:
        THROW_EXCEPTION("invalid type");
        break;
      }
      break;
    case Siconos::SPARSE:
      if(vInNum == Siconos::SPARSE)
        noalias(*vect.Sparse) += *vIn.sparse();
      else
        THROW_EXCEPTION("can not add a dense to a sparse.");
      break;
    default:
      THROW_EXCEPTION("invalid type.");
      break;
    }
  }
  return *this;
}
XSiconosVector& XSiconosVector::operator += (const BlockVector& vIn)
{

  return *this;
}

XSiconosVector& XSiconosVector::operator -= (const XSiconosVector& vIn)
{
  if(&vIn == this)
  {
    this->zero();
    return *this;
  }

  Siconos::UBLAS_TYPE vInNum = vIn.num();
  {
    switch(num())
    {
    case Siconos::DENSE:
      switch(vInNum)
      {
      case Siconos::DENSE:
        noalias(*vect.Dense) -= *vIn.dense();
        break;
      case Siconos::SPARSE:
        noalias(*vect.Dense) -= *vIn.sparse();
        break;
      default:
        THROW_EXCEPTION("invalid type.");
        break;
      }
      break;
    case Siconos::SPARSE:
      if(vInNum == Siconos::SPARSE)
        noalias(*vect.Sparse) -= *vIn.sparse();
      else
        THROW_EXCEPTION("can not sub a dense to a sparse.");
      break;
    default:
      THROW_EXCEPTION("invalid type.");
      break;
    }
  }
  return *this;
}

XSiconosVector& XSiconosVector::operator -= (const BlockVector& vIn)
{
  return *this;
}


//===============
// Comparison
//===============

bool operator == (const XSiconosVector &m, const XSiconosVector &x)
{
  DEBUG_PRINTF("norm = %12.8e \n", (m - x).normInf());
  DEBUG_PRINTF("std::numeric_limits<double>::epsilon() = %12.8e \n", std::numeric_limits<double>::epsilon());
  DEBUG_EXPR(std::cout << std::boolalpha << ((m - x).normInf() <= std::numeric_limits<double>::epsilon()) <<std::endl;);
  double atol = 1e-14;
  double rtol = std::numeric_limits<double>::epsilon();
  return ((m - x).normInf() <= atol + rtol * x.normInf()) ;
}

//==================
// y = scalar * x
//==================

XSiconosVector operator * (const  XSiconosVector&m, double d)
{
  Siconos::UBLAS_TYPE numM = m.num();

  if(numM == Siconos::DENSE)
  {
    // Copy m into p and call siconosBindings::scal(d,p), p = d*p.
    XDenseVect p = *m.dense();
    siconosBindings::scal(d, p);
    return p;
  }
  else// if(numM==Siconos::SPARSE)
  {
    return (SparseVect)(*m.sparse() * d);
  }
}

XSiconosVector operator * (double d, const  XSiconosVector&m)
{
  Siconos::UBLAS_TYPE numM = m.num();

  if(numM == Siconos::DENSE)
  {
    // Copy m into p and call siconosBindings::scal(d,p), p = d*p.
    XDenseVect p = *m.dense();
    siconosBindings::scal(d, p);
    return p;
  }
  else// if(numM==Siconos::SPARSE)
  {
    return (SparseVect)(*m.sparse() * d);
  }
}

XSiconosVector operator / (const XSiconosVector &m, double d)
{
  Siconos::UBLAS_TYPE numM = m.num();

  if(numM == Siconos::DENSE)
  {
    XDenseVect p = *m.dense();
    siconosBindings::scal((1.0 / d), p);
    return p;
  }

  else// if(numM==Siconos::SPARSE){
    return (SparseVect)(*m.sparse() / d);
}

//====================
//  Vectors addition
//====================

XSiconosVector operator + (const  XSiconosVector& x, const  XSiconosVector& y)
{
  if(x.size() != y.size())
    THROW_EXCEPTION("inconsistent sizes");

  Siconos::UBLAS_TYPE numX = x.num();
  Siconos::UBLAS_TYPE numY = y.num();

  if(numX == numY)  // x, y XSiconosVector of the same type
  {
    if(numX == Siconos::DENSE)
    {
      //    siconosBindings::xpy(*x.dense(),p);
      //    return p;
      return (XDenseVect)(*x.dense() + *y.dense());
    }
    else
      return (SparseVect)(*x.sparse() + *y.sparse());
  }

  else // x, y XSiconosVector with y and x of different types
  {
    if(numX == Siconos::DENSE)
      return (XDenseVect)(*x.dense() + *y.sparse());
    else
      return (XDenseVect)(*x.sparse() + *y.dense());
  }

}

void add(const XSiconosVector& x, const XSiconosVector& y, XSiconosVector& z)
{
  // Computes z = x + y in an "optimized" way (in comparison with operator +)

  if(x.size() != y.size() || x.size() != z.size())
    THROW_EXCEPTION("inconsistent sizes");

  Siconos::UBLAS_TYPE numX = x.num();
  Siconos::UBLAS_TYPE numY = y.num();
  Siconos::UBLAS_TYPE numZ = z.num();

  if(&z == &x)  // x, and z are the same object.
  {
    z += y;
  }
  else if(&z == &y)  // y and z are the same object, different from x
  {
    z += x;
  }
  else // No common memory between x,y and z
  {

    if(numZ != 0)  // z is a XSiconosVector
    {
      if(numX == numY && numX != 0)  // x, y XSiconosVector of the same type
      {
        if(numX == Siconos::DENSE)
        {
          if(numZ != Siconos::DENSE)
            THROW_EXCEPTION("Addition of two dense vectors into a sparse.");
          noalias(*z.dense()) = *x.dense() + *y.dense() ;
        }
        else
        {
          if(numZ == Siconos::DENSE)
            noalias(*z.dense()) = *x.sparse() + *y.sparse() ;
          else
            noalias(*z.sparse()) = *x.sparse() + *y.sparse() ;
        }
      }
      else if(numX != 0 && numY != 0)  // x and y of different types => z must be dense.
      {
        if(numZ != Siconos::DENSE)
          THROW_EXCEPTION("z can not be sparse.");
        if(numX == Siconos::DENSE)
          noalias(*z.dense()) = *x.dense() + *y.sparse();
        else
          noalias(*z.dense()) = *x.sparse() + *y.dense() ;
      }
    }
  }
}

//======================
//  Vectors subtraction
//======================

XSiconosVector operator - (const  XSiconosVector& x, const  XSiconosVector& y)
{
  if(x.size() != y.size())
    THROW_EXCEPTION("inconsistent sizes");

  Siconos::UBLAS_TYPE numX = x.num();
  Siconos::UBLAS_TYPE numY = y.num();

  if(numX == numY)  // x, y XSiconosVector of the same type
  {
    if(numX == Siconos::DENSE)
    {
      //    siconosBindings::xpy(*x.dense(),p);
      //    return p;
      return (XDenseVect)(*x.dense() - *y.dense());
    }
    else
      return (SparseVect)(*x.sparse() - *y.sparse());
  }
  else // x, y XSiconosVector with y and x of different types
  {
    if(numX == Siconos::DENSE)
      return (XDenseVect)(*x.dense() - *y.sparse());
    else
      return (XDenseVect)(*x.sparse() - *y.dense());
  }
}

void sub(const XSiconosVector& x, const XSiconosVector& y, XSiconosVector& z)
{
  // Computes z = x - y in an "optimized" way (in comparison with operator +)

  if(x.size() != y.size() || x.size() != z.size())
    THROW_EXCEPTION("inconsistent sizes");

  Siconos::UBLAS_TYPE numX = x.num();
  Siconos::UBLAS_TYPE numY = y.num();
  Siconos::UBLAS_TYPE numZ = z.num();

  if(&z == &x)  // x and z are the same object.
  {
    z -= y;
  }
  else if(&z == &y)  // y and z are the same object
  {
    {
      if(numX == Siconos::DENSE)
      {
        if(numZ != Siconos::DENSE)
          THROW_EXCEPTION("Subtraction of two dense vectors into a sparse.");
        *z.dense() = *x.dense() - *y.dense() ;
      }
      else
      {
        if(numZ == Siconos::DENSE)
          *z.dense() = *x.sparse() - *y.dense() ;
        else
          *z.sparse() = *x.sparse() - *y.sparse() ;
      }
    }
  }
  else // No common memory between x or y and z
  {

    if(numZ != 0)  // z is a XSiconosVector
    {
      if(numX == numY && numX != 0)  // x, y XSiconosVector of the same type
      {
        if(numX == Siconos::DENSE)
        {
          if(numZ != Siconos::DENSE)
            THROW_EXCEPTION("Addition of two dense vectors into a sparse.");
          noalias(*z.dense()) = *x.dense() - *y.dense() ;
        }
        else
        {
          if(numZ == Siconos::DENSE)
            noalias(*z.dense()) = *x.sparse() - *y.sparse() ;
          else
            noalias(*z.sparse()) = *x.sparse() - *y.sparse() ;
        }
      }
      else if(numX != 0 && numY != 0)  // x and y of different types => z must be dense.
      {
        if(numZ != Siconos::DENSE)
          THROW_EXCEPTION("z can not be sparse.");
        if(numX == Siconos::DENSE)
          noalias(*z.dense()) = *x.dense() - *y.sparse();
        else
          noalias(*z.dense()) = *x.sparse() - *y.dense() ;
      }
    }
  }
}

void axpby(double a, const XSiconosVector& x, double b, XSiconosVector& y)
{
  // Computes y = ax + by

  if(x.size() != y.size())
    THROW_EXCEPTION("inconsistent sizes");

  Siconos::UBLAS_TYPE numX = x.num();
  Siconos::UBLAS_TYPE numY = y.num();

  if(numX == numY)  // x and y of the same type
  {
    if(numX == Siconos::DENSE)  // all dense
    {
      siconosBindings::scal(b, *y.dense());
      siconosBindings::axpy(a, *x.dense(), *y.dense());
    }
    else // all sparse
    {
      *y.sparse() *= b;
      if(&y != &x)
        noalias(*y.sparse()) += a**x.sparse();
      else
        *y.sparse() += a**x.sparse();
    }
  }

  else // x and y of different types
  {
    y *= b;
    {
      if(numX == Siconos::DENSE)
        *y.sparse() += a**x.dense();
      else
        *y.dense() +=  a**x.sparse();
    }
  }
}

void axpy(double a, const XSiconosVector& x, XSiconosVector& y)
{
  // Computes y = ax + y

  if(x.size() != y.size())
    THROW_EXCEPTION("nconsistent sizes");

  Siconos::UBLAS_TYPE numX = x.num();
  Siconos::UBLAS_TYPE numY = y.num();

  if(numX == numY)  // x and y of the same type
  {
    if(numX == Siconos::DENSE)  // all dense
      siconosBindings::axpy(a, *x.dense(), *y.dense());

    else // all sparse
    {
      if(&y != &x)
        noalias(*y.sparse()) += a**x.sparse();
      else
        *y.sparse() += a**x.sparse();
    }
  }

  else // x and y of different types
  {
    {
      if(numX == Siconos::DENSE)
        *y.sparse() += a**x.dense();
      else
        *y.dense() +=  a**x.sparse();
    }
  }
}

double inner_prod(const XSiconosVector &x, const XSiconosVector &m)
{
  if(x.size() != m.size())
    THROW_EXCEPTION("inconsistent sizes");

  Siconos::UBLAS_TYPE numM = m.num();
  Siconos::UBLAS_TYPE numX = x.num();

  if(numX == numM)
  {
    if(numM == Siconos::DENSE)
      return siconosBindings::dot(*x.dense(), *m.dense());
    else
      return inner_prod(*x.sparse(), *m.sparse());
  }
  else if(numM == Siconos::DENSE)
    return inner_prod(*x.sparse(), *m.dense());
  else
    return inner_prod(*x.dense(), *m.sparse());
}

// outer_prod(v,w) = trans(v)*w
SimpleMatrix outer_prod(const XSiconosVector &x, const XSiconosVector& m)
{
  Siconos::UBLAS_TYPE numM = m.num();
  Siconos::UBLAS_TYPE numX = x.num();

  if(numM == Siconos::DENSE)
  {
    if(numX == Siconos::DENSE)
      return (DenseMat)(outer_prod(*x.dense(), *m.dense()));

    else// if(numX == Siconos::SPARSE)
      return (DenseMat)(outer_prod(*x.sparse(), *m.dense()));
  }
  else // if(numM == Siconos::SPARSE)
  {
    if(numX == Siconos::DENSE)
      return (DenseMat)(outer_prod(*x.dense(), *m.sparse()));

    else //if(numX == Siconos::SPARSE)
      return (DenseMat)(outer_prod(*x.sparse(), *m.sparse()));
  }
}

void scal(double a, const XSiconosVector & x, XSiconosVector & y, bool init)
{
  // To compute y = a *x (init = true) or y += a*x (init = false)

  if(&x == &y)
  {
    if(init)
      y *= a;
    else
    {
      y *= (1.0 + a);
    }
  }
  else
  {
    unsigned int sizeX = x.size();
    unsigned int sizeY = y.size();

    if(sizeX != sizeY)
      THROW_EXCEPTION("sizes are not consistent.");

    Siconos::UBLAS_TYPE numY = y.num();
    Siconos::UBLAS_TYPE numX = x.num();
    if(numX == numY)
    {

      if(numX == Siconos::DENSE)  // ie if both are Dense
      {
        if(init)
          //siconosBindings::axpby(a,*x.dense(),0.0,*y.dense());
          noalias(*y.dense()) = a * *x.dense();
        else
          noalias(*y.dense()) += a * *x.dense();
      }
      else  // if both are sparse
      {
        if(init)
          noalias(*y.sparse()) = a**x.sparse();
        else
          noalias(*y.sparse()) += a**x.sparse();
      }
    }
    else
    {
      if(numY == 0 || numX == 0)  // if y or x is block
      {
        if(init)
        {
          y = x;
          y *= a;
        }
        else
        {
          XSiconosVector tmp(x);
          tmp *= a;
          y += tmp;
        }
      }
      else
      {
        if(numY == Siconos::DENSE)  // if y is dense
        {
          if(init)
            noalias(*y.dense()) = a**x.sparse();
          else
            noalias(*y.dense()) += a**x.sparse();

        }
        else
          THROW_EXCEPTION("Operation not allowed on non-dense vector.");
      }
    }
  }
}

void subscal(double a, const XSiconosVector & x, XSiconosVector & y, const Index& coord, bool init)
{
  // To compute sub_y = a *sub_x (init = true) or sub_y += a*sub_x (init = false)
  // Coord  = [r0x r1x r0y r1y];
  // subX is the sub-vector of x, for row numbers between r0x and r1x-1.
  // The same for y with riy.


  // Check dimensions
  unsigned int dimX = coord[1] - coord[0];
  unsigned int dimY = coord[3] - coord[2];
  if(dimY != dimX)
    THROW_EXCEPTION("inconsistent sizes between (sub)x and (sub)y.");
  if(dimY > y.size() || dimX > x.size())
    THROW_EXCEPTION("input index too large.");

  Siconos::UBLAS_TYPE numY = y.num();
  Siconos::UBLAS_TYPE numX = x.num();

  if(&x == &y)  // if x and y are the same object
  {
    if(numX == Siconos::DENSE)  // Dense
    {
      ublas::vector_range<XDenseVect> subY(*y.dense(), ublas::range(coord[2], coord[3]));
      if(coord[0] == coord[2])
      {
        if(init)
          subY *= a;
        else
          subY *= (1.0 + a);
      }
      else
      {
        ublas::vector_range<XDenseVect> subX(*x.dense(), ublas::range(coord[0], coord[1]));
        if(init)
          subY = a * subX;
        else
          subY += a * subX;
      }
    }
    else //if (numX == Siconos::SPARSE) // Sparse
    {
      ublas::vector_range<SparseVect> subY(*y.sparse(), ublas::range(coord[2], coord[3]));
      if(coord[0] == coord[2])
      {
        if(init)
          subY *= a;
        else
          subY *= (1.0 + a);
      }
      else
      {
        ublas::vector_range<SparseVect> subX(*x.sparse(), ublas::range(coord[0], coord[1]));
        if(init)
          subY = a * subX;
        else
          subY += a * subX;
      }
    }
  }
  else
  {
    if(numX == numY)
    {
      if(numX == Siconos::DENSE)  // ie if both are Dense
      {
        ublas::vector_range<XDenseVect> subX(*x.dense(), ublas::range(coord[0], coord[1]));
        ublas::vector_range<XDenseVect> subY(*y.dense(), ublas::range(coord[2], coord[3]));

        if(init)
          noalias(subY) = a * subX;
        else
          noalias(subY) += a * subX;
      }
      else  // if both are sparse
      {
        ublas::vector_range<SparseVect> subX(*x.sparse(), ublas::range(coord[0], coord[1]));
        ublas::vector_range<SparseVect> subY(*y.sparse(), ublas::range(coord[2], coord[3]));

        if(init)
          noalias(subY) = a * subX;
        else
          noalias(subY) += a * subX;
      }
    }
    else // x and y of different types ...
    {
      if(numY == Siconos::DENSE)  // y dense, x sparse
      {
        ublas::vector_range<XDenseVect> subY(*y.dense(), ublas::range(coord[2], coord[3]));
        ublas::vector_range<SparseVect> subX(*x.sparse(), ublas::range(coord[0], coord[1]));

        if(init)
          noalias(subY) = a * subX;
        else
          noalias(subY) += a * subX;
      }
      else // y sparse, x dense => fails
        THROW_EXCEPTION("Operation not allowed (y sparse, x dense).");
    }
  }
}
void cross_product(const XSiconosVector& V1, const XSiconosVector& V2, XSiconosVector& VOUT)
{
  if(V1.size() != 3 || V2.size() != 3 || VOUT.size() != 3)
    THROW_EXCEPTION("allowed only with dim 3.");

  double aux = V1.getValue(1) * V2.getValue(2) - V1.getValue(2) * V2.getValue(1);
  VOUT.setValue(0, aux);

  aux = V1.getValue(2) * V2.getValue(0) - V1.getValue(0) * V2.getValue(2);
  VOUT.setValue(1, aux);

  aux = V1.getValue(0) * V2.getValue(1) - V1.getValue(1) * V2.getValue(0);
  VOUT.setValue(2, aux);

}

//

void abs_wise(const XSiconosVector& V, XSiconosVector& Vabs)
{
  for(unsigned int it = 0; it < V.size(); ++it)
  {
    Vabs.setValue(it, std::abs(V.getValue(it)));
  };
}

//

void getMax(const XSiconosVector& V, double& maxvalue, unsigned int& idmax)
{
  maxvalue = V.getValue(0);
  idmax = 0;
  for(unsigned int it = 1; it < V.size(); ++it)
  {
    if(V.getValue(it) > maxvalue)
    {
      maxvalue = V.getValue(it);
      idmax = it;
    };
  };
}

//

void getMin(const XSiconosVector& V, double& minvalue, unsigned int& idmin)
{
  minvalue = V.getValue(0);
  idmin = 0;
  for(unsigned int it = 1; it < V.size(); ++it)
  {
    if(V.getValue(it) < minvalue)
    {
      minvalue = V.getValue(it);
      idmin = it;
    };
  };
}


// struct exp_op
// {
//   double operator()(double d) const
//   {
//     return std::exp(d);
//   }
// };

void setBlock(const XSiconosVector& vIn, SP::XSiconosVector vOut, unsigned int sizeB,
              unsigned int startIn, unsigned int startOut)
{
  unsigned int endOut = startOut + sizeB;
  Siconos::UBLAS_TYPE numIn = vIn.num();
  Siconos::UBLAS_TYPE numOut = vOut->num();
  assert(vOut->size() >= endOut && "The output vector is too small");
  if(numIn == numOut)
  {
    if(numIn == Siconos::DENSE)  // vIn / vOut are Dense
      noalias(ublas::subrange(*vOut->dense(), startOut, endOut)) = ublas::subrange(*vIn.dense(), startIn, startIn + sizeB);
    else // if(numIn == Siconos::SPARSE)// vIn / vOut are Sparse
      noalias(ublas::subrange(*vOut->sparse(), startOut, endOut)) = ublas::subrange(*vIn.sparse(), startIn, startIn + sizeB);
  }
  else // vIn and vout of different types ...
  {
    if(numIn == Siconos::DENSE)  // vIn Dense
      noalias(ublas::subrange(*vOut->sparse(), startOut, endOut)) = ublas::subrange(*vIn.dense(), startIn, startIn + sizeB);
    else // if(numIn == Siconos::SPARSE)// vIn Sparse
      noalias(ublas::subrange(*vOut->dense(), startOut, endOut)) = ublas::subrange(*vIn.sparse(), startIn, startIn + sizeB);
  }
}

unsigned int XSiconosVector::size(void) const
{
  if(!_dense)
  {
    return (vect.Sparse->size());
  }
  else
  {
    return (vect.Dense->size());
  }
}

XSiconosVector& operator *= (XSiconosVector& v, const double& s)
{
  if(v._dense)
    *v.dense() *= s;
  else
    *v.sparse() *= s;
  return v;
}


XSiconosVector& operator /= (XSiconosVector& v, const double& s)
{
  if(v._dense)
    *v.dense() /= s;
  else
    *v.sparse() /= s;
  return v;
}

XSiconosVector::iterator XSiconosVector::begin()
{
//  return XSiconosVector::iterator(*this, 0);
}

XSiconosVector::const_iterator XSiconosVector::begin() const
{
}

XSiconosVector::iterator XSiconosVector::end()
{
}

XSiconosVector::const_iterator XSiconosVector::end() const
{
}

XSiconosVector::operator std::vector<double>()
{
  std::vector<double> v(size());
  std::copy(begin(), end(), v.begin());
  return v;
}
