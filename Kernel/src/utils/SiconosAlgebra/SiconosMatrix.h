/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */

/*! \file SiconosMatrix.h
    \brief SiconosMatrix class

*/

#ifndef __SiconosMatrix__
#define __SiconosMatrix__

#include "SiconosAlgebra.h"

/** Union of DenseMat pointer, TriangMat pointer BandedMat, SparseMat, SymMat, Zero and Identity mat pointers.
 */
union Mat
{
  DenseMat *Dense;
  TriangMat *Triang;
  SymMat *Sym;
  SparseMat *Sparse;
  BandedMat *Banded;
  ZeroMat *Zero;
  IdentityMat *Identity;
};

class SimpleVector;


/** Virtual Interface for matrices in Siconos
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \author NDIAYE Abdel-Aziz
 *  \date (creation) 07/21/2006
 *  Matrices can be either block or Simple.
 *  See Derived classes for details.
 */
class SiconosMatrix
{
protected:
  /**\var isBlockMatrix
   * bool to check the type of the current matrix; true if block else false. */
  bool isBlockMatrix;

  /** default constructor */
  SiconosMatrix(bool = false);

public:

  /** Destructor. */
  virtual ~SiconosMatrix();

  /** true if the matrix is block else false.
  * \return a bool.*/
  inline bool isBlock(void) const
  {
    return isBlockMatrix;
  }

  /** determines if the matrix is square
  *  \return true if the matrix is square
  */
  virtual bool isSquare() const = 0;

  /** determines if the matrix has been inversed
  *  \return true if the matrix is inversed
  */
  virtual bool isInversed() const  = 0;

  /** determines if the matrix has been factorized
  *  \return true if the matrix is factorized
  */
  virtual bool isFactorized() const  = 0;

  /** get the attribute num of current matrix
  * \return an unsigned int.
  */
  virtual unsigned int getNum(void) const = 0;

  /** get the number of block (i=0, row, i=1 col)
  *  \param unsigned int(i=0, row, i=1 col)
  *  \return an unsigned int. 1 as default, for SimpleMatrix.
  */
  virtual inline unsigned int getNumberOfBlocks(unsigned int) const
  {
    return 1;
  };

  /** get DenseMat matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return a DenseMat
  */
  virtual const DenseMat getDense(unsigned int = 0, unsigned int = 0) const = 0;

  /** get TriangMat matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return a TriangMat
  */
  virtual const TriangMat getTriang(unsigned int = 0, unsigned int = 0) const  = 0;

  /** get SymMat matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return a SymMat
  */
  virtual const SymMat getSym(unsigned int = 0, unsigned int = 0) const  = 0;

  /** get BandedMat matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return a BandedMat
  */
  virtual const BandedMat getBanded(unsigned int = 0, unsigned int = 0) const  = 0;

  /** get SparseMat matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return a SparseMat
  */
  virtual const SparseMat getSparse(unsigned int = 0, unsigned int = 0) const = 0;

  /** get ZeroMat matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return a ZeroMat
  */
  virtual const ZeroMat getZero(unsigned int = 0, unsigned int = 0) const = 0;

  /** get  getIdentity matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return an IdentityMat
  */
  virtual const IdentityMat getIdentity(unsigned int = 0, unsigned int = 0) const = 0;

  /** get a pointer on DenseMat matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return a DenseMat*
  */
  virtual  DenseMat* getDensePtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** get a pointer on TriangMat matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return a TriangMat*
  */
  virtual TriangMat* getTriangPtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** get a pointer on SymMat matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return a SymMat*
  */
  virtual SymMat* getSymPtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** get a pointer on BandedMat matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return a BandedMat*
  */
  virtual BandedMat* getBandedPtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** get a pointer on SparseMat matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return a SparseMat*
  */
  virtual SparseMat* getSparsePtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** get a pointer on ZeroMat matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return a ZeroMat*
  */
  virtual ZeroMat* getZeroPtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** get a pointer on Identity matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return an IdentityMat*
  */
  virtual IdentityMat* getIdentityPtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** return the adress of the array of double values of the matrix ( for block(i,j) if this is a block matrix)
  *  \param: row position for the required block
  *  \param: col position for the required block
  *  \return double* : the pointer on the double array
  */
  virtual double* getArray(unsigned int = 0, unsigned int = 0) const = 0;

  /** get BlocksMat matrix
  *  \useless for SimpleMatrix
  *  \return a BlocksMat
  */
  virtual const BlocksMat getAllBlocks(void) const = 0;

  /** get block corresponding to lines given in numRow and columns in numCol
  *  \param unsigned int, row index
  *  \param unsigned int, col index
  *  \param a SiconosMatrix (in-out paramater)
  */
  virtual void getBlock(unsigned int, unsigned int, SiconosMatrix&) const = 0;

  /** get block at position row-col if BlockMatrix, else if SimpleMatrix return this
  *  \param unsigned int row
  *  \param unsigned int col
  */
  virtual SiconosMatrix* getBlockPtr(unsigned int = 0, unsigned int = 0) = 0;

  /** get std::deque of bool
  *   \useless for SimpleMatrix
  *   \return a std::deque<bool>
  */
  virtual const std::deque<bool> getBlockAllocated(void) const = 0;

  /** get row index of current matrix and save it into vOut
  *  \param unsigned int: index of required line
  *  \param ref to SimpleVector: in-out parameter
  */
  virtual void getRow(unsigned int, SimpleVector&) const = 0;

  /** get column index of current matrix and save it into vOut
  *  \param unsigned int: index of required column
  *  \param ref to SimpleVector: in-out parameter
  */
  virtual void getCol(unsigned int, SimpleVector&) const = 0;

  /** set line row of the current matrix with vector v
  *  \param an unsigned int and a SimpleVector
  */
  virtual void setRow(unsigned int, const SimpleVector&) = 0;

  /** set column col of the current matrix with vector v
  *  \param an unsigned int and a SimpleVector
  */
  virtual void setCol(unsigned int, const SimpleVector&) = 0;

  /** get the number of rows or columns of the matrix
   *  \exception SiconosMatrixException
   *  \param : unsigned int, 0 for rows, 1 for columns
   *  \return an int
   */
  virtual unsigned int size(unsigned int)const = 0;

  /** resize the matrix with nbrow rows and nbcol columns, upper and lower are only useful for BandedMatrix .
  *   The existing elements of the matrix are preseved when specified.
  *  \exception SiconosMatrixException
  *  \param 2 unsigned int: number of rows and columns
  *  \param 2 unsigned int: for banded matrices
  */
  virtual void resize(unsigned int, unsigned int, unsigned int lower = 0, unsigned int upper = 0, bool = true) = 0;

  /** compute the infinite norm of the matrix
  *  \return a double
  */
  virtual const double normInf(void) const = 0;

  /** display data on standard output
  */
  virtual void display(void) const = 0;

  /** sets all the values of the matrix to 0.0
  */
  virtual void zero(void) = 0;

  /** set an identity matrix
  */
  virtual void eye(void) = 0;

  // --- MATRICES HANDLING AND OPERATORS ---

  /** copy the matrix "blockMat" into the matrix "mat" at the position (xPos, yPos)
  *  \param SiconosMatrix& : the matrix to be copied in the current matrix
  *  \param unsigned int : the line position to start the copy of the blockmatrix
  *  \param unsigned int : the column position to start the copy of the blockmatrix
  */
  virtual void matrixCopy(const SiconosMatrix&, unsigned int, unsigned int) = 0;

  // Note: in the following functions, row and col are general;
  // that means that for a SimpleMatrix m, m(i,j) is index (i,j) element but
  // for a BlockMatrix w that contains 2 SiconosMatrix of size 3
  // w(1, 4) corresponds to the element (1,1) of the second matrix.


  /** get or set the element matrix[i,j]
  *  \param an unsigned int i
  *  \param an unsigned int j
  *  \exception SiconosMatrixException
  *  \return the element matrix[i,j]
  */
  virtual double& operator()(unsigned int , unsigned int) = 0;

  /** return the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \return a double
   */
  virtual double getValue(unsigned int, unsigned int) = 0;

  /** set the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \param the value
   */
  virtual void setValue(unsigned int, unsigned int, double) = 0;

  /** get or set the element matrix[i,j]
  *  \param an unsigned int i
  *  \param an unsigned int j
  *  \exception SiconosMatrixException
  *  \return the element matrix[i,j]
  */
  virtual double operator()(unsigned int , unsigned int) const = 0;

  /** operator =
   *  \param SiconosMatrix : the matrix to be copied
   */
  virtual SiconosMatrix& operator  =(const SiconosMatrix&) = 0;

  /** operator /=
   *  \param double, a scalar
   */
  virtual SiconosMatrix& operator /=(double) = 0;

  /** operator /=
   *  \param int, a scalar
   */
  virtual SiconosMatrix& operator /=(int) = 0;

  /** operator +=
   *  \param SiconosMatrix : a matrix to add
   */
  virtual SiconosMatrix& operator +=(const SiconosMatrix&) = 0;

  /** operator -=
   *  \param SiconosMatrix : a matrix to subtract
   */
  virtual SiconosMatrix& operator -=(const SiconosMatrix&) = 0;

  /** operator *=
   *  \param double, a scalar
   */
  virtual SiconosMatrix& operator *=(double) = 0;

  /** operator * =
   *  \param int, a scalar
   */
  virtual SiconosMatrix& operator *=(int) = 0;

  /** computes an LU factorization of a general M-by-N matrix using partial pivoting with row interchanges.
  *  The result is returned in this (InPlace). Based on Blas dgetrf function.
  */
  virtual void PLUFactorizationInPlace(void) = 0;

  /**  compute inverse of this thanks to LU factorization with Partial pivoting. This method inverts U and then computes inv(A) by solving the system
  *  inv(A)*L = inv(U) for inv(A). The result is returned in this (InPlace). Based on Blas dgetri function.
  */
  virtual void  PLUInverseInPlace(void) = 0;

  /** solves a system of linear equations A * X = B  (A=this) with a general N-by-N matrix A using the LU factorization computed
  *   by PLUFactorizationInPlace. Based on Blas dgetrs function.
  *  \param input: the RHS matrix b - output: the result x
  */
  virtual void  PLUForwardBackwardInPlace(SiconosMatrix &B) = 0;

  /** solves a system of linear equations A * X = B  (A=this) with a general N-by-N matrix A using the LU factorization computed
  *   by PLUFactorizationInPlace.  Based on Blas dgetrs function.
  *  \param input: the RHS matrix b - output: the result x
  */
  virtual void   PLUForwardBackwardInPlace(SiconosVector &B) = 0;
};

#endif
