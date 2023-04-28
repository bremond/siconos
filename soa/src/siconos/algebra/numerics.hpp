#pragma once

#include "siconos/utils/pattern.hpp"
#include "CSparseMatrix_internal.h"                   // for CSparseMatrix
#include "NumericsSparseMatrix.h"

namespace siconos {

  void fill_csc(CSparseMatrix*csc,
                match::ublas_sparse_matrix auto m,
                match::indice auto row_off,
                match::indice auto col_off,
                match::scalar auto tol)
  {
  }

  void solve(match::ublas_sparse_matrix auto a,
             match::ublas_sparse_matrix auto b,
             match::ublas_sparse_matrix auto c)
  {
  }

  void prod(match::ublas_sparse_matrix auto a,
            match::ublas_sparse_matrix auto b,
            match::ublas_sparse_matrix auto c)
  {
  }
}
