#pragma once

#include "CSparseMatrix_internal.h"  // for CSparseMatrix
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "siconos/algebra/linear_algebra.hpp"
#include "siconos/utils/pattern.hpp"
#include "siconos/utils/traits.hpp"

namespace siconos {

static constexpr auto zero_threshold = 1e-30;
namespace numerics {

namespace match {
template <typename T>
concept any_mat = requires { typename T::any_mat_t; };

template <typename T>
concept mat = requires { typename T::mat_t; };

template <typename T>
concept diag_mat = requires { typename T::diag_mat_t; };
}  // namespace match

struct any_mat {
  using any_mat_t = void;
};

template <typename T>
struct mat : any_mat {
  using mat_t = void;
  static constexpr auto vncols = T::ColsAtCompileTime;
  static constexpr auto vnrows = T::RowsAtCompileTime;
  NumericsMatrix* _m = nullptr;
  bool _inversed = false; // for diagonal format
  NumericsMatrix* _mt = nullptr; // transposed is allocated with Numerics
  constexpr mat() {}

  ~mat()
  {
    if (_m) {
      _m = NM_free(_m);
    }
  }
};

template <typename T>
struct diag_mat : mat<T> {
  using any_mat_t = typename mat<T>::any_mat_t;
  using diag_mat_t = void;

  static constexpr auto vncols = T::ColsAtCompileTime;
  static constexpr auto vnrows = T::RowsAtCompileTime;

};

template <typename T>
struct vec {
  static constexpr auto vnrows = T::RowsAtCompileTime;

  std::vector<T> _data;

  ~vec() {}
};

const auto size0(match::any_mat auto& m)
{
  return m._m->size0 / m.vnrows;
};
const auto size1(match::any_mat auto& m)
{
  return m._m->size1 / m.vncols;
};

void resize(match::any_mat auto& m, siconos::match::indice auto nrows,
            siconos::match::indice auto ncols)
{
  if (m._m) m._m = NM_free(m._m);
  m._mt = nullptr;
  m._inversed = false;

  m._m =
      NM_create(NM_SPARSE, nrows * m.vnrows, ncols * m.vncols);
  NM_triplet_alloc(m._m, 1);
}

  void transpose(match::any_mat auto& m)
  {
    if (!m._mt)
    {
      m._mt = NM_transpose(m._m);
    }
  }

void setup(match::any_mat auto& m) { resize(m, 1, 1); }

template <typename T>
void set_value(match::any_mat auto& m, siconos::match::indice auto i,
               siconos::match::indice auto j, const T& value)
{
  if constexpr (siconos::match::scalar<T>) {
    NM_zentry(m._m, i*m.vnrows, j*m.vncols, value, zero_threshold);
  }
  // diagonal block
  else if constexpr (siconos::match::diagonal_matrix<T>) {
    for (decltype(i) k = 0; k < traits::ncols(T{}); ++k) {
      NM_zentry(m._m, i*m.vnrows + k, j*m.vncols + k, value.diagonal()(k), zero_threshold);
    }
  }
  // full block
  else if constexpr (siconos::match::matrix<T>) {
    for (decltype(i) k = 0; k < traits::nrows(T{}); ++k) {
      for (decltype(j) l = 0; l < traits::ncols(T{}); ++l) {
        NM_zentry(m._m, i*m.vnrows + k, j*m.vncols + l, value(k, l), zero_threshold);
      }
    }
  }
}
  template<siconos::match::diagonal_matrix A>
  void inverse(diag_mat<A>& a)
  {
    if(!a._inversed)
    {
      for (auto i=0; i<NM_triplet(a._m)->nz; ++i)
        a._m->matrix2->triplet->x[i] = 1.0/a._m->matrix2->triplet->x[i];
    }
    a._inversed = true;
  }

  // c <- a b
template <typename A, typename B>
void prod(mat<A>& a, mat<B>& b, mat<prod_t<A, B>>& c)
{
  NM_gemm(1, a._m, b._m, 1, c._m);
}

  // c <- a b^t
  template <typename A, typename B>
  void prodt(diag_mat<A>& a, mat<B>& b, mat<prod_t<A, trans_t<B>>>& c)
{
  assert(b._mt);
  NM_gemm(1, a._m, b._mt, 1, c._m);
}

// c <- a^-1 b^t
  template <siconos::match::diagonal_matrix A, siconos::match::matrix B>
  void solvet(diag_mat<A>& a, mat<B>& b, mat<trans_t<B>>& c)
{
  inverse(a);
  transpose(b);
  prodt(a, b, c);
}


}  // namespace numerics
}  // namespace siconos
