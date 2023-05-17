#pragma once

#include "CSparseMatrix_internal.h"  // for CSparseMatrix
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "siconos/algebra/eigen.hpp"
#include "siconos/algebra/linear_algebra.hpp"
#include "siconos/utils/pattern.hpp"
#include "siconos/utils/traits.hpp"

namespace siconos {

static constexpr auto zero_threshold = 1e-30;
namespace numerics {

namespace match {

template <typename T>
concept vec = requires { typename T::vec_t; };

template <typename T>
concept any_mat = !vec<T> && requires { typename T::any_mat_t; };

template <typename T>
concept mat = any_mat<T> && requires { typename T::mat_t; };

template <typename T>
concept diag_mat = any_mat<T> && requires { typename T::diag_mat_t; };

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
  bool _inversed = false;         // for diagonal format
  NumericsMatrix* _mt = nullptr;  // transposed is allocated with Numerics
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
  using vec_t = void;
  static constexpr auto vncols = 1;
  static constexpr auto vnrows = T::RowsAtCompileTime;

  NumericsMatrix* _v = nullptr;

  constexpr vec() {}

  ~vec()
  {
    if (_v) {
      _v = NM_free(_v);
    }
  }
};

const auto size0(match::any_mat auto& m) { return m._m->size0 / m.vnrows; };
const auto size0(match::vec auto& v) { return v._v->size0 / v.vnrows; };
const auto size1(match::any_mat auto& m) { return m._m->size1 / m.vncols; };

static_assert(vec<vector<int, 2>>::vncols == 1);

void resize(match::any_mat auto& m, siconos::match::indice auto nrows,
            siconos::match::indice auto ncols)
{
  if (m._m) m._m = NM_free(m._m);
  m._mt = nullptr;
  m._inversed = false;

  m._m = NM_create(NM_SPARSE, nrows * m.vnrows, ncols * m.vncols);
  NM_triplet_alloc(m._m, 1);
}

void resize(match::vec auto& v, siconos::match::indice auto nrows)
{
  if (v._v) v._v = NM_free(v._v);

  //  static_assert(m.vncols == 1);  // only vector of vectors
  // dense vector
  v._v = NM_create(NM_DENSE, nrows * v.vnrows, v.vncols);
}

void transpose(match::any_mat auto& m)
{
  if (!m._mt) {
    m._mt = NM_transpose(m._m);
  }
}

void setup(match::any_mat auto& m) { resize(m, 1, 1); }

template <typename T>
void set_value(match::any_mat auto& m, siconos::match::indice auto i,
               siconos::match::indice auto j, const T& value)
{
  if constexpr (siconos::match::scalar<T>) {
    NM_zentry(m._m, i * m.vnrows, j * m.vncols, value, zero_threshold);
  }
  // diagonal block
  else if constexpr (siconos::match::diagonal_matrix<T>) {
    for (decltype(i) k = 0; k < traits::ncols(T{}); ++k) {
      NM_zentry(m._m, i * m.vnrows + k, j * m.vncols + k, value.diagonal()(k),
                zero_threshold);
    }
  }
  // full block
  else if constexpr (siconos::match::matrix<T>) {
    for (decltype(i) k = 0; k < traits::nrows(T{}); ++k) {
      for (decltype(j) l = 0; l < traits::ncols(T{}); ++l) {
        NM_zentry(m._m, i * m.vnrows + k, j * m.vncols + l, value(k, l),
                  zero_threshold);
      }
    }
  }
  else {
    []<bool flag = false>()
    {
      static_assert(flag, "set_value: cannot insert this value");
    }
    ();
  }
}

template <typename T>
void set_value(match::vec auto& m, siconos::match::indice auto i,
               const T& value)
{
  if constexpr (siconos::match::scalar<T>) {
    NM_zentry(m._m, i * m.vnrows, 0, value, zero_threshold);
  }
  // vector block
  else if constexpr (siconos::match::vector<T>) {
    for (decltype(i) k = 0; k < traits::nrows(T{}); ++k) {
      NM_zentry(m._v, i * m.vnrows + k, 0, value(k), zero_threshold);
    }
  }
  // compile time error
  else {
    []<bool flag = false>()
    {
      static_assert(flag, "set_value: cannot insert this value");
    }
    ();
  }
}

template <siconos::match::diagonal_matrix A>
void inverse(diag_mat<A>& a)
{
  if (!a._inversed) {
    for (auto i = 0; i < NM_triplet(a._m)->nz; ++i)
      a._m->matrix2->triplet->x[i] = 1.0 / a._m->matrix2->triplet->x[i];
  }
  a._inversed = true;
}

// b += a
template <typename T>
void add(vec<T>& a, vec<T>& b)
{
// improve Numerics  cblas_daxpy(size0(a), 1, a._v->matrix0, 1, b._v->matrix0);

  for (size_t i=0; i<size0(a)*a.vnrows; ++i)
  {
    b._v->matrix0[i] += a._v->matrix0[i];
  }
}

  template<typename T>
  void scal(siconos::match::scalar auto h, vec<T>& v)
  {
    NM_scal(h, v._v);
  }
// c <- a b
// Matrix Matrix
template <typename A, typename B>
void prod(mat<A>& a, mat<B>& b, mat<prod_t<A, B>>& c)
{
  NM_gemm(1, a._m, b._m, 1, c._m);
}

// Matrix Vector
template <typename A, typename B>
void prod(mat<A>& a, vec<B>& b, vec<prod_t<A, B>>& c)
{
  NM_gemv(1, a._m, b._v->matrix0, 1, c._v->matrix0);
}

// c <- a^t b
template <typename A, typename B>
void prodt1(mat<A>& a, vec<B>& b, vec<prod_t<trans_t<A>, B>>& c)
{
  assert(a._mt);
  NM_gemv(1, a._mt, b._v->matrix0, 1, c._v->matrix0);
}
// c <- a b^t
template <typename A, typename B>
void prodt2(diag_mat<A>& a, mat<B>& b, mat<prod_t<A, trans_t<B>>>& c)
{
  assert(b._mt);
  NM_gemm(1, a._m, b._mt, 1, c._m);
}

// c <- a^-1 b
template <siconos::match::diagonal_matrix A, typename B>
void solve(diag_mat<A>& a, vec<B>& b, vec<B>& c)
{
  inverse(a);
  prod(a, b, c);
}

  template <siconos::match::diagonal_matrix A, siconos::match::matrix B>
void solvet(diag_mat<A>& a, mat<B>& b, mat<trans_t<B>>& c)
{
  inverse(a);
  transpose(b);
  prodt2(a, b, c);
}

}  // namespace numerics
}  // namespace siconos
