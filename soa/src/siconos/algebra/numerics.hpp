#pragma once

#include "siconos/utils/pattern.hpp"
#include "siconos/utils/traits.hpp"
#include "siconos/algebra/linear_algebra.hpp"
#include "CSparseMatrix_internal.h"                   // for CSparseMatrix
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"

namespace siconos {

  static constexpr auto zero_threshold = 1e-30;
  namespace numerics
  {
    template<match::matrix T>
    struct mat
    {
      static constexpr auto value_nb_rows = traits::get_nb_rows(T{});
      static constexpr auto value_nb_cols = traits::get_nb_cols(T{});
      NumericsMatrix* _m = nullptr;

      mat(){}

      ~mat()
      {
        if(_m)
        {
          _m = NM_free(_m);
        }
      }
    };

  template<typename T>
  struct vec
  {
    static constexpr auto value_nb_rows = traits::get_nb_rows(T{});

    std::vector<T> _data;

    ~vec()
    {
    }
  };

  template<typename T>
  struct diag_mat : vec<T>
  {
  };

  template<typename T>
  void resize(mat<T>& m,
              match::indice auto nrows,
              match::indice auto ncols)
  {
    if (m._m) m._m = NM_free(m._m);

    m._m =NM_create(NM_SPARSE,
                    nrows*m.value_nb_rows,
                    ncols*m.value_nb_cols);
    m._m->matrix2->origin = NSM_CSC;
    NM_csc_alloc(m._m, 1);
  }

  template<typename T>
  void set_value(mat<T>& m,
                 match::indice auto i,
                 match::indice auto j,
                 const T& value)
  {
    if constexpr (match::scalar<T>)
    {
      NM_zentry(m._m, i, j, value, zero_threshold);
    }
    else if constexpr (match::matrix<T>)
    {
      for (decltype(i) k=0; k<traits::get_nb_rows(T{}); ++k)
      {
        for (decltype(j) l=0; l<traits::get_nb_cols(T{}); ++l)
        {
          NM_zentry(m._m, k, l, value[k][l], zero_threshold);
        }
      }
    }
  }


// c <- a^-1 b^t
  template<match::matrix A, match::matrix B>
  void solvet(mat<A>& a,
              mat<B>& b,
              mat<trans_t<B>>& c)
  {
  }

    // c <- a b
  template<typename A, typename B>
  void prod(mat<A>& a,
            mat<B>& b,
            mat<prod_t<A, B>>& c)
  {
  }

  } // namespace numerics
} // namespace siconos
