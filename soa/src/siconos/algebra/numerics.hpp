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

    namespace match
    {
      template<typename T>
      concept any_mat = requires { typename T::any_mat_t; };

      template<typename T>
      concept mat = requires { typename T::mat_t; };

      template<typename T>
      concept diag_mat = requires { typename T::diag_mat_t; };
    }

    struct any_mat
    {
      using any_mat_t = void;
    };

    template<typename T>
    struct mat : any_mat
    {
      using mat_t = void;
      static constexpr auto value_nb_rows = traits::get_nb_rows(T{});
      static constexpr auto value_nb_cols = traits::get_nb_cols(T{});
      NumericsMatrix* _m = nullptr;

      constexpr mat(){}

      ~mat()
      {
        if(_m)
        {
          _m = NM_free(_m);
        }
      }
    };

    template<typename T>
    struct diag_mat  : mat<T>
    {
      using any_mat_t = typename mat<T>::any_mat_t;
      using diag_mat_t = void;
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


  void resize(match::any_mat auto& m,
              siconos::match::indice auto nrows,
              siconos::match::indice auto ncols)
  {
    if (m._m) m._m = NM_free(m._m);

    m._m =NM_create(NM_SPARSE,
                    nrows*m.value_nb_rows,
                    ncols*m.value_nb_cols);
    m._m->matrix2->origin = NSM_CSC;
    NM_csc_alloc(m._m, 10);
  }

  template<typename T>
  void set_value(match::any_mat auto& m,
                 siconos::match::indice auto i,
                 siconos::match::indice auto j,
                 const T& value)
  {
    if constexpr (siconos::match::scalar<T>)
    {
      NM_zentry(m._m, i, j, value, zero_threshold);
    }
    else if constexpr (siconos::match::vector<T>)
    {
      for (decltype(i) k=0; k<traits::get_nb_rows(T{}); ++k)
      {
        NM_entry(m._m, i+k, j+k, value[k]);
      }
    }
    else if constexpr (siconos::match::matrix<T>)
    {
      for (decltype(i) k=0; k<traits::get_nb_rows(T{}); ++k)
      {
        for (decltype(j) l=0; l<traits::get_nb_cols(T{}); ++l)
        {
          NM_zentry(m._m, i+k, j+l, value[k][l], zero_threshold);
        }
      }
    }
  }


    template<typename A>
    void factorize(match::diag_mat auto& a)
    {
    }

    // c <- a^-1 b^t
    template<siconos::match::matrix B>
    void solvet(match::diag_mat auto& a,
                mat<B> & b,
                mat<trans_t<B>>& c)
    {
    }

    // c <- a b
  template<typename A, typename B>
  void prod(mat<A>& a,
            mat<B>& b,
            mat<prod_t<A, B>>& c)
  {
    NM_gemm(1, a._m, b._m, 1, c._m);
  }

  } // namespace numerics
} // namespace siconos
