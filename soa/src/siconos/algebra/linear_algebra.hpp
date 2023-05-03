#pragma once

#include "siconos/utils/pattern.hpp"
#include "siconos/utils/traits.hpp"
#include <range/v3/view/join.hpp>
#include <range/v3/view/stride.hpp>
#include <range/v3/view/zip.hpp>
#include <type_traits>

namespace siconos
{

  using size_t = std::size_t;

  template<typename T, size_t N>
  using vector = std::array<T, N>;

  template<typename T, size_t N, size_t M>
  using matrix = std::array<vector <T, N>, M>;

  template<match::matrix A>
  using trans_t = matrix<typename A::value_type::value_type,
                         traits::get_nb_rows(A{}),
                         traits::get_nb_cols(A{})>;

  template<match::matrix A, match::matrix B>
  using prod_t = matrix<typename A::value_type::value_type,
                        traits::get_nb_rows(A{}),
                        traits::get_nb_cols(B{})>;

  static_assert (std::is_same_v<matrix<int,1,2>,trans_t<matrix<int,2,1>>>);

  namespace linear_algebra
  {

    static auto solve_in_place =
      ground::overload(
        []<match::abstract_vector V>(V& m, V&x)
      constexpr -> decltype(auto)
      {
        for (auto [xi,mi] : ranges::views::zip(x, m))
        {
          xi = xi/mi; // diagonal terms
        };
      },
      [](auto& m, auto &x)
      {
         []<bool flag = false>()
         {
            static_assert(flag, "not implemented");
         }();
      });

  }

  static constexpr decltype(auto) operator -
    (const match::abstract_vector auto& a, const match::abstract_vector auto& b)
    {
      return ground::itransform(a,
                                [&b](const auto i,
                                     const auto& ai)
                                { return ai-b[i];});
    }

  static constexpr decltype(auto) operator -
    (const match::abstract_vector auto&a)
    {
      return ground::itransform(a,
                                [](const auto i, const auto& ai)
                                { return -ai;});
    }


  static constexpr decltype(auto) operator +
  (const match::abstract_vector auto& a, const match::abstract_vector auto& b)
    {
      return ground::itransform(a,
                                [&b](const auto i, const auto& ai)
                                { return ai+b[i];});
    }

  static constexpr decltype(auto) operator *
  (const match::scalar auto& v, const match::abstract_vector auto& a)
    {
      return ground::itransform(a,
                                [&v](const auto i, const auto& ai)
                                { return v*ai;});
    }

}
