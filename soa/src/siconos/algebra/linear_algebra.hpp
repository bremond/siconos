#pragma once

#include "siconos/utils/pattern.hpp"

#include <range/v3/view/join.hpp>
#include <range/v3/view/stride.hpp>
#include <range/v3/view/zip.hpp>
#include <type_traits>

namespace siconos
{
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
