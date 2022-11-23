#ifndef SICONOS_LINEAR_ALGEBRA_HPP
#define SICONOS_LINEAR_ALGEBRA_HPP


#include "siconos_pattern.hpp"

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
        []<match::vector V>(V& m, V&x)
      constexpr -> decltype(auto)
      {
        for (auto [xi,mi] : ranges::views::zip(m, x))
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
    (const match::vector auto& a, const match::vector auto& b)
    {
      return ground::itransform(a,
                                [&b](const auto i,
                                     const auto& ai)
                                { return ai-b[i];});
    }

  static constexpr decltype(auto) operator -
    (const match::vector auto&a)
    {
      return ground::itransform(a,
                                [](const auto i, const auto& ai)
                                { return -ai;});
    }


  static constexpr decltype(auto) operator +
    (const auto &a, const auto& b)
    {
      return ground::itransform(a,
                                [&b](const auto i, const auto& ai)
                                { return ai+b[i];});
    }

  static constexpr decltype(auto) operator *
    (const auto& v, const auto& a)
    {
      return ground::itransform(a,
                                [&v](const auto i, const auto& ai)
                                { return v*ai;});
    }

}

#endif
