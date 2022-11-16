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
    static auto solve = []<bool over=true>(const auto& m, auto&x)
      constexpr -> decltype(auto)
    {
      return ground::itransform(
        x,
        [&m,&x](const auto i,
                const match::scalar auto& xi)
        {
          static_assert(match::scalar<std::decay_t<decltype(m[0])>>);
          return xi / m[i*(std::size(x)+1)]; // diagonal terms
        });
    };

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
