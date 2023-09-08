#pragma once

#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/traits/traits.hpp"

#include <range/v3/view/join.hpp>
#include <range/v3/view/stride.hpp>
#include <range/v3/view/zip.hpp>
#include <type_traits>


namespace siconos {

using size_t = std::size_t;

static auto solve_in_place(auto& m, auto& x)
{
  if constexpr (requires(decltype(m) m_) { m_.diagonal()[0]; }) {
    for (auto [xi, mi] : ranges::views::zip(x, m)) {
      xi = xi / mi;  // diagonal terms
    };
  }
  else
    []<bool flag = false>()
    {
      static_assert(flag, "not implemented");
      // storage::ground::type_trace<decltype(m)>();
    }
  ();
}

// static constexpr decltype(auto) operator-(const match::vector auto& a,
//                                           const match::vector auto& b)
// {
//   return ground::itransform(
//       a, [&b](const auto i, const auto& ai) { return ai - b[i]; });
// }

// static constexpr decltype(auto) operator-(const match::vector auto& a)
// {
//   return ground::itransform(a,
//                             [](const auto i, const auto& ai) { return -ai;
//                             });
// }

// static constexpr decltype(auto) operator+(const match::vector auto& a,
//                                           const match::vector auto& b)
// {
//   return ground::itransform(
//       a, [&b](const auto i, const auto& ai) { return ai + b[i]; });
// }

// static constexpr decltype(auto) operator*(const match::scalar auto& v,
//                                           const match::vector auto& a)
// {
//   return ground::itransform(
//       a, [&v](const auto i, const auto& ai) { return v * ai; });
// }

}  // namespace siconos
