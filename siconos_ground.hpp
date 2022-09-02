#ifndef SICONOS_GROUND
#define SICONOS_GROUND

#include <numeric>
#include <boost/hana.hpp>
#include <boost/hana/fwd/fold_left.hpp>
#include <boost/hana/ext/std/tuple.hpp>

namespace siconos
{
  namespace ground
  {
    static auto fold_left = [](const auto& array,
                               const auto& initial_state,
                               const auto&& fun) constexpr -> decltype(auto)
    {
      using array_type = std::decay_t<decltype(array)>;

      if constexpr (boost::hana::Foldable<array_type>::value)
      {
        return
          boost::hana::fold_left(array, initial_state, fun);
      }
      else if constexpr (
        std::copy_constructible<array_type> &&
        std::equality_comparable<decltype(array.begin())> &&
        std::input_iterator<decltype(array.begin())>)
      {
        return std::accumulate(array.begin(), array.end(), initial_state, fun);
      }
      else
      {
        // cf https://stackoverflow.com/questions/38304847/constexpr-if-and-static-assert
        []<bool flag = false>()
          {
            static_assert(flag, "cannot fold_left with these parameters");
          }();
      }
    };
  }
}

#endif
