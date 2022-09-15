#ifndef SICONOS_GROUND
#define SICONOS_GROUND

#include <boost/hana/fwd/map.hpp>
#include <boost/hana/pair.hpp>
#include <numeric>
#include <boost/hana.hpp>
#include <boost/hana/fwd/fold_left.hpp>
#include <boost/hana/ext/std/tuple.hpp>
#include <boost/hana/ext/std/array.hpp>
#include <boost/hana/fwd/find_if.hpp>
#include <boost/hana/fwd/for_each.hpp>



namespace siconos
{
  namespace ground
  {
    namespace hana = boost::hana;

    static auto fold_left = [](const auto& array,
                               const auto& initial_state,
                               const auto&& fun) constexpr -> decltype(auto)
    {
      using array_type = std::decay_t<decltype(array)>;

      /* ~ static */
      if constexpr (hana::Foldable<array_type>::value)
      {
        return
          hana::fold_left(array, initial_state, fun);
      }
      /* ~ dynamic (aka std::vector) */
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

    static auto find_if = hana::find_if;

    static auto for_each = hana::for_each;

    template<typename First, typename Second>
    using pair = hana::pair<decltype(hana::type_c<First>), Second>;

    template<typename ...Pairs>
    using map = hana::map<Pairs...>;

    template<typename T>
    static auto get = [](auto&& data) constexpr -> decltype(auto)
    {
//      if constexpr (hana::contains(data_t{}, hana::type_c<T>))
      {
        return data[hana::type_c<T>];
      }
//      else
//      {
        // cf https://stackoverflow.com/questions/38304847/constexpr-if-and-static-assert
//        []<bool flag = false>()
//        {
//          static_assert(flag, "map: unknown key");
//        }();
//      }
    };

    static auto transform = hana::transform;

    static auto map_transform = hana::demux(hana::to<hana::map_tag>)
        (hana::compose(transform, hana::to<hana::tuple_tag>));

    static auto dup = []<typename F>(F&& f) constexpr -> decltype(auto)
    {
      return
      [&f]<typename X>(X&& x)
      {
        auto&& px = std::forward<X>(x);     // x must be forwarded once!!
        return std::forward<F>(f)(px, px);
      };
    };

    static_assert(dup(hana::plus)(1) == 2);

    static auto transform_map_value = []<typename F>(F&& f)
    {
      return
      dup(hana::lockstep(hana::make_pair)
          (hana::first,
           hana::compose(std::forward<F>(f), hana::second)));
    };


    static auto map_vtransform = []<typename F>(F&& f)
    { return
      hana::curry<2>(hana::flip(map_transform))
      (transform_map_value(std::forward<F>(f)));
    };
  }
}
#endif
