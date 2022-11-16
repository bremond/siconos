#ifndef SICONOS_GROUND
#define SICONOS_GROUND

#include <algorithm>
#include <numeric>


#include <boost/hana/fwd/map.hpp>
#include <boost/hana/pair.hpp>
#include <boost/hana.hpp>
#include <boost/hana/fwd/fold_left.hpp>
#include <boost/hana/ext/std/tuple.hpp>
#include <boost/hana/ext/std/array.hpp>
#include <boost/hana/fwd/find_if.hpp>
#include <boost/hana/fwd/for_each.hpp>
#include <boost/hana/experimental/type_name.hpp>


namespace siconos
{
  namespace ground
  {
    namespace hana = boost::hana;

    template<std::size_t N>
    static auto range = hana::range_c<std::size_t, 0, N>;

    template<std::size_t N>
    static auto iterate = hana::iterate<N>;

    static auto fold_left = []
      <typename Array, typename State, typename Fun>
      (Array&& array,
       State&& initial_state,
       Fun&& fun)
      constexpr -> decltype(auto)
    {
      using array_type = std::decay_t<Array>;

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


    static auto overload = hana::overload;

    static auto apply = hana::apply;

    static auto find_if = hana::find_if;

    static auto for_each = hana::for_each;

    static auto insert = hana::insert;

    static auto partial = hana::partial;

    static auto first = hana::first;

    static auto second = hana::second;

    // f(T{}, ...) -> f<T>(...)
    template<typename T>
    static auto t_arg = []<typename F>(F&& f) { return ground::partial(f, T{}); };

    template<typename First, typename Second>
    using pair = hana::pair<decltype(hana::type_c<First>), Second>;

    static auto make_pair = hana::make_pair;

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

    // map -> tuple -> tranform -> map
    static auto map_transform = hana::demux(hana::to<hana::map_tag>)
        (hana::compose(transform, hana::to<hana::tuple_tag>));

    // dup(f)(x) = f(x, x)
    static auto dup = []<typename F>(F&& f)
      constexpr -> decltype(auto)
    {
      return
      [&f]<typename X>(X&& x)
      {
        auto&& px = std::forward<X>(x);     // x must be forwarded once!!
        return std::forward<F>(f)(px, px);
      };
    };

    static_assert(dup(hana::plus)(1) == 2);
    static_assert(dup(hana::mult)(2) == 4);

    // map_transform pair(first, f(first, second)),
    static auto map_value_transform =
      []<typename M, typename F>(M&& m, F&& f)
      constexpr -> decltype(auto)
    {
      return map_transform(std::forward<M>(m),
                           dup(hana::lockstep(hana::make_pair)
                               (hana::first,
                                dup(hana::lockstep(std::forward<F>(f))
                                    (hana::first,
                                     hana::second)))));

    };


    // compile-time itransform
    static constexpr const auto itransform_ct(const auto& a, auto&& f)
    {
      using array_type = std::decay_t<decltype(a)>;
      using size_type = typename array_type::size_type;
      array_type ta; // 'a' passed as const ref is not a constant expression
      return
        [&f, &a]<size_type... I>(std::index_sequence<I...>)
        {
          return (array_type{f(I, a[I])... });
        }(std::make_integer_sequence<size_type, std::size(ta)>{});
    }

    static constexpr const auto itransform(const auto& array,
                                           auto&& func)
    {
      using array_type = std::decay_t<decltype(array)>;
      if constexpr (hana::Foldable<array_type>::value)
      {
        return itransform_ct(array, func);
      }
      else if constexpr (
        std::equality_comparable<decltype(array.begin())> &&
        std::input_iterator<decltype(array.begin())>)
      {
        array_type res;
        using size_type = typename array_type::size_type;
        std::transform(array.begin(), array.end(),
                       std::back_inserter(res),
                       [&func,&res](const auto& x)
                       {
                         size_type i = std::size(res);
                         return func(i, x);
                       });
        return res;
      }
      else
      {
        // cf https://stackoverflow.com/questions/38304847/constexpr-if-and-static-assert
        []<bool flag = false>()
          {
            static_assert(flag, "cannot transform with these parameters");
          }();
      }
    }

  }
}
#endif
