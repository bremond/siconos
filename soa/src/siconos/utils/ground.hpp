#pragma once

#include <boost/hana/fwd/type.hpp>
#define BOOST_HANA_CONFIG_ENABLE_STRING_UDL 1
#include <algorithm>
#include <boost/hana/functional/overload_linearly.hpp>
#include <numeric>
#include <type_traits>

#include <boost/hana/fwd/map.hpp>
#include <boost/hana/pair.hpp>
#include <boost/hana.hpp>
#include <boost/hana/fwd/fold_left.hpp>
#include <boost/hana/ext/std/tuple.hpp>
#include <boost/hana/integral_constant.hpp>
#include <boost/hana/equal.hpp>
#include <boost/hana/not_equal.hpp>
#include <boost/hana/ext/std/array.hpp>
#include <boost/hana/fwd/find_if.hpp>
#include <boost/hana/fwd/for_each.hpp>
#include <boost/hana/string.hpp>

#include "boost/hana/ext/boost/ublas.hpp"

#include <typeinfo>
// cf https://www.boost.org/doc/libs/1_80_0/libs/hana/doc/html/structboost_1_1hana_1_1string.html#ad77f7afff008c2ce15739ad16a8bf0a8

#include <string_view>


namespace siconos
{

  namespace ground
  {


    template <typename T>
    constexpr auto type_name() {
      std::string_view name, prefix, suffix;
#ifdef __clang__
      name = __PRETTY_FUNCTION__;
      prefix = "auto type_name() [T = ";
      suffix = "]";
#elif defined(__GNUC__)
      name = __PRETTY_FUNCTION__;
      prefix = "constexpr auto type_name() [with T = ";
      suffix = "]";
#elif defined(_MSC_VER)
      name = __FUNCSIG__;
      prefix = "auto __cdecl type_name<";
      suffix = ">(void)";
#endif
      name.remove_prefix(prefix.size());
      name.remove_suffix(suffix.size());
      return name;
    }

    // debug (see_below for clang, gcc above...)
#if defined( __clang__)
    template <typename ...Ts> constexpr bool type_trace() { std::tuple<Ts...> see_messages_below; return true;};
#else
    template <typename ...Ts> constexpr bool type_trace() { std::tuple<Ts...> see_messages_above; return true;};
#endif

    namespace hana = boost::hana;

    static constexpr auto typeid_ = hana::typeid_;

    template<template<typename ...Ts> typename F>
    static constexpr auto trait = hana::trait<F>;
    static constexpr auto compose = hana::compose;

    static constexpr auto lockstep = hana::lockstep;

    static constexpr auto any_of = hana::any_of;

    template<std::size_t N>
    static constexpr auto range = hana::range_c<std::size_t, 0, N>;

    template<std::size_t N>
    static constexpr auto iterate = hana::iterate<N>;

    static auto fold_left = []<typename Array, typename State, typename Fun>(
                                Array &&array, State &&initial_state,
                                Fun &&fun) constexpr -> decltype(auto) {
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
        }
        ();
      }
    };

    static auto overload = hana::overload_linearly;

    static auto apply = hana::apply;

    static auto find_if = hana::find_if;

    static auto filter = hana::filter;

    static auto for_each = hana::for_each;

    static auto insert = hana::insert;

    static auto partial = hana::partial;

    static auto first = hana::first;

    static auto second = hana::second;

    // f(T{}, ...) -> f<T>(...)
    template<typename T>
    static auto t_arg = []<typename F>(F&& f) { return ground::partial(f, T{}); };

    template <typename First, typename Second>
    using pair =
        hana::pair<std::decay_t<decltype(hana::type_c<First>)>, Second>;

    template <typename First, typename Second>
    using ipair = hana::pair<First, Second>;

    template <typename T> static auto type_c = hana::type_c<T>;

    static auto make_pair =
        []<typename First, typename Second>(First, Second &&second) constexpr
        -> decltype(auto) {
      return hana::make_pair(hana::type_c<First>, std::forward<Second>(second));
      };

      template <typename... Pairs> using map = hana::map<Pairs...>;

      template <typename Data, typename Key>
      concept has_key = requires(Data m) {
        m[hana::type_c<Key>];
      };

      template <typename T>
      static auto get = []<has_key<T> D>(D && data) constexpr->decltype(auto) {
        return data[hana::type_c<T>];
      };

      static auto transform = hana::transform;

      // map -> tuple -> tranform -> map
      static auto map_transform = hana::demux(hana::to<hana::map_tag>)(
          hana::compose(transform, hana::to<hana::tuple_tag>));

      // dup(f)(x) = f(x, x)
      static auto dup = []<typename F>(F && f) constexpr->decltype(auto) {
        return [&f]<typename X>(X &&x) {
          auto &&px = std::forward<X>(x); // x must be forwarded once!!
          return std::forward<F>(f)(px, px);
        };
      };

      static_assert(dup(hana::plus)(1) == 2);
      static_assert(dup(hana::mult)(2) == 4);

      // map_transform pair(first, f(first, second)),
      static auto map_value_transform =
          []<typename M, typename F>(M && m, F && f) constexpr->decltype(auto) {
        return map_transform(
            std::forward<M>(m),
            dup(hana::lockstep(hana::make_pair)(
                hana::first, dup(hana::lockstep(std::forward<F>(f))(
                                 hana::first, hana::second)))));
      };

      // compile-time itransform
      static constexpr const auto itransform_ct(const auto &a, auto &&f) {
        using array_type = std::decay_t<decltype(a)>;
        using size_type = typename array_type::size_type;
        array_type ta; // 'a' passed as const ref is not a constant expression
        return [&f, &a ]<size_type... I>(std::index_sequence<I...>) {
          return (array_type{f(I, a[I])... });
        }
        (std::make_integer_sequence<size_type, std::size(ta)>{});
      }

      static constexpr const auto itransform(const auto &array, auto &&func) {
        using array_type = std::decay_t<decltype(array)>;
        if constexpr (hana::Foldable<array_type>::value) {
          return itransform_ct(array, func);
        } else if constexpr (std::equality_comparable<
                                 decltype(array.begin())> &&
                             std::input_iterator<decltype(array.begin())>) {
          array_type res;
          using size_type = typename array_type::size_type;
          std::transform(array.begin(), array.end(), std::back_inserter(res),
                         [&func, &res](const auto &x) {
                           size_type i = std::size(res);
                           return func(i, x);
                         });
          return res;
        } else {
          // cf
          // https://stackoverflow.com/questions/38304847/constexpr-if-and-static-assert
          []<bool flag = false>() {
            static_assert(flag, "cannot transform with these parameters");
          }
          ();
        }
      }
    template <typename Base>
    struct from
    {
      template <typename T>
      struct is_a_derivation
      {
        using type = is_a_derivation<T>;
        static constexpr bool value = std::derived_from<T, Base>;
      };
    };

    template <auto F, typename... Ts>
    using check = std::conditional_t<F.template operator()<Ts...>(),
                                     std::true_type, std::false_type>;

    template <auto F, typename... Ts>
    struct on_concept {
      template <typename T2>
      struct is_a_model {
        using type = is_a_model<T2>;
        static constexpr bool value = check<F, T2, Ts...>::value;
      };
    };

    template <auto F, typename... Ts>
    static constexpr auto is_a_model =
        compose(trait<on_concept<F, Ts...>::template is_a_model>, typeid_);

    static constexpr auto is_integral =
      is_a_model <
      []<typename T>() consteval { return std::is_integral<T>::value; }>;

    template<typename B>
    static constexpr auto derive_from =
      is_a_model<
      []<typename T>() consteval { return std::derived_from<T, B>; }>;
  }
}

