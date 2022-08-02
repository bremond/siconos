#ifndef SICONOS_UTIL_HPP
#define SICONOS_UTIL_HPP

#include <tuple>
#include <functional>
#include <iostream>

//https://stackoverflow.com/questions/81870/is-it-possible-to-print-a-variables-type-in-standard-c/56766138#56766138
#include <string_view>

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

namespace siconos
{

  template<typename Attrs>
  static auto tags = []() constexpr
  {
    return transform([]<typename Attr>(Attr)
                     {
                       return std::tuple<typename Attr::tag>{};
                     }, Attrs{});
  };

  namespace concepts
  {
    // https://stackoverflow.com/questions/68443804/c20-concept-to-check-tuple-like-types
    // https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2021/p2165r2.pdf
    template<class T, std::size_t N>
    concept has_tuple_element =
      requires(T t)
    {
      typename std::tuple_element_t<N, std::remove_const_t<T>>;
      { std::get<N>(t) } -> std::convertible_to<const std::tuple_element_t<N, T>&>;
    };

    template<class T>
    concept tuple_like_ =
      !std::is_reference_v<T> &&
      requires(T t)
    {
      typename std::tuple_size<T>::type;
      requires std::derived_from<
      std::tuple_size<T>,
      std::integral_constant<std::size_t, std::tuple_size_v<T>>
      >;
    } && []<std::size_t... N>(std::index_sequence<N...>)
    {
      return (has_tuple_element<T, N> && ...);
    }(std::make_index_sequence<std::tuple_size_v<T>>());

    template<class T>
    concept tuple_like = tuple_like_<std::decay_t<T>>;

    template<typename T>
    concept array_like = requires (T a)
    { typename T::size_type;
      typename T::value_type;
      { std::size(a) } -> std::convertible_to<std::size_t>;
    } && []<std::size_t... N>(std::index_sequence<N...>)
    {
      return (has_tuple_element<T, N> && ...);
    }(std::make_index_sequence<std::tuple_size_v<T>>());

    template<typename T>
    concept keep = requires (T t) { typename T::tag; { t.size }; };

    template<typename T>
    concept use = requires { typename T::use_t; };

    template<typename T>
    concept use_items = requires { typename T::use_items; };

    template<typename T>
    concept tag = requires { typename T::tag; typename T::tag_t; };

    template<typename T>
    concept structure = requires { typename T::structure_t; };

    template<typename T>
    concept attribute = requires { typename T::tag; typename T::symbol; };

    template<typename T>
    concept has_attributes = requires { typename T::attributes; };

    template<typename T>
    concept terminal_attribute =
      attribute<T> &&
      requires { typename T::structure::type; };

    template<typename ...Ts>
    concept at_least_one_terminal_attribute = (terminal_attribute<Ts> || ...);

    template<typename T>
    concept has_at_least_one_terminal_attribute =
      has_attributes<T> &&
      at_least_one_terminal_attribute<typename T::attributes>;

    template<typename T>
    concept item = has_attributes<typename T::definition>;

    template<typename T>
    concept vertex_item = requires { { T::definition::vertex_item == true }; };

  }

  template<typename ...Args>
  using tuple = std::tuple<Args...>;

  static auto for_each = []<typename ...Args>(auto&& fun, const tuple<Args...> tpl)
    constexpr
  {
    std::apply([&fun](auto&&... args) { ((fun(args)), ...);}, tpl);
  };



  template<typename U>
  struct contains_p
  {
    template<typename T>
    struct make
    {
      static constexpr auto value = contains<U>(T{});
    };
  };




  // https://ctrpeach.io/posts/cpp20-string-literal-template-parameters/
  // https://stackoverflow.com/questions/62266052/c20-string-literal-template-argument-working-example
  // https://www.cppstories.com/2021/constexpr-vecstr-cpp20/

  template<std::size_t N>
  struct string_literal
  {
    constexpr string_literal(const char (&str)[N])
    {
      std::copy_n(str, N, value);
    }
    char value[N];
  };

  template<string_literal symb>
  static auto get_m =
    [](auto&& attributes)
    {
      constexpr auto search = []<typename Attrs>(auto&& loop, Attrs&& attrs)
        constexpr
      {
        using attrs_t = std::decay_t<Attrs>;
        if constexpr (std::tuple_size_v<attrs_t> == 0)
        {
          // not found
          // cf https://stackoverflow.com/questions/38304847/constexpr-if-and-static-assert
          []<bool flag = false, typename T=Attrs>()
            {
              static_assert(flag, "symbol not found");
            }();
        }
        constexpr auto attr = car(std::decay_t<Attrs>{});

        if constexpr (attr.symbol_v.data.value == symb.value)
        {
          return car(attrs);
        }
        else
        {
          return loop(loop, cdr(attrs));
        }
      };
      search(search, attributes);
    };


}
#endif
