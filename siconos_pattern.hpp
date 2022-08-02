#ifndef SICONOS_PATTERN_HPP
#define SICONOS_PATTERN_HPP

#include "siconos_util.hpp"

#include <tuple>
#include <functional>
#include <type_traits>
#include <initializer_list>
namespace siconos
{
  // let rec
  // https://stackoverflow.com/questions/2067988/recursive-lambda-functions-in-c11
  template <typename F>
  struct recursive
  {
    F f;
    template <typename... Ts> constexpr
    decltype(auto) operator()(Ts&&... ts)  const
    {
      return f(std::ref(*this), std::forward<Ts>(ts)...);
    }

    template <typename... Ts> constexpr
    decltype(auto) operator()(Ts&&... ts)
    {
      return f(std::ref(*this), std::forward<Ts>(ts)...);
    }
  };

  template <typename F> recursive(F) -> recursive<F>;

  static auto const rec = [](auto f) constexpr
  {
    return recursive{std::move(f)};
  };

  static auto proj = [] (auto& data)
      constexpr -> decltype(auto)
  {
    return
      ([&data](auto&& fun) constexpr -> decltype(auto)
      {
        return
          ([&data,&fun]<typename ...As>(As&&...args) constexpr -> decltype(auto)
          {
            return fun(std::forward<As>(args)..., data);
          });
      });
  };

  static auto compose = [](auto&& f, auto&& g) constexpr -> decltype(auto)
  {
    return
      [&f,&g]<typename Arg>(Arg arg) constexpr -> decltype(auto)
      {
        return f(g(std::forward<Arg>(arg)));
      };
  };

  template<typename ... Ts>
  struct overload : Ts ... {
    using Ts::operator() ...;
  };

  template<typename ...Args>
  using gather = std::tuple<Args...>;

  static auto car = []<typename Tpl>(Tpl tpl)
    constexpr
  {
    return std::get<0>(tpl);
  };

  static_assert(std::is_same_v<decltype(car(gather<int,float,char>{})), int>);

  static auto cdr =
    []<typename A0, typename... As>(std::tuple<A0, As...> tpl)
    constexpr
  {
    return std::apply(
      [](const A0&, const As&... args)
      {
        return std::make_tuple(args...);
      }, tpl);
  };

  static_assert(std::is_same_v<decltype(cdr(gather<int,float,char>{})), gather<float, char>>);

  static auto cons =
    []<typename A, typename... As>(A, std::tuple<As...>)
    constexpr
  {
    return std::tuple<A, As...>{};
  };

  static_assert(std::is_same_v<decltype(cons(int{}, gather<float,char>{})), gather<int,float, char>>);

  static auto append = []<concepts::tuple_like ...Tpls>(Tpls... tpls)
    constexpr
  {
    return std::tuple_cat(tpls...);
  };

  static auto flatten = []<concepts::tuple_like Tpl>(Tpl tpl)
      constexpr
  {
    if constexpr (concepts::tuple_like<decltype(car(tpl))>)
    {
      return std::apply(append, tpl);
    }
    else
    {
      return tpl;
    }
  };

  static auto transform = []<typename ...Args>(auto&& fun, const std::tuple<Args...>)
    constexpr
  {
    return std::tuple {fun(Args{})...};
  };

  template<typename F>
  static auto filter =
    []<typename Tpl>(const Tpl tpl) constexpr
    {
      auto loop = rec([]<typename Itpl>(auto&& loop, Itpl) constexpr
      {
        auto v = car(Itpl{});
        using v_t = std::decay_t<decltype(v)>;
        if constexpr (F{}.value(v_t{}))
        {
          if constexpr (std::tuple_size_v<Itpl> > 1)
          {
            return cons(v, loop(cdr(Itpl{})));
          }
          else
          {
            return std::tuple<v_t>(v);
          }
        }
        else if constexpr (std::tuple_size_v<Itpl> > 1)
        {
          return loop(cdr(Itpl{}));
        }
        else
        {
          return std::tuple<>{};
        }
      });

      if constexpr (std::tuple_size_v<Tpl> > 0)
      {
        return loop(Tpl{});
      }
      else
      {
        return std::tuple<>{};
      }
    };

  template<typename T>
  static constexpr auto instance = T{};

  template<typename T>
  struct hold
  {
    static constexpr auto value = T{};
  };
  template<typename T>
  static constexpr auto make_instance(T&&) { return instance<T>; };

  template<typename U>
  static auto contains = []<typename ...Attrs>(gather<Attrs...>) constexpr -> bool
  {
    return (std::is_same_v<U, Attrs> || ...);
  };

  namespace must
  {
    template<typename T, typename Tpl>
    concept contains = contains<T>(instance<Tpl>);
  }

  namespace match
  {
    template<typename T>
    concept attribute = requires { typename T::attribute_t; };

    template<typename T>
    concept attributes = requires { typename T::attributes; };

    template<typename T>
    concept vdescriptor = requires { typename T::vdescriptor_t; };

    template<typename T>
    concept item = requires { typename T::item_t; };

    template<typename T>
    concept graph_item =
      item<typename T::item> &&
      requires { typename T::graph_item_t; };

    template<typename T>
    concept vertex_item =
      graph_item<typename T::graph_item> &&
      requires { typename T::vertex_item_t; };

    template<typename T>
    concept edge_item =
      graph_item<typename T::graph_item> &&
      requires { typename T::edge_item_t; };

    template<typename T, typename I>
    concept attribute_of =
      attribute<T> &&
      item<I> &&
      must::contains<T, typename I::attributes>;

    template<typename T>
    concept uses =
      item<T> &&
      requires { typename T::uses; };

    template<typename T>
    concept keeps =
      item<T> &&
      requires { typename T::keeps; };

    template<typename T>
    concept push_back = requires (T a) { { a.push_back(typename T::value_type{}) }; };
  }

  static_assert(match::push_back<std::vector<double>>);
  static_assert(!match::push_back<std::array<double,3>>);

  namespace some
  {
    template<size_t ...Sizes>
    struct attribute
    {
      using attribute_t = void;

      static constexpr auto sizes = std::tuple {Sizes...};
    };

    struct scalar : attribute<>
    {};

    struct indice : attribute<>
    {};

    template<typename T>
    struct vdescriptor : attribute<>
    {
      using vdescriptor_t = void;
      static constexpr bool descriptor = true;
      using type = T;
    };

    template<std::size_t N, std::size_t M>
    struct matrix : attribute<N, M> {};

    template<std::size_t N>
    struct vector : attribute<N> {};

  }

  template<string_literal Text>
  struct text
  {
    static constexpr auto data = Text;
  };

  template<string_literal Symbol>
  struct symbol : text<Symbol>{};;

  template<string_literal Descr>
  struct description : text<Descr>{};

  template<typename T>
  struct tag
  {
    using type = T;
  };

  // concept item
  template<match::item T>
  struct use
  {
    using use_t = void;
    using type = T;
  };


  template<typename T>
  struct attribute_p
  {
    static constexpr auto value = match::attribute<T>;
  };

  template<match::attribute Attr, std::size_t N>
  struct keep
  {
    using attribute = Attr;
    static constexpr std::size_t size = N;
  };

  template<typename ...Args>
  struct item
  {
    using args = gather<Args...>;
    using item_t = void;
  };

  template<typename ...Args>
  struct graph_item : item<Args...>
  {
    using graph_item_t = void;
  };

  template<typename ...Args>
  struct vertex_item : graph_item<Args...>
  {
    using vertex_item_t = void;
  };

  template<typename ...Args>
  struct edge_item : graph_item<Args...>
  {
    using edge_item_t = void;
  };


  template<typename... Ts>
  struct handle
  {
    using type = std::tuple<Ts...>;
    type value = {};
    handle(Ts... ts) : value{std::forward<Ts>(ts)...} {};
    handle() : value{} {};
  };

  template<typename T>
  struct place_holder
  {
    using type = std::array<T, 1>;
  };

  namespace concepts
  {
    // T is a tag
    template<typename T, typename Data>
    concept vertex_item_t = requires (T t) { { static_cast<typename Data::vertex_items>(t) }; };
  }


  template<match::attribute Attr>
  static auto item_attribute = [](concepts::tuple_like auto items) constexpr
  {
    auto loop = rec([]<typename Tpl>(auto&& loop, Tpl tpl)
    {
      using tpl_t = std::decay_t<Tpl>;
      using item_t = std::decay_t<decltype(car(tpl))>;

      static_assert(match::item<item_t>);

      if constexpr(match::attribute_of<Attr, item_t>)
      {
        return item_t{};
      }
      else if constexpr (std::tuple_size_v<tpl_t> > 1)
      {
        return loop(cdr(tpl)) ;
      }
      else
      {
        []<bool flag = false>()
          {
            static_assert(flag, "item not found");
          }();
      }
    });

    return loop(items);
  };

  static auto all_attributes =
  rec(
    [](auto&& all_attributes, match::item auto&& t) constexpr
    {
      using type_t = std::decay_t<decltype(t)>;

      // sub items
      if constexpr (match::uses<type_t>)
      {
        if constexpr (match::attributes<type_t>)
        {
          // attributes + sub items attributes
          return
            flatten(
              append(
                instance<typename type_t::attributes>,
                transform(all_attributes, instance<typename type_t::uses>)));
        }
        else
        {
          return
            flatten(transform(all_attributes, instance<typename type_t::uses>));
        }
      }
      // item
      else
      {
        return instance<typename type_t::attributes>;
      }

    });

  static auto all_items =
    rec(
      [](auto&& all_items, match::item auto root_item)
      {
        using type_t = std::decay_t<decltype(root_item)>;

        if constexpr (match::uses<type_t>)
        {
          return cons(root_item,
                      flatten(transform(all_items,
                                        typename type_t::uses{})));
        }
        else
        {
          return gather<type_t>{};
        }
      });

  static auto all_graph_items = [](match::item auto&& t)
    constexpr -> decltype(auto)
  {
    return filter<hold<decltype([]<typename T>(T) { return match::graph_item<T>; })>>
      (all_items(t));
  };

  static auto all_vertex_items = [](match::item auto&& t)
    constexpr -> decltype(auto)
  {
    return filter<hold<decltype([]<typename T>(T) constexpr { return match::vertex_item<T>; })>>
      (all_items(t));
  };



}

#endif
