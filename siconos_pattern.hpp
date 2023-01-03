#ifndef SICONOS_PATTERN_HPP
#define SICONOS_PATTERN_HPP

#include "siconos_util.hpp"
#include "siconos_ground.hpp"

#include <tuple>
#include <utility>
#include <functional>
#include <type_traits>
#include <initializer_list>

namespace siconos
{

  template<typename ...Args>
  using gather = std::tuple<Args...>;

  struct empty{};
  struct linear{};
  struct time_invariant{};


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


  static auto car = []<typename Tpl>(Tpl tpl)
    constexpr
  {
    static_assert (std::tuple_size_v<Tpl> > 0);
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
    if constexpr (std::tuple_size_v<Tpl> > 0)
    {
      if constexpr (concepts::tuple_like<decltype(car(tpl))>)
      {
        return std::apply(append, tpl);
      }
      else
      {
        return tpl;
      }
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

  template<std::size_t N>
  struct degrees_of_freedom
  {
    using degrees_of_freedom_t = void;
    static constexpr std::size_t value = N;
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
    concept scalar = std::is_scalar_v<T>;

    template<typename T>
    concept indice = std::is_scalar_v<T> &&
      requires (T i) { std::array<double,1>{}[i]; };

    template<typename T>
    concept degrees_of_freedom = requires { typename T::degrees_of_freedom_t; };

    template<typename T>
    concept attributes = requires { typename T::attributes; };

    template<typename T>
    concept vdescriptor = requires { typename T::vdescriptor_t; };

    template<typename T>
    concept item = requires { typename T::item_t; };

    template<typename T>
    concept items =
      item<T> &&
      requires { typename T::items; };

    template<typename T>
    concept properties =
      item<T> &&
      requires { typename T::properties; };

    template<typename T>
    concept size = requires (T a) { { std::size(a) }; };

    template<typename T>
    concept push_back = requires (T a) { { a.push_back(typename T::value_type{}) }; };

    template<typename T, typename I>
    concept handle =
      item<I> && std::derived_from<I, typename T::type> &&
      requires { typename T::handle_t; };

    template<typename T>
    concept any_full_handle =
      requires { typename T::full_handle_t; };

    template<typename T, typename I>
    concept full_handle =
      any_full_handle<T> &&
      item<I> && std::derived_from<I, typename T::type>;



  }

  static_assert(match::size<std::vector<double>>);
  static_assert(match::size<std::array<double,1>>);
  static_assert(match::push_back<std::vector<double>>);
  static_assert(!match::push_back<std::array<double,3>>);

  namespace some
  {
    template<typename ...Args>
    struct attribute
    {
      using args = gather<Args...>;
      using attribute_t = void;
    };

    struct property
    {
      using property_t = void;
    };

    struct attached_storage : property {};

    // not here
    struct time_invariant : property
    {
      using time_invariant_t = void;
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

    struct undefined_matrix : attribute<> {};
    struct undefined_diagonal_matrix : undefined_matrix {};
    struct undefined_vector : attribute<> {};

    template<size_t ...Sizes>
    struct with_sizes
    {
       static constexpr auto sizes = std::tuple {Sizes...};
    };

    template<typename Type>
    struct with_type
    {
      using type = Type;
    };

    template<typename ...Types>
    struct with_types
    {
      using types = std::tuple<Types...>;
    };

    template<std::size_t N, std::size_t M, typename Type = some::scalar>
    struct matrix : undefined_matrix,
      with_sizes<N, M>, with_type<Type> {};

    template<typename M, typename Type = some::scalar>
    struct diagonal_matrix : undefined_diagonal_matrix,
      with_sizes<std::get<0>(M::sizes)>, with_type<Type> {};

    template<std::size_t N, typename Type = some::scalar>
    struct vector : undefined_vector, with_sizes<N>, with_type<Type> {};

    struct undefined_graph : attribute<> {};

    template<typename Edge, typename Vertice>
    struct graph : undefined_graph, with_types<Edge, Vertice> {};

    struct undefined_unbounded_collection : attribute<> {};
    struct undefined_bounded_collection : attribute<> {};

    template<typename Type>
    struct unbounded_collection : undefined_unbounded_collection, with_type<Type> {};

    template<typename Type, std::size_t N>
    struct bounded_collection : undefined_bounded_collection, with_sizes<N>, with_type<Type> {};

    template<match::item T>
    struct item_ref : attribute<>
    {
      using type = T;
    };

  }

  namespace match
  {
    template<typename T>
    concept attribute = requires { typename T::attribute_t; };

    template<typename T>
    concept attribute_or_item = attribute<T> || item<T>;

    template<typename T, typename I>
    concept attribute_of =
      attribute<T> &&
      item<I> &&
      must::contains<T, typename I::attributes>;

    template<typename T, typename I>
    concept attached_storage =
      attribute<T> &&
      item<I> &&
      (std::derived_from<I, typename T::item> ||  std::derived_from<typename T::item, I>);

    template<typename T>
    concept item_ref =
      attribute<T> &&
      item<typename T::type>;

    template<typename T>
    concept vector = std::derived_from<T, some::vector<std::get<0>(T::sizes)>>
      || requires (T a) { a[0] ;};

    template<typename T>
    concept matrix = std::derived_from<T, some::matrix<std::get<0>(T::sizes), std::get<1>(T::sizes)>> || requires (T a) { a[0] ;};

    template<typename T, std::size_t N>
    concept cvector = requires (T a) { a[N-1]; };

    template<typename T>
    concept property = std::derived_from<T, some::property>;

    template<typename T, typename K>
    concept property_of = property<T> && property<K> && std::derived_from<T, K>;

    template<typename T, typename Ks>
    concept any_of_property = ground::any_of(
      Ks{},
      []<match::property K>(K) { return std::derived_from<T, K>; });

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

  template<typename T>
  struct degrees_of_freedom_p
  {
    static constexpr auto value = match::degrees_of_freedom<T>;
  };

  template<typename ...Args>
  struct frame
  {
    using args = gather<Args...>;

    static constexpr std::size_t dof =
      ground::find_if(args{},
                      []<typename T>(T) { return degrees_of_freedom_p<T>{}; }).
      value_or([]<bool flag = false>() { static_assert(flag, "need some dof"); }).value;
  };

  template<typename ...Args>
  struct item
  {
    using args = gather<Args...>;
    using item_t = void;

    friend auto operator<=>(const item<Args...>&, const item<Args...>&) = default;
  };

  template<typename Wrap>
  struct wrap : item<>, Wrap
  {
    static_assert(match::item<typename Wrap::type>);
    using attributes = typename Wrap::type::attributes;
    using type = typename Wrap::type;
  };

  template<match::item I, typename Tag, match::attribute DataSpec>
  struct attached_storage : some::attached_storage, DataSpec
  {
    using item    = I;
    using tag      = Tag;
    using data_spec = DataSpec;
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
      else if constexpr (match::attached_storage<Attr, item_t>)
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

  static auto attributes = []<match::item Item>(Item&& t)
    constexpr -> decltype(auto)
  {
    if constexpr (match::attributes<Item>)
    {
      return typename Item::attributes{};
    }
    else
    {
      return gather<>{};
    }
  };

  // deprecated
  static auto all_attributes =
  rec(
    [](auto&& all_attributes, match::item auto&& t) constexpr
    {
      using type_t = std::decay_t<decltype(t)>;

      auto items_ref = []()
      {
        if constexpr (match::attributes<type_t>)
        {
          return transform([]<typename T>(T) { return typename T::type{}; },
                           filter<hold<decltype([]<typename T>(T) { return match::item_ref<T>; })>>(typename type_t::attributes{}));
        }
        else
        {
          return gather<>{};
        }
      };
      // sub items
      if constexpr (match::items<type_t>)
      {
        if constexpr (match::attributes<type_t>)
        {
          // attributes + sub items attributes
          return
            flatten(
              append(
                instance<typename type_t::attributes>,
                transform(all_attributes, append(items_ref(),
                                                 instance<typename type_t::items>))));
        }
        else
        {
          return
            flatten(transform(all_attributes,
                              append(items_ref(),
                                     instance<typename type_t::items>)));
        }
      }
      // item
      else
      {
        if constexpr (match::attributes<type_t>)
        {
          return flatten(append(instance<gather<typename type_t::attributes>>,
                                transform(all_attributes, items_ref())));
        }
        else
        {
          return flatten(transform(all_attributes, items_ref()));
        }
      }

    });

  static auto all_items =
    rec(
      [](auto&& all_items, match::item auto root_item)
      {
        using type_t = std::decay_t<decltype(root_item)>;

        auto items_ref = []()
        {
          if constexpr (match::attributes<type_t>)
          {
            return transform([]<typename T>(T) { return typename T::type{}; },
                             filter<hold<decltype([]<typename T>(T) { return match::item_ref<T>; })>>(typename type_t::attributes{}));
          }
          else
          {
            return gather<>{};
          }
        };
        if constexpr (match::items<type_t>)
        {
          return cons(root_item,
                      flatten(transform(all_items,
                                        append(items_ref(), typename type_t::items{}))));
        }
        else
        {
          return cons(root_item, flatten(transform(all_items,
                                                   items_ref())));;
        }
      });

//  template<typename K>
//  static auto all_items_of_property = [](match::item auto&& t)
//    constexpr -> decltype(auto)
//  {
//    return filter<hold<decltype([]<typename T>(T) { return match::property<T,K>; })>>
//      (all_items(t));
//  };

  namespace match
  {
    //
    template<typename H, typename A>
    concept handle_attribute =
      attribute<A> &&
      item<typename H::type> &&
      must::contains<A, typename H::type::attributes>;

    // suspect
    template<typename H, typename IA>
    concept handle_intern_attribute =
      attribute<IA> &&
      item<typename H::type> &&
      must::contains<IA, decltype(all_attributes(typename H::type::attributes{}))>;
  }
  namespace type
  {
    template<template<typename T> typename Transform, typename ...Args>
    using transform = decltype(
      transform([]<typename A>(A){ return Transform<A>{};},
                gather<Args...>{}));
    template<match::attribute ...Attrs>
    using attributes = gather<Attrs...>;

    template<match::item ...Items>
    using attributes_of_items = decltype(flatten(append(typename Items::attributes{}...)));
  }

  template<typename Handle>
  struct default_interface
  {
    decltype(auto) self() { return static_cast<Handle*>(this); };
  };
}

#endif
