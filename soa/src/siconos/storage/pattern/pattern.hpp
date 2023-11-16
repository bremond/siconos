#pragma once

#include <functional>
#include <initializer_list>
#include <tuple>
#include <type_traits>
#include <utility>

#include "siconos/storage/ground/ground.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/base_concepts.hpp"
#include "siconos/storage/some/some.hpp"

namespace siconos::storage::pattern {

// namespace some = storage::some;
// namespace ground = storage::ground;

// let rec
// https://stackoverflow.com/questions/2067988/recursive-lambda-functions-in-c11
template <typename F>
struct recursive {
  F f;
  template <typename... Ts>
  constexpr decltype(auto) operator()(Ts&&... ts) const
  {
    return f(std::ref(*this), std::forward<Ts>(ts)...);
  }

  template <typename... Ts>
  constexpr decltype(auto) operator()(Ts&&... ts)
  {
    return f(std::ref(*this), std::forward<Ts>(ts)...);
  }
};

template <typename F>
recursive(F) -> recursive<F>;

static auto const rec = [](auto f) constexpr {
  return recursive{std::move(f)};
};

static auto proj = [](auto& data) constexpr -> decltype(auto) {
  return ([&data](auto&& fun) constexpr -> decltype(auto) {
    return ([&data,
             &fun]<typename... As>(As&&... args) constexpr -> decltype(auto) {
      return fun(std::forward<As>(args)..., data);
    });
  });
};

static auto compose = [](auto&& f, auto&& g) constexpr -> decltype(auto) {
  return [&f, &g]<typename Arg>(Arg arg) constexpr -> decltype(auto) {
    return f(g(std::forward<Arg>(arg)));
  };
};

static auto car = []<typename Tpl>(Tpl tpl) constexpr {
  static_assert(std::tuple_size_v<Tpl> > 0);
  return std::get<0>(tpl);
};

static_assert(std::is_same_v<decltype(car(gather<int, float, char>{})), int>);

static auto cdr = []<typename A0, typename... As>(
                      std::tuple<A0, As...> tpl) constexpr {
  return std::apply(
      [](const A0&, const As&... args) { return std::make_tuple(args...); },
      tpl);
};

static_assert(std::is_same_v<decltype(cdr(gather<int, float, char>{})),
                             gather<float, char>>);

static auto cons = []<typename A, typename... As>(
                       A, std::tuple<As...>) constexpr {
  return std::tuple<A, As...>{};
};

static_assert(std::is_same_v<decltype(cons(int{}, gather<float, char>{})),
                             gather<int, float, char>>);

static auto append = []<concepts::tuple_like... Tpls>(Tpls... tpls) constexpr
{
  return std::tuple_cat(tpls...);
};

static auto flatten = []<concepts::tuple_like Tpl>(Tpl tpl) constexpr {
  if constexpr (std::tuple_size_v < Tpl >> 0) {
    if constexpr (concepts::tuple_like<decltype(car(tpl))>) {
      return std::apply(append, tpl);
    }
    else {
      return tpl;
    }
  }
  else {
    return tpl;
  }
};

static auto transform = []<typename... Args>(
                            auto&& fun, const std::tuple<Args...>) constexpr {
  return std::tuple{fun(Args{})...};
};

template <typename F>
static auto filter = []<typename Tpl>(const Tpl tpl) constexpr {
  auto loop = rec([]<typename Itpl>(auto&& loop, Itpl) constexpr {
    auto v = car(Itpl{});
    using v_t = std::decay_t<decltype(v)>;
    if constexpr (F{}.value(v_t{})) {
      if constexpr (std::tuple_size_v < Itpl >> 1) {
        return cons(v, loop(cdr(Itpl{})));
      }
      else {
        return std::tuple<v_t>(v);
      }
    }
    else if constexpr (std::tuple_size_v<Itpl> > 1) {
      return loop(cdr(Itpl{}));
    }
    else {
      return std::tuple<>{};
    }
  });

  if constexpr (std::tuple_size_v < Tpl >> 0) {
    return loop(Tpl{});
  }
  else {
    return std::tuple<>{};
  }
};

template <std::size_t N>
struct size {
  using size_t = void;
  static constexpr std::size_t value = N;
};

template <std::size_t N>
struct degrees_of_freedom : size<N> {
  using degrees_of_freedom_t = void;
};

template <typename T>
static constexpr auto instance = T{};

template <typename T>
struct hold {
  static constexpr auto value = T{};
};
template <typename T>
static constexpr auto make_instance(T&&)
{
  return instance<T>;
};

template <typename U>
static auto contains =
    []<typename... Attrs>(gather<Attrs...>) constexpr -> bool {
  return (std::is_same_v<U, Attrs> || ...);
};

namespace must {
template <typename T, typename Tpl>
concept contains = contains<T>(instance<Tpl>);
}

static_assert(match::size<std::vector<double>>);
static_assert(match::size<std::array<double, 1>>);
static_assert(match::push_back<std::vector<double>>);
static_assert(!match::push_back<std::array<double, 3>>);

namespace match {

template <typename T>
concept type_t = requires { typename T::type_t; };

template <typename T>
concept attribute = requires { typename T::attribute_t; };

template <typename T>
concept attribute_with_internal_type =
    attribute<T> && requires { typename T::type; };

template <typename T>
concept attribute_or_item = attribute<T> || item<T>;

template <typename T, typename I>
concept attached_storage =
    attribute<T> && item<I> &&
    (std::derived_from<I, typename T::item> ||
     std::derived_from<typename T::item, I>
     // wrap case
     || std::derived_from<typename I::type, typename T::item>);

template <typename T, typename Tag>
concept tag = std::derived_from<typename T::tag, Tag>;

template <typename T>
concept item_ref = attribute<T> && item<typename T::type>;

template <typename T>
concept unbounded_storage = std::derived_from<T, some::unbounded_storage>;

template <typename T>
concept bounded_storage = std::derived_from<T, some::bounded_storage>;

template <typename T>
concept abstract_vector =
    std::derived_from<T, some::undefined_vector> ||
    std::derived_from<T, some::undefined_unbounded_vector>;

template <typename T>
concept abstract_matrix =
    std::derived_from<T, some::undefined_matrix> ||
    std::derived_from<T, some::undefined_unbounded_matrix>;

template <typename T, std::size_t N>
concept cvector = requires(T a) { a[N - 1]; };

template <typename T>
concept property = std::derived_from<T, some::property>;

template <typename T, typename K>
concept property_of = property<T> && property<K> && std::derived_from<T, K>;

template <typename T, typename Ks>
concept any_of_property = ground::any_of(Ks{}, []<match::property K>(K) {
  return std::derived_from<T, K>;
});

}  // namespace match

template <typename T>
struct tag {
  using type = T;
};

// concept item
template <match::item T>
struct use {
  using use_t = void;
  using type = T;
};

template <typename T>
struct attribute_p {
  static constexpr auto value = match::attribute<T>;
};

template <typename T>
struct degrees_of_freedom_p {
  static constexpr auto value = match::degrees_of_freedom<T>;
};

template <typename... Args>
struct frame {
  using args = gather<Args...>;

  static constexpr std::size_t dof =
      ground::find_if(args{},
                      []<typename T>(T) { return degrees_of_freedom_p<T>{}; })
          .value_or([]<bool flag = false>() {
            static_assert(flag, "need some dof");
          })
          .value;
};

template <typename... Args>
struct item {
  using item_t = void;

  using args = gather<Args...>;
  using attributes = decltype(ground::filter(
      args{}, ground::derive_from<some::attribute<>>));

  template <typename H>
  struct interface {
    decltype(auto) self()
    {
      return static_cast<H*>(this);  // handle inherits from default_interface
    }
  };

  friend auto operator<=>(const item<Args...>&,
                          const item<Args...>&) = default;
};

struct any_wrapper {};

struct any_bounded_wrapper : any_wrapper {};
struct any_unbounded_wrapper : any_wrapper {};

template <template <typename... Ts> typename Wrapper, match::item Item,
          typename... Args>
struct wrap : Wrapper<Item, Args...>, Item, any_wrapper {
  using wrap_t = void;
  template <typename T>
  using wrapper = Wrapper<T, Args...>;
  //  using attributes = typename Item::attributes;
  using type = Item;
};

template <typename T>
struct place_holder {
  using type = std::array<T, 1>;
};

namespace concepts {
// T is a tag
template <typename T, typename Data>
concept vertex_item_t = requires(T t) {
                          {
                            static_cast<typename Data::vertex_items>(t)
                          };
                        };
}  // namespace concepts

template <string_literal Name, match::attribute A>
struct attribute : A, symbol<Name> {};

template <match::item Item, match::attribute A>
struct paired : A {
  using item = Item;
};

template <match::item Item, match::attribute A>
static constexpr decltype(auto) named_attribute_maybe(Item, A)
{
  if constexpr (std::derived_from<A, any_symbol>) {
    // for an attribute declared inline inside attributes type, attach
    // Item
    if constexpr (match::wrap<Item>) {
      return paired<typename Item::type, A>{};
    }
    else {
      static_assert(!match::wrap<Item>);
      return paired<Item, A>{};
    }
  }
  else {
    return A{};
  }
}

static auto attributes =
    []<match::item Item>(Item) constexpr -> decltype(auto) {
  if constexpr (match::attributes<Item>) {
    return ground::transform(typename Item::attributes{},
                             []<match::attribute A>(A) {
                               return named_attribute_maybe(Item{}, A{});
                             });
  }
  else {
    return gather<>{};
  }
};

template <typename Attrs, string_literal S>
using get_attr_t = std::decay_t<decltype(std::get<0>(
    ground::filter(Attrs{}, ground::derive_from<symbol<S>>)))>;

template <match::item Item, string_literal S>
using attr_t = std::decay_t<decltype(std::get<0>(
    ground::filter(attributes(Item{}), ground::derive_from<symbol<S>>)))>;

static auto properties =
    []<match::item Item>(Item) constexpr -> decltype(auto) {
  if constexpr (match::properties<Item>) {
    return typename Item::properties{};
  }
  else {
    return gather<>{};
  }
};

static auto all_items = rec([](auto&& all_items, match::item auto root_item) {
  using type_t = std::decay_t<decltype(root_item)>;

  auto items_ref = []() {
    if constexpr (match::attributes<type_t>) {
      return transform([]<typename T>(T) { return typename T::type{}; },
                       filter<hold<decltype([]<typename T>(T) {
                         return match::item_ref<T>;
                       })>>(attributes(type_t{})));
    }
    else {
      return gather<>{};
    }
  };

  auto poly_ref = []() {
    if constexpr (match::attributes<type_t>) {
      return transform(
          []<typename T>(T) { return typename T::type{}; },
          flatten(transform([]<typename T>(T) { return typename T::types{}; },
                            filter<hold<decltype([]<typename T>(T) {
                              return match::polymorphic_type<T>;
                            })>>(attributes(type_t{})))));
    }
    else {
      return gather<>{};
    }
  };

  if constexpr (match::items<type_t>) {
    return cons(
        root_item,
        flatten(transform(all_items,
                          append(items_ref(), typename type_t::items{}))));
  }
  else {
    return cons(root_item, flatten(transform(
                               all_items, append(items_ref(), poly_ref()))));
    ;
  }
});

static auto all_attributes = []<match::item Item>(Item) constexpr {
  return flatten(append(transform(attributes, all_items(Item{}))));
};

static auto all_properties = []<match::item Item>(Item) constexpr {
  return flatten(append(transform(properties, all_items(Item{}))));
};

//  template<typename K>
//  static auto all_items_of_property = [](match::item auto&& t)
//    constexpr -> decltype(auto)
//  {
//    return filter<hold<decltype([]<typename T>(T) { return
//    match::property<T,K>; })>>
//      (all_items(t));
//  };

namespace match {

template <typename H, typename A>
concept handle_attribute =
    attribute<A> && item<typename H::type> &&
    must::contains<A, decltype(siconos::storage::pattern::attributes(
                          typename H::type{}))>;

template <typename H, typename A>
concept handle_attached_storage =
    item<typename H::type> &&
    std::tuple_size_v<
        std::decay_t<decltype(filter<hold<decltype([]<typename T>(T) {
          return (match::attached_storage<T, typename H::type> &&
                  match::tag<T, A>);
        })>>(typename H::info_t::all_properties_t{}))>> >=
        1;  // not an attached storage

template <typename T, typename I>
concept attribute_of =
    attribute<T> && item<I> &&
    must::contains<T, decltype(siconos::storage::pattern::attributes(I{}))>;
}  // namespace match
namespace types {
template <template <typename T> typename Transform, typename... Args>
using transform = decltype(transform(
    []<typename A>(A) { return Transform<A>{}; }, gather<Args...>{}));
template <match::attribute... Attrs>
using attributes = gather<Attrs...>;

template <match::item... Items>
using attributes_of_items =
    decltype(flatten(append(typename Items::attributes{}...)));

template <match::item... Items>
using properties_of_items =
    decltype(flatten(append(typename Items::properties{}...)));

template <size_t N, typename tpl>
using nth_t = std::decay_t<decltype(std::get<N>(tpl{}))>;
}  // namespace types

template <string_literal S>
struct indice_value : symbol<S> {
  std::size_t value;
};

template <string_literal S>
struct param : symbol<S> {};

template <auto V>
struct param_val {
  static constexpr auto value = V;
};

template <typename T>
struct param_type {
  using type = T;
};

template <match::attribute Attr>
static auto item_attribute = [](concepts::tuple_like auto items) constexpr {
  using items_t = std::decay_t<decltype(items)>;

  auto loop = rec([]<typename Tpl>(auto&& loop, Tpl tpl) {
    using tpl_t = std::decay_t<Tpl>;
    using item_t = std::decay_t<decltype(car(tpl))>;

    static_assert(match::item<item_t>);

    if constexpr (match::attribute_of<Attr, item_t>) {
      return item_t{};
    }
    else if constexpr (match::attached_storage<Attr, item_t>) {
      return item_t{};
    }
    else if constexpr (std::tuple_size_v<tpl_t> > 1) {
      return loop(cdr(tpl));
    }
    else {
      []<typename Attribute = Attr, typename LastItem = item_t,
         typename Items = items_t, bool flag = false>()
      {
        static_assert(flag, "item not found");
      }
      ();
    }
  });

  return loop(items);
};

static auto constexpr attribute_name(match::attribute auto a)
{
  return a.str.value;
};

}  // namespace siconos::storage::pattern
